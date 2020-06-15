#include <hellfire.h>
#include <noc.h>

#include "filter.h"
#include "image.h"

#define BUFFER_SIZE 768

_Noreturn void sender(void)
{
    int32_t i;
    uint32_t crc;
    int8_t buf[BUFFER_SIZE];
    uint8_t * ptr;
    int16_t val, channel;

    ptr = image;

    if (hf_comm_create(hf_selfid(), 1000, 0))
        panic(0xff);

    delay_ms(50);

    srand(hf_cpuid());

    // generate a unique channel number for this CPU
    channel = hf_cpuid();
    while (1){
//        for (i = 0; i < sizeof(buf) - 4; i++) {
//            buf[i] = *(ptr++);
//        }
        crc = hf_crc32((int8_t *)ptr, sizeof(buf) - 4);
        memcpy(ptr + BUFFER_SIZE - 4, &crc, 4);
        val = hf_send(1, 5000, (int8_t *) ptr, BUFFER_SIZE, 1);
        printf("sizeof %d\n", sizeof(buf));
        if (val)
            printf("hf_send(): error %d\n", val);
        else
            ptr = ptr + BUFFER_SIZE - 4;
    }
}

_Noreturn void receiver(void)
{
    int8_t buf[1500];
    uint16_t cpu, port, size;
    int16_t val;
    uint32_t crc;
    int32_t i;

    if (hf_comm_create(hf_selfid(), 5000, 0))
        panic(0xff);

    while (1){
        i = hf_recvprobe();
        if (i < 0) continue;

        val = hf_recv(&cpu, &port, buf, &size, i);
        if (val){
            printf("hf_recv(): error %d\n", val);
        } else {
            memcpy(&crc, buf + size - 4, 4);
            printf("cpu %d, port %d, channel %d, size %d, crc %08x [free queue: %d]", cpu, port, i, size, crc,
                   hf_queue_count(pktdrv_queue));

            if (hf_crc32(buf, size - 4) == crc) {
                printf(" (CRC32 pass)\n");
                for (i = 0; i < sizeof(buf) - 4; i++) {
                    printf("%02x ", (uint8_t *)buf[i]);
                }
                printf("\n");
            }
           else {
                printf(" (CRC32 fail)\n");
            }
        }
    }
}

void app_main(void)
{
    switch (hf_cpuid()) {
        case 0:
            hf_spawn(sender, 0, 0, 0, "sender", 4096);
        case 1:
            hf_spawn(receiver, 0, 0, 0, "processor", 4096);
    }
}
