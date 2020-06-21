#ifndef HELLFIREOS_APP_H
#include <hellfire.h>
#include <noc.h>
#define HELLFIREOS_APP_H
#endif //HELLFIREOS_APP_H
#include "filter.h"
#include "image.h"


#define SIDE_IMAGE  256
#define SIDE_SOBEL  5
#define SIDE_GAUSS  7
#define SIZE_CRC    4
#define SIDE_MULTIPLE  1
#define SIZE_BUFFER SIDE_IMAGE * SIDE_MULTIPLE + SIZE_CRC
#define TIMEOUT 500

_Noreturn void master(void)
{
    int32_t i;
    uint32_t crc;
    int8_t buf[SIZE_BUFFER];
    uint8_t * ptr;
    int16_t val, channel;

    ptr = image;

    if (hf_comm_create(hf_selfid(), 1000, 0))
        panic(0xff);

    delay_ms(50);

    srand(hf_cpuid());
    printf("Start\n");

    // generate a unique channel number for this CPU
    channel = hf_cpuid();
    while (1){
//        for (i = 0; i < sizeof(buf) - 4; i++) {
//            buf[i] = *(ptr++);
//        }
        crc = hf_crc32((int8_t *)ptr, sizeof(buf) - 4);
        memcpy(ptr + SIZE_BUFFER - 4, &crc, 4);
        val = hf_send(1, 2000, (int8_t *) ptr, SIZE_BUFFER, 1);
//        val = hf_sendack(1, 2000, (int8_t *) ptr, SIZE_BUFFER, 1, TIMEOUT);
        printf("sizeof %d\n", sizeof(buf));
//        if (val) printf("hf_send(): error %d\n", val);
        if (val) printf("sender, hf_sendack(): error %d\n", val);

        ptr = ptr + SIZE_BUFFER - 4;
        printf("ptr %d\n", ptr);
        delay_ms(10);
    }
}


_Noreturn void receiver(void)
{
    int8_t buf[SIZE_BUFFER];
    uint16_t cpu, port, size;
    int16_t val;
    uint32_t crc;
    int32_t i;

    if (hf_comm_create(hf_selfid(), 2000, 0))
        panic(0xff);

    while (1){
        i = hf_recvprobe();
        if (i < 0) continue;

        val = hf_recv(&cpu, &port, buf, &size, i);
//        val = hf_recvack(&cpu, &port, buf, &size, i);
//        if (val) printf("Master, hf_sendack(): error %d\n", val);

        if (val) {
            printf("hf_recv(): error %d\n", val);
            continue;
        }

        memcpy(&crc, buf + size - 4, 4);
        printf("R cpu %d, port %d, channel %d, size %d, crc %08x [free queue: %d]", cpu, port, i, size, crc,
               hf_queue_count(pktdrv_queue));

        if (hf_crc32(buf, size - 4) != crc) {
            printf(" (CRC32 fail)\n");
            continue;
        }

        printf(" (CRC32 pass)\n");
//        val = hf_send(1, 3000, buf, size, 1);
//        if (val)
//            printf("hf_send(): error %d\n", val);

    }
}


_Noreturn void worker(void)
{
    int8_t buf[SIZE_BUFFER * SIDE_SOBEL];
    uint16_t cpu, port, size, recv_messages=0;
    int16_t val;
    uint32_t crc;
    int32_t channel;

    uint8_t *img_gauss, *img_sobel;
//    uint32_t time;

    img_sobel = (uint8_t *) malloc(SIDE_IMAGE * SIDE_SOBEL);
//    img_gauss = (uint8_t *) malloc(SIDE_IMAGE * SIDE_GAUSS);
//    if (img_gauss == NULL || img_sobel == NULL){
//        printf("\nmalloc() failed!\n");
//        for(;;);
//    }


    if (hf_comm_create(hf_selfid(), 3000, 0))
        panic(0xff);

    while (1){
        channel = hf_recvprobe();
        if (channel < 0) continue;

        val = hf_recv(&cpu, &port, buf, &size, channel);
        if (val){
            printf("hf_recv(): error %d\n", val);
            continue;
        }

        memcpy(&crc, buf + size - 4, 4);
        printf("W cpu %d, port %d, channel %d, size %d, crc %08x [free queue: %d]", cpu, port, channel, size, crc,
               hf_queue_count(pktdrv_queue));

        if (hf_crc32(buf, size - 4) != crc) {
            printf(" (CRC32 fail)\n");
            continue;
        }
        printf(" (CRC32 pass)\n");

        memmove(img_sobel + ((recv_messages % SIDE_SOBEL) * SIZE_BUFFER), buf, sizeof(buf)- 4);

        recv_messages++;
        printf("recv_messages %d\n", recv_messages);

        if (recv_messages < SIDE_SOBEL) continue;
        printf("*** do_sobel ***\n");

//        do_sobel((uint8_t *) buf, img_sobel, SIDE_IMAGE, SIDE_MULTIPLE);
//                for (i = 0; i < sizeof(img_sobel) - 4; i++) {
//                    printf("S %d %x\n", i, img_sobel[i]);
//                }
        printf("\n");

    }
}


void task(void){
    uint32_t i, j, k = 0;
    uint8_t *img_gauss, *img_sobel;
    uint32_t time;

    while(1) {
        img_gauss = (uint8_t *) malloc(height * width);
        img_sobel = (uint8_t *) malloc(height * width);
        if (img_gauss == NULL || img_sobel == NULL){
            printf("\nmalloc() failed!\n");
            for(;;);
        }

        printf("\n\nstart of processing!\n\n");

        time = _readcounter();

        do_gaussian(image, img_gauss, width, height);
        do_sobel(img_gauss, img_sobel, width, height);

        time = _readcounter() - time;

        printf("done in %d clock cycles.\n\n", time);

        printf("\n\nint32_t width = %d, height = %d;\n", width, height);
        printf("uint8_t image[] = {\n");
        for (i = 0; i < height; i++){
            for (j = 0; j < width; j++){
                printf("0x%x", img_sobel[i * width + j]);
                if ((i < height-1) || (j < width-1)) printf(", ");
                if ((++k % 16) == 0) printf("\n");
            }
        }
        printf("};\n");

        free(img_gauss);
        free(img_sobel);

        printf("\n\nend of processing!\n");
        panic(0);
    }

}


_Noreturn void sender(void)
{
    int32_t i;
    uint32_t crc;
    int8_t buf[SIZE_BUFFER];
    uint8_t * ptr;
    int16_t val, channel;

    ptr = image;

    if (hf_comm_create(hf_selfid(), 4000, 0))
        panic(0xff);

    delay_ms(50);

    srand(hf_cpuid());

    // generate a unique channel number for this CPU
    channel = hf_cpuid();
    while (1){
//        for (i = 0; i < sizeof(buf) - 4; i++) {
//            buf[i] = *(ptr++);
//        }
//        crc = hf_crc32((int8_t *)ptr, sizeof(buf) - 4);
//        memcpy(ptr + BUFFER_SIZE - 4, &crc, 4);
        val = hf_send(1, 1000, (int8_t *) ptr, SIZE_BUFFER, 1);
        printf("sizeof %d\n", sizeof(buf));
        if (val)
            printf("hf_send(): error %d\n", val);
        else
            ptr = ptr + SIZE_BUFFER - 4;
    }
}


void app_main(void)
{
    printf(" OI\n");

    switch (hf_cpuid()) {
        case 0:
            hf_spawn(master, 0, 0, 0, "M", 4096);
        case 1:
            hf_spawn(receiver, 0, 0, 0, "R", 4096);
//            hf_spawn(worker, 0, 0, 0, "W", 4096);
    }
}
