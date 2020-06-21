#ifndef HELLFIREOS_APP_H
#include <hellfire.h>
#include <noc.h>
#define HELLFIREOS_APP_H
#endif //HELLFIREOS_APP_H
#include "filter.h"
#include "image.h"


#define WIDTH_FILTER  256
#define HEIGHT_FILTER  5
#define CENTER_LINE 2
#define SIZE_FILTER WIDTH_FILTER * HEIGHT_FILTER
#define SIZE_CRC    4
#define SIZE_BUFFER WIDTH_FILTER + SIZE_CRC
#define TIMEOUT 500
#define PORT_SENDER 1000
#define PORT_RECEIVER 2000
#define PORT_WORKER 3000
#define NUM_CPU 1
#define IDX_SENDED 0
#define IDX_RECEIVED 1


_Noreturn void master_send(void)
{
    uint32_t crc;
    int8_t buf[SIZE_BUFFER];
    uint8_t * ptr;
    int16_t val, channel;
    int16_t control[NUM_CPU][2];

    ptr = image;

    if (hf_comm_create(hf_selfid(), PORT_SENDER, 0))
        panic(0xff);

    delay_ms(50);

    srand(hf_cpuid());
    printf("Start\n");

    // generate a unique channel number for this CPU
    channel = hf_cpuid();
    while (1){
        crc = hf_crc32((int8_t *)ptr, sizeof(buf) - SIZE_CRC);
        memcpy(ptr + SIZE_BUFFER - SIZE_CRC, &crc, SIZE_CRC);
        val = hf_send(1, PORT_RECEIVER, (int8_t *) ptr, SIZE_BUFFER, 1);
        if (val) printf("hf_send(): error %d\n", val);
//        val = hf_sendack(1, 2000, (int8_t *) ptr, SIZE_BUFFER, 1, TIMEOUT);
//        if (val) printf("sender, hf_sendack(): error %d\n", val);

        ptr = ptr + SIZE_BUFFER - SIZE_CRC;
        printf("ptr %d\n", ptr);
        delay_ms(10);
    }
}



_Noreturn void receiver(void)
{
    int8_t buf[SIZE_BUFFER];
    uint16_t cpu, port, size, cpuid;
    int16_t val;
    uint32_t crc;
    int32_t i;

    if (hf_comm_create(hf_selfid(), PORT_RECEIVER, 0))
        panic(0xff);
    cpuid = hf_cpuid();

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

        memcpy(&crc, buf + size - SIZE_CRC, SIZE_CRC);
//        printf("R cpu %d, port %d, channel %d, size %d, crc %08x [free queue: %d]", cpu, port, i, size, crc,
//               hf_queue_count(pktdrv_queue));

        if (hf_crc32(buf, size - SIZE_CRC) != crc) {
            printf(" (CRC32 fail)\n");
            continue;
        }

        printf(" (CRC32 pass)\n");
        val = hf_send(cpuid, PORT_WORKER, buf, size, 1);
        if (val)
            printf("hf_send(): error %d\n", val);

        delay_ms(10);

    }
}


_Noreturn void worker(void)
{
    int8_t buf[SIZE_BUFFER];
    uint16_t cpu, port, size, cpuid, recv_messages=0;
    uint16_t i, j, k;
    int16_t val;
    uint32_t crc;
    int32_t channel;

    uint8_t *img_gauss, *img_sobel, *img;
//    uint32_t time;
    cpuid = hf_cpuid();

    img_sobel = (uint8_t *) malloc(SIZE_FILTER);
    img = (uint8_t *) malloc(SIZE_FILTER);
    img_gauss = (uint8_t *) malloc(SIZE_FILTER);
    if (img_gauss == NULL || img_sobel == NULL){
        printf("\nmalloc() failed!\n");
        for(;;);
    }


    if (hf_comm_create(hf_selfid(), PORT_WORKER, 0))
        panic(0xff);

    while (1){
        channel = hf_recvprobe();
        if (channel < 0) continue;

        val = hf_recv(&cpu, &port, buf, &size, channel);
        if (val){
            printf("hf_recv(): error %d\n", val);
            continue;
        }

        memcpy(&crc, buf + size - SIZE_CRC, SIZE_CRC);
        printf("W cpu %d, port %d, channel %d, size %d, crc %08x [free queue: %d]", cpu, port, channel, size, crc,
               hf_queue_count(pktdrv_queue));

        if (hf_crc32(buf, size - SIZE_CRC) != crc) {
            printf(" (CRC32 fail)\n");
            continue;
        }
        printf(" (CRC32 pass)\n");

        memmove(img, img + WIDTH_FILTER, sizeof(img) - WIDTH_FILTER - SIZE_CRC);

        memmove(img + (SIZE_FILTER - WIDTH_FILTER), buf, SIZE_BUFFER - SIZE_CRC);

        recv_messages++;
        printf("recv_messages %d\n", recv_messages);

        if (recv_messages < HEIGHT_FILTER) continue;

        do_sobel2((uint8_t *) img, img_sobel, WIDTH_FILTER, HEIGHT_FILTER);
        printf("do_sobel\n");

//        for (i = 0; i < SIDE_IMAGE * SIDE_SOBEL; i++) {
//            printf("0x%02x ", img_sobel[i]);
//            if ((i % SIDE_IMAGE) == 0) printf("\n");
//        }
//        for (i = 0; i < HEIGHT_FILTER; i++){
//            for (j = 0; j < WIDTH_FILTER; j++){
//                printf("0x%02x", img_sobel[i * width + j]);
//                if ((i < height-1) || (j < width-1)) printf(", ");
//                k++;
//                if ((k % WIDTH_FILTER) == 0) printf("\n");
//            }
//        }

        printf("\n");

        val = hf_send(cpuid, PORT_SENDER, img_sobel + CENTER_LINE * WIDTH_FILTER, SIZE_BUFFER, 1);
        if (val)
            printf("hf_send(): error %d\n", val);

        delay_ms(10);

        break;

    }
}


_Noreturn void sender(void)
{
    int8_t buf[SIZE_BUFFER];
    uint16_t cpu, port, size, cpuid;
    int16_t val;
    uint32_t crc;
    int32_t i;

    if (hf_comm_create(hf_selfid(), PORT_RECEIVER, 0))
        panic(0xff);

    cpuid = hf_cpuid();

    while (1){
        i = hf_recvprobe();
        if (i < 0) continue;

        val = hf_recv(&cpu, &port, buf, &size, i);

        if (val) {
            printf("hf_recv(): error %d\n", val);
            continue;
        }

        memcpy(&crc, buf + size - SIZE_CRC, SIZE_CRC);
//        printf("R cpu %d, port %d, channel %d, size %d, crc %08x [free queue: %d]", cpu, port, i, size, crc,
//               hf_queue_count(pktdrv_queue));

        if (hf_crc32(buf, size - SIZE_CRC) != crc) {
            printf(" (CRC32 fail)\n");
            continue;
        }

        printf(" (CRC32 pass)\n");
        val = hf_send(0, PORT_RECEIVER, buf, size, cpuid);
        if (val)
            printf("hf_send(): error %d\n", val);

        delay_ms(10);

    }
}


_Noreturn void master_receiver(void)
{
    int8_t buf[SIZE_BUFFER];
    uint16_t cpu, port, size, channel;
    int16_t val;
    uint32_t crc;

    if (hf_comm_create(hf_selfid(), PORT_RECEIVER, 0))
        panic(0xff);

    while (1){
        i = hf_recvprobe();
        if (i < 0) continue;

        val = hf_recv(&cpu, &port, buf, &size, channel);

        if (val) {
            printf("hf_recv(): error %d\n", val);
            continue;
        }

        memcpy(&crc, buf + size - SIZE_CRC, SIZE_CRC);
//        printf("R cpu %d, port %d, channel %d, size %d, crc %08x [free queue: %d]", cpu, port, i, size, crc,
//               hf_queue_count(pktdrv_queue));

        if (hf_crc32(buf, size - SIZE_CRC) != crc) {
            printf(" (CRC32 fail)\n");
            continue;
        }
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



void app_main(void)
{
    printf(" OI\n");

    switch (hf_cpuid()) {
        case 0:
            hf_spawn(master_send, 0, 0, 0, "S", 4096);
            hf_spawn(master_receiver, 0, 0, 0, "R", 4096);
        case 1:
            hf_spawn(receiver, 0, 0, 0, "R", 4096);
            hf_spawn(worker, 0, 0, 0, "W", 4096);
            hf_spawn(sender, 0, 0, 0, "S", 4096);
    }
}
