#ifndef HELLFIREOS_APP_H
#include <hellfire.h>
#include <noc.h>
#define HELLFIREOS_APP_H
#endif //HELLFIREOS_APP_H
#include "filter.h"
#include "image.h"


#define WIDTH_IMAGE  128
#define HEIGHT_KERNEL  5
#define CENTER_LINE 2
#define SIZE_CRC    4
#define SIZE_PROC_BUFFER WIDTH_IMAGE * HEIGHT_KERNEL
#define SIZE_COMM_BUFFER WIDTH_IMAGE + SIZE_CRC
#define TIMEOUT 500
#define PORT_SOURCE 2000
#define PORT_TARGET 3000
#define PORT_WORKER 4000
#define CHANNEL_SOURCE 2
#define CHANNEL_TARGET 3
#define CHANNEL_WORKER 4
#define NUM_CPU 1
#define IDX_SENDED 0
#define IDX_RECEIVED 1


void source(void)
{
    uint32_t crc,  x, z;
    int8_t buf[SIZE_COMM_BUFFER];
    uint8_t * ptr;
    int16_t val;
//    int16_t control[NUM_CPU][2];
    uint32_t sended_messages=0;

    uint16_t cpu, port, size;
    int32_t channel;

    ptr = image;

    if (hf_comm_create(hf_selfid(), PORT_SOURCE, 0))
        panic(0xff);

    delay_ms(50);

//    srand(hf_cpuid());
//    printf("Start\n");

    // generate a unique channel number for this CPU
//    channel = hf_cpuid();
    while (1){
        if (sended_messages > WIDTH_IMAGE - 1) break;

        channel = hf_recvprobe();
        if (channel < 0) {
            printf("hf_recvprobe(): error %d\n", channel);
            hf_recv(&cpu, &port, buf, &size, channel);
            continue;
        }

        val = hf_recv(&cpu, &port, buf, &size, channel);
        if (val){
            printf("hf_recv(): error %d\n", val);
            continue;
        }

        printf("S cpu %d, port %d, channel %d [free queue: %d] ", cpu, port, channel, size, crc,
               hf_queue_count(pktdrv_queue));

        crc = hf_crc32((int8_t *)ptr, sizeof(buf) - SIZE_CRC);
        memcpy(ptr + WIDTH_IMAGE, &crc, SIZE_CRC);
        val = hf_send(cpu, PORT_WORKER, (int8_t *) ptr, SIZE_COMM_BUFFER, cpu);

        if (val){
            printf("hf_send(): error %d\n", val);
            continue;
        }

//        for(x = control * CENTER_LINE; x < WIDTH_IMAGE * (CENTER_LINE + 1); x++){
//            printf("0x%02x", image[x]);
//            printf(", ");
//            if ((++z % 16) == 0) printf("\n");
//        }
//        printf("\n");

        sended_messages++;
        ptr = ptr + WIDTH_IMAGE;
        printf("sended_messages %d ", sended_messages);
        printf("ptr %d", ptr);
        printf("\n");
    }

    printf("\n\nend of processing!\n\n");
//    panic(0);
    while (1);
}


void worker(void)
{
    int8_t buf[SIZE_COMM_BUFFER], buf_dummy[1];
    uint16_t cpu, port, size, cpuid;
    int16_t val, x, z = 0;
    uint32_t crc;
    int32_t channel;
    uint32_t control=0, recv_messages=0;

    uint8_t *img_gauss, *img_sobel, *img;
//    uint32_t time;
    cpuid = hf_cpuid();

    img = (uint8_t *) malloc(SIZE_PROC_BUFFER);
    img_sobel = (uint8_t *) malloc(SIZE_PROC_BUFFER);
    img_gauss = (uint8_t *) malloc(SIZE_PROC_BUFFER);

    if (img_gauss == NULL || img_sobel == NULL){
        printf("\nmalloc() failed!\n");
        for(;;);
    }

    if (hf_comm_create(hf_selfid(), PORT_WORKER, 0))
        panic(0xff);

    while (1){
        if (recv_messages > WIDTH_IMAGE - 1) break;

        // request data to source
        val = hf_send(0, PORT_SOURCE, buf_dummy, 0, hf_cpuid());
        if (val) {
            printf("hf_send(): error %d\n", val);
            delay_ms(5);
            continue;
        }

        // receive data from source
        channel = hf_recvprobe();
        if (channel < 0) {
//            hf_recv(&cpu, &port, buf, &size, channel);
            printf("hf_recvprobe(): error %d\n", channel);
            continue;
        }
//        else if (channel == 1) {
//            hf_recv(&cpu, &port, buf, &size, channel);
//            printf("hf_recvprobe(): error %d\n", channel);
//        }

        val = hf_recv(&cpu, &port, buf, &size, channel);
        if (val){
            printf("hf_recv(): error %d\n", val);
            continue;
        }

        memcpy(&crc, buf + size - SIZE_CRC, SIZE_CRC);
        printf("W cpu %d, port %d, channel %d, size %d, crc %08x [free queue: %d]", cpu, port, channel, size, crc,
               hf_queue_count(pktdrv_queue));

        if (hf_crc32(buf, size - SIZE_CRC) != crc) {
            printf(" (CRC32 fail) \n");
            continue;
        }
        printf(" (CRC32 pass) ");

        memmove(img, img + WIDTH_IMAGE, SIZE_PROC_BUFFER - WIDTH_IMAGE);
        memmove(img + (SIZE_PROC_BUFFER - WIDTH_IMAGE), buf, WIDTH_IMAGE);

        recv_messages++;
        printf("recv_messages %d ", recv_messages);
        printf("\n");

        if (recv_messages < HEIGHT_KERNEL) continue;

        do_sobel0((uint8_t *)img, img_sobel, WIDTH_IMAGE, HEIGHT_KERNEL);
//        do_gaussian(img_sobel, img_gauss, WIDTH_IMAGE, HEIGHT_KERNEL);

//        for(x = WIDTH_IMAGE * CENTER_LINE; x < WIDTH_IMAGE * (CENTER_LINE + 1); x++){
//            printf("0x%02x", img[x]);
//            printf(", ");
//            if ((++z % 16) == 0) printf("\n");
//        }


        // send data to target
//        val = hf_send(2, PORT_TARGET, img_gauss + CENTER_LINE * SIZE_PROC_BUFFER,
//                      SIZE_COMM_BUFFER, cpuid);
//        if (val)
//            printf("hf_send(): error %d\n", val);

//        delay_ms(5);

    }

//    for(x = WIDTH_IMAGE * CENTER_LINE; x < SIZE_PROC_BUFFER - 1; x++){
//        printf("0x%02x", img_gauss[x]);
//        printf(", ");
//        if ((++z % 16) == 0) printf("\n");
//    }
    free(img);
    free(img_gauss);
    free(img_sobel);

    printf("\n\nend of processing!\n\n");
    panic(0);

}


void target(void)
{
    int8_t buf[SIZE_COMM_BUFFER];
    uint32_t i, j, k = 0;

    uint16_t cpu, port, size;
    int16_t val;
    uint32_t crc;
    int32_t channel;
    uint8_t *filter_image;
    uint8_t * ptr;
    uint32_t received_messages=0;

    if (hf_comm_create(hf_selfid(), PORT_TARGET, 0))
        panic(0xff);

    filter_image = (uint8_t *) malloc(height * width);
    ptr = filter_image;
    ptr = ptr + WIDTH_IMAGE * CENTER_LINE;

    while (1){
        if (received_messages > WIDTH_IMAGE - 1) break;

        channel = hf_recvprobe();
        if (channel < 0) continue;

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
        memmove(ptr, buf, WIDTH_IMAGE);

        ptr = ptr + WIDTH_IMAGE;

        printf("control %d\n", received_messages);
        received_messages++;
    }

    printf("\n\nint32_t width = %d, height = %d;\n", width, height);
    printf("uint8_t image[] = {\n");
    for (i = 0; i < height; i++){
        for (j = 0; j < width; j++){
            printf("0x%02x", filter_image[i * width + j]);
            if ((i < height-1) || (j < width-1)) printf(", ");
            if ((++k % 16) == 0) printf("\n");
        }
    }
    printf("};\n");

}


void app_main(void)
{
    uint8_t i;
    switch (hf_cpuid()) {
        case 0:
            hf_spawn(source, 0, 0, 0, "s", 4096);
//            hf_spawn(target, 0, 0, 0, "w", 2048);
        case 1:
            hf_spawn(worker, 0, 0, 0, "w", 4096);
//        case 2:
//            hf_spawn(target, 0, 0, 0, "w", 2048);
    }

//    if (hf_cpuid() == 0)
//        hf_spawn(source, 0, 0, 0, "s", 2048);
//
//    for (i = 1; i < hf_ncores(); i++, msg++) {
//        if (hf_cpuid() == i) {
//            hf_spawn(worker, 0, 0, 0, "w", 2048);
//        }
//    }
}
