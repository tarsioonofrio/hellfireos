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
#define CPU_SOURCE 0
#define CPU_TARGET 2
#define NUM_CPU 1
#define MESSAGE_PER_CPU WIDTH_IMAGE / NUM_CPU
#define IDX_SENDED 0
#define IDX_RECEIVED 1


void source(void)
{
    uint32_t crc,  x, z;
    int8_t buf[SIZE_COMM_BUFFER];
    uint8_t * ptr[NUM_CPU], i;
    int16_t val;
//    int16_t control[NUM_CPU];
    uint32_t sended_messages=0;
    printf("MESSAGE_PER_CPU %d\n", MESSAGE_PER_CPU);

    uint16_t cpu, port, size;
    int32_t channel;

    for (i=0; i < NUM_CPU; i++){
        ptr[i] = &image[MESSAGE_PER_CPU * i];
    }

//    ptr[0] = image;

    if (hf_comm_create(hf_selfid(), PORT_SOURCE, 0))
        panic(0xff);

    delay_ms(10);

//    srand(hf_cpuid());
//    printf("Start\n");

    // generate a unique channel number for this CPU
//    channel = hf_cpuid();
    while (1){
        if (sended_messages > WIDTH_IMAGE - 1) break;

        channel = hf_recvprobe();
        if (channel < 0) {
//            printf("hf_recvprobe(): error %d\n", channel);
//            hf_recv(&cpu, &port, buf, &size, channel);
            continue;
        }
        else if (channel == 0){
//            printf("hf_recvprobe(): error %d\n", channel);
            hf_recv(&cpu, &port, buf, &size, channel);
            continue;
        }

        val = hf_recv(&cpu, &port, buf, &size, channel);
        if (val){
            printf("hf_recv(): error %d\n", val);
            continue;
        }
        if ((cpu == CPU_SOURCE) || (cpu == CPU_TARGET)) {
            printf("hf_recv(): error cpu %d\n", cpu);
            continue;
        }

        printf("S cpu %d, port %d, channel %d, size %d, [free queue: %d]", cpu, port, channel, size,
               hf_queue_count(pktdrv_queue));

        crc = hf_crc32((int8_t *)ptr[cpu - 1], sizeof(buf) - SIZE_CRC);
        memcpy(ptr[cpu - 1] + WIDTH_IMAGE, &crc, SIZE_CRC);
        val = hf_send(cpu, PORT_WORKER, (int8_t *) ptr[cpu - 1], SIZE_COMM_BUFFER, 1);

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
        ptr[cpu - 1] = ptr[cpu - 1] + WIDTH_IMAGE;
        printf("sended_messages %d ", sended_messages);
        printf("ptr %d", ptr[cpu - 1]);
        printf("\n");
    }

    free(ptr);
    printf("\n\nend of processing!\n\n");
    while (1);
}


void worker(void)
{
    int8_t buffer_source[SIZE_COMM_BUFFER], buf_dummy[1], buffer_target[SIZE_COMM_BUFFER];
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
        val = hf_send(CPU_SOURCE, PORT_SOURCE, buf_dummy, 1, 1); // cpuid
        if (val) {
            printf("hf_send(): error %d\n", val);
            continue;
        }

        // receive data from source
        channel = hf_recvprobe();
        if (channel < 0) {
            printf("hf_recvprobe(): error %d\n", channel);
//            hf_recv(&cpu, &port, buf, &size, channel);
            continue;
        }
        else if (channel == 0){
            printf("hf_recvprobe(): error %d\n", channel);
            hf_recv(&cpu, &port, buffer_source, &size, channel);
            continue;
        }

        val = hf_recv(&cpu, &port, buffer_source, &size, channel);
        if (val){
            printf("hf_recv(): error %d\n", val);
            continue;
        }

        memcpy(&crc, buffer_source + size - SIZE_CRC, SIZE_CRC);
        printf("W cpu %d, port %d, channel %d, size %d, crc %08x [free queue: %d]", cpu, port, channel, size, crc,
               hf_queue_count(pktdrv_queue));

        if (hf_crc32(buffer_source, size - SIZE_CRC) != crc) {
            printf(" (CRC32 fail) \n");
            continue;
        }
        printf(" (CRC32 pass) ");

        memmove(img, img + WIDTH_IMAGE, SIZE_PROC_BUFFER - WIDTH_IMAGE);
        memmove(img + (SIZE_PROC_BUFFER - WIDTH_IMAGE), buffer_source, WIDTH_IMAGE);

        recv_messages++;
        printf("received_messages %d ", recv_messages);
        printf("\n");

        if (recv_messages < HEIGHT_KERNEL) continue;

        do_sobel0((uint8_t *)img, img_sobel, WIDTH_IMAGE, HEIGHT_KERNEL);
        do_gaussian(img_sobel, img_gauss, WIDTH_IMAGE, HEIGHT_KERNEL);


        memmove(buffer_target, img_gauss + CENTER_LINE * WIDTH_IMAGE, WIDTH_IMAGE);
        crc = hf_crc32((int8_t *)buffer_target, sizeof(buffer_target) - SIZE_CRC);
        memcpy(buffer_target + WIDTH_IMAGE, &crc, SIZE_CRC);

        // send data to target
        val = hf_send(CPU_TARGET, PORT_TARGET, buffer_target,SIZE_COMM_BUFFER, cpuid);
        if (val)
            printf("hf_send(): error %d\n", val);

    }

    free(img);
    free(img_gauss);
    free(img_sobel);

    printf("\n\nend of processing!\n\n");
    while (1);
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
        if (received_messages > WIDTH_IMAGE - 5) break;

        channel = hf_recvprobe();
        if (channel < 0) continue;

        val = hf_recv(&cpu, &port, buf, &size, channel);

        if (val) {
            printf("hf_recv(): error %d\n", val);
            continue;
        }
        if ((cpu == CPU_SOURCE) || (cpu == CPU_TARGET)) {
            printf("hf_recv(): error cpu %d\n", cpu);
            continue;
        }

        memcpy(&crc, buf + size - SIZE_CRC, SIZE_CRC);
        printf("T cpu %d, port %d, channel %d, size %d, crc %08x [free queue: %d]", cpu, port, channel, size, crc,
               hf_queue_count(pktdrv_queue));

        if (hf_crc32(buf, size - SIZE_CRC) != crc) {
            printf(" (CRC32 fail)\n");
            continue;
        }
        printf(" (CRC32 pass) ");

        memmove(ptr, buf, WIDTH_IMAGE);

        ptr = ptr + WIDTH_IMAGE;

        received_messages++;
        printf("received_messages %d\n", received_messages);
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
    free(filter_image);
    free(ptr);
    panic(0);

}


void app_main(void)
{
    uint8_t i;
    switch (hf_cpuid()) {
        case 0:
            hf_spawn(source, 0, 0, 0, "S", 4096);
        case 1:
            hf_spawn(worker, 0, 0, 0, "W", 4096);
        case 2:
            hf_spawn(target, 0, 0, 0, "T", 4096);
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
