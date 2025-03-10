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
#define CPU_TARGET 5
#define NUM_CPU 3
#define MESSAGE_PER_CPU WIDTH_IMAGE * WIDTH_IMAGE / NUM_CPU
#define LINES_PER_CPU WIDTH_IMAGE / NUM_CPU
#define IDX_SENDED 0
#define IDX_RECEIVED 1


void source(void)
{
    uint8_t i;
    uint8_t * ptr;
    uint16_t cpu, port, size;
    uint32_t crc, sended_messages=0;

    int8_t buf[SIZE_COMM_BUFFER];
    int16_t val;
    int32_t shift;

    printf("cpu %d, name %s, thread %d.\n", hf_cpuid(), hf_selfname(), hf_selfid());


    if (hf_comm_create(hf_selfid(), PORT_SOURCE, 0))
        panic(0xff);

//    delay_ms(10);


    // generate a unique channel number for this CPU
//    channel = hf_cpuid();
    while (1){
//        if (sended_messages > WIDTH_IMAGE - 1) break;

        shift = hf_recvprobe();
        if (shift < 0) {
            printf("hf_recvprobe(): error %d\n", shift);
            continue;
        }
//        else if (channel == 0){
//            printf("hf_recvprobe(): error %d\n", channel);
//            hf_recv(&cpu, &port, buf, &size, channel);
//            continue;
//        }
        val = hf_recv(&cpu, &port, buf, &size, shift);
        if (val){
            printf("hf_recv(): error %d\n", val);
            continue;
        }
        if ((cpu == CPU_SOURCE) || (cpu == CPU_TARGET)) {
            printf("hf_recv(): error cpu %d\n", cpu);
            continue;
        }

        printf("S cpu %d, port %d, channel %d, size %d, [free queue: %d] ", cpu, port, shift, size,
               hf_queue_count(pktdrv_queue));
//        ptr = image + MESSAGE_PER_CPU * cpu + channel;
        ptr = image + shift;
        crc = hf_crc32((int8_t *)ptr, sizeof(buf) - SIZE_CRC);
//        crc = hf_crc32((int8_t *)ptr[cpu - 1], sizeof(buf) - SIZE_CRC);
        memcpy(ptr + WIDTH_IMAGE, &crc, SIZE_CRC);
        val = hf_send(cpu, port, (int8_t *) ptr, SIZE_COMM_BUFFER, 1);

        if (val){
            printf("hf_send(): error %d\n", val);
            continue;
        }

        sended_messages++;
        printf("sended_messages %d ", sended_messages);
        printf("\n");
    }

    free(ptr);
    printf("\n\nend of processing!\n\n");
    while (1);
}


void worker(void)
{
    uint8_t *img_gauss, *img_sobel, *img;
    uint16_t cpu, port, size, cpuid, shift_source, shift_target;
    uint32_t crc, recv_messages=0, i;

    int8_t buffer_source[SIZE_COMM_BUFFER], buf_dummy[1], buffer_target[SIZE_COMM_BUFFER];
    int16_t val;
    int32_t channel;

    cpuid = hf_cpuid();
    shift_source = MESSAGE_PER_CPU * (cpuid - 1);

    img = (uint8_t *) malloc(SIZE_PROC_BUFFER);
    img_sobel = (uint8_t *) malloc(SIZE_PROC_BUFFER);
    img_gauss = (uint8_t *) malloc(SIZE_PROC_BUFFER);

    printf("cpu %d, name %s, thread %d.\n", hf_cpuid(), hf_selfname(), hf_selfid());

    if (img_gauss == NULL || img_sobel == NULL){
        printf("\nmalloc() failed!\n");
        for(;;);
    }
    val = hf_comm_create(hf_selfid(), PORT_WORKER + hf_selfid() + cpuid, 0);
//    val = hf_comm_create(hf_selfid(), PORT_WORKER, 0);
    if (val) {
        printf("hf_comm_create error %d\n", val);
        panic(0xff);
    }

    while (1){
        if (recv_messages > LINES_PER_CPU - 1) break;

        // request data to source
        val = hf_send(CPU_SOURCE, PORT_SOURCE, buf_dummy, 1, shift_source);


//        delay_ms(NUM_CPU * 3);
        if (val) {
            printf("hf_send(): error %d\n", val);
            continue;
        }

        for (i = 0; i < 25000; i++) {
            // receive data from source
            channel = hf_recvprobe();
            if (channel >-1)
                break;
        }

        if (channel == 0){
            printf("hf_recvprobe(): error %d\n", channel);
            hf_recv(&cpu, &port, buffer_source, &size, channel);
            continue;
        }

        val = hf_recv(&cpu, &port, buffer_source, &size, channel);
        if (val){
            printf("hf_recv(): error %d\n", val);
            continue;
        }

        if (cpu != CPU_SOURCE) {
            printf("hf_recv(): error cpu %d\n", cpu);
            continue;
        }

        memcpy(&crc, buffer_source + size - SIZE_CRC, SIZE_CRC);
        printf("%s: cpu %d, port %d, channel %d, size %d, crc %08x [free queue: %d] ", hf_selfname(), cpu, port, channel, size, crc,
               hf_queue_count(pktdrv_queue));

        if (hf_crc32(buffer_source, size - SIZE_CRC) != crc) {
            printf(" (CRC32 fail) \n");
            continue;
        }
        printf(" (CRC32 pass) ");

        memmove(img, img + WIDTH_IMAGE, SIZE_PROC_BUFFER - WIDTH_IMAGE);
        memmove(img + (SIZE_PROC_BUFFER - WIDTH_IMAGE), buffer_source, WIDTH_IMAGE);

        recv_messages++;
        printf("%s received_messages %d ", hf_selfname(), recv_messages);
        shift_target = shift_source;
        shift_source = shift_source + WIDTH_IMAGE;
        if (recv_messages < HEIGHT_KERNEL) {
            printf("\n");
            continue;
        }

        do_gaussian((uint8_t *)img, img_gauss, WIDTH_IMAGE, HEIGHT_KERNEL);
        do_sobel0(img_gauss, img_sobel, WIDTH_IMAGE, HEIGHT_KERNEL);
        printf("process ");


        memmove(buffer_target, img_sobel + CENTER_LINE * WIDTH_IMAGE, WIDTH_IMAGE);
        crc = hf_crc32((int8_t *)buffer_target, sizeof(buffer_target) - SIZE_CRC);
        memcpy(buffer_target + WIDTH_IMAGE, &crc, SIZE_CRC);
        printf("copy ");

        // send data to target
        val = hf_send(CPU_TARGET, PORT_TARGET, buffer_target,SIZE_COMM_BUFFER, shift_target);
        if (val)
            printf("hf_send(): error %d\n", val);
        printf("%s send ", hf_selfname());

        printf("\n");

    }

    free(img);
    free(img_gauss);
    free(img_sobel);

    printf("\n\nend of processing!\n\n");
    while (1);
}



void target(void)
{
    uint8_t * filter_image;
    uint8_t * ptr;
    uint16_t cpu, port, size;
    uint32_t crc, received_messages=0;
    uint32_t i, j, k = 0;
    uint32_t time;

    int8_t buf[SIZE_COMM_BUFFER];
    int16_t val;
    int32_t shift;

    printf("cpu %d, name %s, thread %d.\n", hf_cpuid(), hf_selfname(), hf_selfid());

    if (hf_comm_create(hf_selfid(), PORT_TARGET, 0))
        panic(0xff);

    filter_image = (uint8_t *) malloc(height * width);

    time = _readcounter();

    while (1){
        if (received_messages > WIDTH_IMAGE - 5) break;

        shift = hf_recvprobe();
        if (shift < 0) continue;

        val = hf_recv(&cpu, &port, buf, &size, shift);

        if (val) {
            printf("hf_recv(): error %d\n", val);
            continue;
        }
        if ((cpu == CPU_SOURCE) || (cpu == CPU_TARGET)) {
            printf("hf_recv(): error cpu %d\n", cpu);
            continue;
        }

        memcpy(&crc, buf + size - SIZE_CRC, SIZE_CRC);
        printf("T cpu %d, port %d, channel %d, size %d, crc %08x [free queue: %d] ", cpu, port, shift, size, crc,
               hf_queue_count(pktdrv_queue));

        if (hf_crc32(buf, size - SIZE_CRC) != crc) {
            printf(" (CRC32 fail)\n");
            continue;
        }
        printf(" (CRC32 pass) ");

        ptr = filter_image + shift;
        memmove(ptr, buf, WIDTH_IMAGE);

        received_messages++;

        printf("received_messages %d ", received_messages);
        printf("\n");
    }

    time = _readcounter() - time;
    printf("done in %d clock cycles.\n\n", time);

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

void test(void) {
    int32_t i;
    uint32_t crc;
    int8_t buf[500];
    int16_t val, channel;
    printf("cpu %d, name %s, thread %d.\n", hf_cpuid(), hf_selfname(), hf_selfid());

    val = hf_comm_create(hf_selfid(), 1000 + hf_selfid(), 0);
    if (val) {
        printf("hf_comm_create error %d\n", val);
        panic(0xff);
    }
    while (1);

}

void app_main(void)
{
    switch (hf_cpuid()) {
        case CPU_SOURCE:
            hf_spawn(source, 0, 0, 0, "S", 4096);
        case 1:
            hf_spawn(worker, 0, 0, 0, "W1", 4096);
        case 2:
            hf_spawn(worker, 0, 0, 0, "W2", 4096);
        case 3:
        case CPU_TARGET:
            hf_spawn(target, 0, 0, 0, "T", 4096);
    }


//    uint8_t i;
//    if (hf_cpuid() == 0)
//        hf_spawn(source, 0, 0, 0, "s", 2048);
//
//    for (i = 1; i < hf_ncores(); i++, msg++) {
//        if (hf_cpuid() == i) {
//            hf_spawn(worker, 0, 0, 0, "w", 2048);
//        }
//    }
}
