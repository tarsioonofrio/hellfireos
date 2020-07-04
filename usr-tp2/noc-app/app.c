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
    uint32_t crc;
    int8_t buf[SIZE_COMM_BUFFER];
    uint8_t * ptr;
    int16_t val;
//    int16_t control[NUM_CPU][2];
    int16_t control=0;

    uint16_t cpu, port, size;
    int32_t channel;

    ptr = image;

    if (hf_comm_create(hf_selfid(), PORT_SOURCE, 0))
        panic(0xff);

    delay_ms(50);

    srand(hf_cpuid());
    printf("Start\n");

    // generate a unique channel number for this CPU
//    channel = hf_cpuid();
    while (1){
        if (control > WIDTH_IMAGE - 1) break;

        channel = hf_recvprobe();
        if (channel < 0) continue;

        val = hf_recv(&cpu, &port, buf, &size, channel);
        if (val){
            printf("hf_send(): error %d\n", val);
            continue;
        }
        //        val = hf_sendack(1, 2000, (int8_t *) ptr, SIZE_BUFFER, 1, TIMEOUT);

        crc = hf_crc32((int8_t *)ptr, sizeof(buf) - SIZE_CRC);
        memcpy(ptr + WIDTH_IMAGE, &crc, SIZE_CRC);
        val = hf_send(cpu, PORT_WORKER, (int8_t *) ptr, SIZE_COMM_BUFFER, 1);
        if (val){
            printf("hf_send(): error %d\n", val);
            continue;
        }//        val = hf_sendack(1, 2000, (int8_t *) ptr, SIZE_BUFFER, 1, TIMEOUT);
//        if (val) printf("sender, hf_sendack(): error %d\n", val);
        control++;
        ptr = ptr + WIDTH_IMAGE;
        printf("ptr %d\n", ptr);
        printf("control %d\n", control);
//        delay_ms(20);
    }
}



void worker(void)
{
    int8_t buf[SIZE_COMM_BUFFER], buf_dummy[1];
    uint16_t cpu, port, size, cpuid, recv_messages=0;
    int16_t val, x, z = 0;
    uint32_t crc;
    int32_t channel;
    int16_t control=0;

    uint8_t *img_gauss, *img_sobel, *img;
//    uint32_t time;
    cpuid = hf_cpuid();

    img_sobel = (uint8_t *) malloc(SIZE_PROC_BUFFER);
    img = (uint8_t *) malloc(SIZE_PROC_BUFFER);
    img_gauss = (uint8_t *) malloc(SIZE_PROC_BUFFER);

    if (img_gauss == NULL || img_sobel == NULL){
        printf("\nmalloc() failed!\n");
        for(;;);
    }

    if (hf_comm_create(hf_selfid(), PORT_WORKER, 0))
        panic(0xff);

    while (1){
        if (control > WIDTH_IMAGE - 1) break;

        // request data to source
        val = hf_send(0, PORT_SOURCE, buf_dummy, sizeof(buf_dummy), cpuid);
        if (val) {
            printf("hf_recv(): error %d\n", val);
            continue;
        }

        // receive data from source
        channel = hf_recvprobe();
        if (channel < 0) continue;

        val = hf_recv(&cpu, &port, buf, &size, channel);
        if (val){
//            printf("hf_recv(): error %d\n", val);
            continue;
        }

        memcpy(&crc, buf + size - SIZE_CRC, SIZE_CRC);
        printf("W cpu %d, port %d, channel %d, size %d, crc %08x [free queue: %d]", cpu, port, channel, size, crc,
               hf_queue_count(pktdrv_queue));

        if (hf_crc32(buf, size - SIZE_CRC) != crc) {
//            printf(" (CRC32 fail)\n");
            continue;
        }
//        printf(" (CRC32 pass)\n");

        memmove(img, img + WIDTH_IMAGE, SIZE_PROC_BUFFER - WIDTH_IMAGE);
        memmove(img + (SIZE_PROC_BUFFER - WIDTH_IMAGE), buf, WIDTH_IMAGE);

        recv_messages++;
//        printf("recv_messages %d\n", recv_messages);

        if (recv_messages < HEIGHT_KERNEL) continue;

        do_sobel0((uint8_t *)img, img_sobel, WIDTH_IMAGE, HEIGHT_KERNEL);
        do_gaussian(img_sobel, img_gauss, WIDTH_IMAGE, HEIGHT_KERNEL);

        for(x = WIDTH_IMAGE * CENTER_LINE; x < WIDTH_IMAGE * (CENTER_LINE + 1); x++){
            printf("0x%02x", img_sobel[x]);
            printf(", ");
            if ((++z % 16) == 0) printf("\n");
        }

//        printf("control %d\n", control);
        control++;

        // send data to target
        val = hf_send(1, PORT_TARGET, img_gauss + CENTER_LINE * SIZE_PROC_BUFFER,
                      SIZE_COMM_BUFFER, 1);
        if (val)
            printf("hf_send(): error %d\n", val);
    }

//    for(x = WIDTH_IMAGE * CENTER_LINE; x < SIZE_PROC_BUFFER - 1; x++){
//        printf("0x%02x", img_gauss[x]);
//        printf(", ");
//        if ((++z % 16) == 0) printf("\n");
//    }

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
    int16_t control=0;

    if (hf_comm_create(hf_selfid(), PORT_TARGET, 0))
        panic(0xff);

    filter_image = (uint8_t *) malloc(height * width);
    ptr = filter_image;
    ptr = ptr + WIDTH_IMAGE * CENTER_LINE;

    while (1){
        if (control > WIDTH_IMAGE - 1) break;

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

        printf("control %d\n", control);
        control++;
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


void _source(void)
{
    int8_t buf[SIZE_COMM_BUFFER], buf_dummy[1], ready=0;
    uint16_t cpu, port, size, cpuid;
    int16_t val;
    uint32_t crc;
    int32_t channel;
    int16_t control=0;

    if (hf_comm_create(hf_selfid(), PORT_TARGET, 0))
        panic(0xff);

    cpuid = hf_cpuid();

    while (1){
        if (control > WIDTH_IMAGE - 1)
            break;

        // receive ready from worker
        channel = hf_recvprobe();
        if (channel == CHANNEL_WORKER) {
            ready = 1;
            hf_recv(&cpu, &port, buf, &size, channel);
        }


        val = hf_send(0, PORT_SOURCE, buf_dummy, sizeof(buf_dummy), cpuid);
        if (val) {
            printf("hf_recv(): error %d\n", val);
            continue;
        }

        channel = hf_recvprobe();
        if (channel < 0) continue;

        val = hf_recv(&cpu, &port, buf, &size, channel);
//        val = hf_recvack(&cpu, &port, buf, &size, i);
//        if (val) printf("Master, hf_sendack(): error %d\n", val);

        if (val) {
            printf("hf_recv(): error %d\n", val);
            continue;
        }

        memcpy(&crc, buf + size - SIZE_CRC, SIZE_CRC);
        printf("R cpu %d, port %d, channel %d, size %d, crc %08x [free queue: %d]", cpu, port, channel, size, crc,
               hf_queue_count(pktdrv_queue));

        if (hf_crc32(buf, size - SIZE_CRC) != crc) {
            printf(" (CRC32 fail)\n");
            continue;
        }

        printf(" (CRC32 pass)\n");

        val = hf_send(2, PORT_WORKER, buf, size, 1);
        if (val) {
            printf("hf_send(): error %d\n", val);
        }
        printf("control %d\n", control);
        control++;
        delay_ms(5);

    }
}

void _target(void)
{
    int8_t buf[SIZE_COMM_BUFFER];
    uint16_t cpu, port, size, cpuid;
    int16_t val;
    uint32_t crc;
    int32_t i;
    int16_t control=0;

    if (hf_comm_create(hf_selfid(), PORT_SOURCE, 0))
        panic(0xff);

    cpuid = hf_cpuid();

    while (1){
        if (control > WIDTH_IMAGE - 1)
            break;

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
        printf("R cpu %d, port %d, channel %d, size %d, crc %08x [free queue: %d]", cpu, port, i, size, crc,
               hf_queue_count(pktdrv_queue));

        if (hf_crc32(buf, size - SIZE_CRC) != crc) {
            printf(" (CRC32 fail)\n");
            continue;
        }

        printf(" (CRC32 pass)\n");

        val = hf_send(4, PORT_TARGET, buf, size, 1);
        if (val) {
            printf("hf_send(): error %d\n", val);
        }
        printf("control %d\n", control);
        control++;
//        delay_ms(10);

    }
}



void test(void) {
    int32_t comm;
    comm = hf_comm_create(hf_selfid(), 1000, 0);
    if (comm) {
        printf("ERROR %d", comm);
        panic(0xff);
    }
}


void sender0(void)
{
    int32_t msg = 0;
    int8_t buf[500], i;
    int16_t val, channel;

    if (hf_comm_create(hf_selfid(), 1000, 0))
        panic(0xff);

    delay_ms(50);

    // generate a unique channel number for this CPU
    channel = hf_cpuid();
    while (1){
        for (i = 0; i < hf_ncores(); i++, msg++) {
            if (hf_cpuid() != i) {
                sprintf(buf, "i am cpu %d, thread %d, channel %d: msg %d size: %d\n", hf_cpuid(), hf_selfid(), channel,
                        msg++, sizeof(buf));
                val = hf_send(i, 5000, buf, sizeof(buf), channel);
                if (val)
                    printf("hf_send(): error %d\n", val);
                delay_ms(10);
            }
        }
    }
}


void receiver0(void)
{
    int8_t buf[1500];
    uint16_t cpu, task, size;
    int16_t val;
    int32_t i;

    if (hf_comm_create(hf_selfid(), 5000, 0))
        panic(0xff);

    while (1){
        i = hf_recvprobe();
        if (i >= 0) {
            val = hf_recv(&cpu, &task, buf, &size, i);
            sprintf(buf, "i am cpu %d, task %d, channel %d: size: %d\n", cpu, task, i, sizeof(buf));
            if (val)
                printf("hf_recv(): error %d\n", val);
            else
                printf("%s", buf);
        }
    }
}


void app_main(void)
{
    switch (hf_cpuid()) {
        case 0:
            hf_spawn(source, 0, 0, 0, "s", 2048);
//            hf_spawn(test, 0, 0, 0, "m", 2048);
//            hf_spawn(master_receiver, 0, 0, 0, "R", 4096);
        case 1:
            hf_spawn(worker, 0, 0, 0, "w", 2048);
//        case 2:
//            hf_spawn(test, 0, 0, 0, "w", 2048);
//        case 2:
//            hf_spawn(test, 0, 0, 0, "t", 2048);
//        case 4:
//            hf_spawn(master_receiver, 0, 0, 0, "MR", 4096);
    }
}
