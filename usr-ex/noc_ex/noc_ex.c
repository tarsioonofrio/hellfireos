#include <hellfire.h>
#include <noc.h>

typedef struct{
    uint8_t size;
    uint8_t *task;
} struct_node;

typedef struct{
    uint8_t node; //target_node
    uint8_t task;
    uint8_t size;
    int8_t *buf;
} struct_task;


void sender(void) {
    int32_t i;
    uint8_t n;
    uint32_t crc;
    int8_t buf[500];
    int16_t val, channel;
    int task_num = 3;
    struct_task task_send[3];
    task_send[0].node = 1;
    task_send[0].task = 4;
    task_send[0].size = 75;
    task_send[0].buf = (int8_t *)malloc(task_send[0].size*sizeof(int8_t));
    task_send[1].node = 1;
    task_send[1].task = 5;
    task_send[1].size = 55;
    task_send[1].buf = (int8_t *)malloc(task_send[1].size*sizeof(int8_t));
    task_send[2].node = 1;
    task_send[2].task = 3;
    task_send[2].size = 85;
    task_send[2].buf = (int8_t *)malloc(task_send[2].size*sizeof(int8_t));

    if (hf_comm_create(hf_selfid(), 1000, 0)){
        panic(0xff);
    }

    delay_ms(50);

    srand(hf_cpuid());

    // generate a unique channel number for this CPU
    while (1) {
        for (n = 0; n < task_num; n++) {
            for (i = 0; i < task_send[n].size - 4; i++) {
                task_send[n].buf[i] = random() % 255;
            }
            crc = hf_crc32(task_send[n].buf, task_send[n].size - 4);
            memcpy(task_send[n].buf + task_send[n].size - 4, &crc, 4);
            val = hf_send(task_send[n].node, 5000 + task_send[n].task, task_send[n].buf, task_send[n].size, 0); //task_send[n].task);
            if (val)
                printf("hf_send(): error %d\n", val);
//            delay_ms(10);
        }
    }
}

void receiver(void)
{
	int8_t buf[1500];
	uint16_t send_node, send_port, send_size, send_task;
	int16_t val;
	uint32_t crc;
	int32_t i;
    uint8_t task_id = strtol(hf_selfname(), (char **)NULL, 10);

	printf("hf_selfid %d, hf_selfname %d\n", hf_selfid(), task_id);

	if (hf_comm_create(hf_selfid(), 5000 + task_id, 0))
		panic(0xff);
	
	while (1){
		i = hf_recvprobe();
		if (i >= 0) {
			val = hf_recv(&send_node, &send_port, buf, &send_size, 0);
			if (val){
				printf("hf_recv(): error %d\n", val);
			} else {
				memcpy(&crc, buf + send_size - 4, 4);
				printf("cpu %d, port %d, channel %d, size %d, crc %08x [free queue: %d], task_id %d", send_node, send_port, i, send_size, crc, hf_queue_count(pktdrv_queue), task_id);
				if (hf_crc32(buf, send_size - 4) == crc)
					printf(" (CRC32 pass)\n");
				else
					printf(" (CRC32 fail)\n");
			}
		}
	}
}

void app_main(void)
{
    switch (hf_cpuid()) {
        case 0:
            hf_spawn(sender, 0, 0, 0, "1", 4096);
        case 1:
            hf_spawn(receiver, 0, 0, 0, "4", 4096);
            hf_spawn(receiver, 0, 0, 0, "5", 4096);
            hf_spawn(receiver, 0, 0, 0, "3", 4096);
//        default:
//            hf_spawn(receiver, 0, 0, 0, "receiver", 4096);
    }
}
