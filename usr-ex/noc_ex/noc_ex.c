#include <hellfire.h>
#include <noc.h>

#define EDGE_SOURCE 0
#define EDGE_TARGET 1
#define EDGE_WEIGHT 2
#define EDGE_SIZE 17


int arr_task_edge [EDGE_SIZE][3] = {
        {1, 4, 75},
        {1, 5, 55},
        {1, 3, 85},
        {4, 2, 35},
        {4, 8, 88},
        {5, 8, 780},
        {3, 8, 64},
        {3, 7, 224},
        {2, 6, 120},
        {2, 7, 450},
        {8, 9, 1150},
        {8, 10, 116},
        {6, 9, 155},
        {7, 11, 890},
        {7, 10, 435},
        {9, 11, 220},
        {10, 11, 440}
};

uint8_t arr_task_node [11] = {0, 2, 1, 1, 1, 3, 3, 2, 4, 4, 5};


typedef struct{
    uint8_t size;
    uint8_t *task;
} struct_node;

typedef struct{
    uint8_t target_node; //target_node
    uint8_t target_task;
    uint8_t size;
    int8_t *buf;
} struct_task;


void sender(void) {
    int32_t i;
    uint8_t n;
    uint32_t crc;
//    int8_t buf[500];
    int16_t val;
    int task_num=0;

    uint8_t task_id = strtol(hf_selfname(), (char **)NULL, 10);
    printf("task_id %d, arr_task %d, hf_cpuid %d\n", task_id, arr_task_node[task_id - 1], hf_cpuid());

    for (n = 0; n < EDGE_SIZE; n++) {
        if (task_id == arr_task_edge[n][EDGE_SOURCE])
            task_num++;
    }

    struct_task task_send[task_num];

    for (n = 0; n < EDGE_SIZE; n++) {
        if (task_id == arr_task_edge[n][EDGE_SOURCE]){
            printf("%d %d %d %d\n", arr_task_edge[n][EDGE_SOURCE], arr_task_edge[n][EDGE_TARGET],
                    arr_task_edge[n][EDGE_WEIGHT], arr_task_node[arr_task_edge[n][EDGE_TARGET] -1]);

            task_send[n].target_task = arr_task_edge[n][EDGE_TARGET];
            task_send[n].target_node = arr_task_node[arr_task_edge[n][EDGE_TARGET] -1];
            task_send[n].size = arr_task_edge[n][EDGE_WEIGHT];
            task_send[n].buf = (int8_t *)malloc(task_send[n].size*sizeof(int8_t));
        }
    }

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
            val = hf_send(task_send[n].target_node, 5000 + task_send[n].target_task, task_send[n].buf, task_send[n].size, 0); // task_send[n].target_task);
            if (val)
                printf("hf_send(): error %d\n", val);

//            delay_ms(10);
        }
    }
}

_Noreturn void receiver(void)
{
	int8_t buf[1500];
	uint16_t source_cpu, source_port, source_size, source_channel;
	int16_t val;
	uint32_t crc;
	int32_t i;

    uint8_t task_id = strtol(hf_selfname(), (char **)NULL, 10);
    printf("task_id %d, arr_task %d, hf_cpuid %d\n", task_id, arr_task_node[task_id - 1], hf_cpuid());
    source_channel = arr_task_node[task_id - 1];

	if (hf_comm_create(hf_selfid(), 5000 + task_id, 0))
		panic(0xff);
	
	while (1){
		i = hf_recvprobe();
		if (i >= 0) {
			val = hf_recv(&source_cpu, &source_port, buf, &source_size, 0);//source_channel);
			if (val){
				printf("hf_recv(): error %d\n", val);
			} else {
				memcpy(&crc, buf + source_size - 4, 4);
				printf("cpu %d, port %d, channel %d, size %d, crc %08x [free queue: %d]", source_cpu, source_port, i, source_size, crc, hf_queue_count(pktdrv_queue));
				if (hf_crc32(buf, source_size - 4) == crc)
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
        default:
            printf("\n *** hf_cpuid %d *** \n", hf_cpuid());
    }
}
