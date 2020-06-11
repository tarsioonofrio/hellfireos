#include <hellfire.h>

#define Q_SIZE	40

/*
XTEA encryption algorithm

based on reference code released into the public domain by David Wheeler and Roger Needham
the code takes 64 bits of data in v[0] and v[1] and 128 bits of key in key[0] - key[3]

recommended number of rounds is 32 (2 Feistel-network rounds are performed on each iteration).
*/

void xtea_encrypt(uint32_t v[2], const uint32_t key[4], uint32_t num_rounds)
{
    uint32_t i;
    uint32_t v0 = v[0], v1 = v[1], sum = 0, delta = 0x9E3779B9;

    for (i = 0; i < num_rounds; i++){
        v0 += (((v1 << 4) ^ (v1 >> 5)) + v1) ^ (sum + key[sum & 3]);
        sum += delta;
        v1 += (((v0 << 4) ^ (v0 >> 5)) + v0) ^ (sum + key[(sum >> 11) & 3]);
    }
    v[0] = v0; v[1] = v1;
}

void xtea_decrypt(uint32_t v[2], const uint32_t key[4], uint32_t num_rounds)
{
    uint32_t i;
    uint32_t v0 = v[0], v1 = v[1], delta = 0x9E3779B9, sum = delta * num_rounds;

    for (i = 0; i < num_rounds; i++){
        v1 -= (((v0 << 4) ^ (v0 >> 5)) + v0) ^ (sum + key[(sum >> 11) & 3]);
        sum -= delta;
        v0 -= (((v1 << 4) ^ (v1 >> 5)) + v1) ^ (sum + key[sum & 3]);
    }
    v[0] = v0; v[1] = v1;
}


/* XTEA stream cipher, CBC mode
 * CBC mode based on https://en.wikipedia.org/wiki/Block_cipher_mode_of_operation
 */
#define BLOCKLEN	8		// in bytes

void xtea_cbc_encrypt(uint8_t *out, uint8_t *in, uint32_t len, const uint32_t key[4], const uint32_t iv[2])
{
    uint32_t i, rem, block[2], tiv[2];

    rem = len % BLOCKLEN;
    tiv[0] = iv[0];
    tiv[1] = iv[1];
    for (i = 0; i < len; i += BLOCKLEN) {
        memcpy((char *)block, in, BLOCKLEN);
        block[0] ^= tiv[0];
        block[1] ^= tiv[1];
        xtea_encrypt(block, key, 32);
        tiv[0] = block[0];
        tiv[1] = block[1];
        memcpy(out, (char *)block, BLOCKLEN);
        in += BLOCKLEN;
        out += BLOCKLEN;
    }
    if (rem) {
        memcpy((char *)block, in, BLOCKLEN - rem);
        memset((char *)block + rem, 0, BLOCKLEN - rem);
        block[0] ^= tiv[0];
        block[1] ^= tiv[1];
        xtea_encrypt(block, key, 32);
        memcpy(out, (char *)block, BLOCKLEN - rem);
    }
}

void xtea_cbc_decrypt(uint8_t *out, uint8_t *in, uint32_t len, const uint32_t key[4], const uint32_t iv[2])
{
    uint32_t i, rem, block[2], block2[2], tiv[2];

    rem = len % BLOCKLEN;
    tiv[0] = iv[0];
    tiv[1] = iv[1];
    for (i = 0; i < len; i += BLOCKLEN) {
        memcpy((char *)block, in, BLOCKLEN);
        block2[0] = block[0];
        block2[1] = block[1];
        xtea_decrypt(block, key, 32);
        block[0] ^= tiv[0];
        block[1] ^= tiv[1];
        tiv[0] = block2[0];
        tiv[1] = block2[1];
        memcpy(out, (char *)block, BLOCKLEN);
        in += BLOCKLEN;
        out += BLOCKLEN;
    }
    if (rem) {
        memcpy((char *)block, in, BLOCKLEN - rem);
        memset((char *)block + rem, 0, BLOCKLEN - rem);
        tiv[0] = block[0];
        tiv[1] = block[1];
        xtea_decrypt(block, key, 32);
        block[0] ^= tiv[0];
        block[1] ^= tiv[1];
        memcpy(out, (char *)block, BLOCKLEN - rem);
    }
}


struct queue *q;
mutex_t m;
	
void sender(void)
{
	int32_t i = 0;
	int8_t *buf;

    uint32_t xtea_key[4] = {0xf0e1d2c3, 0xb4a59687, 0x78695a4b, 0x3c2d1e0f};
    uint32_t iv[2] = {0x11223344, 0x55667788};


    for(;;){
		if (hf_queue_count(q) < Q_SIZE){
			buf = malloc(sizeof(int8_t) * 100);
			if (buf){
				hf_mtxlock(&m);
                sprintf(buf, "hello from task %s, counting %d (item at %08x, %d free)", hf_selfname(), i++, (uint32_t)buf, hf_freemem());
                xtea_cbc_encrypt(buf, buf, sizeof(buf), xtea_key, iv);
				if (hf_queue_addtail(q, buf)){
//					printf("queue is full!\n");
					free(buf);
				}
				hf_mtxunlock(&m);
			}else{
//				printf("malloc() failed!\n");
			}
		}
	}
}

void receiver(void)
{
	int8_t *buf;

    uint32_t xtea_key[4] = {0xf0e1d2c3, 0xb4a59687, 0x78695a4b, 0x3c2d1e0f};
    uint32_t iv[2] = {0x11223344, 0x55667788};

    for(;;){
		if (hf_queue_count(q)){
			hf_mtxlock(&m);
			buf = hf_queue_remhead(q);
            xtea_cbc_decrypt(buf, buf, sizeof(buf), xtea_key, iv);
			hf_mtxunlock(&m);
			if (buf){
				printf("task %s -> %s\n", hf_selfname(), buf);
				free(buf);
			}
		}
	}
}

void logthread(void)
{
	for(;;){
//		printf("queue: %d (tick time %dus)\n", hf_queue_count(q), hf_ticktime());
		hf_yield();
	}
}

void app_main(void){
	hf_mtxinit(&m);
	q = hf_queue_create(Q_SIZE);

	hf_spawn(sender, 0, 0, 0, "sender 1", 1024);
	hf_spawn(sender, 0, 0, 0, "sender 2", 1024);
	hf_spawn(sender, 0, 0, 0, "sender 3", 1024);
	hf_spawn(sender, 0, 0, 0, "sender 4", 1024);
	hf_spawn(sender, 0, 0, 0, "sender 5", 1024);
	hf_spawn(sender, 0, 0, 0, "sender 6", 1024);
	hf_spawn(receiver, 0, 0, 0, "receiver 1", 1024);
	hf_spawn(receiver, 0, 0, 0, "receiver 2", 1024);
	hf_spawn(receiver, 0, 0, 0, "receiver 3", 1024);
	hf_spawn(logthread, 100, 1, 100, "log", 1024);
}

