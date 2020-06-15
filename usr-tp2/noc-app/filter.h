//
// Created by tarsio on 14/06/2020.
//

uint32_t isqrt(uint32_t a);
uint8_t gaussian(uint8_t buffer[5][5]);
uint8_t sobel(uint8_t buffer[3][3]);
void do_gaussian(uint8_t *input, uint8_t *output, int32_t width, int32_t height);
void do_sobel(uint8_t *input, uint8_t *output, int32_t width, int32_t height);