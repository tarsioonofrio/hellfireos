//
// Created by tarsio on 09/06/2020.
//
// https://rosettacode.org/wiki/Fast_Fourier_transform#C
//C
//
//        Inplace FFT with O(n) memory usage. Note: array size is assumed to be power of 2 and not checked by code; you can just pad it with 0 otherwise.
//Also, complex is C99 standard.



#include <stdio.h>
#include <math.h>
#include <complex.h>

double PI;
typedef double complex cplx;

void _fft(cplx buf[], cplx out[], int n, int step)
{
    if (step < n) {
        _fft(out, buf, n, step * 2);
        _fft(out + step, buf + step, n, step * 2);

        for (int i = 0; i < n; i += 2 * step) {
            cplx t = cexp(-I * PI * i / n) * out[i + step];
            buf[i / 2]     = out[i] + t;
            buf[(i + n)/2] = out[i] - t;
        }
    }
}

void fft(cplx buf[], int n)
{
    cplx out[n];
    for (int i = 0; i < n; i++) out[i] = buf[i];

    _fft(buf, out, n, 1);
}


void show(const char * s, cplx buf[]) {
    printf("%s", s);
    for (int i = 0; i < 8; i++)
        if (!cimag(buf[i]))
            printf("%g ", creal(buf[i]));
        else
            printf("(%g, %g) ", creal(buf[i]), cimag(buf[i]));
}

int main()
{
    PI = atan2(1, 1) * 4;
    cplx buf[] = {1, 1, 1, 1, 0, 0, 0, 0};

    show("Data: ", buf);
    fft(buf, 8);
    show("\nFFT : ", buf);

    return 0;
}



//Output:
//
//Data: 1 1 1 1 0 0 0 0
//FFT : 4 (1, -2.41421) 0 (1, -0.414214) 0 (1, 0.414214) 0 (1, 2.41421)