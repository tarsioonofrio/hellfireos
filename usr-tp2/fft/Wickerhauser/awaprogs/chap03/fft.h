/*
 * Declare the functions used in factored discrete Fourier transforms.
 * 
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef FFT_HDR_ALREADY_INCLUDED
#define FFT_HDR_ALREADY_INCLUDED

#include "real.h"
#include "complex.h"

extern void
  fftproduct(                  /* Apply sparse matrix product. */
	     complex *f,       /* Input and output vector.     */
	     int q,            /* Length of `f[]' is N=1<<q.   */
	     const complex *W); /* Exponentials: `Omega(N/2)' */

extern complex * 
  fftomega(			/* Return exp(-PI*i*n/M), */
	   int M);		/* for n=0,1,2,...,|M|-1. */

extern void
  fftnormal(			/* Multiply `f[n].Re' and     */
	    complex *f,		/* `f[n].Im' by `1.0/sqrt(N), */
	    int N );		/* for n=0,1,2,...,N.         */

extern complex *
  dft(				/* Allocate, assign and return  */
      const complex *f,		/* a complex vector, the (1<<q) */
      int q);			/* point DFT of the input `f[]' */


#endif /* FFT_HDR_ALREADY_INCLUDED */
