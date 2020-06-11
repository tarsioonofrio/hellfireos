/*
 * Declare functions used to perform the discrete Hartley transform.
 * 
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef HARTLEY_HDR_ALREADY_INCLUDED
#define HARTLEY_HDR_ALREADY_INCLUDED

#include "real.h"

extern void
  dhtproduct(			/* Apply sparse matrix product. */
	     real *f,		/* Input and output vector.     */
	     int q,		/* Length of `f[]' is N=1<<q.   */
	     const real *c,	/* Cosines c_0,...,c_{N/2-1}.   */
	     const real *s);	/* Sines s_0,s_1,...,s_{N/2-1}. */

extern void 
  dhtcossin(			/* Fill the arrays with:  */
	    real *c,		/*   (real)cos(n*PI/N),   */
	    real *s,		/*   (real)sin(n*PI/N),   */
	    int N);		/* for n=0,1,2,...,N/2-1. */

extern void
  dhtnormal(
	    real *f,		/* Multiply `f[n]' for n=0,1, */
	    int N );		/* ...,N  by `1.0/sqrt(N).'   */

extern real *
  dht(				/* Allocate, assign and return  */
      const real *f,		/* a real vector, the (1<<q)   */
      int q);			/* point DFT of the input `f[]' */

#endif /* HARTLEY_HDR_ALREADY_INCLUDED */
