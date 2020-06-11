/*
 * Discrete trigonometric transforms based on FFT.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "real.h"
#include "common.h"		/* For the constants SR2 and PI.   */
#include "bitrev.h"		/* Generic bit-reversal functions. */
#include "fft.h"
#include "dtts.h"

/****************************************************************
 * Fill the preallocated arrays of real's `c[]' and `s[]', each
 * assumed to be of length `N/2', with the respective  values:
 * 	c[n] = cos(n*PI/N), n=0,1,...,N/2-1.
 * 	s[n] = sin(n*PI/N), n=0,1,...,N/2-1.
 */ 
extern void 
  mkcossin(			/* Fill the arrays with:  */
	   real *c,		/*   (real)cos(n*PI/N),   */
	   real *s,		/*   (real)sin(n*PI/N),   */
	   int N)		/* for n=0,1,2,...,N/2-1. */
{
  double factor, theta;
  int n;

  factor = PI/(double)N;

  c[0] = (real)1.0;
  s[0] = (real)0.0;

  theta = 0.0;
  for(n=1; n<(N>>1); n++)
    {
      theta += factor;
      c[n] = (real)cos(theta);
      s[n] = (real)sin(theta);
    }
  return;
}


/****************************************************************
 * dct_i()
 *
 * This function accepts an input array of N + 1 `real' elements,
 * for N=2^q, and replaces its elements with the DCT-I transform.
 * It uses the radix-2 FFT on 2^{q+1} points, using a temporary
 * array of 2^{q+1} extra `complex' data structs which it
 * allocates and frees upon termination.  It is its own inverse.
 */
extern void
  dct_i(
	real *x,		/* Input; length (1<<q)+1. */
	int q)			/* Nonnegative integer. */
{
  complex *f, *W;
  int N, NN, n, m;
  real norm;

  N = (1<<q);			/* number of inputs == (1+N) */
  NN = (N<<1);			/* length of temp array ==2N */

  /* Allocate and zero the temporary array */
  f = (complex *)calloc(NN, sizeof(complex));

  /* Inject the input vector into `f[]' using the U matrix: */
  f[0].Re = x[0]*SR2;
  f[N].Re = x[N]*SR2;
  for(n=1, m=NN-1; n<N; n++, m--)
    {
      f[m].Re = f[n].Re = x[n];
    }

  /* Do a 2N-point non-normalized FFT in place on `f[]' */
  bitrevi(f, q+1, sizeof(complex));
  W = fftomega(N);		/* generate sines/cosines */
  fftproduct(f, q+1, W);	/* apply sparse matrices */
  free(W);			/* clean up. */
  
  /* Project back to the array `x[]' using U*: */
  x[0] = f[0].Re * SR2;
  x[N] = f[N].Re * SR2;
  for(n=1, m=NN-1; n<N; n++, m--)
    {
      x[n] = f[m].Re + f[n].Re;
    }
  free(f);			/* clean up the temp array */

  /* Normalize: divide by 2*sqrt(2N): */
  if(q&1)
    norm = 0.5/(real)(1<<(q+1)/2);
  else
    norm = 0.5/(real)sqrt((double)NN);
  for(n=0; n<=N; n++)  x[n] *= norm;

  return;
}


/****************************************************************
 * dst_i()
 *
 * This function accepts an input array of N - 1 `real' elements,
 * for N=2^q, and replaces its elements with the DST-I transform.
 * It uses the radix-2 FFT on 2^{q+1} points, using a temporary
 * array of 2^{q+1} extra `complex' data structs which it
 * allocates and frees upon termination.  It is its own inverse.
 * Note: we assume that the input array is x[0],x[1],...,x[N-1],
 * with the input values located in x[1],...,x[N-1].
 */
extern void
  dst_i(
	real *x,	/* Input; length (1<<q). */
	int q)		/* Nonnegative integer. */
{
  complex *f, *W;
  int N, NN, n, m;
  real norm;

  N = (1<<q);			/* number of inputs == (N-1) */
  NN = (N<<1);			/* length of temp array ==2N */

  /* Allocate and zero the temporary array */
  f = (complex *)calloc(NN, sizeof(complex));

  /* Inject the input vector into `f[]' using the V matrix: */
  for(n=1, m=NN-1; n<N; n++, m--)
    {
      f[n].Re =  x[n];
      f[m].Re = -x[n];
    }

  /* Do a 2N-point non-normalized FFT in place on `f[]' */
  bitrevi(f, q+1, sizeof(complex));
  W = fftomega(N);		/* generate sines/cosines */
  fftproduct(f, q+1, W);	/* apply sparse matrices */
  free(W);			/* clean up. */
  
  /* Project back to the array `x[]' using iV*: */
  for(n=1, m=NN-1; n<N; n++, m--)
    {
      x[n] = f[m].Im - f[n].Im;
    }
  free(f);			/* clean up the temp array */

  /* Normalize: divide by 2*sqrt(2N): */
  if(q&1)
    norm = 0.5/(real)(1<<(q+1)/2);
  else
    norm = 0.5/(real)sqrt((double)NN);
  for(n=1; n<N; n++)  x[n] *= norm;

  return;
}


/****************************************************************
 * dct_ii()
 *
 * This function accepts an input array of N `real' elements,
 * for N=2^q, and replaces its elements with the DCT-II transform.
 * It uses the radix-2 FFT on 2^{q+1} points, using a temporary
 * array of 2^{q+1} extra `complex' data structs which it
 * allocates and frees upon termination.  Its inverse is the
 * DCT-III function below.
 */
extern void
  dct_ii(
	 real *x,	/* Input; length (1<<q). */
	 int q)		/* Nonnegative integer. */
{
  complex *f, *W;
  int N, NN, n, m;
  real norm, *c, *s;

  N = (1<<q);			/* number of inputs == N */
  NN = (N<<1);			/* length of temp array ==2N */

  /* Allocate and zero the temporary array */
  f = (complex *)calloc(NN, sizeof(complex));

  /* Inject the input vector into `f[]' using the Q matrix: */
  for(n=0, m=NN-1; n<N; n++, m--)
    {
      f[n].Re =  x[n];
      f[m].Re =  x[n];
    }

  /* Do a 2N-point non-normalized FFT in place on `f[]' */
  bitrevi(f, q+1, sizeof(complex));
  W = fftomega(N);		/* generate sines/cosines */
  fftproduct(f, q+1, W);	/* apply sparse matrices */
  free(W);			/* clean up the exponentials. */
  
  /* Project back to the array `x[]' using P*: */
  x[0] = f[0].Re * SR2;
  c = (real *)malloc(NN*sizeof(real));  s = c+N;
  mkcossin(c, s, NN);
  for(n=1, m=NN-1; n<N; n++, m--)
    {
      x[n] =
	CRRMULRE(f[n], c[n], -s[n])
	  + CRRMULRE(f[m], c[n], s[n]);
    }
  free(f);			/* clean up the temp array */
  free(c);			/* clean up the sines/cosines. */

  /* Normalize: divide by 2*sqrt(2N): */
  if(q&1)
    norm = 0.5/(real)(1<<(q+1)/2);
  else
    norm = 0.5/(real)sqrt((double)NN);
  for(n=0; n<N; n++)  x[n] *= norm;

  return;
}

/****************************************************************
 * dst_ii()
 *
 * This function accepts an input array of N `real' elements,
 * for N=2^q, and replaces its elements with the DST-II transform.
 * It uses the radix-2 FFT on 2^{q+1} points, using a temporary
 * array of 2^{q+1} extra `complex' data structs which it
 * allocates and frees upon termination.  Its inverse is the
 * DST-III function below.
 */
extern void
  dst_ii(
	 real *x,	/* Input; length (1<<q). */
	 int q)		/* Nonnegative integer. */
{
  complex *f, *W;
  int N, NN, n, m;
  real norm, *c, *s;

  N = (1<<q);			/* number of inputs == N */
  NN = (N<<1);			/* length of temp array ==2N */

  /* Allocate and zero the temporary array */
  f = (complex *)calloc(NN, sizeof(complex));

  /* Inject the input vector into `f[]' using iQ~: */
  for(n=0, m=NN-1; n<N; n++, m--)
    {
      f[n].Im =  x[n];
      f[m].Im = -x[n];
    }

  /* Do a 2N-point non-normalized FFT in place on `f[]' */
  bitrevi(f, q+1, sizeof(complex));
  W = fftomega(N);		/* generate sines/cosines */
  fftproduct(f, q+1, W);	/* apply sparse matrices */
  free(W);			/* clean up the exponentials. */
  
  /* Project back to the array `x[]' using P~*: */
  c = (real *)malloc(NN*sizeof(real)); s=c+N;
  mkcossin(c, s, NN);
  for(n=1, m=NN-1; n<N; n++, m--)
    {
      x[n-1] = /* x[n-1] = Re{ f[n]*Wbar[n] - f[2N-n]*W[n] } */
	 CRRMULRE(f[n], c[n], -s[n])
	  - CRRMULRE(f[m], c[n], s[n]);
    }
  x[N-1] = f[N].Im * SR2;  /* x[N-1]= Re{ f[N]*(-i)*sqrt(2) } */
  free(f);			/* clean up the temp array */
  free(c);			/* clean up the exponentials. */

  /* Normalize: divide by 2*sqrt(2N): */
  if(q&1)
    norm = 0.5/(real)(1<<(q+1)/2);
  else
    norm = 0.5/(real)sqrt((double)NN);
  for(n=0; n<N; n++)  x[n] *= norm;

  return;
}


/****************************************************************
 * dct_iii()
 *
 * This function accepts an input array of N `real' elements,
 * for N=2^q, and replaces its elements with the DCT-III transform.
 * It uses the radix-2 FFT on 2^{q+1} points, using a temporary
 * array of 2^{q+1} extra `complex' data structs which it
 * allocates and frees upon termination.  Its inverse is the
 * DCT-III function below.
 */
extern void
  dct_iii(
	 real *x,	/* Input; length (1<<q). */
	 int q)		/* Nonnegative integer. */
{
  complex *f, *W;
  int N, NN, n, m;
  real norm, *c, *s;

  N = (1<<q);			/* number of inputs == N */
  NN = (N<<1);			/* length of temp array ==2N */

  /* Allocate and zero the temporary array */
  f = (complex *)calloc(NN, sizeof(complex));

  /* Inject the input vector into `f[]' using the Pbar matrix: */
  c = (real *)malloc(NN*sizeof(real)); s = c+N;
  mkcossin(c, s, NN);
  f[0].Re = x[0] * SR2;		/* f[0] = x[0]*sqrt(2) */
  for(n=1, m=NN-1; n<N; n++, m--)
    {
      f[m].Re = x[n] * c[n];	/* f[2N-n] = x[n]*W[n] */
      f[m].Im = x[n] * s[n];
      f[n].Re =  f[m].Re;	/* f[n] = x[n]*Wbar[n] */
      f[n].Im = -f[m].Im;
    }
  free(c);			/* clean up the sines/cosines. */

  /* Do a 2N-point non-normalized FFT in place on `f[]' */
  bitrevi(f, q+1, sizeof(complex));
  W = fftomega(N);		/* generate sines/cosines */
  fftproduct(f, q+1, W);	/* apply sparse matrices */
  free(W);			/* clean up the exponentials. */
  
  /* Project back to the array `x[]' using Q*: */
  for(n=0, m=NN-1; n<N; n++, m--)
    {
      x[n] = f[m].Re + f[n].Re;
    }
  free(f);			/* clean up the temp array */

  /* Normalize: divide by 2*sqrt(2N): */
  if(q&1)
    norm = 0.5/(real)(1<<(q+1)/2);
  else
    norm = 0.5/(real)sqrt((double)NN);
  for(n=0; n<N; n++)  x[n] *= norm;

  return;
}

/****************************************************************
 * dst_iii()
 *
 * This function accepts an input array of N `real' elements, for
 * N=2^q, and replaces its elements with the DST-III transform.
 * It uses the radix-2 FFT on 2^{q+1} points, using a temporary
 * array of 2^{q+1} extra `complex' data structs which it
 * allocates and frees upon termination.  Its inverse is the
 * DCT-III function below.
 */
extern void
  dst_iii(
	  real *x,	/* Input; length (1<<q). */
	  int q)	/* Nonnegative integer. */
{
  complex *f, *W;
  int N, NN, n, m;
  real norm, *c, *s;

  N = (1<<q);			/* number of inputs == N */
  NN = (N<<1);			/* length of temp array ==2N */

  /* Allocate and zero the temporary array */
  f = (complex *)calloc(NN, sizeof(complex));

  /* Inject the input into `f[]' using the P~bar matrix: */
  c = (real *)malloc(NN*sizeof(real)); s = c+N;
  mkcossin(c, s, NN);
  for(n=1, m=NN-1; n<N; n++, m--)
    {
      f[n].Re =  x[n-1] * c[n];	/* f[n] = x[n-1]*Wbar[n] */
      f[n].Im = -x[n-1] * s[n];
      f[m].Re = -f[n].Re;	/* f[2N-n] = -x[n-1]*W[n] */
      f[m].Im =  f[n].Im;
    }
  f[N].Im = -x[N-1] * SR2;	/* f[N] = x[N-1]*(-i)*sqrt(2) */
  free(c);			/* clean up the sines/cosines. */

  /* Do a 2N-point non-normalized FFT in place on `f[]' */
  bitrevi(f, q+1, sizeof(complex));
  W = fftomega(N);		/* generate sines/cosines */
  fftproduct(f, q+1, W);	/* apply sparse matrices */
  free(W);			/* clean up the exponentials. */
  
  /* Project back to the array `x[]' using iQ~*: */
  for(n=0, m=NN-1; n<N; n++, m--)
    {
      x[n] =  f[m].Im - f[n].Im;
    }
  free(f);			/* clean up the temp array */

  /* Normalize: divide by 2*sqrt(2N): */
  if(q&1)
    norm = 0.5/(real)(1<<(q+1)/2);
  else
    norm = 0.5/(real)sqrt((double)NN);
  for(n=0; n<N; n++)  x[n] *= norm;

  return;
}

/****************************************************************
 * dct_iv()
 *
 * This function accepts an input array of N `real' elements,
 * for N=2^q, and replaces its elements with the DCT-IV transform.
 * It uses the radix-2 FFT on 2^{q+1} points, using a temporary
 * array of 2^{q+1} extra `complex' data structs which it
 * allocates and frees upon termination.  It is its own inverse.
 */
extern void
  dct_iv(
	 real *x,	/* Input; length (1<<q). */
	 int q)		/* Nonnegative integer. */
{
  complex *f, *W, w, tmp;
  int N, NN, n, m;
  real norm, *c, *s;

  N = (1<<q);			/* number of inputs == N */
  NN = (N<<1);			/* length of temp array ==2N */

  /* Allocate and zero the temporary array */
  f = (complex *)calloc(NN, sizeof(complex));

  /* Inject the input vector into `f[]' using the R matrix: */
  c = (real *)malloc(NN*sizeof(real)); s = c+N;
  mkcossin(c, s, NN);
  f[0].Re = x[0];		/* f[0] = x[0]    */
  f[N].Im = x[N-1];		/* f[N] = ix[N-1] */
  for(n=1, m=NN-1; n<N; n++, m--)
    {
      f[n].Re = x[n] * c[n];	/* f[n] = x[n]*Wbar[n] */
      f[n].Im = - x[n] * s[n];
      f[m].Re = x[n-1] * c[n];	/* f[2N-n] = x[n-1]*W[n] */
      f[m].Im = x[n-1] * s[n];
    }
  /* keep the sines and cosines in `c[]' and `s[]' */

  /* Do a 2N-point non-normalized FFT in place on `f[]' */
  bitrevi(f, q+1, sizeof(complex));
  W = fftomega(N);		/* generate sines/cosines */
  fftproduct(f, q+1, W);	/* apply sparse matrices */
  free(W);			/* clean up the exponentials. */
  
  /* Project back to the array `x[]' using Rt: */
  w.Re = (real)cos(-PI/(N<<2));
  w.Im = (real)sin(-PI/(N<<2));

  /* x[0] = Re{ w*(f[0] + W[1]*f[2N-1]) } */
  tmp.Re = f[0].Re + CRRMULRE( f[NN-1], c[1], s[1] );
  tmp.Im = f[0].Im + CRRMULIM( f[NN-1], c[1], s[1] );
  x[0] = CCMULRE(tmp, w );

  /* x[N-1] = Re{ w*(Wbar[N-1]*f[N-1] + i*f[N]) } */
  tmp.Re = CRRMULRE(f[N-1], c[N-1], -s[N-1]) - f[N].Im;
  tmp.Im = CRRMULIM(f[N-1], c[N-1], -s[N-1]) + f[N].Re;
  x[N-1] = CCMULRE(tmp, w );

  for(n=1, m=NN-2; n<N-1; n++, m--)
    {
      /* x[n] = Re{ w*(Wbar[n]*f[n] + W[n+1]*f[2N-n-1]) } */
      tmp.Re = 
	CRRMULRE(f[n], c[n], -s[n]) 
	  + CRRMULRE(f[m], c[n+1], s[n+1]);
      tmp.Im = 
	CRRMULIM(f[n], c[n], -s[n])
	  + CRRMULIM(f[m], c[n+1], s[n+1]);
      x[n] = CCMULRE( tmp, w );
    }
  free(f);			/* clean up the temp array */
  free(c);			/* clean up the exponentials. */

  /* Normalize: divide by 2*sqrt(2N): */
  if(q&1)
    norm = 0.5/(real)(1<<(q+1)/2);
  else
    norm = 0.5/(real)sqrt((double)NN);
  for(n=0; n<N; n++)  x[n] *= norm;

  return;
}


/****************************************************************
 * dst_iv()
 *
 * This function accepts an input array of N `real' elements,
 * for N=2^q, and replaces its elements with the DST-IV transform.
 * It uses the radix-2 FFT on 2^{q+1} points, using a temporary
 * array of 2^{q+1} extra `complex' data structs which it
 * allocates and frees upon termination.  It is its own inverse.
 */
extern void
  dst_iv(
	 real *x,	/* Input; length (1<<q). */
	 int q)		/* Nonnegative integer. */
{
  complex *f, *W, w, tmp;
  int N, NN, n, m;
  real norm, *c, *s;

  N = (1<<q);			/* number of inputs == N */
  NN = (N<<1);			/* length of temp array ==2N */

  /* Allocate and zero the temporary array */
  f = (complex *)calloc(NN, sizeof(complex));

  /* Inject the input vector into `f[]' using  R~: */
  c = (real *)malloc(NN*sizeof(real)); s = c+N;
  mkcossin(c, s, NN);
  f[0].Re = x[0];		/* f[0] = x[0]     */
  f[N].Im = -x[N-1];		/* f[N] = -ix[N-1] */
  for(n=1, m=NN-1; n<N; n++, m--)
    {
      f[n].Re = x[n] * c[n];	/* f[n] = x[n]*Wbar[n] */
      f[n].Im = -x[n] * s[n];
      f[m].Re = -x[n-1] * c[n];	/* f[2N-n] = -x[n-1]*W[n] */
      f[m].Im = -x[n-1] * s[n];
    }
  /* keep the sines and cosines in `c[]' and `s[]' */

  /* Do a 2N-point non-normalized FFT in place on `f[]' */
  bitrevi(f, q+1, sizeof(complex));
  W = fftomega(N);		/* generate sines/cosines */
  fftproduct(f, q+1, W);	/* apply sparse matrices */
  free(W);			/* clean up the exponentials. */
  
  /* Project back to the array `x[]' using iw*R~t: */
  w.Re = (real)cos(-PI/(N<<2));
  w.Im = (real)sin(-PI/(N<<2));

  /* x[0] = -Im{ w*(f[0] - W[1]*f[2N-1]) } */
  tmp.Re = f[0].Re - CRRMULRE( f[NN-1], c[1], s[1] );
  tmp.Im = f[0].Im - CRRMULIM( f[NN-1], c[1], s[1] );
  x[0] = -CCMULIM(tmp, w );

  /* x[N-1] = -Im{ w*(Wbar[N-1]*f[N-1] - i*f[N]) } */
  tmp.Re = CRRMULRE(f[N-1], c[N-1], -s[N-1]) + f[N].Im;
  tmp.Im = CRRMULIM(f[N-1], c[N-1], -s[N-1]) - f[N].Re;
  x[N-1] = -CCMULIM(tmp, w );

  for(n=1, m=NN-2; n<N-1; n++, m--)
    {
      /* x[n] = -Im{ w*(Wbar[n]*f[n] - W[n+1]*f[2N-n-1]) } */
      tmp.Re = 
	CRRMULRE(f[n], c[n], -s[n]) 
	  - CRRMULRE(f[m], c[n+1], s[n+1]);
      tmp.Im = 
	CRRMULIM(f[n], c[n], -s[n])
	  - CRRMULIM(f[m], c[n+1], s[n+1]);
      x[n] = -CCMULIM( tmp, w );
    }
  free(f);			/* clean up the temp array */
  free(c);			/* clean up the exponentials. */

  /* Normalize: divide by 2*sqrt(2N): */
  if(q&1)
    norm = 0.5/(real)(1<<(q+1)/2);
  else
    norm = 0.5/(real)sqrt((double)NN);
  for(n=0; n<N; n++)  x[n] *= norm;

  return;
}

