/*
 * 
 * Functions to perform the discrete Hartley transform.
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
#include "bitrev.h"
#include "hartley.h"

#define SQH	(0.707106781186547) /* Square root of 1/2 */
#define PI	(3.1415926535897932)

#define dhtbitrev(out,in,q)  bitrevd(out,in,q,sizeof(real));
#define dhtbrinpl(x,q)       bitrevi(x,q,sizeof(real));

/****************************************************************
 * This function applies the sparse matrices G_{q-1}, G_{q-2},
 *  ... G_1, and finally G_0 to the real input vector `f[]' in
 *  place.  We assume that q>2.
 */
extern void
  dhtproduct(                   /* Apply sparse matrix product. */
             real *f,		/* Input and output vector.     */
	     int q,		/* Length of `f[]' is N=1<<q.   */
	     const real *c,	/* Cosines c_0,...,c_{N/2-1}.   */
	     const real *s)	/* Sines s_0,s_1,...,s_{N/2-1}. */
{
  int b, j, k, n, M2, M, N, N1;
  real *fptr1, *fptr2, *fptr3, *fptr4, tmp1, tmp2;

  assert(q>2);			/* Test our assumption.        */
  N= 1<<q;			/* Length of the input signal. */

  /* The first matrix G_{q-1}= A2+A2+...+A2 is a special case:
   * We apply the block diagonal matrix:
   *          1  1
   *          1 -1
   *               1  1
   *               1 -1
   *                    .
   *                      .
   *                        1  1
   *                        1 -1
   */
  for( b=0; b<N; b+=2 )
    {
      fptr1 = f + b;
      fptr2 = fptr1 + 1;

      tmp1 = *fptr1 - *fptr2;
      *fptr1 += *fptr2;
      *fptr2 = tmp1;
    }


  /* Matrices G_{q-2}, G_{q-3},..., G_0 have butterflies: */
  for(k=q-2; k>=0; k--)
    {
      N1 = N>>k;		/* size of blocks in G_k */
      M  = N1>>1;		/* butterfly size in G_k */
      M2 = M>>1;		/* midpoint of butterfly */

      /* Each G_k is the block diagonal direct sum of
       * 2^{k} blocks A_{2M}, each of size N1 x N1, where 
       * N1 = 2M = 2^{q-k} and the matrix A_{2M} is
       * defined by:
       *             I_M   B_M
       *             I_M  -B_M
       *
       * The matrix B_M is the "butterfly" below:
       *
       *      1 0       . . .       0
       *      0 c1                  s1
       *           c2            s2
       *             .         .
       *      .        cK   sK
       *      .           1
       *      .        sK  -cK
       *             .         .
       *           s2           -c2
       *      0 s1                 -c1
       *
       * for K=(M/2)-1, ck=cos(PI*k/M), sk=sin(PI*k/M).
       *
       */
      for( b=0; b<N; b+=N1 )
	{
	  /* Multiply `f[]' by a direct summand `A_{2M}':
	   *   We use temporary variables `tmp1' and `tmp2'
	   *      to allow multiplication in place.
	   *   We keep 4 pointers `fptr1' to `fptr4' into
	   *     the `f[]' array since 4 output values are 
	   *     computed simultaneously.
	   *   We use the identities for n=0,1,2,...,M/2-1:
	   *
	   *	     cos(n*PI/M) = cos(n*(N/2M)*P/(N/2))
	   *	     sin(n*PI/M) = sin(n*(N/2M)*P/(N/2))
	   *
	   *     so we only need one table of length N/4 
	   *     for each of sines and cosines. Note that 
	   *     (N/2M) = (1<<k).
	   */

	  /* Do the ones: */
	  fptr1 = f + b;
	  fptr2 = fptr1 + M2;
	  fptr3 = fptr1 + M;
	  fptr4 = fptr3 + M2;
	  
	  tmp1 = *fptr1 - *fptr3;
	  *fptr1 += *fptr3;
	  *fptr3 = tmp1;

	  tmp2 = *fptr2 - *fptr4;
	  *fptr2 += *fptr4;
	  *fptr4 = tmp2;

	  /* Do the sines/cosines: */
	  fptr4 += M2;
	  fptr2 += M2;
	  n = 0;
	  for( j=1; j<M2; j++ )
	    {
	      ++fptr1; --fptr2;
	      ++fptr3; --fptr4;
	      n += (1<<k);

	      tmp1 = *fptr3*c[n] + *fptr4*s[n];
	      tmp2 = *fptr3*s[n] - *fptr4*c[n];
	      *fptr3 = *fptr1 - tmp1;
	      *fptr4 = *fptr2 - tmp2;
	      *fptr1 += tmp1;
	      *fptr2 += tmp2;
	    }
	}
    }
  return;
}

/****************************************************************
 * Fill the preallocated arrays of real's `c[]' and `s[]', each
 * assumed to be of length `N/2', with the respective  values:
 * 	c[n] = cos(n*PI/N), n=0,1,...,N/2-1.
 * 	s[n] = sin(n*PI/N), n=0,1,...,N/2-1.
 */ 
extern void 
  dhtcossin(			/* Fill the arrays with:  */
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
 * Normalize `f[]', a vector of length `N', by dividing each
 *  component by sqrt(N).
 */
extern void
  dhtnormal(
	    real *f,		/* Multiply `f[n]' for n=0,1, */
	    int N )		/* ...,N  by `1.0/sqrt(N).'   */
{
  double norm;
  int n;
  
  norm = sqrt(1.0/(double)N);
  for(n=0; n<N; n++)
    {
      f[n] = norm * f[n];
    }
  return;
}

/****************************************************************
 * dht()
 *
 * This function allocates, computes and returns an array of
 * real's which is the discrete Hartley transform of the real
 * input vector `f[]'.
 *
 * Calling sequence:
 *	dht( f, q )
 *
 * Inputs:
 *	(const real *) f	This is a preallocated array of 
 *				  `real's.  It is not changed.
 *
 *	(int) q			This nonnegative integer must
 *				  satisfy  q == log_2(length f)
 *
 * Output:
 *	(real *) dht		Function `dht()' allocates a
 *				  vector of length `1<<q'
 *				  and fills it with the DHT
 *				  of `f[]'.
 *
 * External functions called:
 * 	malloc(),free()	From <stdlib.h>; to allocate arrays.
 *	dhtbitrev()	Permute `f[]' to `fhat[]' by bit-reversal.
 *	dhtnormal()	Divide an `N'-component vector by sqrt(N).
 *	dhtcossin()	Fill arrays with sines and cosines..
 *	dhtproduct()	Apply the sparse matrices F_{q-1}...F_0.
 */
extern real *
  dht(				/* Allocate, assign and return  */
      const real *f,		/* a real vector, the (1<<q)   */
      int q)			/* point DFT of the input `f[]' */
{
  int N;
  real *fh, *c, *s;

  assert(q>=0);			/* Accept only positive `q'. */
  N = 1<<q;			/* Compute length of `f[]'.  */

  /* Allocate an output vector of same length, `N', as `f[]': */
  fh = (real *)malloc(N*sizeof(real));
  assert(fh!=0);		/* Test `calloc()' success. */

  /* Compute three trivial cases directly: */
  if( q == 0 )
    {
      fh[0] = f[0]; 
      return(fh);
    }
  if( q == 1 )
    {
      fh[0] = (f[0] + f[1])*SQH; 
      fh[1] = (f[0] - f[1])*SQH;
      return(fh);
    }
  if( q == 2 )
    {
      fh[0] = (f[0] + f[1] + f[2] + f[3])*0.5; 
      fh[1] = (f[0] + f[1] - f[2] - f[3])*0.5;
      fh[2] = (f[0] - f[1] + f[2] - f[3])*0.5; 
      fh[3] = (f[0] - f[1] - f[2] + f[3])*0.5;
      return(fh);
    }
  /* Use the DHT factorization for longer vectors `f[]'. */

  /* Permute by bit-reversal `fh[] = P f[]': */
  dhtbitrev(fh, f, q);

  /* Generate the table of sines and cosines: */
  c = (real *)malloc((N/2)*sizeof(real));
  s = c + (N/4);
  dhtcossin( c, s, N/2);

  /* Apply sparse matrices `G_{q-1} G_{q-2}...G_1 G_0 fh[]': */
  dhtproduct( fh, q, c, s );

  /* Normalize `fh[]': divide by sqrt(N): */
  dhtnormal(fh, N);

  /* Clean up by deallocating the sine/cosine table. */
  free(c);

  return(fh);
}


