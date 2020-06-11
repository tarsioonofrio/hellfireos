/*
 * 
 * These functions can be assembled into various kinds of factored
 * discrete Fourier transforms. 
 * 
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "complex.h"
#include "common.h"
#include "bitrev.h"
#include "fft.h"

#define fftbitrev(out,in,q)  bitrevd(out,in,q,sizeof(complex))
#define fftbrinpl(x,q)       bitrevi(x,q,sizeof(complex))

/****************************************************************
 * This function successively applies the sparse matrices
 * F_{q-1}, F_{q-2}, ..., F_1, and finally F_0 to the complex
 * input vector `f[]', transforming it in place.
 */
extern void
  fftproduct(                   /* Apply sparse matrix product. */
	     complex *f,        /* Input and output vector.     */
	     int q,             /* Length of `f[]' is N=1<<q.   */
	     const complex *W)	/* Exponentials: `Omega(N/2)' */
{
  int b, j, k, N, N1, M;
  complex  *fptr1, *fptr2, tmp;
  const complex *Wptr;

  N  = 1<<q;
  for(k=q-1; k>=0; k--)
    {
      N1 = N>>k;		/* block size */
      M = N1>>1;		/* butterfly size */
      /* Each F_k has 2^{k} blocks E_{2M}, where 2M=2^{q-k},
       * and the matrix E_M is defined by:
       *	I_M   W_M
       *	I_M  -W_M
       */
      for( b=0; b<N; b+=N1 )
	{
	  /* Multiply `f[]' by a direct summand `E_{2M}':
	   *   We use a temporary variable `tmp' to allow
	   *     multiplication in place.
	   *   We keep 2 pointers `fptr1' and `fptr2' into
	   *     the `f' array and increment them across
	   *     blocks to avoid some arithmetic.
	   *   We use the identities
	   *
	   *	     W_M[j] = W_{N/2}[j*(N/2M)]
	   *
	   *     for j=0,1,2,... so we only need
	   *     the single array Omega(N/2).
	   */
	  fptr1 = f + b;
	  fptr2 = fptr1 + M;

	  /* Do the ones: */
	  tmp.Re = fptr2->Re;	/* Mult. f[M] by +/-W[0]= +/-1 */
	  tmp.Im = fptr2->Im;
	  fptr2->Re = fptr1->Re - tmp.Re;
	  fptr2->Im = fptr1->Im - tmp.Im;
	  fptr1->Re += tmp.Re;
	  fptr1->Im += tmp.Im;

	  /* Do the cosine/sine terms */
	  Wptr = W;
	  for( j=1; j<M; j++ )
	    {
	      ++fptr1; ++fptr2;
	      Wptr += (1<<k);	/* increment W by N/2M */
	      tmp.Re = CCMULRE(*fptr2, *Wptr ); 
	      tmp.Im = CCMULIM(*fptr2, *Wptr ); 
	      fptr2->Re = fptr1->Re - tmp.Re;
	      fptr2->Im = fptr1->Im - tmp.Im;
	      fptr1->Re += tmp.Re;
	      fptr1->Im += tmp.Im;
	    }
	}
    }
  return;
}

/****************************************************************
 * Return a complex vector `Omega[]' of length `|M|' containing
 * the values:
 * 	W[n] = exp(PI*i*n/M), n = 0, 1,...,|M|-1.
 * If M<0, then the returned vector is the complex conjugate
 * of Omega(|M|).
 */ 
extern complex * 
  fftomega(			/* Return exp(-PI*i*n/M), */
	   int M)		/* for n=0,1,2,...,|M|-1. */
{
  complex *W;
  double factor, theta;
  int n;

  factor = -PI/(double)M;

  M = absval(M);
  W = (complex *)malloc(M*sizeof(complex));

  theta = 0.0;
  for(n=0; n<M; n++)
    {
      W[n].Re = (real)cos(theta);
      W[n].Im = (real)sin(theta);
      theta += factor;
    }
  return(W);
}

/****************************************************************
 * Normalize `f[]', a vector of length `N', by dividing each
 *  component by sqrt(N).
 */
extern void
  fftnormal(			/* Multiply `f[n].Re' and     */
	    complex *f,		/* `f[n].Im' by `1.0/sqrt(N), */
	    int N )		/* for n=0,1,2,...,N.         */
{
  double norm;
  int n;

  norm = sqrt(1.0/(double)N);
  for(n=0; n<N; n++)
    {
      f[n].Re = norm * f[n].Re;
      f[n].Im = norm * f[n].Im;
    }
  return;
}

/****************************************************************
 * dft()
 *
 * This function allocates, computes and returns a vector of
 * `complex' data structs which is the  discrete Fourier
 * transform of the complex input vector `f[]'.  If the
 * log_2(length) parameter `q' is negative, then it computes
 * the inverse discrete Fourier transform.
 *
 * Calling sequence:
 *	dft( f, q )
 *
 * Inputs:
 *	(const complex *) f	This is a preallocated array of 
 *				  `complex' data structs.  It is
 *				  not changed.
 *
 *	(int) q			This integer must satisfy
 *					|q| =  log_2(length f)
 *
 * Output:
 *	(complex *) dft		Function `dft()' allocates a
 *				  vector of length `1<<absval(q)'
 *				  and fills it with either:
 *				    - the DFT of `f[]'
 *					(if q >= 0 ), or
 *				    - the inverse DFT of `f[]'
 *					(if q < 0).
 *
 * External functions called:
 * 	calloc(),free()	From <stdlib.h>; allocate output vector.
 *	fftbitrev()	Permute `f[]' to `fhat[]' by bit-reversal.
 *	fftnormal()	Divide an `N'-component vector by sqrt(N).
 *	fftomega()	Return a vector of complex exponentials.
 *	fftproduct()	Apply the sparse matrices F_{q-1}...F_0.
 */
extern complex *
  dft(				/* Allocate, assign and return  */
      const complex *f,		/* a complex vector, the (1<<q) */
      int q)			/* point DFT of the input `f[]' */
{
  int N, inverse;
  complex *W, *fhat;

  /* If q<0, then the inverse DFT is to be computed: */
  inverse = 0;
  if( q < 0 )
    {
      q = -q ;
      inverse = 1;
    } 
  N = 1<<q;			/* Compute length of `f[]'. */

  /* Allocate an output vector of same length, `N', as `f[]': */
  fhat = (complex *)calloc(N, sizeof(complex));
  assert(fhat!=0);		/* Test `calloc()' success. */

  /* Compute two trivial cases directly: */
  if( q == 0 )
    {
      fhat[0].Re = f[0].Re; 
      fhat[0].Im = f[0].Im; 
      return(fhat);
    }
  if( q == 1 )
    {
      fhat[0].Re = (f[0].Re + f[1].Re)*SQH; 
      fhat[0].Im = (f[0].Im + f[1].Im)*SQH; 
      fhat[1].Re = (f[0].Re - f[1].Re)*SQH;
      fhat[1].Im = (f[0].Im - f[1].Im)*SQH;
      return(fhat);
    }
  /* Use the FFT factorization for longer vectors `f[]' */

  /* Permute by bit-reversal `fhat[] = P f[]': */
  fftbitrev( fhat, f, q);
  
  /* Generate the complex exponential vector `Omega_{N/2}': */
  if(inverse)
    W = fftomega(-N/2);		/* conjugate, for inverse DFT */
  else
    W = fftomega( N/2);		/* for DFT */

  /* Apply sparse matrices `F_{q-1} F_{q-2}...F_1 F_0 fhat[]': */
  fftproduct( fhat, q, W );

  /* Normalize `fhat[]': divide by sqrt(N): */
  fftnormal( fhat, N );

  free(W);			/* Clean up */

  return(fhat);
}


