/* 
 *  Discrete wavelet transform functions for one-dimensional signals.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 * 
 */

#include <assert.h>
#include <stdlib.h>
#include "real.h"
#include "interval.h"
#include "dwt.h"
#include "common.h"		/* for `max()' and `min()'. */
#include "cd.h"

/*********************************************************
 * dwtpd0()
 *
 * [D]iscrete [W]avelet [T]ransform on [P]eriodic
 * [D]isjoint arrays, down to level [0]; Recursive.
 *
 * Calling sequence and basic algorithm:
 *
 *  dwtpd0( DIFS, SUMS, IN, N, H, G ):
 *     If N > 1 then
 *        cdpi( DIFS+N/2, 1, IN, N, G )
 *        cdpi( SUMS+N/2, 1, IN, N, H )
 *        dwtpd0( DIFS, SUMS, SUMS+N/2, N/2, H, G )
 *     Else
 *        Let DIFS[0] += IN[0]
 *
 * Inputs:
 *	(real *)difs		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(real *)sums		This preallocated array must have
 *				  at least `N' zeroes.
 *
 *	(const real *)in	This input array must have at least
 *				  `N' elements.  It is not changed.
 *
 *	(int)N			This must be 2^L for some L>=0.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)difs		This array is filled by reference with
 *				  wavelet coefficients.
 *
 *	(real *)sums		This array is filled by reference with
 *				  scaling function coefficients.
 *
 * Assumptions:
 *	1. `difs[]' and `sums[]' are disjoint and disjoint from `in[]'.
 *	2. `N' is a nonnegative integer power of 2.
 *	3. `sums[0]...sums[N-1]' are 0 at the outset.
 *
 * External functions called:
 *	cdpi()
 */
extern void
  dwtpd0( 
	 real *difs,     /* Output array with N elements.  */
	 real *sums,     /* Working array with N zeroes.   */
	 const real *in, /* Input array with N elements.   */
	 int N,          /* This is the number of inputs.  */
	 const pqf *H,   /* Low-pass QF data structure.    */
	 const pqf *G)   /* High-pass QF data structure.   */
{
  int N2;

  if( N>1 )
    {
      N2 = N/2;
      cdpi( difs+N2, 1, in, N, G );
      cdpi( sums+N2, 1, in, N, H );
      dwtpd0( difs, sums, sums+N2, N2, H, G );
    }
  else
    difs[0] += in[0];

  return;
}

/*********************************************************
 * dwtpd0n():
 *
 * [D]iscrete [W]avelet [T]ransform on [P]eriodic
 * [D]isjoint arrays, down to level [0]; [N]on-recursive.
 *
 *  Calling sequence and basic algorithm:
 *
 *  dwtpd0n( DIFS, SUMS, IN, N, H, G ):
 *     While N > 1
 *        cdpi( DIFS+(N/2), 1, IN, N, G )
 *        cdpi( SUMS+(N/2), 1, IN, N, H )
 *        Let IN = SUMS+(N/2)
 *        N /= 2
 *     Let DIFS[0] += IN[0]
 *
 * Inputs:
 *	(real *)difs		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(real *)sums		This preallocated array must have
 *				  at least `N' zeroes.
 *
 *	(const real *)in	This input array must have at least
 *				  `N' elements.  It is not changed.
 *
 *	(int)N			This must be 2^L for some L>=0.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)difs		This array is filled by reference with
 *				  wavelet coefficients.
 *
 *	(real *)sums		This array is filled by reference with
 *				  scaling function coefficients.
 *
 * Assumptions:
 *	1. `difs[]' and `sums[]' are disjoint and disjoint from `in[]'.
 *	2. `N' is a nonnegative integer power of 2.
 *	3. `sums[0]...sums[N-1]' are 0 at the outset.
 *
 * External functions called:
 *	cdpi()
 */
extern void
  dwtpd0n( 
	  real *difs,     /* Output array with N elements.  */
	  real *sums,     /* Working array with N zeroes.   */
	  const real *in, /* Input array with N elements.   */
	  int N,          /* This is the number of inputs.  */
	  const pqf *H,   /* Low-pass QF data structure.    */
	  const pqf *G)   /* High-pass QF data structure.   */
{
  int N2;

  while( N > 1 )
    {
      N2 = N/2;
      cdpi( difs+N2, 1, in, N, G );
      cdpi( sums+N2, 1, in, N, H );
      in = sums + N2;
      N = N2;
    }
  difs[0] += in[0];
  return;
}

/*********************************************************
 * dwtpi0():
 *
 * [D]iscrete [W]avelet [T]ransform on [P]eriodic
 * [I]n-place arrays, down to level [0]; Recursive.
 *
 * Calling sequence and basic algorithm:
 *
 *  dwtpi0( DATA, WORK, N, H, G ):
 *     If N > 1 then
 *        cdpe( WORK+N/2, 1, DATA, N, G )
 *        cdpe( WORK, 1, DATA, N, H )
 *        For I = 0 to N-1
 *           Let DATA[I] = WORK[I]
 *        dwtpi0( DATA, WORK, N/2, H, G )
 *
 * Inputs:
 *	(real *)data		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(real *)work		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(int)N			This must be 2^L for some L>=0.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)data		This array is overwritten by reference with
 *				  wavelet coefficients.
 *
 *	(real *)work		This array is overwritten by reference with
 *				  scaling function coefficients.
 *
 * Assumptions:
 *	1. `data[]' and `work[]' are disjoint arrays.
 *	2. `N' is a nonnegative integer power of 2.
 *
 * External functions called:
 *	cdpe()
 */
extern void
  dwtpi0( 
	 real *data,     /* Input/Output array, length N.  */
	 real *work,     /* Working array, length N.       */
	 int N,          /* This is the number of inputs.  */
	 const pqf *H,   /* Low-pass QF data structure.    */
	 const pqf *G)   /* High-pass QF data structure.   */
{
  int i;

  if( N > 1 )
    {
      cdpe( work+N/2, 1, data, N, G );
      cdpe(  work,   1, data, N, H );
      for(i=0; i<N; i++) data[i] = work[i];
      dwtpi0( data, work, N/2, H, G );
    }
  return;
}

/*********************************************************
 * dwtpi0n():
 *
 * [D]iscrete [W]avelet [T]ransform on [P]eriodic
 * [I]n-place arrays, down to level [0]; [N]on-recursive.
 *
 *  Calling sequence and basic algorithm:
 *
 *  dwtpi0n( DATA, WORK, N, H, G ):
 *     While N > 1
 *        cdpe( WORK+N/2, 1, DATA, N, G )
 *        cdpe( WORK, 1, DATA, N, H )
 *        For I = 0 to N-1
 *           Let DATA[I] = WORK[I]
 *        N /= 2
 *
 *
 * Inputs:
 *	(real *)data		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(real *)work		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(int)N			This must be 2^L for some L>=0.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)data		This array is overwritten by reference with
 *				  wavelet coefficients.
 *
 *	(real *)work		This array is overwritten by reference with
 *				  scaling function coefficients.
 *
 * Assumptions:
 *	1. `data[]' and `work[]' are disjoint arrays.
 *	2. `N' is a nonnegative integer power of 2.
 *
 * External functions called:
 *	cdpe()
 */
extern void
  dwtpi0n( 
	  real *data,     /* Input/Output array, length N.  */
	  real *work,     /* Working array, length N.       */
	  int N,          /* This is the number of inputs.  */
	  const pqf *H,   /* Low-pass QF data structure.    */
	  const pqf *G)   /* High-pass QF data structure.   */
{
  int i;

  while( N > 1 )
    {
      cdpe( work+(N/2), 1, data, N, G );
      cdpe(    work,    1, data, N, H );
      for(i=0; i<N; i++) data[i]=work[i];
      N /= 2;
    }
  return;
}

/*********************************************************
 * dwtpd():
 *
 * [D]iscrete [W]avelet [T]ransform on [P]eriodic
 * [D]isjoint arrays to a specified level; Recursive.
 *
 *  Calling sequence and basic algorithm:
 *
 *  dwtpd( DIFS, SUMS, IN, N, L, H, G ):
 *     If L > 0 then
 *        cdpo( DIFS+N/2, 1, IN, N, G )
 *        cdpo( SUMS+N/2, 1, IN, N, H )
 *        dwtpd( DIFS, SUMS, SUMS+N/2, N/2, L-1, H, G )
 *     Else
 *        For K = 0 to N-1
 *           DIFS[K] += IN[K]
 *
 * Inputs:
 *	(real *)difs		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(real *)sums		This preallocated array must have
 *				  at least `N' zeroes.
 *
 *	(const real *)in	This input array must have at least
 *				  `N' elements.  It is not changed.
 *
 *	(int)N			This must be `M*2^L' for some M>0.
 *
 *	(int)L			This must be a nonnegative integer.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)difs		This array is filled by reference with
 *				  wavelet coefficients.
 *
 *	(real *)sums		This array is filled by reference with
 *				  scaling function coefficients.
 *
 * Assumptions:
 *	1. `difs[]' and `sums[]' are disjoint and disjoint from `in[]'.
 *	2. `N' is divisible by `1<<L'.
 *	3. `L' is a nonnegative integer.
 *	4. `sums[0]...sums[N-1]' are 0 at the outset.
 *
 * External functions called:
 *	cdpo()
 */
extern void
  dwtpd( 
	real *difs,       /* Output array with N elements.  */
	real *sums,       /* Working array with N zeroes.   */
	const real *in,   /* Input array with N elements.   */
	int N,            /* This is the number of inputs.  */
	int L,            /* This is the number of levels.  */
	const pqf *H,     /* Low-pass QF data structure.    */
	const pqf *G)     /* High-pass QF data structure.   */
{
  int N2;

  if( L>0 )
    {
      N2 = N/2;
      cdpo( difs+N2, 1, in, N, G );
      cdpo( sums+N2, 1, in, N, H );
      dwtpd( difs, sums, sums+N2, N2, L-1, H, G );
    }
  else
    while( N-- > 0 ) *difs++ += *in++;
  return;
}

/*********************************************************
 * dwtpdn():
 *
 * [D]iscrete [W]avelet [T]ransform on [P]eriodic
 * [D]isjoint arrays to a specified level; [N]on-recursive.
 *
 *  Calling sequence and basic algorithm:
 *
 *  dwtpdn( DIFS, SUMS, IN, N, L, H, G ):
 *     While L > 0
 *        cdpi( DIFS+N/2, 1, IN, N, G )
 *        cdpi( SUMS+N/2, 1, IN, N, H )
 *        L -= 1
 *        N /= 2
 *        Let IN = SUMS + N
 *     For K = 0 to N-1
 *        DIFS[K] += IN[K]
 *
 * Calling sequence:
 *	dwtpdn( difs, sums, in, N, L, H, G )
 *
 * Inputs:
 *	(real *)difs		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(real *)sums		This preallocated array must have
 *				  at least `N' zeroes.
 *
 *	(const real *)in	This input array must have at least
 *				  `N' elements.  It is not changed.
 *
 *	(int)N			This must be `M*2^L' for some M>0.
 *
 *	(int)L			This must be a nonnegative integer.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)difs		This array is filled by reference with
 *				  wavelet coefficients.
 *
 *	(real *)sums		This array is filled by reference with
 *				  scaling function coefficients.
 *
 * Assumptions:
 *	1. `difs[]' and `sums[]' are disjoint and disjoint from `in[]'.
 *	2. `N' is divisible by `1<<L'.
 *	3. `L' is a nonnegative integer.
 *	4. `sums[0]...sums[N-1]' are 0 at the outset.
 *
 * External functions called:
 *	cdpi()
 */
extern void
  dwtpdn( 
	 real *difs,      /* Output array with N elements.  */
	 real *sums,      /* Working array with N zeroes.   */
	 const real *in,  /* Input array with N elements.   */
	 int N,           /* This is the number of inputs.  */
	 int L,           /* This is the number of levels.  */
	 const pqf *H,    /* Low-pass QF data structure.    */
	 const pqf *G)    /* High-pass QF data structure.   */
{
  int N2;

  while( L-- >0 )
    {
      N2 = N/2;
      cdpi( difs+N2, 1, in, N, G );
      cdpi( sums+N2, 1, in, N, H );
      in = sums+N2;
      N = N2;
    }
  while( N-- > 0 ) *difs++ += *in++;
  return;
}

/*********************************************************
 * dwtpi():
 *
 * [D]iscrete [W]avelet [T]ransform on [P]eriodic
 * [I]n-place arrays, to a specified level; Recursive.
 *
 * Calling sequence and basic algorithm:
 *
 *  dwtpi( DATA, WORK, N, L, H, G ):
 *     If L > 0 then
 *        cdpe( WORK+N/2, 1, DATA, N, G )
 *        cdpe( WORK, 1, DATA, N, H )
 *        For I = 0 to N-1
 *           Let DATA[I] = WORK[I]
 *        dwtpi( DATA, WORK, N/2, L-1, H, G )
 *
 *
 * Inputs:
 *	(real *)data		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(real *)work		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(int)N			This must be `M*2^L' for some M>0.
 *
 *	(int)L			This must be a nonnegative integer.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)data		This array is overwritten by reference with
 *				  wavelet coefficients.
 *
 *	(real *)work		This array is overwritten by reference with
 *				  scaling function coefficients.
 *
 * Assumptions:
 *	1. `data[]' and `work[]' are disjoint arrays.
 *	2. `N' is divisible by `1<<L'.
 *	3. `L' is a nonnegative integer.
 *
 * External functions called:
 *	cdpe()
 */
extern void
  dwtpi( 
	real *data,     /* Input/Output array, length N.  */
	real *work,     /* Working array, length N.       */
	int N,          /* This is the number of inputs.  */
	int L,          /* This is the number of levels.  */
	const pqf *H,   /* Low-pass QF data structure.    */
	const pqf *G)   /* High-pass QF data structure.   */
{
  int i;

  if( L > 0 )
    {
      cdpe( work+(N/2), 1, data, N, G );
      cdpe(    work,    1, data, N, H );
      for(i=0; i<N; i++) data[i] = work[i];
      dwtpi( data, work, N/2, L-1, H, G );
    }
  return;
}

/*********************************************************
 * dwtpin():
 *
 * [D]iscrete [W]avelet [T]ransform on [P]eriodic
 * [I]n-place arrays, to a specified level; [N]on-recursive.
 *
 * Calling sequence and basic algorithm:
 *
 *  dwtpin( DATA, WORK, N, L, H, G ):
 *     While L > 0
 *        cdpe( WORK+N/2, 1, DATA, N, G )
 *        cdpe( WORK, 1, DATA, N, H )
 *        For I = 0 to N-1
 *           Let DATA[I] = WORK[I]
 *        N /= 2
 *        L -= 1
 *
 *
 * Inputs:
 *	(real *)data		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(real *)work		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(int)N			This must be `M*2^L' for some M>0.
 *
 *	(int)L			This must be a nonnegative integer.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)data		This array is overwritten by reference with
 *				  wavelet coefficients.
 *
 *	(real *)work		This array is overwritten by reference with
 *				  scaling function coefficients.
 *
 * Assumptions:
 *	1. `data[]' and `work[]' are disjoint arrays.
 *	2. `N' is divisible by `1<<L'.
 *	3. `L' is a nonnegative integer.
 *
 * External functions called:
 *	cdpe()
 */
extern void
  dwtpin( 
	 real *data,     /* Input/Output array, length N.  */
	 real *work,     /* Working array, length N.       */
	 int N,          /* This is the number of inputs.  */
	 int L,          /* This is the number of levels.  */
	 const pqf *H,   /* Low-pass QF data structure.    */
	 const pqf *G)   /* High-pass QF data structure.   */
{
  int i;

  while( L-- > 0 )
    {
      cdpe( work+(N/2), 1, data, N, G );
      cdpe(    work,    1, data, N, H );
      for(i=0; i<N; i++) data[i] = work[i];
      N /= 2;
    }
  return;
}

/*********************************************************
 * idwtpd0():
 *
 * [I]nverse [D]iscrete [W]avelet [T]ransform on [P]eriodic
 * [D]isjoint arrays, down to level [0]; Recursive.
 *
 * Calling sequence and basic algorithm:
 *
 *  idwtpd0( OUT, SUMS, IN, N, H, G ):
 *     If N > 1 then
 *        N /=2
 *        idwtpd0( SUMS+N, SUMS, IN, N, H, G)
 *        acdpi( OUT, 1, SUMS+N, N, H )
 *        acdpi( OUT, 1,  IN +N, N, G )
 *     Else
 *        Let OUT[0] += IN[0]
 *
 *
 * Inputs:
 *	(real *)out		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(real *)sums		This preallocated array must have
 *				  at least `N' zeroes.
 *
 *	(const real *)in	This input array must have at least
 *				  `N' elements.  It is not changed.
 *
 *	(int)N			This must be 2^L for some L>=0.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)out		This array is filled by reference with
 *				  the reconstruction from wavelet coefficients.
 *
 *	(real *)sums		This array is filled by reference with
 *				  scaling function coefficients.
 *
 * Assumptions:
 *	1. `out[]' and `sums[]' are disjoint and disjoint from `in[]'.
 *	2. `N' is a nonnegative integer power of 2.
 *	3. `sums[0]...sums[N-1]' are 0 at the outset.
 *
 * External functions called:
 *	acdpi()
 */
extern void
  idwtpd0( 
	  real *out,       /* Preallocated, zeroed output.  */
	  real *sums,      /* Working array of N zeroes.    */
	  const real *in,  /* Input array with N elements.  */
	  int N,           /* This is the number of inputs. */
	  const pqf *H,    /* Low-pass QF data structure.   */
	  const pqf *G)    /* High-pass QF data structure.  */
{
  if( N > 1 )
    {
      N /= 2;
      idwtpd0( sums+N, sums, in, N, H, G);
      acdpi( out, 1, sums+N, N, H );
      acdpi( out, 1,  in +N, N, G );
    }
  else
    *out += *in;
  return;
}

/*********************************************************
 * idwtpd0n():
 *
 * [I]nverse [D]iscrete [W]avelet [T]ransform on [P]eriodic
 * [D]isjoint arrays, down to level [0]; [N]on-recursive.
 *
 * Calling sequence and basic algorithm:
 *
 *  idwtpd0n( OUT, SUMS, IN, N, H, G ):
 *     If N > 1 then
 *        SUMS += 1
 *        Let M = 1
 *        Let SUMS[0] = IN[0]
 *        IN += 1
 *        While M < N/2
 *           acdpo( SUMS+M, 1, SUMS, M, H )
 *           acdpo( SUMS+M, 1,  IN,  M, G )
 *           SUMS += M
 *           IN   += M
 *           M *= 2
 *        acdpo( OUT, 1, SUMS, N/2, H )
 *        acdpo( OUT, 1,  IN,  N/2, G )
 *     Else
 *        OUT[0] += IN[0]
 *
 *
 * Inputs:
 *	(real *)out		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(real *)sums		This preallocated array must have
 *				  at least `N' zeroes.
 *
 *	(const real *)in	This input array must have at least
 *				  `N' elements.  It is not changed.
 *
 *	(int)N			This must be 2^L for some L>=0.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)out		This array is filled by reference with
 *				  the reconstruction from wavelet coefficients.
 *
 *	(real *)sums		This array is filled by reference with
 *				  scaling function coefficients.
 *
 * Assumptions:
 *	1. `out[]' and `sums[]' are disjoint and disjoint from `in[]'.
 *	2. `N' is a nonnegative integer power of 2.
 *	3. `sums[0]...sums[N-1]' are 0 at the outset.
 *
 * External functions called:
 *	acdpo()
 */
extern void
  idwtpd0n( 
	   real *out,      /* Preallocated, zeroed output.  */
	   real *sums,     /* Working array of N zeroes.    */
	   const real *in, /* Input array with N elements.  */
	   int N,          /* This is the number of inputs. */
	   const pqf *H,   /* Low-pass QF data structure.   */
	   const pqf *G)   /* High-pass QF data structure.  */
{
  int M;

  if( N > 1 )
    {
      M = 1;
      sums += M;
      sums[0] = in[0];
      in += 1;
      while( M < N/2 )
	{
	  acdpo( sums+M, 1, sums, M, H );
	  acdpo( sums+M, 1,  in,  M, G );
	  sums += M;
	  in += M;
	  M *= 2;
	}
      acdpo( out, 1, sums, M, H );
      acdpo( out, 1,  in,  M, G );
    }
  else
    *out += *in;
  return;
}

/*********************************************************
 * idwtpi0():
 *
 * [I]nverse [D]iscrete [W]avelet [T]ransform on [P]eriodic
 * [I]n-place arrays, to level [0]; Recursive.
 *
 * Calling sequence and basic algorithm:
 *
 *  idwtpi0( DATA, WORK, N, L, H, G ):
 *     If N > 1 then
 *        idwtpi0( DATA, WORK, N/2, H, G )
 *        acdpe( WORK, 1,    DATA,    N/2, H )
 *        acdpo( WORK, 1, DATA+(N/2), N/2, G )
 *        For I = 0 to N-1
 *           Let DATA[I] = WORK[I]
 *
 *
 * Inputs:
 *	(real *)data		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(real *)work		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(int)N			This must be 2^L for some L>=0.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)data		This array is filled by reference with
 *				  the reconstruction from wavelet coefficients.
 *
 *	(real *)work		This array is filled by reference with
 *				  scaling function coefficients.
 *
 * Assumptions:
 *	1. `work[]' and `data[]' are disjoint.
 *	2. `N == 1<<L' for some `L>=0'.
 *
 * External functions called:
 *	acdpe(), acdpo()
 */
extern void
  idwtpi0( 
	  real *data,     /* Input and output, length N.    */
	  real *work,     /* Working array of N elements.   */
	  int N,          /* This is the number of inputs.  */
	  const pqf *H,   /* Low-pass QF data structure.    */
	  const pqf *G)   /* High-pass QF data structure.   */
{
  int N2, i;

  if( N>1 )
    {
      N2 = N/2;
      idwtpi0 ( data, work, N2, H, G );
      acdpe ( work, 1,   data,  N2, H );
      acdpo ( work, 1, data+N2, N2, G );
      for(i=0; i<N; i++) data[i] = work[i];
    }
  return;
}

/*********************************************************
 * idwtpi0n():
 *
 * [I]nverse [D]iscrete [W]avelet [T]ransform on [P]eriodic
 * [I]n-place arrays, to level [0]; [N]on-recursive.
 *
 * Calling sequence and basic algorithm:
 *
 *  idwtpi0n( DATA, WORK, N, H, G ):
 *     Let M = 1
 *     While M < N
 *        acdpe( WORK, 1,  DATA,  M, H )
 *        acdpo( WORK, 1, DATA+M, M, G )
 *        For I = 0 to 2*M-1
 *           Let DATA[I] = WORK[I]
 *        M *= 2
 *
 *
 * Inputs:
 *	(real *)data		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(real *)work		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(int)N			This must be 2^L for some L>=0.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)data		This array is filled by reference with
 *				  the reconstruction from wavelet coefficients.
 *
 *	(real *)work		This array is filled by reference with
 *				  scaling function coefficients.
 *
 * Assumptions:
 *	1. `work[]' and `data[]' are disjoint.
 *	2. `N' is divisible by `1<<L'.
 *	3. `L' is a nonnegative integer.
 *
 * External functions called:
 *	acdpe(), acdpo()
 */
extern void
  idwtpi0n( 
	   real *data,     /* Input and output, length N.    */
	   real *work,     /* Working array of N elements.   */
	   int N,          /* This is the number of inputs.  */
	   const pqf *H,   /* Low-pass QF data structure.    */
	   const pqf *G)   /* High-pass QF data structure.   */
{
  int  M, i;

  M = 1;
  while( M<N )
    {
      acdpe ( work, 1,  data,  M, H );
      acdpo ( work, 1, data+M, M, G );
      M *= 2;
      for(i=0; i<M; i++) data[i] = work[i];
    }
  return;
}

/*********************************************************
 * idwtpd():
 *
 * [I]nverse [D]iscrete [W]avelet [T]ransform on [P]eriodic
 * [D]isjoint arrays, to specified level; Recursive.
 *
 * Calling sequence and basic algorithm:
 *
 *  idwtpd( OUT, SUMS, IN, N, L, H, G ):
 *     If L > 0 then
 *        N /= 2
 *        idwtpd( SUMS+N, SUMS, IN, N, L-1, H, G )
 *        acdpo( OUT, 1, SUMS+N, N, H )
 *        acdpo( OUT, 1,  IN +N, N, G )
 *     Else
 *        For K = 0 to N-1
 *           Let OUT[K] += IN[K]
 *
 *
 * Inputs:
 *	(real *)out		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(real *)sums		This preallocated array must have
 *				  at least `N' zeroes.
 *
 *	(const real *)in	This input array must have at least
 *				  `N' elements.  It is not changed.
 *
 *	(int)N			This must be `M*2^L' for some M>0.
 *
 *	(int)L			This must be a nonnegative integer.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)out		This array is filled by reference with
 *				  the reconstruction from wavelet coefficients.
 *
 *	(real *)sums		This array is filled by reference with
 *				  scaling function coefficients.
 *
 * Assumptions:
 *	1. `out[]' and `sums[]' are disjoint and disjoint from `in[]'.
 *	2. `N' is divisible by `1<<L'.
 *	3. `L' is a nonnegative integer.
 *	4. `sums[0]...sums[N-1]' are 0 at the outset.
 *
 * External functions called:
 *	acdpo()
 */
extern void
  idwtpd( 
	 real *out,        /* Preallocated, zeroed output.  */
	 real *sums,       /* Working array of N zeroes.    */
	 const real *in,   /* Input array with N elements.  */
	 int N,            /* This is the number of inputs. */
	 int L,            /* This is the number of levels. */
	 const pqf *H,     /* Low-pass QF data structure.   */
	 const pqf *G)     /* High-pass QF data structure.  */
{
  if( L > 0 )
    {
      N /= 2;
      idwtpd( sums+N, sums, in, N, L-1, H, G);
      acdpo( out, 1, sums+N, N, H );
      acdpo( out, 1,  in +N, N, G );
    }
  else
    while( N-- > 0 ) *out++ += *in++;
  return;
}

/*********************************************************
 * idwtpdn():
 *
 * [I]nverse [D]iscrete [W]avelet [T]ransform on [P]eriodic
 * [D]isjoint arrays, to specified level; [N]on-recursive.
 *
 * Calling sequence and basic algorithm:
 *
 *  idwtpdn( OUT, SUMS, IN, N, L, H, G ):
 *     If L > 0 then
 *        M = N>>L
 *        SUMS += M
 *        For K=0 to M-1
 *           Let SUMS[K] = IN[K]
 *        IN += M
 *        While M < N/2
 *           acdpi( SUMS+M, 1, SUMS, M, H )
 *           acdpi( SUMS+M, 1,  IN,  M, G )
 *           SUMS += M
 *           IN   += M
 *           M  *= 2
 *        acdpi( OUT, 1, SUMS, N/2, H )
 *        acdpi( OUT, 1, IN,   N/2, G )
 *     Else
 *        For K=0 to N-1
 *           Let OUT[K] += IN[K]
 *
 *
 * Inputs:
 *	(real *)out		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(real *)sums		This preallocated array must have
 *				  at least `N' zeroes.
 *
 *	(const real *)in	This input array must have at least
 *				  `N' elements.  It is not changed.
 *
 *	(int)N			This must be `M*2^L' for some M>0.
 *
 *	(int)L			This must be a nonnegative integer.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)out		This array is filled by reference with
 *				  the reconstruction from wavelet coefficients.
 *
 *	(real *)sums		This array is filled by reference with
 *				  scaling function coefficients.
 *
 * Assumptions:
 *	1. `out[]' and `sums[]' are disjoint and disjoint from `in[]'.
 *	2. `N' is divisible by `1<<L'.
 *	3. `L' is a nonnegative integer.
 *	4. `sums[0]...sums[N-1]' are 0 at the outset.
 *
 * External functions called:
 *	acdpi(),memcpy()
 */
extern void
  idwtpdn( 
	  real *out,      /* Output, N preallocated zeroes. */
	  real *sums,     /* Working array of N zeroes.     */
	  const real *in, /* Input array with N elements.   */
	  int N,          /* This is the number of inputs.  */
	  int L,          /* This is the number of levels.  */
	  const pqf *H,   /* Low-pass QF data structure.    */
	  const pqf *G)   /* High-pass QF data structure.   */
{
  int M;

  if( L>0 )
    {
      M = N>>L;
      sums += M;
      memcpy( sums, in, M*sizeof(real) );
      in += M;
      while( M < N/2 )
	{
	  acdpi( sums+M, 1, sums, M, H );
	  acdpi( sums+M, 1,  in,  M, G );
	  sums += M;
	  in += M;
	  M *= 2;
	}
      acdpi( out, 1, sums, M, H );
      acdpi( out, 1,  in,  M, G );
    }
  else
    while( N-- > 0 ) *out++ += *in++;
  return;
}

/*********************************************************
 * idwtpi():
 *
 * [I]nverse [D]iscrete [W]avelet [T]ransform on [P]eriodic
 * [I]n-place arrays, to specified level; Recursive.
 *
 * Calling sequence and basic algorithm:
 *
 *  idwtpi( DATA, WORK, N, L, H, G ):
 *     If L > 0 then
 *        Let M = N/2
 *        idwtpi( DATA, WORK, M, L-1, H, G )
 *        acdpe( WORK, 1,  DATA,  M, H )
 *        acdpo( WORK, 1, DATA+M, M, G )
 *        For I = 0 to N-1
 *           Let DATA[I] = WORK[I]
 *
 *
 * Inputs:
 *	(real *)data		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(real *)work		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(int)N			This must be `M*2^L' for some M>0.
 *
 *	(int)L			This must be a nonnegative integer.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)data		This array is filled by reference with
 *				  the reconstruction from wavelet coefficients.
 *
 *	(real *)work		This array is filled by reference with
 *				  scaling function coefficients.
 *
 * Assumptions:
 *	1. `work[]' and `data[]' are disjoint.
 *	2. `N' is divisible by `1<<L'.
 *	3. `L' is a nonnegative integer.
 *
 * External functions called:
 *	acdpe(), acdpo()
 */
extern void
  idwtpi( 
	 real *data,     /* Input and output, length N.    */
	 real *work,     /* Working array of N elements.   */
	 int N,          /* This is the number of inputs.  */
	 int L,          /* This is the number of levels.  */
	 const pqf *H,   /* Low-pass QF data structure.    */
	 const pqf *G)   /* High-pass QF data structure.   */
{
  int  M, i;

  if( L>0 )
    {
      M = N/2;
      idwtpi ( data, work, M, L-1, H, G );
      acdpe ( work, 1,  data,  M, H );
      acdpo ( work, 1, data+M, M, G );
      for(i=0; i<N; i++) data[i] = work[i];
    }
  return;
}

/*********************************************************
 * idwtpin():
 *
 * [I]nverse [D]iscrete [W]avelet [T]ransform on [P]eriodic
 * [I]n-place arrays, to specified level; [N]on-recursive.
 *
 * Calling sequence and basic algorithm:
 *
 *  idwtpin( DATA, WORK, N, L, H, G ):
 *     Let M = N>>L
 *     While M < N
 *        acdpi( WORK, 1,  DATA,  M, H )
 *        acdpi( WORK, 1, DATA+M, M, G )
 *        For I = 0 to 2*M-1
 *           Let DATA[I] = WORK[I]
 *        M *= 2
 *
 *
 * Inputs:
 *	(real *)data		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(real *)work		This preallocated array must have
 *				  at least `N' elements.
 *
 *	(int)N			This must be `M*2^L' for some M>0.
 *
 *	(int)L			This must be a nonnegative integer.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)data		This array is filled by reference with
 *				  the reconstruction from wavelet coefficients.
 *
 *	(real *)work		This array is filled by reference with
 *				  scaling function coefficients.
 *
 * Assumptions:
 *	1. `work[]' and `data[]' are disjoint.
 *	2. `N' is divisible by `1<<L'.
 *	3. `L' is a nonnegative integer.
 *
 * External functions called:
 *	acdpe(), acdpo()
 */
extern void
  idwtpin( 
	  real *data,     /* Input and output, length N.    */
	  real *work,     /* Working array of N elements.   */
	  int N,          /* This is the number of inputs.  */
	  int L,          /* This is the number of levels.  */
	  const pqf *H,   /* Low-pass QF data structure.    */
	  const pqf *G)   /* High-pass QF data structure.   */
{
  int  M, i;

  M = N>>L;
  while( M < N )
    {
      acdpe ( work, 1,  data,  M, H );
      acdpo ( work, 1, data+M, M, G );
      M *= 2;
      for(i=0; i<M; i++) data[i] = work[i];
    }
  return;
}

/*********************************************************
 * dwtaintervals():
 *
 * Compute the endpoints of INTERVALs to be filled by the
 * aperiodic discrete wavelet transform down to a specified
 * level.
 *
 * Calling sequence and basic algorithm:
 *
 *   dwtaintervals( V, W, IN, L, H, G ):
 *      If L>0 then
 *         Let V.LEAST = cdaleast( IN, H )
 *         Let V.FINAL = cdafinal( IN, H )
 *         Let W.LEAST = cdaleast( IN, G )
 *         Let W.FINAL = cdafinal( IN, G )
 *         dwtaintervals( V+1, W+1, V, L-1, H, G )
 *      Else 
 *         Let W.LEAST = IN.LEAST
 *         Let W.FINAL = IN.FINAL  
 *     Return
 *
 *
 * Inputs:
 *	(interval *)V		This preallocated array must have
 *				  at least `L' INTERVAL elements.
 *
 *	(interval *)W		This preallocated array must have
 *				  at least `L+1' INTERVAL elements.
 *
 *	(interval *)in		This must point to a preallocated
 *				  and assigned INTERVAL
 *
 *	(int)L			This must be a nonnegative integer.
 *
 *	(const pqf *)H		These are the QF structs to be used
 *	(const pqf *)G		  for filter convolution-decimation.
 *
 * Outputs:
 *	(interval *)V		The `least' and `final' members of these
 *	(interval *)W		  arrays are filled by reference with
 *				  their correct endpoints for `dwta()'.
 *
 */
extern void
  dwtaintervals( 
		interval *V,	/* The `L' scaling subspaces.    */
		interval *W,    /* The `L+1' wavelet subspaces.  */
		interval *in,   /* This is the input interval.   */
		int L,          /* This is the number of levels. */
		const pqf *H,   /* Low-pass QF data structure.   */
		const pqf *G)   /* High-pass QF data structure.  */
{

  if( L>0 )
    {
      V->least = cdaleast( in, H );
      V->final = cdafinal( in, H );
      W->least = cdaleast( in, G );
      W->final = cdafinal( in, G );
      dwtaintervals( V+1, W+1, V, L-1, H, G );
    }
  else
    { 
      W->least = in->least;
      W->final = in->final;
    }
  return;
}

/*********************************************************
 * dwtaorigins():
 *
 * Compute the endpoints of INTERVALs to be filled by the
 * aperiodic discrete wavelet transform down to a specified
 * level.
 *
 * Calling sequence and basic algorithm:
 *
 *   dwtaorigins( V, W, SUMS, DIFS, L ):
 *      If L>0 then
 *         Let V.ORIGIN = SUMS + shifttoorigin( V )
 *         Let W.ORIGIN = DIFS + shifttoorigin( W )
 *         SUMS += shifttonextinterval( V )
 *         DIFS += shifttonextinterval( W )
 *         dwtaorigins( V+1, W+1, SUMS, DIFS, L-1 )
 *      Else
 *         Let W.ORIGIN = DIFS + shifttoorigin( V )
 *      Return
 *
 *
 * Inputs:
 *	(interval *)V		This preallocated array must have
 *				  at least `L' INTERVAL elements.
 *
 *	(interval *)W		This preallocated array must have
 *				  at least `L+1' INTERVAL elements.
 *
 *	(real *)sums		These must point to the beginnings of
 *	(real *)difs		 preallocated, sufficiently long arrays.
 *
 *	(int)L			This must be a nonnegative integer.
 *
 * Outputs:
 *	(interval *)V		The `origin' members of these arrays
 *	(interval *)W		  are set by reference to addresses
 *				  within `sums[]' and `difs[]'.
 *
 */
extern void
  dwtaorigins( 
	      interval *V,	/* The `L' scaling subspaces.       */
	      interval *W,	/* The `L+1' wavelet subspaces.     */
	      real *sums,	/* Concatenated scaling amplitudes. */
	      real *difs,	/* Concatenated wavelet amplitudes. */
	      int L)		/* This is the number of levels.    */
{
  if( L>0 )
    {
      V->origin = sums + shifttoorigin( V );
      W->origin = difs + shifttoorigin( W );
      sums += shifttonextinterval( V );
      difs += shifttonextinterval( W );
      dwtaorigins( V+1, W+1, sums, difs, L-1 );
    }
  else
    {
      W->origin = difs + shifttoorigin( V );
    }
  return;
}

/*********************************************************
 * dwta():
 *
 * [D]iscrete [W]avelet [T]ransform, [A]periodic.
 *
 * Calling sequence and basic algorithm:
 *
 *   dwta( V, W, IN, L, H, G ):
 *      If L > 0
 *         cdao( V.ORIGIN, 1, IN.ORIGIN, IN.LEAST, IN.FINAL, H )
 *         cdao( W.ORIGIN, 1, IN.ORIGIN, IN.LEAST, IN.FINAL, G )
 *         dwta( V+1, W+1, V, L-1, H, G )
 *      Else
 *         For N = IN.LEAST to IN.FINAL
 *            Let W.ORIGIN[N] = IN.ORIGIN[N]
 *      Return
 *
 *
 * Inputs:
 *	(interval *)V		This array must contain at least
 *				  `L' intervals.
 *
 *	(interval *)W		This array must contain at least
 *				  `L+1' intervals.
 *
 *	(interval *)in		This interval holds the input signal.
 *				  It is not changed by this function.
 *
 *	(int)L			This is the nonnegative number of levels.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(interval *)V		These intervals have the aperiodic scaling
 *				  amplitudes superposed into their arrays.
 *
 *	(interval *)W		These intervals have the aperiodic wavelet
 *				  amplitudes superposed into their arrays.
 *
 *
 * External functions called:
 *	cdao()
 */
extern void
  dwta( 
       interval *V,    /* Array of `L' scaling subspaces.   */
       interval *W,    /* Array of `L+1' wavelet subspaces. */
       interval *in,   /* Preallocated input interval.      */
       int L,          /* This is the number of levels.     */
       const pqf *H,   /* Low-pass QF data structure.       */
       const pqf *G)   /* High-pass QF data structure.      */
{
  if( L > 0 )
    {
      cdao( V->origin, 1, in->origin, in->least, in->final, H );
      cdao( W->origin, 1, in->origin, in->least, in->final, G );
      dwta( V+1, W+1, V, L-1, H, G );
    }
  else
    {
      int n;
      for( n = in->least; n <= in->final; n++ )
	 W->origin[n] = in->origin[n];
    }
  return;
}

/****************************************************************
 * dwtan():
 *
 * [D]iscrete [W]avelet [T]ransform, [A]periodic, [N]onrecursive.
 *
 * Calling sequence and basic algorithm:
 *
 *   dwtan( V, W, IN, L, H, G ):
 *      For K = 0 to L-1
 *         cdao( V[K].ORIGIN, 1, IN.ORIGIN, IN.LEAST, IN.FINAL, H )
 *         cdao( W[K].ORIGIN, 1, IN.ORIGIN, IN.LEAST, IN.FINAL, G )
 *         Let IN = V[K]
 *      For N = IN.LEAST to IN.FINAL
 *         Let W[L].ORIGIN[N] = IN.ORIGIN[N]
 *
 *
 * Inputs:
 *	(interval *)V		This array must contain at least
 *				  `L' intervals.
 *
 *	(interval *)W		This array must contain at least
 *				  `L+1' intervals.
 *
 *	(interval *)in		This interval holds the input signal.
 *				  It is not changed by this function.
 *
 *	(int)L			This is the nonnegative number of levels.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(interval *)V		These intervals have the aperiodic scaling
 *				  amplitudes superposed into their arrays.
 *
 *	(interval *)W		These intervals have the aperiodic wavelet
 *				  amplitudes superposed into their arrays.
 *
 * External functions called:
 *	cdao()
 */
extern void
  dwtan( 
	interval *V,		/* Array of `L' scaling subspaces.   */
	interval *W,		/* Array of `L+1' wavelet subspaces. */
	interval *in,		/* Preallocated input interval.      */
	int L,			/* This is the number of levels.     */
	const pqf *H,		/* Low-pass QF data structure.       */
	const pqf *G)		/* High-pass QF data structure.      */
{
  int k, n;

  for( k=0; k<L; k++ )
    {
      cdao( V[k].origin, 1, in->origin, in->least, in->final, H );
      cdao( W[k].origin, 1, in->origin, in->least, in->final, G );
      in = V+k;
    }
  for( n=in->least; n<=in->final; n++ )
    W[L].origin[n] = in->origin[n];
  return;
}

/****************************************************************
 * dwtacomplete():
 *
 * Complete [D]iscrete [W]avelet [T]ransform, [A]periodic
 * from an array.
 *
 * Calling sequence and basic algorithm:
 *
 *  dwtacomplete( DATA, LENGTH, L, H, G ):
 *     Allocate L empty INTERVALs at V
 *     Allocate L+1 empty INTERVALs at W
 *     Let IN = makeinterval( DATA, 0, LENGTH-1 )
 *     dwtaintervals( V, W, IN, L, H, G )
 *     Let NS = intervalstotal( V, L )
 *     Allocate an array of NS zeroes at SUMS
 *     Let ND = intervalstotal( W, L+1 )
 *     Allocate an array of NS zeroes at DIFS
 *     dwtaorigins( V, W, SUMS+NS, DIFS+ND, L )
 *     dwta( V, W, IN, L, H, G )
 *     Deallocate SUMS[]
 *     Deallocate IN, but leave DATA[] alone
 *     Deallocate V[]
 *     Return W
 *
 *
 * Inputs:
 *	(real *)data		This array contains `length' input values.
 *
 *	(int)length		This must be a positive integer
 *
 *	(int)L			This is the nonnegative number of levels.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(interval *)dwtacomplete	These intervals have the aperiodic
 *				  scaling amplitudes superposed into their
 *				  data arrays.
 *
 * External functions called:
 *	cdao(), dwtaintervals(), dwtaorigins(),
 *	intervalstotal(),
 *	calloc(), free()
 */
extern interval *
  dwtacomplete( 
	       real *data,	/* Array of `length' input values.     */
	       int length,	/* Positive length of the input array. */
	       int L,		/* This is the number of levels.       */
	       const pqf *H,	/* Low-pass QF data structure.         */
	       const pqf *G)	/* High-pass QF data structure.        */
{
  interval *V, *W, *in;
  real *sums, *difs;
  int ns, nd;

  V = (interval *)calloc( L, sizeof(interval)); assert(V);
  W = (interval *)calloc( L+1, sizeof(interval)); assert(W);
  in = makeinterval( data, 0, length-1 );
  dwtaintervals( V, W, in, L, H, G );
  ns = intervalstotal( V, L );
  sums = (real *)calloc(ns, sizeof(real)); assert(sums);
  nd = intervalstotal( W, L+1 );
  difs = (real *)calloc(nd, sizeof(real)); assert(difs);
  dwtaorigins( V, W, sums+ns, difs+nd, L ); /* cat from end */
  dwta( V, W, in, L, H, G );
  free(sums);
  in->origin = 0; free(in);
  freeinterval(V);
  return(W);
}

/****************************************************************
 * dwtalocal():
 *
 * [D]iscrete [W]avelet [T]ransform, [A]periodic from
 * an INTERVAL, with local allocation of subspace INTERVALs.
 *
 * Calling sequence and basic algorithm:
 *
 *  dwtalocal( V, W, IN, L, H, G ):
 *     If L > 0
 *        Let LEAST = cdaleast( IN, H )
 *        Let FINAL = cdafinal( IN, H )
 *        Let V = makeinterval( NULL, LEAST, FINAL )
 *        cdai( V.ORIGIN, 1, IN.ORIGIN, IN.LEAST, IN.FINAL, H )
 *        Let LEAST = cdaleast( IN, G )
 *        Let FINAL = cdafinal( IN, G )
 *        Let W = makeinterval( NULL, LEAST, FINAL )
 *        cdai( W.ORIGIN, 1, IN.ORIGIN, IN.LEAST, IN.FINAL, G )
 *        dwtalocal( V+1, W+1, V, L-1, H, G )
 *     Else
 *        Let W = makeinterval( NULL, IN.LEAST, IN.FINAL )
 *        For N = IN.LEAST to IN.FINAL
 *           Let W.ORIGIN[N] = IN.ORIGIN[N]
 *     Return
 *
 *
 * Inputs:
 *	(interval *)V		This array must contain at least
 *				  `L' intervals.
 *
 *	(interval *)W		This array must contain at least
 *				  `L+1' intervals.
 *
 *	(interval *)in		This interval holds the input signal.
 *				  It is not changed by this function.
 *
 *	(int)L			This is the nonnegative number of levels.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(interval *)V		These intervals have the aperiodic wavelet
 *	(interval *)W		  and scaling amplitudes written into their
 *				  newly allocated data arrays.
 *
 * External functions called:
 *	cdai(), makeinterval(), cdaleast(), cdafinal()
 */
extern void
  dwtalocal( 
	    interval *V,	/* Array of `L' empty INTERVALs.     */
	    interval *W,	/* Array of `L+1' empty INTERVALs.   */
	    interval *in,	/* Input INTERVAL.                   */
	    int L,		/* This is the number of levels.     */
	    const pqf *H,	/* Low-pass QF data structure.       */
	    const pqf *G)	/* High-pass QF data structure.      */
{
  int least, final, n;
  if( L>0 )
    {
      least = cdaleast( in, H );
      final = cdafinal( in, H );
      V = makeinterval( 0, least, final );
      cdai( V->origin, 1, in->origin, in->least, in->final, H );
      least = cdaleast( in, G );
      final = cdafinal( in, G );
      W = makeinterval( 0, least, final );
      cdai( W->origin, 1, in->origin, in->least, in->final, G );
      dwtalocal( V+1, W+1, V, L-1, H, G );
    }
  else
    {
      W = makeinterval( 0, in->least, in->final );
      for( n=in->least; n<=in->final; n++ )
	W->origin[n] = in->origin[n];
    }
  return;
}

/*********************************************************
 * idwtaintervals():
 *
 * Compute the endpoints of INTERVALs to be filled by the
 * aperiodic inverse discrete wavelet transform down to a
 * specified level.
 *
 * Calling sequence and basic algorithm:
 *
 *   idwtaintervals( OUT, V, W, L, H, G ):
 *      If L>0 then
 *         idwtaintervals( V, V+1, W+1, L-1, H, G )
 *         Let OUT.LEAST = min( acdaleast(V,H), acdaleast(W,G) )
 *         Let OUT.FINAL = max( acdafinal(V,H), acdafinal(W,G) )
 *      Else
 *         Let OUT.LEAST = W.LEAST
 *         Let OUT.FINAL = W.FINAL
 *      Return
 *
 *
 * Inputs:
 *	(interval *)out		This must point to a preallocated
 *				  INTERVAL for the output.
 *
 *	(interval *)V		This preallocated array must have
 *				  at least `L' INTERVAL elements.
 *
 *	(interval *)W		This preallocated array must have
 *				  at least `L+1' INTERVAL elements.
 *
 *	(int)L			This must be a nonnegative integer.
 *
 *	(const pqf *)H		These are the QF structs to be used
 *	(const pqf *)G		  for adjoint convolution-decimation.
 *
 * Outputs:
 *	(interval *)V		The `least' and `final' members of these
 *	(interval *)out		  arrays are filled by reference with
 *				  their correct endpoints for `idwta()'.
 *
 */
extern void
  idwtaintervals( 
		 interval *out, /* This is the output INTERVAL.  */
		 interval *V,	/* The `L' scaling subspaces.    */
		 interval *W,   /* The `L+1' wavelet subspaces.  */
		 int L,         /* This is the number of levels. */
		 const pqf *H,  /* Low-pass QF data structure.   */
		 const pqf *G)  /* High-pass QF data structure.  */
{

  if( L>0 )
    {
      idwtaintervals( V, V+1, W+1, L-1, H, G );
      out->least = min( acdaleast(V,H), acdaleast(W,G) );
      out->final = max( acdafinal(V,H), acdafinal(W,G) );
    }
  else
    {
      out->least = W->least;
      out->final = W->final;
    }
  return;
}

/*********************************************************
 * idwtaorigins():
 *
 * Compute the endpoints of INTERVALs to be filled by the
 * aperiodic adjoint discrete wavelet transform down to a
 * specified level.
 *
 * Calling sequence and basic algorithm:
 *
 *   idwtaorigins( V, SUMS, L ):
 *      If L>0 then
 *         Let V.ORIGIN = SUMS + shifttoorigin( V )
 *         SUMS += shifttonextinterval( V )
 *         idwtaorigins( V+1, SUMS, L-1 )
 *      Return
 *
 *
 * Inputs:
 *	(interval *)V		This preallocated array must have
 *				  at least `L' INTERVAL elements.
 *
 *	(real *)sums		This must point to the beginnings of a
 *				 preallocated, sufficiently long array.
 *
 *	(int)L			This must be a nonnegative integer.
 *
 * Outputs:
 *	(interval *)V		The `origin' members of this array
 *				  are set by reference to addresses
 *				  within `sums[]'.
 *
 */
extern void
  idwtaorigins( 
	       interval *V,	/* The `L' scaling subspaces.         */
	       real *sums,	/* This is the array of amplitudes.   */
	       int L)		/* This is the number of levels.      */
{
  if( L>0 )
    {
      V->origin = sums + shifttoorigin( V );
      sums += shifttonextinterval( V );
      idwtaorigins( V, sums, L-1 );
    }
  return;
}

/*********************************************************
* idwta():
 *
 * [I]nverse [D]iscrete [W]avelet [T]ransform, [A]periodic.
 *
 * Calling sequence and basic algorithm:
 *
 *   idwta( OUT, V, W, L, H, G ):
 *      If L>0 then
 *         idwta( V, V+1, W+1, L-1, H, G )
 *         acdao( OUT.ORIGIN, 1, V.ORIGIN, V.LEAST, V.FINAL, H )
 *         acdao( OUT.ORIGIN, 1, W.ORIGIN, W.LEAST, W.FINAL, G )
 *      Else
 *         For N = W.LEAST to W.FINAL
 *            Let OUT.ORIGIN[N] = W.ORIGIN[N]
 *      Return
 *
 *
 * Inputs:
 *	(interval *)out		This interval must be preallocated and
 *				  assigned an array at its `origin' member.
 *
 *	(interval *)V		This array of `L' intervals must have all 0s
 *				  in its members' data arrays.  Those arrays
 *				  will be filled with scaling amplitudes.
 *
 *	(interval *)W		This array of `L+1' intervals must have the
 *				  wavelet amplitudes in its members' data
 *				  arrays.
 *
 *	(int)L			This is the nonnegative number of levels.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(interval *)out		This interval has the reconstructed sequence
 *				  superposed onto its data array by reference.
 *
 * External functions called:
 *	acdao()
 *
 */
extern void
  idwta( 
	interval *out, /* Preallocated, zeroed output interval.     */
	interval *V,   /* Working intervals for scaling amplitudes. */
	interval *W,   /* Input intervals of wavelet amplitudes.    */
	int L,         /* This is the number of levels. */
	const pqf *H,  /* Low-pass QF data structure.   */
	const pqf *G)  /* High-pass QF data structure.  */
{
  if( L>0 )
    {
      idwta( V, V+1, W+1, L-1, H, G );
      acdao( out->origin, 1, V->origin, V->least, V->final, H );
      acdao( out->origin, 1, W->origin, W->least, W->final, H );
    }
  else
    {
      int n;
      for( n=W->least; n<= W->final; n++)
	out->origin[n] = W->origin[n];
    }
  return;
}

/*********************************************************
* idwtacomplete():
 *
 * Complete [I]nverse [D]iscrete [W]avelet [T]ransform, 
 * [A]periodic, from an input array of INTERVALs.
 *
 * Calling sequence and basic algorithm:
 *
 *   idwtacomplete( W, L, H, G ):
 *      Allocate L empty INTERVALs at V
 *      Allocate 1 empty INTERVAL at OUT
 *      idwtaintervals( OUT, V, W, L, H, G )
 *      Let OUT = enlargeinterval( OUT, OUT.LEAST, OUT.FINAL )
 *      Let NS = intervalstotal( V, L )
 *      Allocate an array of NS zeroes at SUMS
 *      idwtaorigins( V, SUMS+NS, L )
 *      idwta( OUT, V, W, L, H, G )
 *      Deallocate SUMS[]
 *      Deallocate V[]
 *      Return OUT
 *
 *
 * Inputs:
 *	(interval *)W		This array of `L+1' intervals must have the
 *				  wavelet amplitudes in its members' data
 *				  arrays.
 *
 *	(int)L			This is the nonnegative number of levels.
 *
 *	(const pqf *)H		These are QF structs used for low-pass and
 *	(const pqf *)G		  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)idwtacomplete	The return value is a pointer to a newly  
 *				  allocated interval containing the 
 *				  reconstructed signal.
 *
 * External functions called:
 *	acdao()
 *	enlargeinterval(), intervalstotal()
 *	idwtaintervals(), idwtaorigins()
 *	calloc(), free()
 *
 */
extern interval *
  idwtacomplete( 
		interval *W,	/* Input wavelet amplitudes.     */
		int L,		/* This is the number of levels. */
		const pqf *H,	/* Low-pass QF data structure.   */
		const pqf *G)	/* High-pass QF data structure.  */
{
  interval *V, *out;
  int ns;
  real *sums;

  V   = (interval *)calloc( L, sizeof(interval));  assert(V);
  out = (interval *)calloc( 1, sizeof(interval));  assert(out);
  idwtaintervals( out, V, W, L, H, G );
  out = enlargeinterval( out, out->least, out->final );
  ns = intervalstotal( V, L );
  sums = (real *)calloc(ns, sizeof(real)); assert(sums);
  idwtaorigins( V, sums+ns, L ); /* cat from end */
  idwta( out, V, W, L, H, G );
  free(sums);
  free(V);
  return(out);
}

