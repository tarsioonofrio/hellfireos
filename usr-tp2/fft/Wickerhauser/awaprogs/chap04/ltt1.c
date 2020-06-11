/* 
 * Basic local trigonometric transform functions.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include "ltt1.h"

#include <stdlib.h>		/* for calloc(), malloc() */
#include <assert.h>		/* We abort if malloc() returns NULL. */
#include "real.h"
#include "interval.h"		/* Rising cutoff arrays are INTERVALs. */
#include "fold.h"		/* Declare `fips()' and `fipc()' functions. */
#include "dtts.h"		/* Declare various DCT and DST functions. */


/***********************************************************************
 * lct()
 *
 *  [L]ocal [C]osine [T]ransform (in place).
 *
 *  Calling sequence and basic algorithm:
 *
 *    lct( SIG, N, RISE ):
 *       fipc( SIG, SIG, RISE )
 *       fipc( SIG+N, SIG+N, RISE )
 *       dctiv( SIG, N )
 *
 *  Inputs:
 *	(real *)sig		This is the input and output array.
 *
 *	(int)n			This is a nonnegative power of 2.
 *
 *	(const interval *)rise	This is a sampled rising cutoff function.
 *
 * Outputs:
 *	(real *)sig		Two applications of `fipc()' and the
 *				 DCT-IV transform are performed in place
 *				 on `sig[]'.
 *
 * External functions called:
 *	fipc()		Declared in "fold.h"
 *	dctiv()		Declared in "dtts.h"
 */
extern void
  lct(
      real *sig,		/* Input and output array.           */
      int n,			/* Length of the input/output array. */
      const interval *rise)	/* Sampled rising cutoff function.   */
{
  int log2n;

  /* Compute the logarithm base 2 of the signal length `n': */
  assert(n>0);
  log2n = 0;
  while( (1<<log2n) < n ) ++log2n;
  assert((1<<log2n)==n);

  /* Fold at the left and right endpoints: */
  fipc( sig, sig, rise );
  fipc( sig+n, sig+n, rise );

  /* Apply DCT-IV: */
  dct_iv( sig, log2n );

  return;
}

/***********************************************************************
 * lst()
 *
 *  [L]ocal [S]ine [T]ransform (in place).
 *
 *  Calling sequence and basic algorithm:
 *
 *    lst( SIG, N, RISE ):
 *       fips( SIG, SIG, RISE )
 *       fips( SIG+N, SIG+N, RISE )
 *       dstiv( SIG, N )
 *
 *  Inputs:
 *	(real *)sig		This is the input and output array.
 *
 *	(int)n			This is a nonnegative power of 2.
 *
 *	(const interval *)rise	This is a sampled rising cutoff function.
 *
 * Outputs:
 *	(real *)sig		Two applications of `fipc()' and the
 *				 DCT-IV transform are performed in place
 *				 on `sig[]'.
 *
 * External functions called:
 *	fips()		Declared in "fold.h"
 *	dstiv()		Declared in "dtts.h"
 */
extern void
  lst(
      real *sig,		/* Input and output array.           */
      int n,			/* Length of the input/output array. */
      const interval *rise)	/* Sampled rising cutoff function.   */
{
  int log2n;

  /* Compute the logarithm base 2 of the signal length `n': */
  assert(n>0);
  log2n = 0;
  while( (1<<log2n) < n ) ++log2n;
  assert((1<<log2n)==n);

  /* Fold at the left and right endpoints: */
  fips( sig, sig, rise );
  fips( sig+n, sig+n, rise );

  /* Apply DST-IV: */
  dst_iv( sig, log2n );

  return;
}

/***********************************************************************
 * ilct()
 *
 *  [I]nverse [L]ocal [C]osine [T]ransform (in place).
 *
 *  Calling sequence and basic algorithm:
 *
 *    ilct( SIG, N, RISE ):
 *       dctiv( SIG, N )
 *       uipc( SIG, SIG, RISE )
 *       uipc( SIG+N, SIG+N, RISE )
 *
 *  Inputs:
 *	(real *)sig		This is the input and output array.
 *
 *	(int)n			This is a nonnegative power of 2.
 *
 *	(const interval *)rise	This is a sampled rising cutoff function.
 *
 * Outputs:
 *	(real *)sig		The DCT-IV transform and two applications
 *				of `uipc()' are performed in place on `sig[]'.
 *
 * External functions called:
 *	uipc()		Declared in "fold.h"
 *	dctiv()		Declared in "dtts.h"
 */
extern void
  ilct(
       real *sig,		/* Input and output array.           */
       int n,			/* Length of the input/output array. */
       const interval *rise)	/* Sampled rising cutoff function.   */
{
  int log2n;

  /* Compute the logarithm base 2 of the signal length `n': */
  assert(n>0);
  log2n = 0;
  while( (1<<log2n) < n ) ++log2n;
  assert((1<<log2n)==n);

  /* Apply DCT-IV, which is its own inverse: */
  dct_iv( sig, log2n );

  /* Unfold at the left and right endpoints: */
  uipc( sig, sig, rise );
  uipc( sig+n, sig+n, rise );

  return;
}

/***********************************************************************
 * ilst()
 *
 *  [I]nverse [L]ocal [S]ine [T]ransform (in place).
 *
 *  Calling sequence and basic algorithm:
 *
 *    ilst( SIG, N, RISE ):
 *       dstiv( SIG, N )
 *       uips( SIG, SIG, RISE )
 *       uips( SIG+N, SIG+N, RISE )
 *
 *  Inputs:
 *	(real *)sig		This is the input and output array.
 *
 *	(int)n			This is a nonnegative power of 2.
 *
 *	(const interval *)rise	This is a sampled rising cutoff function.
 *
 * Outputs:
 *	(real *)sig		The DST-IV transform and two applications
 *				of `uips()' are performed in place on `sig[]'.
 *
 * External functions called:
 *	uips()		Declared in "fold.h"
 *	dstiv()		Declared in "dtts.h"
 */
extern void
  ilst(
       real *sig,		/* Input and output array.           */
       int n,			/* Length of the input/output array. */
       const interval *rise)	/* Sampled rising cutoff function.   */
{
  int log2n;

  /* Compute the logarithm base 2 of the signal length `n': */
  assert(n>0);
  log2n = 0;
  while( (1<<log2n) < n ) ++log2n;
  assert((1<<log2n)==n);

  /* Apply DST-IV, which is own inverse: */
  dst_iv( sig, log2n );

  /* Unfold at the left and right endpoints: */
  uips( sig, sig, rise );
  uips( sig+n, sig+n, rise );

  return;
}

/***********************************************************************
 * lpdc()
 *
 *  [L]ocal [P]eriodization [D]isjoint, [C]osine polarity.
 *
 *  Calling sequence and basic algorithm:
 *
 *    lpdc( OUT, IN, N, RISE ):
 *       fdcp(  OUT,  1,  IN,   IN,  N/2, RISE )
 *       fdcn( OUT+N, 1, IN+N, IN+N, N/2, RISE )
 *       uipc( OUT+N, OUT, RISE )
 *
 *  Inputs:
 *	(real *)out		This is the output array.
 *
 *	(const real *)out	This is the unchanged input array.
 *
 *	(int)n			This is greater than `2*rise->final'.
 *
 *	(const interval *)rise	This is a sampled rising cutoff function.
 *
 * Outputs:
 *	(real *)out		The periodized signal is written to this array.
 *
 *
 * Assumptions:
 *	1. rise != NULL
 *	2. n>2*rise->final
 *
 * External functions called:
 *	fdcn(), fdcp(), uipc(), assert()
 */
extern void
  lpdc(
       real *out,		/* Preallocated output array.      */
       const real *in,		/* Preallocated input array.       */
       int n,			/* Length to which we periodize.   */
       const interval *rise)	/* Sampled rising cutoff function. */
{
  assert(rise);
  assert(n>2*rise->final);

  /* Fold at the left and right endpoints: */
  fdcp(  out,  1,  in,   in,  n/2, rise );
  fdcn( out+n, 1, in+n, in+n, n/2, rise );

  /* Unfold the right endpoint against the left: */
  uipc( out+n, out, rise );

  return;
}

/***********************************************************************
 * lpic()
 *
 *  [L]ocal [P]eriodization [I]n-place, [C]osine polarity.
 *
 *  Calling sequence and basic algorithm:
 *
 *    lpic( SIG, N, RISE ):
 *       fipc( SIG, SIG, RISE )
 *       fipc( SIG+N, SIG+N, RISE )
 *       uipc( SIG+N, SIG, RISE )
 *
 *  Inputs:
 *	(real *)sig		This is the input and output array.
 *
 *	(int)n			This is greater than `2*rise->final'.
 *
 *	(const interval *)rise	This is a sampled rising cutoff function.
 *
 * Outputs:
 *	(real *)sig		Two applications of `fipc()' and one
 *				 of `uipc()' are performed in place on `sig[]'.
 *
 * External functions called:
 *	fipc(),uipc()		Declared in "fold.h"
 */
extern void
  lpic(
       real *sig,		/* Input and output array.           */
       int n,			/* Length of the input/output array. */
       const interval *rise)	/* Sampled rising cutoff function.   */
{
  assert(rise);
  assert(n>2*rise->final);
  
  /* Fold at the left and right endpoints: */
  fipc( sig, sig, rise );
  fipc( sig+n, sig+n, rise );

  /* Unfold the right endpoint against the left: */
  uipc( sig+n, sig, rise );

  return;
}

/***********************************************************************
 * lpds()
 *
 *  [L]ocal [P]eriodization [D]isjoint, [S]ine polarity.
 *
 *  Calling sequence and basic algorithm:
 *
 *    lpds( OUT, IN, N, RISE ):
 *       fdsp(  OUT,  1,  IN,   IN,  N/2, RISE )
 *       fdsn( OUT+N, 1, IN+N, IN+N, N/2, RISE )
 *       uips( OUT+N, OUT, RISE )
 *
 *  Inputs:
 *	(real *)out		This is the output array.
 *
 *	(const real *)out	This is the unchanged input array.
 *
 *	(int)n			This is greater than `2*rise->final'.
 *
 *	(const interval *)rise	This is a sampled rising cutoff function.
 *
 * Outputs:
 *	(real *)out		The periodized signal is written to this array.
 *
 *
 * Assumptions:
 *	1. rise != NULL
 *	2. n>2*rise->final
 *
 * External functions called:
 *	fdsn(), fdsp(), uips(), assert()
 */
extern void
  lpds(
       real *out,		/* Preallocated output array.      */
       const real *in,		/* Preallocated input array.       */
       int n,			/* Length to which we periodize.   */
       const interval *rise)	/* Sampled rising cutoff function. */
{
  assert(rise);
  assert(n>2*rise->final);

  /* Fold at the left and right endpoints: */
  fdsp(  out,  1,  in,   in,  n/2, rise );
  fdsn( out+n, 1, in+n, in+n, n/2, rise );

  /* Unfold the right endpoint against the left: */
  uips( out+n, out, rise );

  return;
}

/***********************************************************************
 * lpis()
 *
 *  [L]ocal [P]eriodization [I]n-place, [S]ine polarity.
 *
 *  Calling sequence and basic algorithm:
 *
 *    lpis( SIG, N, RISE ):
 *       fips( SIG, SIG, RISE )
 *       fips( SIG+N, SIG+N, RISE )
 *       uips( SIG+N, SIG, RISE )
 *
 *  Inputs:
 *	(real *)sig		This is the input and output array.
 *
 *	(int)n			This is greater than `2*rise->final'.
 *
 *	(const interval *)rise	This is a sampled rising cutoff function.
 *
 * Outputs:
 *	(real *)sig		Two applications of `fips()' and one
 *				 of `uips()' are performed in place on `sig[]'.
 *
 * External functions called:
 *	fips(),uips()		Declared in "fold.h"
 */
extern void
  lpis(
       real *sig,		/* Input and output array.           */
       int n,			/* Length of the input/output array. */
       const interval *rise)	/* Sampled rising cutoff function.   */
{
  assert(rise);
  assert(n>2*rise->final);

  /* Fold at the left and right endpoints: */
  fips( sig, sig, rise );
  fips( sig+n, sig+n, rise );

  /* Unfold the right endpoint against the left: */
  uips( sig+n, sig, rise );

  return;
}

/***********************************************************************
 * ilpic()
 *
 *  [I]nverse [L]ocal [P]eriodization [I]n-place, [C]osine polarity.
 *
 *  Calling sequence and basic algorithm:
 *
 *    ilpic( SIG, N, RISE ):
 *       fipc( SIG+N, SIG, RISE )
 *       uipc( SIG, SIG, RISE )
 *       uipc( SIG+N, SIG+N, RISE )
 *
 *  Inputs:
 *	(real *)sig		This is the input and output array.
 *
 *	(int)n			This is greater than `2*rise->final'.
 *
 *	(const interval *)rise	This is a sampled rising cutoff function.
 *
 * Outputs:
 *	(real *)sig		One applications of `fipc()' and two
 *				 of `uipc()' are performed in place on `sig[]'.
 *
 * External functions called:
 *	fipc(),uipc()		Declared in "fold.h"
 */
extern void
  ilpic(
	real *sig,		/* Input and output array.           */
	int n,			/* Length of the input/output array. */
	const interval *rise)	/* Sampled rising cutoff function.   */
{
  assert(rise);
  assert(n>2*rise->final);

  /* Fold the right endpoint against the left: */
  fipc( sig+n, sig, rise );

  /* Unfold at the left and right endpoints: */
  uipc( sig, sig, rise );
  uipc( sig+n, sig+n, rise );

  return;
}

/***********************************************************************
 * ilpis()
 *
 *  [I]nverse [L]ocal [P]eriodization [I]n-place, [S]ine polarity.
 *
 *  Calling sequence and basic algorithm:
 *
 *    ilpis( SIG, N, RISE ):
 *       fips( SIG+N, SIG, RISE )
 *       uips( SIG, SIG, RISE )
 *       uips( SIG+N, SIG+N, RISE )
 *
 *  Inputs:
 *	(real *)sig		This is the input and output array.
 *
 *	(int)n			This is greater than `2*rise->final'.
 *
 *	(const interval *)rise	This is a sampled rising cutoff function.
 *
 * Outputs:
 *	(real *)sig		Two applications of `uips()' and one
 *				 of `fips()' are performed in place on `sig[]'.
 *
 * External functions called:
 *	fips(),uips()		Declared in "fold.h"
 */
extern void
  ilpis(
       real *sig,		/* Input and output array.           */
       int n,			/* Length of the input/output array. */
       const interval *rise)	/* Sampled rising cutoff function.   */
{
  assert(rise);
  assert(n>2*rise->final);

  /* Fold the right endpoint against the left: */
  fips( sig+n, sig, rise );

  /* Unfold at the left and right endpoints: */
  uips( sig, sig, rise );
  uips( sig+n, sig+n, rise );

  return;
}

/***********************************************************************
 * lpica()
 *
 *  [L]ocal [P]eriodization [I]n-place, [C]osine polarity, [A]djacent
 *  intervals.
 *
 *  Calling sequence and basic algorithm:
 *
 *    lpica( SIG, LENGTHS, NUM, RISE ):
 *       fipc( SIG, SIG, RISE )
 *       For I = 0 to NUM-1
 *          fipc( SIG+LENGTHS[I], SIG+LENGTHS[I], RISE )
 *          uipc( SIG+LENGTHS[I], SIG, RISE )
 *          SIG += LENGTHS[I]
 *
 *  Inputs:
 *	(real *)sig		This is the input and output array.
 *
 *	(int *)lengths		This is an array of `num' integers,
 *				  each greater than `2*rise->final+1'.
 *
 *	(int)num		This is the number of intervals.
 *
 *	(const interval *)rise	This is a sampled rising cutoff function.
 *
 * Outputs:
 *	(real *)sig		Two applications of `fipc()' and one
 *				 of `uipc()' are performed in place on 
 *				 each subinterval of `sig[]'.
 *
 * External functions called:
 *	fipc(),uipc()		Declared in "fold.h"
 */
extern void
  lpica(
	real *sig,		/* Input and output array.           */
	int *lengths,		/* Lengths of the subintervals.      */
	int  num,		/* Number of subintervals.           */
	const interval *rise)	/* Sampled rising cutoff function.   */
{
  int i;

  assert(rise);

  /* Fold at the leftmost endpoint: */
  fipc( sig, sig, rise );

  /* Loop over the intervals, folding at the right endpoint then unfolding
   * the left endpoint against the right endpoint: */
  for( i=0; i<=num; i++)
    {
      assert(lengths[i]>2*rise->final);

      fipc( sig+lengths[i], sig+lengths[i], rise );
      uipc( sig+lengths[i], sig, rise );
      sig += lengths[i];
    }
  return;
}

/***********************************************************************
 * lpisa()
 *
 *  [L]ocal [P]eriodization [I]n-place, [S]ine polarity, [A]djacent
 *  intervals.
 *
 *  Calling sequence and basic algorithm:
 *
 *    lpisa( SIG, LENGTHS, NUM, RISE ):
 *       fips( SIG, SIG, RISE )
 *       For I = 0 to NUM-1
 *          fips( SIG+LENGTHS[I], SIG+LENGTHS[I], RISE )
 *          uips( SIG+LENGTHS[I], SIG, RISE )
 *          SIG += LENGTHS[I]
 *
 *  Inputs:
 *	(real *)sig		This is the input and output array.
 *
 *	(int *)lengths		This is an array of `num' integers,
 *				  each greater than `2*rise->final+1'.
 *
 *	(int)num		This is the number of intervals.
 *
 *	(const interval *)rise	This is a sampled rising cutoff function.
 *
 * Outputs:
 *	(real *)sig		Two applications of `fips()' and one
 *				 of `uips()' are performed in place on 
 *				 each subinterval of `sig[]'.
 *
 * External functions called:
 *	fips(),uips()		Declared in "fold.h"
 */
extern void
  lpisa(
	real *sig,		/* Input and output array.           */
	int *lengths,		/* Lengths of the subintervals.      */
	int  num,		/* Number of subintervals.           */
	const interval *rise)	/* Sampled rising cutoff function.   */
{
  int i;

  assert(rise);

  /* Fold at the leftmost endpoint: */
  fips( sig, sig, rise );

  /* Loop over the intervals, folding at the right endpoint then unfolding
   * the left endpoint against the right endpoint: */
  for( i=0; i<=num; i++)
    {
      assert(lengths[i]>2*rise->final);

      fips( sig+lengths[i], sig+lengths[i], rise );
      uips( sig+lengths[i], sig, rise );
      sig += lengths[i];
    }
  return;
}

/***********************************************************************
 * ilpica()
 *
 *  [I]nverse [L]ocal [P]eriodization [I]n-place, [C]osine polarity,
 *  [A]djacent intervals.
 *
 *  Calling sequence and basic algorithm:
 *
 *     ilpica( SIG, LENGTHS, NUM, RISE ):
 *        For I = 0 to NUM-1
 *           fipc( SIG+LENGTHS[I], SIG, RISE )
 *           uipc( SIG, SIG, RISE )
 *           SIG += LENGTHS[I]
 *        uipc( SIG, SIG, RISE )
 *
 *  Inputs:
 *	(real *)sig		This is the input and output array.
 *
 *	(int *)lengths		This is an array of `num' integers,
 *				  each greater than `2*rise->final+1'.
 *
 *	(int)num		This is the number of intervals.
 *
 *	(const interval *)rise	This is a sampled rising cutoff function.
 *
 * Outputs:
 *	(real *)sig		Two applications of `uipc()' and one
 *				 of `fipc()' are performed in place on 
 *				 each subinterval of `sig[]'.
 *
 * External functions called:
 *	fipc(),uipc()		Declared in "fold.h"
 */
extern void
  ilpica(
	 real *sig,		/* Input and output array.           */
	 int *lengths,		/* Lengths of the subintervals.      */
	 int  num,		/* Number of subintervals.           */
	 const interval *rise)	/* Sampled rising cutoff function.   */
{
  int i;

  assert(rise);

  /* Loop over the intervals, folding the left endpoint against the
   * right endpoint and then unfolding at the left endpoint: */
  for( i=0; i<=num; i++)
    {
      assert(lengths[i]>2*rise->final);

      fipc( sig+lengths[i], sig, rise );
      uipc( sig, sig, rise );
      sig += lengths[i];
    }

  /* Unfold at the rightmost endpoint: */
  uipc( sig, sig, rise );

  return;
}

/***********************************************************************
 * ilpisa()
 *
 *  [I]nverse [L]ocal [P]eriodization [I]n-place, [S]ine polarity,
 *  [A]djacent intervals.
 *
 *  Calling sequence and basic algorithm:
 *
 *     ilpisa( SIG, LENGTHS, NUM, RISE ):
 *        For I = 0 to NUM-1
 *           fips( SIG+LENGTHS[I], SIG, RISE )
 *           uips( SIG, SIG, RISE )
 *           SIG += LENGTHS[I]
 *        uips( SIG, SIG, RISE )
 *
 *  Inputs:
 *	(real *)sig		This is the input and output array.
 *
 *	(int *)lengths		This is an array of `num' integers,
 *				  each greater than `2*rise->final+1'.
 *
 *	(int)num		This is the number of intervals.
 *
 *	(const interval *)rise	This is a sampled rising cutoff function.
 *
 * Outputs:
 *	(real *)sig		Two applications of `uips()' and one
 *				 of `fips()' are performed in place on 
 *				 each subinterval of `sig[]'.
 *
 * External functions called:
 *	fips(),uips()		Declared in "fold.h"
 */
extern void
  ilpisa(
	 real *sig,		/* Input and output array.           */
	 int *lengths,		/* Lengths of the subintervals.      */
	 int  num,		/* Number of subintervals.           */
	 const interval *rise)	/* Sampled rising cutoff function.   */
{
  int i;

  assert(rise);

  /* Loop over the intervals, folding the left endpoint against the
   * right endpoint and then unfolding at the left endpoint: */
  for( i=0; i<=num; i++)
    {
      assert(lengths[i]>2*rise->final);

      fips( sig+lengths[i], sig, rise );
      uips( sig, sig, rise );
      sig += lengths[i];
    }

  /* Unfold at the rightmost endpoint: */
  uips( sig, sig, rise );

  return;
}
