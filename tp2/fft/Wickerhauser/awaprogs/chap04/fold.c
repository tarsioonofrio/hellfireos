/*
 * 
 * These functions perform ``folding'' and ``unfolding,'' the
 * preprocessing step used before DCT-IV and DST-IV for the local
 * trigonometric transforms and their inverses.
 * 
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <assert.h>
#include "real.h"
#include "interval.h"
#include "rcf.h"
#include "fold.h"

/*********************************************************************
 * fdcn()
 *
 * [F]old a given array into a [D]isjoint output array, with [C]osine
 * polarity (odd on the left and even on the right), into the [N]egative
 * part of the output array.  Elements `outneg[-n*step],...,outneg[-1*step]'
 * of the output array are assigned.
 *
 * Calling sequence and basic algorithm:
 *
 *    fdcn( ONEG, STEP, INEG, IPOS, N, RISE ):
 *      For K = -N to RISE.LEAST-1
 *        Let ONEG[K*STEP] = INEG[K]
 *      For K = RISE.LEAST to -1
 *        Let ONEG[K*STEP] = RISE.ORIGIN[-1-K] * INEG[K]
 *                            - RISE.ORIGIN[K] * IPOS[-1-K]
 *
 * The output array locations `outneg[-n*step]...outneg[-1*step]' must
 * not overlap with any of the input arrays.
 * 
 * Input:
 * 	(real *)outneg		Storage for `outneg[-n*step]', ...,
 *				  `outneg[-2*step]', `outneg[-1*step]'
 *				  must be allocated.
 *
 *	(int)step		This positive integer is the increment to
 *				  use between output array elements.
 *
 *	(const real *)inneg	Values `inneg[-n], ..., inneg[-1]' must
 *				  be allocated and defined.
 *
 *	(const real *)inpos	Values `inpos[0]...inpos[rise->final]' must
 *				  be allocated and defined.
 *
 *	(int)n			This positive integer gives the length of
 *				  the output and one of the input arrays.
 *
 *	(const interval *)rise	This data structure must contain a valid
 *				  rising cutoff function.
 *
 * Assumptions:
 *	(1)  n >= rise->final >= 0
 *	(2)  step > 0
 *	
 * Output:
 *   Array values `outneg[-n*step], ..., outneg[-1*step]' are replaced
 *      with linear combinations of `inneg[]' and `inpos[]'.
 */
extern void
  fdcn(
       real *outneg,		/* Negative or left-half output array. */
       int   step,		/* Output array increment: step > 0. */
       const real *inneg,	/* Negative or left-half input array. */
       const real *inpos,	/* Positive or right-half input array. */
       int   n,			/* Length of `inneg[]' and `outneg[]'. */
       const interval *rise)	/* Sampled rising cutoff function. */
{
  const real *rpos, *rneg;
  int k;

  assert(rise);
  assert(rise->final>=0);
  assert(n >= rise->final);
  assert(step>0);

  rpos = rise->origin;
  rneg = rise->origin;

  k=0;
  while( k <= rise->final )
    {
      outneg -= step;
      *outneg = (*rpos++)*(*--inneg) - (*--rneg)*(*inpos++);
      ++k;
    }
  while(  k < n )
    {
      outneg -= step;
      *outneg = *--inneg;
      ++k;
    }     
  return;
}

/*********************************************************************
 * fdcp()
 *
 * [F]old a given array into a [D]isjoint output array, with [C]osine
 * polarity (odd on the left and even on the right), into the [P]ositive
 * part of the output array.  Elements `outpos[0],...,outpos[(n-1)*step]'
 * of the output array are assigned. The output array locations `outpos[1]',
 * `outpos[1]', ..., `outpos[(n-1)*step]' must not overlap with any of the
 * input arrays.
 *
 * Calling sequence and basic algorithm:
 *
 *      fdcp( OPOS, STEP, INEG, IPOS, N, RISE ):
 *        For K = 0 to RISE.FINAL
 *          Let OPOS[K*STEP] = RISE.ORIGIN[K] * IPOS[K]
 *                              + RISE.ORIGIN[-1-K] * INEG[-1-K]
 *        For K = RISE.FINAL+1 to N-1
 *          OPOS[K*STEP]=IPOS[K]
 * 
 * Input:
 * 	(real *)outpos		Storage for `outpos[0]',`outpos[1*step]',
 *				   ..., `outpos[(n-1)*step]' must be 
 *				  pre-allocated.
 *
 *	(int)step		This positive integer is the increment to
 *				  use between output array elements.
 *
 *	(const real *)inneg	Values `inneg[rise->least], ..., inneg[-1]'
 *				  must be allocated and defined.
 *
 *	(const real *)inpos	Values `inpos[0]...inpos[n-1]' must be
 *				  allocated and defined.
 *
 *	(int)n			This positive integer gives the length of
 *				  the output and one of the input arrays.
 *
 *	(const interval *)rise	This data structure must contain a valid
 *				  rising cutoff function.
 *
 * Assumptions:
 *	(1)  n >= rise->final >= 0
 *	(2)  step > 0
 *	
 * Output:
 *   Array values `outpos[0]', ..., `outpos[(n-1)*step]' are replaced
 *     with linear combinations of `inneg[]' and `inpos[]'.
 */
extern void
  fdcp(
       real *outpos,		/* Positive or right-half output array. */
       int   step,		/* Output array increment: step > 0. */
       const real *inneg,	/* Negative or left-half input array. */
       const real *inpos,	/* Positive or right-half input array. */
       int   n,			/* Length of `inpos[]' and `outpos[]'. */
       const interval *rise)	/* Sampled rising cutoff function. */
{
  const real *rpos, *rneg;
  int k;

  assert(rise);
  assert(rise->final>=0);
  assert(n >= rise->final);
  assert(step>0);

  rpos = rise->origin;
  rneg = rise->origin;

  k=0;
  while( k <= rise->final )
    {
      *outpos = (*rpos++)*(*inpos++) + (*--rneg)*(*--inneg);
      outpos += step;
      ++k;
    }
  while(  k < n )
    {
      *outpos = *inpos++;
      outpos += step;
      ++k;
    }     
  return;
}

/*********************************************************************
 * fdsn()		
 *
 * [F]old a given array into a [D]isjoint output array, with [S]ine
 * polarity (odd on the left and even on the right), into the [N]egative
 * part of the output array.  Elements `outneg[-n*step],...,outneg[-1*step]'
 * of the output array are assigned.
 *
 * Calling sequence and basic algorithm:
 *
 *    fdsn( ONEG, STEP, INEG, IPOS, N, RISE ):
 *      For K = -N to RISE.LEAST-1
 *        Let ONEG[K*STEP] = INEG[K]
 *      For K = RISE.LEAST to -1
 *        Let ONEG[K*STEP] = RISE.ORIGIN[-1-K] * INEG[K]
 *                            + RISE.ORIGIN[K] * IPOS[-1-K]
 *
 * The output array locations `outneg[-n*step]...outneg[-1*step]' must
 * not overlap with any of the input arrays.
 * 
 * Input:
 * 	(real *)outneg		Storage for `outneg[-n*step]', ...,
 *				  `outneg[-2*step]', `outneg[-1*step]'
 *				  must be allocated.
 *
 *	(int)step		This positive integer is the increment to
 *				  use between output array elements.
 *
 *	(const real *)inneg	Values `inneg[-n], ..., inneg[-1]' must
 *				  be allocated and defined.
 *
 *	(const real *)inpos	Values `inpos[0]...inpos[rise->final]' must
 *				  be allocated and defined.
 *
 *	(int)n			This positive integer gives the length of
 *				  the output and one of the input arrays.
 *
 *	(const interval *)rise	This data structure must contain a valid
 *				  rising cutoff function.
 *
 * Assumptions:
 *	(1)  n >= rise->final >= 0
 *	(2)  step > 0
 *	
 * Output:
 *   Array values `outneg[-n*step], ..., outneg[-1*step]' are replaced
 *      with linear combinations of `inneg[]' and `inpos[]'.
 */
extern void
  fdsn(
       real *outneg,		/* Negative or left-half output array. */
       int   step,		/* Output array increment: step > 0. */
       const real *inneg,	/* Negative or left-half input array. */
       const real *inpos,	/* Positive or right-half input array. */
       int   n,			/* Length of `inneg[]' and `outneg[]'. */
       const interval *rise)	/* Sampled rising cutoff function. */
{
  const real *rpos, *rneg;
  int k;

  assert(rise);
  assert(rise->final>=0);
  assert(n >= rise->final);
  assert(step>0);

  rpos = rise->origin;
  rneg = rise->origin;

  k=0;
  while( k <= rise->final )
    {
      outneg -= step;
      *outneg = (*rpos++)*(*--inneg) + (*--rneg)*(*inpos++);
      ++k;
    }
  while(  k < n )
    {
      outneg -= step;
      *outneg = *--inneg;
      ++k;
    }     
  return;
}

/*********************************************************************
 * fdsp()		
 *
 * [F]old a given array into a [D]isjoint output array, with [S]ine
 * polarity (odd on the left and even on the right), into the [P]ositive
 * part of the output array.  Elements `outpos[0],...,outpos[(n-1)*step]'
 * of the output array are assigned.
 *
 * Calling sequence and basic algorithm:
 *
 *      fdsp( OPOS, STEP, INEG, IPOS, N, RISE ):
 *        For K = 0 to RISE.FINAL
 *          Let OPOS[K*STEP] = RISE.ORIGIN[K] * IPOS[K]
 *                              - RISE.ORIGIN[-1-K] * INEG[-1-K]
 *        For K = RISE.FINAL+1 to N-1
 *          OPOS[K*STEP]=IPOS[K]
 * 
 * Input:
 * 	(real *)outpos		Storage for `outpos[0]',`outpos[1*step]',
 *				   ..., `outpos[(n-1)*step]' must be 
 *				  pre-allocated.
 *
 *	(int)step		This positive integer is the increment to
 *				  use between output array elements.
 *
 *	(const real *)inneg	Values `inneg[rise->least], ..., inneg[-1]'
 *				  must be allocated and defined.
 *
 *	(const real *)inpos	Values `inpos[0]...inpos[n-1]' must be
 *				  allocated and defined.
 *
 *	(int)n			This positive integer gives the length of
 *				  the output and one of the input arrays.
 *
 *	(const interval *)rise	This data structure must contain a valid
 *				  rising cutoff function.
 *
 * Assumptions:
 *	(1)  n >= rise->final >= 0
 *	(2)  step > 0
 *	
 * Output:
 *   Array values `outpos[0]', ..., `outpos[(n-1)*step]' are replaced
 *     with linear combinations of `inneg[]' and `inpos[]'.
 */
extern void
  fdsp(
       real *outpos,		/* Positive or right-half output array. */
       int   step,		/* Output array increment: step > 0. */
       const real *inneg,	/* Negative or left-half input array. */
       const real *inpos,	/* Positive or right-half input array. */
       int   n,			/* Length of `inpos[]' and `outpos[]'. */
       const interval *rise)	/* Sampled rising cutoff function. */
{
  const real *rpos, *rneg;
  int k;

  assert(rise);
  assert(rise->final>=0);
  assert(n >= rise->final);
  assert(step>0);

  rpos = rise->origin;
  rneg = rise->origin;

  k=0;
  while( k <= rise->final )
    {
      *outpos = (*rpos++)*(*inpos++) - (*--rneg)*(*--inneg);
      outpos += step;
      ++k;
    }
  while(  k < n )
    {
      *outpos = *inpos++;
      outpos += step;
      ++k;
    }     
  return;
}

/*********************************************************************
 * fipc()		
 *
 *  [F]old [I]n-[P]lace, [C]osine polarity.
 *
 *  Given an array, smoothly fold the odd part into the left half
 *  and the even part into the right half.  This is implemented as
 *  a transformation in place of two arrays: the negative or left half
 *  `oneg[]', and the positive or right half `opos[]'.  The indexing
 *  is chosen so that `opos' and `oneg' are typically identical
 *  pointers to the first element of a  block of the given array.
 *  The function then folds the leading edge of the  block into the
 *  trailing edge of the previous block.  However, the arrays
 *  `oneg[]' and `opos[]' must not overlap except possibly at 0.
 *  Temporary variables are used to allow transformation in place. 
 *
 * Calling sequence and basic algorithm:
 *
 *      fipc( ONEG, OPOS, RISE ):
 *        For K = 0 to RISE.FINAL
 *          Let TEMP = RISE.ORIGIN[K] * OPOS[K]
 *                      + RISE.ORIGIN[-K-1] * ONEG[-K-1]
 *          Let ONEG[-K-1] = RISE.ORIGIN[K] * ONEG[-K-1]
 *                            - RISE.ORIGIN[-K-1] * OPOS[K]
 *          Let OPOS[K] = TEMP
 * 
 * Input:
 *	(real *)oneg		Values `oneg[rise->least],...,oneg[-1]'
 *				  must be allocated and defined.
 *
 * 	(real *)opos		Values `opos[0],...,opos[rise->final]'
 *				  must be  allocated and defined.
 *
 *	(const interval *)rise	This data structure must contain a valid
 *				  rising cutoff function.
 *	
 * Output:
 *   Array values `oneg[rise->least]', ..., `oneg[-1]' and `opos[0]',...,
 *   `opos[rise->final]' are replaced with orthogonal linear combinations.
 *
 * Assumptions:
 *	(1)  rise->final >= 0
 */
extern void
  fipc(
       real *oneg,		/* Negative or left-half input/output. */
       real *opos,		/* Positive or right-half input/output. */
       const interval *rise)	/* Sampled rising cutoff function. */
{
  int k;
  real temp;
  const real *rneg, *rpos;

  assert(rise);
  assert(rise->final>=0);

  rneg = rise->origin;
  rpos = rise->origin;

  for( k=0; k<=rise->final; k++)
    {
      temp  = (*rpos)*(*opos) + (*--rneg)*(*--oneg);
      *oneg  = (*rpos++)*(*oneg) - (*rneg)*(*opos);
      *opos++ = temp;
    }
  return;
}

/***************************************************************************
 * fips()		
 *
 *  [F]old [I]n-[P]lace, [S]ine polarity.
 *
 *  Given an array, smoothly fold the odd part into the left half
 *  and the even part into the right half.  This is implemented as
 *  a transformation in place of two arrays: the negative or left half
 *  `oneg[]', and the positive or right half `opos[]'.  The indexing
 *  is chosen so that `opos' and `oneg' are typically identical
 *  pointers to the first element of a  block of the given array.
 *  The function then folds the leading edge of the  block into the
 *  trailing edge of the previous block.  However, the arrays
 *  `oneg[]' and `opos[]' must not overlap except possibly at 0.
 *  Temporary variables are used to allow transformation in place. 
 *
 * Calling sequence and basic algorithm:
 *
 *      fips( ONEG, OPOS, RISE ):
 *        For K = 0 to RISE.FINAL
 *          Let TEMP = RISE.ORIGIN[K] * OPOS[K]
 *                      - RISE.ORIGIN[-K-1] * ONEG[-K-1]
 *          Let ONEG[-K-1] = RISE.ORIGIN[K] * ONEG[-K-1]
 *                            + RISE.ORIGIN[-K-1] * OPOS[K]
 *          Let OPOS[K] = TEMP
 * 
 * Input:
 *	(real *)oneg		Values `oneg[rise->least],...,oneg[-1]'
 *				  must be allocated and defined.
 *
 * 	(real *)opos		Values `opos[0],...,opos[rise->final]'
 *				  must be  allocated and defined.
 *
 *	(const interval *)rise	This data structure must contain a valid
 *				  rising cutoff function.
 *	
 * Output:
 *   Array values `oneg[rise->least]', ..., `oneg[-1]' and `opos[0]',...,
 *   `opos[rise->final]' are replaced with orthogonal linear combinations.
 *
 * Assumptions:
 *	(1)  rise->final >= 0
 */
extern void
  fips(
       real *oneg,		/* Negative or left-half input/output. */
       real *opos,		/* Positive or right-half input/output. */
       const interval *rise)	/* Sampled rising cutoff function. */
{
  int k;
  real temp;
  const real *rneg, *rpos;

  assert(rise);
  assert(rise->final>=0);

  rneg = rise->origin;
  rpos = rise->origin;

  for( k=0; k<=rise->final; k++)
    {
      temp  = (*rpos)*(*opos) - (*--rneg)*(*--oneg);
      *oneg  = (*rpos++)*(*oneg) + (*rneg)*(*opos);
      *opos++ = temp;
    }
  return;
}
