/*
 * 
 * This header file declares basic 1-dimensional local trigonometric
 *  transform and smooth local periodization functions.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef LTT1_HDR_ALREADY_INCLUDED
#define LTT1_HDR_ALREADY_INCLUDED

#include "real.h"
#include "interval.h"		/* Rising cutoff arrays are INTERVALs. */

extern void
  lct(
      real *sig,		/* Input and output array.           */
      int n,			/* Length of the input/output array. */
      const interval *rise);	/* Sampled rising cutoff function.   */

extern void
  lst(
      real *sig,		/* Input and output array.           */
      int n,			/* Length of the input/output array. */
      const interval *rise);	/* Sampled rising cutoff function.   */

extern void
  ilct(
       real *sig,		/* Input and output array.           */
       int n,			/* Length of the input/output array. */
       const interval *rise);	/* Sampled rising cutoff function.   */

extern void
  ilst(
       real *sig,		/* Input and output array.           */
       int n,			/* Length of the input/output array. */
       const interval *rise);	/* Sampled rising cutoff function.   */

extern void
  lpds(
       real *out,		/* Preallocated output array.      */
       const real *in,		/* Preallocated input array.       */
       int n,			/* Length to which we periodize.   */
       const interval *rise);	/* Sampled rising cutoff function. */

extern void
  lpic(
       real *sig,		/* Input and output array.           */
       int n,			/* Length of the input/output array. */
       const interval *rise);	/* Sampled rising cutoff function.   */

extern void
  lpds(
       real *out,		/* Preallocated output array.      */
       const real *in,		/* Preallocated input array.       */
       int n,			/* Length to which we periodize.   */
       const interval *rise);	/* Sampled rising cutoff function. */

extern void
  lpis(
       real *sig,		/* Input and output array.           */
       int n,			/* Length of the input/output array. */
       const interval *rise);	/* Sampled rising cutoff function.   */

extern void
  ilpic(
	real *sig,		/* Input and output array.           */
	int n,			/* Length of the input/output array. */
	const interval *rise);	/* Sampled rising cutoff function.   */

extern void
  ilpis(
       real *sig,		/* Input and output array.           */
       int n,			/* Length of the input/output array. */
       const interval *rise);	/* Sampled rising cutoff function.   */

extern void
  lpica(
	real *sig,		/* Input and output array.           */
	int *lengths,		/* Lengths of the subintervals.      */
	int  num,		/* Number of subintervals.           */
	const interval *rise);	/* Sampled rising cutoff function.   */

extern void
  lpisa(
	real *sig,		/* Input and output array.           */
	int *lengths,		/* Lengths of the subintervals.      */
	int  num,		/* Number of subintervals.           */
	const interval *rise);	/* Sampled rising cutoff function.   */

extern void
  ilpica(
	 real *sig,		/* Input and output array.           */
	 int *lengths,		/* Lengths of the subintervals.      */
	 int  num,		/* Number of subintervals.           */
	 const interval *rise);	/* Sampled rising cutoff function.   */

extern void
  ilpisa(
	 real *sig,		/* Input and output array.           */
	 int *lengths,		/* Lengths of the subintervals.      */
	 int  num,		/* Number of subintervals.           */
	 const interval *rise);	/* Sampled rising cutoff function.   */

#endif /* LTT1_HDR_ALREADY_INCLUDED */

