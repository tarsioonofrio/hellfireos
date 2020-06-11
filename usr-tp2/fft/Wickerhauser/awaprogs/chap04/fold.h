/* 
 * Declare functions to perform folding and unfolding in one dimension.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */


#ifndef FOLD_HDR_ALREADY_INCLUDED
# define FOLD_HDR_ALREADY_INCLUDED

#include "real.h"
#include "interval.h"

/* Function to perform "in-place" folding: */
typedef void
  (*fi_fn)(
	   real * ,		/* Negative (left) side: -reach <= x < 0. */
	   real * ,		/* Positive (right) side: 0 <= x < reach. */
	   const interval * );	/* Rising cutoff function. */

/* Function to perform folding and copying into a disjoint array: */
typedef void
  (*fd_fn)(
	   real * ,		/* Output array: -N<=x<0 or 0<=x<N. */
	   int ,		/* Output array increment: step > 0. */
	   const real * ,	/* Negative (left) input: -N <= x < 0. */
	   const real * ,	/* Positive (right) input: 0 <= x < reach. */
	   int ,		/* Length of the longer input array. */
	   const interval * );	/* Rise increases on -reach <= x < reach. */

/* Disjoint folding and copying functions: */
extern void
  fdcn(
       real *outneg,		/* Negative or left-half output array. */
       int   step,		/* Output array increment: step > 0. */
       const real *inneg,	/* Negative or left-half input array. */
       const real *inpos,	/* Positive or right-half input array. */
       int   n,			/* Length of `inneg[]' and `outneg[]'. */
       const interval *rise);	/* Sampled rising cutoff function. */

extern void
  fdcp(
       real *outpos,		/* Positive or right-half output array. */
       int   step,		/* Output array increment: step > 0. */
       const real *inneg,	/* Negative or left-half input array. */
       const real *inpos,	/* Positive or right-half input array. */
       int   n,			/* Length of `inpos[]' and `outpos[]'. */
       const interval *rise);	/* Sampled rising cutoff function. */

extern void
  fdsn(
       real *outneg,		/* Negative or left-half output array. */
       int   step,		/* Output array increment: step > 0. */
       const real *inneg,	/* Negative or left-half input array. */
       const real *inpos,	/* Positive or right-half input array. */
       int   n,			/* Length of `inneg[]' and `outneg[]'. */
       const interval *rise);	/* Sampled rising cutoff function. */

extern void
  fdsp(
       real *outpos,		/* Positive or right-half output array. */
       int   step,		/* Output array increment: step > 0. */
       const real *inneg,	/* Negative or left-half input array. */
       const real *inpos,	/* Positive or right-half input array. */
       int   n,			/* Length of `inpos[]' and `outpos[]'. */
       const interval *rise);	/* Sampled rising cutoff function. */

/* The adjoint of sine-polarity folding is cosine-polarity folding: */
#define udcn fdsn		/* [U]nfold [D]isjoint [C]osine [N]egative. */
#define udcp fdsp		/* [U]nfold [D]isjoint [C]osine [P]ositive. */
#define udsn fdcn		/* [U]nfold [D]isjoint [S]ine [N]egative. */
#define udsp fdcp		/* [U]nfold [D]isjoint [S]ine [P]ositive. */

/* Functions to perform folding in place: */
extern void
  fipc(
       real *oneg,		/* Negative or left-half input/output. */
       real *opos,		/* Positive or right-half input/output. */
       const interval *rise);	/* Sampled rising cutoff function. */

extern void
  fips(
       real *oneg,		/* Negative or left-half input/output. */
       real *opos,		/* Positive or right-half input/output. */
       const interval *rise);	/* Sampled rising cutoff function. */

/* The adjoint of fips() is fipc(): */
#define uipc fips  /* [U]nfold [I]n-[P]lace, [C]osine polarity */
#define uips fipc  /* [U]nfold [I]n-[P]lace, [S]ine polarity */

#endif /* FOLD_HDR_ALREADY_INCLUDED */
