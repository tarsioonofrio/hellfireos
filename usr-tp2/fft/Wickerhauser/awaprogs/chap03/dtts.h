/*
 * Discrete sine and cosine transforms based on FFT.
 * 
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef DTTS_HDR_ALREADY_INCLUDED
#define DTTS_HDR_ALREADY_INCLUDED

#include "real.h"

extern void
  dct_i(
	real *x,	/* Input; length (1<<q)+1. */
	int q);		/* Nonnegative integer. */

extern void
  dst_i(
	real *x,	/* Input; length (1<<q). */
	int q);		/* Nonnegative integer. */

extern void
  dct_ii(
	 real *x,	/* Input; length (1<<q). */
	 int q);	/* Nonnegative integer. */

extern void
  dst_ii(
	 real *x,	/* Input; length (1<<q). */
	 int q);	/* Nonnegative integer. */

extern void
  dct_iii(
	  real *x,	/* Input; length (1<<q). */
	  int q);	/* Nonnegative integer. */

extern void
  dst_iii(
	  real *x,	/* Input; length (1<<q). */
	  int q);	/* Nonnegative integer. */

extern void
  dct_iv(
	 real *x,	/* Input; length (1<<q). */
	 int q);	/* Nonnegative integer. */

extern void
  dst_iv(
	 real *x,	/* Input; length (1<<q). */
	 int q);	/* Nonnegative integer. */

#endif /* DTTS_HDR_ALREADY_INCLUDED */
