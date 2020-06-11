/*
 * 
 * Header for adapted local trigonometric transform functions.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef ALTT1_HDR_ALREADY_INCLUDED
#define ALTT1_HDR_ALREADY_INCLUDED

#include "real.h"
#include "interval.h"		/* Rising cutoff arrays are INTERVALs. */
#include "btn.h"
#include "hedge.h"

extern interval *
  initrcf(
	  int e);		/* Range of midpoint folding. */

extern void
  lcadf(
	real *parent,		/* Input and output array binary tree. */
	int n,			/* Length of one row in the tree.      */
	int L,			/* Number of levels in the tree.       */
	interval *rise);	/* Sampled rising cutoff function.   */

extern interval **
  initrcfs(
	   int n,		/* Length of signal to decompose. */
	   int L);		/* Bottom decomposition level.    */

extern btn *
  initlcabtn(
	     const real *in,	/* Input signal */
	     int n,		/* Length of the input signal. */
	     interval *rise);	/* Sampled rising cutoff function. */

extern void
  lcadm(
	btn *root,		/* Current root node in the BTN tree. */
	int s,			/* Current level in the tree. */
	int L,			/* Maximum level in the tree. */
	interval **rs);		/* Sampled rising cutoff functions. */

extern real *
  lcadf2hedge(
	      hedge *graph,	/* Partially assigned HEDGE for the output. */
	      const real *in,	/* Input signal to analyze. */
	      int n,		/* Length of the input signal. */
	      int L);		/* Maximum level in the analysis tree. */

extern btn *
  lcadm2hedge(
	      hedge *graph,	/* Partially assigned HEDGE for the output. */
	      const real *in,	/* Input signal to analyze. */
	      int n,		/* Length of the input signal. */
	      int L);		/* Maximum level in the analysis tree. */

extern void
  lcsdf(
	hedge *graph,		/* Pointers to subintervals in an array. */
	int n,			/* Length of the original signal.        */
	interval *rise);	/* Sampled rising cutoff function.       */

#endif /* ALTT1_HDR_ALREADY_INCLUDED */
