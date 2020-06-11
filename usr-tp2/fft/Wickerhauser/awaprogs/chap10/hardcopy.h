/*
 * Declare functions to prepare PostScript files depicting 1-dimensional
 * signals and their idealized time-frequency plane density plots.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef HARDCOPY_HDR_ALREADY_INCLUDED
# define HARDCOPY_HDR_ALREADY_INCLUDED

#include <stdio.h>
#include "real.h"
#include "tfa.h"

/* The following is a system-dependent minimum displayable gray-level. */
#define		MINGRAY		(0.001)

/* The numbers below are the coordinates (in PostScript points) of the 
 * lower left (LL*) and upper right (UR*) corners of the BoundingBox
 * which contains the plotted signal.
 */
#define		LLXS	72
#define		LLYS	72
#define		URXS	528
#define		URYS	254

/* The following BoundingBox coordinates correspond to a square centered 
 * in the upper part of an 8.5" x 11" sheet of paper.  It represents the
 * analysis by time-frequency atoms.
 */
#define		LLXA	 72
#define		LLYA	264
#define		URXA	528
#define		URYA	720

extern void
  plotsig(
	  FILE *psfile,		/* Write into this file stream. */
	  const real  *signal,	/* Array of `length' Y-values.  */
	  int length);		/* Positive integer. */

extern void
  tfa1s2ps(
	   FILE *psfile,	/* Name of the output file to use.  */
	   int samples,		/* Number of samples in the signal. */
	   tfa1 *tfarray,	/* Array of TFA1 data structures.   */
	   int num);		/* Number of TFA1 data structures.  */

#endif /* HARDCOPY_HDR_ALREADY_INCLUDED */
