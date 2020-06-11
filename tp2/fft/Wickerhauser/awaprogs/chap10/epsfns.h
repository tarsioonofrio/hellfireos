/*
 * Declare functions to write PostScript drawing primitives to a stream.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef EPSFNS_HDR_ALREADY_INCLUDED
# define EPSFNS_HDR_ALREADY_INCLUDED

#include <stdio.h>		/* Declaration of FILE */
#include "real.h"

extern void
  epsepilogue(
	      FILE  *psfile);	/* File to write into. */

extern void
  epsfrect(
	   FILE  *psfile,	/* File to write into. */
	   real xmin,		/* Least x-coord (in [0.0, 1.0]). */
	   real ymin,		/* Least y-coord (in [0.0, 1.0]). */
	   real xmax,		/* Greatest x-coord (in [0.0, 1.0]). */
	   real ymax,		/* Greatest y-coord (in [0.0, 1.0]). */
	   real gray);		/* Gray level: 0.0=black, 1.0=white. */

extern void
  epslineto(
	    FILE  *psfile,	/* File into which we will write. */
	    real xval,		/* x-coord (in [0.0, 1.0]).       */
	    real yval);		/* y-coord (in [0.0, 1.0]).       */

extern void
  epsmoveto(
	    FILE  *psfile,	/* File into which we will write. */
	    real xval,		/* x-coord (in [0.0, 1.0]).       */
	    real yval);		/* y-coord (in [0.0, 1.0]).       */


extern void
  epsprologue(
	      FILE *psfile,	/* Print to this file.  */
	      int   bbxmin,	/* X-coord of lower left BoundingBox corner. */
	      int   bbymin,	/* Y-coord of lower left BoundingBox corner. */
	      int   bbxmax,	/* X-coord of upper left BoundingBox corner. */
	      int   bbymax);	/* Y-coord of upper left BoundingBox corner. */

extern void
  epsstroke(
	    FILE  *psfile);	/* File into which we will write. */

extern void
  epstranslate(
	       FILE *psfile,	/* Print to this file. */
	       int   xptval,	/* X-coord of translated origin. */
	       int   yptval);	/* Y-coord of translated origin. */

#endif /* EPSFNS_HDR_ALREADY_INCLUDED */
