/* 
 * Declare various rising cutoff functions for folding and unfolding.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef RCF_HDR_ALREADY_INCLUDED
#define RCF_HDR_ALREADY_INCLUDED

#include "real.h"
#include "interval.h"
#include "common.h"		/* Defines PI */

#define rcf(t)		rcfis(1,t) /* Default: one iteration of sine. */

#define theta1(t)	(0.25*PI*(1.0 + (t)))
#define theta3(t)	(0.125*PI*(2.0+(t)*(3.0 - (t)*(t))))
#define th		theta1	/* Default: one iteration of sine. */
#define rcfth0(t)	(t>-1.0?(t<1.0?sin(th(t)):1.0):0.0)
#define rcfth(t)	(t>-1.0?(t<0?cos(th(-(t))):(t<1.0?sin(th(t)):1.0)):0)

extern real
  rcfis(			/* Iterated-sine rising cutoff function. */
	int  n,			/* How many iterations of sin() to use.  */
	real t);		/* Point of evaluation. */

extern void
  rcfgrid(			/* Rising cutoff sampled at gridpoints. */
	  interval *r);		/* Put samples in this data structure. */

extern void
  rcfmidp(			/* Rising cutoff sampled at midpoints. */
	  interval *r);		/* Put samples in this data structure. */


#endif /* RCF_HDR_ALREADY_INCLUDED */
