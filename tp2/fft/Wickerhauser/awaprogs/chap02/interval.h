/* 
 * This header file declares the `interval' data type and some functions
 * which manipulate it.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef INTERVAL_HDR_ALREADY_INCLUDED
#define INTERVAL_HDR_ALREADY_INCLUDED

#include "real.h"

typedef struct {
  real *origin;			/* Pointer to zero-indexed element */
  int   least;			/* Least valid index of `origin[]' */
  int   final;			/* Final valid index of `origin[]' */
} interval;

extern interval *
  makeinterval( 
	       const real *data,/* NULL, or a valid array starting at 0. */
	       int least,	/* Least valid index of output `origin'. */
	       int final );	/* Final valid index of output `origin'. */

extern interval *
  freeinterval( 
	       interval *seg);	/* `interval' data structure to free. */

extern interval *
  enlargeinterval( 
		  interval *old, /* NULL, or a valid data array. */
		  int least,	/* New least `origin[]' index.   */
		  int final );	/* New final `origin[]' index.   */

extern int
  ininterval( 
	     interval *segment,	/* `interval' data structure to test. */
	     int offset);	/* Index to validate. */

extern int
  intervalstotal( 
		 interval *in,	/* Array of `interval' data structures. */
		 int num);	/* Number of `interval's in `in[]'. */

#endif /* INTERVAL_HDR_ALREADY_INCLUDED */

