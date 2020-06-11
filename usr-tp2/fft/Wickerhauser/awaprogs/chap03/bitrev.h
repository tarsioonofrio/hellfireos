/* 
 * Declare bit-reversal functions.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef BITREV_HDR_ALREADY_INCLUDED
# define BITREV_HDR_ALREADY_INCLUDED

extern int
  br( 
     int n,			/* Nonnegative integer to bit-reverse. */
     int log2len);		/* Reverse this many bits.  */

extern void
  bitrevd(
	  void *out,		/* Pointer to base of the output array. */
	  const void *in,	/* Pointer to base of the input array. */
	  int   q,		/* # of index bits, or elements of `in[]'. */
	  int   size);		/* Number of bytes in an `in[]' element.  */

extern void
  bitrevi(
	  void *x,		/* Pointer to the input/output array. */
	  int   q,		/* # of index bits, or elements of `x[]'. */
	  int   size);		/* Number of bytes in an `x[]' element.  */


#endif /* BITREV_HDR_ALREADY_INCLUDED */
