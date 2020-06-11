/*
 * 
 * Bit-reversal functions.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <assert.h>
#include <stdlib.h>		/* For malloc() and free() */
#include <string.h>		/* For memcpy() */
#include "bitrev.h"

/********************************************************************
 * br()
 *
 *	  Return the input integer bit-reversed.
 *
 *  Calling sequence and basic algorithm:
 *
 *   br( N, LOG2LEN ):
 *      Let U = N&1
 *      For J = 1 to LOG2LEN 
 *	  N >>= 1
 *	  U <<= 1
 *	  U += N&1
 *      Return U
 *
 *
 *  Inputs:
 *	(int)n		This is the nonnegative input value.
 *	(int)log2len	This is the positive number of bits to reverse.
 *
 * Outputs:
 *	(int)br		The  return value is `n' with its bits reversed.
 *
 * Assumptions:
 *	1.  n >=0
 *	2.  n < (1<<log2len)
 */
extern int
  br( 
     int n,			/* Nonnegative integer to bit-reverse. */
     int log2len)		/* Reverse this many bits.  */
{
  int u;

  assert(n>=0);
  assert(n<(1<<log2len));

  u = n&1;
  while( --log2len )
    {
      n >>= 1;
      u <<= 1;
      u += n&1;
    }
  return(u);
}

/********************************************************************
 * bitrevd()
 *
 *	Permute to a disjoint array by index bit-reversal.
 *
 *  Calling sequence and basic algorithm:
 *
 *  bitrevd( OUT, IN, Q ):
 *     Let M = 1<<Q
 *     For N = 0 to M-1
 *        Let U = br(N, Q)
 *        Let OUT[U] = IN[N]
 *
 *  Inputs:
 *	(void *)out	This is the start of the output array.
 *	(void *)in	This is the start of the input array.
 *	(int)q		This is the positive number of bits in the indices.
 *	(int)size	This is the size (in bytes) of elements of `in[]'
 *			  or `out[]'.
 *
 *  Outputs:
 *	Elements of `in[]' are assigned to `out[]' by side effect.
 *
 *  Assumptions:
 *	1.  q >=0
 *	2.  size > 0
 *	3.  in != NULL
 *	4.  out != NULL
 *
 *  External functions called:
 *	br(), memcpy();
 */
extern void
  bitrevd(
	  void *out,		/* Pointer to base of the output array. */
	  const void *in,	/* Pointer to base of the input array. */
	  int   q,		/* # of index bits, or elements of `in[]'. */
	  int   size)		/* Number of bytes in an `in[]' element.  */
{
  int u, n;

  assert(in);
  assert(out);
  assert(q>=0);
  assert(size>0);

  for(n=0; n<(1<<q); n++)
    {
      u = br(n,q);
      memcpy( out+n*size, in+u*size, size );
    }
  return;
}

/********************************************************************
 * bitrevi()
 *
 *	Permute an array in place by index bit-reversal.
 *
 *  Calling sequence and basic algorithm:
 *
 *  bitrevi( X, Q ):
 *     Let M = 1<<Q
 *     For N = 0 to M-1
 *        Let U = br(N, Q)
 *        If U > N then
 *           Let TEMP = X[N]
 *           Let X[N] = X[U]
 *           Let X[U] = TEMP
 *
 *  Inputs:
 *	(void *)x	This is the joint input and output array.
 *	(int)q		This is the positive number of bits in the indices.
 *	(int)size	This is the size (in bytes) of elements of `x[]'
 *
 *  Outputs:
 *	Elements of `x[]' are permuted by side effect exchanges.
 *
 *  Assumptions:
 *	1.  q >=0
 *	2.  size > 0
 *	3.  x != NULL
 *
 *  External functions called:
 *	br(), memcpy(), malloc(), free();
 */
extern void
  bitrevi(
	  void *x,		/* Pointer to the input/output array. */
	  int   q,		/* # of index bits, or elements of `x[]'. */
	  int   size)		/* Number of bytes in an `x[]' element.  */
{
  int u, n;
  void *temp, *xn, *xu;

  assert(x);
  assert(q>=0);
  assert(size>0);

  temp = malloc(size);
  for(n=0; n<(1<<q); n++)
    {
      u = br(n,q);
      if( u > n )
	{
	  xu  = x+u*size;
	  xn =  x+n*size;
	  memcpy( temp, xu, size );
	  memcpy( xu, xn, size );
	  memcpy( xn, temp, size );
	}
    }
  free(temp);
  return;
}

