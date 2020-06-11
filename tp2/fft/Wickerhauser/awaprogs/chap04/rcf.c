/*
 *
 * These utilties produce rising cutoff functions needed for folding.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <stdlib.h>		/* For malloc() */
#include <assert.h>
#include <math.h>		/* For sin() */
#include "interval.h"
#include "rcf.h"
#include "common.h"

/***********************************************************************
 *  rcfis()
 *
 *  Return evaluations of a rising cutoff function defined on [-1,1] by
 *  the iterated sine formula.
 *
 *  Calling sequence and basic algorithm:
 *    rcfis( N, T ):
 *       If T > -1.0 then
 *          If T < 1.0 then
 *             For I = 0 to N-1
 *                Let R = sin( 0.5 * PI * T )
 *             Let R = sin( 0.25 * PI * (1.0 + R) )
 *          Else
 *             Let R = 1.0
 *       Else
 *          Let R = 0.0
 *       Return R
 *
 *  Input:
 *	(int)n		This is the small nonnegative number of iterations.
 *	(real)t		This is any real number.
 *
 *  Output:
 *	(real)rcfis	Return an evaluation of a rising cutoff function.
 */
extern real
  rcfis(			/* Iterated-sine rising cutoff function. */
	int  n,			/* How many iterations of sin() to use.  */
	real t)			/* Point of evaluation. */
{
  int i;

  assert(n>=0);
  assert(n<=32);

  if(t > -1.0)
    if(t < 1.0)
      {
	for(i=0; i<n; i++) 
	  t = sin(0.5*PI*t);	/* iterate the sine function. */
	t  = sin(0.25*PI*(1.0+t));
      }
    else
      t = 1.0;
  else
    t = 0.0;
  return(t);
}


/***********************************************************************
 *  rcfgrid()
 *
 *  Sample a rising cutoff function at gridpoints.
 *
 *  Calling sequence and basic algorithm:
 *
 *       rcfgrid( R ):
 *          Let X  = 0.0
 *          Let DX = 1.0/(R.FINAL+1.0)
 *          Let R.ORIGIN[0] = sqrt(0.5)
 *          For J = 1 to R.FINAL
 *             X += DX
 *             Let R.ORIGIN[J]  = rcf(X)
 *             Let R.ORIGIN[-J] = rcf(-X)
 *
 *  Input:
 *	(interval *)r	This is a preallocated INTERVAL.
 *
 *  Output:
 *	(interval *)r	This INTERVAL is filled with samples.
 *
 *  Assumptions:
 *	1. r != NULL
 *	2. r->origin != NULL
 *	3. r->final >= 0
 *	4. r->least == -(r->final)
 */
extern void
  rcfgrid(			/* Rising cutoff sampled at gridpoints. */
	  interval *r)		/* Put samples in this data structure. */
{
  real x, dx;
  int j;

  assert(r);
  assert(r->origin);
  assert(r->final>=0);
  assert(r->least== -(r->final));

  x = 0.0;
  dx = 1.0/(r->final + 1.0);
  r->origin[0] = (real)SRH;
  for(j=1; j <= r->final; j++)
    {
      x += dx;
      r->origin[j] = rcf(x);
      r->origin[-j] = rcf(-x);
    }
  return;
}


/***********************************************************************
 *  rcfmidp()
 *
 *  Sample a rising cutoff function between gridpoints.
 *
 *  Calling sequence and basic algorithm:
 *
 *    rcfmidp( R ):
 *       Let X =  0.5/(R.FINAL+1.0)
 *       Let DX = 1.0/(R.FINAL+1.0)
 *       For J = 0 to R.FINAL
 *          Let R.ORIGIN[J]    = rcf(X)
 *          Let R.ORIGIN[-J-1] = rcf(-X)
 *          X += DX
 *
 *  Input:
 *	(interval *)r	This is a preallocated INTERVAL.
 *
 *  Output:
 *	(interval *)r	This INTERVAL is filled with samples.
 *
 *  Assumptions:
 *	1. r != NULL
 *	2. r->origin != NULL
 *	3. r->final >= 0
 *	4. r->least == -(r->final+1)
 */
extern void
  rcfmidp(			/* Rising cutoff sampled at gridpoints. */
	  interval *r)		/* Put samples in this data structure. */
{
  real x, dx;
  int j;

  assert(r);
  assert(r->origin);
  assert(r->final>=0);
  assert(r->least== -(r->final+1));

  x = 0.5/(r->final + 1.0);
  dx = 1.0/(r->final + 1.0);
  for(j=0; j <= r->final; j++)
    {
      r->origin[j] = rcf(x);
      r->origin[-j-1] = rcf(-x);
      x += dx;
    }
  return;
}
