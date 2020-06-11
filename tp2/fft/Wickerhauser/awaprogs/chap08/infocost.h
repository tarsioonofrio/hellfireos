/*
 *  Declare various information cost functions.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 *
 */

#ifndef INFOCOST_HDR_ALREADY_INCLUDED
# define INFOCOST_HDR_ALREADY_INCLUDED

#include "real.h"

/* Use the following to define a cost function with a parameter: */
/* #ifndef infocost
   # define infocost(u, least, final)	lpnormp(u, least, final, 0.25)
   #endif
*/

/* Define a default cost function: */
#ifndef infocost
# define infocost	ml2logl2
#endif

extern real
  thresh(			/* Threshold cost function.            */
	 const real *u,		/* Preallocated array of coefficients. */
	 int least,		/* Least valid `u[]' index.            */
	 int final,		/* Greatest valid `u[]' index.         */
	 real epsilon);		/* Absolute value threshold.           */


extern real
  l1norm(			/* Absolute summation cost function.   */
	 const real *u,		/* Preallocated array of coefficients. */
	 int least,		/* Least valid `u[]' index.            */
	 int final);		/* Greatest valid `u[]' index.         */

extern real
  lpnormp(			/* p-th power summation cost function. */
	  const real *u,	/* Preallocated array of coefficients. */
	  int least,		/* Least valid `u[]' index.            */
	  int final,		/* Greatest valid `u[]' index.         */
	  double p);		/* Power to use in summation.          */

extern real
  ml2logl2(			/* ``Entropy'' cost function.          */
	  const real *u,	/* Preallocated array of coefficients. */
	  int least,		/* Least valid `u[]' index.            */
	  int final);		/* Greatest valid `u[]' index.         */

extern real
  logl2(			/* Gauss-Markov cost function.         */
	const real *u,		/* Preallocated array of coefficients. */
	int least,		/* Least valid `u[]' index.            */
	int final);		/* Greatest valid `u[]' index.         */

extern real
  tdim(				/* Theoretical dimension.              */
	const real *u,		/* Preallocated array of coefficients. */
	int least,		/* Least valid `u[]' index.            */
	int final);		/* Greatest valid `u[]' index.         */


#endif /* INFOCOST_HDR_ALREADY_INCLUDED */
