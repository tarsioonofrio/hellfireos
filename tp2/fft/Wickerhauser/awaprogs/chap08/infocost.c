/*
 *
 * These are various kinds of information cost functions useable
 * for best-basis searches.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <assert.h>
#include <math.h>
#include "real.h"
#include "common.h"
#include "infocost.h"

/*********************************************************************
 * thresh()
 *
 * Threshold information cost function.
 *
 * Calling sequence and basic algorithm:
 *   thresh( U, LEAST, FINAL, EPSILON ):
 *      Let COST = 0
 *      For K = LEAST to FINAL
 *	  If absval(U[K]) > EPSILON then
 *	    COST += 1
 *      Return COST
 *
 * Inputs:
 *	(const real *)u		This must point to a preallocated array.
 *
 *	(int)least		These are the least and final valid
 *	(int)final		  indices in `u[]'.
 *
 *	(real)epsilon		This is a positive, absolute threshold.
 *
 * Output:
 *	(real)thresh		The return value is the number of elements
 *				  of `u[]' whose absolute values are
 *				  greater than `epsilon'.
 *
 *  Assumptions:
 *	1. epsilon>=0
 *
 * External functions called:
 *	absval()
 */
extern real
  thresh(			/* Threshold cost function.            */
	 const real *u,		/* Preallocated array of coefficients. */
	 int least,		/* Least valid `u[]' index.            */
	 int final,		/* Greatest valid `u[]' index.         */
	 real epsilon)		/* Absolute value threshold.           */
{
  real cost;
  int k;

  assert(epsilon>=0);

  cost = 0;
  for( k=least; k<=final; k++ )
    if( absval(u[k])>epsilon ) cost += 1.0;
  return(cost);
}

/*********************************************************************
 * l1norm()
 *
 * Absolute summation information cost function.
 *
 * Calling sequence and basic algorithm:
 *   l1norm( U, LEAST, FINAL ):
 *      Let COST = 0
 *      For K = LEAST to FINAL
 *	   COST += absval( U[K] )
 *      Return COST
 *
 * Inputs:
 *	(const real *)u		This must point to a preallocated array.
 *
 *	(int)least		These are the least and final valid
 *	(int)final		  indices in `u[]'.
 *
 * Output:
 *	(real)l1norm		The return value is the sum of the absolute
 *				  values of the elements of `u[]'.
 *
 * External functions called:
 *	absval()
 */
extern real
  l1norm(			/* Absolute summation cost function.   */
	 const real *u,		/* Preallocated array of coefficients. */
	 int least,		/* Least valid `u[]' index.            */
	 int final)		/* Greatest valid `u[]' index.         */
{
  real cost;
  int k;

  cost = 0;
  for( k=least; k<=final; k++ )  cost += absval(u[k]);
  return(cost);
}

/*********************************************************************
 * lpnormp()
 *
 * p-th power summation information cost function.
 *
 * Calling sequence and basic algorithm:
 *   lpnormp( U, LEAST, FINAL, P ):
 *      Let COST = 0
 *      For K = LEAST to FINAL
 *         Let ABSU = absval( U[K] )
 *         If ABSU > 0 then
 *           COST += exp( P*log( ABSU ) )
 *      Return COST
 *
 * Inputs:
 *	(const real *)u		This must point to a preallocated array.
 *
 *	(int)least		These are the least and final valid
 *	(int)final		  indices in `u[]'.
 *
 *	(double)p		Use this power in the summation.
 *
 * Output:
 *	(real)lpnormp		The return value is the sum of the `p'-th
 *				 power of the absolute values of the
 *				 elements of `u[]'.
 *
 * External functions called:
 *	absval(), exp(), log()
 */
extern real
  lpnormp(			/* p-th power summation cost function. */
	  const real *u,	/* Preallocated array of coefficients. */
	  int least,		/* Least valid `u[]' index.            */
	  int final,		/* Greatest valid `u[]' index.         */
	  double p)		/* Power to use in summation.          */
{
  double cost, absu;
  int k;

  cost = 0;
  for( k=least; k<=final; k++ )
    {
      absu = (double)absval(u[k]);
      if( absu>0 ) cost += exp( p * log(absu) );
    }
  return((real)cost);
}

/*********************************************************************
 * ml2logl2()
 *
 * ``Entropy'' information cost function.
 *
 * Calling sequence and basic algorithm:
 *   ml2logl2( U, LEAST, FINAL ):
 *      Let COST = 0
 *      For K = LEAST to FINAL
 *         Let USQ =  U[K] * U[K]
 *         If USQ > 0 then
 *            COST -= USQ * log( USQ )
 *      Return COST
 *
 * Inputs:
 *	(const real *)u		This must point to a preallocated array.
 *
 *	(int)least		These are the least and final valid
 *	(int)final		  indices in `u[]'.
 *
 * Output:
 *	(real)ml2logl2		The return value is the sum of `u*u*log(u*u)'
 *				 for all nonzero values of `u*u'.
 *
 * External functions called:
 *	 log()
 */
extern real
  ml2logl2(			/* ``Entropy'' cost function.          */
	  const real *u,	/* Preallocated array of coefficients. */
	  int least,		/* Least valid `u[]' index.            */
	  int final)		/* Greatest valid `u[]' index.         */
{
  double cost, usq;
  int k;

  cost = 0;
  for( k=least; k<=final; k++ )
    {
      usq = (double)(u[k]*u[k]);
      if( usq>0 ) cost -= usq * log(usq);
    }
  return((real)cost);
}

/*********************************************************************
 * logl2()
 *
 * Gauss-Markov information cost function.
 *
 * Calling sequence and basic algorithm:
 *   logl2( U, LEAST, FINAL ):
 *      Let COST = 0
 *      For K = LEAST to FINAL
 *         COST += log( U[K]*U[K] )
 *      Return COST
 *
 * Inputs:
 *	(const real *)u		This must point to a preallocated array.
 *
 *	(int)least		These are the least and final valid
 *	(int)final		  indices in `u[]'.
 *
 * Output:
 *	(real)logl2		The return value is the sum of `log(u*u)'
 *				 for all values of `u*u'.
 *
 * External functions called:
 *	log()
 */
extern real
  logl2(			/* Gauss-Markov cost function.         */
	const real *u,		/* Preallocated array of coefficients. */
	int least,		/* Least valid `u[]' index.            */
	int final)		/* Greatest valid `u[]' index.         */
{
  double cost;
  int k;

  cost = 0;
  for( k=least; k<=final; k++ )
    cost += log((double)(u[k]*u[k]));
  return((real)cost);
}

/*********************************************************************
 * tdim()
 *
 * Theoretical dimension of a sequence
 *
 * Calling sequence and basic algorithm:
 *   tdim( U, LEAST, FINAL ):
 *      Let MU2LOGU2 = 0
 *      Let ENERGY = 0
 *      For K = LEAST to FINAL
 *         Let USQ =  U[K] * U[K]
 *         If USQ > 0 then
 *            MU2LOGU2 -= USQ * log( USQ )
 *            ENERGY += USQ
 *      If ENERGY > 0 then
 *         Let THEODIM =  ENERGY * exp( MU2LOGU2 / ENERGY )
 *      Else
 *         Let THEODIM = 1
 *      Return THEODIM
 *
 * Inputs:
 *	(const real *)u		This must point to a preallocated array.
 *
 *	(int)least		These are the least and final valid
 *	(int)final		  indices in `u[]'.
 *
 * Output:
 *	(real)tdim		The return value is exp of th sum of
 *				 `u*u*log(u*u)' for all nonzero values
 *				 of `u*u'.
 *
 * External functions called:
 *	log(), exp()
 */
extern real
  tdim(				/* Theoretical dimension.              */
	const real *u,		/* Preallocated array of coefficients. */
	int least,		/* Least valid `u[]' index.            */
	int final)		/* Greatest valid `u[]' index.         */
{
  double theodim, mu2logu2, energy, usq;
  int k;

  mu2logu2 = 0;
  energy = 0;
  for( k=least; k<=final; k++ )
    {
      usq =  (double)( u[k] * u[k] );
      if( usq>0 )
	{
	  mu2logu2 -= usq * log( usq );
	  energy += usq;
	}
    }
   if(energy>0)
     theodim =  energy * exp( mu2logu2 / energy );
   else
     theodim = 0.0;
  return(theodim);
}

