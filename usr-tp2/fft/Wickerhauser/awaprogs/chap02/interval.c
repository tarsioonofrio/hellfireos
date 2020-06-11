/*
 * Utilities to manipulate `interval' data structures.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <stdlib.h>		/* for calloc(), malloc() */
#include <assert.h>		/* We abort if malloc() returns NULL. */
#include "real.h"
#include "common.h"
#include "interval.h"


/********************************************************************
 * makeinterval()
 *
 *  Allocate an `interval' and assign its data array.
 *
 *  Calling sequence and basic algorithm:
 *
 *    makeinterval( DATA, LEAST, FINAL ):
 *       Allocate an INTERVAL at SEG with all members 0
 *       Let LENGTH = 1+FINAL-LEAST
 *       If LENGTH>0 then
 *          Allocate an array of LENGTH real's at SEG.ORIGIN
 *          Shift SEG.ORIGIN -= LEAST
 *          If DATA != NULL then
 *             For K = LEAST to FINAL
 *                Let SEG.ORIGIN[K] = DATA[K-LEAST]
 *       Let SEG.LEAST = LEAST
 *       Let SEG.FINAL = FINAL
 *       Return SEG
 *
 *
 *  Inputs:
 *	(const real *)data	If this is non-NULL, we copy its contents
 *			into the newly allocated `origin' array, with a
 *			shift so that `origin[least]==data[0]'
 *
 *	(int)least	This must be a valid index for non-NULL `data[]'.
 *
 *	(int)final	This must be a valid index for non-NULL `data[]'.
 *
 * Outputs:
 *	(interval *)makeinterval	The  return value is a pointer
 *			to a newly-allocated `interval' data structure
 *			with `data[0]', ..., `data[final-least]' copied into
 *			its newly-allocated `origin' member at indices
 *			`origin[least]', ..., `origin[final]', or else
 *			with an array of zeroes with valid indices `least'
 *			through `final'.   If `least > final', any `data[]'
 *			array is ignored and NULL is assigned into the
 *			`interval's `origin' array.
 *
 * External functions called:
 *	calloc(),malloc()	Declared in <stdlib.h>
 *
 * Assumption:
 *	If `data != NULL' then `data[n]' is a valid array element for
 *		each `0 <= n <= final-least'.
 */
extern interval *
  makeinterval( 
	       const real *data,/* NULL, or a valid array starting at 0. */
	       int least,	/* Least valid index of output `origin'. */
	       int final )	/* Final valid index of output `origin'. */
{
  interval *seg;
  int length, k;

  seg = (interval *)calloc(1, sizeof(interval));  assert(seg);
  length = 1 + final - least;
  if(length>0)
    {
      seg->origin = (real *)calloc(length, sizeof(real)); assert(seg->origin);
      seg->origin -= least;
      if( data )
	for( k=least; k<=final; k++ )
	  seg->origin[k] = *data++;
    }	
  seg->least = least;
  seg->final = final;
  return(seg);
}

/********************************************************************
 * freeinterval()
 *
 *  Deallocate an `interval' and its data array.
 *
 *  Calling sequence:
 *	freeinterval( seg )
 *
 *  Basic algorithm:
 *
 *   If SEG != NULL then
 *      If SEG.ORIGIN != NULL then
 *         Shift SEG.ORIGIN += LEAST
 *         Free SEG.ORIGIN
 *      Free SEG
 *   Return NULL
 *
 *
 *  Inputs:
 *	(interval *)seg      If non-NULL, we deallocate it and any non-NULL
 *				data array at `seg->origin'.
 *
 * Outputs:
 *	(interval *)freeinterval	The  return value is always NULL.
 *
 * External functions called:
 *	free()			Declared in <stdlib.h>
 *
 * Assumption:
 *	If `seg != NULL' and `seg->origin+seg->least != NULL', then 
 *	`seg->origin+seg->least' points to the start of a previously-
 *	allocated array.
 */
extern interval *
  freeinterval( 
	       interval *seg)	/* `interval' data structure to free. */
{
  if( seg )
    {
      if( seg->origin )
	free(seg->origin + seg->least);
      free( seg );
    }
  return(NULL);      
}

/********************************************************************
 * enlargeinterval()
 *
 *  Enlarge an `interval', preserving any contents in its
 *	data array.
 *
 *  Calling sequence:
 *	enlargeinterval( old, least, final )
 *
 *  Basic algorithm:
 *
 *  If OLD.ORIGIN == NULL then
 *     Let LENGTH = 1+FINAL-LEAST
 *     If LENGTH>0 then
 *        Allocate an array of LENGTH real's at OLD.ORIGIN
 *        Shift OLD.ORIGIN -= LEAST
 *        Let OLD.LEAST = LEAST
 *        Let OLD.FINAL = FINAL
 *  Else
 *     If OLD.LEAST<LEAST || OLD.FINAL>FINAL then
 *        Let LEAST = min( OLD.LEAST, LEAST )
 *        Let FINAL = max( OLD.FINAL, FINAL )
 *        Let LENGTH = 1 + FINAL-LEAST
 *        If LENGTH>0 then
 *           Allocate an array of LENGTH real's at NEWDATA
 *           Shift NEWDATA -= LEAST
 *           For J = OLD.LEAST to OLD.FINAL
 *              Let NEWDATA[J] = OLD.ORIGIN[J]
 *           Shift OLD.ORIGIN += OLD.LEAST
 *           Deallocate OLD.ORIGIN[]
 *           Let OLD.ORIGIN = NEWDATA
 *           Let OLD.LEAST = LEAST
 *           Let OLD.FINAL = FINAL
 *  Return OLD
 *
 *
 *  Inputs:
 *	(interval *)old	  If `old->origin' is NULL or `least < old->least'
 *				or `final > old->final', we allocate a new
 *				array, copy `old->origin[]'s contents into
 *				it, shift it and assign it to `old->origin'.
 *
 *	(int)least	This is the new least valid index.
 *
 *	(int)final	This is the new final valid index.
 *
 * Outputs:
 *	(interval *)enlargeinterval	The  return value is a pointer
 *			to `old', but with its `origin' member possibly
 *			pointing to a longer array.
 *
 * External functions called:
 *	calloc(),free()		Declared in <stdlib.h>
 *	makeinterval()
 *
 * Assumptions:
 *	We must be allowed to deallocate `old->origin[]' if it is too short.
 */
extern interval *
  enlargeinterval( 
		  interval *old, /* NULL, or a valid data array. */
		  int least,	/* New least `origin[]' index.   */
		  int final )	/* New final `origin[]' index.   */
{
  int i, length;
  real *newdata;

  if( old )			/* Fix the `old' interval */
    {
      if( old->origin )		/* Non-NULL data array. */
	{
	  if( least < old->least || final > old->final )
	    {
	      least = min( least, old->least );
	      final = max( final, old->final );
	      length = 1+final-least;
	      if( length > 0 )
		{
		  newdata = (real *)calloc(length, sizeof(real));
		  assert(newdata);
		  newdata -= least;
		  for( i=old->least; i<=old->final; i++)
		    newdata[i] = old->origin[i];
		  free( old->origin - old->least );
		  old->origin = newdata;
		  old->least = least;
		  old->final = final;
		}
	      /* ...otherwise do nothing. */
	    }
	  /* ...otherwise do nothing. */
	}
      else			/* NULL data array. */
	{
	  length = 1+final-least;
	  if( length > 0 )
	    {
	      newdata = (real *)calloc(length, sizeof(real));
	      assert(newdata);
	      old->least = least;
	      old->final = final;
	      old->origin = newdata - least;
	    }
	  /* ...otherwise do nothing. */
	}
    }
  else				/* NULL input interval. */
    {
      old = makeinterval( 0, least, final );
    }
  return(old);
}

/********************************************************************
 * ininterval()
 *
 *  Check if an offset is in an interval.
 *
 *  Calling sequence:
 *	ininterval( segment, offset )
 *
 *  Basic algorithm:
 *
 *   If SEGMENT.ORIGIN != NULL then
 *      If OFFSET>=SEGMENT.LEAST && OFFSET<=SEGMENT.FINAL then
 *         Return TRUE
 *   Return FALSE
 *
 *
 *  Inputs:
 *	(interval *)segment   This must a valid, preallocated `interval'.
 *
 *	(int)offset	      This can be any index.
 *
 * Outputs:
 *	(int)ininterval	      The  return value is nonzero if 
 *				`segment->origin[offset]' is a valid
 *				array element.  Otherwise, it is 0.
 *
 * Assumption:
 *	If `segment != NULL' and `segment->origin != NULL', then 
 *	`segment->least' through `segment->final' are valid indices in the
 *      array `segment->origin[]'.
 */
extern int
  ininterval( 
	     interval *segment,	/* `interval' data structure to test. */
	     int offset)	/* Index to validate. */
{
  if( segment )
    if( segment->origin )
      if( offset>=segment->least && offset<=segment->final)
	return(TRUE);
  return(FALSE);      
}

/********************************************************************
 * intervalstotal()
 *
 *  Sum the lengths in a list of nonoverlapping intervals.
 *
 *  Calling sequence:
 *	intervalstotal( in, num )
 *
 *  Basic algorithm:
 *
 *   Let TOTAL = 0
 *   For K = 0 to NUM-1
 *      TOTAL += 1 + IN[K].FINAL - IN[K].LEAST
 *   Return TOTAL
 *
 *
 *  Inputs:
 *	(interval *)in   This must be an array of preallocated `interval's.
 *
 *	(int)n	      `in[0]' through `in[n-1] must be `interval's.
 *
 * Outputs:
 *	(int)intervalstotal	      The  return value is the total
 *				number of storage elements needed to
 *			       hold all the data of all the intervals.
 *
 */
extern int
  intervalstotal( 
		 interval *in,	/* Array of `interval' data structures. */
		 int num)	/* Number of `interval's in `in[]'. */
{
  int total,i;

  total = 0;
  if( in )
    for(i=0; i<num; i++)
      total += 1 + in[i].final - in[i].least;
  return(total); 
}
