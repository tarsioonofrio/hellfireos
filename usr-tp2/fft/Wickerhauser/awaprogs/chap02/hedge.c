/*
 * Utilities to manipulate `hedge' data structures.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 *   by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <stdlib.h>		/* For calloc(), malloc(), and free(). */
#include <assert.h>		/* We abort if malloc() returns NULL. */
#include "hedge.h"


/********************************************************************
 * makehedge()
 *
 *	  Allocate a `hedge' and assign it a data array.
 *
 *  Calling sequence and basic algorithm:
 *
 *	makehedge( BLOCKS, CONTENTS, LEVELS, TAG ):
 *	   Allocate a HEDGE data structure at OUT
 *	   Let OUT.BLOCKS = BLOCKS
 *	   If CONTENTS==NULL then
 *	      Allocate an array of BLOCKS pointers at OUT.CONTENTS
 *	   Else
 *	      Let OUT.CONTENTS = CONTENTS
 *	   If LEVELS==NULL then
 *	      Allocate an array of BLOCKS bytes at OUT.LEVELS
 *	   Else
 *	      Let OUT.LEVELS = LEVELS
 *	   Let OUT.TAG = TAG
 *	   Return OUT
 *
 *
 *  Inputs:
 *	(int)blocks	This is the number of elements in the levels array,
 *			 and the number of elements in the contents array,
 *			 or else 0.
 *
 *	(void *)contents  This should be a valid array of data structures
 *			 with `blocks' elements, or else NULL.
 *
 *	(unsigned char *)levels	  This must be a valid array of `blocks'
 *				levels in encounter order, or else NULL
 *
 *	(void *)tag	This should point to a data structure, or
 *			 else be NULL.
 *
 *
 * Outputs:
 *	(hedge *)makehedge	The  return value is a pointer
 *			to a newly-allocated `hedge' data structure
 *			with members `contents', `levels' and `tag', or
 *
 *
 * External functions called:
 *	calloc()	Declared in <stdlib.h>
 *
 * Assumption:
 *	If `data != NULL' then `data[n]' is a valid array element for
 *		each `least <= n <= final'.
 */
extern hedge *
  makehedge( 
	    int blocks,		/* Length of `contents[]' and `levels[]'. */
	    void **contents,	/* Array of length `blocks', or NULL.  */
	    unsigned char *levels, /* Array in encounter order, or NULL. */
	    void *tag)		/* Data structure pointer, or NULL. */
{
  hedge *out;

  out = (hedge *)calloc(1, sizeof(hedge)); assert(out);
  out->blocks = blocks;
  if( blocks>0 )
    {
      if( contents==0 )
	{
	  out->contents = (void **)calloc(blocks, sizeof(void *));
	  assert(out->contents);
	}
      else
	out->contents = contents;

      if( levels==0 )
	{
	  out->levels = (unsigned char *)calloc(blocks, sizeof(unsigned char));
	  assert(out->contents);
	}
      else
	out->levels = levels;
      
      out->tag = tag;
    }
  return(out);
}

/********************************************************************
 * freehedge()
 *
 *	  Deallocate a `hedge' and its arrays.
 *
 *  Calling sequence and basic algorithm:
 *
 *	freehedge( IN, FREECONTENT, FREETAG ):
 *	   If IN.CONTENTS != NULL then
 *            For B = 0 to IN.BLOCKS-1
 *               Deallocate IN.CONTENTS[B] with FREECONTENT()
 *	      Deallocate IN.CONTENTS
 *	   If IN.LEVELS != NULL then
 *	      Deallocate IN.LEVELS
 *	   If IN.TAG != NULL then
 *	      Deallocate IN.TAG with FREETAG()
 *	   Deallocate IN
 *	   Return NULL
 *
 *  Inputs:
 *	(interval *)in      If non-NULL, we deallocate it and any non-NULL
 *			      members at `in->contents', `in->levels', or
 *			      `in->tag'.
 *
 *	(freetype)freecontent	This points to a function which deallocates
 *				  members of the contents array.
 *
 *	(freetype)freetag	This points to a function which deallocates
 *				  the tag members of the input hedge.
 *
 * Outputs:
 *	(hedge *)freehedge	The  return value is always NULL.
 *
 * External functions called:
 *	free()			Declared in <stdlib.h>
 *
 */
extern hedge *
  freehedge( 
	    hedge *in,		/* `hedge' data structure to free.  */
	    freetype freecontent, /* Used to free contents members. */
	    freetype freetag)	/* Used to free the tag member.     */
{
  if( in )
    {
      if( in->contents )
	{
	  int j;
	  for(j=0; j<in->blocks; j++ ) 
	    in->contents[j] = (*freecontent)(in->contents[j]);
	  free(in->contents);
	}
      if( in->levels ) free(in->levels);
      if( in->tag )
	in->tag = (*freetag)(in->tag);
      free( in );
    }
  return(0);      
}
