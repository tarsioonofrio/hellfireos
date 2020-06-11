
/*
 * This header file declares the `hedge' data type and some functions
 * which manipulate it.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef HEDGE_HDR_ALREADY_INCLUDED
# define HEDGE_HDR_ALREADY_INCLUDED

#include "fntype.h"		/* for `freetype' type definition. */

typedef struct {
  int   blocks;			/* # of `contents[]' and `levels[] elements. */
  void **contents;		/* Start of array of `real's or `void *'.  */
  unsigned char *levels;	/* Levels list in encounter order. */
  void *tag;			/* Information data structure. */
} hedge;

extern hedge *
  makehedge( 
	    int blocks,		/* Length of `contents[]' and `levels[]'. */
	    void **contents,	/* Array of length `blocks', or NULL.  */
	    unsigned char *levels, /* Array in encounter order, or NULL. */
	    void *tag);		/* Data structure pointer, or NULL. */

extern hedge *
  freehedge( 
	    hedge *in,		/* `hedge' data structure to free.  */
	    freetype freecontent, /* Used to free contents members. */
	    freetype freetag);	/* Used to free the tag member.     */

#endif /* HEDGE_HDR_ALREADY_INCLUDED */

