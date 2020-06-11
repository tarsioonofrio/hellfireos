/*
 * This header file declares the  `btn' data type and the functions
 *  which allocate and deallocate structures of `btn' nodes.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef BTN_HDR_ALREADY_INCLUDED
#define BTN_HDR_ALREADY_INCLUDED

#include "real.h"
#include "fntype.h"

typedef struct btn {
  void       *content;		/* Interval or whatnot. */
  struct btn *left;		/* Pointer to left descendent. */
  struct btn *right;		/* Pointer to right descendent. */
  void       *tag;		/* Data structure with more info. */
}
btn;				/* [B]inary [T]ree [N]ode */


extern btn *
  makebtn(
	  void *content,	/* Content data structure, or NULL. */
	  btn  *left,		/* Left descendent node, or NULL. */
	  btn  *right,		/* Right descendent node, or NULL. */
	  void *tag);		/* Extra information, or NULL. */

extern btn *
  freebtn(
	  btn  *node,		/* Previously-allocated BTN structure. */
	  freetype freecontent, /* Used to free the content member. */
	  freetype freetag);	/* Used to free the tag member.     */

extern btn *
  makebtnt(
	   int level);		/* Number of levels in the tree. */

extern btn *
  freebtnt(
	   btn  *root,		/* Root of pre-allocated BTN tree. */
	   freetype freecontent, /* Used to free the content member. */
	   freetype freetag);	/* Used to free the tag member.      */

extern btn *
  freebtns(
	   btn  *root,		/* Root of pre-allocated BTN tree. */
	   freetype freetag);	/* Used to free the tag member.      */

extern btn *
  btn2branch(
	     btn  *self,	/* Current root of the BTN tree. */
	     int  level,	/* Target level with respect to `self'. */
	     int  block);	/* Block index with respect to `self'. */
extern btn *
  btnt2btn(
	   btn  *root,		/* Current root of the BTN tree. */
	   int  level,		/* Target level with respect to `root'. */
	   int  block);		/* Block index with respect to `root'. */


#endif /* BTN_HDR_ALREADY_INCLUDED */

