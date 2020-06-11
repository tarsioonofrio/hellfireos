
/*
 * This header file declares functions used to convert between the HEDGE,
 * TFA1, array binary tree, and BTN tree data structures.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef UTILITY_HDR_ALREADY_INCLUDED
#define UTILITY_HDR_ALREADY_INCLUDED

#include "real.h"
#include "btn.h"
#include "hedge.h"
#include "interval.h"
#include "tfa.h"

extern void
  abt2hedge( 
	    hedge *graph,	/* Has `levels[]' but not `contents[]'. */
	    real  *data,	/* Source ABT for `graph->contents[]'.  */
	    int    length);	/* Length of one level of `data[]'. */

extern void
  hedge2abt(
	    real  *data,	/* Start of the array binary tree. */
	    hedge *graph,	/* Fully-assigned `hedge' structure. */
	    int    length);	/* Columns per level in `data[]'. */

extern int
  hedge2btnt(
	     btn   *root,	/* NULL or pre-allocated root of a tree */
	     hedge *graph);	/* Pre-allocated, consistent data structure. */

extern int
  btnt2hedge(
	    hedge *graph,	/* Partly-assigned `hedge' structure. */
	    btn   *root);	/* Preassigned `btn' binary tree. */

extern void
  tfa1s2abt(
	    real *data,		/* Preallocated array binary tree */
	    int   n,		/* Length of one row of `data[]'. */
	    tfa1 *atoms,	/* List to superpose into `data[]'. */
	    int   num);		/* Number of elements in `atoms[]'. */

extern int
  tfa1inabt(
	    tfa1 *atom,		/* Test if this fits in the tree. */
	    int   n,		/* Length of one row of the tree. */
	    int   maxlevel);	/* Maximum level in the tree. */

extern int
  tfa1sinabt(
	     tfa1 *atoms,	/* Test if these fit into the tree. */
	     int   num,		/* Test this many TFA1s. */
	     int   length,	/* Length of one row of the tree. */
	     int   maxlevel);	/* Maximum level in the tree. */

extern void
  abt2tfa1(
	   tfa1 *atom,		/* Partially assigned with tags. */
	   real *data,		/* Source of target amplitudes. */
	   int   length);	/* Length of one row of the tree. */

extern void
  abt2tfa1s(
	    tfa1 *atoms,	/* Partially assigned with tags. */
	    int   num,		/* Length of `atoms[]'.*/
	    real *data,		/* Source of target amplitudes. */
	    int   length);	/* Length of one row of the tree. */

extern void
  tfa12btnt(
	    btn  *root,		/* Root of the target tree. */
	    tfa1 *atom);	/* Fully assigned, with tags. */

extern void
  tfa1s2btnt(
	     btn  *root,	/* Root of the target tree. */
	     tfa1 *atoms,	/* Fully assigned, with tags. */
	     int   num);	/* Elements in `atoms[]'. */

extern void
  btnt2tfa1(
	    tfa1 *atom,		/* Partially assigned, with tags. */
	    btn *root);		/* Root of the source tree. */

extern void
  btnt2tfa1s(
	     tfa1 *atoms,	/* Partially assigned, with tags. */
	     int   num,		/* Elements in `atoms[]'. */
	     btn  *root);	/* Root of the source tree. */

extern void
  array2tfa1s(
	      tfa1 *atoms,	/* `num' preallocated TFA1s to assign. */
	      int   num,	/* Elements in `atoms[]', `amplitudes[]' */
	      real *amplitudes,	/* Source amplitudes to assign */
	      int   block,	/* Block index of the assigned TFA1s */
	      int   level);	/* Level index of the assigned TFA1s */

extern void
  hedgeabt2tfa1s(
		 tfa1  *atoms,	/* Preallocated target array. */
		 hedge *graph,	/* Preallocated source HEDGE. */
		 int    length,	/* Number of columns in tree. */
		 int  maxlevel);/* Depth of the underlying tree. */

extern int
  abthedge2tfa1s(
		 tfa1  *atoms,	/* Preallocated target array. */
		 hedge *graph,	/* Preallocated source HEDGE. */
		 int    length);/* Number of columns in tree. */

extern void
  interval2tfa1s(
		 tfa1 *atoms,	/* The preallocated TFA1s to assign. */
		 interval *segment,	/* Source amplitudes to assign */
		 int   block,	/* Block index of the assigned TFA1s */
		 int   level);	/* Level index of the assigned TFA1s */

extern int
  intervalhedge2tfa1s(
		      tfa1  *atoms,	/* Preallocated target array. */
		      hedge *graph);	/* Preallocated source HEDGE. */

extern btn *
  abt2btnt(
	   real *data,		/* Preallocated array binary tree. */
	   int   length,	/* One row of the array binary tree. */
	   int   maxlevel,	/* Depth of the array binary tree. */
	   int   level);	/* `data' points to start of this level. */

#endif /* UTILITY_HDR_ALREADY_INCLUDED */

