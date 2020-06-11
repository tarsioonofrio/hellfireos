/*
 * Declare functions which extract graph bases from 1-dimensional
 * wavelet packet and local trigonometric analyses.
 * 
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */


#ifndef GETBASIS_HDR_ALREADY_INCLUDED
# define GETBASIS_HDR_ALREADY_INCLUDED

#include "real.h"
#include "hedge.h"
#include "btn.h"

extern real
  costs2bbasis(			/* Return the best information cost.  */
	       hedge *graph,	/* Empty, to receive the best basis.  */
	       btn *root,	/* Root of the cost-tagged subtree.   */
	       int level);	/* Level index of `root' in the tree. */

extern void
  btnt2costs(			/* Return the best information cost.  */
	     btn *root);	/* Root of the current subtree.   */

extern hedge *
  btnt2bbasis(			/* Return the best graph basis.   */
	      btn *root,	/* Root of the current subtree.   */
	      int maxlevel);	/* Search the tree to this depth. */

extern btn *
  abt2costs(			/* Return a costs-tagged BTN tree.   */
	    real *data,		/* The underlying array binary tree. */
	    int length,		/* Length of one row of `data[]'.    */
	    int maxlevel);	/* This is the depth of `data[]'.    */

extern hedge *
  abt2bbasis(			/* Return the best graph basis.    */
	     real *data,	/* Start of the array binary tree. */
	     int length,	/* Length of one row of the tree.  */
	     int maxlevel);	/* Depth of the array binary tree. */

extern real
  levelcost(			/* Return the cost of a level.     */
	    btn *root,		/* Start of the array binary tree. */
	    int level);		/* Level whose cost is desired.    */

extern real
  costs2blevel(			/* Return the best level basis cost.     */
	       hedge *graph,	/* Empty, to receive the best level.     */
	       btn *root,	/* Root of the cost-tagged subtree.      */
	       int minlevel,	/* Shallowest level of `root' to search. */
	       int maxlevel);	/* Deepest level of `root' to search.    */

extern real
  btnt2bbhedge(			/* Return the best graph basis cost.   */
	       hedge *graph,	/* Partially allocated output struct.  */
	       btn *root,	/* Root of the current subtree.        */
	       int s,		/* Current level in `root'.            */
	       int l);		/* Search the tree to this depth. */

extern hedge *
  btnt2blevel(			/* Return the best graph basis.     */
	      btn *root,	/* Root of the current subtree.     */
	      int minlevel,	/* Search the tree from this depth. */
	      int maxlevel);	/* Search the tree to this depth.   */

extern real
  btntlevelcost(		/* Return the cost of a level.  */
		btn *root,	/* Root of the BTN tree.        */
		int level);	/* Level whose cost is desired. */

extern real
  btnt2blhedge(			/* Return the best graph basis.      */
	       hedge *graph,	/* Partially allocated output hedge. */
	       btn *root,	/* Root of the current subtree.      */
	       int minlevel,	/* Search the tree from this depth.  */
	       int maxlevel);	/* Search the tree to this depth.    */

extern hedge *
  abt2blevel(			/* Return the best level basis.      */
	     real *data,	/* Start of the array binary tree.   */
	     int length,	/* Length of one row of the tree.    */
	     int minlevel,	/* Minimum search depth in the tree. */
	     int maxlevel);	/* Maximum search depth in the tree. */

extern real
  abt2blhedge(			/* Return the best level basis cost. */
	      hedge *graph,	/* Partial HEDGE for the output.     */
	      real *data,	/* Start of the array binary tree.   */
	      int length,	/* Length of one row of the tree.    */
	      int minlevel,	/* Minimum search depth in the tree. */
	      int maxlevel);	/* Maximum search depth in the tree. */

#endif /* GETBASIS_HDR_ALREADY_INCLUDED */
