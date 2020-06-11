/*
 * Utilities to manipulate `interval', `hedge', `btn', and `tfa1' data
 * structures.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <stdlib.h>		/* For calloc(), malloc(), and free(). */
#include <assert.h>		/* We abort if malloc() returns NULL. */
#include "real.h"
#include "abt.h"
#include "btn.h"
#include "common.h"		/* For TRUE,FALSE and min(),max() */
#include "hedge.h"
#include "interval.h"
#include "tfa.h"
#include "utility.h"

/********************************************************************
 * abt2hedge()
 *
 *  Complete a levels-specified `hedge' from an array binary tree.
 *
 *  Calling sequence and basic algorithm:
 *
 *   abt2hedge( GRAPH, DATA, LENGTH ):
 *     Let COLUMN = 0
 *     For I = 0 to GRAPH.BLOCKS-1
 *        Let GRAPH.CONTENTS[I] = DATA+COLUMN+LENGTH*GRAPH.LEVELS[I]
 *        COLUMN += LENGTH>>GRAPH.LEVELS[I]
 *
 *
 *  Inputs:
 *	(hedge *)graph	This points to a data structure. `graph->levels[]'
 *			 must be an array of `graph->blocks' numbers, and
 *			 `graph->contents[]' must be an array of 
 *			 `graph->blocks' pointers.
 *
 *	(real *)data	This must point to the first element of an array
 *			 binary tree with `length' columns per level.
 *
 *	(int)length	This is the length of one level of the tree.
 *
 * Outputs:
 *	The array `graph->contents[]' is filled with pointers into `data[]'
 *
 * Assumption:
 *	There are enough levels in the array binary tree as will be
 *	requested by the elements of `graph->levels[]'.
 */
extern void
  abt2hedge( 
	    hedge *graph,	/* Has `levels[]' but not `contents[]'. */
	    real  *data,	/* Source ABT for `graph->contents[]'.  */
	    int    length)	/* Length of one level of `data[]'. */
{
  int i, column;

  column = 0;
  for( i=0; i<graph->blocks; i++ )
    {
      graph->contents[i] = (void *)(data + column+length*graph->levels[i]);
      column += length >> graph->levels[i];
    }
}

/********************************************************************
 * hedge2abt()
 *
 *   Superpose amplitudes from a hedge into an array binary tree.
 *
 *  Calling sequence and basic algorithm:
 *
 *   hedge2abt( DATA, GRAPH, LENGTH ):
 *     Let COLUMN = 0
 *     For I = 0 to GRAPH.BLOCKS-1
 *        Let BLENGTH = LENGTH>>GRAPH.LEVELS[I]
 *        Let BLOCK = DATA + LENGTH*GRAPH.LEVELS[I] + COLUMN
 *        For J = 0 to BLENGTH-1
 *           BLOCK[J] += GRAPH.CONTENTS[J]
 *        COLUMN += BLENGTH
 *
 *
 *  Inputs:
 *	(real *)data	This points to the first element of an array
 *			 binary tree.
 *
 *	(hedge *)graph	This must be a fully-assigned `hedge' with 
 *			 pointers to `real' as its contents.
 *
 *	(int)length	This must be a positive integer.
 *
 * Outputs:
 *	Amplitudes from `graph->contents[]' are superposed into the 
 *	array binary tree at `data[]' by side effect.
 *
 * Assumptions:
 *	The array binary tree at `data[]' has `length' columns per
 *	level and at least as many levels as the maximum element
 *	of `graph->levels[]'.
 *
 */
extern void
  hedge2abt(
	    real  *data,	/* Start of the array binary tree. */
	    hedge *graph,	/* Fully-assigned `hedge' structure. */
	    int    length)	/* Columns per level in `data[]'. */
{
  int column, i, j, blength;
  real *block, *content;

  column = 0;
  for( i=0; i<graph->blocks; i++ )
    {
      content = (real *)(graph->contents[i]);
      block = data + length*graph->levels[i] + column;
      blength = length >> graph->levels[i];
      for( j=0; j<blength; j++ )
	block[j] += content[j];
      column += blength;
    }
}

/********************************************************************
 * btnt2hedge()
 *
 *	Fill a levels-specified hedge from a `btn' tree.
 *
 *  Calling sequence and basic algorithm:
 *
 *	btnt2hedge( GRAPH, ROOT ):
 *	   Let MAXLEVEL = 0
 *	   Let FRACTION = 0
 *	   For I = 0 to GRAPH.BLOCKS-1
 *	      Let LEVEL = GRAPH.LEVELS[I]
 *	      If LEVEL>MAXLEVEL then
 *	         FRACTION <<= LEVEL-MAXLEVEL
 *	         Let MAXLEVEL = LEVEL
 *	         Let BLOCK = FRACTION
 *	      Else
 *	         Let BLOCK = FRACTION>>(MAXLEVEL-LEVEL)
 *	      Let NODE = btnt2btn( ROOT, LEVEL, BLOCK )
 *	      Let GRAPH.CONTENTS[I] = NODE.CONTENT
 *	      FRACTION += 1<<LEVEL
 *	   Return MAXLEVEL
 *
 *
 *  Inputs:
 *	(hedge *)graph	This is a partially-allocated `hedge'.
 *
 *	(btn *)root	This is the root of a pre-allocated `btn'
 *			 tree deep enough to fill `graph'.
 *
 * Outputs:
 *	(int)btnt2hedge	The return value is the maximum level in
 *			 the levels member of `graph'.
 *			 
 *	Content members from nodes in the tree at `root' which are
 *	 specified by the `graph->levels[]' are assigned 
 *	array binary tree at `data[]' by side effect.
 *
 * Assumptions:
 *	1.   `root' is non-NULL
 *	2.   `graph' is non-NULL
 *	3.   `graph->contents[]' has `graph->blocks' elements.
 *	4.   `graph->levels[]' has `graph->blocks' elements.
 *
 */
extern int
  btnt2hedge(			/* Return the depth of `graph->levels'. */
	    hedge *graph,	/* Partly-assigned `hedge' structure. */
	    btn   *root)	/* Preassigned `btn' binary tree. */
{
  int  maxlevel, fraction, i, level, block;
  btn *node;

  assert(graph);
  assert(root);

  maxlevel = 0;
  fraction = 0;
  for(i=0; i<graph->blocks; i++)
    {
      level = graph->levels[i];
      if (level>maxlevel)
	{
	  fraction <<= (level-maxlevel);
	  maxlevel = level;
	  block = fraction;
	}
      else
	block = fraction >> (maxlevel-level);
      node = btnt2btn( root, level, block );
      graph->contents[i] = node->content;
      fraction += 1<<level;
    }
  return(maxlevel);
}

/******************************************************************
 * hedge2btnt()
 *
 *	Make a partial BTN tree from a hedge
 *
 *  Calling sequence and basic algorithm:
 *
 *	hedge2btnt( ROOT, GRAPH ):
 *	   Let MAXLEVEL = 0
 *	   Let FRACTION = 0
 *	   For I = 0 to GRAPH.BLOCKS-1
 *	      Let LEVEL = GRAPH.LEVELS[I]
 *	      If LEVEL>MAXLEVEL then
 *	         FRACTION <<= LEVEL-MAXLEVEL
 *	         Let MAXLEVEL = LEVEL
 *	         Let BLOCK = FRACTION
 *	      Else
 *	         Let BLOCK = FRACTION>>(MAXLEVEL-LEVEL)
 *	      Let NODE = btn2branch( ROOT, LEVEL, BLOCK )
 *	      Let NODE.CONTENT = GRAPH.CONTENTS[I]
 *	      FRACTION += 1<<LEVEL
 *	   Return MAXLEVEL
 *
 *  Inputs:
 *	(btn *)root	This must be either NULL or preallocated.
 *
 *	(hedge *)graph	This must be non-NULL.
 *
 *  Outputs:
 *	(int)hedge2btnt	The returned value is the maximum level
 *			prduced in the tree.
 *
 *  External functions called:
 *	btn2branch()
 *
 *  Assumptions:
 *	1.  `graph->levels[]' has `graph->blocks' members.
 *	2.  `graph->contents[]' has `graph->blocks' members.
 *	3.  `graph->levels[]' has only nonnegative elements.
 */
extern int
  hedge2btnt(
	     btn   *root,	/* NULL or pre-allocated root of a tree */
	     hedge *graph)	/* Pre-allocated, consistent data structure. */
{
  int maxlevel, i, fraction, block, level;
  btn  *node;

  assert(graph);

  maxlevel = 0;
  fraction = 0;
  for(i=0; i<graph->blocks; i++)
    {
      level = graph->levels[i];
      assert(level>=0);
      if(level>maxlevel)
	{
	  fraction <<= (level-maxlevel);
	  maxlevel = level;
	  block = fraction;
	}
      else
	block = fraction >>(maxlevel-level);
      node = btn2branch( root, level, block );
      node->content = graph->contents[i];
      fraction += 1<<level;
    }
  return(maxlevel);
}

/****************************************************************
 * tfa1s2abt()
 *
 *	Superpose TFA1's into an array binary tree.
 *
 *  Calling sequence and basic algorithm:
 *
 *	tfa1s2abt( DATA, N, ATOMS, NUM ):
 *	   For K = 0 to NUM-1
 *	      Let START = abtblock( N, ATOMS[K].LEVEL, ATOMS[K].BLOCK )
 *	      DATA[START + ATOMS[K].OFFSET] += ATOMS[K].AMPLITUDE
 *
 *  Inputs:
 *	(real *)data	This array should be pre-allocated.
 *
 *	(int)n		This is the positive length of one row of `data[]'.
 *
 *	(taf1 *)atoms	This is a preallocated array.
 *
 *	(int)num	This the nonnegative length of `atoms[]'.
 *
 *  Outputs:
 *	Amplitudes from `atoms' are superposed into the array binary tree.
 *
 *  Assumptions:
 *	1.  `data != NULL'
 *	2.  `data[]' is long enough to hold the amplitudes from `atoms[]'.
 */
extern void
  tfa1s2abt(
	    real *data,		/* Preallocated array binary tree */
	    int   n,		/* Length of one row of `data[]'. */
	    tfa1 *atoms,	/* List to superpose into `data[]'. */
	    int   num)		/* Number of elements in `atoms[]'. */
{
  int k;

  assert( data );
  assert( atoms || !num );
 
  for(k=0; k<num; k++)
    {
      data[ abtblock(n,atoms[k].level,atoms[k].block) + atoms[k].offset ]
	+= atoms[k].amplitude;
    }
}

/****************************************************************
 * tfa1inabt()
 *
 *	Verify that a TFA1 fits into an array binary tree.
 *
 *  Calling sequence and basic algorithm:
 *
 *	tfa1inabt( ATOM, N, MAXLEVEL ):
 *	   If ATOM.LEVEL>=0 && ATOM.LEVEL<=MAXLEVEL then
 *	     If ATOM.BLOCK>=0 && ATOM.BLOCK<(1<<ATOM.LEVEL) then
 *	       If ATOM.OFFSET>=0 && ATOM.OFFSET<(N>>ATOM.LEVEL) then
 *	          Return TRUE
 *	   Return FALSE
 *
 *  Inputs:
 *	(tfa1 *)atom	This atom should have all its members assigned.
 *
 *	(int)n		This is the positive length of one row of the tree.
 *
 *	(int)maxlevel	This is the nonnegative maximum level in the tree.
 *
 *  Outputs:
 *	(int)tfa1inabt	The return value is TRUE (nonzero) if it fits
 *			or FALSE (zero) if it doesn't.
 */
extern int
  tfa1inabt(
	    tfa1 *atom,		/* Test if this fits in the tree. */
	    int   n,		/* Length of one row of the tree. */
	    int   maxlevel)	/* Maximum level in the tree. */
{
  if( atom->level >= 0 && atom->level <= maxlevel )
    if( atom->block >= 0 && atom->block < (1<<atom->level) )
      if( atom->offset >= 0 && atom->offset < abtblength(n,atom->level) )
	return(TRUE);
  return(FALSE);
}

/****************************************************************
 * tfa1sinabt()
 *
 *	Verify that a list of TFA1's fits into an array binary tree.
 *
 *  Calling sequence and basic algorithm:
 *
 *	tfa1sinabt( atoms, num, lengthx, maxlevel )
 *	   For K = 0 to NUM-1
 *	      If !tfa1inabt( ATOMS[K], LENGTH, MAXLEVEL ) then
 *		 Return FALSE
 *	   Return TRUE
 *
 *  Inputs:
 *	(tfa1 *)atoms	This should be a pre-assigned array of TFA1s.
 *
 *	(int)num	This is the nonnegative length of `atoms[]'.
 *
 *	(int)length	This is the positive length of one row of the tree.
 *
 *	(int)maxlevel	This is the nonnegative maximum level in the tree.
 *
 *  Outputs:
 *	(int)tfa1sinabt	The return value is TRUE (nonzero) if they fit
 *			or FALSE (zero) if they don't.
 *
 *  External functions called:
 *	tfa1inabt()
 */
extern int
  tfa1sinabt(
	     tfa1 *atoms,	/* Test if these fit into the tree. */
	     int   num,		/* Test this many TFA1s. */
	     int   length,	/* Length of one row of the tree. */
	     int   maxlevel)	/* Maximum level in the tree. */
{
  int k;
  for(k=0; k<num; k++)
    if( !tfa1inabt(atoms+k, length, maxlevel) ) return(FALSE);
  return(TRUE);
}

/****************************************************************
 * abt2tfa1()
 *
 *	Superpose an amplitude from an array binary tree onto the
 *	AMPLITUDE member of a partially assigned TFA1, using its
 *	BLOCK, LEVEL and OFFSET members to locate the amplitude
 *	within the tree.
 *
 *  Calling sequence and basic algorithm:
 *
 *   abt2tfa1( ATOM, DATA, LENGTH ):
 *       Let BLOCK = abtblock( LENGTH, ATOM.LEVEL, ATOM.BLOCK )
 *       ATOM.AMPLITUDE += DATA[BLOCK+ATOM.OFFSET]
 *
 *  Inputs:
 *	(tfa1 *)atom	This TFA1 should have its BLOCK, LEVEL and OFFSET
 *			  members assigned.
 *
 *	(real *)data	This is a preassigned array binary tree.
 *
 *	(int)length	This is the positive length of one row of the tree.
 *
 *  Outputs:
 *	The amplitude at (`atom->level', `atom->block', `atom->offset')
 *	  in the array binary tree `data[]' is superposed onto
 *	  `atom->amplitude' by side effect.
 *
 *  External functions called:
 *	abtblock()
 */
extern void
  abt2tfa1(
	   tfa1 *atom,		/* Partially assigned with tags. */
	   real *data,		/* Source of target amplitudes. */
	   int   length)	/* Length of one row of the tree. */
{
  int block;
  block = abtblock(length, atom->level, atom->block);
  atom->amplitude += data[block+atom->offset];
}

/****************************************************************
 * abt2tfa1s()
 *
 *	Superpose amplitudes from an array binary tree onto the
 *	AMPLITUDE members of an array of partially assigned TFA1,
 *	using the BLOCK, LEVEL and OFFSET members to locate the
 *	amplitude within the tree.
 *
 *  Calling sequence and basic algorithm:
 *
 *   abt2tfa1s( ATOMS, NUM, DATA, LENGTH ):
 *      For K = 0 to NUM-1
 *	   abt2tfa1( ATOMS+K, DATA, LENGTH )
 *
 *  Inputs:
 *	(tfa1 *)atoms	These TFA1s should have their BLOCK, LEVEL
 *			  and OFFSET members assigned.
 *
 *	(int)num	This is the nonnegative length of `atoms[]'.
 *
 *	(real *)data	This is a preassigned array binary tree.
 *
 *	(int)length	This is the positive length of one row of the tree.
 *
 *  Outputs:
 *	Amplitudes at (`atoms[k].level', `atoms[k].block', `atoms[k].offset')
 *	  in the array binary tree `data[]' are superposed onto
 *	  `atoms[k].amplitude' by side effect.
 *
 *  External functions called:
 *	abt2tfa1()
 */
extern void
  abt2tfa1s(
	    tfa1 *atoms,	/* Partially assigned with tags. */
	    int   num,		/* Length of `atoms[]'.*/
	    real *data,		/* Source of target amplitudes. */
	    int   length)	/* Length of one row of the tree. */
{
  int k;

  assert(num>=0);
  assert(data);
  assert(atoms);

  for(k=0; k<num; k++)
    abt2tfa1(atoms+k, data, length);
}

/****************************************************************
 * tfa12btnt()
 *
 *	Superpose the AMPLITUDE from one TFA1 into a partial BTN
 *	tree.  If the target node does not exist, create it.  If
 *	the INTERVAL content of the target node is too small,
 *	enlarge it.
 *
 *  Calling sequence and basic algorithm:
 *
 *  tfa12btnt( ROOT, ATOM ):
 *     Let NODE = btn2branch( ROOT, ATOM.LEVEL, ATOM.BLOCK )
 *     If NODE.CONTENT is NULL then
 *	  Let NODE.CONTENT = makeinterval(0, ATOM.OFFSET, ATOM.OFFSET)
 *     Else
 *	  Let LEAST = min( NODE.CONTENT.LEAST, ATOM.OFFSET )
 *	  Let FINAL = max( NODE.CONTENT.FINAL, ATOM.OFFSET )
 *	  enlargeinterval( NODE.CONTENT, LEAST, FINAL )
 *     NODE.CONTENT.ORIGIN[ATOM.OFFSET] += ATOM.AMPLITUDE
 *
 *  Inputs:
 *	(btn *)root	This must be the pre-allocated root of a tree.
 *
 *	(tfa1 *)atom	This must be a fully assigned TFA1 data structure.
 *
 *  Outputs:
 *	`atom->amplitude' is superposed at index `atom->offset' of the
 *	  INTERVAL content of the BTN at (`atom->level', `atom->block'),
 *	  with respect to the tree at `root', by side effect.  If
 *	  necessary, the BTN tree and the INTERVAL are enlarged.
 *
 *  External functions called:
 *	btn2branch(), enlargeinterval(), makeinterval()
 */
extern void
  tfa12btnt(
	    btn  *root,		/* Root of the target tree. */
	    tfa1 *atom)		/* Fully assigned, with tags. */
{
  btn *node;
  interval *content;
  int least, final;

  assert(root);
  assert(atom);

  node = btn2branch( root, atom->level, atom->block );
  if(!node->content)
    {
      node->content = content = makeinterval(0, atom->offset, atom->offset);
    }
  else
    {
      content = (interval *)(node->content);
      least = min( content->least, atom->offset );
      final = max( content->final, atom->offset );
      enlargeinterval( content, least, final );
    }
  content->origin[atom->offset] += atom->amplitude;
}

/****************************************************************
 * tfa1s2btnt()
 *
 *	Superpose the AMPLITUDEs from a list of TFA1s into a
 *	partial BTN tree.  If the target node does not exist,
 *	create it.  If the INTERVAL content of the target node
 *	is too small, enlarge it.
 *
 *  Calling sequence and basic algorithm:
 *
 *  tfa1s2btnt( ROOT, ATOMS, NUM ):
 *     For K = 0 to NUM -1
 *        tfa12btnt( ROOT, ATOMS[K] )
 *
 *  Inputs:
 *	(btn *)root	This must be the pre-allocated root of a tree.
 *
 *	(tfa1 *)atoms	This must be an array of fully assigned TFA1
 *			  data structures.
 *
 *	(int)num	This is the nonnegative length of `atoms[]'.
 *
 *  Outputs:
 *	`atoms[k].amplitude' is superposed at index `atoms[k].offset'
 *	  of the INTERVAL content of the BTN at (`atoms[k]level',
 *	 `atom->block'), with respect to the tree at `root' by side
 *	 effect.  If necessary, the BTN tree and INTERVALs are enlarged.
 *
 *  External functions called:
 *	tfa12btnt()
 */
extern void
  tfa1s2btnt(
	     btn  *root,	/* Root of the target tree. */
	     tfa1 *atoms,	/* Fully assigned, with tags. */
	     int   num)		/* Elements in `atoms[]'. */
{
  int k;

  assert(root);
  assert(atoms);

  for(k=0; k<num; k++)
    tfa12btnt( root, atoms+k );
}

/****************************************************************
 * btnt2tfa1()
 *
 *	Superpose an amplitude from a BTN tree into a partially
 *	assigned TFA1, using the tags to locate the amplitude.
 *	If the target node or its INTERVAL don't exist, do nothing.
 *
 *  Calling sequence and basic algorithm:
 *
 *   btnt2tfa1( ATOM, ROOT ):
 *     Let NODE = btnt2btn( ROOT, ATOM.LEVEL, ATOM.BLOCK )
 *     If NODE != NULL
 *        If NODE.CONTENT != NULL
 *	  If ininterval( NODE.CONTENT, ATOM.OFFSET ) then
 *	     ATOM.AMPLITUDE += NODE.CONTENT.ORIGIN[ATOM.OFFSET]
 *
 *  Inputs:
 *	(tfa1 *)atom	This must be a partially assigned TFA1 structure.
 *
 *	(btn *)root	This must be the pre-allocated root of a tree.
 *
 *  Outputs:
 *	The amplitude at (`atom->level', `atom->block', `atom->offset')
 *	  in the BTN tree of INTERVALs at `root', if it exists,  is
 *	  superposed onto `atom->amplitude'.
 *
 *  External functions called:
 *	btnt2btn(), ininterval()
 */
extern void
  btnt2tfa1(
	    tfa1 *atom,		/* Partially assigned, with tags. */
	    btn *root)		/* Root of the source tree. */
{
  btn *node;
  interval *content;

  assert(atom);
  
  node = btnt2btn( root, atom->level, atom->block );
  if(node)
    if(node->content)
      {
	content = (interval *)(node->content);
	if( ininterval( content, atom->offset ) )
	  atom->amplitude += content->origin[atom->offset];
      }
}

/****************************************************************
 * btnt2tfa1s()
 *
 *	Superpose amplitudes from a BTN tree into a list of partially
 *	assigned TFA1s, using their tags to locate the amplitudes.
 *	If a target node or its INTERVAL doesn't exist, skip it.
 *
 *  Calling sequence and basic algorithm:
 *
 *   btnt2tfa1s( ATOMS, NUM, ROOT ):
 *     For K = 0 to NUM-1
 *        btnt2tfa1( ATOMS[K], ROOT )
 *
 *  Inputs:
 *	(tfa1 *)atoms	This must be an array of partially assigned
 *			  TFA1 data structures.
 *
 *	(int)num	This the nonnegative length of `atoms[]'.
 *
 *	(btn *)root	This must be the pre-allocated root of a tree.
 *
 *  Outputs:
 *	The amplitudes at (`atoms[k].level', `atoms[k].block',
 *	  `atoms[k].offset') in the BTN tree of INTERVALs at `root',
 *	  if they exist,  are superposed onto `atoms[k].amplitude'.
 *
 *  External functions called:
 *	btnt2tfa1()
 */
extern void
  btnt2tfa1s(
	     tfa1 *atoms,	/* Partially assigned, with tags. */
	     int   num,		/* Elements in `atoms[]'. */
	     btn  *root)	/* Root of the source tree. */
{
  int k;

  assert(atoms);

  for(k=0; k<num; k++)
    btnt2tfa1( atoms+k, root );
}

/***********************************************************************
 * array2tfa1s()
 *
 *	Assign TFA1s with a fixed block and level from an array
 *	of amplitudes.
 *
 *  Calling sequence and basic algorithm:
 *
 *    array2tfa1s( ATOMS, NUM, AMPLITUDES, BLOCK, LEVEL ):
 *       For I = 0 to NUM-1
 *	    Let ATOMS[I].AMPLITUDE = AMPLITUDES[I]
 *	    Let ATOMS[I].BLOCK = BLOCK
 *	    Let ATOMS[I].LEVEL = LEVEL
 *	    Let ATOMS[I].OFFSET = I
 *
 *  Inputs:
 *
 *	(tfa1 *)atoms	This must be preallocated to length `num'.
 *
 *	(int)num	This is the positive number of elements of
 *			  `atoms[]' and `amplitudes[]'.
 *
 *	(real *)amplitudes	Assign these `num' amplitudes.
 *
 *	(int)block	Use this nonnegative block index.
 *
 *	(int)level	Use this as the small nonnegative level index.
 *
 *  Outputs:
 *	Amplitudes are assigned to the first `num' TFA1s in `atoms[]'.
 */
extern void
  array2tfa1s(
	      tfa1 *atoms,	/* `num' preallocated TFA1s to assign. */
	      int   num,	/* Elements in `atoms[]', `amplitudes[]' */
	      real *amplitudes,	/* Source amplitudes to assign */
	      int   block,	/* Block index of the assigned TFA1s */
	      int   level)	/* Level index of the assigned TFA1s */
{
  int i;

  assert( (atoms && amplitudes) || !num );

  for(i=0; i<num; i++)
    {
      atoms->block = block;
      atoms->level = level;
      atoms->offset = i;
      atoms->amplitude = *amplitudes++;
      ++atoms;
    }
}

/***********************************************************************
 * hedgeabt2tfa1s()
 *
 * Convert an array binary tree HEDGE of known depth to a list of TFA1s.
 *
 *  Calling sequence and basic algorithm:
 *
 *  hedgeabt2tfa1s( ATOMS, GRAPH, LENGTH, MAXLEVEL ):
 *    Let START = 0
 *    For J = 0 to GRAPH.BLOCKS-1
 *      Let LEVEL = GRAPH.LEVELS[J]
 *      Let BLOCK = START>>(MAXLEVEL-LEVEL)
 *      Let NUM = abtblength( LENGTH, LEVEL )
 *      START += NUM
 *      array2tfa1s( ATOMS, NUM, GRAPH.CONTENTS+J, BLOCK, LEVEL )
 *      ATOMS += NUM
 *
 *  Inputs:
 *	(tfa1 *)atoms	This must be a preallocated array.
 *
 *	(hedge *)graph	This must be a preallocated and assigned
 *			  data structure.
 *
 *	(int)length	This is the positive number of columns in one
 *			  row of the underlying array binary tree.
 *
 *	(int)maxlevel	Use this small nonnegative integer as the 
 *			  maximum level of the array binary tree.
 *
 *  Output:
 *	Amplitudes and their level, block and offset tags are
 *	  assigned into the array of TFA1s by side effect.
 *
 *  External functions called:
 *	array2tfa1s(), abtblength()
 */
extern void
  hedgeabt2tfa1s(
		 tfa1  *atoms,	/* Preallocated target array. */
		 hedge *graph,	/* Preallocated source HEDGE. */
		 int    length,	/* Number of columns in tree. */
		 int  maxlevel)	/* Depth of the underlying tree. */
{
  int j, start, level, block, num;
  real *amplitudes;

  assert(graph);
  assert( atoms || !graph->blocks);

  start = 0;
  for(j=0; j<graph->blocks; j++)
    {
      level = graph->levels[j];
      block = start>>(maxlevel-level);
      num = abtblength(length,level);
      start += num;
      amplitudes = (real *)(graph->contents[j]);
      array2tfa1s( atoms, num, amplitudes, block, level );
      atoms += num;
    }
}

/***********************************************************************
 * abthedge2tfa1s()
 *
 * Convert an array binary tree HEDGE to a list of TFA1s.
 *
 *  Calling sequence and basic algorithm:
 *
 *  abthedge2tfa1s( ATOMS, GRAPH, LENGTH ):
 *    Let MAXLEVEL = 0
 *    Let FRACTION = 0
 *    For J = 0 to GRAPH.BLOCKS-1
 *       Let LEVEL = GRAPH.LEVELS[J]
 *       If LEVEL>MAXLEVEL then
 *	    FRACTION <<= LEVEL-MAXLEVEL
 *	    Let MAXLEVEL = LEVEL
 *	    Let BLOCK = FRACTION
 *       Else
 *	    Let BLOCK = FRACTION>>(MAXLEVEL-LEVEL)
 *       Let NUM = abtblength( LENGTH, LEVEL )
 *       array2tfa1s( ATOMS, NUM, GRAPH.CONTENTS+J, BLOCK, LEVEL )
 *       FRACTION += 1<<LEVEL
 *       ATOMS += NUM
 *    Return MAXLEVEL
 *
 *  Inputs:
 *	(tfa1 *)atoms	This must be a preallocated array.
 *
 *	(hedge *)graph	This must be a preallocated and assigned
 *			  data structure.
 *
 *	(int)length	This is the positive number of columns in one
 *			  row of the underlying array binary tree.
 *
 *  Output:
 *	Amplitudes and their level, block and offset tags are
 *	  assigned into the array of TFA1s by side effect.
 *
 *  External functions called:
 *	array2tfa1s(), abtblength()
 */
extern int
  abthedge2tfa1s(
		 tfa1  *atoms,	/* Preallocated target array. */
		 hedge *graph,	/* Preallocated source HEDGE. */
		 int    length)	/* Number of columns in tree. */
{
  int  maxlevel, fraction, i, level, block, num;
  real *amplitudes;

  assert(graph);
  assert( atoms || !graph->blocks);

  maxlevel = 0;
  fraction = 0;
  for(i=0; i<graph->blocks; i++)
    {
      level = graph->levels[i];
      if (level>maxlevel)
	{
	  fraction <<= (level-maxlevel);
	  maxlevel = level;
	  block = fraction;
	}
      else
	block = fraction >> (maxlevel-level);
      num = abtblength( length, level );
      amplitudes = (real *)(graph->contents[i]);
      array2tfa1s( atoms, num, amplitudes, block, level );
      fraction += 1<<level;
      atoms += num;
    }
  return(maxlevel);
}

/***********************************************************************
 * interval2tfa1s()
 *
 *	Assign TFA1s with a fixed block and level from an INTERVAL
 *	data structure of amplitudes.
 *
 *  Calling sequence and basic algorithm:
 *
 *    interval2tfa1s( ATOMS, SEGMENT, BLOCK, LEVEL ):
 *       For I = 0 to SEGMENT.FINAL-SEGMENT.LEAST
 *	   Let ATOMS[I].AMPLITUDE = SEGMENT.ORIGIN[I+SEGMENT.LEAST]
 *	   Let ATOMS[I].BLOCK = BLOCK
 *	   Let ATOMS[I].LEVEL = LEVEL
 *	   Let ATOMS[I].OFFSET = I
 *
 *  Inputs:
 *
 *	(tfa1 *)atoms	This must be preallocated to length `num'.
 *
 *	(interval *)segment	We assign from this source.
 *
 *	(int)block	This is the nonnegative index of the block.
 *
 *	(int)level	This is the small nonnegative level index.
 *
 *  Outputs:
 *	Amplitudes are assigned to the TFA1s in `atoms[]'.
 */
extern void
  interval2tfa1s(
		 tfa1 *atoms,	/* The preallocated TFA1s to assign. */
		 interval *segment,	/* Source amplitudes to assign */
		 int   block,	/* Block index of the assigned TFA1s */
		 int   level)	/* Level index of the assigned TFA1s */
{
  int i;

  assert( atoms );
  assert( segment );

  for(i=segment->least; i<segment->final; i++)
    {
      atoms->block = block;
      atoms->level = level;
      atoms->offset = i;
      atoms->amplitude = segment->origin[i];
      ++atoms;
    }
}

/***********************************************************************
 * intervalhedge2tfa1s()
 *
 * Convert a HEDGE containing INTERVALs to a list of TFA1s.
 *
 *  Calling sequence and basic algorithm:
 *
 *    intervalhedge2tfa1s( ATOMS, GRAPH ):
 *       Let MAXLEVEL = 0
 *       Let FRACTION = 0
 *       For J = 0 to GRAPH.BLOCKS-1
 *          Let LEVEL = GRAPH.LEVELS[J]
 *          If LEVEL>MAXLEVEL then
 *             FRACTION <<= LEVEL-MAXLEVEL
 *             Let MAXLEVEL = LEVEL
 *             Let BLOCK = FRACTION
 *          Else
 *             Let BLOCK = FRACTION>>(MAXLEVEL-LEVEL)
 *          Let SEGMENT = GRAPH.CONTENTS + J
 *          interval2tfa1s( ATOMS, SEGMENT, BLOCK, LEVEL )
 *          ATOMS += 1 + SEGMENT.FINAL - SEGMENT.LEAST
 *          FRACTION += 1<<LEVEL
 *       Return MAXLEVEL
 *
 *  Inputs:
 *	(tfa1 *)atoms	This must be a preallocated array.
 *
 *	(hedge *)graph	This must be a preallocated and assigned
 *			  data structure.
 *
 *  Output:
 *	(int)intervalhedge2tfa1s	The return value is the depth
 *					of the underlying tree.
 *
 *	Amplitudes and their level, block and offset tags are
 *	  assigned into the array of TFA1s by side effect.
 *
 *  External functions called:
 *	interval2tfa1s()
 */
extern int
  intervalhedge2tfa1s(
		      tfa1  *atoms,	/* Preallocated target array. */
		      hedge *graph)	/* Preallocated source HEDGE. */
{
  int  maxlevel, fraction, i, level, block;
  interval *segment;

  assert( graph );
  assert( atoms || !graph->blocks);

  maxlevel = 0;
  fraction = 0;
  for(i=0; i<graph->blocks; i++)
    {
      level = graph->levels[i];
      if ( level>maxlevel )
	{
	  fraction <<= (level-maxlevel);
	  maxlevel = level;
	  block = fraction;
	}
      else
	block = fraction >> (maxlevel-level);
      segment = (interval *)(graph->contents[i]);
      interval2tfa1s( atoms, segment, block, level );
      atoms += 1 + segment->final - segment->least;
      fraction += 1<<level;
    }
  return(maxlevel);
}

/***********************************************************************
 * abt2btnt()
 *
 *   Convert an array binary tree to a BTN tree whose contents are
 *   pointers into the array.
 *
 *  Calling sequence and basic algorithm:
 *
 *   abt2btnt( DATA, LENGTH, MAXLEVEL, LEVEL ):
 *     Let ROOT = makebtn( DATA, NULL, NULL, NULL )
 *     If LEVEL<MAXLEVEL then
 *       Let CHILD = DATA+LENGTH
 *       Let ROOT.LEFT = abt2btnt(CHILD, LENGTH, MAXLEVEL, LEVEL+1)
 *       Let CHILD = DATA+LENGTH+(LENGTH>>LEVEL)/2
 *       Let ROOT.RIGHT = abt2btnt(CHILD, LENGTH, MAXLEVEL, LEVEL+1)
 *     Return ROOT
 *
 *  Inputs:
 *	(real *)data	This must be a pointer into a preallocated
 *			  and assigned array binary tree.
 *
 *	(int)length	This is the positive number of columns in
 *			  one row of the array binary tree.
 *
 *	(int)maxlevel	This is the nonnegative depth of the
 *			  array binary tree.
 *
 *	(int)level	This is the nonnegative current level
 *			  in the array binary tree.
 *
 *  Output:
 *	(int)abt2btnt	The return value is a BTN data structure
 *			  which is the root of a tree.
 *
 *	A complete BTN tree is allocated and linked, with the content
 *	  members of its nodes pointing into the array binary tree.
 *
 *  External functions called:
 *	makebtn()
 */
extern btn *
  abt2btnt(
	   real *data,		/* Preallocated array binary tree. */
	   int   length,	/* One row of the array binary tree. */
	   int   maxlevel,	/* Depth of the array binary tree. */
	   int   level)		/* `data' points to start of this level. */
{
  btn *root;
  real *child;

  root = makebtn( data, 0, 0, 0 );
  if(level<maxlevel)
    {
      child = data + length;
      root->left = abt2btnt( child, length, maxlevel, level+1);
      child = data + length + (length>>level)/2;
      root->right = abt2btnt( child, length, maxlevel, level+1);
    }
  return(root);
}


