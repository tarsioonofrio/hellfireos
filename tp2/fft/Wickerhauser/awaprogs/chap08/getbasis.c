/*
 * 
 * These functions extract various basis subsets from library trees.
 * 
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <assert.h>
#include <stdlib.h>
#include "abt.h"
#include "infocost.h"
#include "getbasis.h"
#include "utility.h"

/*********************************************************************
 * costs2bbasis()
 *
 * This function extracts a best-basis from a BTN tree
 * with costs tags.
 *
 * Calling order and basic algorithm:
 *   costs2bbasis( GRAPH, ROOT, LEVEL ):
 *      If ROOT.LEFT==NULL && ROOT.RIGHT==NULL then
 *         Let GRAPH.LEVELS[GRAPH.BLOCKS] = LEVEL
 *         Let GRAPH.CONTENTS[GRAPH.BLOCKS] = ROOT.CONTENT
 *         GRAPH.BLOCKS += 1
 *         Let BESTCOST = ROOT.TAG
 *      Else
 *         Let BLOCKS = GRAPH.BLOCKS
 *         Let COST = 0
 *         If ROOT.LEFT != NULL then 
 *            COST += costs2bbasis( GRAPH, ROOT.LEFT, LEVEL+1 )
 *         If ROOT.RIGHT != NULL then
 *            COST += costs2bbasis( GRAPH, ROOT.RIGHT, LEVEL+1 )
 *         If ROOT.TAG>COST then
 *            Let BESTCOST = COST
 *         Else
 *            Let BESTCOST = ROOT.TAG
 *            Let GRAPH.BLOCKS = BLOCKS
 *            Let GRAPH.LEVELS[GRAPH.BLOCKS] = LEVEL
 *            Let GRAPH.CONTENTS[GRAPH.BLOCKS] = ROOT.CONTENT
 *            GRAPH.BLOCKS += 1
 *   Return BESTCOST
 *
 *
 * Input:
 *	(hedge *)graph		This points to a partially assigned data
 *				  structure with enough contents members
 *				  to hold the longest possible basis.
 *
 *	(btn *)root		This is the root of the current BTN tree,
 *				  all of whose `tag' members must contain
 *				  information costs.
 *
 *	(int)level		This is the current level in the BTN tree.
 *
 * Output:
 *	(real)costs2bbasis	The return value is the information cost of
 *				  of the best graph basis below `root'.  The
 *				  levels description of this basis and the
 *				  contents of its nodes are assigned to 
 *				  `graph' by side effect.
 *
 * Assumptions:
 *	1. `graph' is non-NULL
 *	2. `root' is non-NULL
 *	3. level>=0
 *
 * External functions called:
 *	assert()
 */
extern real
  costs2bbasis(			/* Return the best information cost.  */
	       hedge *graph,	/* Empty, to receive the best basis.  */
	       btn *root,	/* Root of the cost-tagged subtree.   */
	       int level)	/* Level index of `root' in the tree. */
{
  real bestcost, cost;
  int blocks;

  assert(graph);
  assert(root);
  assert(level>=0);
  
  if( !root->left && !root->right )
    {
      graph->levels[graph->blocks] = level;
      graph->contents[graph->blocks] = root->content;
      graph->blocks += 1;
      bestcost = *(real *)root->tag;
    }
  else
    {
      blocks = graph->blocks;
      cost = 0;
      if( root->left ) cost += costs2bbasis( graph, root->left, level+1 );
      if( root->right ) cost += costs2bbasis( graph, root->right, level+1 );
      if( *(real *)root->tag > cost )
	{
          bestcost = cost;
	}
      else
	{
	  bestcost = *(real *)root->tag;
	  graph->blocks = blocks;
	  graph->levels[graph->blocks] = level;
	  graph->contents[graph->blocks] = root->content;
	  graph->blocks += 1;
	}
    }
  return(bestcost);
}

/*********************************************************************
 * btnt2costs()
 *
 * This function puts costs tags into a BTN tree of INTERVALs.
 *
 * Calling order and basic algorithm:
 *   btnt2costs( ROOT ):
 *      If ROOT != NULL
 *         Let I = NODE.CONTENT
 *         Let ROOT.TAG = infocost(I.ORIGIN, I.LEAST, I.FINAL)
 *         btnt2costs( ROOT.LEFT )
 *         btnt2costs( ROOT.RIGHT )
 *      Return
 *
 * Input:
 *	(btn *)root		This is the root of the current BTN tree.
 *
 * Output:
 *	The INTERVAL contents of the BTN data structures in the tree
 *	 below `root' have their information costs evaluated with the
 *	 function `infocost()' and the result placed into newly allocated
 *	 tag members by side effect.
 *
 * External functions called:
 *	malloc(), assert(), infocost()
 */
extern void
  btnt2costs(			/* Return the best information cost.  */
	     btn *root)		/* Root of the current subtree.   */
{
  if( root )
    {
      interval *i;

      i = (interval *)root->content;
      if(!root->tag)
	{
	  root->tag = malloc(sizeof(real));    assert(root->tag);
	  assert(root->tag);
	}
      *(real *)root->tag = infocost(i->origin, i->least, i->final);
      btnt2costs( root->left );
      btnt2costs( root->right );
    }
   return;
}

/*********************************************************************
 * btnt2bbasis()
 *
 * This function returns the best graph basis of INTERVALs.
 *
 * Calling order and basic algorithm:
 *   btnt2bbasis( ROOT, MAXLEVEL ):
 *      btnt2costs( ROOT )
 *      Let GRAPH = makehedge( 1<<MAXLEVEL, NULL, NULL, NULL )
 *      Let GRAPH.BLOCKS = 0
 *      costs2bbasis( GRAPH, ROOT, 0 )
 *      Return GRAPH
 *
 * Input:
 *	(btn *)root		This is the root of the current BTN tree.
 *
 *	(int)maxlevel		Search to this depth for the best basis.
 *
 * Output:
 *	(hedge *)btnt2bbasis	The return value contains the INTERVAL
 *				  contents of the best basis BTN data
 *				  structures in the tree at `root'.
 *
 *  Assumptions:
 *	1. maxlevel>=0
 *	2. `root' is non-NULL
 *
 * External functions called:
 *	btnt2costs(), makehedge(), costs2bbasis(), assert()
 */
extern hedge *
  btnt2bbasis(			/* Return the best graph basis.   */
	      btn *root,	/* Root of the current subtree.   */
	      int maxlevel)	/* Search the tree to this depth. */
{
  hedge *graph;

  assert(maxlevel>=0);
  assert(root);

  btnt2costs( root );

  /* There can be no more than `1<<maxlevel' blocks in the graph basis: */
  graph = makehedge( 1<<maxlevel, 0, 0, 0 );
  graph->blocks = 0;
  costs2bbasis( graph, root, 0 );
  return( graph );
}

/*********************************************************************
 * abt2costs()
 *
 * Build a costs BTN tree from an array binary tree
 *
 * Calling order and basic algorithm:
 *   abt2costs( DATA, LENGTH, MAXLEVEL ):
 *     Let COSTS = makebtnt( MAXLEVEL )
 *     For LEVEL = 0 to MAXLEVEL
 *       For BLOCK = 0 to (1<<LEVEL)-1
 *         Let ORIGIN = DATA + abtblock( LENGTH, LEVEL, BLOCK )
 *         Let BLENGTH = abtblength( LENGTH, LEVEL )
 *         Let CNODE = btnt2btn( COSTS, LEVEL, BLOCK )
 *         Let CNODE.TAG = infocost( ORIGIN, 0, BLENGTH-1 )
 *         Let CNODE.CONTENT = ORIGIN
 *     Return COSTS
 *
 * Input:
 *	(real *)data		This is an array binary tree of amplitudes.
 *
 *	(int)length		Each row of `data[]' is this long.
 *
 *	(int)maxlevel		Search to this depth for the best basis.
 *
 * Output:
 *	(btn *)abt2costs	The return value is a pointer to the root
 *				  of the newly allocated costs BTN tree.
 *
 *  Assumptions:
 *	1. maxlevel>=0
 *	2. `length' is divisible by `1<<maxlevel'
 *	3. `data[]' is non-NULL
 *
 * External functions called:
 *	abtblock(), abtblength(), infocost(), btnt2btn(), makebtnt(), assert()
 */
extern btn *
  abt2costs(			/* Return a costs-tagged BTN tree.   */
	    real *data,		/* The underlying array binary tree. */
	    int length,		/* Length of one row of `data[]'.    */
	    int maxlevel)	/* This is the depth of `data[]'.    */
{
  btn *costs, *cnode;
  int level, block, blength;
  real *origin;

  assert(maxlevel>=0);
  assert(data);
  assert( (length>>maxlevel)<<maxlevel == length );

  costs = makebtnt( maxlevel );
  for( level=0; level<=maxlevel; level++ )
    for( block=0; block<(1<<level); block++ )
      {
	origin = data + abtblock( length, level, block );
	blength = abtblength( length, level );
	cnode = btnt2btn( costs, level, block );
	cnode->tag = malloc(sizeof(real));  assert(cnode->tag);
	*(real *)cnode->tag = infocost( origin, 0, blength-1 );
	cnode->content = (void *)origin;
      }
  return(costs);
}

/*********************************************************************
 * abt2bbasis()
 *
 * This function returns the best graph basis in an array binary tree.
 *
 * Calling order and basic algorithm:
 *   abt2bbasis( DATA, LENGTH, MAXLEVEL ):
 *      Let COSTS = abt2costs( DATA, LENGTH, MAXLEVEL )
 *      Let GRAPH = makehedge( 1<<MAXLEVEL, NULL, NULL, NULL )
 *      Let GRAPH.BLOCKS = 0
 *      costs2bbasis( GRAPH, COSTS, 0 )
 *      freebtnt( COSTS, NULL, free ) 
 *      Return GRAPH
 *
 * Input:
 *	(real *)data		This points to the array binary tree.
 *
 *	(int)length		One row of `data[]' is this long.
 *
 *	(int)maxlevel		This is the depth of `data[]'.
 *
 * Output:
 *	(hedge *)abt2bbasis	The return value contains the array contents
 *				  of the best graph basis in `data[]'.
 *
 *  Assumptions:
 *	1. maxlevel>=0
 *	2. `data[]' is non-NULL
 *	3. `1<<maxlevel' divides `length'
 *
 * External functions called:
 *	abt2costs(), makehedge(), costs2bbasis(), freebtnt(), freevoid()
 */
extern hedge *
  abt2bbasis(			/* Return the best graph basis.    */
	     real *data,	/* Start of the array binary tree. */
	     int length,	/* Length of one row of the tree.  */
	     int maxlevel)	/* Depth of the array binary tree. */
{
  btn *costs;
  hedge *graph;

  assert(maxlevel>=0);
  assert(data);
  assert( (length>>maxlevel)<<maxlevel == length );

  costs = abt2costs( data, length, maxlevel );
  graph = makehedge( 1<<maxlevel, 0, 0, 0 );
  graph->blocks = 0;
  costs2bbasis( graph, costs, 0 );
  freebtnt( costs, 0, freevoid ); 
  return(graph);
}

/*********************************************************************
 * levelcost()
 *
 * This function returns the cost of one level from a costs BTN tree.
 *
 * Calling order and basic algorithm:
 *   levelcost( ROOT, LEVEL ):
 *     Let COST = 0
 *     For BLOCK = 0 to (1<<LEVEL)-1
 *       Let NODE = btnt2btn( ROOT, LEVEL, BLOCK )
 *       COST += NODE.TAG
 *     Return COST
 *
 * Input:
 *	(btn *)root		This points to the costs BTN tree.
 *
 *	(int)level		This is the level whose cost is desired.
 *
 * Output:
 *	(real)levelcost		The return value is the cost of `level'
 *
 *  Assumptions:
 *	1. The BTN tree at `root' contains a complete `level'
 *	2. `root' is non-NULL
 *
 * External functions called:
 *	btnt2btn(), assert()
 */
extern real
  levelcost(			/* Return the cost of a level.   */
	    btn *root,		/* Root BTN tree.                */
	    int level)		/* Level whose cost is desired.  */
{
  real cost;
  int block;
  btn *node;

  assert(root);
  assert( level>=0 );

  cost = 0;
  for( block=0; block<(1<<level); block++ )
    {
      node = btnt2btn( root, level, block );
      assert(node);
      cost += *(real *)node->tag;
    }
  return(cost);
}

/*********************************************************************
 * costs2blevel()
 *
 * This function extracts a best level basis from a BTN tree
 * with costs tags.
 *
 * Calling order and basic algorithm:
 *   costs2blevel( GRAPH, ROOT, MINLEVEL, MAXLEVEL ):
 *     Let BESTLEVEL = MINLEVEL
 *     Let BESTCOST = levelcost( ROOT, MINLEVEL )
 *     For LEVEL = MINLEVEL+1 to MAXLEVEL
 *       Let COST = levelcost( ROOT, LEVEL )
 *       If COST<BESTCOST then
 *         Let BESTCOST = COST
 *         Let BESTLEVEL = LEVEL
 *     Let GRAPH.BLOCKS = 1<<BESTLEVEL
 *     For BLOCK = 0 to GRAPH.BLOCKS-1
 *       Let GRAPH.LEVELS[BLOCK] = BESTLEVEL
 *       Let GRAPH.CONTENTS[BLOCK] = btnt2btn(ROOT,BESTLEVEL,BLOCK)
 *     Return BESTCOST
 *
 *
 * Input:
 *	(hedge *)graph		This points to a partially assigned data
 *				  structure with enough contents members
 *				  to hold the longest possible basis.
 *
 *	(btn *)root		This is the root of the current BTN tree,
 *				  all of whose `tag' members must contain
 *				  information costs.
 *
 *	(int)minlevel		These are the shallowest and deepest levels
 *	(int)maxlevel		  of the costs BTN tree to search.
 *
 * Output:
 *	(real)costs2blevel	The return value is the information cost of
 *				  of the best level basis below `root'.  The
 *				  levels description of this basis and the
 *				  contents of its nodes are assigned to 
 *				  `graph' by side effect.
 *
 * Assumptions:
 *	1. minlevel<=maxlevel
 *	2. `graph' is non-NULL
 *	3. `root' is non-NULL
 *
 * External functions called:
 *	btnt2btn(), levelcost(), assert()
 */
extern real
  costs2blevel(			/* Return the best level basis cost.     */
	       hedge *graph,	/* Empty, to receive the best level.     */
	       btn *root,	/* Root of the cost-tagged subtree.      */
	       int minlevel,	/* Shallowest level of `root' to search. */
	       int maxlevel)	/* Deepest level of `root' to search.    */
{
  int bestlevel, level, block;
  real bestcost, cost;

  assert( minlevel<=maxlevel );
  assert( graph );
  assert( root );

  bestlevel = minlevel;
  bestcost = levelcost( root, minlevel );
  for(level=minlevel+1; level<=maxlevel; level++ )
    {
      cost = levelcost( root, level );
      if( cost<bestcost )
	{
	  bestcost = cost;
	  bestlevel = level;
	}
    }
  graph->blocks = 1<<bestlevel;
  for( block=0; block<graph->blocks; block++ )
    {
      graph->levels[block] = bestlevel;
      graph->contents[block] = btnt2btn( root, bestlevel, block );
    }
  return(bestcost);
}

/*********************************************************************
 * btnt2bbhedge()
 *
 * This function returns the best graph basis of INTERVALs.
 *
 * Calling order and basic algorithm:
 *   btnt2bbhedge( GRAPH, ROOT, S, L ):
 *     Let MYCOST = infocost( ROOT.CONTENT.ORIGIN,
 *                      ROOT.CONTENT.LEAST, ROOT.CONTENT.FINAL )
 *     If S == L then
 *       Let GRAPH.CONTENT[GRAPH.BLOCKS] = ROOT.CONTENT
 *       Let GRAPH.LEVELS[GRAPH.BLOCKS] = L
 *       GRAPH.BLOCKS += 1
 *       Let BESTCOST = MYCOST
 *     Else
 *       Let BLOCKS = GRAPH.BLOCKS
 *       Let LCOST = btnt2bbhedge( GRAPH, ROOT.LEFT, S+1, L )
 *       Let RCOST = btnt2bbhedge( GRAPH, ROOT.RIGHT, S+1, L )
 *       If MYCOST > LCOST + RCOST then
 *          Let BESTCOST = LCOST + RCOST
 *       Else
 *          Let BESTCOST = MYCOST
 *          Let GRAPH.BLOCKS = BLOCKS
 *          Let GRAPH.CONTENTS[GRAPH.BLOCKS] = ROOT.CONTENT
 *          Let GRAPH.LEVELS[GRAPH.BLOCKS] = S
 *          GRAPH.BLOCKS += 1
 *     Return BESTCOST
 *
 * Input:
 *	(hedge *)graph		This is a partially allocated output HEDGE.
 *
 *	(btn *)root		This is the root of the current BTN tree.
 *
 *	(int)s			This is the current level in `root'.
 *
 *	(int)l			Search to this depth for the best basis.
 *
 * Output:
 *	(real)btnt2bbhedge	The return value is the information cost
 *				  of the best graph basis in `root'.
 *
 *  Assumptions:
 *	1. s<=l
 *	2. `root' is non-NULL
 *	3. `graph' is non-NULL
 *
 * External functions called:
 *	assert(), infocost()
 */
extern real
  btnt2bbhedge(			/* Return the best graph basis cost.   */
	       hedge *graph,	/* Partially allocated output struct.  */
	       btn *root,	/* Root of the current subtree.        */
	       int s,		/* Current level in `root'.            */
	       int l)		/* Search the tree to this depth. */
{
  interval *seg;
  real cost, bestcost;
  int blocks;

  assert(s<=l);
  assert(root);
  assert(graph);

  seg = (interval *)root->content;
  bestcost = infocost( seg->origin, seg->least, seg->final );

  if( s==l )
    {
      graph->contents[graph->blocks] = root->content;
      graph->levels[graph->blocks] = l;
      graph->blocks += 1;
    }
  else
    {
      blocks = graph->blocks;
      cost = 0;
      cost += btnt2bbhedge( graph, root->left, s+1, l );
      cost += btnt2bbhedge( graph, root->right, s+1, l );
      if( bestcost > cost )
	 bestcost = cost;
      else
	{
	  graph->blocks = blocks;
	  graph->contents[graph->blocks] = root->content;
	  graph->levels[graph->blocks] = s;
	  graph->blocks += 1;
	}
    }
  return(bestcost);
}

/*********************************************************************
 * btnt2blevel()
 *
 * This function returns the best level basis of INTERVALs.
 *
 * Calling order and basic algorithm:
 *   btnt2blevel( ROOT, MINLEVEL, MAXLEVEL ):
 *     btnt2costs( ROOT )
 *     Let GRAPH = makehedge( 1<<MAXLEVEL, NULL, NULL, NULL )
 *     Let GRAPH.BLOCKS = 0
 *     costs2blevel( GRAPH, ROOT, MINLEVEL, MAXLEVEL )
 *     Return GRAPH
 *
 * Input:
 *	(btn *)root		This is the root of the current BTN tree.
 *
 *	(int)minlevel		Search from this depth for the best level.
 *
 *	(int)maxlevel		Search to this depth for the best level.
 *
 * Output:
 *	(hedge *)btnt2blevel	The return value contains the INTERVAL
 *				  contents of the best basis BTN data
 *				  structures in the tree at `root'.
 *
 *  Assumptions:
 *	1. maxlevel>=minlevel
 *	2. `root' is non-NULL
 *
 * External functions called:
 *	btnt2costs(), makehedge(), costs2blevel(), assert()
 */
extern hedge *
  btnt2blevel(			/* Return the best graph basis.     */
	      btn *root,	/* Root of the current subtree.     */
	      int minlevel,	/* Search the tree from this depth. */
	      int maxlevel)	/* Search the tree to this depth.   */
{
  hedge *graph;

  assert( maxlevel>=minlevel );
  assert( root );

  btnt2costs( root );
  graph = makehedge( 1<<maxlevel, 0, 0, 0 );
  graph->blocks = 0;
  costs2blevel( graph, root, minlevel, maxlevel );
  return(graph);
}

/*********************************************************************
 * btntlevelcost()
 *
 * This function returns the cost of one level from an untagged BTN tree.
 *
 * Calling order and basic algorithm:
 *   btntlevelcost( ROOT, LEVEL ):
 *     Let COST = 0
 *     For BLOCK = 0 to (1<<LEVEL)-1
 *       Let NODE = btnt2btn( ROOT, LEVEL, BLOCK )
 *       Let I = NODE.CONTENT
 *       COST += infocost( I.ORIGIN, I.LEAST, I.FINAL )
 *     Return COST
 *
 * Input:
 *	(btn *)root		This points to the BTN tree.
 *
 *	(int)level		This is the level whose cost is desired.
 *
 * Output:
 *	(real)btntlevelcost	The return value is the cost of `level'
 *
 *  Assumptions:
 *	1. The BTN tree at `root' contains a complete `level'
 *	2. `root' is non-NULL
 *	3. level>=0
 *
 * External functions called:
 *	btnt2btn(), assert()
 */
extern real
  btntlevelcost(		/* Return the cost of a level.  */
		btn *root,	/* Root of the BTN tree.        */
		int level)	/* Level whose cost is desired. */
{
  real cost;
  interval *seg;
  btn *node;
  int block;

  assert(root);
  assert(level>=0);

  cost = 0;
  for( block=0; block<(1<<level); block++ )
    {
      node = btnt2btn( root, level, block );
      seg = (interval *)node->content;
      cost += infocost( seg->origin, seg->least, seg->final );
    }
  return(cost);
}

/*********************************************************************
 * btnt2blhedge()
 *
 * This function returns the best level basis cost.
 *
 * Calling order and basic algorithm:
 *   btnt2blhedge( GRAPH, ROOT, MINLEVEL, MAXLEVEL ):
 *     Let BESTLEVEL = MINLEVEL
 *     Let BESTCOST = btntlevelcost( ROOT, MINLEVEL )
 *     For LEVEL = MINLEVEL+1 to MAXLEVEL
 *       Let COST = btntlevelcost( ROOT, LEVEL )
 *       If COST<BESTCOST then
 *         Let BESTCOST = COST
 *         Let BESTLEVEL = LEVEL
 *     Let GRAPH.BLOCKS = 1<<BESTLEVEL
 *     For BLOCK = 0 to GRAPH.BLOCKS-1
 *         Let NODE = btnt2btn( ROOT, LEVEL, BLOCK )
 *         Let GRAPH.CONTENTS[BLOCK] = NODE.CONTENT
 *         Let GRAPH.LEVELS[BLOCK] = BESTLEVEL
 *     Return BESTCOST
 *
 * Input:
 *	(hedge *)graph		Output goes into this partially filled struct.
 *
 *	(btn *)root		This is the root of the current BTN tree.
 *
 *	(int)minlevel		Search from this depth for the best level.
 *
 *	(int)maxlevel		Search to this depth for the best level.
 *
 * Output:
 *	(real)btnt2blhedge	The return value is the information cost
 *				  of the best level, described by the levels
 *				  list in `graph'.
 *
 *  Assumptions:
 *	1. maxlevel>=minlevel
 *	2. `root' is non-NULL
 *	3. `graph' is non-NULL
 *
 * External functions called:
 *	btnt2btn(), btntlevelcost(), assert()
 */
extern real
  btnt2blhedge(			/* Return the best graph basis.      */
	       hedge *graph,	/* Partially allocated output hedge. */
	       btn *root,	/* Root of the current subtree.      */
	       int minlevel,	/* Search the tree from this depth.  */
	       int maxlevel)	/* Search the tree to this depth.    */
{
  btn *node;
  real cost, bestcost;
  int level, block, bestlevel;

  assert( maxlevel>=minlevel );
  assert( root );
  assert( graph );

  bestlevel = minlevel;
  bestcost = btntlevelcost( root, minlevel );

  /* Bubble up the best level: */
  for( level=minlevel+1; level<=maxlevel; level++ )
    {
      cost = btntlevelcost( root, level );
      if( cost<bestcost )
	{
	  bestcost = cost;
	  bestlevel = level;
	}
    }

  /* Assign the output HEDGE data structure: */
  graph->blocks = 1<<bestlevel;
  for( block=0; block<graph->blocks; block++ )
    {
      node = btnt2btn( root, level, block );
      graph->contents[block] = node->content;
      graph->levels[block] = bestlevel;
    }
  return(bestcost);
}

/*********************************************************************
 * abt2blevel()
 *
 * This function returns the best level basis in an array binary tree.
 *
 * Calling order and basic algorithm:
 *   abt2blevel( DATA, LENGTH, MINLEVEL, MAXLEVEL ):
 *     Let COSTS = abt2costs( DATA, LENGTH, MAXLEVEL )
 *     Let GRAPH = makehedge( 1<<MAXLEVEL, NULL, NULL, NULL )
 *     Let GRAPH.BLOCKS = 0
 *     costs2blevel( GRAPH, COSTS, MINLEVEL, MAXLEVEL )
 *     freebtnt( COSTS, NULL, free ) 
 *     Return GRAPH
 *
 * Input:
 *	(real *)data		This points to the array binary tree.
 *
 *	(int)length		One row of `data[]' is this long.
 *
 *	(int)minlevel		Search from this depth in `data[]'.
 *
 *	(int)maxlevel		Stop at this depth in `data[]'.
 *
 * Output:
 *	(hedge *)abt2blevel	The return value contains the array contents
 *				  of the best level basis in `data[]'.
 *
 *  Assumptions:
 *	1. maxlevel>=minlevel
 *	2. `data[]' is non-NULL
 *	3. `1<<maxlevel' divides `length'
 *
 * External functions called:
 *	abt2costs(), makehedge(), costs2bbasis(),
 *	freebtnt(), freevoid(), assert()
 */
extern hedge *
  abt2blevel(			/* Return the best level basis.      */
	     real *data,	/* Start of the array binary tree.   */
	     int length,	/* Length of one row of the tree.    */
	     int minlevel,	/* Minimum search depth in the tree. */
	     int maxlevel)	/* Maximum search depth in the tree. */
{
  btn *costs;
  hedge *graph;

  assert( maxlevel>=minlevel );
  assert( data );
  assert( (length>>maxlevel)<<maxlevel == length );

  costs = abt2costs( data, length, maxlevel );
  graph = makehedge( 1<<maxlevel, 0, 0, 0 );
  graph->blocks = 0;
  costs2blevel( graph, costs, minlevel, maxlevel );
  freebtnt( costs, 0, freevoid );
  return(graph);
}

/*********************************************************************
 * abt2blhedge()
 *
 * This function returns the best level basis in an array binary tree.
 * Since the amplitudes in an array binary tree are stored contiguously by
 * level, we can search for the level with least information cost without
 * introducing any additional data structures.
 *
 * Calling order and basic algorithm:
 *   abt2blhedge( GRAPH, DATA, LENGTH, MINLEVEL, MAXLEVEL ):
 *     Let BESTLEVEL = MINLEVEL
 *     Let BESTCOST = infocost(DATA+LENGTH*MINLEVEL, 0, LENGTH-1)
 *     For LEVEL = MINLEVEL+1 to MAXLEVEL
 *       Let COST =  infocost(DATA+LENGTH*LEVEL, 0, LENGTH-1)
 *       If COST<BESTCOST then
 *         Let BESTCOST = COST
 *         Let BESTLEVEL = LEVEL
 *     Let GRAPH.BLOCKS = 1<<BESTLEVEL
 *     For BLOCK = 0 to GRAPH.BLOCKS-1
 *         Let GRAPH.CONTENTS[BLOCK] = 
 *                DATA + abtblock( LENGTH, BESTLEVEL, BLOCK )
 *         Let GRAPH.LEVELS[BLOCK] = BESTLEVEL
 *     Return BESTCOST
 *
 * Input:
 *	(hedge *)graph		Output goes to this data structure.
 *
 *	(real *)data		This points to the array binary tree.
 *
 *	(int)length		One row of `data[]' is this long.
 *
 *	(int)minlevel		Search from this depth in `data[]'.
 *
 *	(int)maxlevel		Stop at this depth in `data[]'.
 *
 * Output:
 *	(real)abt2blhedge	The return value is the cost of the
 *				  best level basis in `data[]'.
 *
 *  Assumptions:
 *	1. maxlevel>=minlevel
 *	2. `data[]' is non-NULL
 *	3. `1<<maxlevel' divides `length'
 *	4. `graph' is non-NULL
 *
 * External functions called:
 *	infocost(), abtblock(), abtbllength(),
 *	assert()
 */
extern real
  abt2blhedge(			/* Return the best level basis cost. */
	      hedge *graph,	/* Partial HEDGE for the output.     */
	      real *data,	/* Start of the array binary tree.   */
	      int length,	/* Length of one row of the tree.    */
	      int minlevel,	/* Minimum search depth in the tree. */
	      int maxlevel)	/* Maximum search depth in the tree. */
{
  real bestcost, cost;
  int level, block, bestlevel;

  assert( maxlevel>=minlevel );
  assert( data );
  assert( graph );
  assert( (length>>maxlevel)<<maxlevel == length );

  /* Find the best level between `minlevel' and `maxlevel': */
  bestlevel = minlevel;
  bestcost = infocost(data+length*minlevel, 0, length-1);
  for( level=minlevel+1; level<=maxlevel; level++ )
    {
      cost =  infocost(data+length*level, 0, length-1);
      if( cost<bestcost )
	{
	  bestcost = cost;
	  bestlevel = level;
	}
    }

  /* Write to the output HEDGE: */
  graph->blocks = 1<<bestlevel;
  for( block=0; block<graph->blocks; block++ )
    {
      graph->contents[block] = 
	data + abtblock( length, bestlevel, block );
      graph->levels[block] = bestlevel;
    }
  return(bestcost);
}
