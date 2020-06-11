/*
 *
 * Two-dimensional best basis and custom basis development.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <assert.h>
#include "real.h"
#include "qf.h"
#include "hedge.h"
#include "infocost.h"
#include "cd.h"
#include "cdp2d.h"
#include "wpt2.h"

/*******************************************************************
 * bbwp2()
 *
 *  Compute and write the best basis of 2-dimensional isotropic
 *  wavelet packets to an output hedge, using the space-saving
 *  periodic discrete wavelet packet analysis algorithm.
 *
 *  Calling sequence and basic algorithm:
 *
 *   bbwp2( GRAPH, KD, IN, X,Y, S, L, WK, H, G ):
 *     Let XY = X*Y
 *     Let COST = infocost( IN, 0, XY-1 )
 *     If  S<L then
 *       Let BLOCK = GRAPH.BLOCKS
 *       For J = 0 to 3
 *         Let K[J] = KD + J*XY/4
 *       Let KD = K[3]  + XY/4
 *       scdpe2( K[0], K[1], K[2], K[3], IN, X,Y, WK, H, G)
 *       Let KCOST = 0
 *       For J = 0 to 3
 *         KCOST += bbwp2(GRAPH,KD, K[J], X/2,Y/2, S+1, L, WK,H,G)
 *       If KCOST < COST then
 *         Let COST = KCOST
 *       Else
 *         Let GRAPH.LEVELS[BLOCK] = S
 *         Let OUT = GRAPH.CONTENTS[BLOCK]
 *         For I = 0 to XY-1
 *           Let OUT[I] = IN[I]
 *         Let GRAPH.BLOCKS = BLOCK+1
 *         Let GRAPH.CONTENTS[GRAPH.BLOCKS] = OUT + XY
 *     Else
 *       Let GRAPH.LEVELS[GRAPH.BLOCKS] = L
 *       Let OUT = GRAPH.CONTENTS[GRAPH.BLOCKS]
 *       For I = 0 to Y*X-1
 *         Let OUT[I] = IN[I]
 *       GRAPH.BLOCKS += 1
 *       Let GRAPH.CONTENTS[GRAPH.BLOCKS] = OUT + XY
 *     Return COST
 *
 *
 *  Inputs:
 *
 *	(hedge *)graph		This preallocated HEDGE must contain enough
 *				  space for any possible graph basis, i.e.,
 *				  it should have `1<<(2*L)' elements in its
 *				  `contents' and `levels' arrays.
 *
 *	(real *)kd		Entire branches of child subspaces will be
 *				  developed in this array, which must have
 *				  length `(4*X*Y)/3'
 *
 *	(real *)in		This points to the input array, which is left
 *				  unchanged.
 *
 *	(int)x			These positive integers are the numbers of
 *	(int)y			  rows and columns of `in[]', respectively.
 *
 *	(int)s			This small integer is the current level.
 *
 *	(int)L			This small integer is the deepest level.
 *
 *	(real *)wk		This is a scratch array for transposition.
 *
 *	(pqf *)H		Use these low-pass and high-pass periodized
 *	(pqf *)G		  filters for convolution-decimation.
 *
 *  Output:
 *	(real *)bbwp2		The return value is the cost of the best 2-D
 *				  periodic graph basis transform of `in[]'.
 *
 *	The array at `graph->contents[0]' is replaced with its coefficients
 *	in the 2-dimensional wavelet packet best basis defined by the cost
 *	function `infocost()', and the other pointers in `graph->contents[]'
 *	are set to their corresponding wavelet packet subspaces.
 *
 *  Assumptions:
 *	1. x>0; y>0
 *	2. L >= s >= 0
 *	3. in != NULL
 *	4. wk != NULL
 *	5. kd != NULL
 *	6. graph != NULL
 *
 *  External functions called:
 *	scdpe2(), infocost(), assert()
 *
 */
extern real
  bbwp2(
	hedge *graph,		/* Preallocated struct for the best basis  */
	real  *kd,		/* Array to hold one branch of 2-D DWPAP   */
	const real *in,		/* Input array, will not be changed        */
	int    x,		/* Length of one column of `in[]'          */
	int    y,		/* Length of one row of `in[]'             */
	int    s,		/* Current level in the DWPA tree          */
	int    L,		/* Maximum depth of the decomposition      */
	real  *wk,		/* Scratch array for transposition         */
	pqf   *H,		/* Low-pass periodic quadrature filter     */
	pqf   *G)		/* High-pass periodic quadrature filter    */
{
  int block, i, j, xy;
  real *out, *k[4];
  real  cost, kcost;

  assert(x>0);
  assert(y>0);
  assert(s<=L);
  assert(s>=0);
  assert(in);
  assert(wk);
  assert(kd);
  assert(graph);

  xy = x*y;
  cost = infocost( in, 0, xy-1 );
  if(  s<L )
    {
      block = graph->blocks;
      for( j=0; j<4; j++ )
	 k[j] = kd + j*xy/4;
      kd = k[3]  + xy/4;
      scdpe2( k[0], k[1], k[2], k[3], in, x,y, wk, H, G);
      kcost = 0;
      for( j=0; j<4; j++ )
	kcost += bbwp2( graph, kd, k[j], x/2, y/2, s+1, L, wk, H, G);
      if( kcost<cost )
	cost = kcost;
      else
	{
	  graph->levels[block] = s;
	  out = (real *)graph->contents[block];
	  for( i=0; i<xy; i++ )
	    out[i] = in[i];
	  graph->blocks = block+1; 
	  graph->contents[graph->blocks] = out + xy;
	}
    }
  else
    {
      graph->levels[graph->blocks] = L;
      out = graph->contents[graph->blocks];
      for( i=0; i<xy; i++ )
	out[i] = in[i];
      graph->blocks += 1;
      graph->contents[graph->blocks] = out + xy;
    }
  return(cost);
}

/*******************************************************************
 * cbwp2()
 *
 *  Compute and write a custom basis of 2-dimensional isotropic
 *  wavelet packets to an output hedge, doing all development
 *  into periodic discrete wavelet packets in place.
 *
 *  Calling sequence and basic algorithm:
 *
 *   cbwp2( GRAPH, LEVEL, IX, IY, WORK, H, G ):
 *     If LEVEL < GRAPH.LEVELS[GRAPH.BLOCKS]  then
 *       scdpi2( GRAPH.CONTENTS[GRAPH.BLOCKS], IX,IY, WORK, H, G )
 *       For K = 0 to 3
 *          cbwp2( GRAPH, LEVEL+1, IX/2, IY/2, WORK, H, G )
 *     Else
 *       Let GRAPH.CONTENTS[GRAPH.BLOCKS+1] =
 *                      GRAPH.CONTENTS[GRAPH.BLOCKS] + IX*IY
 *       GRAPH.BLOCKS += 1
 *     Return
 *
 *
 *  Inputs:
 *
 *	(hedge *)graph		This preallocated HEDGE must contain enough
 *				  space for the desired graph basis, up to
 *				  `1<<(2*L)' elements in its
 *				  `contents' and `levels' arrays.
 *
 *	(int)level		This small integer is the current level.
 *
 *	(int)ix			These positive integers are the numbers of
 *	(int)iy			  rows and columns of the input.
 *
 *	(real *)work		This is a scratch array for transposition.
 *
 *	(pqf *)H		Use these low-pass and high-pass periodized
 *	(pqf *)G		  filters for convolution-decimation.
 *
 *  Output:
 *	The array at `graph->contents[0]' is replaced with its coefficients
 *	in the 2-dimensional wavelet packet basis defined by the encounter
 *	list `graph->levels[]', and the other pointers in `graph->contents[]'
 *	are set to their corresponding wavelet packet subspaces.
 *
 *  Assumptions:
 *	1. ix>0; iy>0
 *	2. level >= 0
 *	3. graph != NULL
 *	4. work != NULL
 *
 *  External functions called:
 *	assert(), scdpi2()
 *
 */
extern void
  cbwp2(
	hedge *graph,		/* Preallocated struct for the best basis  */
	int    level,		/* Current level in the DWPA tree          */
	int    ix,		/* Length of one column of `in[]'          */
	int    iy,		/* Length of one row of `in[]'             */
	real  *work,		/* Scratch array for transposition         */
	pqf   *H,		/* Low-pass periodic quadrature filter     */
	pqf   *G)		/* High-pass periodic quadrature filter    */
{
  real *out;
  int k;

  assert(ix>0);
  assert(iy>0);
  assert(level>=0);
  assert(work);
  assert(graph);

  if( level<graph->levels[graph->blocks] )
    {
      out = (real *)graph->contents[graph->blocks];
      scdpi2( out, ix,iy, work, H, G );
      for( k=0; k<4; k++ )
	cbwp2( graph, level+1, ix/2, iy/2, work, H, G );
    }
  else
    {
      graph->contents[graph->blocks+1]
	= (real *)graph->contents[graph->blocks] + ix*iy;
      graph->blocks += 1;
    }
  return;
}
