/*
 * Discrete wavelet packet analysis routines in one dimension.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 * 
 */

#include <assert.h>
#include <stdlib.h>
#include "real.h"
#include "common.h"		/* for max() */
#include "utility.h"
#include "abt.h"
#include "btn.h"
#include "cd.h"
#include "fntype.h"
#include "hedge.h"
#include "qf.h"
#include "tfa.h"
#include "dwpa.h"

/*********************************************************
 * dwpap2abt0()
 *
 * [D]iscrete [W]avelet [P]acket [A]nalysis, [P]eriodic
 * to a preallocated [A]rray [B]inary [T]ree
 *
 * Calling sequence and basic algorithm:
 *
 *   dwpap2abt0( DATA, N, MAXLEVEL, HQF, GQF ):
 *      For L = 0 to MAXLEVEL-1
 *         Let NPARENT = abtblength( N, L )
 *         For B = 0 to (1<<L)-1
 *            Let PARENT = DATA + abtblock( N, L, B )
 *            Let CHILD = DATA + abtblock( N, L+1, 2*B )
 *            cdpe( CHILD, 1, PARENT, NPARENT, HQF )
 *            Let CHILD = DATA + abtblock( N, L+1, 2*B+1 )
 *            cdpe( CHILD, 1, PARENT, NPARENT, GQF )
 *
 * Inputs:
 *	(real *)data		This preallocated array binary tree has
 *				  `n' elements per row and `maxlevel' rows.
 *
 *	(int)n			This is the positive row length of `data[]'.
 *
 *	(int)maxlevel		This is the depth of `data[]'.
 *
 *	(const pqf *)HQF	These are QF structs used for low-pass and
 *	(const pqf *)GQF	  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)data		This array is filled by reference with
 *				  wavelet packet coefficients.
 *
 * Assumptions:
 *	1. `n' is divisible by `1<<maxlevel'.
 *	2. `data[]' is preallocated with `n*(1+maxlevel)' elements.
 *
 * External functions called:
 *	cdpe(), abtblock(), abtblength()
 */
extern void
  dwpap2abt0(
	     real *data,	/* Input/output array binary tree. */
	     int n,		/* The row length of `data[]'.     */
	     int maxlevel,	/* `data[]' has `1+maxlevel' rows. */
	     const pqf *HQF,	/* Low-pass QF data structure.     */
	     const pqf *GQF)	/* High-pass QF data structure.    */
{
  real *parent, *child;
  int l, b, nparent;

  assert( (n>>maxlevel)<<maxlevel == n );

  for( l=0; l<maxlevel; l++ )
    {
      nparent = abtblength( n, l );
      for( b=0; b<(1<<l); b++ )
	{
         parent = data + abtblock( n, l, b );
	 child  = data + abtblock( n, l+1, 2*b );
         cdpe( child, 1, parent, nparent, HQF );
	 child  = data + abtblock( n, l+1, 2*b+1 );
         cdpe( child, 1, parent, nparent, GQF );
       }
    }
  return;
}


/*********************************************************
 * dwpap2abt()
 *
 * [D]iscrete [W]avelet [P]acket [A]nalysis, [P]eriodic
 * to an [A]rray [B]inary [T]ree, with preparation.
 *
 * Calling sequence and basic algorithm:
 *
 *   dwpap2abt( IN, N, MAXLEVEL, HQF, GQF ):
 *      Allocate an array of (MAXLEVEL+1)*N REALs at DATA
 *      For I = 0 to N-1
 *         Let DATA[I] = IN[I]
 *      dwpap2abt0( DATA, N, MAXLEVEL, HQF, GQF )
 *      Return DATA
 *
 *
 * Inputs:
 *	(real *)in		This array contains only the `n' inputs.
 *
 *	(int)n			This is the positive row length of `data[]'.
 *
 *	(int)maxlevel		This is the depth of `data[]'.
 *
 *	(const pqf *)HQF	These are QF structs used for low-pass and
 *	(const pqf *)GQF	  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)dwpap2abt	The function returns a pointer to an array
 *				  binary tree filled by reference with
 *				  wavelet packet coefficients.
 *
 * Assumptions:
 *	1. `n' is divisible by `1<<maxlevel'.
 *
 * External functions called:
 *	dwpap2abt0(), calloc()
 */
extern real *
  dwpap2abt(
	     const real *in,	/* Input array, not changed.       */
	     int n,		/* The length of the array `in[]'. */
	     int maxlevel,	/* Expand to depth `maxlevel'.     */
	     const pqf *HQF,	/* Low-pass QF data structure.     */
	     const pqf *GQF)	/* High-pass QF data structure.    */
{
  real *data;
  int i;

  assert( (n>>maxlevel)<<maxlevel == n );

  data = (real *)calloc( (maxlevel+1)*n, sizeof(real));  assert(data);
  for( i=0; i<n; i++ ) data[i] = in[i];
  dwpap2abt0( data, n, maxlevel, HQF, GQF );
  return(data);
}

/*********************************************************
 * dwpap2hedger()
 *
 * In place [D]iscrete [W]avelet [P]acket [A]nalysis,
 * [P]eriodic, to a [HEDGE], [R]ecursive core.
 *
 * Calling sequence and basic algorithm:
 *
 *   dwpap2hedger( GRAPH, J, N, S, HQF, GQF, WORK ):
 *     If GRAPH.LEVELS[J]==S then
 *       J += 1
 *     Else
 *       Let PARENT = GRAPH.CONTENTS[J]
 *       cdpe(  WORK, 1, PARENT, N, HQF )
 *       cdpe( WORK+N/2, 1, PARENT, N, GQF )
 *       For I = 0 to N-1
 *         Let PARENT[I] = WORK[I]
 *       Let J = dwpap2hedger( GRAPH, J, N/2, S+1, HQF, GQF, WORK )
 *       Let GRAPH.CONTENTS[J] = PARENT + N/2
 *       Let J = dwpap2hedger( GRAPH, J, N/2, S+1, HQF, GQF, WORK )
 *     Return J
 *
 * Inputs:
 *	(hedge *)graph		This data structure should be preallocated
 *				  and assigned a valid levels list.
 *
 *	(int)j			This is the current block index.
 *
 *	(int)n			This is the positive signal length.
 *
 *	(int)s			This is the current scale index.
 *
 *	(const pqf *)HQF	These are QF structs used for low-pass and
 *	(const pqf *)GQF	  high-pass convolution-decimation.
 *
 *	(real *)work		This is a scratch array to hold intermediate
 *				  wavelet packet amplitudes.
 *
 * Outputs:
 *	(int)dwpap2hedger	The return value is the next unoccupied
 *				  block index.
 *
 * Assumptions:
 *	1. `n' is divisible by `1<<s'.
 *
 * External functions called:
 *	cdpe()
 */
extern int
  dwpap2hedger(
	       hedge *graph,	/* Partially specified data struct. */
	       int j,		/* The current block index.         */
	       int n,		/* The current block length.        */
	       int s,		/* The current scale index.         */
	       const pqf *HQF,	/* Low-pass QF data structure.      */
	       const pqf *GQF,	/* High-pass QF data structure.     */
	       real *work)	/* Intermediate wavelet amplitudes. */
{
  int i;
  real *parent;

  assert( (n>>s)<<s == n );

  if(graph->levels[j]==s)
    {
      j += 1;
    }
  else
    {
      parent = (real *)graph->contents[j];
      cdpe(  work, 1, parent, n, HQF );
      cdpe( work+n/2, 1, parent, n, GQF );
      for(i=0; i<n; i++) parent[i] = work[i];
      j = dwpap2hedger( graph, j, n/2, s+1, HQF, GQF, work );
      graph->contents[j] = (void *)(parent + n/2);
      j = dwpap2hedger( graph, j, n/2, s+1, HQF, GQF, work );
    }
  return(j);
}

/*********************************************************
 * dwpap2hedge()
 *
 * In place [D]iscrete [W]avelet [P]acket [A]nalysis,
 * [P]eriodic, to a [HEDGE], ab initio
 *
 * Calling sequence and basic algorithm:
 *
 *   dwpap2hedge( IN, LENGTH, LEVELS, BLOCKS, HQF, GQF ):
 *     Let GRAPH = makehedge( BLOCKS, NULL, LEVELS, NULL )
 *     Let GRAPH.CONTENTS[0] = IN
 *     Allocate an array of LENGTH REALs at WORK
 *     dwpap2hedger( GRAPH, 0, LENGTH, 0, HQF, GQF, WORK )
 *     Deallocate WORK[]
 *     Return GRAPH
 *
 * Inputs:
 *	(real *)in		This is the input array
 *
 *	(int)length		This is the positive length of `in[]'
 *
 *	(unsigned char *)levels		This is the list of levels defining
 *				  the graph basis.
 *
 *	(int)blocks		This is the number of blocks in the graph
 *				  basis defined by `levels[]'.
 *
 *	(const pqf *)HQF	These are QF structs used for low-pass and
 *	(const pqf *)GQF	  high-pass convolution-decimation.
 *
 * Outputs:
 *	(hedge *)dwpap2hedge	The return value is a HEDGE data structure
 *				  pointing to various offsets of `in[]',
 *				  which has been been replaced with wavelet
 *				  packet amplitudes.
 *
 * Assumptions:
 *	1. `levels[]' is a valid graph basis description.
 *
 * External functions called:
 *	makehedge(), assert(), calloc(), dwpap2hedger(), free()
 */
extern hedge *
  dwpap2hedge(
	      real *in,		/* Joint input and output array.    */
	      int length,	/* Number of elements of `in[]'.    */
	      unsigned char *levels, /* Graph basis levels list.    */
	      int blocks,	/* Number of blocks in the basis.   */
	      const pqf *HQF,	/* Low-pass QF data structure.      */
	      const pqf *GQF)	/* High-pass QF data structure.     */
{
  hedge *graph;
  real *work;

  graph = makehedge( blocks, 0, levels, 0 );
  graph->contents[0] = (void *)in;
  work = (real *)calloc(length, sizeof(real));  assert(work);
  dwpap2hedger( graph, 0, length, 0, HQF, GQF, work );
  free(work);
  return(graph);
}

/*********************************************************
 * abt2dwpsp()
 *
 * [A]rray [B]inary [T]ree to [D]iscrete [W]avelet [P]acket
 * [S]ynthesis, [P]eriodic
 *
 * Calling sequence and basic algorithm:
 *
 *   abt2dwpsp( DATA, N, MAXLEVEL, HQF, GQF ):
 *      Let L = MAXLEVEL
 *      While L>0
 *         Let NCHILD = abtblength( N, L )
 *         L -= 1
 *         For B = 0 to (1<<L)-1
 *            Let DPARENT = DATA + abtblock( N, L, B )
 *            Let DCHILD = DATA + abtblock( N, L+1, 2*B )
 *            acdpi( DPARENT, 1, DCHILD, NCHILD, HQF )
 *            Let DCHILD = DATA + abtblock( N, L+1, 2*B+1 )
 *            acdpi( DPARENT, 1, DCHILD, NCHILD, GQF )
 *
 * Inputs:
 *	(real *)data		This array binary contains the amplitudes of
 *				  the wavelet packets to superpose.
 *
 *	(int)n			This is the positive row length of the ABT.
 *
 *	(int)maxlevel		This is depth of the ABT.
 *
 *	(const pqf *)HQF	These are QF structs used for low-pass and
 *	(const pqf *)GQF	  high-pass convolution-decimation.
 *
 * Outputs:
 *	The first `n'-element row of the array binary tree at `data[]'
 *	has the wavelet packets constructed from all deeper rows
 *	superposed into itself.  Intermediate depth rows also have partial
 *	reconstructions superposed onto themselves.
 *
 * Assumptions:
 *	1. `n' is divisible by `1<<maxlevel'.
 *
 * External functions called:
 *	acdpi(), abtblock(), abtblength(), assert()
 */
extern void
  abt2dwpsp(
	    real *data,		/* Input/output array binary tree.  */
	    int n,		/* Row length of the ABT `data[]'.  */
	    int maxlevel,	/* Depth of the ABT `data[]'.       */
	    const pqf *HQF,	/* Low-pass QF data structure.      */
	    const pqf *GQF)	/* High-pass QF data structure.     */
{
  int level, block, nchild;
  real *dchild, *dparent;

  assert( (n>>maxlevel)<<maxlevel == n );

  level = maxlevel;
  while( level>0 )
    {
      nchild = abtblength( n, level );
      --level;
      for( block=0; block<(1<<level); block++ )
	{
	  dparent = data + abtblock( n, level, block );
	  /* Left child: */
	  dchild  = data + abtblock( n, level+1, 2*block );
	  acdpi( dparent, 1, dchild, nchild, HQF );
	  /* Right child: */
	  dchild  = data + abtblock( n, level+1, 2*block+1 );
	  acdpi( dparent, 1, dchild, nchild, GQF );
	}
    }
  return;
}

/*********************************************************
 * tfa1s2dwpsp()
 *
 * From [T]ime [F]requency [A]toms to the [D]iscrete [W]avelet
 * [P]acket [S]ynthesis, [P]eriodic
 *
 * Calling sequence and basic algorithm:
 *
 *   tfa1s2dwpsp( ATOMS, NUM, N, HQF, GQF ):
 *      Let MAXLEVEL = ATOMS[0].LEVEL
 *      For I = 1 to NUM-1
 *         Let MAXLEVEL = max( MAXLEVEL, ATOMS[I].LEVEL )
 *      Allocate an array of N*(MAXLEVEL+1) 0s at DATA
 *      tfa1s2abt( DATA, N, ATOMS, NUM )
 *      abt2dwpsp( DATA, N, MAXLEVEL, HQF, GQF )
 *      Return DATA
 *
 * Inputs:
 *	(tfa1 *)atoms		This is the input array
 *
 *	(int)num		This is the number of elements in `atoms[]'
 *
 *	(int)n			This is the period of the original signal
 *				  represented by the atoms.
 *
 *	(const pqf *)HQF	These are QF structs used for low-pass and
 *	(const pqf *)GQF	  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)dwpap2hedge	The return value is an array binary tree
 *				  whose first row contains the superposition
 *				  of the wavelet packets described by the
 *				  array `atoms[]'.
 *
 * Assumptions:
 *	1. `atoms[]' contains a valid list of TFA1s.
 *
 * External functions called:
 *	max(), calloc(), assert(), tfa1s2abt(), abt2dwpsp()
 */
extern real *
  tfa1s2dwpsp(
	      tfa1 *atoms,	/* Input array of time-frequency atoms.    */
	      int num,		/* Number of elements of `atoms[]'.        */
	      int n,		/* Period of the signal to be synthesized. */
	      const pqf *HQF,	/* Periodized low-pass QF data structure.  */
	      const pqf *GQF)	/* Periodized high-pass QF data structure. */
{
  int maxlevel, i;
  real *data;

  /* Find the depth of the analysis tree: */
  maxlevel = atoms[0].level;
  for( i=1; i<num; i++ )
    maxlevel = max( maxlevel, atoms[i].level );

  /* Allocate an array binary tree, fill it with amplitudes, and superpose: */
  data = (real *)calloc(n*(1+maxlevel), sizeof(real));  assert(data);
  tfa1s2abt( data, n, atoms, num );
  abt2dwpsp( data, n, maxlevel, HQF, GQF );

  return(data);
}

/*********************************************************
 * hedge2dwpspr()
 *
 * From [HEDGE] to [D]iscrete [W]avelet [P]acket [S]ynthesis,
 * [P]eriodic, in place; [R]ecursive core.
 *
 * Calling sequence and basic algorithm:
 *
 *   hedge2dwpspr( GRAPH, J, N, S, HQF, GQF, WORK ):
 *     If S < GRAPH.LEVELS[J] then
 *       Let LEFT = GRAPH.CONTENTS[J]
 *       Let J = hedge2dwpspr( GRAPH, J, N/2, S+1, HQF, GQF, WORK )
 *       Let RIGHT = GRAPH.CONTENTS[J]
 *       Let J = hedge2dwpspr( GRAPH, J, N/2, S+1, HQF, GQF, WORK )
 *       acdpe( WORK, 1, LEFT, N/2, HQF ) 
 *       acdpo( WORK, 1, RIGHT, N/2, GQF ) 
 *       For I = 0 to N-1
 *          Let LEFT[I] = WORK[I]
 *     Else
 *       J += 1
 *     Return J
 *
 * Inputs:
 *	(hedge *)graph		This data structure should be a properly
 *				  allocated and assigned graph basis.
 *
 *	(int)j			This is the current block index.
 *
 *	(int)n			This is the positive signal length.
 *
 *	(int)s			This is the current scale index.
 *
 *	(const pqf *)HQF	These are QF structs used for low-pass and
 *	(const pqf *)GQF	  high-pass convolution-decimation.
 *
 *	(real *)work		This is a scratch array to hold intermediate
 *				  wavelet packet amplitudes.
 *
 * Outputs:
 *	(int)hedge2dwpspr	The return value is the next unreconstructed
 *				  block index.
 *
 * External functions called:
 *	acdpe(), acdpo()
 */
extern int
  hedge2dwpspr(
	       hedge *graph,	/* Input data structure.            */
	       int j,		/* The current block index.         */
	       int n,		/* The current block length.        */
	       int s,		/* The current scale index.         */
	       const pqf *HQF,	/* Low-pass QF data structure.      */
	       const pqf *GQF,	/* High-pass QF data structure.     */
	       real *work)	/* Intermediate wavelet amplitudes. */
{
  int i;
  real *left, *right;

  if( s < graph->levels[j] )
    {
      left = (real *)graph->contents[j];
      j = hedge2dwpspr( graph, j, n/2, s+1, HQF, GQF, work );
      right = (real *)graph->contents[j];
      j = hedge2dwpspr( graph, j, n/2, s+1, HQF, GQF, work );
      acdpe( work, 1, left, n/2, HQF );
      acdpo( work, 1, right, n/2, GQF );
      for( i=0; i<n; i++ ) left[i] = work[i];
    }
  else
    {
      j += 1;
    }
  return(j);
}

/*********************************************************
 * hedge2dwpsp()
 *
 * From [HEDGE] to [D]iscrete [W]avelet [P]acket [S]ynthesis,
 * [P]eriodic, in place; ab initio
 *
 * Calling sequence and basic algorithm:
 *
 *   hedge2dwpsp( GRAPH, LENGTH, HQF, GQF ):
 *     Allocate an array of LENGTH REALs at WORK
 *     hedge2dwpspr( GRAPH, 0, LENGTH, 0, HQF, GQF, WORK )
 *     Deallocate WORK[]
 *     Return GRAPH.CONTENTS[0]
 *
 * Inputs:
 *	(hedge *)graph		This data structure should be a properly
 *				  allocated and assigned graph basis.
 *
 *	(int)length		This is the period of the signal.
 *
 *	(const pqf *)HQF	These are QF structs used for low-pass and
 *	(const pqf *)GQF	  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)hedge2dwpsp	The return value is a pointer to the start
 *				  of the reconstructed signal array.
 *
 * Assumptions:
 *	1. `graph' is a valid graph basis for period `length'.
 *
 * External functions called:
 *	hedge2dwpspr(), calloc(), assert(), free()
 */
extern real *
  hedge2dwpsp(
	      hedge *graph,	/* Input data structure.                 */
	      int length,	/* Period of the synthesized signal.     */
	      const pqf *HQF,	/* Low-pass periodic QF data structure.  */
	      const pqf *GQF)	/* High-pass periodic QF data structure. */
{
  real *work;

  work = (real *)calloc( length, sizeof(real));  assert(work);
  hedge2dwpspr( graph, 0, length, 0, HQF, GQF, work );
  free(work);
  return((real *)graph->contents[0]);
}

/*********************************************************
 * hedge2dwpspabt()
 *
 * From a [HEDGE] to a [D]iscrete [W]avelet [P]acket
 * [S]ythesis, [P]eriodic [A]rray [B]inary [T]ree, in natural order
 *
 * Calling sequence and basic algorithm:
 *
 *   hedge2dwpspabt( GRAPH, N, HQF, GQF ):
 *      Let MAXLEVEL = GRAPH.LEVELS[0]
 *      For I = 1 to GRAPH.BLOCKS-1
 *         Let MAXLEVEL = max( MAXLEVEL, GRAPH.LEVELS[I] )
 *      Allocate an array of N*(MAXLEVEL+1) 0s at DATA
 *      hedge2abt( DATA, GRAPH, N )
 *      abt2dwpsp( DATA, N, MAXLEVEL, HQF, GQF )
 *      Return DATA
 *
 * Inputs:
 *	(hedge *)graph		This defines the  wavelet packet basis.
 *
 *	(int)n			This is the period of the original signal.
 *
 *	(const pqf *)HQF	These are QF structs used for low-pass and
 *	(const pqf *)GQF	  high-pass convolution-decimation.
 *
 * Outputs:
 *	(real *)hedge2dwpspabt	Return an array binary tree containing the
 *				  superposed wavelet packets defined by
 *				  `graph' in its first row.
 *
 * Assumptions:
 *	1. `graph' is a valid graph basis for signals of length `n'.
 *
 * External functions called:
 *	max(), hedge2abt(), abt2dwpsp(), assert(), calloc()
 */
extern real *
  hedge2dwpspabt(
		 hedge *graph,	/* Basis amplitudes to reconstruct.    */
		 int n,		/* Period of the reconstructed signal. */
		 const pqf *HQF, /* Low-pass QF data structure.        */
		 const pqf *GQF) /* High-pass QF data structure.       */
{
  real *data;
  int maxlevel, i;

  /* Find the depth of the tree necessary to hold the atoms: */
  maxlevel = graph->levels[0];
  for( i=1; i<graph->blocks; i++ )
    maxlevel = max( maxlevel, graph->levels[i] );

  data = (real *)calloc(n*(maxlevel+1), sizeof(real));  assert(data);
  hedge2abt( data, graph, n );
  abt2dwpsp( data, n, maxlevel, HQF, GQF );
  return(data);  
}

/*********************************************************
 * cdachild()
 *
 * Allocate and superpose amplitudes into a child block by
 * aperiodic convolution-decimation.
 *
 * Calling sequence and basic algorithm:
 *
 *   cdachild( CHILD, PARENT, F ):
 *      If PARENT != NULL then
 *         If PARENT.ORIGIN != NULL 
 *            Let CLEAST = cdaleast( PARENT, F )
 *            Let CFINAL = cdafinal( PARENT, F )
 *            Let CHILD = enlargeinterval( CHILD, CLEAST, CFINAL )
 *            cdai( CHILD.ORIGIN, 1, PARENT.ORIGIN,
 *                                 PARENT.LEAST, PARENT.FINAL, F )
 *      Return CHILD
 *
 * Inputs:
 *	(interval *)child	This data structure should be properly
 *				  allocated to receive the output.
 *
 *	(interval *)parent	This data structure should be a properly
 *				  allocated and assigned INTERVAL.
 *
 *	(const pqf *)F		This is the PQF data structure used for
 *				  aperiodic convolution-decimation.
 *
 * Outputs:
 *	(interval *)child	The return value is a pointer to the possibly
 *				  enlarged child interval.
 *
 * External functions called:
 *	enlargeinterval(), cdaleast(), cdafinal(), cdai()
 */
extern interval *
  cdachild(
	   interval *child,	/* Output INTERVAL data structure. */
	   interval *parent,	/* Input INTERVAL data structure.  */
	   const pqf *F)	/* Aperiodic QF data structure.    */
{
  int cleast, cfinal;

  if( parent )
    if( parent->origin )
      {
	cleast = cdaleast( parent, F );
	cfinal = cdafinal( parent, F );
	child = enlargeinterval( child, cleast, cfinal );
	cdai(child->origin, 1, parent->origin, parent->least,parent->final, F);
      }

  return(child);
}

/*********************************************************
 * dwpaa2btntr()
 *
 * [D]iscrete [W]avelet [P]acket [A]nalysis,
 * [A]periodic, to a [BTN] [T]ree; [R]ecursive core.
 *
 * Calling sequence and basic algorithm:
 *
 *   dwpaa2btntr( NODE, LEVEL, HQF, GQF ):
 *      If LEVEL>0 then
 *         Let CHILD = cdachild( NULL, NODE.CONTENT, HQF )
 *         Let NODE.LEFT = makebtn( CHILD, NULL, NULL, NULL )
 *         dwpaa2btntr( NODE.LEFT, LEVEL-1, HQF, GQF )
 *         Let CHILD = cdachild( NULL, NODE.CONTENT, GQF )
 *         Let NODE.RIGHT = makebtn( CHILD, NULL, NULL, NULL )
 *         dwpaa2btntr( NODE.RIGHT, LEVEL-1, HQF, GQF )
 *      Return
 *
 * Inputs:
 *	(btn *)node		This is the binary tree node to develop
 *
 *	(int)level		This is the depth of the subtree to develop.
 *
 *	(const pqf *)HQF	These are QF structs used for low-pass and
 *	(const pqf *)GQF	  high-pass convolution-decimation.
 *
 * Outputs:
 *	A binary tree of BTN data structures is built below `node' to
 *	  a depth of `level'.
 *
 * Assumptions:
 *	1. level >= 0
 *	2. node != NULL
 *
 * External functions called:
 *	cdachild(), makebtn(), assert()
 */
extern void
  dwpaa2btntr(
	      btn *node,	/* Root of the current BTN tree. */
	      int level,	/* Depth of the tree to develop. */
	      const pqf *HQF,	/* Low-pass QF data structure.   */
	      const pqf *GQF)	/* High-pass QF data structure.  */
{
  interval *parent, *child;

  assert(level>=0);
  assert(node);

  if( level>0 )
    {
      parent = (interval *)node->content;
      child = cdachild( 0, parent, HQF );
      node->left = makebtn( child, 0, 0, 0 );
      dwpaa2btntr( node->left, level-1, HQF, GQF );
      child = cdachild( 0, parent, GQF );
      node->right = makebtn( child, 0, 0, 0 );
      dwpaa2btntr( node->right, level-1, HQF, GQF );
    }
  return;
}

/*********************************************************
 * dwpaa2btnt()
 *
 * Complete [D]iscrete [W]avelet [P]acket [A]nalysis,
 * [A]periodic, to a [BTN] [T]ree, from an array.
 *
 * Calling sequence and basic algorithm:
 *
 *   dwpaa2btnt( IN, LEAST, FINAL, MAXLEVEL, HQF, GQF ):
 *      Let ROOT = makebtn( NULL, NULL, NULL, NULL )
 *      Let ROOT.CONTENT = makeinterval( IN, LEAST, FINAL )
 *      dwpaa2btntr( ROOT, MAXLEVEL, HQF, GQF )
 *      Return ROOT
 *
 * Inputs:
 *	(real *)in		This is the input signal.
 *
 *	(int)least		This is the least valid index of `in[]'.
 *
 *	(int)final		This is the greatest valid index of `in[]'.
 *
 *	(int)maxlevel		This is the depth of the wavelet packet tree.
 *
 *	(const pqf *)HQF	These are QF structs used for low-pass and
 *	(const pqf *)GQF	  high-pass convolution-decimation.
 *
 * Outputs:
 *	(btn *)dwpaa2btnt	A binary tree of BTN data structures is built
 *				  by side effect; the return value is its root.
 *
 * Assumptions:
 *	1. maxlevel >= 0
 *
 * External functions called:
 *	dwpaa2btntr(), makebtn(), makeinterval(), assert()
 */
extern btn *
  dwpaa2btnt(
	     real *in,		/* Input signal array. */
	     int least,		/* Least valid index of `in[]'. */
	     int final,		/* Greatest valid index of `in[]'. */
	     int maxlevel,	/* Depth of the tree to develop. */
	     const pqf *HQF,	/* Low-pass QF data structure.   */
	     const pqf *GQF)	/* High-pass QF data structure.  */
{
  btn *root;

  assert(maxlevel>=0);
  
  root = makebtn(0, 0, 0, 0);
  root->content = (void *)makeinterval( in, least, final );
  dwpaa2btntr( root, maxlevel, HQF, GQF );
  return(root);
}

/*********************************************************
 * dwpaa2hedger()
 *
 * [D]iscrete [W]avelet [P]acket [A]nalysis,
 * [A]periodic, to a [HEDGE]; [R]ecursive core.
 *
 * Calling sequence and basic algorithm:
 *
 *   dwpaa2hedger( GRAPH, NODE, J, S, HQF, GQF ):
 *     If GRAPH.LEVELS[J]==S then
 *       Let GRAPH.CONTENTS[J] = NODE.CONTENT
 *       Let NODE.CONTENT = NULL
 *       J += 1
 *     Else
 *       Let CHILD = cdachild( NULL, NODE.CONTENT, HQF )
 *       Let NODE.LEFT = makebtn( CHILD, NULL, NULL, NULL )
 *       Let J = dwpaa2hedger( GRAPH, NODE.LEFT, J, S+1, HQF, GQF )
 *       Let CHILD = cdachild( NULL, NODE.CONTENT, GQF )
 *       Let NODE.RIGHT = makebtn( CHILD, NULL, NULL, NULL )
 *       Let J = dwpaa2hedger( GRAPH, NODE.RIGHT, J, S+1, HQF, GQF )
 *     Return J
 *
 * Inputs:
 *	(hedge *)graph		This contains a valid levels list.
 *
 *	(btn *)node		This is the root of the current BTN tree.
 *
 *	(int)j			This is the current block number.
 *
 *	(int)s			This is the current scale index.
 *
 *	(const pqf *)HQF	These are QF structs used for low-pass and
 *	(const pqf *)GQF	  high-pass convolution-decimation.
 *
 * Outputs:
 *	(int)dwpaa2hedger	The return value is the number of the next
 *				  block to decompose.
 *
 * External functions called:
 *	cdachild(), makebtn(), assert()
 */
extern int
  dwpaa2hedger(
	       hedge *graph,	/* Level-specified graph basis.   */
	       btn *node,	/* Root of the current BTN tree.  */
	       int j,		/* Current block number.          */
	       int s,		/* Current level in the BTN tree. */
	       const pqf *HQF,	/* Low-pass QF data structure.    */
	       const pqf *GQF)	/* High-pass QF data structure.   */
{
  interval *parent, *child;

  if( graph->levels[j]==s )
    {
      graph->contents[j] = node->content;
      node->content = 0;	/* prevent this INTERVAL from being free'd */
      ++j;
    }
  else
    {
      parent = (interval *)node->content;
      child = cdachild( 0, parent, HQF );
      node->left = makebtn( child, 0, 0, 0 );
      j = dwpaa2hedger( graph, node->left, j, s+1, HQF, GQF );
      child = cdachild( 0, parent, GQF );
      node->right = makebtn( child, 0, 0, 0 );
      j = dwpaa2hedger( graph, node->right, j, s+1, HQF, GQF );
    }
  return(j);
}

/*********************************************************
 * dwpaa2hedge()
 *
 * Complete [D]iscrete [W]avelet [P]acket [A]nalysis,
 * [A]periodic, to a [HEDGE], from an array.
 *
 * Calling sequence and basic algorithm:
 *
 *   dwpaa2hedge( LEVELS, BLOCKS, DATA, LEAST, FINAL, HQF, GQF ):
 *      Let ROOT = makebtn( NULL, NULL, NULL, NULL )
 *      Let ROOT.CONTENT = makeinterval( DATA, LEAST, FINAL )
 *      Let GRAPH = makehedge( BLOCKS, NULL, LEVELS, NULL )
 *      dwpaa2hedger( GRAPH, ROOT, 0, HQF, GQF )
 *      Let ROOT.CONTENT = NULL
 *      freebtnt( ROOT, freeinterval, free )
 *      Return GRAPH
 *
 * Inputs:
 *	(unsigned char *)levels	This specifies the graph basis to develop.
 *
 *	(int)blocks		This is the positive length of `levels[]'
 *
 *	(real *)data		This is the input signal.
 *
 *	(int)least		This is the least valid index of `data[]'.
 *
 *	(int)final		This is the greatest valid index of `data[]'.
 *
 *	(const pqf *)HQF	These are QF structs used for low-pass and
 *	(const pqf *)GQF	  high-pass convolution-decimation.
 *
 * Outputs:
 *	(hedge *)dwpaa2hedge	A binary tree of BTN data structures is built
 *				  and deallocated by side effect, but the
 *				  INTERVAL content of the nodes specified by
 *				  `levels[]' are put into the output HEDGE.
 *
 * Assumptions:
 *	1. `levels[]' is a valid graph basis description.
 *
 * External functions called:
 *	
 */
extern hedge *
  dwpaa2hedge(
	      unsigned char *levels, /* Levels specifying the basis. */
	      int blocks,	/* Number of blocks in the basis.    */
	      real *data,	/* Input signal array.               */
	      int least,	/* Least valid index of `data[]'.    */
	      int final,	/* Greatest valid index of `data[]'. */
	      const pqf *HQF,	/* Low-pass QF data structure.       */
	      const pqf *GQF)	/* High-pass QF data structure.      */
{
  hedge *graph;
  btn *root;

  root = makebtn( 0, 0, 0, 0 );
  root->content = makeinterval( data, least, final );
  graph = makehedge( blocks, 0, levels, 0 );
  dwpaa2hedger( graph, root, 0, 0, HQF, GQF );
  root->content = 0;
  freebtnt( root, (freetype)freeinterval, freevoid );
  return(graph);
}

/*********************************************************
 * acdaparent()
 *
 * Allocate and superpose amplitudes into a parent block by
 * adjoint aperiodic convolution-decimation.
 *
 * Calling sequence and basic algorithm:
 *
 *   acdaparent( PARENT, CHILD, F ):
 *      If CHILD != NULL then
 *         If CHILD.ORIGIN != NULL then
 *            Let LEAST = acdaleast( CHILD, F )
 *            Let LEAST = min( PARENT.LEAST, LEAST )
 *            Let FINAL = acdafinal( CHILD, F )
 *            Let FINAL = max( PARENT.FINAL, FINAL )
 *            Let PARENT = enlargeinterval( PARENT, LEAST, FINAL )
 *            acdai( PARENT.ORIGIN, 1, CHILD.ORIGIN, 
 *                                    CHILD.LEAST, CHILD.FINAL, F )
 *      Return PARENT
 *
 * Inputs:
 *	(interval *)parent	This data structure should be properly
 *				  allocated to receive the output.
 *
 *	(interval *)child	This data structure should be a properly
 *				  allocated and assigned INTERVAL.
 *
 *	(const pqf *)F		This is the PQF data structure used for
 *				  aperiodic adjoint convolution-decimation.
 *
 * Outputs:
 *	(interval *)parent	The return value is a pointer to the possibly
 *				  enlarged parent interval.
 *
 * External functions called:
 *	enlargeinterval(), acdaleast(), acdafinal(), acdai(), max(), min()
 */
extern interval *
  acdaparent(
	     interval *parent,	/* Output INTERVAL data structure. */
	     interval *child,	/* Input INTERVAL data structure.  */
	     const pqf *F)	/* Aperiodic QF data structure.    */
{
  if( child )
    if( child->origin )
      {
	int least, final;
	least = acdaleast( child, F );
	least = min( parent->least, least );
	final = acdafinal( child, F );
	final = max( parent->final, final );
	parent = enlargeinterval( parent, least, final );
	acdai(parent->origin, 1, child->origin, child->least,child->final, F);
      }
  return(parent);
}

/*********************************************************
 * btnt2dwpsa()
 *
 * From a [BTN] [T]ree to a [D]iscrete [W]avelet [P]acket
 * [S]ythesis, [A]periodic, natural order
 *
 * Calling sequence and basic algorithm:
 *
 *   btnt2dwpsa( ROOT, HQF, GQF ):
 *      If ROOT != NULL then
 *         btnt2dwpsa( ROOT.LEFT, HQF, GQF )
 *         If ROOT.LEFT != NULL then
 *            Let ROOT.CONTENT = acdaparent( ROOT.CONTENT,
 *                                        ROOT.LEFT.CONTENT, HQF )
 *         btnt2dwpsa( ROOT.RIGHT, HQF, GQF )
 *         If ROOT.RIGHT != NULL then
 *            Let ROOT.CONTENT = acdaparent( ROOT.CONTENT,
 *                                       ROOT.RIGHT.CONTENT, GQF )
 *     Return
 *
 * Inputs:
 *	(btn *)root		This is the root of the tree to synthesize.
 *
 *	(const pqf *)HQF	These are QF structs used for low-pass and
 *	(const pqf *)GQF	  high-pass convolution-decimation.
 *
 * Outputs:
 *	The binary tree of BTN data structures at `root' is traversed and
 *	wavelet packets are superposed into its nodes by side effect.
 *
 * Assumptions:
 *	1. Leaf nodes have NULL left and right members.
 *
 * External functions called:
 *	
 */
extern void
  btnt2dwpsa(
	     btn *root,		/* Root of the tree to reconstruct. */
	     const pqf *HQF,	/* Low-pass QF data structure.      */
	     const pqf *GQF)	/* High-pass QF data structure.     */
{
  if( root )
    {
      btnt2dwpsa( root->left, HQF, GQF );
      if( root->left )
	root->content = 
	  (void *)acdaparent((interval *)root->content,
			     (interval *)root->left->content, HQF );
      btnt2dwpsa( root->right, HQF, GQF );
      if( root->right )
	root->content =
	  (void *)acdaparent((interval *)root->content,
                             (interval *)root->right->content, GQF );
    }
  return;
}

/*********************************************************
 * btnt2dwpsa0()
 *
 * From a [BTN] [T]ree to a [D]iscrete [W]avelet [P]acket
 * [S]ythesis, [A]periodic, with immediate deallocation of
 * intermediate nodes, in natural order
 *
 * Calling sequence and basic algorithm:
 *
 *   btnt2dwpsa0( ROOT, HQF, GQF ):
 *     If ROOT != NULL then
 *       btnt2dwpsa0( ROOT.LEFT, HQF, GQF )
 *       If ROOT.LEFT != NULL then
 *         Let ROOT.CONTENT = acdaparent( ROOT.CONTENT,
 *                                     ROOT.LEFT.CONTENT, HQF )
 *         Let ROOT.LEFT = freebtn( ROOT.LEFT, freeinterval, free )
 *       btnt2dwpsa0( ROOT.RIGHT, HQF, GQF )
 *       If ROOT.RIGHT != NULL then
 *         Let ROOT.CONTENT = acdaparent( ROOT.CONTENT,
 *                                     ROOT.RIGHT.CONTENT, GQF )
 *         Let ROOT.RIGHT = freebtn( ROOT.RIGHT, freeinterval, free )
 *     Return
 *
 * Inputs:
 *	(btn *)root		This is the root of the tree to synthesize.
 *
 *	(const pqf *)HQF	These are QF structs used for low-pass and
 *	(const pqf *)GQF	  high-pass convolution-decimation.
 *
 * Outputs:
 *	The binary tree of BTN data structures at `root' is traversed and
 *	wavelet packets are superposed into its nodes by side effect.
 *
 * Assumptions:
 *	1. Leaf nodes have NULL left and right members.
 *
 * External functions called:
 *	acdaparent(), freebtn()
 */
extern void
  btnt2dwpsa0(
	      btn *root,	/* Root of the tree to reconstruct. */
	      const pqf *HQF,	/* Low-pass QF data structure.      */
	      const pqf *GQF)	/* High-pass QF data structure.     */
{
  if( root )
    {
      btnt2dwpsa0( root->left, HQF, GQF );
      if( root->left )
	{
	  root->content = 
	    (void *)acdaparent((interval *)root->content,
			       (interval *)root->left->content, HQF );
	  root->left = freebtn(root->left, (freetype)freeinterval, freevoid);
	}
      btnt2dwpsa0( root->right, HQF, GQF );
      if( root->right )
	{
	  root->content =
	    (void *)acdaparent((interval *)root->content,
			       (interval *)root->right->content, GQF );
	  root->right = freebtn(root->right, (freetype)freeinterval, freevoid);
	}
    }
  return;
}

/*********************************************************
 * tfa1s2dwpsa()
 *
 * From a list of [TFA1s] to a [D]iscrete [W]avelet [P]acket
 * [S]ythesis, [A]periodic, in natural order
 *
 * Calling sequence and basic algorithm:
 *
 *   tfa1s2dwpsa( ATOMS, NUM, HQF, GQF ):
 *      Let ROOT = makebtn( NULL,NULL, NULL, NULL )
 *      tfa1s2btnt( ROOT, ATOMS, NUM )
 *      btnt2dwpsa( ROOT, HQF, GQF )
 *      Let OUT = ROOT.CONTENT
 *      Let ROOT.CONTENT = NULL
 *      freebtnt( ROOT, FREEINTERVAL, FREE )
 *      Return OUT
 *
 * Inputs:
 *	(tfa1 *)atoms		These are the wavelet packets to superpose.
 *
 *	(int)num		This is the length of the list `atoms[]'.
 *
 *	(const pqf *)HQF	These are QF structs used for low-pass and
 *	(const pqf *)GQF	  high-pass convolution-decimation.
 *
 * Outputs:
 *	(interval *)tfa1s2dwpsa	Return an interval containing the superposed
 *				  wavelet packets defined by `atoms[]'.
 *
 * Assumptions:
 *	1. Elements of `atoms[]' have valid indices.
 *
 * External functions called:
 *	makebtn(), tfa1s2btnt(), btnt2dwpsa(), 
 *	freebtnt(), freeinterval(), freevoid()
 */
extern interval *
  tfa1s2dwpsa(
	      tfa1 *atoms,	/* List of wavelet packets to reconstruct. */
	      int num,		/* Number of wavelet packets in `atoms[]'. */
	      const pqf *HQF,	/* Low-pass QF data structure.      */
	      const pqf *GQF)	/* High-pass QF data structure.     */
{
  btn *root;
  interval *out;

  root = makebtn(0, 0, 0, 0);
  tfa1s2btnt( root, atoms, num );
  btnt2dwpsa( root, HQF, GQF );
  out = (interval *)root->content;
  root->content = 0;
  freebtnt( root, (freetype)freeinterval, freevoid );
  return(out);
}

/*********************************************************
 * hedge2dwpsa()
 *
 * From a [HEDGE] to a [D]iscrete [W]avelet [P]acket
 * [S]ythesis, [A]periodic, in natural order
 *
 * Calling sequence and basic algorithm:
 *
 *   hedge2dwpsa( GRAPH, HQF, GQF ):
 *      Let ROOT = makebtn( NULL,NULL, NULL, NULL )
 *      hedge2btnt( ROOT, GRAPH )
 *      btnt2dwpsa( ROOT, HQF, GQF )
 *      Let OUT = ROOT.CONTENT
 *      Let ROOT.CONTENT = NULL
 *      freebtnt( ROOT, freeinterval, free )
 *      Return OUT
 *
 * Inputs:
 *	(hedge *)graph		This defines the  wavelet packet basis.
 *
 *	(const pqf *)HQF	These are QF structs used for low-pass and
 *	(const pqf *)GQF	  high-pass convolution-decimation.
 *
 * Outputs:
 *	(interval *)hedge2dwpsa	Return an interval containing the superposed
 *				  wavelet packets defined by `graph'.
 *
 * Assumptions:
 *	1. `graph' is a valid graph basis.
 *
 * External functions called:
 *	makebtn(), hedge2btnt(), btnt2dwpsa(), 
 *	freebtnt(), freeinterval(), freevoid()
 */
extern interval *
  hedge2dwpsa(
	      hedge *graph,	/* Basis amplitudes to reconstruct. */
	      const pqf *HQF,	/* Low-pass QF data structure.      */
	      const pqf *GQF)	/* High-pass QF data structure.     */
{
  btn *root;
  interval *out;

  root = makebtn(0, 0, 0, 0);
  hedge2btnt( root, graph );
  btnt2dwpsa( root, HQF, GQF );
  out = (interval *)root->content;
  root->content = 0;
  freebtnt( root, (freetype)freeinterval, freevoid );
  return(out);
}
