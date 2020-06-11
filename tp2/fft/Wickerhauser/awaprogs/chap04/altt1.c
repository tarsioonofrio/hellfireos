/* 
 * Adapted local trigonometric transform functions.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include "altt1.h"

#include <stdlib.h>		/* for calloc(), malloc() */
#include <assert.h>		/* We abort if malloc() returns NULL. */
#include "real.h"
#include "interval.h"		/* Rising cutoff arrays are INTERVALs.       */
#include "fold.h"		/* Declare `fips()' and `fipc()' functions.  */
#include "dtts.h"		/* Declare various DCT and DST functions.    */
#include "utility.h"		/* Declare `abt2hedge()' and `btnt2hedge()'. */



/***********************************************************************
 * initrcf()
 *
 *  [INIT]ialize [R]ising [C]utoff [F]unction.
 *
 *  Calling sequence and basic algorithm:
 *
 *    initrcf( E ):
 *       Let RISE = makeinterval(NULL, -E, E-1 )
 *       rcfmidp( RISE )
 *       Return RISE
 *
 *  Inputs:
 *	(int)e		This positive integer is the range of folding.
 *
 * Outputs:
 *	(interval *)initrcf	The return value is a pointer to a
 *				 midpoint folding function.
 *
 * External functions called:
 *	rcfmidp()	Declared in "rcf.h"
 */
extern interval *
  initrcf(
	  int e)		/* Range of midpoint folding. */
{
  interval *rise;

  assert(e>0);

  rise = makeinterval(0, -e, e-1);
  return(rise);
}

/***********************************************************************
 * lcadf()
 *
 *  [L]ocal [C]osine [A]nalysis, [D]yadic, [F]ixed folding.
 *  This transformation works in place on a preallocated array
 *  binary tree, with the signal in the first row.
 *
 *  Calling sequence and basic algorithm:
 *
 *    lcadf( PARENT, N, L, RISE ):
 *       Let NP = N
 *       For LEVEL = 0 to L-1
 *          Let NC = NP/2
 *          For PBLOCK = 1 to 1<<LEVEL
 *             Let MIDP = PARENT + NC
 *             Let CHILD = MIDP + N
 *             fdcn( CHILD, 1, MIDP, MIDP, NC, RISE )
 *             fdcp( CHILD, 1, MIDP, MIDP, NC, RISE )
 *             dctiv( PARENT, NP )
 *             PARENT += NP
 *          Let NP = NC
 *       For PBLOCK = 1 to 1<<L
 *          dctiv( PARENT, NP )
 *          PARENT += NP
 *
 *  Inputs:
 *	(real *)parent		This is the input and output array.
 *
 *	(int)n			This is a nonnegative power of 2.
 *
 *	(int)L			This is a nonnegative integer.
 *
 *	(interval *)rise	This is a sampled rising cutoff function.
 *
 * Outputs:
 *	The array binary tree at `parent' is filled with a complete
 *	dyadic local cosine analysis of the signal in the first row.
 *
 * External functions called:
 *	fdcn(),fdcp()	Declared in "fold.h"
 *	dct_iv()	Declared in "dtts.h"
 *
 * Assumptions:
 *	1. L>=0
 *	2. n>0 is a power of 2 and divisible by (1<<L).
 *	3. rise != 0
 *	4. rise->final < (n>>L)
 */
extern void
  lcadf(
	real *parent,		/* Input and output array binary tree. */
	int n,			/* Length of one row in the tree.      */
	int L,			/* Number of levels in the tree.       */
	interval *rise)		/* Sampled rising cutoff function.   */
{
  int log2n, lp, pblock, nc, np;
  real *child, *midp;

  assert( n>0 );
  assert( L>=0 );
  assert( (n>>L)<<L == n );
  assert( rise );
  assert( rise->final < (n>>L) );

  /* Compute the logarithm base 2 of the signal length `n': */
  log2n = 0;
  while( (1<<log2n) < n ) ++log2n;
  assert((1<<log2n)==n);

  /* Fill the array binary tree: */
  np = n;
  for( lp=0; lp<L; lp++ )
    {
      nc = np/2;
      for( pblock=0; pblock < (1<<lp); pblock++ )
	{
	  midp = parent + nc;
	  child = midp + n;
	  fdcn( child, 1, midp, midp, nc, rise );
	  fdcp( child, 1, midp, midp, nc, rise );
	  dct_iv( parent, log2n-lp );
	}
      parent += np;
      np = nc;
    }
  for( pblock=0; pblock<(1<<L); pblock++ )
    {
      dct_iv( parent, log2n-L );
      parent += np;
    }

  return;
}

/***********************************************************************
 * initrcfs()
 *
 *  [INIT]ialize [R]ising [C]utoff [F]unction[S], one for each level.
 *
 *  Calling sequence and basic algorithm:
 *
 *    initrcfs( N, L ):
 *       Allocate an array of L INTERVALs at RS
 *       For I = 0 to L-1
 *          Let RS[I] = initrcf( (N/2)>>I )
 *       Return RS
 *
 *  Inputs:
 *	(int)n		This positive integer is the length of the signal.
 *	(int)L		This positive integer is the number of levels.
 *
 * Outputs:
 *	(interval **)initrcfs	The return value is a pointer to the first
 *				element of an array of `L' sampled RCFs.
 *
 * External functions called:
 *	initrcf()
 */
extern interval **
  initrcfs(
	   int n,		/* Length of signal to decompose. */
	   int L)		/* Bottom decomposition level.    */
{
  interval **rs;
  int i;

  assert(n>0);
  assert(L>=0);

  rs = (interval **)calloc(1+L, sizeof(interval *));
  for( i=0; i<=L; i++ )
    {
      rs[i] = initrcf(n/2);
      n /= 2;
    }
  return(rs);
}

/***********************************************************************
 * initlcabtn()
 *
 *  [INIT]ialize for [L]ocal [C]osine [A]nalysis to a [BTN] tree.
 *
 *  Calling sequence and basic algorithm:
 *
 *   initlcabtn( IN, N, RISE )
 *      Let ROOT = makebtn( NULL, NULL, NULL, NULL )
 *      Let ROOT.CONTENT = makeinterval( IN, 0, N-1 )
 *      fipc( ROOT.CONTENT.ORIGIN+N, ROOT.CONTENT.ORIGIN, RISE )
 *      Return ROOT
 *
 *  Inputs:
 *	(const real *)in	The input signal in this array is not changed.
 *	(int)n			Expect this many input samples.
 *	(interval *)rise	Use this to fold input to the root BTN.
 *
 * Outputs:
 *	(btn *)initlcabtn	The return value is a pointer to the newly
 *			allocated root of the analysis BTN tree.
 *
 * External functions called:
 *	fipc()		Declared in "fold.h"
 *	makeinterval()	Declared in "interval.h"
 *	makebtn()	Declared in "btn.h"
 */
extern btn *
  initlcabtn(
	     const real *in,	/* Input signal */
	     int n,		/* Length of the input signal. */
	     interval *rise)	/* Sampled rising cutoff function. */
{
  btn *root;
  real *data;

  assert(n>0);
  assert( in );
  assert( rise );

  root = makebtn(0,0,0,0);
  root->content = (void *)makeinterval( in, 0, n-1 );
  data = ((interval *)root->content)->origin;
  fipc( data + n, data, rise );
  return(root);
}

/***********************************************************************
 * lcadm()
 *
 *  [L]ocal [C]osine [A]nalysis, [D]yadic, [M]ultiple folding.
 *  This transformation builds a BTN tree from a preallocated
 *  root BTN with the signal as its contents.
 *
 *  Calling sequence and basic algorithm:
 *
 *   lcadm( ROOT, S, L, RS ):
 *     Let LENGTH = 1 + ROOT.CONTENT.FINAL
 *     If S < L then
 *       Let NC = LENGTH/2
 *       Let MIDP = ROOT.CONTENT.ORIGIN + NC
 *       Let LCHILD = makeinterval( NULL, 0, NC-1 )
 *       Let ROOT.LEFT = makebtn( LCHILD, NULL, NULL, NULL )
 *       fdcn( LCHILD.ORIGIN+NC, 1, MIDP, MIDP, NC, RS[S] )
 *       lcadm( ROOT.LEFT, S+1, L, RS )
 *       Let RCHILD = makeinterval( NULL, 0, NC-1 )
 *       Let ROOT.RIGHT = makebtn( RCHILD, NULL, NULL, NULL )
 *       fdcp( RCHILD.ORIGIN, 1, MIDP, MIDP, NC, RS[S] )
 *       lcadm( ROOT.RIGHT, S+1, L, RS )
 *     dctiv( ROOT.CONTENT.ORIGIN, LENGTH )
 *     Return
 *
 *  Inputs:
 *	(btn *)root		This must be preallocated and assigned.
 *	(int)s			This nonnegative number is the current level.
 *	(int)L			This is a small nonnegative maximum level.
 *	(interval **)rs		This array holds pointers to a sampled RCFs.
 *
 * Outputs:
 *	The BTN tree at `root' is filled with a complete
 *	dyadic local cosine analysis of the signal in `root->origin[]'.
 *
 * External functions called:
 *	fdcn(),fdcp()	Declared in "fold.h"
 *	dct_iv()	Declared in "dtts.h"
 *
 * Assumptions:
 *	1. 0 <= s <= L
 *	2. rs != 0
 *	3. root->content->least == 0
 *	4. root->content->final >= 0
 */
extern void
  lcadm(
	btn *root,		/* Current root node in the BTN tree. */
	int s,			/* Current level in the tree. */
	int L,			/* Maximum level in the tree. */
	interval **rs)		/* Sampled rising cutoff functions. */
{
  int length, log2n, nc;
  real *midp;
  interval *child;

  assert( s >= 0);
  assert( s <= L);
  assert( rs );
  assert( root );
  assert( root->content );
  assert( ((interval *)root->content)->least == 0 );
  assert( ((interval *)root->content)->final >= 0 );

  length = 1+((interval *)root->content)->final;
  if( s<L )
    {
      nc = length/2;
      midp = ((interval *)root->content)->origin + nc;

      child = makeinterval( 0, 0, nc-1 );
      root->left = makebtn( child, 0, 0, 0 );
      fdcn( child->origin + nc, 1, midp, midp, nc, rs[s] );
      lcadm( root->left, s+1, L, rs );

      child = makeinterval( 0, 0, nc-1 );
      root->right = makebtn( child, 0, 0, 0 );
      fdcp( child->origin, 1, midp, midp, nc, rs[s] );
      lcadm( root->right, s+1, L, rs );
    }

  /* Compute the logarithm base 2 of the signal length: */
  log2n = 0;
  while( (1<<log2n) < length ) ++log2n;
  assert((1<<log2n)==length);

  dct_iv( ((interval *)root->content)->origin, log2n );

  return;
}

/***********************************************************************
 * lcadf2hedge()
 *
 *  Initialize a few variables, then use `lcadf()' to develop
 *  a complete array binary tree, then extract a basis into a 
 *  partially assigned HEDGE data structure.
 *
 *  Calling sequence and basic algorithm:
 *
 *   lcadf2hedge( GRAPH, IN, N, L )
 *     Let RISE = initrcf( (N>>L)/2 )
 *     Allocate an array binary tree OUT of length (1+L)*N
 *     fdcp(  OUT,  1, IN+N, IN, N, RISE )
 *     fdcn( OUT+N, 1, IN+N, IN, N, RISE )
 *     lcadf( OUT, N, L, RISE )
 *     abt2hedge( GRAPH, OUT, N )
 *     Return OUT
 *
 *  Inputs:
 *	(hedge *)graph		This must be preallocated with its `levels'
 *				  array assigned.
 *	(const real *)in	Expect a signal in the first `n' locations.
 *	(int)n			This is the positive length of the signal.
 *	(int)L			This is the small nonnegative maximum level.
 *
 * Outputs:
 *	An array binary tree is allocated and filled with a complete
 *	fixed-folding dyadic local cosine analysis of the signal in `in[]'.
 *
 *	(real *)lcadf2hedge	The return value is a pointer to the start
 *				of the array binary tree holding the analysis.
 *
 * External functions called:
 *	initrcf()
 *	calloc()	Declared in <stdlib.h>
 *	fdcn(),fdcp()	Declared in "fold.h"
 *	abt2hedge()	Declared in "utility.h"
 *	freeinterval()	Declared in "interval.h"
 *
 * Assumptions:
 *	1. `graph' is allocated and assigned `levels' and `contents' arrays.
 *	2. `in[]' is preallocated.
 *	3. n>0 is a power of 2.
 *	4. L>=0.
 */
extern real *
  lcadf2hedge(
	      hedge *graph,	/* Partially assigned HEDGE for the output. */
	      const real *in,	/* Input signal to analyze. */
	      int n,		/* Length of the input signal. */
	      int L)		/* Maximum level in the analysis tree. */
{
  interval *rise;
  int log2n;
  real *out;

  assert( graph );
  assert( graph->levels );
  assert( graph->contents );
  assert( in );
  assert( n > 0 );
  assert( L >= 0 );

  /* Compute the logarithm base 2 of the signal length: */
  log2n = 0;
  while( (1<<log2n) < n ) ++log2n;
  assert( (1<<log2n) == n );

  /* Allocate an array binary tree: */
  out = (real *)calloc( (1+L)*n, sizeof(real) );  assert( out );

  /* Allocate and assign a fixed-length sampled RCF: */
  rise = initrcf( (n>>L)/2 );

  /* Fold the input into the first `out[]' row: */
  fdcp(  out,  1, in+n, in, n/2, rise );
  fdcn( out+n, 1, in+n, in, n/2, rise );

  lcadf( out, n, L, rise );
  abt2hedge( graph, out, n );

  /* Deallocate temporary array. */
  freeinterval(rise);

  return(out);
}

/***********************************************************************
 * lcadm2hedge()
 *
 *  Initialize a few variables, then use `lcadm()' to develop
 *  a complete BTN tree, then extract a basis subset into a 
 *  partially assigned HEDGE data structure.
 *
 *  Calling sequence and basic algorithm:
 *
 *   lcadm2hedge( GRAPH, IN, N, L )
 *      Let RS = initrcfs( N, L )
 *      Let ROOT = initlcabtn( IN, N, RS[0] )
 *      lcadm( ROOT, 0, L, RS )
 *      btnt2hedge( GRAPH, ROOT )
 *      For I = 0 to L-1
 *         freeinterval( RS[I] )
 *      Deallocate the array at RS
 *      Return ROOT
 *
 *  Inputs:
 *	(hedge *)graph		This must be preallocated with its `levels'
 *				  array assigned.
 *	(const real *)in	Expect a signal in the first `n' locations.
 *	(int)n			This is the positive length of the signal.
 *	(int)L			This is the small nonnegative maximum level.
 *
 * Outputs:
 *	A complete BTN tree is allocated and filled with a complete dyadic
 *	multiple-folding local cosine analysis of the signal `in[]'.
 *
 *	(btn *)lcadm2hedge	The return value is a pointer to the root
 *				of the BTN tree holding the analysis.
 *
 * External functions called:
 *	initrcfs(),initlcabtn()
 *	btnt2hedge()		Declared in "utility.h"
 *	freeinterval()		Declared in "interval.h"
 *
 * Assumptions:
 *	1. `graph' is allocated and assigned `levels' and `contents' arrays.
 *	2. `in[]' is preallocated.
 *	3. n>0 is a power of 2.
 *	4. L>=0.
 */
extern btn *
  lcadm2hedge(
	      hedge *graph,	/* Partially assigned HEDGE for the output. */
	      const real *in,	/* Input signal to analyze. */
	      int n,		/* Length of the input signal. */
	      int L)		/* Maximum level in the analysis tree. */
{
  interval **rs;
  int i;
  btn *root;

  assert( graph );
  assert( graph->levels );
  assert( graph->contents );
  assert( in );
  assert( n > 0 );
  assert( L >= 0 );

  /* Allocate and assign a fixed-length sampled RCF: */
  rs = initrcfs( n, L );

  /* Allocate a root BTN and fold the input into it: */
  root = initlcabtn( in, n, rs[0] );

  lcadm( root, 0, L, rs );
  btnt2hedge( graph, root );

  /* Deallocate temporary arrays. */
  for(i=0; i<L; i++)
    freeinterval(rs[i]);
  free(rs);

  return(root);
}

/***********************************************************************
 * lcsdf()
 *
 *  [L]ocal [C]osine [S]ynthesis, [D]yadic, [F]ixed folding.
 *
 *  Calling sequence and basic algorithm:
 *
 *    lcsdf( GRAPH, N, RISE ):
 *       Let DATA = GRAPH.CONTENTS[0]
 *       Let NSEG = N>>GRAPH.LEVELS[0]
 *       dctiv( DATA, NSEG )
 *       For BLOCK = 1 to GRAPH.BLOCKS-1
 *          Let SEG = GRAPH.CONTENTS[BLOCK]
 *          Let NSEG = N>>GRAPH.LEVELS[BLOCK]
 *          dctiv( SEG, NSEG )
 *          uipc( SEG, SEG, RISE )
 *       uipc( DATA+N, DATA, RISE )
 *
 *
 *  Inputs:
 *	(hedge *)graph		This is the input and output array.
 *
 *	(int)n			This is the period of the original signal
 *
 *	(interval *)rise	This is a sampled rising cutoff function.
 *
 * Outputs:
 *	The array binary tree at `parent' is filled with a complete
 *	dyadic local cosine analysis of the signal in the first row.
 *
 * External functions called:
 *	uipc()		Declared in "fold.h"
 *	dct_iv()	Declared in "dtts.h"
 *
 * Assumptions:
 *	1. n>0 is a power of 2.
 *	2. rise != 0
 */
extern void
  lcsdf(
	hedge *graph,		/* Pointers to subintervals in an array. */
	int n,			/* Length of the original signal.        */
	interval *rise)		/* Sampled rising cutoff function.       */
{
  int log2n, log2nseg, block;
  real  *seg;

  assert( n>0 );
  assert( rise );

  /* Compute the logarithm base 2 of the signal length `n': */
  log2n = 0;
  while( (1<<log2n) < n ) ++log2n;
  assert((1<<log2n)==n);

  dct_iv( graph->contents[0], log2n - graph->levels[0] );
  for( block = 1; block<graph->blocks; block++ )
    {
      seg = graph->contents[block];
      log2nseg = log2n - graph->levels[block];
      dct_iv( seg, log2nseg );
      uipc( seg, seg, rise );
    }
  uipc( graph->contents[0], graph->contents[0]+n, rise );

  return;
}
