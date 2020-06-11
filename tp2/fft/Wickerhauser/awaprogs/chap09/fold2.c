/* 
 * Separable 2-D folding/unfolding routines.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <assert.h>
#include "real.h"
#include "rcf.h"
#include "interval.h"
#include "fold.h"
#include "fold2.h"

/**********************************************************************
 * fdc2()
 *
 * [F]old [D]isjoint [C]osine, [2]-dimensional, separable
 *
 *  Scheme:
 *
 *    `in[]'         `work[]'                `out?[]'
 *    ------         --------                --------
 *  0  1  2  3      0  4  8  c  ____     0  1'        
 *  4  5  6  7 -->  1' 5' 9' d'     \    4' 5"        2' 3   
 *  8  9  a  b			     ->               6" 7' 
 *  c  d  e  f                 	         8' 9"    ->       
 *            \                          c  d'   /    a" b'  
 *	       \		                /     e' f  
 *              ->  2' 6' a' e'    ____________/            
 *                  3  7  b  f 	               
 *
 *
 *  Calling sequence and basic algorithm:
 *
 *    fdc2( OUT0,OUT1,OUT2,OUT3, IN, X0, X1, Y0, Y1, WORK, RX, RY ):
 *       Let IX = X0 + X1
 *       Let IY = Y0 + Y1
 *       Let IPTR = IN + Y0
 *       Let WPTR = WORK + Y0*IX
 *       For I = 0 to IX-1
 *         fdcn( WPTR+I, IX, IPTR, IPTR, Y0, RY )
 *         IPTR += IY
 *       Let WPTR = WORK + X0
 *       OUT0 += Y0*X0
 *       For I = 0 to Y0-1
 *         fdcn( OUT0+I, Y0, WPTR, WPTR, X0, RX )
 *         fdcp( OUT2+I, Y0, WPTR, WPTR, X1, RX )
 *         WPTR += IX
 *       Let IPTR = IN + Y0
 *       Let WPTR = WORK
 *       For I = 0 to IX-1
 *         fdcp( WPTR+I, IX, IPTR, IPTR, Y1, RY )
 *         IPTR += IY
 *       Let WPTR = WORK + X0
 *       OUT1 += X0*Y1
 *       For I=0 to Y1-1
 *         fdcn( OUT1+I, Y1, WPTR, WPTR, X0, RX )
 *         fdcp( OUT3+I, Y1, WPTR, WPTR, X1, RX )
 *         WPTR += IX
 *
 *  Inputs:
 *	(real *)out0		These arrays must be preallocated and 
 *	(real *)out1		  assigned. They will be overwritten,
 *	(real *)out2		  so they must not overlap.
 *	(real *)out3
 *
 *	(const real *)in	This array must be preallocated and defined
 *				  in at least `(x0+x1)*(y0+y1)' locations.
 *				  It is not changed by this function.
 *
 *	(int)x0			These positive integers determine the number
 *	(int)x1			  of rows and columns in the arrays `in[]'
 *	(int)y0			  (respectively `r0+r1' and `c0+c1') and the
 *	(int)y1			  number of rows and columns of the 4 output
 *				  arrays `out0[],out1[],out2[],out3[]',
 *				  respectively `x0*y0, x0*y1, x1*y0, x1*y1'.
 *
 *	(real *)work		This temporary array must be preallocated
 *				  with at least `(r0+r1)*max(c0,c1)' memory
 *				  locations.  It will be overwritten.
 *
 *      (const interval *)rx,	These must be preallocated and assigned 
 *	(const interval *)ry,	  arrays defining rising cutoff functions
 *				  to be used by the folding functions.
 *
 *  Outputs:
 *	(real *)out0		These arrays are filled by side effect.
 *	(real *)out1
 *	(real *)out2
 *	(real *)out3
 *
 *	(real *)work		This array is trashed by side effect.
 *
 *  External functions called:
 *	fdcp(), fdcn(), assert()
 *
 *  Assumptions:
 *	1. `x0,x1,y0,y1' are positive integers.
 *	2. rx != NULL;  ry != NULL
 *	3. 0 < 1 + rx->final - rx->least <= min(x0,x1).
 *	4. 0 < 1 + ry->final - ry->least <= min(y0,y1).
 */
extern void
  fdc2(
       real *out0,		/* (X=0, Y=0) output array pointer. */
       real *out1,		/* (X=0, Y=1) output array pointer. */
       real *out2,		/* (X=1, Y=0) output array pointer. */
       real *out3,		/* (X=1, Y=1) output array pointer. */
       const real *in,		/* Pointer to input (parent) array. */
       int x0,			/* Rows in arrays `out0[],out1[]' */
       int x1,			/* Rows in arrays `out2[],out3[]' */
       int y0,			/* Columns in arrays `out0[],out2[]' */
       int y1,			/* Columns in arrays `out1[],out3[]' */
       real *work,		/* Preallocated scratch array. */
       const interval *rx,	/* Rising cutoff function for X-folding. */
       const interval *ry)	/* Rising cutoff function for Y-folding. */
{
  int i, ix, iy;
  const real *iptr;
  real *wptr;

  assert(x0>0);
  assert(x1>0);
  assert(y0>0);
  assert(y1>0);
  assert(rx);
  assert(ry);
  assert( 0 < 1 + rx->final - rx->least );
  assert( 1 + rx->final - rx->least <= min(x0,x1) );
  assert( 0 < 1 + ry->final - ry->least );
  assert( 1 + ry->final - ry->least <= min(y0,y1) );

  ix = x0+x1;  iy = y0+y1;  /* Number of rows and columns of `in[]'. */

  /* Y-DIRECTION, left half:     -then-    X-DIRECTION:
   *
   *       `in + y0'     `work + x0'
   *          /               /
   * in[]:   /    work[]:    /            out0[]:
   *  0  1 *2  3       0  4 *8  c              0  1'
   *  4  5  6  7  -->  1' 5' 9' d'      -->    4' 5"
   *  8  9  a  b       *                       *--- `out0 + x0*y0'
   *  c  d  e  f        \                 out2[]:
   *               `work + ix*y0'              8' 9"
   *                                           c  d'
   */
  iptr = in + y0;		/* Move to midpoint of first row of `in[]'. */
  wptr = work + ix*y0;		/* Move to end of first column of `work[]'. */
  for( i=0; i<ix; ++i )
    {
      fdcn( wptr++, ix, iptr, iptr, y0, ry );
      iptr += iy;
    }
  wptr = work + x0;		/* Move to midpoint of 1st row of `work[]'. */
  out0 += x0*y0;		/* Move to end of first column of `out0[]'. */
  for( i=0; i<y0; ++i )
    {
      fdcn( out0++, y0, wptr, wptr, x0, rx );
      fdcp( out2++, y0, wptr, wptr, x1, rx );
      wptr += ix;
    }

  /* Y-DIRECTION, right half:       -then-    X-DIRECTION:
   *
   *      `in + y0'
   *          /
   * in[]:   /      work[]:                 out1[]:  
   *  0  1 *2  3           `work + x0'          2' 3
   *  4  5  6  7             /                  6" 7'
   *  8  9  a  b  -->  2' 6'*a' e'   -->        *----`out1 + x0*y1'
   *  c  d  e  f       3  7  b  f          out3[]:         
   *                                            a" b'       
   *                        		        e' f
   */
  iptr = in + y0;		/* Move to midpoint of first row of `in[]'. */
  wptr = work;			/* Start at top of first column of `work[]'. */
  for( i=0; i<ix; ++i )
    {
      fdcp( wptr++, ix, iptr, iptr, y1, ry );
      iptr += iy;
    }
  wptr = work + x0;		/* Move to midpoint of 1st row of `work[]'. */
  out1 += x0*y1;		/* Move to end of first column of `out1[]'. */
  for( i=0; i<y1; ++i )
    {
      fdcn( out1++, y1, wptr, wptr, x0, rx );
      fdcp( out3++, y1, wptr, wptr, x1, rx );
      wptr += ix;
    }
  return;
}


/**********************************************************************
 * fds2()
 *
 * [F]old [D]isjoint [S]ine, [2]-dimensional, separable
 *
 *  Scheme:
 *
 *    `in[]'         `work[]'                `out?[]'
 *    ------         --------                --------
 *  0  1  2  3      0  4  8  c  ____     0  1'        
 *  4  5  6  7 -->  1' 5' 9' d'     \    4' 5"        2' 3   
 *  8  9  a  b			     ->               6" 7' 
 *  c  d  e  f                 	         8' 9"    ->       
 *            \                          c  d'   /    a" b'  
 *	       \		                /     e' f  
 *              ->  2' 6' a' e'    ____________/            
 *                  3  7  b  f 	               
 *
 *
 *  Calling sequence and basic algorithm:
 *
 *    fds2( OUT0,OUT1,OUT2,OUT3, IN, X0, X1, Y0, Y1, WORK, RX, RY ):
 *       Let IX = X0 + X1
 *       Let IY = Y0 + Y1
 *       Let IPTR = IN + Y0
 *       Let WPTR = WORK + Y0*IX
 *       For I = 0 to IX-1
 *         fdsn( WPTR+I, IX, IPTR, IPTR, Y0, RY )
 *         IPTR += IY
 *       Let WPTR = WORK + X0
 *       OUT0 += Y0*X0
 *       For I = 0 to Y0-1
 *         fdsn( OUT0+I, Y0, WPTR, WPTR, X0, RX )
 *         fdsp( OUT2+I, Y0, WPTR, WPTR, X1, RX )
 *         WPTR += IX
 *       Let IPTR = IN + Y0
 *       Let WPTR = WORK
 *       For I = 0 to IX-1
 *         fdsp( WPTR+I, IX, IPTR, IPTR, Y1, RY )
 *         IPTR += IY
 *       Let WPTR = WORK + X0
 *       OUT1 += X0*Y1
 *       For I=0 to Y1-1
 *         fdsn( OUT1+I, Y1, WPTR, WPTR, X0, RX )
 *         fdsp( OUT3+I, Y1, WPTR, WPTR, X1, RX )
 *         WPTR += IX
 *
 *  Inputs:
 *	(real *)out0		These arrays must be preallocated and 
 *	(real *)out1		  assigned. They will be overwritten,
 *	(real *)out2		  so they must not overlap.
 *	(real *)out3
 *
 *	(const real *)in	This array must be preallocated and defined
 *				  in at least `(x0+x1)*(y0+y1)' locations.
 *				  It is not changed by this function.
 *
 *	(int)x0			These positive integers determine the number
 *	(int)x1			  of rows and columns in the arrays `in[]'
 *	(int)y0			  (respectively `r0+r1' and `c0+c1') and the
 *	(int)y1			  number of rows and columns of the 4 output
 *				  arrays `out0[],out1[],out2[],out3[]',
 *				  respectively `x0*y0, x0*y1, x1*y0, x1*y1'.
 *
 *	(real *)work		This temporary array must be preallocated
 *				  with at least `(r0+r1)*max(c0,c1)' memory
 *				  locations.  It will be overwritten.
 *
 *      (const interval *)rx,	These must be preallocated and assigned 
 *	(const interval *)ry,	  arrays defining rising cutoff functions
 *				  to be used by the folding functions.
 *
 *  Outputs:
 *	(real *)out0		These arrays are filled by side effect.
 *	(real *)out1
 *	(real *)out2
 *	(real *)out3
 *
 *	(real *)work		This array is trashed by side effect.
 *
 *  External functions called:
 *	fdsp(), fdsn(), assert()
 *
 *  Assumptions:
 *	1. `x0,x1,y0,y1' are positive integers.
 *	2. rx != NULL;  ry != NULL
 *	3. 0 < 1 + rx->final - rx->least <= min(x0,x1).
 *	4. 0 < 1 + ry->final - ry->least <= min(y0,y1).
 */
extern void
  fds2(
       real *out0,		/* (X=0, Y=0) output array pointer. */
       real *out1,		/* (X=0, Y=1) output array pointer. */
       real *out2,		/* (X=1, Y=0) output array pointer. */
       real *out3,		/* (X=1, Y=1) output array pointer. */
       const real *in,		/* Pointer to input (parent) array. */
       int x0,			/* Rows in arrays `out0[],out1[]' */
       int x1,			/* Rows in arrays `out2[],out3[]' */
       int y0,			/* Columns in arrays `out0[],out2[]' */
       int y1,			/* Columns in arrays `out1[],out3[]' */
       real *work,		/* Preallocated scratch array. */
       const interval *rx,	/* Rising cutoff function for X-folding. */
       const interval *ry)	/* Rising cutoff function for Y-folding. */
{
  int i, ix, iy;
  const real *iptr;
  real *wptr;

  assert(x0>0);
  assert(x1>0);
  assert(y0>0);
  assert(y1>0);
  assert(rx);
  assert(ry);
  assert( 0 < 1 + rx->final - rx->least );
  assert( 1 + rx->final - rx->least <= min(x0,x1) );
  assert( 0 < 1 + ry->final - ry->least );
  assert( 1 + ry->final - ry->least <= min(y0,y1) );

  ix = x0+x1;  iy = y0+y1;  /* Number of rows and columns of `in[]'. */

  /* Y-DIRECTION, left half:     -then-    X-DIRECTION:
   *
   *       `in + y0'     `work + x0'
   *          /               /
   * in[]:   /    work[]:    /            out0[]:
   *  0  1 *2  3       0  4 *8  c              0  1'
   *  4  5  6  7  -->  1' 5' 9' d'      -->    4' 5"
   *  8  9  a  b       *                       *--- `out0 + x0*y0'
   *  c  d  e  f        \                 out2[]:
   *               `work + ix*y0'              8' 9"
   *                                           c  d'
   */
  iptr = in + y0;		/* Move to midpoint of first row of `in[]'. */
  wptr = work + ix*y0;		/* Move to end of first column of `work[]'. */
  for( i=0; i<ix; ++i )
    {
      fdsn( wptr++, ix, iptr, iptr, y0, ry );
      iptr += iy;
    }
  wptr = work + x0;		/* Move to midpoint of 1st row of `work[]'. */
  out0 += x0*y0;		/* Move to end of first column of `out0[]'. */
  for( i=0; i<y0; ++i )
    {
      fdsn( out0++, y0, wptr, wptr, x0, rx );
      fdsp( out2++, y0, wptr, wptr, x1, rx );
      wptr += ix;
    }

  /* Y-DIRECTION, right half:       -then-    X-DIRECTION:
   *
   *      `in + y0'
   *          /
   * in[]:   /      work[]:                 out1[]:  
   *  0  1 *2  3           `work + x0'          2' 3
   *  4  5  6  7             /                  6" 7'
   *  8  9  a  b  -->  2' 6'*a' e'   -->        *----`out1 + x0*y1'
   *  c  d  e  f       3  7  b  f          out3[]:         
   *                                            a" b'       
   *                        		        e' f
   */
  iptr = in + y0;		/* Move to midpoint of first row of `in[]'. */
  wptr = work;			/* Start at top of first column of `work[]'. */
  for( i=0; i<ix; ++i )
    {
      fdsp( wptr++, ix, iptr, iptr, y1, ry );
      iptr += iy;
    }
  wptr = work + x0;		/* Move to midpoint of 1st row of `work[]'. */
  out1 += x0*y1;		/* Move to end of first column of `out1[]'. */
  for( i=0; i<y1; ++i )
    {
      fdsn( out1++, y1, wptr, wptr, x0, rx );
      fdsp( out3++, y1, wptr, wptr, x1, rx );
      wptr += ix;
    }
  return;
}

/**********************************************************************
 * udc2()
 *
 *  [U]nfolding, [D]isjoint [C]osine polarity; separable [2]-dimensional.
 *
 *
 *  Scheme:
 *
 *      `in?[]'          `work[]'          `out[]'              
 *      ---------          --------         -------
 *           __________  0  4  8  c    
 *  0  1'   /            1' 5' 9' d'\     0  1  2  3
 *  4' 5"\ /                         \__  4  5  6  7
 *        /  2' 3                    /    8  9  a  b
 *  8' 9"/   6" 7'\      2' 6' a' e'/     c  d  e  f
 *  c  d'          >--   3  7  b  f    
 *           a" b'/   
 *           e' f   
 *                                    
 *
 *  Calling sequence and basic algorithm:
 *
 *    udc2( OUT, IN0, IN1, IN2, IN3, X0, X1, Y0, Y1, WORK, RX, RY ):
 *       Let OX = X0 + X1
 *       Let OY = Y0 + Y1
 *       Let WPTR = WORK + Y0*OX
 *       Let INNEG = IN0 + Y0
 *       Let INPOS = IN1
 *       For I = 0 to X0-1
 *         udcn( WPTR+I, OX, INNEG, INPOS, Y0, RY )
 *         INNEG += Y0
 *         INPOS += Y1
 *       Let WPTR = WORK + Y0*OX + X0
 *       Let INNEG = IN2 + Y0
 *       Let INPOS = IN3
 *       For I = 0 to X0-1
 *         udcn( WPTR+I, OX, INNEG, INPOS, Y0, RY )
 *         INNEG += Y0
 *         INPOS += Y1
 *       Let WPTR = WORK + X0
 *       OUT += OY*X0
 *       For I = 0 to Y0-1
 *         udcn( OUT+I, OY, WPTR, WPTR, X0, RX )
 *         udcp( OUT+I, OY, WPTR, WPTR, X1, RX )
 *         WPTR += OX
 *       Let WPTR = WORK
 *       Let INNEG = IN0 + Y0
 *       Let INPOS = IN1
 *       For I = 0 to X0-1
 *         udcp( WPTR+I, OX, INNEG, INPOS, Y1, RY )
 *         INNEG += Y0
 *         INPOS += Y1
 *       Let WPTR = WORK + X0
 *       Let INNEG = IN2 + Y0
 *       Let INPOS = IN3
 *       For I = 0 to X1-1
 *         udcp( WPTR+I, OX, INNEG, INPOS, Y1, RY )
 *         INNEG += Y0
 *         INPOS += Y1
 *       Let WPTR = WORK + X0
 *       OUT += Y0
 *       For I = 0 to Y1-1
 *         udcn( OUT+I, OY, WPTR, WPTR, X0, RX )
 *         udcp( OUT+I, OY, WPTR, WPTR, X1, RX )
 *         WPTR += OX
 *
 *
 *  Inputs:
 *	(real *)out		This array must be preallocated with at 
 *				  least `(y0+y1)*(x0+x1)' elements.
 *				  It will be overwritten.
 *
 *	(const real *)in0	These arrays must be preallocated and defined
 *	(const real *)in1	  with at least `x0*y0, x0*y1, x1*y0, x1*y1'
 *	(const real *)in2	  locations, respectively, each.  They are
 *	(const real *)in3	  not changed by this function.
 *
 *	(int)x0			Array `in0[]' has `x0' rows and `y0' columns,
 *	(int)x1			  `in1[]' has `x0' rows and `y1' columns,
 *	(int)y0			  `in2[]' has `x1' rows and `y0' columns,
 *	(int)y1			  `in3[]' has `x1' rows and `y1' columns,
 *
 *	(real *)work		This temporary array must be preallocated
 *				  with at least `(x0+x1)*max(y0,y1)' memory
 *				  locations.  It will be overwritten.
 *
 *      (const interval *)rx,	These must be preallocated and assigned 
 *	(const interval *)ry,	  arrays defining rising cutoff functions
 *				  to be used by the unfolding functions.
 *
 *  Outputs:
 *	(real *)out		This array is filled by side effect.
 *
 *	(real *)work		This array is trashed by side effect.
 *
 *  External functions called:
 *	udcn(), udcp(), assert()
 *
 *  Assumptions:
 *	1. `y0,y1,x0,x1' are positive integers.
 *	2. rx != NULL;  ry != NULL
 *	3. 0 < 1 + rx->final - rx->least <= min(x0,x1).
 *	4. 0 < 1 + ry->final - ry->least <= min(y0,y1).
 */
extern void
  udc2(				/* Separable 2-D adjoint folding. */
       real *out,		/* Pointer to output array. */
       const real *in0,		/* Pointer to (x=0, y=0) input array. */
       const real *in1,		/* Pointer to (x=0, y=1) input array. */
       const real *in2,		/* Pointer to (x=1, y=0) input array. */
       const real *in3,		/* Pointer to (x=1, y=1) input array. */
       int x0,			/* Rows in arrays `in0[],in1[]' */
       int x1,			/* Rows in arrays `in2[],in3[]' */
       int y0,			/* Columns in arrays `in0[],in2[]' */
       int y1,			/* Columns in arrays `in1[],in3[]' */
       real *work,		/* Preallocated scratch array. */
       const interval *rx,	/* Rising cutoff function for X-unfolding. */
       const interval *ry)	/* Rising cutoff function for Y-unfolding. */
{
  int i, oy, ox;
  real *wptr;
  const real *inpos, *inneg;

  assert(y0>0);
  assert(y1>0);
  assert(x0>0);
  assert(x1>0);
  assert(rx);
  assert(ry);
  assert( 0 < 1 + rx->final - rx->least );
  assert( 1 + rx->final - rx->least <= min(x0,x1) );
  assert( 0 < 1 + ry->final - ry->least );
  assert( 1 + ry->final - ry->least <= min(y0,y1) );

  oy = y0+y1;  ox = x0+x1;	/* Number of columns and rows in `out[]' */

  /* Y-DIRECTION, left half:     -then-     X-DIRECTION:
   *
   *  in0[]:    in1[]:         work[]:
   *   0  1' *   *2' 3       0  4'   *8' c   out[]:
   *   4' 5"      6" 7'      1  5'    9' d     0  1  .  .
   *     \___________\_____/ *      / *   \__  4  5  .  .
   *                               /          *8  9  .  .
   *   in2[]:    in3[]:           /            c  d  .  .
   *    8' 9" *   *a" b'         /
   *    c  d'      e' f         /
   *        \__________\_______/
   */
  wptr = work + y0*ox;		/* Move past end of first `work[]' column. */
  inneg = in0 + y0;		/* Move past end of first row of `in0[]'.  */
  inpos = in1;			/* Move to start of first row of `in1[]'.  */
  for( i=0; i<x0; ++i )
    {
      udcn( wptr++, ox, inneg, inpos, y0, ry );
      inneg += y0;  inpos += y1;
    }
  inneg = in2 + y0;		/* Move past end of first row of `in2[]'. */
  inpos = in3;			/* Move to start of first row of `in3[]'. */
  for( i=0; i<x1; ++i )
    {
      udcn( wptr++, ox, inneg, inpos, y0, ry );
      inneg += y0;  inpos += y1;
    }
  wptr = work + x0;		/* Move to middle of first `work[]' row.   */
  out += oy*x0;			/* Move to middle of first `out[]' column. */
  for( i=0; i<y0; ++i )
    {
      udcn(  out,  oy, wptr, wptr, x0, rx );
      udcp( out++, oy, wptr, wptr, x1, rx );
      wptr += ox;
    }

  /* Y-DIRECTION, right half:     -then-     X-DIRECTION:
   *
   *  in0[]:    in1[]:         work[]:
   *   0  1' *   *2' 3                       out[]:
   *   4' 5"      6" 7'                         .  .  2  3
   *     \___________\______                __  .  .  6  7
   *                        \              /    .  . *a  b
   *   in2[]:    in3[]:     *2  6'   *a' e      .  .  e  f
   *    8' 9" *   *a" b'     3  7'    b' f
   *    c  d'      e' f             /
   *        \__________\___________/
   */
  wptr = work;			/* Move to start of first `work[]' column. */
  inneg = in0 + y0;		/* Move past end of first row of `in0[]'.  */
  inpos = in1;			/* Move to start of first row of `in1[]'.  */
  for( i=0; i<x0; ++i )
    {
      udcp( wptr++, ox, inneg, inpos, y1, ry );
      inneg += y0;  inpos += y1;
    }
  inneg = in2 + y0;		/* Move past end of first row of `in2[]'. */
  inpos = in3;			/* Move to start of first row of `in3[]'. */
  for( i=0; i<x1; ++i )
    {
      udcp( wptr++, ox, inneg, inpos, y1, ry );
      inneg += y0;  inpos += y1;
    }
  wptr = work + x0;		/* Move to middle of first `work[]' row. */
  /* `out[]' currently points to the middle of the middle `out[]' column. */
  for( i=0; i<y1; ++i )
    {
      udcn(  out,  oy, wptr, wptr, x0, rx );
      udcp( out++, oy, wptr, wptr, x1, rx );
      wptr += ox;
    }
  return;
}

/**********************************************************************
 * uds2()
 *
 *  [U]nfolding, [D]isjoint [S]ine polarity; separable [2]-dimensional.
 *
 *
 *  Scheme:
 *
 *      `in?[]'          `work[]'          `out[]'              
 *      ---------          --------         -------
 *           __________  0  4  8  c    
 *  0  1'   /            1' 5' 9' d'\     0  1  2  3
 *  4' 5"\ /                         \__  4  5  6  7
 *        /  2' 3                    /    8  9  a  b
 *  8' 9"/   6" 7'\      2' 6' a' e'/     c  d  e  f
 *  c  d'          >--   3  7  b  f    
 *           a" b'/   
 *           e' f   
 *                                    
 *
 *  Calling sequence and basic algorithm:
 *
 *    uds2( OUT, IN0, IN1, IN2, IN3, X0, X1, Y0, Y1, WORK, RX, RY ):
 *       Let OX = X0 + X1
 *       Let OY = Y0 + Y1
 *       Let WPTR = WORK + Y0*OX
 *       Let INNEG = IN0 + Y0
 *       Let INPOS = IN1
 *       For I = 0 to X0-1
 *         udsn( WPTR+I, OX, INNEG, INPOS, Y0, RY )
 *         INNEG += Y0
 *         INPOS += Y1
 *       Let WPTR = WORK + Y0*OX + X0
 *       Let INNEG = IN2 + Y0
 *       Let INPOS = IN3
 *       For I = 0 to X0-1
 *         udsn( WPTR+I, OX, INNEG, INPOS, Y0, RY )
 *         INNEG += Y0
 *         INPOS += Y1
 *       Let WPTR = WORK + X0
 *       OUT += OY*X0
 *       For I = 0 to Y0-1
 *         udsn( OUT+I, OY, WPTR, WPTR, X0, RX )
 *         udsp( OUT+I, OY, WPTR, WPTR, X1, RX )
 *         WPTR += OX
 *       Let WPTR = WORK
 *       Let INNEG = IN0 + Y0
 *       Let INPOS = IN1
 *       For I = 0 to X0-1
 *         udsp( WPTR+I, OX, INNEG, INPOS, Y1, RY )
 *         INNEG += Y0
 *         INPOS += Y1
 *       Let WPTR = WORK + X0
 *       Let INNEG = IN2 + Y0
 *       Let INPOS = IN3
 *       For I = 0 to X1-1
 *         udsp( WPTR+I, OX, INNEG, INPOS, Y1, RY )
 *         INNEG += Y0
 *         INPOS += Y1
 *       Let WPTR = WORK + X0
 *       OUT += Y0
 *       For I = 0 to Y1-1
 *         udsn( OUT+I, OY, WPTR, WPTR, X0, RX )
 *         udsp( OUT+I, OY, WPTR, WPTR, X1, RX )
 *         WPTR += OX
 *
 *
 *  Inputs:
 *	(real *)out		This array must be preallocated with at 
 *				  least `(y0+y1)*(x0+x1)' elements.
 *				  It will be overwritten.
 *
 *	(const real *)in0	These arrays must be preallocated and defined
 *	(const real *)in1	  with at least `x0*y0, x0*y1, x1*y0, x1*y1'
 *	(const real *)in2	  locations, respectively, each.  They are
 *	(const real *)in3	  not changed by this function.
 *
 *	(int)x0			Array `in0[]' has `x0' rows and `y0' columns,
 *	(int)x1			  `in1[]' has `x0' rows and `y1' columns,
 *	(int)y0			  `in2[]' has `x1' rows and `y0' columns,
 *	(int)y1			  `in3[]' has `x1' rows and `y1' columns,
 *
 *	(real *)work		This temporary array must be preallocated
 *				  with at least `(x0+x1)*max(y0,y1)' memory
 *				  locations.  It will be overwritten.
 *
 *      (const interval *)rx,	These must be preallocated and assigned 
 *	(const interval *)ry,	  arrays defining rising cutoff functions
 *				  to be used by the unfolding functions.
 *
 *  Outputs:
 *	(real *)out		This array is filled by side effect.
 *
 *	(real *)work		This array is trashed by side effect.
 *
 *  External functions called:
 *	udsn(), udsp(), assert()
 *
 *  Assumptions:
 *	1. `y0,y1,x0,x1' are positive integers.
 *	2. rx != NULL;  ry != NULL
 *	3. 0 < 1 + rx->final - rx->least <= min(x0,x1).
 *	4. 0 < 1 + ry->final - ry->least <= min(y0,y1).
 */
extern void
  uds2(				/* Separable 2-D adjoint folding. */
       real *out,		/* Pointer to output array. */
       const real *in0,		/* Pointer to (x=0, y=0) input array. */
       const real *in1,		/* Pointer to (x=0, y=1) input array. */
       const real *in2,		/* Pointer to (x=1, y=0) input array. */
       const real *in3,		/* Pointer to (x=1, y=1) input array. */
       int x0,			/* Rows in arrays `in0[],in1[]' */
       int x1,			/* Rows in arrays `in2[],in3[]' */
       int y0,			/* Columns in arrays `in0[],in2[]' */
       int y1,			/* Columns in arrays `in1[],in3[]' */
       real *work,		/* Preallocated scratch array. */
       const interval *rx,	/* Rising cutoff function for X-unfolding. */
       const interval *ry)	/* Rising cutoff function for Y-unfolding. */
{
  int i, oy, ox;
  real *wptr;
  const real *inpos, *inneg;

  assert(y0>0);
  assert(y1>0);
  assert(x0>0);
  assert(x1>0);
  assert(rx);
  assert(ry);
  assert( 0 < 1 + rx->final - rx->least );
  assert( 1 + rx->final - rx->least <= min(x0,x1) );
  assert( 0 < 1 + ry->final - ry->least );
  assert( 1 + ry->final - ry->least <= min(y0,y1) );

  oy = y0+y1;  ox = x0+x1;	/* Number of columns and rows in `out[]' */

  /* Y-DIRECTION, left half:     -then-     X-DIRECTION:
   *
   *  in0[]:    in1[]:         work[]:
   *   0  1' *   *2' 3       0  4'   *8' c   out[]:
   *   4' 5"      6" 7'      1  5'    9' d     0  1  .  .
   *     \___________\_____/ *      / *   \__  4  5  .  .
   *                               /          *8  9  .  .
   *   in2[]:    in3[]:           /            c  d  .  .
   *    8' 9" *   *a" b'         /
   *    c  d'      e' f         /
   *        \__________\_______/
   */
  wptr = work + y0*ox;		/* Move past end of first `work[]' column. */
  inneg = in0 + y0;		/* Move past end of first row of `in0[]'.  */
  inpos = in1;			/* Move to start of first row of `in1[]'.  */
  for( i=0; i<x0; ++i )
    {
      udsn( wptr++, ox, inneg, inpos, y0, ry );
      inneg += y0;  inpos += y1;
    }
  inneg = in2 + y0;		/* Move past end of first row of `in2[]'. */
  inpos = in3;			/* Move to start of first row of `in3[]'. */
  for( i=0; i<x1; ++i )
    {
      udsn( wptr++, ox, inneg, inpos, y0, ry );
      inneg += y0;  inpos += y1;
    }
  wptr = work + x0;		/* Move to middle of first `work[]' row.   */
  out += oy*x0;			/* Move to middle of first `out[]' column. */
  for( i=0; i<y0; ++i )
    {
      udsn(  out,  oy, wptr, wptr, x0, rx );
      udsp( out++, oy, wptr, wptr, x1, rx );
      wptr += ox;
    }

  /* Y-DIRECTION, right half:     -then-     X-DIRECTION:
   *
   *  in0[]:    in1[]:         work[]:
   *   0  1' *   *2' 3                       out[]:
   *   4' 5"      6" 7'                         .  .  2  3
   *     \___________\______                __  .  .  6  7
   *                        \              /    .  . *a  b
   *   in2[]:    in3[]:     *2  6'   *a' e      .  .  e  f
   *    8' 9" *   *a" b'     3  7'    b' f
   *    c  d'      e' f             /
   *        \__________\___________/
   */
  wptr = work;			/* Move to start of first `work[]' column. */
  inneg = in0 + y0;		/* Move past end of first row of `in0[]'.  */
  inpos = in1;			/* Move to start of first row of `in1[]'.  */
  for( i=0; i<x0; ++i )
    {
      udsp( wptr++, ox, inneg, inpos, y1, ry );
      inneg += y0;  inpos += y1;
    }
  inneg = in2 + y0;		/* Move past end of first row of `in2[]'. */
  inpos = in3;			/* Move to start of first row of `in3[]'. */
  for( i=0; i<x1; ++i )
    {
      udsp( wptr++, ox, inneg, inpos, y1, ry );
      inneg += y0;  inpos += y1;
    }
  wptr = work + x0;		/* Move to middle of first `work[]' row. */
  /* `out[]' currently points to the middle of the middle `out[]' column. */
  for( i=0; i<y1; ++i )
    {
      udsn(  out,  oy, wptr, wptr, x0, rx );
      udsp( out++, oy, wptr, wptr, x1, rx );
      wptr += ox;
    }
  return;
}
