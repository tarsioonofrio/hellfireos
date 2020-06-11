/*
 * 
 * Declare separable, 2-dimensional folding and unfolding functions.
 * 
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef FOLD2_HDR_ALREADY_INCLUDED
# define FOLD2_HDR_ALREADY_INCLUDED

#include "real.h"
#include "interval.h"

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
       const interval *ry);	/* Rising cutoff function for Y-folding. */

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
       const interval *ry);	/* Rising cutoff function for Y-folding. */

extern void
  udc2(
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
       const interval *ry);	/* Rising cutoff function for Y-unfolding. */

extern void
  uds2(
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
       const interval *ry);	/* Rising cutoff function for Y-unfolding. */

#endif /* FOLD_HDR_ALREADY_INCLUDED */
