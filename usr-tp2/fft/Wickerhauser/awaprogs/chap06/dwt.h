/*
 * Declare functions for one-dimensional discrete wavelet transforms.
 * 
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 * 
 */

#ifndef DWT_HDR_ALREADY_INCLUDED
# define DWT_HDR_ALREADY_INCLUDED

#include "real.h"
#include "interval.h"
#include "qf.h"

/* We concatenate from the end by default, unless overridden:  */
#ifndef DWTA_CONCATENATE_FROM_START
 /* Shifts to concatenate intervals from the end of the long array */
 #define shifttoorigin( I )         ( -1 - (I)->final )
 #define shifttonextinterval( I )   ( (I)->least - (I)->final - 1 )
#else
 /* Shifts to concatenate intervals from the start of the long array */
 #define shifttoorigin( I )         ( -(I)->least )
 #define shifttonextinterval( I )   ( 1 + (I)->final -  (I)->least )
#endif

extern void
  dwtpd0( 
	 real *difs,     /* Output array with N elements. */
	 real *sums,     /* Working array with N zeroes.  */
	 const real *in, /* Input array with N elements.  */
	 int N,          /* This is the number of inputs. */
	 const pqf *H,   /* Low-pass QF data structure.   */
	 const pqf *G);  /* High-pass QF data structure.  */

extern void
  dwtpd0n( 
	  real *difs,     /* Output array with N elements. */
	  real *sums,     /* Working array with N zeroes.  */
	  const real *in, /* Input array with N elements.  */
	  int N,          /* This is the number of inputs. */
	  const pqf *H,   /* Low-pass QF data structure.   */
	  const pqf *G);  /* High-pass QF data structure.  */

extern void
  dwtpi0( 
	 real *data,     /* Input/Output array, length N. */
	 real *work,     /* Working array, length N.      */
	 int N,          /* This is the number of inputs. */
	 const pqf *H,   /* Low-pass QF data structure.   */
	 const pqf *G);  /* High-pass QF data structure.  */

extern void
  dwtpi0n( 
	  real *data,     /* Input/Output array, length N. */
	  real *work,     /* Working array, length N.      */
	  int N,          /* This is the number of inputs. */
	  const pqf *H,   /* Low-pass QF data structure.   */
	  const pqf *G);  /* High-pass QF data structure.  */

extern void
  dwtpd( 
	real *difs,       /* Output array with N elements. */
	real *sums,       /* Working array with N zeroes.  */
	const real *in,   /* Input array with N elements.  */
	int N,            /* This is the number of inputs. */
	int L,            /* This is the number of levels. */
	const pqf *H,     /* Low-pass QF data structure.   */
	const pqf *G);    /* High-pass QF data structure.  */

extern void
  dwtpdn( 
	 real *difs,      /* Output array with N elements. */
	 real *sums,      /* Working array with N zeroes.  */
	 const real *in,  /* Input array with N elements.  */
	 int N,           /* This is the number of inputs. */
	 int L,           /* This is the number of levels. */
	 const pqf *H,    /* Low-pass QF data structure.   */
	 const pqf *G);   /* High-pass QF data structure.  */

extern void
  dwtpi( 
	real *data,     /* Input/Output array, length N.  */
	real *work,     /* Working array, length N.       */
	int N,          /* This is the number of inputs.  */
	int L,          /* This is the number of levels.  */
	const pqf *H,   /* Low-pass QF data structure.    */
	const pqf *G);  /* High-pass QF data structure.   */

extern void
  dwtpin( 
	 real *data,     /* Input/Output array, length N.  */
	 real *work,     /* Working array, length N.       */
	 int N,          /* This is the number of inputs.  */
	 int L,          /* This is the number of levels.  */
	 const pqf *H,   /* Low-pass QF data structure.    */
	 const pqf *G);  /* High-pass QF data structure.   */

extern void
  idwtpd0( 
	  real *out,       /* Preallocated, zeroed output.  */
	  real *sums,      /* Working array of N zeroes.    */
	  const real *in,  /* Input array with N elements.  */
	  int N,           /* This is the number of inputs. */
	  const pqf *H,    /* Low-pass QF data structure.   */
	  const pqf *G);   /* High-pass QF data structure.  */

extern void
  idwtpd0n( 
	   real *out,      /* Preallocated, zeroed output.  */
	   real *sums,     /* Working array of N zeroes.    */
	   const real *in, /* Input array with N elements.  */
	   int N,          /* This is the number of inputs. */
	   const pqf *H,   /* Low-pass QF data structure.   */
	   const pqf *G);  /* High-pass QF data structure.  */

extern void
  idwtpi0( 
	  real *data,     /* Input and output, length N.    */
	  real *work,     /* Working array of N elements.   */
	  int N,          /* This is the number of inputs.  */
	  const pqf *H,   /* Low-pass QF data structure.    */
	  const pqf *G);  /* High-pass QF data structure.   */

extern void
  idwtpi0n( 
	   real *data,     /* Input and output, length N.    */
	   real *work,     /* Working array of N elements.   */
	   int N,          /* This is the number of inputs.  */
	   const pqf *H,   /* Low-pass QF data structure.    */
	   const pqf *G);  /* High-pass QF data structure.   */

extern void
  idwtpd( 
	 real *out,        /* Preallocated, zeroed output.  */
	 real *sums,       /* Working array of N zeroes.    */
	 const real *in,   /* Input array with N elements.  */
	 int N,            /* This is the number of inputs. */
	 int L,            /* This is the number of levels. */
	 const pqf *H,     /* Low-pass QF data structure.   */
	 const pqf *G);    /* High-pass QF data structure.  */

extern void
  idwtpdn( 
	  real *out,      /* Output, N preallocated zeroes. */
	  real *sums,     /* Working array of N zeroes.     */
	  const real *in, /* Input array with N elements.   */
	  int N,          /* This is the number of inputs.  */
	  int L,          /* This is the number of levels.  */
	  const pqf *H,   /* Low-pass QF data structure.    */
	  const pqf *G);  /* High-pass QF data structure.   */

extern void
  idwtpi( 
	 real *data,     /* Input and output, length N.    */
	 real *work,     /* Working array of N elements.   */
	 int N,          /* This is the number of inputs.  */
	 int L,          /* This is the number of levels.  */
	 const pqf *H,   /* Low-pass QF data structure.    */
	 const pqf *G);  /* High-pass QF data structure.   */

extern void
  idwtpin( 
	  real *data,     /* Input and output, length N.    */
	  real *work,     /* Working array of N elements.   */
	  int N,          /* This is the number of inputs.  */
	  int L,          /* This is the number of levels.  */
	  const pqf *H,   /* Low-pass QF data structure.    */
	  const pqf *G);  /* High-pass QF data structure.   */

extern void
  dwtaintervals( 
		interval *V,	/* The `L' scaling subspaces.    */
		interval *W,    /* The `L+1' wavelet subspaces.  */
		interval *in,   /* This is the input interval.   */
		int L,          /* This is the number of levels. */
		const pqf *H,   /* Low-pass QF data structure.   */
		const pqf *G);  /* High-pass QF data structure.  */

extern void
  dwtaorigins( 
	      interval *V,	/* The `L' scaling subspaces.    */
	      interval *W,	/* The `L+1' wavelet subspaces.  */
	      real *sums,	/* This is the input interval.   */
	      real *difs,	/* This is the input interval.   */
	      int L);		/* This is the number of levels. */

extern void
  dwta( 
       interval *V,    /* Array of `L' scaling subspaces.   */
       interval *W,    /* Array of `L+1' wavelet subspaces. */
       interval *in,   /* Preallocated input interval.      */
       int L,          /* This is the number of levels.     */
       const pqf *H,   /* Low-pass QF data structure.       */
       const pqf *G);  /* High-pass QF data structure.      */

extern void
  dwtan( 
	interval *V,		/* Array of `L' scaling subspaces.   */
	interval *W,		/* Array of `L+1' wavelet subspaces. */
	interval *in,		/* Preallocated input interval.      */
	int L,			/* This is the number of levels.     */
	const pqf *H,		/* Low-pass QF data structure.       */
	const pqf *G);		/* High-pass QF data structure.      */

extern interval *
  dwtacomplete( 
	       real *data,	/* Array of `length' input values.     */
	       int length,	/* Positive length of the input array. */
	       int L,		/* This is the number of levels.       */
	       const pqf *H,	/* Low-pass QF data structure.         */
	       const pqf *G);	/* High-pass QF data structure.        */

extern void
  dwtalocal( 
	    interval *V,	/* Array of `L' empty INTERVALs.     */
	    interval *W,	/* Array of `L+1' empty INTERVALs.   */
	    interval *in,	/* Input INTERVAL.                   */
	    int L,		/* This is the number of levels.     */
	    const pqf *H,	/* Low-pass QF data structure.       */
	    const pqf *G);	/* High-pass QF data structure.      */

extern void
  idwtaintervals( 
		 interval *out, /* This is the output INTERVAL.  */
		 interval *V,	/* The `L' scaling subspaces.    */
		 interval *W,   /* The `L+1' wavelet subspaces.  */
		 int L,         /* This is the number of levels. */
		 const pqf *H,  /* Low-pass QF data structure.   */
		 const pqf *G); /* High-pass QF data structure.  */

extern void
  idwtaorigins( 
	       interval *V,	/* The `L' scaling subspaces.         */
	       real *sums,	/* This is the array of amplitudes.   */
	       int L);		/* This is the number of levels.      */

extern void
  idwta( 
	interval *out, /* Preallocated, zeroed output interval.     */
	interval *V,   /* Working intervals for scaling amplitudes. */
	interval *W,   /* Input intervals of wavelet amplitudes.    */
	int L,         /* This is the number of levels. */
	const pqf *H,  /* Low-pass QF data structure.   */
	const pqf *G); /* High-pass QF data structure.  */

#endif    /* DWT_HDR_ALREADY_INCLUDED */
