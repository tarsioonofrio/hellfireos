/*
 * Declare separable 2-dimensional convolution-decimation functions.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef CDP2D_HDR_ALREADY_INCLUDED
# define CDP2D_HDR_ALREADY_INCLUDED

#include "real.h"
#include "qf.h"		/* Declaration of PQF data structure. */

extern void
  sacdpd2(
	  real *out,		/* Pointer to a preallocated output array    */
	  const real *in0,	/* Pointer to the (HxHy) input array or NULL */
	  const real *in1,	/* Pointer to the (HxGy) input array or NULL */
	  const real *in2,	/* Pointer to the (GxHy) input array or NULL */
	  const real *in3,	/* Pointer to the (GxGy) input array or NULL */
	  int    ix,		/* Rows in each of `in0[]...in3[]'           */
	  int    iy,		/* Columns in each of `in0[]...in3[]'        */
	  real *work,		/* Scratch array of size `2*ix*iy'           */
	  const pqf *HQF,	/* 1-D Adjoint low-pass quadrature filter    */
	  const pqf *GQF);	/* 1-D Adjoint high-pass quadrature filter   */

extern void
  sacdpe2(
	  real *out,		/* Pointer to a preallocated output array    */
	  const real *in0,	/* Pointer to the (HxHy) input array or NULL */
	  const real *in1,	/* Pointer to the (HxGy) input array or NULL */
	  const real *in2,	/* Pointer to the (GxHy) input array or NULL */
	  const real *in3,	/* Pointer to the (GxGy) input array or NULL */
	  int    ix,		/* Rows in each of `in0[]...in3[]'           */
	  int    iy,		/* Columns in each of `in0[]...in3[]'        */
	  real *work,		/* Scratch array of size `2*ix*iy'           */
	  const pqf *HQF,	/* 1-D Adjoint low-pass quadrature filter    */
	  const pqf *GQF);	/* 1-D Adjoint high-pass quadrature filter   */

extern void
  sacdpi2(
	  real *data,		/* Pointer to the input and output array   */
	  int   ix,		/* Rows in the 4 parts of `data[]'         */
	  int   iy,		/* Columns in the 4 parts of `data[]'      */
	  real *work,		/* Scratch array of size `4*ix*iy'         */
	  const pqf *HQF,	/* 1-D Adjoint low-pass quadrature filter  */
	  const pqf *GQF);	/* 1-D Adjoint high-pass quadrature filter */

extern void
  scdpd2(
	 real *out0,		/* Pointer to output (HxHy) array, or NULL  */
	 real *out1,		/* Pointer to output (GxHy) array, or NULL  */
	 real *out2,		/* Pointer to output (HxGy) array, or NULL  */
	 real *out3,		/* Pointer to output (GxGy) array, or NULL  */
	 const real *in,	/* Pointer to the input array: rows first   */
	 int  ix,		/* Number of rows in the array `in[]'       */
	 int  iy,		/* Number of columns in the array `in[]'    */
	 real *work,		/* Scratch array of size `ix*iy/2'          */
	 const pqf *HQF,	/* Low-pass quadrature filter data struct   */
	 const pqf *GQF);	/* High-pass quadrature filter data struct  */

extern void
  scdpe2(
	 real *out0,		/* Pointer to output (HxHy) array, or NULL  */
	 real *out1,		/* Pointer to output (HxGy) array, or NULL  */
	 real *out2,		/* Pointer to output (GxHy) array, or NULL  */
	 real *out3,		/* Pointer to output (GxGy) array, or NULL  */
	 const real *in,	/* Pointer to the input array: rows first   */
	 int  ix,		/* Number of rows in the array `in[]'       */
	 int  iy,		/* Number of columns in the array `in[]'    */
	 real *work,		/* Scratch array of size `ix*iy/2'          */
	 const pqf *HQF,	/* Low-pass quadrature filter data struct   */
	 const pqf *GQF);	/* High-pass quadrature filter data struct  */

extern void
  scdpi2(
	 real *data,		/* Pointer to joint input and output array  */
	 int  ix,		/* Number of rows in the array `data[]'     */
	 int  iy,		/* Number of columns in the array `data[]'  */
	 real *work,		/* Scratch array of size `ix*iy'            */
	 const pqf *HQF,	/* Low-pass quadrature filter data struct   */
	 const pqf *GQF);	/* High-pass quadrature filter data struct  */

#endif /* CDP2D_HDR_ALREADY_INCLUDED */
