/*
 * Declare functions for tensor product D-dimensional discrete wavelet 
 * transforms and their inverses.
 * 
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 * 
 */

#ifndef WTENSOR_HDR_ALREADY_INCLUDED
# define WTENSOR_HDR_ALREADY_INCLUDED

#include "real.h"
#include "qf.h"

extern void
  dwtpt2(
	 real *data,		/* Input and output array.          */
	 int   ix,		/* Number of rows in `data[]'       */
	 int   iy,		/* Number of columns in `data[]'    */
	 int   xlevel,		/* Recursion depth in X-direction   */
	 int   ylevel,		/* Recursion depth in Y-direction   */
	 pqf  *HX,		/* Low-pass filter for X-direction  */
	 pqf  *GX,		/* High-pass filter for X-direction */
	 pqf  *HY,		/* Low-pass filter for Y-direction  */
	 pqf  *GY,		/* High-pass filter for Y-direction */
	 real *work);		/* Array with `max(ix,iy)' elements */

extern void
  idwtpt2(
	  real *data,		/* Input and output array.          */
	  int   ix,		/* Number of rows in `data[]'       */
	  int   iy,		/* Number of columns in `data[]'    */
	  int   xlevel,		/* Recursion depth in X-direction   */
	  int   ylevel,		/* Recursion depth in Y-direction   */
	  pqf  *HX,		/* Low-pass filter for X-direction  */
	  pqf  *GX,		/* High-pass filter for X-direction */
	  pqf  *HY,		/* Low-pass filter for Y-direction  */
	  pqf  *GY,		/* High-pass filter for Y-direction */
	  real *work);		/* Array with `max(ix,iy)' elements */

extern void
  dwtptd(
	 real *data,		/* Input and output array.                 */
	 int   d,		/* Number of dimensions in `data[]'        */
	 int  *lengths,		/* Length of `data[]' in each dimension    */
	 int  *levels,		/* Recursion depth in each dimension       */
	 pqf  *HQF,		/* Low-pass filter for each dimension      */
	 pqf  *GQF,		/* High-pass filter for each dimension     */
	 real *work);		/* Array of `max(lengths[k]: k=0,...,d-1)' */

extern void
  idwtptd(
	  real *data,		/* Input and output array.                 */
	  int   d,		/* Number of dimensions in `data[]'        */
	  int  *lengths,	/* Length of `data[]' in each dimension    */
	  int  *levels,		/* Recursion depth in each dimension       */
	  pqf  *HQF,		/* Low-pass filter for each dimension      */
	  pqf  *GQF,		/* High-pass filter for each dimension     */
	  real *work);		/* Array of `max(lengths[k]: k=0,...,d-1)' */

extern void
  dwtpvd(
	 real *data,		/* Input and output array.                 */
	 int   d,		/* Number of dimensions in `data[]'        */
	 int  *len,		/* Length of `data[]' in each dimension    */
	 int  *lvl,		/* Recursion depth in each dimension       */
	 pqf **QFS,		/* Low,high-pass filter for each dimension */
	 real *work);		/* Array of `max(len[k]: k=0,...,d-1)'     */

extern void
  idwtpvd(
	  real *data,		/* Input and output array.                 */
	  int   d,		/* Number of dimensions in `data[]'        */
	  int  *len,		/* Length of `data[]' in each dimension    */
	  int  *lvl,		/* Recursion depth in each dimension       */
	  pqf **QFS,		/* Low,high-pass filter for each dimension */
	  real *work);		/* Array of `max(len[k]: k=0,...,d-1)'     */

extern void
  dwtpi2(
	 real *data,		/* Input and output array.           */
	 int   ix,		/* Number of rows in `data[]'        */
	 int   iy,		/* Number of columns in `data[]'     */
	 int   lx,		/* Recursion depth in X-direction    */
	 int   ly,		/* Recursion depth in Y-direction    */
	 pqf  *HX,		/* Low-pass filter for X-direction   */
	 pqf  *GX,		/* High-pass filter for X-direction  */
	 pqf  *HY,		/* Low-pass filter for Y-direction   */
	 pqf  *GY,		/* High-pass filter for Y-direction  */
	 real *work);		/* Array with `max(ix,iy)' elements  */

extern void
  idwtpi2(
	  real *data,		/* Input and output array.          */
	  int   ix,		/* Number of rows in `data[]'       */
	  int   iy,		/* Number of columns in `data[]'    */
	  int   lx,		/* Recursion depth in X-direction   */
	  int   ly,		/* Recursion depth in Y-direction   */
	  pqf  *HX,		/* Low-pass filter for X-direction  */
	  pqf  *GX,		/* High-pass filter for X-direction */
	  pqf  *HY,		/* Low-pass filter for Y-direction  */
	  pqf  *GY,		/* High-pass filter for Y-direction */
	  real *work);		/* Array with `max(ix,iy)' elements */

#endif  /* WTENSOR_HDR_ALREADY_INCLUDED */
