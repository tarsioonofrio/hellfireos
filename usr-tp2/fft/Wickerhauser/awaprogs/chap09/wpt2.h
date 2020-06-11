/*
 * Declare functions to perform the memory-saving 2-dimensional discrete
 * wavelet packet transform to a custom graph basis or the best graph basis.
 * 
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 * 
 */

#ifndef WPT2_HDR_ALREADY_INCLUDED
# define WPT2_HDR_ALREADY_INCLUDED

#include "real.h"
#include "qf.h"
#include "hedge.h"

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
	pqf   *G);		/* High-pass periodic quadrature filter    */

extern void
  cbwp2(
	hedge *graph,		/* Preallocated struct for the best basis  */
	int    level,		/* Current level in the DWPA tree          */
	int    ix,		/* Length of one column of `in[]'          */
	int    iy,		/* Length of one row of `in[]'             */
	real  *work,		/* Scratch array for transposition         */
	pqf   *H,		/* Low-pass periodic quadrature filter     */
	pqf   *G);		/* High-pass periodic quadrature filter    */

#endif /* WPT2_HDR_ALREADY_INCLUDED */

