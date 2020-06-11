/*
 *
 * Compute the 2-dimensional periodized tensor wavelet transform.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <assert.h>
#include "real.h"
#include "qf.h"
#include "xp.h"
#include "dwt.h"
#include "wtensor.h"

/**********************************************************************
 * dwtpt2()
 *
 * Tensor product of 2 1-dimensional periodized discrete wavelet bases.
 *
 * This function reads as input REALs representing a 2-dimensional
 * signal in the form of a 1-dimensional array of concatenated rows, and
 * writes its transform into coefficients with respect to the tensor product
 * of two 1-dimensional wavelet bases into a disjoint output array. The
 * output coefficients are not labeled, so they must be recognized by their
 * position in the output array.  For example, in a 3x,3y-level decomposition
 * the output is structured as follows:
 *         0 1 44 5555
 *         2 3 88 9999
 *         6 a cc dddd
 *         6 a cc dddd
 *         7 b ee ffff
 *         7 b ee ffff
 *         7 b ee ffff
 *         7 b ee ffff
 * Each number 0,1,2,3,4,5,6,7,8,9,a,b,c,d,e,f represents a different
 * shape of wavelet, from the product of the largest-scale scaling
 * functions (0) to the product of the smallest scale wavelets (f).  The
 * output is presented in the same order as the input, scanned in rows
 * from top to bottom and left to right.
 *
 * There can be as many levels in each dimension (X or Y) as the number
 * of times that dimension can be divided by 2.  It is not necessary to
 * use the same number of levels in both the X and Y, nor even to use the
 * same filter in both X and Y, but it is necessary to use the periodic
 * DWT to preserve the number of elements in each row and column.
 *
 * The QF's to use in the development are passed as the parameters `HX',
 * `GX', `HY' and `GY', which are pointers to predefined data structures
 * of type PQF used by the convolution-decimation function `cdpo()'.  They
 * should be periodized, perfect reconstruction,
 * finite impulse response conjugate orthogonal or biorthogonal quadrature
 * filters of the type described by I. Daubechies, in the paper
 * "Orthonormal bases of compactly supported wavelets", Communications
 * on Pure and Applied Mathematics, XLI(1988),909-996 and the book
 * "Ten Lectures on Wavelets", ISBN 0-89871-274-2, SIAM Press, 1993. 
 *
 *  Calling sequence and basic algorithm:
 *
 *   dwtpt2(DATA, IX,IY, XLEVEL,YLEVEL, HX,GX, HY,GY, WORK):
 *      Let I = 0
 *      While I < IY*IX
 *         dwtpi( DATA+I, WORK, IY, YLEVEL, HY, GY )
 *         I += IY
 *      xpi2( DATA, IX, IY )
 *      Let I = 0
 *      While I <  IY*IX
 *         dwtpi( DATA+I, WORK, IX, XLEVEL, HX, GX )
 *         I += IX
 *      xpi2( DATA, IY, IX )
 *
 *  Inputs:
 *	(real *)data		This must reference a pointer to a 
 *				  preallocated array with at least
 *				  `iy*ix' elements.
 *
 *	(int)ix			These positive integers give the number of
 *	(int)iy			  rows and columns in the array.
 *
 *	(int)xlevel		These must be nonnegative integers.  They
 *	(int)ylevel		   define the expansion depth in X and Y.
 *
 *	(pqf *)HX		These are pointers to the low-pass (s) and
 *	(pqf *)GX		  high-pass (d) quadrature filter structs
 *	(pqf *)HY		  used for convolution-decimation in the
 *	(pqf *)GY		  x-direction and y-direction, respectively.
 *
 *	(real *)work		This is a preallocated temporary array
 *				  with at least `max(iy,ix)' elements.
 *
 * Outputs:
 *	(real *)data		The transform overwrites the input.
 *
 * Assumptions:
 *	1. data != NULL
 *	2. work != NULL
 *	3. xlevel>=0
 *	4. ylevel>=0
 *	5. ix>0
 *	6. iy>0
 *
 * External functions called:
 *	dwtpi()	 Periodized 1-D in place discrete wavelet transform.
 *	xpi2()	 Transpose a 2-dimensional array in place.
 *	assert()
 */
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
	 real *work)		/* Array with `max(ix,iy)' elements */
{
  int i;

  assert(work);
  assert(data);
  assert(xlevel>=0);
  assert(ylevel>=0);
  assert(ix>0);
  assert(iy>0);

  /* Wavelet transform in place along rows: */
  for( i=0; i<iy*ix; i+=iy )
    dwtpi( data+i, work, iy, ylevel, HY, GY );

  /* Transposition in place: */
   xpi2( data, ix, iy, sizeof(real) );

  /* Wavelet transform in place along columns: */
  for( i=0; i<iy*ix; i+=ix )
      dwtpi( data+i, work, ix, xlevel, HX, GX );

  /* Transpose in place back to lists of rows: */
  xpi2( data, iy, ix, sizeof(real) );

  return;
}

/**********************************************************************
 * idwtpt2()
 *
 * Reconstruction from tensor product of 2 1-dimensional periodized DWTs.
 *
 * This function inverts `dwtpt2()' by reconstructing a periodic
 * signal from a 1-dimensional array of concatenated rows, considered 
 * to be the wavelet coefficients of the signal with respect to a tensor 
 * product of two 1-dimensional wavelet bases.  It writes the reconstructed
 * signal into a disjoint output array. The input coefficients must be
 * recognized by their position in the input array.  For example, in a
 * 3x,3y-level decomposition the input is structured as follows:
 *         0 1 44 5555
 *         2 3 88 9999
 *         6 a cc dddd
 *         6 a cc dddd
 *         7 b ee ffff
 *         7 b ee ffff
 *         7 b ee ffff
 *         7 b ee ffff
 * Each number 0,1,2,3,4,5,6,7,8,9,a,b,c,d,e,f represents a different
 * shape of wavelet, from the product of the largest-scale scaling
 * functions (0) to the product of the smallest scale wavelets (f).  The
 * output is presented in the same order as the input, scanned in rows
 * from top to bottom and left to right.
 *
 * There can be as many levels in each dimension (X or Y) as the number
 * of times that dimension can be divided by 2.  It is not necessary to
 * use the same number of levels in both the X and Y, nor even to use the
 * same filter in both X and Y, but it is necessary to use the periodic
 * DWT to preserve the number of elements in each row and column.
 *
 * The QF's to use in the development are passed as the parameters `HX',
 * `GX', `HY' and `GY', which are pointers to predefined data structures
 * of type PQF used by the adjoint convolution-decimation function `acdpo()'.
 * They should be periodized, perfect reconstruction,
 * finite impulse response conjugate orthogonal or biorthogonal quadrature
 * filters of the type described by I. Daubechies, in the paper
 * "Orthonormal bases of compactly supported wavelets", Communications
 * on Pure and Applied Mathematics, XLI(1988),909-996 and the book
 * "Ten Lectures on Wavelets", ISBN 0-89871-274-2, SIAM Press, 1993. 
 *
 *  Calling sequence and basic algorithm:
 *
 *   idwtpt2(DATA, IX,IY, XLEVEL,YLEVEL, HX,GX, HY,GY, WORK):
 *      Let I = 0
 *      While I < IY*IX
 *         idwtpi( DATA+I, WORK, IY, YLEVEL, HY, GY )
 *         I += IY
 *      xpi2( DATA, IX, IY )
 *      Let I = 0
 *      While I <  IY*IX
 *         idwtpi( DATA+I, WORK, IX, XLEVEL, HX, GX )
 *         I += IX
 *      xpi2( DATA, IY, IX )
 *
 *  Inputs:
 *	(real *)data		This must reference a pointer to a 
 *				  preallocated array with at least
 *				  `iy*ix' elements.
 *
 *	(int)ix			These positive integers give the number of
 *	(int)iy			  rows and columns in `data[]'.
 *
 *	(int)xlevel		These must be nonnegative integers.  They
 *	(int)ylevel		   define the expansion depth in X and Y.
 *
 *	(pqf *)HX		These are pointers to the low-pass and
 *	(pqf *)GX		  high-pass quadrature filter structs
 *	(pqf *)HY		  used for convolution-decimation in the
 *	(pqf *)GY		  x-direction and y-direction, respectively.
 *
 *	(real *)work		This is a preallocated temporary array
 *				  with at least `max(iy,ix)' elements.
 *
 * Outputs:
 *	(real *)data		The transform overwrites the input.
 *
 * Assumptions:
 *	1. data != NULL
 *	2. work != NULL
 *	3. xlevel>=0
 *	4. ylevel>=0
 *	5. ix>0
 *	6. iy>0
 *
 * External functions called:
 *	idwtpi() Periodized 1-D in place inverse DWT.
 *	xpi2()	 Transpose a 2-dimensional array in place.
 *	assert()
 */
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
	  real *work)		/* Array with `max(ix,iy)' elements */
{
  int i;

  assert(work);
  assert(data);
  assert(xlevel>=0);
  assert(ylevel>=0);
  assert(ix>0);
  assert(iy>0);

  /* Inverse wavelet transform in place along rows: */
  for( i=0; i<iy*ix; i+=iy )
    idwtpi( data+i, work, iy, ylevel, HY, GY );

  /* Transposition in place: */
   xpi2( data, ix, iy, sizeof(real) );

  /* Inverse wavelet transform in place along columns: */
  for( i=0; i<iy*ix; i+=ix )
    idwtpi( data+i, work, ix, xlevel, HX, GX );

  /* Transpose in place back to lists of rows: */
  xpi2( data, iy, ix, sizeof(real) );

  return;
}

/**********************************************************************
 * dwtptd()
 *
 * Tensor product of `d' 1-dimensional periodized discrete wavelet bases.
 *
 * This function reads as input REALs representing a `d'-dimensional
 * signal in the form of a 1-dimensional array of concatenated rows, and
 * writes its transform into coefficients with respect to the tensor product
 * of `d' 1-dimensional wavelet bases into a disjoint output array. The
 * output coefficients are not labeled, so they must be recognized by their
 * position in the output array.  This is a generalization of the
 * 2-dimensional algorithm `dwtptd()' defined above.   The
 * output is presented in the same order as the input, scanned in rows
 * from top to bottom and left to right.
 *
 * There can be as many levels in each dimension `0,1,...,d-1' as the
 * number of times that dimension can be divided by 2.  Different 
 * depths of decomposition can be used for each dimension.  It is,
 * however, necessary to use the periodic DWT to preserve the number
 * of elements in each row and column.
 *
 * The same filters `HQF' and `GQF' are used in the development of each
 * dimension.  They should be  periodized, perfect reconstruction,
 * finite impulse response conjugate orthogonal or biorthogonal quadrature
 * filters of the type described by I. Daubechies, in the paper
 * "Orthonormal bases of compactly supported wavelets", Communications
 * on Pure and Applied Mathematics, XLI(1988),909-996 and the book
 * "Ten Lectures on Wavelets", ISBN 0-89871-274-2, SIAM Press, 1993. 
 *
 *  Calling sequence and basic algorithm:
 *
 *   dwtptd( DATA, D, LENGTHS, LEVELS, HQF,GQF, WORK ):
 *     Let VOLUME = LENGTHS[0]
 *     For K = 1 to D-1
 *       VOLUME *= LENGTHS[K]
 *     For K = 1 to D
 *       Let I = 0
 *       While I < VOLUME
 *         dwtpi( DATA+I, WORK, LENGTHS[D-1], LEVELS[D-K], HQF,GQF )
 *         I += LENGTHS[D-1]
 *       xpid( DATA, LENGTHS, D )
 *
 *  Inputs:
 *	(real *)data		This must reference a pointer to a 
 *				  preallocated array with at least
 *				  `lengths[0]*...*lengths[d-1]' elements.
 *
 *	(int)d			This is the dimension of the array `data[]'
 *
 *	(int *)lengths		These positive integers give the dimensions
 *				  of the array of coefficients.
 *
 *	(int *)levels		This is an array of `d' nonnegative integers;
 *				  `levels[k]' defines the expansion depth
 *				   in dimension `k', for `k=0,...,d-1'.
 *
 *	(pqf *)HQF		These are pointers to the low-pass and
 *	(pqf *)GQF		  high-pass quadrature filter structs
 *				  used for convolution-decimation in the
 *				  various directions.
 *
 *	(real *)work		This is a preallocated temporary array
 *				  with `max(lengths[0],...,lengths[d-1])'
 *				  elements, used by `dwtpi()'.
 *
 * Outputs:
 *	(real *)data		The transform overwrites the input.
 *
 * Assumptions:
 *	1. 0<d<32
 *	2. data != NULL
 *	3. lengths != NULL
 *	4. levels != NULL
 *	5. work != NULL
 *	6. HQF != NULL
 *	7. GQF != NULL
 *
 * External functions called:
 *	dwtpi()	 Periodized 1-D in place discrete wavelet transform.
 *	xpid()	 Cyclic d-dimensional transposition in place.
 *	assert()
 */
extern void
  dwtptd(
	 real *data,		/* Input and output array.                 */
	 int   d,		/* Number of dimensions in `data[]'        */
	 int  *lengths,		/* Length of `data[]' in each dimension    */
	 int  *levels,		/* Recursion depth in each dimension       */
	 pqf  *HQF,		/* Low-pass filter for each dimension      */
	 pqf  *GQF,		/* High-pass filter for each dimension     */
	 real *work)		/* Array of `max(lengths[k]: k=0,...,d-1)' */
{
  long int i, volume;
  int k, level, row;

  assert(d>0);
  assert(d<32);
  assert(lengths);
  assert(levels);
  assert(data);
  assert(work);
  assert(HQF);
  assert(GQF);

  /* Compute the total length of `data[]': */
  volume = lengths[0];
  for( k=1;   k<d;  k++ )
    volume *= lengths[k];

  /* Loop over the dimensions: */
  for( k=d-1; k>=0; k-- )
    {
      row = lengths[d-1];
      level =  levels[k];

      for( i=0; i<volume; i+=row )
	dwtpi( data+i, work, row, level, HQF, GQF );

      /* Cyclic transposition: */
      xpid( data, lengths, d, sizeof(real) );
    }
  return;
}

/**********************************************************************
 * idwtptd()
 *
 * Inverse of tensor product of `d' 1-dimensional periodized DWTs.
 *
 * This function inverts `dwtptd()' by reconstructing a periodic
 * signal from a 1-dimensional array of concatenated rows, considered 
 * to be the wavelet coefficients of the signal with respect to a tensor 
 * product of `d' 1-dimensional wavelet bases.  It writes the reconstructed
 * signal into a disjoint output array. The input coefficients must be
 * recognized by their position in the input array.  This is a generalization
 * of the 2-dimensional algorithm `idwtpt2()' defined above.   The
 * output is presented in the same order as the input, scanned in rows
 * from top to bottom and left to right.
 *
 * There can be as many levels in each dimension `0,1,...,d-1' as the
 * number of times that dimension can be divided by 2.  Different 
 * depths of decomposition can be used for each dimension.  It is,
 * however, necessary to use the periodic DWT to preserve the number
 * of elements in each row and column.
 *
 * The same filters `HQF' and `GQF' are used in the reconstruction of each
 * dimension.  They should be  periodized, perfect reconstruction,
 * finite impulse response conjugate orthogonal or biorthogonal quadrature
 * filters of the type described by I. Daubechies, in the paper
 * "Orthonormal bases of compactly supported wavelets", Communications
 * on Pure and Applied Mathematics, XLI(1988),909-996 and the book
 * "Ten Lectures on Wavelets", ISBN 0-89871-274-2, SIAM Press, 1993. 
 *
 *  Calling sequence and basic algorithm:
 *
 *   idwtptd( DATA, D, LENGTHS, LEVELS, HQF,GQF, WORK ):
 *     Let VOLUME = LENGTHS[0]
 *     For K = 1 to D-1
 *       VOLUME *= LENGTHS[K]
 *     For K = 1 to D
 *       Let I = 0
 *       While I < VOLUME
 *         idwtpi( DATA+I, WORK, LENGTHS[D-1], LEVELS[D-K], HQF,GQF )
 *         I += LENGTHS[D-1]
 *       xpid( DATA, LENGTHS, D )
 *
 *  Inputs:
 *	(real *)data		This must reference a pointer to a 
 *				  preallocated array with at least
 *				  `lengths[0]*...*lengths[d-1]' elements.
 *
 *	(int)d			This is the dimension of the array `data[]'
 *
 *	(int *)lengths		These positive integers give the dimensions
 *				  of the array of coefficients.
 *
 *	(int *)levels		This is an array of `d' nonnegative integers;
 *				  `levels[k]' defines the expansion depth
 *				   in dimension `k', for `k=0,...,d-1'.
 *
 *	(pqf *)HQF		These are pointers to the low-pass and
 *	(pqf *)GQF		  high-pass quadrature filter structs
 *				  used for convolution-decimation in the
 *				  various directions.
 *
 *	(real *)work		This is a preallocated temporary array
 *				  with `max(lengths[0],...,lengths[d-1])'
 *				  elements, used by `idwtpi()'.
 *
 * Outputs:
 *	(real *)data		The transform overwrites the input.
 *
 * Assumptions:
 *	1. 0<d<32
 *	2. data != NULL
 *	3. lengths != NULL
 *	4. levels != NULL
 *	5. work != NULL
 *	6. HQF != NULL
 *	7. GQF != NULL
 *
 * External functions called:
 *	idwtpi() Periodized 1-D in place inverse discrete wavelet transform.
 *	xpid()	 Cyclic d-dimensional transposition in place.
 *	assert()
 */
extern void
  idwtptd(
	  real *data,		/* Input and output array.                 */
	  int   d,		/* Number of dimensions in `data[]'        */
	  int  *lengths,	/* Length of `data[]' in each dimension    */
	  int  *levels,		/* Recursion depth in each dimension       */
	  pqf  *HQF,		/* Low-pass filter for each dimension      */
	  pqf  *GQF,		/* High-pass filter for each dimension     */
	  real *work)		/* Array of `max(lengths[k]: k=0,...,d-1)' */
{
  long int i, volume;
  int k, level, row;

  assert(d>0);
  assert(d<32);
  assert(lengths);
  assert(levels);
  assert(data);
  assert(work);
  assert(HQF);
  assert(GQF);

  /* Compute the total length of `data[]': */
  volume = lengths[0];
  for( k=1;   k<d;  k++ )
    volume *= lengths[k];

  /* Loop over the dimensions: */
  for( k=d-1; k>=0; k-- )
    {
      row = lengths[d-1];
      level =  levels[k];

      for( i=0; i<volume; i+=row )
	idwtpi( data+i, work, row, level, HQF, GQF );

      /* Cyclic transposition: */
      xpid( data, lengths, d, sizeof(real) );
    }
  return;
}

/**********************************************************************
 * dwtpvd()
 *
 * Tensor product of `d' 1-dimensional periodized discrete wavelet bases,
 * with varying filters and levels of decomposition in each dimention.
 *
 * This function reads as input REALs representing a d-dimensional
 * signal in the form of a 1-dimensional array of concatenated rows, and
 * writes its transform into coefficients with respect to the tensor product
 * of `d' 1-dimensional wavelet bases into a disjoint output array. The
 * output coefficients are not labeled, so they must be recognized by their
 * position in the output array.  This is a generalization of the
 * algorithm `dwtptd()' defined above.   The output is presented in the
 * same order as the input, scanned in rows from top to bottom and left
 * to right.
 *
 * The QF's to use in the development are passed as adjacent elements
 * in an array of pointers, two for each dimension.  They
 * should be  periodized, perfect reconstruction,
 * finite impulse response conjugate orthogonal or biorthogonal quadrature
 * filters of the type described by I. Daubechies, in the paper
 * "Orthonormal bases of compactly supported wavelets", Communications
 * on Pure and Applied Mathematics, XLI(1988),909-996 and the book
 * "Ten Lectures on Wavelets", ISBN 0-89871-274-2, SIAM Press, 1993. 
 *
 * There can be as many levels in each dimension 0,1,...,d-1 as the number
 * of times that dimension can be divided by 2.  Different filters and
 * depths of decomposition can be used for each dimension.  It is,
 * however, necessary to use the periodic DWT to preserve the number
 * of elements in each row and column.
 *
 *  Calling sequence and basic algorithm:
 *
 *   dwtpvd( DATA, D, LEN, LVL, QFS, WORK ):
 *     Let VOLUME = LEN[0]
 *     For K = 1 to D-1
 *       VOLUME *= LEN[K]
 *     Let K = D-1
 *     While K >= 0
 *       Let I = 0
 *       While I < VOLUME
 *         dwtpi(DATA+I,WORK,LEN[D-1],LVL[K],QFS[2*K],QFS[2*K+1])
 *         I += LEN[D-1]
 *       xpid( DATA, LEN, D )
 *       K -= 1
 *
 *  Inputs:
 *	(real *)data		This must reference a pointer to a 
 *				  preallocated array with at least
 *				  `len[0]*...*len[d-1]' elements.
 *
 *	(int)d			This is the dimension of the array `data[]'
 *
 *	(int *)len		These `d' positive integers give the length
 *				  of the array in each dimension.
 *
 *	(int *)lvl		These `d' small positive integers give the
 *				  depth of decomposition in each dimension.
 *
 *	(pqf **)QFS		These are pointers to the low-pass and
 *				  high-pass quadrature filter structs
 *				  used for convolution-decimation in the
 *				  various directions.  `QFS[2*k]' is the
 *				  low-pass QF to use in dimension `k', and 
 *				  `QFS[2*k+1]' is the high pass filter
 *				  for that dimension.
 *
 *	(real *)work		This is a preallocated scratch array
 *				  with `max(len[0],...,len[d-1])'
 *				  elements, used by `dwtpi()'.
 *
 * Outputs:
 *	(real *)data		The transform overwrites the input.
 *
 * Assumptions:
 *	1. 0<d<32
 *	2. data != NULL
 *	3. len != NULL
 *	4. lvl != NULL
 *	5. work != NULL
 *	6. QFS != NULL
 *
 * External functions called:
 *	dwtpi()	 Periodized 1-D in place discrete wavelet transform.
 *	xpid()	 Cyclic d-dimensional transposition in place.
 *	assert()
 */
extern void
  dwtpvd(
	 real *data,		/* Input and output array.                 */
	 int   d,		/* Number of dimensions in `data[]'        */
	 int  *len,		/* Length of `data[]' in each dimension    */
	 int  *lvl,		/* Recursion depth in each dimension       */
	 pqf **QFS,		/* Low,high-pass filter for each dimension */
	 real *work)		/* Array of `max(len[k]: k=0,...,d-1)'     */
{
  long int i, volume;
  int level, row, k;
  pqf *HQF, *GQF;

  assert(d>0);
  assert(d<32);
  assert(len);
  assert(lvl);
  assert(data);
  assert(work);
  assert(QFS);

  /* Compute the total length of `data[]': */
  volume = len[0];
  for( k=1;   k<d;  k++ )
    volume *= len[k];

  /* Loop over the dimensions: */
  for( k=d-1; k>=0; k-- )
    {
      row = len[d-1];
      level =  lvl[k];
      HQF =  QFS[2*k];
      GQF =  QFS[2*k+1];

      for( i=0; i<volume; i+=row )
	dwtpi( data+i, work, row, level, HQF, GQF );

      /* Cyclic transposition: */
      xpid( data, len, d, sizeof(real) );
    }
  return;
}

/**********************************************************************
 * idwtpvd()
 *
 * Inverse of tensor product of `d' 1-dimensional periodized DWTs,
 * with varying filters and levels of decomposition in each dimention.
 *
 * This function inverts `dwtptd()' by reconstructing a periodic
 * signal from a 1-dimensional array of concatenated rows, considered 
 * to be the wavelet coefficients of the signal with respect to a tensor 
 * product of `d' 1-dimensional wavelet bases.  It writes the reconstructed
 * signal into a disjoint output array. The input coefficients must be
 * recognized by their position in the input array.  This is a generalization
 * of the fixed-filters algorithm `idwtptd()' defined above.   The
 * output is presented in the same order as the input, scanned in rows
 * from top to bottom and left to right.
 *
 * The QF's to use in the reconstruction are passed as adjacent elements
 * in an array of pointers, two for each dimension.  They
 * should be  periodized, perfect reconstruction,
 * finite impulse response conjugate orthogonal or biorthogonal quadrature
 * filters of the type described by I. Daubechies, in the paper
 * "Orthonormal bases of compactly supported wavelets", Communications
 * on Pure and Applied Mathematics, XLI(1988),909-996 and the book
 * "Ten Lectures on Wavelets", ISBN 0-89871-274-2, SIAM Press, 1993. 
 *
 * There can be as many levels in each dimension 0,1,...,d-1 as the number
 * of times that dimension can be divided by 2.  Different filters and
 * depths of decomposition can be used for each dimension.  It is,
 * however, necessary to use the periodic DWT to preserve the number
 * of elements in each row and column.
 *
 *  Calling sequence and basic algorithm:
 *
 *   idwtpvd( DATA, D, LEN, LVL, QFS, WORK ):
 *     Let VOLUME = LEN[0]
 *     For K = 1 to D-1
 *       VOLUME *= LEN[K]
 *     Let K = D-1
 *     While K >= 0
 *       Let I = 0
 *       While I < VOLUME
 *         idwtpi(DATA+I,WORK,LEN[D-1],LVL[K],QFS[2*K],QFS[2*K+1])
 *         I += LEN[D-1]
 *       xpid( DATA, LEN, D )
 *       K -= 1
 *
 *  Inputs:
 *	(real *)data		This must reference a pointer to a 
 *				  preallocated array with at least
 *				  `len[0]*...*len[d-1]' elements.
 *
 *	(int)d			This is the dimension of the array `data[]'
 *
 *	(int *)len		These `d' positive integers give the length
 *				  of the array in each dimension.
 *
 *	(int *)lvl		These `d' small positive integers give the
 *				  depth of decomposition in each dimension.
 *
 *	(pqf **)QFS		These are pointers to the low-pass and
 *				  high-pass quadrature filter structs
 *				  used for convolution-decimation in the
 *				  various directions.  `QFS[2*k]' is the
 *				  low-pass QF to use in dimension `k', and 
 *				  `QFS[2*k+1]' is the high pass filter
 *				  for that dimension.
 *
 *	(real *)work		This is a preallocated scratch array
 *				  with `max(len[0],...,len[d-1])'
 *				  elements, used by `idwtpi()'.
 *
 * Outputs:
 *	(real *)data		The transform overwrites the input.
 *
 * Assumptions:
 *	1. 0<d<32
 *	2. data != NULL
 *	3. len != NULL
 *	4. lvl != NULL
 *	5. work != NULL
 *	6. QFS != NULL
 *
 * External functions called:
 *	idwtpi() Periodized 1-D in place inverse discrete wavelet transform.
 *	xpid()	 Cyclic d-dimensional transposition in place.
 *	assert()
 */
extern void
  idwtpvd(
	  real *data,		/* Input and output array.                 */
	  int   d,		/* Number of dimensions in `data[]'        */
	  int  *len,		/* Length of `data[]' in each dimension    */
	  int  *lvl,		/* Recursion depth in each dimension       */
	  pqf **QFS,		/* Low,high-pass filter for each dimension */
	  real *work)		/* Array of `max(len[k]: k=0,...,d-1)'     */
{
  long int i, volume;
  int level, row, k;
  pqf *HQF, *GQF;

  assert(d>0);
  assert(d<32);
  assert(len);
  assert(lvl);
  assert(data);
  assert(work);
  assert(QFS);

  /* Compute the total length of `data[]': */
  volume = len[0];
  for( k=1;   k<d;  k++ )
    volume *= len[k];

  /* Loop over the dimensions: */
  for( k=d-1; k>=0; k-- )
    {
      row = len[d-1];
      level =  lvl[k];
      HQF =  QFS[2*k];
      GQF =  QFS[2*k+1];

      for( i=0; i<volume; i+=row )
	idwtpi( data+i, work, row, level, HQF, GQF );

      /* Cyclic transposition: */
      xpid( data, len, d, sizeof(real) );
    }
  return;
}

/**********************************************************************
 * dwtpi2()
 *
 * Isotropic in place two-dimensional periodized discrete wavelet transform.
 *
 * This function reads as input REALs representing a 2-dimensional
 * signal in the form of a 1-dimensional array of concatenated rows, and
 * replaces it with its coefficients with respect to the isotropic
 * two-dimensional wavelet bases.  The output coefficients are not labeled,
 * so they must be recognized by their position in the output array.  The
 * output is presented in the same order as the input, scanned in rows
 * from top to bottom and left to right.
 *
 * There can be as many levels in each dimension (X or Y) as the number
 * of times that dimension can be divided by 2.  It is not necessary to
 * use the same number of levels in both the X and Y, nor even to use the
 * same filter in both X and Y, but it is necessary to use the periodic
 * DWT to preserve the number of elements in each row and column.
 *
 * The QF's to use in the development are passed as the parameters `HX',
 * `GX', `HY' and `GY', which are pointers to predefined data structures
 * of type PQF used by the convolution-decimation function `cdpo()'.  They
 * should be periodized, perfect reconstruction,
 * finite impulse response conjugate orthogonal or biorthogonal quadrature
 * filters of the type described by I. Daubechies, in the paper
 * "Orthonormal bases of compactly supported wavelets", Communications
 * on Pure and Applied Mathematics, XLI(1988),909-996 and the book
 * "Ten Lectures on Wavelets", ISBN 0-89871-274-2, SIAM Press, 1993. 
 *
 *  Calling sequence and basic algorithm:
 *
 *   dwtpi2(DATA, IX,IY, LX,LY, HX,GX, HY,GY, WORK):
 *      Let M = IX>>LX
 *      Let I = 0
 *      While I < M*IY
 *         dwtpi( DATA, WORK, IY, LY, HY, GY )
 *         I += IY
 *      While M < IX
 *         While I < 2*M*IY
 *            dwtpi( DATA, WORK, IY, LY, HY, GY )
 *            I += IY
 *         M *= 2
 *      xpi2( DATA, IX, IY )
 *      Let M = IY>>LY
 *      Let I = 0
 *      While I < M*IX
 *         dwtpi( DATA, WORK, IX, LX, HX, GX )
 *         I += IX
 *      While M < IY
 *         While I < 2*M*IX
 *            dwtpi( DATA, WORK, IX, LX, HX, GX )
 *            I += IX
 *         M *= 2
 *      xpi2( DATA, IY, IX )
 *
 *  Inputs:
 *	(real *)data		This must reference a pointer to a 
 *				  preallocated array with at least
 *				  `iy*ix' elements.
 *
 *	(int)ix			These positive integers give the number of rows
 *	(int)iy			  and columns in the array of coefficients.
 *
 *	(int)lx			These must be nonnegative integers.  They
 *	(int)ly			   define the expansion depth in X and Y.
 *
 *	(pqf *)HX		These are pointers to the low-pass and
 *	(pqf *)GX		  high-pass quadrature filter structs
 *	(pqf *)HY		  used for convolution-decimation in the
 *	(pqf *)GY		  x-direction and y-direction, respectively.
 *
 *	(real *)work		This is a preallocated temporary array
 *				  with at least `max(iy,ix)' elements.
 *
 * Outputs:
 *	(real *)data		The transform overwrites the input.
 *
 * Assumptions:
 *	1. data != NULL
 *	2. work != NULL
 *	3. lx>=0
 *	4. ly>=0
 *	5. `ix' is divisible by `1<<lx'
 *	6. `iy' is divisible by `1<<ly'
 *
 * External functions called:
 *	dwtpi()	 Periodized 1-D in place discrete wavelet transform.
 *	xpi2()	 Transpose a 2-dimensional array in place.
 *	assert()
 */
extern void
  dwtpi2(
	 real *data,		/* Input and output array.          */
	 int   ix,		/* Number of rows in `data[]'       */
	 int   iy,		/* Number of columns in `data[]'    */
	 int   lx,		/* Recursion depth in X-direction   */
	 int   ly,		/* Recursion depth in Y-direction   */
	 pqf  *HX,		/* Low-pass filter for X-direction  */
	 pqf  *GX,		/* High-pass filter for X-direction */
	 pqf  *HY,		/* Low-pass filter for Y-direction  */
	 pqf  *GY,		/* High-pass filter for Y-direction */
	 real *work)		/* Array with `max(ix,iy)' elements */
{
  int i, m;

  assert(work);
  assert(data);
  assert(lx>=0);
  assert(ly>=0);
  assert(ix>0);
  assert(iy>0);

  /* DWTP along the rows (Y-direction): */
  m = ix>>lx;
  for( i=0; i<m*iy; i+=iy )
    dwtpi( data, work, iy, ly, HY, GY );

  for( ; m<ix; m*=2 )
    for( ; i<2*m*iy; i+=iy )
      dwtpi( data, work, iy, ly, HY, GY );

  /* Tranpose in place: */
  xpi2( data, ix, iy, sizeof(real) );

  /* DWTP along the columns (X-direction): */
  m = iy>>ly;
  for( i=0;  i<m*ix; i+=ix )
    dwtpi( data, work, ix, lx, HX, GX );

  for( ; m<iy; m*=2 )
    for( ; i<2*m*ix; i+=ix )
      dwtpi( data, work, ix, lx, HX, GX );

  /* Transpose in place back to original orientation: */
   xpi2( data, iy, ix, sizeof(real) );

  return;
}

/**********************************************************************
 * idwtpi2()
 *
 * Isotropic in place two-dimensional inverse periodized DWT.
 *
 * This function inverts `dwtpi2()' by reconstructing a periodic
 * signal from a 1-dimensional array of concatenated rows, considered 
 * to be the wavelet coefficients of the signal with respect to an
 * isotropic two-dimensional wavelet transform.  It writes the reconstructed
 * signal into a disjoint output array. The input coefficients must be
 * recognized by their position in the input array.  The output is
 * presented in the same order as the input, scanned in rows from top
 * to bottom and left to right.
 *
 * There can be as many levels in each dimension (X or Y) as the number
 * of times that dimension can be divided by 2.  It is not necessary to
 * use the same number of levels in both the X and Y, nor even to use the
 * same filter in both X and Y, but it is necessary to use the periodic
 * iDWT to preserve the number of elements in each row and column.
 *
 * The QF's to use in the development are passed as the parameters `HX',
 * `GX', `HY' and `GY', which are pointers to predefined data structures
 * of type PQF used by the adjoint convolution-decimation function `acdpo()'.
 * They should be periodized, perfect reconstruction,
 * finite impulse response conjugate orthogonal or biorthogonal quadrature
 * filters of the type described by I. Daubechies, in the paper
 * "Orthonormal bases of compactly supported wavelets", Communications
 * on Pure and Applied Mathematics, XLI(1988),909-996 and the book
 * "Ten Lectures on Wavelets", ISBN 0-89871-274-2, SIAM Press, 1993. 
 *
 *  Calling sequence and basic algorithm:
 *
 *   idwtpi2(DATA, IX,IY, LX,LY, HX,GX, HY,GY, WORK):
 *      Let M = IX>>LX
 *      Let I = 0
 *      While I < M*IY
 *         idwtpi( DATA, WORK, IY, LY, HY, GY )
 *         I += IY
 *      While M < IX
 *         While I < 2*M*IY
 *            idwtpi( DATA, WORK, IY, LY, HY, GY )
 *            I += IY
 *         M *= 2
 *      xpi2( DATA, IX, IY )
 *      Let M = IY>>LY
 *      Let I = 0
 *      While I < M*IX
 *         idwtpi( DATA, WORK, IX, LX, HX, GX )
 *         I += IX
 *      While M < IY
 *         While I < 2*M*IX
 *            idwtpi( DATA, WORK, IX, LX, HX, GX )
 *            I += IX
 *         M *= 2
 *      xpi2( DATA, IY, IX )
 *
 *  Inputs:
 *	(real *)data		This must reference a pointer to a 
 *				  preallocated array with at least
 *				  `iy*ix' elements.
 *
 *	(int)ix			These positive integers give the number of rows
 *	(int)iy			  and columns in the array of coefficients.
 *
 *	(int)lx			These must be nonnegative integers.  They
 *	(int)ly			   define the expansion depth in X and Y.
 *
 *	(pqf *)HX		These are pointers to the low-pass and
 *	(pqf *)GX		  high-pass quadrature filter structs
 *	(pqf *)HY		  used for convolution-decimation in the
 *	(pqf *)GY		  x-direction and y-direction, respectively.
 *
 *	(real *)work		This is a preallocated temporary array
 *				  with at least `max(iy,ix)' elements.
 *
 * Outputs:
 *	(real *)data		The transform overwrites the input.
 *
 * Assumptions:
 *	1. data != NULL
 *	2. work != NULL
 *	3. lx>=0
 *	4. ly>=0
 *	5. `ix' is divisible by `1<<lx'
 *	6. `iy' is divisible by `1<<ly'
 *
 * External functions called:
 *	idwtpi() Periodized 1-D in place inverse discrete wavelet transform.
 *	xpi2()	 Transpose a 2-dimensional array in place.
 *	assert()
 */
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
	  real *work)		/* Array with `max(ix,iy)' elements */
{
  int i, m;

  assert(work);
  assert(data);
  assert(lx>=0);
  assert(ly>=0);
  assert(ix>0);
  assert(iy>0);

  /* iDWTP along the rows (Y-direction): */
  m = ix>>lx;
  for( i=0; i<m*iy; i+=iy )
    idwtpi( data, work, iy, ly, HY, GY );

  for( ; m<ix; m*=2 )
    for( ; i<2*m*iy; i+=iy )
      idwtpi( data, work, iy, ly, HY, GY );

  /* Tranpose in place: */
  xpi2( data, ix, iy, sizeof(real) );

  /* iDWTP along the columns (X-direction): */
  m = iy>>ly;
  for( i=0;  i<m*ix; i+=ix )
    idwtpi( data, work, ix, lx, HX, GX );

  for( ; m<iy; m*=2 )
    for( ; i<2*m*ix; i+=ix )
      idwtpi( data, work, ix, lx, HX, GX );

  /* Transpose in place back to original orientation: */
   xpi2( data, iy, ix, sizeof(real) );

  return;
}
