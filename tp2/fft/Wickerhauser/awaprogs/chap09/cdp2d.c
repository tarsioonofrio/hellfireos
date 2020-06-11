/*
 * 
 * Separable 2-D convolution-decimation routines.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <assert.h>
#include "real.h"
#include "qf.h"
#include "cd.h"
#include "cdp2d.h"

/**********************************************************************
 * sacdpd2()
 *
 *  [S]eparable [A]djoint [C]onvolution-[D]ecimation;
 *  [P]eriodic, [D]isjoint, [2]-dimensional.
 *
 *  This function inverts the separable 2-dimensional periodic 
 *  convolution-decimation performed by `scdpd2()' or `scdpi2()',
 *  using the quadrature mirror filters applied by the function
 *  `acdpo()'.  It superposes the results into the array `out[]'.
 *  It takes 4 input arrays as arguments.  It assumes that the
 *  input and output are disjoint.
 *
 *  Note: if any of the 4 output arrays are NULL, then they are
 *  not included in the reconstruction.  This provides a means of
 *  cheaply including totally zero descendents.
 *
 *  Calling sequence and basic algorithm:
 *
 *   sacdpd2( OUT, IN0,IN1,IN2,IN3, IX,IY, WORK, HQF,GQF ):
 *      Let OY = 2*IY
 *      For I=0 to IX-1
 *         acdpe( WORK+I, IX, IN0+I*IY, IY, HQF )
 *         acdpo( WORK+I, IX, IN1+I*IY, IY, GQF )
 *      For I=0 to OY-1
 *         acdpo( OUT+I, OY, WORK+I*IX, IX, HQF )
 *      For I=0 to IX-1
 *         acdpe( WORK+I, IX, IN2+I*IY, IY, HQF )
 *         acdpo( WORK+I, IX, IN3+I*IY, IY, GQF )
 *      For I=0 to OY-1
 *         acdpo( OUT+I, OY, WORK+I*IX, IX, GQF )
 *
 *
 *  Inputs:
 *	(real *)out		This array must be preallocated with at 
 *				  least `4*iy*ix' elements.
 *				  It will be overwritten even before all
 *				  the arrays `in0[]...in3[]' have been read.
 *
 *	(const real *)in0	These arrays must be preallocated and defined
 *	(const real *)in1	  with at least `iy*ix' locations each,
 *	(const real *)in2	  or else be NULL.  They are not changed by
 *	(const real *)in3	  this function.
 *
 *	(int)ix			These positive integers are the number of rows
 *	(int)iy			  and columns in the arrays `in0[]...in3[]'
 *
 *	(real *)work		This temporary array must be preallocated
 *				  with at least `2*ix*iy' elements.  
 *				  It will be overwritten even before all the
 *				  arrays `in0[]...in3[],' have been read.
 *
 *	(const pqf *)HQF	These point to structs containing filter
 *	(const pqf *)GQF	  coefficients for adjoint convolution-
 *				  decimation, used in reconstruction.
 *
 *  Outputs:
 *	(real *)out		This array is filled by side effect.
 *
 *	(real *)work		This array is trashed by side effect.
 *
 *  External functions called:
 *	acdpo(), acdpe(), assert
 *
 *  Assumptions:
 *	1. ix>0
 *	2. iy>0
 *	3. HQF != NULL
 *	4. GQF != NULL
 *	5. `out[0,...,4*ix*iy-1]' is disjoint from each of 
 *		`in0[0,...,ix*iy-1]',...,`in3[0,...,ix*iy-1]'.
 *	6. work!=NULL
 *	7. out!=NULL 
 *	8. in0!=NULL || in1!=NULL || in2!=NULL || in3!=NULL
 */
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
	  const pqf *GQF)	/* 1-D Adjoint high-pass quadrature filter   */
{
  int i, oy;
  real *wptr;

  assert(ix>0);
  assert(iy>0);
  assert(HQF);
  assert(GQF);
  assert(work);
  assert(out);
  assert( in0 || in1 || in2 || in3 );

  oy = 2*iy;

  if( in0 || in1 )
    {
      /* Y-direction 1: Loop over the rows of `in0[]' and/or `in1[]': */
      if( in0 )
	{
	  for( i=0; i<ix; i++ )
	    {
	      acdpe ( work+i, ix, in0, iy, HQF );
	      in0 += iy;
	    }
	  if( in1 )
	    for( i=0; i<ix; i++ )
	      {
		acdpo ( work+i, ix, in1, iy, GQF );
		in1 += iy;
	      }
	}
      else			/* `in0==NULL' but `in1 != NULL' */
	for( i=0; i<ix; i++ )
	  {
	    acdpe ( work+i, ix, in1, iy, GQF );
	    in1 += iy;
	  }

      /* X-direction 1: Loop over the rows of `work[]': */
      wptr = work;
      for( i=0; i<oy; i++ )
	{
	  acdpo ( out+i, oy, wptr, ix, HQF );
	  wptr += ix;
	}
    }

  if( in2 || in3 )
    {
      /* Y-direction 2: Loop over the rows of `in2[]' and/or `in3[]': */
      if( in2 )
	{
	  for( i=0; i<ix; i++ )
	    {
	      acdpe ( work+i, ix, in2, iy, HQF );
	      in2 += iy;
	    }
	  if( in3 )
	    for( i=0; i<ix; i++ )
	      {
		acdpo ( work+i, ix, in3, iy, GQF );
		in3 += iy;
	      }
	}
      else			/* `in2==NULL' but `in3 != NULL' */
	for( i=0; i<ix; i++ )
	  {
	    acdpe ( work+i, ix, in3, iy, GQF );
	    in3 += iy;
	  }
      
      /* X-direction 2: Loop over the rows of `work[]': */
      wptr = work;
      for(i=0; i<oy; i++)
	{
	  acdpo ( out+i, oy, wptr, ix, GQF );
	  wptr += ix;
	}
    }
  return;
}



/**********************************************************************
 * sacdpe2()
 *
 *  [S]eparable [A]djoint [C]onvolution-[D]ecimation;
 *  [P]eriodic, [E]quals, [2]-dimensional.
 *
 *  This function inverts the separable 2-dimensional periodic 
 *  convolution-decimation performed by `scdpd2()' or `scdpi2()',
 *  using the quadrature mirror filters applied by the function
 *  `acdpo()'.  It assigns the results into the array `out[]'.
 *  It takes 4 input arrays as arguments.  It assumes that the
 *  input and output are disjoint.
 *
 *  Note: if any of the 4 output arrays are NULL, then they are
 *  not included in the reconstruction.  This provides a means of
 *  cheaply including totally zero descendents.
 *
 *  Calling sequence and basic algorithm:
 *
 *   sacdpe2( OUT, IN0,IN1,IN2,IN3, IX,IY, WORK, HQF,GQF ):
 *      Let OY = 2*IY
 *      For I=0 to IX-1
 *         acdpe( WORK+I, IX, IN0+I*IY, IY, HQF )
 *         acdpo( WORK+I, IX, IN1+I*IY, IY, GQF )
 *      For I=0 to OY-1
 *         acdpe( OUT+I, OY, WORK+I*IX, IX, HQF )
 *      For I=0 to IX-1
 *         acdpe( WORK+I, IX, IN2+I*IY, IY, HQF )
 *         acdpo( WORK+I, IX, IN3+I*IY, IY, GQF )
 *      For I=0 to OY-1
 *         acdpo( OUT+I, OY, WORK+I*IX, IX, GQF )
 *
 *
 *  Inputs:
 *	(real *)out		This array must be preallocated with at 
 *				  least `4*iy*ix' elements.
 *				  It will be overwritten even before all
 *				  the arrays `in0[]...in3[]' have been read.
 *
 *	(const real *)in0	These arrays must be preallocated and defined
 *	(const real *)in1	  with at least `iy*ix' locations each,
 *	(const real *)in2	  or else be NULL.  They are not changed by
 *	(const real *)in3	  this function.
 *
 *	(int)ix			These positive integers are the number of rows
 *	(int)iy			  and columns in the arrays `in0[]...in3[]'
 *
 *	(real *)work		This temporary array must be preallocated
 *				  with at least `2*ix*iy' elements.  
 *				  It will be overwritten even before all the
 *				  arrays `in0[]...in3[],' have been read.
 *
 *	(const pqf *)HQF	These point to structs containing filter
 *	(const pqf *)GQF	  coefficients for adjoint convolution-
 *				  decimation, used in reconstruction.
 *
 *  Outputs:
 *	(real *)out		This array is filled by side effect.
 *
 *	(real *)work		This array is trashed by side effect.
 *
 *  External functions called:
 *	acdpe(), acdpo(), assert()
 *
 *  Assumptions:
 *	1. ix>0
 *	2. iy>0
 *	3. HQF != NULL
 *	4. GQF != NULL
 *	5. `out[0,...,4*ix*iy-1]' is disjoint from each of 
 *		`in0[0,...,ix*iy-1]',...,`in3[0,...,ix*iy-1]'.
 *	6. work!=NULL
 *	7. out!=NULL 
 *	8. in0!=NULL || in1!=NULL || in2!=NULL || in3!=NULL
 */
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
	  const pqf *GQF)	/* 1-D Adjoint high-pass quadrature filter   */
{
  int i, oy;
  real *wptr;

  assert(ix>0);
  assert(iy>0);
  assert(HQF);
  assert(GQF);
  assert(work);
  assert(out);
  assert( in0 || in1 || in2 || in3 );

  oy = 2*iy;

  if( in0 || in1 )
    {
      /* Y-direction 1: Loop over the rows of `in0[]' and/or `in1[]': */
      if( in0 )
	{
	  for( i=0; i<ix; i++ )
	    {
	      acdpe ( work+i, ix, in0, iy, HQF );
	      in0 += iy;
	    }
	  if( in1 )
	    for( i=0; i<ix; i++ )
	      {
		acdpo ( work+i, ix, in1, iy, GQF );
		in1 += iy;
	      }
	}
      else			/* `in0==NULL' but `in1 != NULL' */
	for( i=0; i<ix; i++ )
	  {
	    acdpe ( work+i, ix, in1, iy, GQF );
	    in1 += iy;
	  }

      /* X-direction 1: Loop over the rows of `work[]': */
      wptr = work;
      for( i=0; i<oy; i++ )
	{
	  acdpe ( out+i, oy, wptr, ix, HQF );
	  wptr += ix;
	}
    }
  else				/*  The output array should be all zero: */
    for( i=0; i<4*ix*iy; i++ )
      out[i] = 0;

  if( in2 || in3 )
    {
      /* Y-direction 2: Loop over the rows of `in2[]' and/or `in3[]': */
      if( in2 )
	{
	  for( i=0; i<ix; i++ )
	    {
	      acdpe ( work+i, ix, in2, iy, HQF );
	      in2 += iy;
	    }
	  if( in3 )
	    for( i=0; i<ix; i++ )
	      {
		acdpo ( work+i, ix, in3, iy, GQF );
		in3 += iy;
	      }
	}
      else			/* `in2==NULL' but `in3 != NULL' */
	for( i=0; i<ix; i++ )
	  {
	    acdpe ( work+i, ix, in3, iy, GQF );
	    in3 += iy;
	  }
      
      /* X-direction 2: Loop over the rows of `work[]': */
      wptr = work;
      for(i=0; i<oy; i++)
	{
	  acdpo ( out+i, oy, wptr, ix, GQF );
	  wptr += ix;
	}
    }
  return;
}

/**********************************************************************
 * sacdpi2()
 *
 *  [S]eparable [A]djoint [C]onvolution-[D]ecimation;
 *  [P]eriodic, [I]n-place, [2]-dimensional.
 *
 *  This function inverts the separable 2-dimensional periodic 
 *  convolution-decimation performed by `scdpi2()', using the
 *  quadrature mirror filters applied by the functions `acdpo()'
 *  and `acdpe()' to a single input/output array.   It reads all
 *  of the input before writing any output, which permits it to perform
 *  its transformation in place.  After reading all input, it zeros the
 *  input/output array and then assigns it the computed filtered values.
 *
 *  Calling sequence and basic algorithm:
 *
 *   sacdpi2( DATA, IX,IY, WORK, HQF,GQF ):
 *      Let OY = 2*IY
 *      Let N = IY*IX
 *      Let WORK1 = WORK
 *      Let WORK2 = WORK + OY*IX
 *      For I=0 to IX-1
 *         acdpe( WORK1+I, IX, DATA  +  I*IY, IY, HQF )
 *         acdpo( WORK1+I, IX, DATA+ N +I*IY, IY, GQF )
 *         acdpe( WORK2+I, IX, DATA+2*N+I*IY, IY, HQF )
 *         acdpo( WORK2+I, IX, DATA+3*N+I*IY, IY, GQF )
 *      For I=0 to OY-1
 *         acdpe( DATA+I, OY, WORK1+I*IX, IX, HQF )
 *         acdpo( DATA+I, OY, WORK2+I*IX, IX, GQF )
 *
 *
 *  Inputs:
 *	(real *)data		This array must be preallocated with at 
 *				  least `(2*iy)*(2*ix)' elements.
 *
 *	(int)ix			These positive integers are the number of rows
 *	(int)iy			  and columns in the 4 subarrays of `data[]'
 *
 *	(real *)work		This temporary array must be preallocated
 *				  with at least `4*ix*iy' elements.  
 *
 *	(const pqf *)HQF	These point to structs containing filter
 *	(const pqf *)GQF	  coefficients for adjoint convolution-
 *				  decimation, used in reconstruction.
 *
 *  Outputs:
 *	(real *)data		This array is replaced by side effect.
 *
 *	(real *)work		This array is trashed by side effect.
 *
 *  Assumptions:
 *	1. data != NULL
 *	2. `iy' and `ix' are positive.
 *	3. HQF != NULL
 *	4. GQF != NULL
 *	5. work != NULL
 *
 *  External functions called:
 *	acdpe(), acdpo(), assert()
 */
extern void
  sacdpi2(
	  real *data,		/* Pointer to the input and output array   */
	  int   ix,		/* Rows in the 4 parts of `data[]'         */
	  int   iy,		/* Columns in the 4 parts of `data[]'      */
	  real *work,		/* Scratch array of size `4*ix*iy'         */
	  const pqf *HQF,	/* 1-D Adjoint low-pass quadrature filter  */
	  const pqf *GQF)	/* 1-D Adjoint high-pass quadrature filter */
{
  int i, oy, n;
  real *work1, *work2, *in0, *in1, *in2, *in3;

  assert(data);
  assert(ix>0);
  assert(iy>0);
  assert(work);
  assert(HQF);
  assert(GQF);
  
  oy = 2*iy;
  n = iy*ix;
  
  work1 = work;
  work2 = work1 + 2*n;
  in0 = data;
  in1 = in0 + n;
  in2 = in1 + n;
  in3 = in2 + n;

  for( i=0; i<ix; i++ )
    {
      acdpe ( work1+i, ix, in0, iy, HQF );
      in0 += iy;
      acdpo ( work1+i, ix, in1, iy, GQF );
      in1 += iy;

      acdpe ( work2+i, ix, in2, iy, HQF );
      in2 += iy;
      acdpo ( work2+i, ix, in3, iy, GQF );
      in3 += iy;
    }

  for(i=0; i<oy; i++)
    {
      acdpe ( data+i, oy, work1, ix, HQF );
      work1 += ix;
      acdpo ( data+i, oy, work2, ix, GQF );
      work2 += ix;
    }

  return;
}



/**********************************************************************
 * scdpd2()
 *
 *  [S]eparable [C]onvolution-[D]ecimation;
 *  [P]eriodic, [D]isjoint, [2]-dimensional.
 *
 *  This function performs a separable 2-dimensional periodic convolution-
 *  decimation, using quadrature mirror filter coefficients defined
 *  by some passed parameters.  It superposes the results onto 4 preallocated
 *  arrays.  It assumes that the output and input arrays are disjoint.
 *
 *  Note: if any of the 4 output array pointers are NULL, then they will
 *  not be written.  This provides a mechanism to get a single descendent.
 *
 *  Calling sequence and basic algorithm:
 *
 *    scdpd2( OUT0,OUT1,OUT2,OUT3, IN, IX,IY, WORK, HQF,GQF ):
 *      Let OY = IY/2
 *      For I=0 to IX-1
 *         cdpe( WORK+I, IX, IN+I*IY, IY, HQF )
 *      For I=0 to OY-1
 *         cdpo( OUT0+I, OY, WORK+I*IX, IX, HQF )
 *         cdpo( OUT2+I, OY, WORK+I*IX, IX, GQF )
 *      For I=0 to IX-1
 *         cdpe( WORK+I, IX, IN+I*IY, IY, GQF )
 *      For I=0 to OY-1
 *         cdpo( OUT1+I, OY, WORK+I*IX, IX, HQF )
 *         cdpo( OUT3+I, OY, WORK+I*IX, IX, GQF )
 *
 *  Inputs:
 *	(real *)out0		These arrays must be preallocated with at 
 *	(real *)out1		  least `(iy/2)*(ix/2)' elements each.
 *	(real *)out2		  They will be overwritten, even before all
 *	(real *)out3		  the input data has been read.
 *
 *	(const real *)in	This array must be preallocated and defined
 *				  in at least `iy*ix' locations.
 *				  It is not changed by this function.
 *
 *	(int)ix			These positive integers are the number of rows
 *	(int)iy			  and columns in the array `in[]'
 *
 *	(real *)work		This temporary array must be preallocated
 *				  with at least `iy*ix/2' elements.  
 *				  It will be overwritten.
 *
 *	(const pqf *)HQF	This must be a predefined low-pass QF.
 *	(const pqf *)GQF	This must be a predefined high-pass QF.
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
 *	cdpe(), cdpo()
 *
 *  Assumptions:
 *	1. in != NULL
 *	2. out0!=NULL || out1!=NULL || out2!=NULL || out3!=NULL
 *	3. work != NULL
 *	4. `in[0,...,ix*iy-1]' is disjoint from each of 
 *	   `out0[0,...,ix*iy/4-1]',...,`out3[0,...,ix*iy/4-1]'.
 *	5. HQF != NULL
 *	6. GQF != NULL
 *	7. `ix' is even and positive
 *	8. `iy' is even and positive
 */
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
	 const pqf *GQF)	/* High-pass quadrature filter data struct  */
{
  int i, oy;
  const real *iptr;
  real *wptr;

  assert(in);
  assert(work);
  assert(out0||out1||out2||out3);
  assert(HQF);
  assert(GQF);
  assert(ix>0);  assert( !(ix&1) );
  assert(iy>0);  assert( !(iy&1) );

  oy = iy/2;

  /* First step (`out0[]' and `out2[]'): */
  if( out0 || out2 )
    {
      /* Y-direction 1: Low-pass filter the rows of `in[]': */
      iptr = in;
      for( i=0; i<ix; i++ )
	{
	  cdpe ( work+i, ix, iptr, iy, HQF );
	  iptr += iy;
	}
      
      /* X-direction 1: Loop over the rows of `work[]': */
      wptr = work;
      for( i=0; i<oy; i++ )
	{
	  if( out0 ) cdpo ( out0+i, oy, wptr, ix, HQF );
	  if( out2 ) cdpo ( out2+i, oy, wptr, ix, GQF );
	  wptr += ix;
	}
    }

  /*  Second step (`out1[]' and `out3[]'): */
  if( out1 || out3 )
    {
      /* Y-direction 2: High-pass filter the rows of `in[]': */
      iptr = in;
      for( i=0; i<ix; i++ )
	{
	  cdpe ( work+i, ix, iptr, iy, GQF );
	  iptr += iy;
	}

      /* X-direction 2: Loop over the rows of `work[]': */
      wptr = work;
      for( i=0; i<oy; i++ )
	{
	  if( out1 ) cdpo ( out1+i, oy, wptr, ix, HQF );
	  if( out3 ) cdpo ( out3+i, oy, wptr, ix, GQF );
	  wptr += ix;
	}
    }
  return;
}

/**********************************************************************
 * scdpe2()
 *
 *  [S]eparable [C]onvolution-[D]ecimation;
 *  [P]eriodic, [E]quals, [2]-dimensional.
 *
 *  This function performs a separable 2-dimensional periodic convolution-
 *  decimation, using quadrature mirror filter coefficients defined
 *  by some passed parameters.  It assigns the results into 4 preallocated
 *  arrays.  It assumes that the output and input arrays are disjoint.
 *
 *  Note: if any of the 4 output array pointers are NULL, then they will
 *  not be written.  This provides a mechanism to get a single descendent.
 *
 *  Calling sequence and basic algorithm:
 *
 *    scdpe2( OUT0,OUT1,OUT2,OUT3, IN, IX,IY, WORK, HQF,GQF ):
 *      Let OY = IY/2
 *      For I=0 to IX-1
 *         cdpe( WORK+I, IX, IN+I*IY, IY, HQF )
 *      For I=0 to OY-1
 *         cdpe( OUT0+I, OY, WORK+I*IX, IX, HQF )
 *         cdpe( OUT2+I, OY, WORK+I*IX, IX, GQF )
 *      For I=0 to IX-1
 *         cdpe( WORK+I, IX, IN+I*IY, IY, GQF )
 *      For I=0 to OY-1
 *         cdpe( OUT1+I, OY, WORK+I*IX, IX, HQF )
 *         cdpe( OUT3+I, OY, WORK+I*IX, IX, GQF )
 *
 *  Inputs:
 *	(real *)out0		These arrays must be preallocated with at 
 *	(real *)out1		  least `(iy/2)*(ix/2)' elements each.
 *	(real *)out2		  They will be overwritten, even before all
 *	(real *)out3		  the input data has been read.
 *
 *	(const real *)in	This array must be preallocated and defined
 *				  in at least `iy*ix' locations.
 *				  It is not changed by this function.
 *
 *	(int)ix			These positive integers are the number of rows
 *	(int)iy			  and columns in the array `in[]'
 *
 *	(real *)work		This temporary array must be preallocated
 *				  with at least `iy*ix/2' elements.  
 *				  It will be overwritten.
 *
 *	(const pqf *)HQF	This must be a predefined low-pass QF.
 *	(const pqf *)GQF	This must be a predefined high-pass QF.
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
 *	cdpe(), assert()
 *
 *  Assumptions:
 *	1. in != NULL
 *	2. out0!=NULL || out1!=NULL || out2!=NULL || out3!=NULL
 *	3. work != NULL
 *	4. `in[0,...,ix*iy-1]' is disjoint from each of 
 *	   `out0[0,...,ix*iy/4-1]',...,`out3[0,...,ix*iy/4-1]'.
 *	5. HQF != NULL
 *	6. GQF != NULL
 *	7. `ix' is even and positive
 *	8. `iy' is even and positive
 */
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
	 const pqf *GQF)	/* High-pass quadrature filter data struct  */
{
  int i, oy;
  const real *iptr;
  real *wptr;

  assert(in);
  assert(work);
  assert(out0||out1||out2||out3);
  assert(HQF);
  assert(GQF);
  assert(ix>0);  assert( !(ix&1) );
  assert(iy>0);  assert( !(iy&1) );

  oy = iy/2;

  /* First step (`out0[]' and `out2[]'): */
  if( out0 || out2 )
    {
      /* Y-direction 1: Low-pass filter the rows of `in[]': */
      iptr = in;
      for( i=0; i<ix; i++ )
	{
	  cdpe ( work+i, ix, iptr, iy, HQF );
	  iptr += iy;
	}
      
      /* X-direction 1: Loop over the rows of `work[]': */
      wptr = work;
      for( i=0; i<oy; i++ )
	{
	  if( out0 ) cdpe ( out0+i, oy, wptr, ix, HQF );
	  if( out2 ) cdpe ( out2+i, oy, wptr, ix, GQF );
	  wptr += ix;
	}
    }

  /*  Second step (`out1[]' and `out3[]'): */
  if( out1 || out3 )
    {
      /* Y-direction 2: High-pass filter the rows of `in[]': */
      iptr = in;
      for( i=0; i<ix; i++ )
	{
	  cdpe ( work+i, ix, iptr, iy, GQF );
	  iptr += iy;
	}

      /* X-direction 2: Loop over the rows of `work[]': */
      wptr = work;
      for( i=0; i<oy; i++ )
	{
	  if( out1 ) cdpe ( out1+i, oy, wptr, ix, HQF );
	  if( out3 ) cdpe ( out3+i, oy, wptr, ix, GQF );
	  wptr += ix;
	}
    }
  return;
}

/**********************************************************************
 * scdpi2()
 *
 *  [S]eparable [C]onvolution-[D]ecimation;
 *  [P]eriodic, [I]n-place, [2]-dimensional.
 *
 *  This function performs a separable 2-dimensional periodic convolution-
 *  decimation, using quadrature mirror filter coefficients defined
 *  by some passed parameters.  It assigns the results back into the
 *  input array, which is completely read before any output is written.
 *
 *  Calling sequence basic algorithm:
 *
 *   scdpi2( DATA, IX, IY, WORK, HQF,GQF ):
 *      Let OY = IY/2
 *      Let N  = OY*IX/2
 *      Let WORK1 = WORK
 *      Let WORK2 = WORK + OY*IX
 *      For I=0 to IX-1
 *         cdpe( WORK1+I, IX, DATA+I*IY, IY, HQF )
 *         cdpe( WORK2+I, IX, DATA+I*IY, IY, GQF )
 *      For I=0 to OY-1
 *         cdpe( DATA+ I,    OY, WORK1+I*IX, IX, HQF )
 *         cdpe( DATA+2*N+I, OY, WORK1+I*IX, IX, GQF )
 *         cdpe( DATA+ N +I, OY, WORK2+I*IX, IX, HQF )
 *         cdpe( DATA+3*N+I, OY, WORK2+I*IX, IX, GQF )
 *
 *  Inputs:
 *	(real *)data		This array must be preallocated with at 
 *				  least `iy*ix' elements.  It will
 *				  be overwritten.
 *
 *	(int)ix			These even positive integers are the number
 *	(int)iy			  of rows and columns in the array `data[]'
 *
 *	(real *)work		This temporary array must be preallocated
 *				  with at least `iy*ix' elements.  
 *				  It will be overwritten.
 *
 *	(const pqf *)HQF	This must be a predefined low-pass QF.
 *	(const pqf *)GQF	This must be a predefined high-pass QF.
 *
 *  Outputs:
 *	(real *)data		This arrays is replaced by side effect.
 *
 *	(real *)work		This array is trashed by side effect.
 *
 *  Assumptions:
 *	1. `iy' and `ix' are positive and divisible by 2
 *	2. `data' is non-NULL
 *	3. `work' is non-NULL
 *	4. HQF != NULL
 *	5. GQF != NULL
 *
 *  External functions called:
 *	cdpo(), cdpe(), assert()
 */
extern void
  scdpi2(
	 real *data,		/* Pointer to joint input and output array  */
	 int  ix,		/* Number of rows in the array `data[]'     */
	 int  iy,		/* Number of columns in the array `data[]'  */
	 real *work,		/* Scratch array of size `ix*iy'            */
	 const pqf *HQF,	/* Low-pass quadrature filter data struct   */
	 const pqf *GQF)	/* High-pass quadrature filter data struct  */
{
  int oy, n, i;
  real *work1, *work2;

  assert(iy>0);
  assert(ix>0);
  assert( !(iy&1) );
  assert( !(ix&1) );
  assert( data );
  assert( work );
  assert(HQF); assert(GQF);

  oy = iy/2;
  n  = oy*ix/2;

  work1 = work;
  work2 = work + oy*ix;
  for( i=0; i<ix; i++ )
    {
      cdpe( work1+i, ix, data+i*iy, iy, HQF );
      cdpe( work2+i, ix, data+i*iy, iy, GQF );
    }
  for( i=0; i<oy; i++ )
    {
      cdpe( data+ i,    oy, work1+i*ix, ix, HQF );
      cdpe( data+2*n+i, oy, work1+i*ix, ix, GQF );
      cdpe( data+ n +i, oy, work2+i*ix, ix, HQF );
      cdpe( data+3*n+i, oy, work2+i*ix, ix, GQF );
    }
  return;
}
