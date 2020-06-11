/*
 * 
 * These functions cause PostScript commands to be written which 
 * plot a 1-D signal and its nominal time-frequency-plane image.
 * 
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <stdio.h>
#include "real.h"
#include "common.h"
#include "tfa.h"
#include "epsfns.h"
#include "hardcopy.h"

/***********************************************************************
 * plotsig()
 * 
 * Write instructions in PostScript to a specified file.  The instructions
 * will plot the graph of { (i,signal[i]) : 0<= i < length }, with lines
 * between the points, near the bottom of an 8.5x11" page.
 * 
 * Calling Sequence and Basic Algorithm:
 *
 *  plotsig( PSFILE, SIGNAL, LENGTH ):
 *     epsprologue( PSFILE,  LLXS, LLYS, URXS, URYS )
 *     epstranslate( PSFILE, 0.0, (URYS - LLYS)/2.0 )
 *     Let NORM = 0.0
 *     For N = 0 to LENGTH-1
 *        Let NORM = max( absval(SIGNAL[N]), NORM )
 *     If NORM == 0.0
 *         epsmoveto( PSFILE, 0.0, 0.0 )
 *         epslineto( PSFILE, 1.0, 0.0 )
 *     Else
 *         Let NORM = 0.45 / NORM
 *         Let YVAL = SIGNAL[0]*NORM
 *         epsmoveto( PSFILE, 0.0, YVAL )
 *         If LENGTH == 1
 *            epslineto( PSFILE, 1.0, YVAL )
 *         Else
 *            Let XVAL = 0.0
 *            Let INCR = 1.0 / (LENGTH - 1.0)
 *            For N = 1 to LENGTH-1
 *               XVAL += INCR
 *               Let YVAL = SIGNAL[N]*NORM
 *               epslineto( PSFILE, XVAL, YVAL )
 *         epsstroke( PSFILE )
 *     epsepilogue( PSFILE )
 *
 *
 * Input:
 *	(FILE *)psfile		This is the output file stream, which must
 *				  be opened and writable.
 *
 *	(real *)signal		This array holds the Y-values of a uniformly-
 *				  sampled 1-dimensional signal which will be
 *				  converted to PostScript plotting commands.
 *
 *	(int)length		Plot this many elements of `signal[]'.
 *
 * Output:
 *	The function causes text to be written to an output stream.
 *
 * External Functions Called:
 *   epsprologue()
 *   epsepilogue()
 *   epsmoveto()
 *   epslineto()
 *   epsstroke()
 *   epstranslate()
 */
extern void
  plotsig(
	  FILE *psfile,		/* Write into this file stream. */
	  const real  *signal,	/* Array of `length' Y-values.  */
	  int length)		/* Positive integer. */
{
  real norm, val, xval, yval, incr;
  int n;
  
  /* Write encapsulated PostScript prologue.  */
  epsprologue(psfile, LLXS, LLYS, URXS, URYS);

  /* Translate the origin to left middleheight within the BoundingBox.  */
  epstranslate( psfile, 0, (int)(( URYS - LLYS )/2) );

  /* Find the maximum amplitude for normalization. */
  norm = 0;
  for( n=0; n<length; n++ )
    {
      val = absval(signal[n]);
      norm = max( norm, val );
    }

  if( norm == 0.0 )		/* Signal is 0: plot a straight 0 line  */
    {
      epsmoveto( psfile, 0.0, 0.0 );
      epslineto( psfile, 1.0, 0.0);
    }
  else				/* There exist nonzero points */
    {
      norm = 0.45 / norm;
      yval = (*signal++)*norm;

      /* Start plot at first signal[] value. */
      epsmoveto( psfile, 0.0, yval );

      /* One-point signals are a special case: */
      if( length == 1 )
	{
	  epslineto( psfile, 1.0, yval );
	}
      else			/* Multiple signal points */
	{
	  xval = 0.0;
	  incr = 1.0 /(length - 1);
	  while ( --length > 0 )
	    {
	      xval += incr;
	      yval  = (*signal++)*norm;
	      epslineto( psfile, xval, yval );
	    }
	}

      /* Ink the plot drawn by epslineto() calls. */
      epsstroke( psfile );
    }
  epsepilogue( psfile );	/* EPS epilogue; draws page.  */

  return;
}

/***********************************************************************
 * tfa1s2ps()
 * 
 * Produce PostScript commands to draw a nominal time-frequency density 
 * plot (information cells in the time-frequency plane) from an array of
 * TFA1s.  Write the results to the specified file, with a prologue
 * that centers the density plot in the top part of an 8.5x11" page.
 *
 * Calling sequence:
 *
 *  tfa1s2ps( PSFILE, SAMPLES, ATOMS, NUM ):
 *     epsprologue( PSFILE, LLXA, LLYA, URXA, URYA )
 *     Let ANORM = 0.0
 *     For N = 0 to NUM-1
 *        Let ANORM = max( absval(ATOMS[N].AMPLITUDE), ANORM )
 *     If ANORM > 0.0
 *        Let ANORM  = 1.0 / ANORM
 *        For N = 0 to NUM-1
 *           Let AMPL  = ATOMS[N].AMPLITUDE
 *           Let GRAY  = ANORM * AMPL
 *           If GRAY > MINGRAY then
 *              Let WIDTH = ( 1<<ATOMS[N].LEVEL ) / SAMPLES
 *              Let XMIN = WIDTH * ATOMS[N].OFFSET
 *              Let XMAX = XMIN + WIDTH
 *              Let HEIGHT = 1.0 / WIDTH
 *              Let YMIN = HEIGHT * ATOMS[N].BLOCK
 *              Let YMAX = YMIN + HEIGHT
 *              epsfrect( PSFILE, XMIN, YMIN, XMAX, YMAX, 1.0-GRAY )
 *     epsepilogue( PSFILE )
 *
 * Input:
 * 	(FILE *)psfile		This file stream must be open for writing.
 * 
 *	(int)samples		The original signal has this many samples.
 * 
 *	(tfa1 *)atoms		This must be preallocated and assigned.
 *
 *	(int)num		Use this many elements from `atoms[]'.
 *
 *
 * External functions called:
 *   epsprologue(psfile, bbxmin, bbymin, bbxmax, bbymax)
 *   epsepilogue(psfile)
 *   epsfrect(psfile, xmin, ymin, xmax, ymax, gray)
 *   maxabsq(head)
 * 
 * Explicit side effects:
 *   Writes PostScript commands, depicting a time-frequency decomposition
 *   given by `atoms[]', to the file `psfile'
 */
extern void
  tfa1s2ps(
	   FILE *psfile,	/* Name of the output file to use.  */
	   int samples,		/* Number of samples in the signal. */
	   tfa1 *atoms,		/* Array of TFA1 structs.           */
	   int num)		/* Number of TFA1 structs.          */
{
  real ymin, ymax, xmin, xmax, gray, anorm, ampl, width, height;
  int n;

  /* The EPS prologue specifies the BoundingBox.  These coordinates
   * leave about a 1" margin on an 8.5x11" sheet of paper.
   */
  epsprologue(psfile, LLXA, LLYA, URXA, URYA );

  /* Find the maximum square of all the amplitudes in `atoms[]'. */
  anorm = 0;
  for( n=0; n<num; n++ )
    {
      ampl = atoms[n].amplitude;
      ampl *= ampl;
      anorm = max( ampl, anorm);
    }

  /* Plot only if there are some nonzero amplitudes. */
  if( anorm > 0.0 )
    {
      anorm  = 1.0 / anorm;
      /* Loop through the tf_atom's in the array. */
      for( n=0; n<num; n++ )
	{
	  ampl  = atoms[n].amplitude;
	  gray  = anorm * ampl * ampl;
	  /* Only draw information cells which are dark enough to plot. */
	  if( gray > MINGRAY )
	    {
	      height = 1.0 /  ( 1 << atoms[n].level );
	      ymin = height * atoms[n].block ;
	      ymax = ymin + height;

	      width  = 1.0 / ( samples >> atoms[n].level );
	      xmin  = width * atoms[n].offset;
	      xmax  = xmin + width;

	      epsfrect( psfile, xmin, ymin, xmax, ymax, 1.0-gray );
	    }
	}
    }				/* end "if(anorm==0)" */
  epsepilogue( psfile );	/* Write an epilogue to the output file. */
  return;
}
