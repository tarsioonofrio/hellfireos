/*
 * 
 * Functions to write PostScript commands to a stream.
 * 
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <stdio.h>
#include "real.h"
#include "epsfns.h"

/*********************************************************************
 * epsepilogue()
 * 
 *  This function writes a showpage command to the end of a specified file,
 *  which is intended to contain an Encapsulated PostScript program.
 * 
 * Calling sequence:
 *	epsepilogue(psfile)
 * 
 * Input:
 * 	(FILE *)psfile		This struct must represent a file opened
 *				  for writing.
 * Output: <none>
 * 
 * External Functions Called:
 *   int  fprintf(FILE *file, const char *format, ...)    From <stdio.h>
 * 
 * Explicit Side Effects:
 *   Writes 3 constant strings to `psfile.'
 */
extern void
  epsepilogue(
	      FILE  *psfile)	/* File to write into. */
{
  fprintf( psfile, "%% end of EPS commands\n");
  fprintf( psfile, "showpage\n");
  fprintf( psfile, "%% all done\n");
  return;
}  



/*********************************************************************
 * epsfrect()
 * 
 * Print Encapsulated PostScript commands to draw and fill a
 * rectangle (to specified gray level)into the specified file.
 *
 * Calling sequence:
 *	epsfrect(psfile, xmin, ymin, xmax, ymax, gray)
 *
 * Input:
 *      (FILE *)psfile		This file must be opened for writing.
 * 
 *	(real)xmin		(xmin, ymin) is the lower left corner, and
 * 	(real)ymin		  (xmax, ymax) is the upper right corner of
 * 	(real)xmax		  the rectangle to be drawn.  All coordinates
 * 	(real)ymax		  must be in the interval [0.0, 1.0].
 *
 *	(real)gray		Gray shading of the rectangle: 0.0 is black,
 *				  1.0 is white.
 *
 * Output: <none>
 * 
 * External Functions Called:
 *   int  fprintf(FILE *file, const char *format, ...)    From <stdio.h>
 * 
 * Explicit Side Effects:
 *   Writes 7 strings to `psfile.'
 * 
 */
extern void
  epsfrect(
	   FILE  *psfile,	/* File to write into. */
	   real xmin,		/* Least x-coord (in [0.0, 1.0]). */
	   real ymin,		/* Least y-coord (in [0.0, 1.0]). */
	   real xmax,		/* Greatest x-coord (in [0.0, 1.0]). */
	   real ymax,		/* Greatest y-coord (in [0.0, 1.0]). */
	   real gray)		/* Gray level: 0.0=black, 1.0=white. */
{
  fprintf( psfile, "%% begin new rectangle:\n");
  fprintf( psfile, "%f xunits %f yunits moveto\n", xmin, ymin);
  fprintf( psfile, "%f xunits %f yunits lineto\n", xmin, ymax);
  fprintf( psfile, "%f xunits %f yunits lineto\n", xmax, ymax);
  fprintf( psfile, "%f xunits %f yunits lineto\n", xmax, ymin);
  fprintf( psfile, "closepath\n");
  fprintf( psfile, "%f setgray fill\n", gray);
  return;
}

/*********************************************************************
 * epslineto()
 *
 * This function prints the PostScript command `lineto' into a 
 * specified file, with arguments followed by predefined macros `yunits'
 * and `xunits,' which are defined by `epsprologue()'.
 *
 * Calling sequence:
 *	epslineto(psfile, xval, yval)
 *
 * Input:
 *	(FILE *)psfile		This file must be opened for writing.
 *
 *	(real)xval		These arguments to lineto must be in the
 *	(real)yval		  interval [0.0, 1.0]
 *
 * Output: <none>
 *
 * External Functions Called:
 *   int  fprintf(FILE *file, const char *format, ...)    From <stdio.h>
 *
 * Explicit Side Effects:
 *   Writes 1 string to `psfile.' 
 */
extern void
  epslineto(
	    FILE  *psfile,	/* File into which we will write. */
	    real xval,		/* x-coord (in [0.0, 1.0]).       */
	    real yval)		/* y-coord (in [0.0, 1.0]).       */
{
  fprintf( psfile, "%f xunits %f yunits lineto\n", xval, yval);
  return;
}

/*********************************************************************
 * epsmoveto()
 *
 * This function prints the PostScript command `moveto' into a 
 * specified file, with arguments followed by predefined macros `yunits'
 * and `xunits,' which are defined by `epsprologue()'.
 *
 * Calling sequence:
 *	epsmoveto(psfile, xval, yval)
 *
 * Input:
 *	(FILE *)psfile		This file must be opened for writing.
 *
 *	(real)xval		These arguments to lineto must be in the
 *	(real)yval		  interval [0.0, 1.0]
 *
 * Output: <none>
 *
 * External Functions Called:
 *   int  fprintf(FILE *file, const char *format, ...)    From <stdio.h>
 *
 * Explicit Side Effects:
 *   Writes 1 string to `psfile.' 
 */
extern void
  epsmoveto(
	    FILE  *psfile,	/* File into which we will write. */
	    real xval,		/* x-coord (in [0.0, 1.0]).       */
	    real yval)		/* y-coord (in [0.0, 1.0]).       */
{
  fprintf( psfile, "%f xunits %f yunits moveto\n", xval, yval);
  return;
}

/*********************************************************************
 * epsprologue()
 *
 * This function writes an Encapsulated PostScript prologue to a specified
 * file.  The prologue defines a BoundingBox, draws an outline just inside 
 * it, then sets two scaling factors `yunits' and `xunits' which transforms
 * coordinates in [0,1]x[0,1] to locations within the BoundingBox.
 *
 * Calling sequence:
 *	epsprologue(psfile, bbxmin, bbymin, bbxmax, bbymax)
 *
 * Input:
 *	(FILE *)psfile		This file must be opened for writing.
 * 
 *	(int)bbxmin		(bbxmin, bbymin) is the lower left corner,
 *	(int)bbymin		  and (bbxmax, bbymax) is the upper right
 *	(int)bbxmax		  corner of the BoundingBox, expressed in
 *	(int)bbymax		  standard PostScript coordinates: 72 points
 *				  to the inch, origin at the lower left side
 *				  of the page.
 * 
 * Output:  <none>
 * 
 * External funtions called:
 *    int fprintf(FILE *file, const char *format, ...)   From <stdio.h>
 * 
 * Explicit side effects:
 *    Prints EPS prologue lines, specifying a BoundingBox, to `psfile.'
 */
extern void
  epsprologue(
	      FILE *psfile,	/* Print to this file.  */
	      int   bbxmin,	/* X-coord of lower left BoundingBox corner. */
	      int   bbymin,	/* Y-coord of lower left BoundingBox corner. */
	      int   bbxmax,	/* X-coord of upper left BoundingBox corner. */
	      int   bbymax)	/* Y-coord of upper left BoundingBox corner. */
{
  fprintf( psfile, "%%!PS-Adobe-2.0 EPSF-1.2\n");
  fprintf( psfile, "%%%%Title: \n");
  fprintf( psfile, "%%%%Creator: \n");
  fprintf( psfile, "%%%%CreationDate: \n");
  fprintf( psfile, "%%%%For: \n");
  fprintf( psfile, "%%%%DocumentFonts: (atend)\n");
  fprintf( psfile, "%%%%Pages: 0 1\n");
  fprintf( psfile, "%%%%BoundingBox: %d %d %d %d\n",
	  bbxmin, bbymin, bbxmax, bbymax);
  fprintf( psfile, "%%%%EndComments\n");

  fprintf( psfile, "%%%%BeginSetup\n");
  fprintf( psfile, "/xunits {\n");
  fprintf( psfile, " %d mul\n", bbxmax-bbxmin-4);
  fprintf( psfile, "} def\n");
  fprintf( psfile, "/yunits {\n");
  fprintf( psfile, " %d mul\n", bbymax-bbymin-4);
  fprintf( psfile, "} def\n");
  fprintf( psfile, "%% Outline the rectangle\n");
  fprintf( psfile, "%d %d moveto\n", bbxmin+1, bbymin+1);
  fprintf( psfile, "%d %d lineto\n", bbxmin+1, bbymax-1);
  fprintf( psfile, "%d %d lineto\n", bbxmax-1, bbymax-1);
  fprintf( psfile, "%d %d lineto\n", bbxmax-1, bbymin+1);
  fprintf( psfile, "closepath stroke\n");
  fprintf( psfile, "%d %d translate\n", bbxmin+2, bbymin+2);
  fprintf( psfile, "%%%%EndSetup\n");
  return;
}

/*********************************************************************
 * epsstroke()
 *
 * Print the Encapsulated PostScript command `stroke' into a specified file.
 *
 * Calling sequence:
 *	epsstroke(psfile)
 *
 * Input:
 *	(FILE *)psfile		This file must be opened for writing.
 *
 * Output:  <none>
 * 
 * External Functions Called:
 *   int  fprintf(FILE *file, const char *format, ...)    From <stdio.h>
 *
 * Explicit Side Effects:
 *   Writes 1 string to `psfile.'
 * 
 */
extern void
  epsstroke(
	    FILE  *psfile)	/* File into which we will write. */
{
  fprintf( psfile, "stroke\n");
  return;
}

/*********************************************************************
 * epstranslate()
 *
 * Write the Encapsulated PostScript command `translate' to a file.
 *
 * Calling sequence:
 *	epstranslate(psfile, xptval, yptval)
 *
 * Input:
 *	(FILE *)psfile		This file must be opened for writing.
 * 
 *	(int)xptval		The origin, after translation, will be at
 *	(int)yptval		  current coordinates (xptval, yptval).
 *
 * Output:  <none>
 * 
 * External funtions called:
 *    int fprintf(FILE *file, const char *format, ...)   From <stdio.h>
 * 
 * Explicit side effects:
 *    Prints one line to `psfile.'
 */
extern void
  epstranslate(
	       FILE *psfile,	/* Print to this file. */
	       int   xptval,	/* X-coord of translated origin. */
	       int   yptval)	/* Y-coord of translated origin. */
{
  fprintf( psfile, "%d %d translate\n", xptval, yptval);
  return;
}
