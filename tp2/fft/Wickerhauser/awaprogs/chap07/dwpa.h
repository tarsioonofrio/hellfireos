/*
 * Declare functions to perform 1-D discrete wavelet packet analysis.
 * 
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 * 
 */

#ifndef DWPA_HDR_ALREADY_INCLUDED
# define DWPA_HDR_ALREADY_INCLUDED

#include "real.h"
#include "btn.h"
#include "hedge.h"
#include "qf.h"
#include "tfa.h"

extern void
  dwpap2abt0(
	     real *data,	/* Input/output array binary tree. */
	     int n,		/* The row length of `data[]'.     */
	     int maxlevel,	/* `data[]' has `1+maxlevel' rows. */
	     const pqf *HQF,	/* Low-pass QF data structure.     */
	     const pqf *GQF);	/* High-pass QF data structure.    */

extern real *
  dwpap2abt(
	     const real *in,	/* Input array, not changed.       */
	     int n,		/* The length of the array `in[]'. */
	     int maxlevel,	/* Expand to depth `maxlevel'.     */
	     const pqf *HQF,	/* Low-pass QF data structure.     */
	     const pqf *GQF);	/* High-pass QF data structure.    */

extern int
  dwpap2hedger(
	       hedge *graph,	/* Partially specified data struct. */
	       int j,		/* The current block index.         */
	       int n,		/* The current block length.        */
	       int s,		/* The current scale index.         */
	       const pqf *HQF,	/* Low-pass QF data structure.      */
	       const pqf *GQF,	/* High-pass QF data structure.     */
	       real *work);	/* Intermediate wavelet amplitudes. */

extern hedge *
  dwpap2hedge(
	      real *in,		/* Joint input and output array.    */
	      int length,	/* Number of elements of `in[]'.    */
	      unsigned char *levels, /* Graph basis levels list.    */
	      int blocks,	/* Number of blocks in the basis.   */
	      const pqf *HQF,	/* Low-pass QF data structure.      */
	      const pqf *GQF);	/* High-pass QF data structure.     */

extern real *
  tfa1s2dwpsp(
	      tfa1 *atoms,	/* Input array of time-frequency atoms.    */
	      int num,		/* Number of elements of `atoms[]'.        */
	      int n,		/* Period of the signal to be synthesized. */
	      const pqf *HQF,	/* Periodized low-pass QF data structure.  */
	      const pqf *GQF);	/* Periodized high-pass QF data structure. */

extern int
  hedge2dwpspr(
	       hedge *graph,	/* Input data structure.            */
	       int j,		/* The current block index.         */
	       int n,		/* The current block length.        */
	       int s,		/* The current scale index.         */
	       const pqf *HQF,	/* Low-pass QF data structure.      */
	       const pqf *GQF,	/* High-pass QF data structure.     */
	       real *work);	/* Intermediate wavelet amplitudes. */

extern real *
  hedge2dwpsp(
	      hedge *graph,	/* Input data structure.                 */
	      int length,	/* Period of the synthesized signal.     */
	      const pqf *HQF,	/* Low-pass periodic QF data structure.  */
	      const pqf *GQF);	/* High-pass periodic QF data structure. */

extern real *
  hedge2dwpspabt(
		 hedge *graph,	/* Basis amplitudes to reconstruct.    */
		 int n,		/* Period of the reconstructed signal. */
		 const pqf *HQF, /* Low-pass QF data structure.        */
		 const pqf *GQF); /* High-pass QF data structure.       */

extern interval *
  cdachild(
	   interval *child,	/* Output INTERVAL data structure. */
	   interval *parent,	/* Input INTERVAL data structure.  */
	   const pqf *F);	/* Aperiodic QF data structure.    */

extern void
  dwpaa2btntr(
	      btn *node,	/* Root of the current BTN tree. */
	      int level,	/* Depth of the tree to develop. */
	      const pqf *HQF,	/* Low-pass QF data structure.   */
	      const pqf *GQF);	/* High-pass QF data structure.  */

extern btn *
  dwpaa2btnt(
	     real *in,		/* Input signal array. */
	     int least,		/* Least valid index of `in[]'. */
	     int final,		/* Greatest valid index of `in[]'. */
	     int maxlevel,	/* Depth of the tree to develop. */
	     const pqf *HQF,	/* Low-pass QF data structure.   */
	     const pqf *GQF);	/* High-pass QF data structure.  */

extern int
  dwpaa2hedger(
	       hedge *graph,	/* Level-specified graph basis.   */
	       btn *node,	/* Root of the current BTN tree.  */
	       int j,		/* Current block number.          */
	       int s,		/* Current level in the BTN tree. */
	       const pqf *HQF,	/* Low-pass QF data structure.    */
	       const pqf *GQF);	/* High-pass QF data structure.   */

extern hedge *
  dwpaa2hedge(
	      unsigned char *levels, /* Levels specifying the basis. */
	      int blocks,	/* Number of blocks in the basis.    */
	      real *data,	/* Input signal array.               */
	      int least,	/* Least valid index of `data[]'.    */
	      int final,	/* Greatest valid index of `data[]'. */
	      const pqf *HQF,	/* Low-pass QF data structure.       */
	      const pqf *GQF);	/* High-pass QF data structure.      */

extern interval *
  acdaparent(
	     interval *parent,	/* Output INTERVAL data structure. */
	     interval *child,	/* Input INTERVAL data structure.  */
	     const pqf *F);	/* Aperiodic QF data structure.    */

extern void
  btnt2dwpsa(
	     btn *root,		/* Root of the tree to reconstruct. */
	     const pqf *HQF,	/* Low-pass QF data structure.      */
	     const pqf *GQF);	/* High-pass QF data structure.     */

extern void
  btnt2dwpsa0(
	      btn *root,	/* Root of the tree to reconstruct. */
	      const pqf *HQF,	/* Low-pass QF data structure.      */
	      const pqf *GQF);	/* High-pass QF data structure.     */

extern interval *
  tfa1s2dwpsa(
	      tfa1 *atoms,	/* List of wavelet packets to reconstruct. */
	      int num,		/* Number of wavelet packets in `atoms[]'. */
	      const pqf *HQF,	/* Low-pass QF data structure.      */
	      const pqf *GQF);	/* High-pass QF data structure.     */

extern interval *
  hedge2dwpsa(
	      hedge *graph,	/* Basis amplitudes to reconstruct. */
	      const pqf *HQF,	/* Low-pass QF data structure.      */
	      const pqf *GQF);	/* High-pass QF data structure.     */

#endif    /* DWPA_HDR_ALREADY_INCLUDED */
