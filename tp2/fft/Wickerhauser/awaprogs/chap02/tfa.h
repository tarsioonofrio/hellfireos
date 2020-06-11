/*
 * This header file declares the  `tfa1', `tfa2' and  `tfad' data types.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef TFA_HDR_ALREADY_INCLUDED
#define TFA_HDR_ALREADY_INCLUDED

#include "real.h"

typedef struct {
  real amplitude;		/* Coefficient amplitude. */
  int  level;			/* 0 is the original signal. */
  int  block;			/* Read from left to right. */
  int  offset;			/* Position within the block. */
}
tfa1;				/* 1-dimensional time-frequency atom. */

typedef struct {
  real amplitude;		/* Coefficient amplitude. */
  int  xlevel;			/* 0 is the original signal. */
  int  ylevel;			/* 0 is the original signal. */
  int  xblock;			/* Read from left to right. */
  int  yblock;			/* Read from top to bottom. */
  int  xoffset;			/* Position within the X block. */
  int  yoffset;			/* Position within the Y block. */
}
tfa2;				/* 2-dimensional time-frequency atom. */

typedef struct {
  real amplitude;		/* Coefficient amplitude. */
  int  dimension;		/* The number of signal dimensions D. */
  int  levels;			/* Coded form of the D level indices */
  int  blocks;			/* Coded form of the D block indices. */
  int  offsets;			/* Coded form of the D position indices. */
}
tfad;				/* D-dimensional time-frequency atom. */


#endif /* TFA_HDR_ALREADY_INCLUDED */

