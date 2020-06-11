/*
 * Type definitions and functions for complex arithmetic in Standard C.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef COMPLEX_HDR_ALREADY_INCLUDED
#define COMPLEX_HDR_ALREADY_INCLUDED

#include "real.h"

#define CRRMULRE(z, yre, yim)        ((z).Re*(yre)-(z).Im*(yim))
#define CRRMULIM(z, yre, yim)        ((z).Re*(yim)+(z).Im*(yre))
#define CCMULRE(z1, z2)        ((z1).Re*(z2).Re-(z1).Im*(z2).Im)
#define CCMULIM(z1, z2)        ((z1).Re*(z2).Im+(z1).Im*(z2).Re)

typedef struct {
  real Re;
  real Im;
} complex;

#endif /* COMPLEX_HDR_ALREADY_INCLUDED */
