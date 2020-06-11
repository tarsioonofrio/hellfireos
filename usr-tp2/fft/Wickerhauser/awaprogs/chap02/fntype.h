/*
 * These are function type definitions.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */


#ifndef FNTYPE_HDR_ALREADY_INCLUDED
# define FNTYPE_HDR_ALREADY_INCLUDED

#include "real.h"

/* Function to deallocate a general type: */
typedef void *
  (* freetype)(void *);

/* Do-nothing function, when no deallocation is needed: */
extern void *
  freevoid(void *discard);

#endif /* FNTYPE_HDR_ALREADY_INCLUDED */
