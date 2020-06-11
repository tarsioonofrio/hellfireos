/*
 * Some trivial deallocation or do-nothing functions.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#include <stdlib.h>
#include "fntype.h"

extern void *
  freevoid(
	   void *discard)
{
  if(discard) free(discard);
  return(0);
}
