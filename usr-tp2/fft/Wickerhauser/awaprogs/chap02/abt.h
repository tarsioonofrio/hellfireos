/*
 * This header file declares macros used to index array binary trees.
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef ABT_HDR_ALREADY_INCLUDED
#define ABT_HDR_ALREADY_INCLUDED

#define abtblock(N,L,B) ( (L)*(N) + (B) * ((N)>>(L)) )
#define abtblength(N,L) ((N)>>(L))

#endif /* ABT_HDR_ALREADY_INCLUDED */

