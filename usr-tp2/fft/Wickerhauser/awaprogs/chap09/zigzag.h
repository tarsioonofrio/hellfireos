/*
 * 
 * This header file defines zig-zag scan order in two dimensions
 *
 * Copyright (C) 1991--94 Wickerhauser Consulting.  All Rights Reserved.
 * May be freely copied for personal, noncommercial use by owners of 
 * ``Adapted Wavelet Analysis from Theory to Software'' ISBN 1-56881-041-5
 * by Mladen Victor Wickerhauser [AK Peters, Ltd., Wellesley, Mass., 1994]
 */

#ifndef ZIGZAG_HDR_ALREADY_INCLUDED
# define ZIGZAG_HDR_ALREADY_INCLUDED

typedef struct {
  int i;
  int j;
} ipoint;			/* Integer 2d lattice point */

static ipoint zigzag8x8[64] = {
  {0,0},
  {0,1}, {1,0}, 
  {2,0}, {1,1}, {0,2},
  {0,3}, {1,2}, {2,1}, {3,0},
  {4,0}, {3,1}, {2,2}, {1,3}, {0,4},
  {0,5}, {1,4}, {2,3}, {3,2}, {4,1}, {5,0},
  {6,0}, {5,1}, {4,2}, {3,3}, {2,4}, {1,5}, {0,6},
  {0,7}, {1,6}, {2,5}, {3,4}, {4,3}, {5,2}, {6,1}, {7,0},
  {7,1}, {6,2}, {5,3}, {4,4}, {3,5}, {2,6}, {1,7},
  {2,7}, {3,6}, {4,5}, {5,4}, {6,3}, {7,2},
  {7,3}, {6,4}, {5,5}, {4,6}, {3,7},
  {4,7}, {5,6}, {6,5}, {7,4},
  {7,5}, {6,6}, {5,7},
  {6,7}, {7,6},
  {7,7}
};


#endif /* ZIGZAG_HDR_ALREADY_INCLUDED */
