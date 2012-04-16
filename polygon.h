/*
 * polygon.h
 * ---------
 */

#ifndef __polygon__
#define __polygon__

/* polygonise a grid cell (p, val) into up to five triangles (number returned);
 * zero is returned if the grid cell is either above or below the specified isolevel */
int polygonise (REAL p[8][3], REAL val [8], REAL isolevel, REAL cutoff, REAL triangles [5][3][3]);

#endif
