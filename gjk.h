/*
 * gjk.h
 * ------
 */

#ifndef __gjk__
#define __gjk__

/* (a,na) and (b,nb) are the two input tables of polyhedrons vertices;
 * 'p' and 'q' are the two outputed closest points, respectively in
 * polyhedron (a,na) and polyhedron (b,nb); the distance is returned */
REAL gjk (REAL epsilon, REAL *a, int na, REAL *b, int nb, REAL *p, REAL *q);

#endif
