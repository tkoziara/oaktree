/*
 * shape.c
 * -------
 */

#include "oaktree.h"
#include "error.h"
#include "alg.h"

/* return distance to shape at given point, together with normal and color */
REAL shape_evaluate (struct shape *shape, REAL *point)
{
  struct halfplane *halfplane;
  REAL a, b, c, v, z [3];

  switch (shape->what)
  {
  case ADD:
    a = shape_evaluate (shape->left, point);
    b = shape_evaluate (shape->right, point);
    c = a*a + b*b;
    v = (a + b + sqrt (c)) * c;
    break;
  case MUL:
    a = shape_evaluate (shape->left, point);
    b = shape_evaluate (shape->right, point);
    c = a*a + b*b;
    v = (a + b - sqrt (c)) * c;
    break;
  case SUB:
    a = shape_evaluate (shape->left, point);
    b = -shape_evaluate (shape->right, point);
    c = a*a + b*b;
    v = (a + b - sqrt (c)) * c;
    break;
  case HPL:
    halfplane = shape->data;
    SUB (point, halfplane->p, z);
    v = DOT (z, halfplane->n);
    break;
  }

  return v;
}
