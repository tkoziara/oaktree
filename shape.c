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
  struct superellipsoid *super;
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
  case ELP:
    super = shape->data;
    SUB (super->c, point, z);
    a = DOT (super->u, z);
    b = DOT (super->v, z);
    c = DOT (super->w, z);
    a = fabs (a);
    b = fabs (b);
    c = fabs (c);
    v = pow (pow (a, super->p) + pow (b, super->p), super->q / super->p) * pow (c, super->q) - 1.0;
    break;
  }

  return v;
}
