/*
 * shape.c
 * -------
 */

#include "oaktree.h"
#include "alg.h"
#include "err.h"

/* return distance to shape at given point, together with normal and color */
REAL shape_evaluate (struct shape *shp, REAL *point)
{
  struct superellipsoid *super;
  REAL a, b, c, v, z [3];

  switch (shp->what)
  {
  case ADD:
    a = shape_evaluate (shp->left, point);
    b = shape_evaluate (shp->right, point);
    c = a*a + b*b;
    v = (a + b + sqrt (c)) * c;
    break;
  case MUL:
    a = shape_evaluate (shp->left, point);
    b = shape_evaluate (shp->right, point);
    c = a*a + b*b;
    v = (a + b - sqrt (c)) * c;
    break;
  case SUB:
    a = shape_evaluate (shp->left, point);
    b = -shape_evaluate (shp->right, point);
    c = a*a + b*b;
    v = (a + b - sqrt (c)) * c;
    break;
  case ELP:
    super = shp->data;
    SUB (super->c, point, z);
    a = DOT (super->v[0], z);
    b = DOT (super->v[1], z);
    c = DOT (super->v[2], z);
    a = fabs (a);
    b = fabs (b);
    c = fabs (c);
    v = pow (pow (a, super->r) + pow (b, super->r), super->t / super->r) * pow (c, super->t) - 1.0;
    break;
  }

  return v;
}
