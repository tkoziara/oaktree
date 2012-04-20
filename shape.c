/*
 * shape.c
 * -------
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "oaktree.h"
#include "error.h"
#include "alg.h"

/* recursive copy */
static struct shape* recursive_copy (struct shape *shape)
{
  struct shape *copy;

  ERRMEM (copy = calloc (1, sizeof (struct shape*)));

  if (shape->extents)
  {
    ERRMEM (copy->extents = malloc (sizeof (REAL [6])));
    COPY6 (shape->extents, copy->extents);
  }

  copy->what = shape->what;

  switch (shape->what)
  {
  case ADD:
  case MUL:
  case SUB:
    copy->left = recursive_copy (shape->left);
    copy->right = recursive_copy (shape->right);
    break;
  case HPL:
    {
      struct halfplane *data;

      ERRMEM (data = malloc (sizeof (struct halfplane)));

      memcpy (data, shape->data, sizeof (struct halfplane));

      copy->data = data;
    }
    break;
  }

  return copy;
}

/* copy and label shape */
struct shape* shape_copy (struct shape *shape, char *label)
{
  struct shape *copy;

  copy = recursive_copy (shape);

  ERRMEM (copy->data = malloc (sizeof (char [strlen (label)])));

  strcpy (copy->data, label);

  return copy;
}

/* combine two shapes */
struct shape* shape_combine (struct shape *left, short what, struct shape *right)
{
  struct shape *shape;

  ERRMEM (shape = calloc (1, sizeof (struct shape)));
  shape->what = what;
  shape->left = left;
  shape->right = right;

  return shape;
}

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
