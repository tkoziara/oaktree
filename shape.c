/*
 * shape.c
 * -------
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "oaktree.h"
#include "error.h"
#include "alg.h"

/* recursive copy */
static struct shape* recursive_copy (struct shape *shape)
{
  struct shape *copy;

  ERRMEM (copy = calloc (1, sizeof (struct shape)));

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
  case SPH:
    {
      struct sphere *data;

      ERRMEM (data = malloc (sizeof (struct sphere)));

      memcpy (data, shape->data, sizeof (struct sphere));
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

  if (label)
  {
    ERRMEM (copy->label = malloc (strlen (label) + 1));

    strcpy (copy->label, label);
  }

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
  struct sphere *sphere;
  REAL a, b, v, z [3];

  switch (shape->what)
  {
  case ADD:
    a = shape_evaluate (shape->left, point);
    b = shape_evaluate (shape->right, point);
    v = MIN (a, b);
    break;
  case MUL:
    a = shape_evaluate (shape->left, point);
    b = shape_evaluate (shape->right, point);
    v = MAX (a, b);
    break;
  case SUB:
    a = shape_evaluate (shape->left, point);
    b = -shape_evaluate (shape->right, point);
    v = MAX (a, b);
    break;
  case HPL:
    halfplane = shape->data;
    SUB (point, halfplane->p, z);
    v = DOT (z, halfplane->n);
    break;
  case SPH:
    sphere = shape->data;
    SUB (point, sphere->c, z);
    v = LEN (z) - sphere->r;
    break;
  }

  return v;
}

/* rotate shape about a point using a rotation matrix */
void shape_rotate (struct shape *shape, REAL *point, REAL *matrix)
{
  REAL v [3];

  if (shape->extents)
  {
    REAL p [8][3], *e = shape->extents;
    int i;

    VECTOR (p[0], e[0], e[1], e[2]);
    VECTOR (p[1], e[0], e[4], e[2]);
    VECTOR (p[2], e[3], e[4], e[2]);
    VECTOR (p[3], e[3], e[1], e[2]);
    VECTOR (p[4], e[0], e[1], e[5]);
    VECTOR (p[5], e[0], e[4], e[5]);
    VECTOR (p[6], e[3], e[4], e[5]);
    VECTOR (p[7], e[3], e[1], e[5]);

    e [0] = e [1] = e [2] = FLT_MAX;
    e [3] = e [4] = e [5] = -FLT_MAX;

    for (i = 0; i < 8; i ++)
    {
      SUB (p[i], point, v);
      NVADDMUL (point, matrix, v, p[i]);
      if (p[i][0] < e[0]) e[0] = p[i][0];
      if (p[i][1] < e[1]) e[1] = p[i][1];
      if (p[i][2] < e[2]) e[2] = p[i][2];
      if (p[i][0] > e[3]) e[3] = p[i][0];
      if (p[i][1] > e[4]) e[4] = p[i][1];
      if (p[i][2] > e[5]) e[5] = p[i][2];
    }
  }

  switch (shape->what)
  {
  case ADD:
  case MUL:
  case SUB:
    shape_rotate (shape->left, point, matrix);
    shape_rotate (shape->right, point, matrix);
    break;
  case HPL:
    {
      struct halfplane *data = shape->data;

      SUB (data->p, point, v);
      NVADDMUL (point, matrix, v, data->p);
      COPY (data->n, v);
      NVMUL (matrix, v, data->n);
    }
    break;
  case SPH:
    {
      struct sphere *data = shape->data;

      SUB (data->c, point, v);
      NVADDMUL (point, matrix, v, data->c);
    }
    break;

  }

}

/* free shape memory */
void shape_destroy (struct shape *shape)
{
  switch (shape->what)
  {
  case ADD:
  case MUL:
  case SUB:
    shape_destroy (shape->left);
    shape_destroy (shape->right);
    break;
  case HPL:
  case SPH:
    free (shape->data);
    break;
  }

  if (shape->extents) free (shape->extents);

  if (shape->label) free (shape->label);

  free (shape);
}
