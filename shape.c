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

/* return shape leaves with counter */
static void shape_leaves_with_counter (struct shape *shape, struct shape **leaves, int *i)
{
  switch (shape->what)
  {
  case ADD:
  case MUL:
  case SUB:
    shape_leaves_with_counter (shape->left, leaves, i);
    shape_leaves_with_counter (shape->right, leaves, i);
    break;
  case HPL:
  case SPH:
    leaves [*i] = shape;
    (*i) ++;
    break;
  }
}

/* oobb rotate */
inline static void oobb_rotate (REAL oobb [8][3], REAL *point, REAL *matrix)
{
  REAL v [3];
  int i;

  for (i = 0; i < 8; i++)
  {
    SUB (oobb[i], point, v);
    NVADDMUL (point, matrix, v, oobb[i]);
  }
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

/* count shape leaves */
int shape_leaves_count (struct shape *shape)
{
  switch (shape->what)
  {
  case ADD:
  case MUL:
  case SUB:
    return shape_leaves_count (shape->left) + shape_leaves_count (shape->right);
    break;
  case HPL:
  case SPH:
    return 1;
    break;
  }

  return 0;
}

/* return shape leaves */
void shape_leaves (struct shape *shape, struct shape **leaves)
{
  int i = 0;

  shape_leaves_with_counter (shape, leaves, &i);
}

/* rotate shape about a point using a rotation matrix */
void shape_rotate (struct shape *shape, REAL *point, REAL *matrix)
{
  REAL v [3];

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
      oobb_rotate (data->oobb, point, matrix);
    }
    break;
  case SPH:
    {
      struct sphere *data = shape->data;

      SUB (data->c, point, v);
      NVADDMUL (point, matrix, v, data->c);
      oobb_rotate (data->oobb, point, matrix);
    }
    break;
  }
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

  if (shape->label) free (shape->label);

  free (shape);
}
