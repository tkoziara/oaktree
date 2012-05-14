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

/* count shape leaves */
static int shape_leaves_count (struct shape *shape)
{
  switch (shape->what)
  {
  case ADD:
  case MUL:
    return shape_leaves_count (shape->left) + shape_leaves_count (shape->right);
    break;
  case HSP:
  case SPH:
  case CYL:
    return 1;
    break;
  }

  return 0;
}

/* output shape leaves overlapping (c,r) sphere and return their count */
static void shape_leaves_with_counter (struct shape *shape, REAL c [3], REAL r, struct shape **leaves, int *i)
{
  struct halfspace *halfspace;
  struct cylinder *cylinder;
  struct sphere *sphere;
  REAL d [4];

  switch (shape->what)
  {
  case ADD:

    if (shape->left->what == MUL)
    {
      d [0] = shape_evaluate (shape->left, c);
      if (fabs (d[0]) <= r) shape_leaves_with_counter (shape->left, c, r, leaves, i);
    }
    else shape_leaves_with_counter (shape->left, c, r, leaves, i);

    if (shape->right->what == MUL)
    {
      d [0] = shape_evaluate (shape->right, c);
      if (fabs (d[0]) <= r) shape_leaves_with_counter (shape->right, c, r, leaves, i);
    }
    else shape_leaves_with_counter (shape->right, c, r, leaves, i);

    break;
  case MUL:

    if (shape->left->what == ADD)
    {
      d [0] = shape_evaluate (shape->left, c);
      if (fabs (d[0]) <= r) shape_leaves_with_counter (shape->left, c, r, leaves, i);
    }
    else shape_leaves_with_counter (shape->left, c, r, leaves, i);

    if (shape->right->what == ADD)
    {
      d [0] = shape_evaluate (shape->right, c);
      if (fabs (d[0]) <= r) shape_leaves_with_counter (shape->right, c, r, leaves, i);
    }
    else shape_leaves_with_counter (shape->right, c, r, leaves, i);

    break;
  case HSP:
    halfspace = shape->data;
    SUB (c, halfspace->p, d);
    d [3] = r + halfspace->r;

    if (DOT (d, d) <= d[3]*d[3])
    {
      d [3] = DOT (halfspace->n, d);

      if (fabs (d[3]) <= r)
      {
	leaves [*i] = shape;
	(*i) ++;
      }
    }
    break;
  case SPH:
    sphere = shape->data;
    SUB (c, sphere->c, d);
    d [3] = r + sphere->r;

    if (DOT (d, d) <= d[3]*d[3])
    {
      leaves [*i] = shape;
      (*i) ++;
    }
    break;
  case CYL:
    cylinder = shape->data;
    d [3] = shape_evaluate (shape, c);

    if (fabs (d[3]) <= r)
    {
      leaves [*i] = shape;
      (*i) ++;
    }
    break;
  }
}

/* compare leaves */
static int compare_leaves (struct shape **ll, struct shape **rr)
{
  REAL u, v, w, a [3];

  if ((*ll)->what < (*rr)->what) return -1;
  else if ((*ll)->what > (*rr)->what) return 1;
  else switch ((*ll)->what)
  {
    case HSP:
    {
      struct halfspace *l = (*ll)->data, *r = (*rr)->data;

      u = l->n[0] - r->n[0];
      if (fabs (u) < 1E-10) /* fine <= normalized */
      {
         u = l->n[1] - r->n[1];
	if (fabs (u) < 1E-10)
	{
           u = l->n[2] - r->n[2];
	  if (fabs (u) < 1E-10)
	  {
	    u = -DOT (l->p, l->n);
	    v = -DOT (r->p, r->n);
	    w = u - v; /* difference of values at (0, 0, 0) */

	    if (fabs (w) < 1E-10) return 0; /* this and below are sensitive to user scale XXX */
	  }
	}
      }

      if (l < r) return -1;
      else return 1;
    }
    break;
    case SPH:
    {
      struct sphere *l = (*ll)->data, *r = (*rr)->data;

      u = l->c[0] - r->c[0];
      if (fabs (u) < 1E-10)
      {
        u = l->c[1] - r->c[1];
	if (fabs (u) < 1E-10)
	{
          u = l->c[2] - r->c[2];
	  if (fabs (u) < 1E-10)
	  {
            u = l->r - r->r;
	    if (fabs (u) < 1E-10) return 0;
	  }
	}
      }

      if (l < r) return -1;
      else return 1;
    }
    break;
    case CYL:
    {
      struct cylinder *l = (*ll)->data, *r = (*rr)->data;

      u = l->d[0] - r->d[0];
      if (fabs (u) < 1E-10)
      {
        u = l->d[1] - r->d[1];
	if (fabs (u) < 1E-10)
	{
          u = l->d[2] - r->d[2];
	  if (fabs (u) < 1E-10)
	  {
	    u = l->r - r->r;
	    if (fabs (u) < 1E-10)
	    {
	      SUB (l->p, r->p, a);

	      u = a[1]*l->d[2] - a[2]*l->d[1]; /* (l->p-r->p) x l->d == 0 */
	      if (fabs (u) < 1E-10)
	      {
	        v = a[2]*l->d[0] - a[0]*l->d[2];
	        if (fabs (v) < 1E-10)
		{
	          w = a[0]*l->d[1] - a[1]*l->d[0];
	          if (fabs (w) < 1E-10) return 0;
		}
	      }
	    }
	  }
	}
      }

      if (l < r) return -1;
      else return 1;
    }
    break;
    default:
    break;
  }

  return 0;
}

/* check whether entier shape has been subtracted */
static int subtracted (struct shape *shape)
{
  switch (shape->what)
  {
  case ADD:
  case MUL:
    return subtracted (shape->left) && subtracted (shape->right);
  case HSP:
    {
      struct halfspace *data = shape->data;

      return data->s == -1;
    }
    break;
  case SPH:
    {
      struct sphere *data = shape->data;

      return data->s == -1;
    }
    break;
  case CYL:
    {
      struct cylinder *data = shape->data;

      return data->s == -1;
    }
    break;
  }

  return 0;
}

/* copy shape */
struct shape* shape_copy (struct shape *shape)
{
  struct shape *copy;

  ERRMEM (copy = calloc (1, sizeof (struct shape)));

  copy->what = shape->what;

  switch (shape->what)
  {
  case ADD:
  case MUL:
    copy->left = shape_copy (shape->left);
    copy->right = shape_copy (shape->right);
    break;
  case HSP:
    {
      struct halfspace *data;

      ERRMEM (data = malloc (sizeof (struct halfspace)));

      memcpy (data, shape->data, sizeof (struct halfspace));
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
  case CYL:
    {
      struct cylinder *data;

      ERRMEM (data = malloc (sizeof (struct cylinder)));

      memcpy (data, shape->data, sizeof (struct cylinder));
      copy->data = data;
    }
    break;
  }

  return copy;
}

/* return the same inverted shape */
struct shape* shape_invert (struct shape *shape)
{
  switch (shape->what)
  {
  case ADD:
    shape->what = MUL;

    shape_invert (shape->left);
    shape_invert (shape->right);
    break;
  case MUL:
    shape->what = ADD;

    shape_invert (shape->left);
    shape_invert (shape->right);
    break;
  case HSP:
    {
      struct halfspace *data = shape->data;

      data->s *= -1.0;
    }
    break;
  case SPH:
    {
      struct sphere *data = shape->data;

      data->s *= -1.0;
    }
    break;
  case CYL:
    {
      struct cylinder *data = shape->data;

      data->s *= -1.0;
    }
    break;
  }

  return shape;
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

/* move shape */
void shape_move (struct shape *shape, REAL *vector)
{
  switch (shape->what)
  {
  case ADD:
  case MUL:
    shape_move (shape->left, vector);
    shape_move (shape->right, vector);
    break;
  case HSP:
    {
      struct halfspace *data = shape->data;

      ACC (vector, data->p);
    }
    break;
  case SPH:
    {
      struct sphere *data = shape->data;

      ACC (vector, data->c);
    }
    break;
  case CYL:
    {
      struct cylinder *data = shape->data;

      ACC (vector, data->p);
    }
    break;
  }
}

/* rotate shape about a point using a rotation matrix */
void shape_rotate (struct shape *shape, REAL *point, REAL *matrix)
{
  REAL v [3];

  switch (shape->what)
  {
  case ADD:
  case MUL:
    shape_rotate (shape->left, point, matrix);
    shape_rotate (shape->right, point, matrix);
    break;
  case HSP:
    {
      struct halfspace *data = shape->data;

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
  case CYL:
    {
      struct cylinder *data = shape->data;

      SUB (data->p, point, v);
      NVADDMUL (point, matrix, v, data->p);
      COPY (data->d, v);
      NVMUL (matrix, v, data->d);
    }
    break;
  }
}

/* return distance to shape at given point, together with normal and color */
REAL shape_evaluate (struct shape *shape, REAL *point)
{
  struct halfspace *halfspace;
  struct cylinder *cylinder;
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
  case HSP:
    halfspace = shape->data;
    SUB (point, halfspace->p, z);
    v = halfspace->s * DOT (z, halfspace->n);
    break;
  case SPH:
    sphere = shape->data;
    SUB (point, sphere->c, z);
    v = sphere->s * (LEN (z) - sphere->r);
    break;
  case CYL:
    cylinder = shape->data;
    SUB (point, cylinder->p, z);
    a = DOT (z, cylinder->d);
    SUBMUL (z, a, cylinder->d, z);
    b = LEN (z);
    a = cylinder->r;
    if (b < a) b = 0.5*((b*b)/a + a); /* smooth out inside with y = x**2/(2r) + r/2 */
    v = cylinder->s * (b - a);
    break;
  }

  return v;
}

/* return value and compute shape normal at given point */
REAL shape_normal (struct shape *shape, REAL *point, REAL *normal)
{
  REAL l [3], r [3], z [3], a, b, sq, v;
  struct halfspace *halfspace;
  struct cylinder *cylinder;
  struct sphere *sphere;

  switch (shape->what)
  {
  case ADD:
    a = shape_normal (shape->left, point, l);
    b = shape_normal (shape->right, point, r);
    v = MIN (a, b);
    if (a < b) { COPY (l, normal); }
    else { COPY (r, normal); }
    break;
  case MUL:
    a = shape_normal (shape->left, point, l);
    b = shape_normal (shape->right, point, r);
    v = MAX (a, b);
    if (a > b) { COPY (l, normal); }
    else { COPY (r, normal); }
    break;
  case HSP:
    halfspace = shape->data;
    SUB (point, halfspace->p, z);
    v = halfspace->s * DOT (z, halfspace->n);
    MUL (halfspace->n, halfspace->s, normal);
    break;
  case SPH:
    sphere = shape->data;
    SUB (point, sphere->c, z);
    sq = LEN (z);
    v = sphere->s * (sq - sphere->r);
    sq *= sphere->s;
    DIV (z, sq, normal);
    break;
  case CYL:
    cylinder = shape->data;
    SUB (point, cylinder->p, z);
    a = DOT (z, cylinder->d);
    SUBMUL (z, a, cylinder->d, z);
    b = LEN (z);
    sq = b * cylinder->s;
    DIV (z, sq, normal);
    a = cylinder->r;
    if (b < a) b = 0.5*((b*b)/a + a); /* smooth out inside with y = x**2/(2r) + r/2 */
    v = cylinder->s * (b - a);
    break;
  }

  return v;
}

/* output unique shape leaves overlapping (c,r) sphere and return their count or inside flag if count is zero */
int shape_unique_leaves (struct shape *shape, REAL c [3], REAL r, struct shape ***leaves, int *inside)
{
  struct shape **leaf;
  int i, j, k, n;
  REAL v;
 
  v = shape_evaluate (shape, c);

  *inside = v < 0.0;

  if (v > r) return 0; /* outward distance filter */

  if (v < -r) return 0; /* inward distance filter */

  n = shape_leaves_count (shape);

  ERRMEM ((*leaves) = leaf = malloc (n * sizeof (struct shape*)));

  i = 0;

  shape_leaves_with_counter (shape, c, r, leaf, &i);

  if (i <= 1) return i;

  qsort (leaf, i, sizeof (struct shape*), (int (*) (const void*, const void*)) compare_leaves);

  for (j = 0, k = 1; k < i; )
  {
    while (k < i && compare_leaves (&leaf[j], &leaf[k]) == 0) k ++;

    if (k < i)
    {
      leaf [j+1] = leaf [k];
      j ++;
      k ++;
    }
  }

  return j+1;
}

/* compute shape extents */
void shape_extents (struct shape *shape, REAL *extents)
{
  struct halfspace *halfspace;
  struct sphere *sphere;
  REAL l [6], r [6];

  switch (shape->what)
  {
  case ADD:
  case MUL:
    shape_extents (shape->left, l);
    shape_extents (shape->right, r);

    if (shape->what == MUL && subtracted (shape->right))
    {
      COPY6 (l, extents); /* no contrubution from B in A - B */
    }
    else
    {
      if (l[0] < r[0]) extents[0] = l[0];
      else extents[0] = r[0];
      if (l[1] < r[1]) extents[1] = l[1];
      else extents[1] = r[1];
      if (l[2] < r[2]) extents[2] = l[2];
      else extents[2] = r[2];
      if (l[3] > r[3]) extents[3] = l[3];
      else extents[3] = r[3];
      if (l[4] > r[4]) extents[4] = l[4];
      else extents[4] = r[4];
      if (l[5] > r[5]) extents[5] = l[5];
      else extents[5] = r[5];
    }
    break;
  case HSP:
    halfspace = shape->data;
    extents [0] = halfspace->p[0] - halfspace->r;
    extents [1] = halfspace->p[1] - halfspace->r;
    extents [2] = halfspace->p[2] - halfspace->r;
    extents [3] = halfspace->p[0] + halfspace->r;
    extents [4] = halfspace->p[1] + halfspace->r;
    extents [5] = halfspace->p[2] + halfspace->r;
    break;
  case SPH:
    sphere = shape->data;
    extents [0] = sphere->c[0] - sphere->r;
    extents [1] = sphere->c[1] - sphere->r;
    extents [2] = sphere->c[2] - sphere->r;
    extents [3] = sphere->c[0] + sphere->r;
    extents [4] = sphere->c[1] + sphere->r;
    extents [5] = sphere->c[2] + sphere->r;
    break;
  case CYL:
    /* cylinder is skipped */
    extents [0] = FLT_MAX;
    extents [1] = FLT_MAX;
    extents [2] = FLT_MAX;
    extents [3] = -FLT_MAX;
    extents [4] = -FLT_MAX;
    extents [5] = -FLT_MAX;
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
    shape_destroy (shape->left);
    shape_destroy (shape->right);
    break;
  case HSP:
  case SPH:
  case CYL:
    free (shape->data);
    break;
  }

  free (shape);
}
