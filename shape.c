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
#include "gjk.h"

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

/* count shape leaves */
static int shape_leaves_count (struct shape *shape)
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
  case CYL:
    return 1;
    break;
  }

  return 0;
}

/* return shape leaves with counter */
static void shape_leaves_with_counter (struct shape *shape, REAL p[8][3], REAL cutoff, struct shape **leaves, int *i)
{
  switch (shape->what)
  {
  case ADD:
  case MUL:
  case SUB:
    shape_leaves_with_counter (shape->left, p, cutoff, leaves, i);
    shape_leaves_with_counter (shape->right, p, cutoff, leaves, i);
    break;
  case HPL:
  case SPH:
  case CYL:
    {
      REAL a [3], b [3];

      if (gjk (0.01*cutoff, (REAL*)p, 8, shape->data, 8, a, b) < cutoff)
      {
	leaves [*i] = shape;
	(*i) ++;
      }
    }
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

/* compare leaves */
int compare_leaves (struct shape **sl, struct shape **sr)
{
  REAL v [6];

  if ((*sl)->what < (*sr)->what) return -1;
  else if ((*sl)->what > (*sr)->what) return 1;
  else switch ((*sl)->what)
  {
    case HPL:
    {
      struct halfplane *l = (*sl)->data,
		       *r = (*sr)->data;

      SUB (l->n, r->n, v);
      v [3] = -DOT (l->p, l->n);
      v [4] = -DOT (r->p, r->n);
      v [5] = v[3] - v[4];

      if (fabs (v[0]) < 1E-10) /* XXX */
      {
	if (fabs (v[1]) < 1E-10)
	{
	  if (fabs (v[2]) < 1E-10)
	  {
	    if (fabs (v[5]) < 1E-10) return 0;
	    else if (v[3] < v[4]) return -1;
	    else return 1;
	  }
	  else if (l->n[2] < r->n[2]) return -1;
	  else return 1;
	}
	else if (l->n[1] < r->n[1]) return -1;
	else return 1;
      }
      else if (l->n[0] < r->n[0]) return -1;
      else return 1;
    }
    break;
    case SPH:
    {
      struct sphere *l = (*sl)->data,
		    *r = (*sr)->data;

      SUB (l->c, r->c, v);
      v [3] = l->r - r->r;

      if (fabs (v[0]) < 1E-10)
      {
	if (fabs (v[1]) < 1E-10)
	{
	  if (fabs (v[2]) < 1E-10)
	  {
	    if (fabs (v[3]) < 1E-10) return 0;
	    else if (l->r < r->r) return -1;
	    else return 1;
	  }
	  else if (l->c[2] < r->c[2]) return -1;
	  else return 1;
	}
	else if (l->c[1] < r->c[1]) return -1;
	else return 1;
      }
      else if (l->c[0] < r->c[0]) return -1;
      else return 1;
    }
    break;
    case CYL:
    {
      struct cylinder *l = (*sl)->data,
		      *r = (*sr)->data;

      PRODUCT (l->d, r->d, v);
      v [3] = shape_evaluate (*sl, r->p);
      v [4] = l->r - r->r;

      if (fabs (v[0]) < 1E-10)
      {
	if (fabs (v[1]) < 1E-10)
	{
	  if (fabs (v[2]) < 1E-10)
	  {
	    if (fabs (v[3]) < 1E-10)
	    {
	      if (fabs (v[4]) < 1E-10) return 0;
	      else if (l->r < r->r) return -1;
	      else return 1;
	    }
	    else if (l->r < r->r) return -1;
	    else return 1;
	  }
	  else if (l->d[2] < r->d[2]) return -1;
	  else return 1;
	}
	else if (l->d[1] < r->d[1]) return -1;
	else return 1;
      }
      else if (l->d[0] < r->d[0]) return -1;
      else return 1;
    }
    break;
    default:
    break;
  }

  return 0;
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

/* output unique shape leaves inside of a grid cell and return their count */
int shape_unique_leaves (struct shape *shape, REAL p [8][3], REAL cutoff, struct shape ***leaves)
{
  int i = 0, j, k, n;
  struct shape **leaf;

  n = shape_leaves_count (shape);

  ERRMEM ((*leaves) = leaf = malloc (n * sizeof (struct shape*)));

  shape_leaves_with_counter (shape, p, cutoff, leaf, &i);

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
  case CYL:
    {
      struct cylinder *data = shape->data;

      SUB (data->p, point, v);
      NVADDMUL (point, matrix, v, data->p);
      COPY (data->d, v);
      NVMUL (matrix, v, data->d);
      oobb_rotate (data->oobb, point, matrix);
    }
    break;
  }
}

/* return distance to shape at given point, together with normal and color */
REAL shape_evaluate (struct shape *shape, REAL *point)
{
  struct halfplane *halfplane;
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
  case CYL:
    cylinder = shape->data;
    a = DOT (cylinder->d, cylinder->d);
    SUB (point, cylinder->p, z);
    b = DOT (z, cylinder->d);
    v = b / a;
    ADDMUL (cylinder->p, v, cylinder->d, z);
    SUB (point, z, z);
    v = LEN (z) - cylinder->r;
    break;
  }

  return v;
}

/* return value and compute shape normal at given point */
REAL shape_normal (struct shape *shape, REAL *point, REAL *normal)
{
  REAL l [3], r [3], z [3], a, b, sq, v;
  struct halfplane *halfplane;
  struct cylinder *cylinder;
  struct sphere *sphere;

  switch (shape->what)
  {
  case ADD:
    a = shape_normal (shape->left, point, l);
    b = shape_normal (shape->right, point, r);
    v = MIN (a, b);
#if 0
    sq = sqrt (a*a + b*b);
    a = (1.0 - a / sq);
    b = (1.0 - b / sq);
    normal [0] = a*l[0] + b*r[0];
    normal [1] = a*l[1] + b*r[1];
    normal [2] = a*l[2] + b*r[2];
#else
    if (a < b) { COPY (l, normal); }
    else { COPY (r, normal); }
#endif
    break;
  case MUL:
    a = shape_normal (shape->left, point, l);
    b = shape_normal (shape->right, point, r);
    v = MAX (a, b);
#if 0
    sq = sqrt (a*a + b*b);
    a = (1.0 + a / sq);
    b = (1.0 + b / sq);
    normal [0] = a*l[0] + b*r[0];
    normal [1] = a*l[1] + b*r[1];
    normal [2] = a*l[2] + b*r[2];
#else
    if (a > b) { COPY (l, normal); }
    else { COPY (r, normal); }
#endif
    break;
  case SUB:
    a = shape_normal (shape->left, point, l);
    b = -shape_normal (shape->right, point, r);
    v = MAX (a, b);
#if 0
    sq = sqrt (a*a + b*b);
    a = (1.0 + a / sq);
    b = (1.0 + b / sq);
    normal [0] = a*l[0] - b*r[0];
    normal [1] = a*l[1] - b*r[1];
    normal [2] = a*l[2] - b*r[2];
#else
    if (a > b) { COPY (l, normal); }
    else { MUL (r, -1, normal); }
#endif
    break;
  case HPL:
    halfplane = shape->data;
    SUB (point, halfplane->p, z);
    v = DOT (z, halfplane->n);
    COPY (halfplane->n, normal);
    break;
  case SPH:
    sphere = shape->data;
    SUB (point, sphere->c, z);
    sq = LEN (z);
    v = sq - sphere->r;
    DIV (z, sq, normal);
    break;
  case CYL:
    cylinder = shape->data;
    a = DOT (cylinder->d, cylinder->d);
    SUB (point, cylinder->p, z);
    b = DOT (z, cylinder->d);
    v = b / a;
    ADDMUL (cylinder->p, v, cylinder->d, z);
    SUB (point, z, z);
    sq = LEN (z);
    v = sq - cylinder->r;
    DIV (z, sq, normal);
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
  case CYL:
    free (shape->data);
    break;
  }

  if (shape->label) free (shape->label);

  free (shape);
}
