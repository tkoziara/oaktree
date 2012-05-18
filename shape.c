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
static int leaves_count (struct shape *shape)
{
  switch (shape->what)
  {
  case ADD:
  case MUL:
    return leaves_count (shape->left) + leaves_count (shape->right);
    break;
  case HSP:
  case SPH:
  case CYL:
  case FLT:
    return 1;
    break;
  }

  return 0;
}

/* output shape leaves overlapping (c,r) sphere and return their count */
static void leaves_within_sphere (struct shape *shape, REAL c [3], REAL r, struct shape **leaves, int *i)
{
  struct halfspace *halfspace;
  struct cylinder *cylinder;
  struct sphere *sphere;
  struct fillet *fillet;
  REAL d [4];

  switch (shape->what)
  {
  case ADD:

    if (shape->left->what == MUL)
    {
      d [0] = shape_evaluate (shape->left, c);
      if (fabs (d[0]) <= r) leaves_within_sphere (shape->left, c, r, leaves, i);
    }
    else leaves_within_sphere (shape->left, c, r, leaves, i);

    if (shape->right->what == MUL)
    {
      d [0] = shape_evaluate (shape->right, c);
      if (fabs (d[0]) <= r) leaves_within_sphere (shape->right, c, r, leaves, i);
    }
    else leaves_within_sphere (shape->right, c, r, leaves, i);

    break;
  case MUL:

    if (shape->left->what == ADD)
    {
      d [0] = shape_evaluate (shape->left, c);
      if (fabs (d[0]) <= r) leaves_within_sphere (shape->left, c, r, leaves, i);
    }
    else leaves_within_sphere (shape->left, c, r, leaves, i);

    if (shape->right->what == ADD)
    {
      d [0] = shape_evaluate (shape->right, c);
      if (fabs (d[0]) <= r) leaves_within_sphere (shape->right, c, r, leaves, i);
    }
    else leaves_within_sphere (shape->right, c, r, leaves, i);

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
  case FLT:
    fillet = shape->data;
    d[0] = shape_evaluate (shape->left, c);
    d[1] = shape_evaluate (shape->right, c);
    d[2] = fillet->r;

    if (d[2] > 0)
    {
      if (d[0] + d[1] < d[2] + r) /* if below line y = -x + fillet + r */
      {
         d[3] = d[2] - sqrt ((d[0]-d[2])*(d[0]-d[2])+(d[1]-d[2])*(d[1]-d[2]));
	if (fabs (d[3]) <= r) /* if within outer fillet cylinder */
	{
	  leaves [*i] = shape;
	  (*i) ++;
	}
      }
    }
    else
    {
      if (d[0] + d[1] > d[2] - r) /* if above libe y = -x - fillet - r */
      {
        d[3] = sqrt ((d[0]-d[2])*(d[0]-d[2])+(d[1]-d[2])*(d[1]-d[2])) + d[2];
	if (fabs (d[3]) <= r) /* if within inner fillet cylinder */
	{
	  leaves [*i] = shape;
	  (*i) ++;
	}
      }
    }

    break;
  }
}

/* compare leaves */
static int compare_leaves (struct shape **ll, struct shape **rr)
{
  REAL u, v, w, a [3];

  /* XXX => sensitive to model scale */

#if REAL == float
  #define EPS 1E-6
#else
  #define EPS 1E-10
#endif

  if ((*ll)->what < (*rr)->what) return -1;
  else if ((*ll)->what > (*rr)->what) return 1;
  else switch ((*ll)->what)
  {
    case HSP:
    {
      struct halfspace *l = (*ll)->data, *r = (*rr)->data;

      u = l->n[0] - r->n[0];
      if (fabs (u) < EPS)
      {
         u = l->n[1] - r->n[1];
	if (fabs (u) < EPS)
	{
           u = l->n[2] - r->n[2];
	  if (fabs (u) < EPS)
	  {
	    u = -DOT (l->p, l->n);
	    v = -DOT (r->p, r->n);
	    w = u - v; /* difference of values at (0, 0, 0) */

	    if (fabs (w) < EPS) return 0;
	    else return w < 0 ? -1 : 1;
	  }
          else return u < 0 ? -1 : 1;
	}
        else return u < 0 ? -1 : 1;
      }
      else return u < 0 ? -1 : 1;
    }
    break;
    case SPH:
    {
      struct sphere *l = (*ll)->data, *r = (*rr)->data;

      u = l->c[0] - r->c[0];
      if (fabs (u) < EPS)
      {
        u = l->c[1] - r->c[1];
	if (fabs (u) < EPS)
	{
          u = l->c[2] - r->c[2];
	  if (fabs (u) < EPS)
	  {
            u = l->r - r->r;
	    if (fabs (u) < EPS) return 0;
            else return u < 0 ? -1 : 1;
	  }
          else return u < 0 ? -1 : 1;
	}
        else return u < 0 ? -1 : 1;
      }
      else return u < 0 ? -1 : 1;
    }
    break;
    case CYL:
    {
      struct cylinder *l = (*ll)->data, *r = (*rr)->data;

      u = l->d[0] - r->d[0];
      if (fabs (u) < EPS)
      {
        u = l->d[1] - r->d[1];
	if (fabs (u) < EPS)
	{
          u = l->d[2] - r->d[2];
	  if (fabs (u) < EPS)
	  {
	    u = l->r - r->r;
	    if (fabs (u) < EPS)
	    {
	      SUB (l->p, r->p, a);

	      u = a[1]*l->d[2] - a[2]*l->d[1]; /* (l->p-r->p) x l->d == 0 */
	      if (fabs (u) < EPS)
	      {
	        v = a[2]*l->d[0] - a[0]*l->d[2];
	        if (fabs (v) < EPS)
		{
	          w = a[0]*l->d[1] - a[1]*l->d[0];
	          if (fabs (w) < EPS) return 0;
                  else return w < 0 ? -1 : 1;
		}
                else return v < 0 ? -1 : 1;
	      }
              else return u < 0 ? -1 : 1;
	    }
            else return u < 0 ? -1 : 1;
	  }
          else return u < 0 ? -1 : 1;
	}
        else return u < 0 ? -1 : 1;
      }
      else return u < 0 ? -1 : 1;
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
  case FLT:
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

/* find out relation between two leaves */
static short relation (struct shape *l, struct shape *r)
{
  struct shape **lpath, **rpath, *il, *jr;
  short i, j;

  for (i = 0, il = l; il; il = il->up) i ++;
  for (j = 0, jr = r; jr; jr = jr->up) j ++;

  ERRMEM (lpath = malloc ((i+j) * sizeof (struct shape*)));
  rpath = lpath + i;

  i = j = 0;

  while (l) lpath [i] = l, i ++, l = l->up;
  while (r) rpath [j] = r, j ++, r = r->up;

  while (lpath [--i] == rpath [--j]);

  i = lpath [i+1]->what;

  free (lpath);

  return i;
}

/* delete leaf */
static void delete (struct shape *leaf)
{
  struct shape *up = leaf->up,
	       *upup = up->up;

  if (upup)
  {
    if (up == upup->left)
    {
      if (leaf == up->left)
      {
	upup->left = up->right;
	up->right->up = upup;
      }
      else
      {
	upup->left = up->left;
	up->left->up = upup;
      }
    }
    else
    {
      if (leaf == up->left)
      {
	upup->right = up->right;
	up->right->up = upup;
      }
      else
      {
	upup->right = up->left;
	up->left->up = upup;
      }
    }

    free (leaf->data);
    free (leaf);
    free (up);
  }
}

/* remove duplicated leaves */
static void remove_duplicated_leaves (struct shape *shape)
{
  REAL c[3] = {0, 0, 0};
  struct shape **leaf;
  int i, j, k, n;
 
  n = leaves_count (shape);

  ERRMEM (leaf = malloc (n * sizeof (struct shape*)));

  i = 0;

  leaves_within_sphere (shape, c, FLT_MAX, leaf, &i);

  qsort (leaf, i, sizeof (struct shape*), (int (*) (const void*, const void*)) compare_leaves);

  for (j = 0, k = 1; k < i; )
  {
    while (k < i && compare_leaves (&leaf[j], &leaf[k]) == 0)
    {
      if (!subtracted (leaf[j]) && subtracted (leaf [k]) && relation (leaf[j], leaf[k]) == MUL)
      {
	/* leaves are collected in the left-first order
	 * while inversions happen in subtractions on the right;
	 * in this case we can remove the subtracted conincident leaf */

        delete (leaf [k]);
      }
      else
      {
	if (leaf[j]->what == HSP)
	{
	  struct halfspace *hj = leaf[j]->data,
			   *hk = leaf[k]->data;
	  printf ("duplicated halfspace: (%g,%g,%g) (%d, %d)\n", hj->n[0], hj->n[1], hj->n[2], j, k);
	  printf ("their relation is %s\n", relation (leaf[j], leaf[k]) == ADD ? "ADD" : "MUL");
	  if (hj->s -1) printf ("ealier is inverted\n");
	  if (hk->s -1) printf ("later is inverted\n");
	}

	/* when duplicates are removed after each combination there
	 * should never be more than two identical leaves at a time */

	/* TODO: handle ADD case first => this will fix the union.py example
	 *      (excessive subdivisions on flat face due to limited radii of
	 *       participating coincident halfspaces);
	 *       ---
	 *       figure out how to generally handle ADD and MUL
	 */
      }

      k ++;
    }

    if (k < i)
    {
      leaf [j+1] = leaf [k];
      j ++;
      k ++;
    }
  }

  free (leaf);
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
    copy->left->up = copy;
    copy->right->up = copy;
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
  case FLT:
    {
      struct fillet *data;

      ERRMEM (data = malloc (sizeof (struct fillet)));

      memcpy (data, shape->data, sizeof (struct fillet));
      copy->data = data;

      copy->left = shape_copy (shape->left);
      copy->right = shape_copy (shape->right);
      copy->left->up = copy;
      copy->right->up = copy;
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
  case FLT:
    {
      struct fillet *data = shape->data;

      data->r *= -1.0;

      shape_invert (shape->left);
      shape_invert (shape->right);
    }
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
  left->up = shape;
  right->up = shape;

  remove_duplicated_leaves (shape);

  return shape;
}

/* move shape */
void shape_move (struct shape *shape, REAL *vector)
{
  switch (shape->what)
  {
  case ADD:
  case MUL:
  case FLT:
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
  case FLT:
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

/* insert fillet between surfaces overlapping (c, r) sphere */
void shape_fillet (struct shape *shape, REAL c [3], REAL r, REAL filletr, short scolor)
{
  struct fillet *fillet;
  struct shape *copy;
  struct shape **leaf;
  int i, j, n, rel;

  n = leaves_count (shape);

  ERRMEM (leaf = malloc (n * sizeof (struct shape*)));

  n = 0;

  leaves_within_sphere (shape, c, r, leaf, &n);

  rel = -1;

  for (i = 0; i < n; i ++)
  {
    if (leaf[i]->what == FLT) continue;

    for (j = i+1; j < n; j ++)
    {
      if (leaf[j]->what == FLT) continue;

      if (rel < 0) rel = relation (leaf [i], leaf[j]);
      else if (rel != relation (leaf [i], leaf [j]))
      {
	fprintf (stderr, "#############################################\n");
	fprintf (stderr, "# Cannot create a fillet on a saddle point! #\n");
	fprintf (stderr, "#############################################\n");
	return;
      }
    }
  }

  /* TODO: in fact the saddle point is sort of handled when you
   *       forget about (@@@); but note that in fillets of fillets
   *       we get distorted radii even for right angle faces;
   *       -------
   *       in shape_unique_leaves the corners (duples of triples of fillets)
   *       should be handled by creation of temporary rounding fillets:
   *       - spherical MUL for 3-MULs
   *       - sphereical ADD for 3-ADDs
   *       - ADD of ADDs for 2-ADDs
   *       - MUL of MULs for 2-MULs
   *       - ADD of MUL of MULs for ADD-2-MULs saddle
   *       - MUL of ADD of ADDs for MUL-2-ADDs saddle (is this possible?)
   */

  for (i = 0; i < n; i ++)
  {
    if (leaf[i]->what == FLT) continue;

    for (j = i+1; j < n; j ++)
    {
      if (leaf[j]->what == FLT) continue;

#if 0
      copy = shape_copy (leaf [i]);
      copy->up = leaf[i];
      leaf[i]->right = copy;
      ERRMEM (copy = calloc (1, sizeof (struct shape)));
      copy->up = leaf[i];
      leaf[i]->left = copy;
      copy->what = FLT;
      copy->left = shape_copy (leaf[i]->right);
      copy->right = shape_copy (leaf[j]);
      ERRMEM (fillet = malloc (sizeof (struct fillet)));
      fillet->scolor = scolor;
      copy->data = fillet;
      leaf[i]->what = relation (leaf[i], leaf[j]);
      if (leaf[i]->what == ADD) fillet->r = filletr;
      else fillet->r = -filletr;
      free (leaf[i]->data);
      leaf[i]->data = NULL;
#else
      /* change leaf[j] into the fillet;
       * put the original leaf[j] on the right;
       * put the other leaf[i] on the left */

      copy = shape_copy (leaf [j]);
      copy->up = leaf[j];
      leaf[j]->right = copy;
      ERRMEM (copy = calloc (1, sizeof (struct shape)));
      copy->up = leaf[j];
      leaf[j]->left = copy;
      copy->what = FLT;
      copy->left = shape_copy (leaf[j]->right);
      copy->right = shape_copy (leaf[i]);
      ERRMEM (fillet = malloc (sizeof (struct fillet)));
      fillet->scolor = scolor;
      copy->data = fillet;
      leaf[j]->what = relation (leaf[i], leaf[j]);
      if (leaf[j]->what == ADD) fillet->r = filletr;
      else fillet->r = -filletr;
      free (leaf[j]->data);
      leaf[j]->data = NULL;

      /* not to create fillets of fillets
       * come back to the original */
      leaf[j] = leaf[j]->right; /* (@@@) */

      /* because of subtractions on the right it is more
       * universal to always add fillets on the right;
       * yet when added on the left they can be combined
       * already here (leaf[i] becomes a fillet while it still iterates over leav[j]s);
       * this produces better rendering of corners but is errornous in subtractions */
#endif
    }
  }
}

/* return distance to shape at given point, together with normal and color */
REAL shape_evaluate (struct shape *shape, REAL *point)
{
  struct halfspace *halfspace;
  struct cylinder *cylinder;
  struct sphere *sphere;
  struct fillet *fillet;
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
  case FLT:
    fillet = shape->data;
    v = fillet->r;
    if (v > 0)
    {
      a = shape_evaluate (shape->left, point);
      b = shape_evaluate (shape->right, point);
      if (a > v || b > v) v = MIN (a, b);
      else v = v - sqrt((a-v)*(a-v)+(b-v)*(b-v));
    }
    else
    {
      a = shape_evaluate (shape->left, point);
      b = shape_evaluate (shape->right, point);
      if (a < v || b < v) v = MAX (a, b);
      else v = sqrt((a-v)*(a-v)+(b-v)*(b-v)) + v;
    }
    break;
  }

  return v;
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
  case FLT:
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

/* output unique shape leaves overlapping (c,r) sphere and return their count or inside flag if count is zero */
int shape_unique_leaves (struct shape *shape, REAL c [3], REAL r, struct shape ***leaves, char *inside)
{
  struct shape **leaf;
  int i, j, k, n;
  REAL v;
 
  v = shape_evaluate (shape, c);

  *inside = v < 0.0;

  if (v > r) return 0; /* outward distance filter */

  if (v < -r) return 0; /* inward distance filter */

  n = leaves_count (shape);

  ERRMEM ((*leaves) = leaf = malloc (n * sizeof (struct shape*)));

  i = 0;

  leaves_within_sphere (shape, c, r, leaf, &i);

  if (i == 0) free (leaf);

  if (i <= 1) return i;

  qsort (leaf, i, sizeof (struct shape*), (int (*) (const void*, const void*)) compare_leaves);

  for (j = 0, k = 1; k < i; )
  {
    while (k < i && compare_leaves (&leaf[j], &leaf[k]) == 0)
    {
      if (leaf [j]->what == FLT) break; /* TODO => perhaps temporarily MUL-combine all: leaf[j] = shape_combine (leaf[j], MUL, leaf [k]); */
      k ++;
    }

    if (k < i)
    {
      leaf [j+1] = leaf [k];
      j ++;
      k ++;
    }
  }

  return j+1;
}

/* test whether the leaf is in a union of shapes */
int shape_leaf_in_union (struct shape *leaf)
{
  if (subtracted (leaf)) /* skip subtracted branch */
  {
    do {leaf = leaf->up; } while (leaf && leaf->what == ADD);
  }

  while (leaf && leaf->what != ADD) leaf = leaf->up;

  return leaf ? 1 : 0;
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
  case FLT:
    shape_destroy (shape->left);
    shape_destroy (shape->right);
    free (shape->data);
    break;
  case HSP:
  case SPH:
  case CYL:
    free (shape->data);
    break;
  }

  free (shape);
}
