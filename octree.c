/*
 * octree.c
 * --------
 */

#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>
#include "polygon.h"
#include "oaktree.h"
#include "error.h"
#include "sort.h"
#include "alg.h"

#define PRIMITIVES_PER_OCTANT 4

/* accuracy test */
static int accurate (REAL q [3], REAL d [8], struct shape *shape, REAL cutoff)
{
  REAL u, v, w;

  u = 0.125 * (d[0]+d[1]+d[2]+d[3]+d[4]+d[5]+d[6]+d[7]);
  v = shape_evaluate (shape, q);
  w = u - v;
  if (fabs (w) > cutoff) return 0;
  else return 1;
}

/* find zero point for u * v < 0 */
inline static void zeropoint (REAL a [3], REAL b [3], REAL u, REAL v, REAL z [3])
{
  REAL ba [3], mu;

  SUB (b, a, ba);
  mu = -u / (v - u);
  ADDMUL (a, mu, ba, z);
}

/* split triangle */
static void split (struct shape *src, REAL t [3][3], struct shape **leaf, int k, struct shape *shape, REAL cutoff, REAL (**out) [3][3], int *m, int *size)
{
  REAL s [3][3][3], d [3], n [3], v;
  int i, j;

  if (k == 0)
  {
    if ((*m)+1 >= (*size))
    {
      (*size) *= 2;
      ERRMEM ((*out) = realloc ((*out), (*size) * sizeof (REAL [3][3])));
    }

    MID3 (t[0], t[1], t[2], d);

    if (src == NULL) /* called from trim */
    {
      v = shape_evaluate (shape, d);
      if (v < -cutoff) /* include inner bits only */
      {
	COPY (t[0], (*out) [*m][0]);
	COPY (t[1], (*out) [*m][1]);
	COPY (t[2], (*out) [*m][2]);
	(*m) ++;
      }
      return;
    }

    NORMAL (t[0], t[1], t[2], n);

    NORMALIZE (n);

    /* difference of coincident surfaces will produce zero-measure zero-isosets */
    /* we try to eliminate those by perturbing the mid-point towards the insidide */

    if (src->what == HSP) cutoff *= 0.1; /* stricter test for flat surfaces */
    else cutoff *= ALG_SQR2; /* relaxed test for curved surfaces (cutoff**3 cube has sqrt(2)*cutoff diameter) XXX */

    SUBMUL (d, cutoff, n, d);
    v = shape_evaluate (shape, d);
    if (v < 0.0 && fabs (cutoff + v) < cutoff) /* inside but not too deep */
    {
      i = 1;

      if (shape_leaf_in_union (src)) /* unions may produce internal boundaries (0-levels) */
      {                              /* in this case we also test positive perturbation */
	ADDMUL (d, 2.0*cutoff, n, d);
	v = shape_evaluate (shape, d);
	i = (v > 0.0 && fabs (cutoff - v) < cutoff); /* outside but not too far */
      }

      if (i)
      {
	COPY (t[0], (*out) [*m][0]);
	COPY (t[1], (*out) [*m][1]);
	COPY (t[2], (*out) [*m][2]);
	(*m) ++;
      }
    }
  }
  else
  {
    d [0] = shape_evaluate (leaf[0], t[0]);
    d [1] = shape_evaluate (leaf[0], t[1]);
    d [2] = shape_evaluate (leaf[0], t[2]);

    if (fabs (d[0]) < cutoff)
    {
      if (fabs (d[1]) < cutoff || fabs (d[2]) < cutoff || d[1]*d[2] > 0.0)
      {
	COPY (t[0], s[0][0]);
	COPY (t[1], s[0][1]);
	COPY (t[2], s[0][2]);
	j = 1; goto done;
      }
      else
      {
	COPY (t[0], s[0][0]);
	COPY (t[1], s[0][1]);
	zeropoint (t[1], t[2], d[1], d[2], s[0][2]);
	COPY (t[0], s[1][0]);
	COPY (s[0][2], s[1][1]);
	COPY (t[2], s[1][2]);
	j = 2; goto done;
      }
    }

    if (fabs (d[1]) < cutoff)
    {
      if (fabs (d[2]) < cutoff || d[0]*d[2] > 0.0)
      {
	COPY (t[0], s[0][0]);
	COPY (t[1], s[0][1]);
	COPY (t[2], s[0][2]);
	j = 1; goto done;
      }
      else
      {
	COPY (t[0], s[0][0]);
	COPY (t[1], s[0][1]);
	zeropoint (t[2], t[0], d[2], d[0], s[0][2]);
	COPY (s[0][2], s[1][0]);
	COPY (t[1], s[1][1]);
	COPY (t[2], s[1][2]);
	j = 2; goto done;
      }
    }

    if (fabs (d[2]) < cutoff)
    {
      if (d[0]*d[1] > 0.0)
      {
	COPY (t[0], s[0][0]);
	COPY (t[1], s[0][1]);
	COPY (t[2], s[0][2]);
	j = 1; goto done;
      }
      else
      {
	COPY (t[0], s[0][0]);
	zeropoint (t[0], t[1], d[0], d[1], s[0][1]);
	COPY (t[2], s[0][2]);
	COPY (s[0][1], s[1][0]);
	COPY (t[1], s[1][1]);
	COPY (t[2], s[1][2]);
	j = 2; goto done;
      }
    }

    if ((d[0] > 0.0 && d[1] > 0.0 && d[2] > 0.0) ||
	(d[0] < 0.0 && d[1] < 0.0 && d[2] < 0.0))
    {
	COPY (t[0], s[0][0]);
	COPY (t[1], s[0][1]);
	COPY (t[2], s[0][2]);
	j = 1; goto done;
    }

    if (d[1] * d[2] > 0.0)
    {
      COPY (t[0], s[0][0]);
      zeropoint (t[0], t[1], d[0], d[1], s[0][1]);
      zeropoint (t[2], t[0], d[2], d[0], s[0][2]);
      COPY (s[0][2], s[1][0]);
      COPY (s[0][1], s[1][1]);
      COPY (t[1], s[1][2]);
      COPY (s[0][2], s[2][0]);
      COPY (t[1], s[2][1]);
      COPY (t[2], s[2][2]);
      j = 3; goto done;
    }

    if (d[0] * d[2] > 0.0)
    {
      COPY (t[1], s[0][0]);
      zeropoint (t[1], t[2], d[1], d[2], s[0][1]);
      zeropoint (t[0], t[1], d[0], d[1], s[0][2]);
      COPY (s[0][2], s[1][0]);
      COPY (s[0][1], s[1][1]);
      COPY (t[2], s[1][2]);
      COPY (s[0][2], s[2][0]);
      COPY (t[2], s[2][1]);
      COPY (t[0], s[2][2]);
      j = 3; goto done;
    }

    if (d[0] * d[1] > 0.0)
    {
      COPY (t[2], s[0][0]);
      zeropoint (t[2], t[0], d[2], d[0], s[0][1]);
      zeropoint (t[1], t[2], d[1], d[2], s[0][2]);
      COPY (s[0][2], s[1][0]);
      COPY (s[0][1], s[1][1]);
      COPY (t[0], s[1][2]);
      COPY (s[0][2], s[2][0]);
      COPY (t[0], s[2][1]);
      COPY (t[1], s[2][2]);
      j = 3; goto done;
    }

done:
    for (i = 0; i < j; i ++) split (src, s [i], leaf+1, k-1, shape, cutoff, out, m, size);
  }
}

/* trim internal face of a boundary cell and return its area */
static REAL trim (struct cell *cell, REAL *x, int type, REAL cutoff, struct face *face)
{
  struct shape *leaf [PRIMITIVES_PER_OCTANT];
  REAL t0 [3][3], t1 [3][3], (*s) [3][3];
  REAL a, area = 0.0;
  int n, m, size;
  struct face *f;

  switch (type)
  {
  case -3: /* -xy */
    t0[0][0] = x[0]; t0[0][1] = x[1]; t0[0][2] = x[2];
    t0[1][0] = x[0]; t0[1][1] = x[4]; t0[1][2] = x[2];
    t0[2][0] = x[3]; t0[2][1] = x[1]; t0[2][2] = x[2];
    t1[0][0] = x[0]; t1[0][1] = x[4]; t1[0][2] = x[2];
    t1[1][0] = x[3]; t1[1][1] = x[4]; t1[1][2] = x[2];
    t1[2][0] = x[3]; t1[2][1] = x[1]; t1[2][2] = x[2];
  break;
  case -2: /* -xz */
    t0[0][0] = x[0]; t0[0][1] = x[1]; t0[0][2] = x[2];
    t0[1][0] = x[3]; t0[1][1] = x[1]; t0[1][2] = x[5];
    t0[2][0] = x[0]; t0[2][1] = x[1]; t0[2][2] = x[5];
    t1[0][0] = x[0]; t1[0][1] = x[1]; t1[0][2] = x[2];
    t1[1][0] = x[3]; t1[1][1] = x[1]; t1[1][2] = x[2];
    t1[2][0] = x[3]; t1[2][1] = x[1]; t1[2][2] = x[5];
  break;
  case -1: /* -yz */
    t0[0][0] = x[0]; t0[0][1] = x[4]; t0[0][2] = x[2];
    t0[1][0] = x[0]; t0[1][1] = x[1]; t0[1][2] = x[5];
    t0[2][0] = x[0]; t0[2][1] = x[4]; t0[2][2] = x[5];
    t1[0][0] = x[0]; t1[0][1] = x[4]; t1[0][2] = x[2];
    t1[1][0] = x[0]; t1[1][1] = x[1]; t1[1][2] = x[2];
    t1[2][0] = x[0]; t1[2][1] = x[1]; t1[2][2] = x[5];
  break;
  case 1: /* yz */
    t0[0][0] = x[3]; t0[0][1] = x[4]; t0[0][2] = x[2];
    t0[1][0] = x[3]; t0[1][1] = x[4]; t0[1][2] = x[5];
    t0[2][0] = x[3]; t0[2][1] = x[1]; t0[2][2] = x[5];
    t1[0][0] = x[3]; t1[0][1] = x[4]; t1[0][2] = x[2];
    t1[1][0] = x[3]; t1[1][1] = x[1]; t1[1][2] = x[5];
    t1[2][0] = x[3]; t1[2][1] = x[1]; t1[2][2] = x[2];
  break;
  case 2: /* xz */
    t0[0][0] = x[0]; t0[0][1] = x[4]; t0[0][2] = x[2];
    t0[1][0] = x[0]; t0[1][1] = x[4]; t0[1][2] = x[5];
    t0[2][0] = x[3]; t0[2][1] = x[4]; t0[2][2] = x[5];
    t1[0][0] = x[0]; t1[0][1] = x[4]; t1[0][2] = x[2];
    t1[1][0] = x[3]; t1[1][1] = x[4]; t1[1][2] = x[5];
    t1[2][0] = x[3]; t1[2][1] = x[4]; t1[2][2] = x[2];
  break;
  case 3: /* xy */
    t0[0][0] = x[0]; t0[0][1] = x[1]; t0[0][2] = x[5];
    t0[1][0] = x[3]; t0[1][1] = x[1]; t0[1][2] = x[5];
    t0[2][0] = x[0]; t0[2][1] = x[4]; t0[2][2] = x[5];
    t1[0][0] = x[0]; t1[0][1] = x[4]; t1[0][2] = x[5];
    t1[1][0] = x[3]; t1[1][1] = x[1]; t1[1][2] = x[5];
    t1[2][0] = x[3]; t1[2][1] = x[4]; t1[2][2] = x[5];
  break;
  }

  for (n = 0, f = cell->face; f; f = f->next)
  {
    if (f->leaf) leaf [n] = f->leaf, n ++;
  }

  if (n != 0)
  {
    m = 0;
    size = 8;
    ERRMEM (s = malloc (size * sizeof (REAL [3][3])));

    split (NULL, t0, leaf, n, cell->domain->shape, cutoff, &s, &m, &size);
    split (NULL, t1, leaf, n, cell->domain->shape, cutoff, &s, &m, &size);

    if (m != 0)
    {
      ERRMEM (face->t = malloc (m * sizeof (REAL [3][3])))

      for (n = 0; n < m; n ++)
      {
	COPY (s [n][0], face->t [n][0]);
	COPY (s [n][1], face->t [n][1]);
	COPY (s [n][2], face->t [n][2]);
	TRIANGLE_AREA (s[n][0], s[n][1], s[n][2], a);
	area += a;
      }

      face->n = m;
    }

    free (s);
  }
  else
  {
    ERRMEM (face->t = malloc (2 * sizeof (REAL [3][3])))

    COPY (t0 [0], face->t [0][0]);
    COPY (t0 [1], face->t [0][1]);
    COPY (t0 [2], face->t [0][2]);
    TRIANGLE_AREA (t0[0], t0[1], t0[2], a);
    area += a;
    COPY (t1 [0], face->t [1][0]);
    COPY (t1 [1], face->t [1][1]);
    COPY (t1 [2], face->t [1][2]);
    TRIANGLE_AREA (t1[0], t1[1], t1[2], a);
    area += a;

    face->n = 2;
  }

  return area;
}

/* invert input face into output face and return its area */
static REAL invert (struct face *in, struct face *out)
{
  int i;

  out->normal [0] = -in->normal[0];
  out->normal [1] = -in->normal[1];
  out->normal [2] = -in->normal[2];

  if (in->t)
  {
    ERRMEM (out->t = malloc (in->n * sizeof (REAL [3][3])));

    for (i = 0; i < in->n; i ++)
    {
      COPY (in->t[i][0], out->t[i][2]);
      COPY (in->t[i][1], out->t[i][1]);
      COPY (in->t[i][2], out->t[i][0]);
    }

    out->n = in->n;
  }

  return in->area;
}

/* drop (c, y) down the tree and complete adjacency */
static void drop (struct octree *octree, struct domain *domain, REAL cutoff, struct cell *c, REAL *y)
{
  struct cell *cell = octree->cell; /* current domain cells can only be the heads of octree cell lists  */
  REAL *x = octree->extents, p [3], n [3];
  struct face *face;
  int i, type;

  if (y[3] < x[0] || y[4] < x[1] || y[5] < x[2] || y[0] > x[3] || y[1] > x[4] || y[2] > x[5]) return; /* doesn't overlap */

  if (cell == c) return; /* slef */

  if (cell && cell->domain == domain) /* c and cell overlap in the same domain */
  {
    MID (x, x+3, p);

    i = 0;

    if (p[0] > y[0] && p[0] < y[3]) i |= 0x1;
    if (p[1] > y[1] && p[1] < y[4]) i |= 0x2;
    if (p[2] > y[2] && p[2] < y[5]) i |= 0x4;

    switch (i)
    {
    case 3: /* xy */
    if (p[2] < y[2]) { VECTOR (n, 0, 0, 1); type = 3; }
    else { VECTOR (n, 0, 0, -1); type = -3; }
    break;
    case 5: /* xz */
    if (p[1] < y[1]) { VECTOR (n, 0, 1, 0); type = 2; }
    else { VECTOR (n, 0, -1, 0); type = -2; }
    break;
    case 6: /* yz */
    if (p[0] < y[0]) { VECTOR (n, 1, 0, 0); type = 1; }
    else { VECTOR (n, -1, 0, 0); type = -1; }
    break;
    default:
    return; /* c and cell don't overlap through face */
    }

    ERRMEM (face = calloc (1, sizeof (struct face)));

    face->area = trim (cell, x, type, cutoff, face);

    if (face->area > 0.0)
    {
      COPY (n, face->normal);
      face->leaf = NULL;
      face->adj = c;
      face->next = cell->face;
      cell->face = face;

      ERRMEM (face = calloc  (1, sizeof (struct face)));
      face->area = invert (cell->face, face);
      face->leaf = NULL;
      face->adj = cell;
      face->next = c->face;
      c->face = face;
    }
    else free (face);
  }
  else if (octree->down [0])
  {
    for (i = 0; i < 8; i ++) drop (octree->down [i], domain, cutoff, c, y);
  }
}

/* auxiliary item */
struct item
{
  struct cell *cell;
  REAL *extents;
  short level;
  struct item *next;
};

/* collect cells from given domain into items list */
static void collect_items (struct octree *octree, short level, struct domain *domain, struct item **list)
{
  struct cell *cell = octree->cell;
  struct item *item;
  int i;

  if (cell && cell->domain == domain)
  {
    ERRMEM (item = malloc (sizeof (struct item)));

    item->extents = octree->extents;
    item->cell = cell;
    item->level = level;

    item->next = *list;
    *list = item;
  }
  else if (octree->down [0])
  {
    for (i = 0; i < 8; i ++) collect_items (octree->down [i], level+1, domain, list);
  }
}

/* sort items by levels */
#define LE(i, j) (i)->level <= (j)->level
IMPLEMENT_LIST_SORT (SINGLE_LINKED, sort_items, struct item, prev, next, LE)

/* create cell adjacency */
static void create_cell_adjacency (struct octree *octree, struct domain *domain, REAL cutoff)
{
  struct item *item = NULL, *next;

  collect_items (octree, 0, domain, &item);

  item = sort_items (item);

  for (; item; item = next)
  {
    next = item->next;

    drop (octree, domain, cutoff, item->cell, item->extents);

    free (item);
  }
}

/* create octree down to a cutoff edge length */
struct octree* octree_create (REAL extents [6])
{
  struct octree *octree;

  ERRMEM (octree = calloc (1, sizeof (struct octree)));

  COPY6 (extents, octree->extents);

  return octree;
}

/* insert domain and refine octree down to a cutoff edge length */
void octree_insert_domain (struct octree *octree, struct domain *domain, REAL cutoff)
{
  REAL t [5][3][3], p [8][3], q [2][3], (*d) [8], (*s) [3][3], *x = octree->extents, a;
  char allaccurate, inside, *flagged;
  int i, j, k, l, n, m, o, size;
  struct shape **leaf, **tmp;
  struct face *list, *face;
  struct cell *cell;

  VECTOR (p[0], x[0], x[1], x[2]);
  VECTOR (p[1], x[0], x[4], x[2]);
  VECTOR (p[2], x[3], x[4], x[2]);
  VECTOR (p[3], x[3], x[1], x[2]);
  VECTOR (p[4], x[0], x[1], x[5]);
  VECTOR (p[5], x[0], x[4], x[5]);
  VECTOR (p[6], x[3], x[4], x[5]);
  VECTOR (p[7], x[3], x[1], x[5]);

  MID (p[0], p[6], q[0]);
  SUB (q[0], p[0], q[1]);

  n = shape_unique_leaves (domain->shape, q[0], LEN (q[1]), &leaf, &inside);
  if (n == 0)
  {
    if (inside)
    {
      if (q[1][0] > domain->grid) goto recurse; /* assumption of cubic octants */

      ERRMEM (cell = calloc (1, sizeof (struct cell)));
      cell->octree = octree;
      cell->domain = domain;
      cell->face = NULL;
      cell->next = octree->cell;
      octree->cell = cell;
    }

    goto done;
  }
  else if (q[1][0] > domain->grid) goto recurse; /* assumption of cubic octants */

  size = 128;
  ERRMEM (flagged = calloc (n, 1))
  ERRMEM (tmp = malloc (n * sizeof (struct shape*)));
  ERRMEM (d = malloc (n * sizeof (REAL [8])));
  ERRMEM (s = malloc (size * sizeof (REAL [3][3])));

  allaccurate = 1;
  list = NULL;

  for (l = i = 0; i < n; i ++)
  {
    for (j = 0; j < 8; j ++) d [i][j] = shape_evaluate (leaf[i], p[j]);

    if (!accurate (q[0], d[i], leaf[i], cutoff))  /* but not accurate enough */
    {
      allaccurate = 0;
    }
    else for (j = 1; j < 8; j ++)
    {
      if (d [i][0] * d [i][j] <= 0.0) /* contains 0-isosurface */
      {
	flagged [i] = 1;
	l ++;
	break;
      }
    }
  }

  if (l > 1) /* test if just one leaf would do */
  {
    x = (REAL*)t;

    for (j = 0; j < 8; j ++)  x [j] = shape_evaluate (domain->shape, p[j]); /* sample shape */

    for (i = 0; i < n; i ++)
    {
      if (flagged [i]) /* for flagged leaves */
      {
	for (x[8] = 0, j = 0; j < 8; j ++)
	{
          x[9] = d[i][j] - x[j];
	  x[8] += fabs (x[9]);  /* compute difference between leaf and shape */
	}

	if (x[8] < cutoff) /* if small enough, use only this leaf */
	{
	  for (j = 0; j < 8; j ++) d[0][j] = d[i][j];
	  allaccurate = l = n = 1;
	  leaf [0] = leaf [i];
	  flagged [0] = 1;
	  break;
	}
      }
    }
  }

  /* recurse down the tree if too many leaves */

  if (allaccurate && l > PRIMITIVES_PER_OCTANT) allaccurate = 0; /* XXX: arbitrary threshold */

  /* if all leaves are accorate extract triangles */

  if (allaccurate)
  {
    for (i = 0; i < n; i ++)
    {
      if (flagged [i])
      {
	for (m = j = k = 0; j < n; j ++)
	{
	  if (flagged [j] && j != i)
	  {
	    tmp [k] = leaf [j];
	    k ++;
	  }
	}

	l = polygonise (p, d[i], 0.0, 0.01*cutoff, t);

	for (j = 0; j < l; j ++)
	{
	  split (leaf[i], t[j], tmp, k, domain->shape, cutoff, &s, &m, &size); /* split against all other flagged leaves */
        }

	if (m)
	{
	  ERRMEM (face = calloc (1, sizeof (struct face)));

	  ERRMEM (face->t = malloc (m * sizeof (REAL [3][3])));

	  face->area = 0;

	  for (o = 0; o < m; o ++)
	  {
	    COPY (s [o][0], face->t [o][0]);
	    COPY (s [o][1], face->t [o][1]);
	    COPY (s [o][2], face->t [o][2]);
	    TRIANGLE_AREA (s[o][0], s[o][1], s[o][2], a);
	    face->area += a;
	  }

	  leaf_normal (leaf[i], q[0], face->normal);
	  NORMALIZE (face->normal);
	  face->leaf = leaf[i];
	  face->adj = NULL;
	  face->n = m;
	  face->next = list;
	  list = face;
	}
      }
    }
  }

  free (flagged);
  free (leaf);
  free (tmp);
  free (d);
  free (s);

  if (list || (allaccurate && inside)) /* triangulation was created or inner octant */
  {
    ERRMEM (cell = calloc (1, sizeof (struct cell)));
    cell->octree = octree;
    cell->face = list;
    cell->domain = domain;
    cell->next = octree->cell;
    octree->cell = cell;
  }
  else if (!allaccurate) /* not enough accuracy */
  {
recurse:
    if (!octree->down [0])
    {
      if (q[1][0] <= cutoff) goto done; /* assumption of cubic octants */

      x = (REAL*) t;

      VECTOR (x, p[0][0], p[0][1], p[0][2]);
      VECTOR (x+3, q[0][0], q[0][1], q[0][2]);
      octree->down [0] = octree_create (x);
      octree->down [0]->up = octree;

      VECTOR (x, p[0][0], q[0][1], p[0][2]);
      VECTOR (x+3, q[0][0], p[6][1], q[0][2]);
      octree->down [1] = octree_create (x);
      octree->down [1]->up = octree;

      VECTOR (x, q[0][0], q[0][1], p[0][2]);
      VECTOR (x+3, p[6][0], p[6][1], q[0][2]);
      octree->down [2] = octree_create (x);
      octree->down [2]->up = octree;

      VECTOR (x, q[0][0], p[0][1], p[0][2]);
      VECTOR (x+3, p[6][0], q[0][1], q[0][2]);
      octree->down [3] = octree_create (x);
      octree->down [3]->up = octree;

      VECTOR (x, p[0][0], p[0][1], q[0][2]);
      VECTOR (x+3, q[0][0], q[0][1], p[6][2]);
      octree->down [4] = octree_create (x);
      octree->down [4]->up = octree;

      VECTOR (x, p[0][0], q[0][1], q[0][2]);
      VECTOR (x+3, q[0][0], p[6][1], p[6][2]);
      octree->down [5] = octree_create (x);
      octree->down [5]->up = octree;

      VECTOR (x, q[0][0], q[0][1], q[0][2]);
      VECTOR (x+3, p[6][0], p[6][1], p[6][2]);
      octree->down [6] = octree_create (x);
      octree->down [6]->up = octree;

      VECTOR (x, q[0][0], p[0][1], q[0][2]);
      VECTOR (x+3, p[6][0], q[0][1], p[6][2]);
      octree->down [7] = octree_create (x);
      octree->down [7]->up = octree;
    }

    for (i = 0; i < 8; i ++) octree_insert_domain (octree->down [i], domain, cutoff);
  }

done:
  if (!octree->up) create_cell_adjacency (octree, domain, cutoff);
}

/* free octree memory */
void octree_destroy (struct octree *octree)
{
  struct cell *cell;
  struct face *face;
  void *next;
  int i;

  if (octree->down [0]) for (i = 0; i < 8; i ++) octree_destroy (octree->down [i]);

  for (cell = octree->cell; cell; cell = next)
  {
    for (face = cell->face; face; face = next)
    {
      next = face->next;
      free (face->t);
      free (face);
    }
    next = cell->next;
    free (cell);
  }

  free (octree);
}
