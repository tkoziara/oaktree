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
static void split (struct shape *src, REAL t [3][3], struct shape **leaf, int k, struct shape *shape, REAL cutoff, REAL (**out) [4][3], int *m, int *size)
{
  REAL s [3][3][3], d [3], n [3], v;
  int i, j;

  if (k == 0)
  {
    if ((*m)+1 >= (*size))
    {
      (*size) *= 2;
      ERRMEM ((*out) = realloc ((*out), (*size) * sizeof (REAL [4][3])));
    }

    MID3 (t[0], t[1], t[2], d);

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
	COPY (n, (*out) [*m][3]);
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

/* drop (c, y, l) down the tree and complete adjacency */
static void drop (struct octree *octree, short level, struct domain *domain, struct cell *c, REAL *y, short l)
{
  struct cell *cell = octree->cell; /* current domain cells can only be the heads of octree cell lists  */
  REAL *x = octree->extents, p [3];
  int i;

  if (y[3] < x[0] || y[4] < x[1] || y[5] < x[2] || y[0] > x[3] || y[1] > x[4] || y[2] > x[5]) return; /* doesn't overlap */

  if (cell == c) return; /* slef */

  if (cell && cell->domain == domain) /* c and cell overlap in the same domain */
  {
    MID (x, x+3, p);

    i = 0;

    if (p[0] > y[0] && p[0] < y[3]) i ++;
    if (p[1] > y[1] && p[1] < y[4]) i ++;
    if (p[2] > y[2] && p[2] < y[5]) i ++;

    if (i != 2) return; /* c and cell don't overlap through face */

    ERRMEM (cell->adj = realloc (cell->adj, ((cell->nadj+1) * sizeof (struct cell*))));
    ERRMEM (c->adj = realloc (c->adj, ((c->nadj+1) * sizeof (struct cell*))));
    cell->adj [cell->nadj] = c;
    c->adj [c->nadj] = cell;
    cell->nadj ++;
    c->nadj ++;
  }
  else if (octree->down [0])
  {
    for (i = 0; i < 8; i ++) drop (octree->down [i], level+1, domain, c, y, l);
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
static void create_cell_adjacency (struct octree *octree, struct domain *domain)
{
  struct item *item = NULL, *next;

  collect_items (octree, 0, domain, &item);

  item = sort_items (item);

  for (; item; item = next)
  {
    next = item->next;

    drop (octree, 0, domain, item->cell, item->extents, item->level);

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
  REAL t [5][3][3], p [8][3], q [2][3], (*d) [8], (*s) [4][3], *x = octree->extents;
  char allaccurate, inside, *flagged;
  int i, j, k, l, n, m, o, size;
  struct shape **leaf, **tmp;
  struct cell *cell;
  struct triang *triang;

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
      cell->triang = NULL;
      cell->domain = domain;
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
  ERRMEM (s = malloc (size * sizeof (REAL [4][3])));

  allaccurate = 1;
  triang = NULL;

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

  if (allaccurate && l > 4) allaccurate = 0; /* 4 is arbitrary XXX */

  /* if all leaves are accorate extract triangles */

  if (allaccurate)
  {
    for (i = 0; i < n; i ++)
    {
      if (flagged [i])
      {
	for (j = k = 0; j < n; j ++)
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
	  m = 0;
	  split (leaf[i], t[j], tmp, k, domain->shape, cutoff, &s, &m, &size); /* split against all other flagged leaves */

	  if (m)
	  {
	    if (!triang) ERRMEM (triang = calloc (1, sizeof (struct triang)));

	    ERRMEM (triang->t = realloc (triang->t, (triang->n+m)* sizeof (REAL [4][3])));

	    for (o = 0; o < m; o ++)
	    {
	      COPY (s [o][0], triang->t [triang->n+o][0]);
	      COPY (s [o][1], triang->t [triang->n+o][1]);
	      COPY (s [o][2], triang->t [triang->n+o][2]);
	      COPY (s [o][3], triang->t [triang->n+o][3]);
	    }

	    triang->n += m;
	  }
	}
      }
    }
  }

  free (flagged);
  free (leaf);
  free (tmp);
  free (d);
  free (s);

  if (triang || (allaccurate && inside)) /* triangulation was created or inner octant */
  {
    ERRMEM (cell = calloc (1, sizeof (struct cell)));
    cell->triang = triang;
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
  if (!octree->up) create_cell_adjacency (octree, domain);
}

/* free octree memory */
void octree_destroy (struct octree *octree)
{
  struct cell *cell, *next;
  int i;

  if (octree->down [0]) for (i = 0; i < 8; i ++) octree_destroy (octree->down [i]);

  for (cell = octree->cell; cell; cell = next)
  {
    next = cell->next;
    if (cell->triang)
    {
      free (cell->triang->t);
      free (cell->triang);
      free (cell->adj);
    }
    free (cell);
  }

  free (octree);
}
