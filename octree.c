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
      COPY (t[0], (*out) [*m][0]);
      COPY (t[1], (*out) [*m][1]);
      COPY (t[2], (*out) [*m][2]);
      COPY (n, (*out) [*m][3]);
      (*m) ++;
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

/* create node */
static struct node* node (REAL *p, REAL *y, struct element *e, REAL cutoff, struct node **list)
{
  struct node *node;
  REAL u;

  if (e)
  {
    u = p[0]-y[0];
    if (fabs (u) < cutoff)
    {
      u = p[1]-y[1];
      if (fabs (u) < cutoff)
      {
	u = p[2]-y[2];
	if (fabs (u) < cutoff) return e->node [0];
	else
	{
	  u = p[2]-y[5];
	  if (fabs (u) < cutoff) return e->node [4];
	}
      }
      else
      {
	u = p[1]-y[4];
	if (fabs (u) < cutoff)
	{
	  u = p[2]-y[2];
	  if (fabs (u) < cutoff) return e->node [1];
	  else
	  {
	    u = p[2]-y[5];
	    if (fabs (u) < cutoff) return e->node [5];
	  }
	}
      }
    }
    else
    {
      u = p[0]-y[3];
      if (fabs (u) < cutoff)
      {
	u = p[1]-y[1];
	if (fabs (u) < cutoff)
	{
	  u = p[2]-y[2];
	  if (fabs (u) < cutoff) return e->node [3];
	  else
	  {
	    u = p[2]-y[5];
	    if (fabs (u) < cutoff) return e->node [7];
	  }
	}
	else
	{
	  u = p[1]-y[4];
	  if (fabs (u) < cutoff)
	  {
	    u = p[2]-y[2];
	    if (fabs (u) < cutoff) return e->node [2];
	    else
	    {
	      u = p[2]-y[5];
	      if (fabs (u) < cutoff) return e->node [6];
	    }
	  }
	}
      }
    }

    ERRMEM  (node = calloc (1, sizeof (struct node)));

    node->hanging = 1;

    /* TODO: compute involved element nodes and coefficients */
  }
  else ERRMEM  (node = calloc (1, sizeof (struct node)));

  node->next = *list;
  *list = node;

  return node;
}

/* propagate element */
static void propagate (struct octree *octree, struct solid *solid, REAL *y, struct element *e, REAL cutoff, struct node **list)
{
  static short j [8][3] = {{0, 1, 2}, {0, 4, 2}, {3, 4, 2}, {3, 1, 2}, {0, 1, 5}, {0, 4, 5}, {3, 4, 5}, {3, 1, 5}};
  struct element *element = octree->element;
  REAL p [3], *x = octree->extents;
  int i;

  if (y[3] < x[0] ||
      y[4] < x[1] ||
      y[5] < x[2] ||
      y[0] > x[3] ||
      y[1] > x[4] ||
      y[2] > x[5] ||
      element == e) return;

  if (element && element->solid == solid)
  {
    for (i = 0; i < 8; i ++)
    {
      if (element->node [i]) continue;

      p [0] = x[j[i][0]];
      p [1] = x[j[i][1]];
      p [2] = x[j[i][2]];

      if (p[0] < y [0] ||
	  p[1] < y [1] ||
	  p[2] < y [2] ||
	  p[0] > y [3] ||
	  p[1] > y [4] ||
	  p[2] > y [5]) continue;

      element->node [i] = node (p, y, e, cutoff, list);
    }
  }
  else if (octree->down [0])
  {
    for (i = 0; i < 8; i ++) propagate (octree->down [i], solid, y, e, cutoff, list);
  }
}

/* auxiliary item */
struct item
{
  struct element *element;
  REAL *extents;
  short level;
  struct item *next;
};

/* collect items */
static void collect_items (short level, struct octree *octree, struct solid *solid, struct item **list)
{
  struct element *element = octree->element;
  struct item *item;
  int i;

  if (element && element->solid == solid)
  {
    ERRMEM (item = malloc (sizeof (struct item)));

    item->extents = octree->extents;
    item->element = element;
    item->level = level;

    item->next = *list;
    *list = item;
  }
  else if (octree->down [0])
  {
    for (i = 0; i < 8; i ++) collect_items (level+1, octree->down [i], solid, list);
  }
}

/* sort items by levels */
#define LE(i, j) (i)->level <= (j)->level
IMPLEMENT_LIST_SORT (SINGLE_LINKED, sort_items, struct item, prev, next, LE)

/* create nodes for solid elements */
static struct node* create_nodes (struct octree *octree, struct solid *solid, REAL cutoff)
{
  struct item *item = NULL, *next;
  struct node *list = NULL;
  int i;

  collect_items (0, octree, solid, &item);

  item = sort_items (item);

  for (; item; item = next)
  {
    next = item->next;

    for (i = 0; i < 8; i ++)
    {
      if (!item->element->node [i]) item->element->node [i] = node (NULL, NULL, NULL, cutoff, &list);
    }

    propagate (octree, solid, item->extents, item->element, cutoff, &list);

    free (item);
  }

  return list;
}

/* create octree down to a cutoff edge length */
struct octree* octree_create (REAL extents [6])
{
  struct octree *octree;

  ERRMEM (octree = calloc (1, sizeof (struct octree)));

  COPY6 (extents, octree->extents);

  return octree;
}

/* insert solid and refine octree down to a cutoff edge length */
struct node* octree_insert_solid (struct octree *octree, struct solid *solid, REAL cutoff)
{
  REAL t [5][3][3], p [8][3], q [2][3], (*d) [8], (*s) [4][3], *x = octree->extents;
  int i, j, k, l, n, m, o, size, inside;
  char allaccurate, *flagged;
  struct shape **leaf, **tmp;
  struct element *element;
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

  n = shape_unique_leaves (solid->shape, q[0], LEN (q[1]), &leaf, &inside);
  if (n == 0)
  {
    if (inside)
    {
      ERRMEM (element = calloc (1, sizeof (struct element)));
      element->triang = NULL;
      element->solid = solid;
      element->next = octree->element;
      octree->element = element;
    }

    goto done;
  }
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

  /* will just one primitive do ? */

  if (l > 1)
  {
    REAL z [8], w, v;

    for (j = 0; j < 8; j ++)  z [j] = shape_evaluate (solid->shape, p[j]); /* sample shape */

    for (i = 0; i < n; i ++)
    {
      if (flagged [i]) /* for flagged primitives */
      {
	for (w = 0, j = 0; j < 8; j ++)
	{
          v = d[i][j] - z[j];
	  w += fabs (v);  /* compute difference between primitive and shape */
	}

	if (w < cutoff) /* if small enough, use only this primitive */
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

  if (allaccurate && l > 8) allaccurate = 0; /* 8 is arbitrary XXX */

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
	  split (leaf[i], t[j], tmp, k, solid->shape, cutoff, &s, &m, &size); /* split against all other flagged primitives */

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
    ERRMEM (element = calloc (1, sizeof (struct element)));
    element->triang = triang;
    element->solid = solid;
    element->next = octree->element;
    octree->element = element;
  }
  else if (!allaccurate) /* not enough accuracy */
  {
    if (!octree->down [0])
    {
      if (q[1][0] <= cutoff && q[1][1] <= cutoff && q[1][2] <= cutoff) goto done;

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

    for (i = 0; i < 8; i ++) octree_insert_solid (octree->down [i], solid, cutoff);
  }

done:
  if (!octree->up) return create_nodes (octree, solid, cutoff);
  else return NULL;
}

/* free octree memory */
void octree_destroy (struct octree *octree)
{
  struct element *element, *next;
  int i;

  if (octree->down [0]) for (i = 0; i < 8; i ++) octree_destroy (octree->down [i]);

  for (element = octree->element; element; element = next)
  {
    next = element->next;
    if (element->triang)
    {
      free (element->triang->t);
      free (element->triang);
    }
    free (element);
  }

  free (octree);
}
