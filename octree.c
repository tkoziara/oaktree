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

    if (src->what == HPL) cutoff *= 0.1; /* stricter test for flat surfaces */
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

/* create single node */
static struct node* node (struct node *n0, struct node *n1, struct node *n2, struct node *n3, struct node **list)
{
  struct node *node;

  ERRMEM  (node = malloc (sizeof (struct node)));

  node->node [0] = n0; 
  node->node [1] = n1;
  node->node [2] = n2;
  node->node [3] = n3;

  node->kind = 0;

  if (n0) node->kind ++;
  if (n1) node->kind ++;
  if (n2) node->kind ++;
  if (n3) node->kind ++;

#if DEBUG
  if (!(node->kind == 0 || node->kind == 2 || node->kind == 4))
  {
    fprintf (stderr, "Incorrect node kind, %d, !\n", node->kind);
  }
#endif

  node->next = *list;
  *list = node;

  return node;
}

/* insert nodes */
static void insert_nodes (struct octree *octree, struct solid *solid,
  struct node *n0, struct node *n1, struct node *n2, struct node *n3,
  struct node *n4, struct node *n5, struct node *n6, struct node *n7,
  struct node **list)
{
  struct element *element = octree->element;

  if (element && element->solid == solid)
  {
    if (!n0) n0 = node (NULL, NULL, NULL, NULL, list); /* can only happen from 0-level */
    if (!n1) n1 = node (NULL, NULL, NULL, NULL, list);
    if (!n2) n2 = node (NULL, NULL, NULL, NULL, list);
    if (!n3) n3 = node (NULL, NULL, NULL, NULL, list);
    if (!n4) n4 = node (NULL, NULL, NULL, NULL, list);
    if (!n5) n5 = node (NULL, NULL, NULL, NULL, list);
    if (!n6) n6 = node (NULL, NULL, NULL, NULL, list);
    if (!n7) n7 = node (NULL, NULL, NULL, NULL, list);

    element->node [0] = n0;
    element->node [1] = n1;
    element->node [2] = n2;
    element->node [3] = n3;
    element->node [4] = n4;
    element->node [5] = n5;
    element->node [6] = n6;
    element->node [7] = n7;
  }
  else if (octree->down [0])
  {

    enum {O0=0x01, O1=0x02, O2=0x04, O3=0x08, O4=0x10, O5=0x20, O6=0x40, O7=0x80} code = 0x00;

    if (octree->down [0]->element && octree->down [0]->element->solid == solid) code |= O0;
    if (octree->down [1]->element && octree->down [1]->element->solid == solid) code |= O1;
    if (octree->down [2]->element && octree->down [2]->element->solid == solid) code |= O2;
    if (octree->down [3]->element && octree->down [3]->element->solid == solid) code |= O3;
    if (octree->down [4]->element && octree->down [4]->element->solid == solid) code |= O4;
    if (octree->down [5]->element && octree->down [5]->element->solid == solid) code |= O5;
    if (octree->down [6]->element && octree->down [6]->element->solid == solid) code |= O6;
    if (octree->down [7]->element && octree->down [7]->element->solid == solid) code |= O7;

    struct node *n01 = (code & (O0|O1)) || (n0 && n1) ? node (n0, n1, NULL, NULL, list) : NULL,
		*n12 = (code & (O1|O2)) || (n1 && n2) ? node (n1, n2, NULL, NULL, list) : NULL,
		*n23 = (code & (O2|O3)) || (n2 && n3) ? node (n2, n3, NULL, NULL, list) : NULL,
		*n30 = (code & (O0|O3)) || (n3 && n0) ? node (n3, n0, NULL, NULL, list) : NULL,
		*n45 = (code & (O4|O5)) || (n4 && n5) ? node (n4, n5, NULL, NULL, list) : NULL,
		*n56 = (code & (O5|O6)) || (n5 && n6) ? node (n5, n6, NULL, NULL, list) : NULL,
		*n67 = (code & (O6|O7)) || (n6 && n7) ? node (n6, n7, NULL, NULL, list) : NULL,
		*n74 = (code & (O7|O4)) || (n7 && n4) ? node (n7, n4, NULL, NULL, list) : NULL,
		*n04 = (code & (O0|O4)) || (n0 && n4) ? node (n0, n4, NULL, NULL, list) : NULL,
		*n15 = (code & (O1|O5)) || (n1 && n5) ? node (n1, n5, NULL, NULL, list) : NULL,
		*n26 = (code & (O2|O6)) || (n2 && n6) ? node (n2, n6, NULL, NULL, list) : NULL,
		*n37 = (code & (O3|O7)) || (n3 && n7) ? node (n3, n7, NULL, NULL, list) : NULL,
		*n0123 = (code & (O0|O1|O2|O3)) || (n0 && n1 && n2 && n3) ? node (n0, n1, n2, n3, list) : NULL,
		*n4567 = (code & (O4|O5|O6|O7)) || (n4 && n5 && n6 && n7) ? node (n4, n5, n6, n7, list) : NULL,
		*n0145 = (code & (O0|O1|O4|O5)) || (n0 && n1 && n4 && n5) ? node (n0, n1, n4, n5, list) : NULL,
		*n1256 = (code & (O1|O2|O5|O6)) || (n1 && n2 && n5 && n6) ? node (n1, n2, n5, n6, list) : NULL,
		*n2367 = (code & (O2|O3|O6|O7)) || (n2 && n3 && n6 && n7) ? node (n2, n3, n6, n7, list) : NULL,
		*n3074 = (code & (O3|O0|O7|O4)) || (n3 && n0 && n7 && n4) ? node (n3, n0, n7, n4, list) : NULL,
		*mid = code ? node (NULL, NULL, NULL, NULL, list) : NULL;

    if (!n0 && (code & O0)) n0 = node (NULL, NULL, NULL, NULL, list);
    if (!n1 && (code & O1)) n1 = node (NULL, NULL, NULL, NULL, list);
    if (!n2 && (code & O2)) n2 = node (NULL, NULL, NULL, NULL, list);
    if (!n3 && (code & O3)) n3 = node (NULL, NULL, NULL, NULL, list);
    if (!n4 && (code & O4)) n4 = node (NULL, NULL, NULL, NULL, list);
    if (!n5 && (code & O5)) n5 = node (NULL, NULL, NULL, NULL, list);
    if (!n6 && (code & O6)) n6 = node (NULL, NULL, NULL, NULL, list);
    if (!n7 && (code & O7)) n7 = node (NULL, NULL, NULL, NULL, list);

    insert_nodes (octree->down [0], solid, n0, n01, n0123, n30, n04, n0145, mid, n3074, list);
    insert_nodes (octree->down [1], solid, n01, n1, n12, n0123, n0145, n15, n1256, mid, list);
    insert_nodes (octree->down [2], solid, n0123, n12, n2, n23, mid, n1256, n26, n2367, list);
    insert_nodes (octree->down [3], solid, n30, n0123, n23, n3, n3074, mid, n2367, n37, list);
    insert_nodes (octree->down [4], solid, n04, n0145, mid, n3074, n4, n45, n4567, n74, list);
    insert_nodes (octree->down [5], solid, n0145, n15, n1256, mid, n45, n5, n56, n4567, list);
    insert_nodes (octree->down [6], solid, mid, n1256, n26, n2367, n4567, n56, n6, n67, list);
    insert_nodes (octree->down [7], solid, n3074, mid, n2367, n37, n74, n4567, n67, n7, list);
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

/* insert solid and refine octree down to a cutoff edge length */
struct node* octree_insert_solid (struct octree *octree, struct solid *solid, REAL cutoff)
{
  REAL t [5][3][3], p [8][3], q [2][3], (*d) [8], (*s) [4][3], *x = octree->extents;
  int i, j, k, l, n, m, o, size;
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

  n = shape_unique_leaves (solid->shape, q[0], LEN (q[1]), &leaf, &i);
  if (n == 0)
  {
    if (i)
    {
      ERRMEM (element = malloc (sizeof (struct element)));
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

  if (triang) /* triangulation was created */
  {
    ERRMEM (element = malloc (sizeof (struct element)));
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
  if (!octree->up)
  {
    struct node *list = NULL;

    insert_nodes (octree, solid, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &list);

    return list;
  }
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
