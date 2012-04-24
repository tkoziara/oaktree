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

/* linear shape functions for hexahedron */
#define HEX0(x,y,z) (0.125*(1.0-(x))*(1.0-(y))*(1.0-(z)))
#define HEX1(x,y,z) (0.125*(1.0+(x))*(1.0-(y))*(1.0-(z)))
#define HEX2(x,y,z) (0.125*(1.0+(x))*(1.0+(y))*(1.0-(z)))
#define HEX3(x,y,z) (0.125*(1.0-(x))*(1.0+(y))*(1.0-(z)))
#define HEX4(x,y,z) (0.125*(1.0-(x))*(1.0-(y))*(1.0+(z)))
#define HEX5(x,y,z) (0.125*(1.0+(x))*(1.0-(y))*(1.0+(z)))
#define HEX6(x,y,z) (0.125*(1.0+(x))*(1.0+(y))*(1.0+(z)))
#define HEX7(x,y,z) (0.125*(1.0-(x))*(1.0+(y))*(1.0+(z)))

/* accuracy test */
static int accurate (REAL p [8][3], REAL d [8], struct shape *shape, REAL cutoff)
{
  REAL l [][3] = {{-0.5, -0.5, -0.5}, {0.5, -0.5, -0.5}, {0.5, 0.5, -0.5}, {-0.5, 0.5, -0.5},
                  {-0.5, -0.5, 0.5}, {0.5, -0.5, 0.5}, {0.5, 0.5, 0.5}, {-0.5, 0.5, 0.5}, {0, 0, 0}};
  REAL q [3], u, v, w;
  int i;

#define point(l0, l1, l2, i) (p[0][i]*HEX0(l0, l1, l2) + p[1][i]*HEX1(l0, l1, l2) +\
                              p[2][i]*HEX2(l0, l1, l2) + p[3][i]*HEX3(l0, l1, l2) +\
		              p[4][i]*HEX4(l0, l1, l2) + p[5][i]*HEX5(l0, l1, l2) +\
		              p[6][i]*HEX6(l0, l1, l2) + p[7][i]*HEX7(l0, l1, l2))

#define interpolate(l0, l1, l2) (d[0]*HEX0(l0, l1, l2) + d[1]*HEX1(l0, l1, l2) +\
                                 d[2]*HEX2(l0, l1, l2) + d[3]*HEX3(l0, l1, l2) +\
				 d[4]*HEX4(l0, l1, l2) + d[5]*HEX5(l0, l1, l2) +\
				 d[6]*HEX6(l0, l1, l2) + d[7]*HEX7(l0, l1, l2))

  for (i = 0; i < 9; i ++)
  {
    q [0] = point (l[i][0], l[i][1], l[i][2], 0);
    q [1] = point (l[i][0], l[i][1], l[i][2], 1);
    q [2] = point (l[i][0], l[i][1], l[i][2], 2);
    u = interpolate (l[i][0], l[i][1], l[i][2]);
    v = shape_evaluate (shape, q);
    w = u - v;
    if (fabs (w) > cutoff) return 0;
  }

  return 1;
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
  REAL d [3], s [3][3][3];
  int i, j;

  if (k == 0)
  {
    if ((*m)+1 >= (*size))
    {
      (*size) *= 2;
      ERRMEM ((*out) = realloc ((*out), (*size) * sizeof (REAL [3][3])));
    }

    MID3 (t[0], t[1], t[2], d);

    short pass = 0;

    REAL n0 [3], n1 [3], v0 = shape_normal (shape, d, n0);

    if (src->what == HPL)
    {
      pass = (fabs (v0) <= 0.1*cutoff); /* mid-point and stringent tolerance for planar pieces */
    }
    else /* vertex points and relaxed tolerance for curved pieces */
    {
      REAL v[3] = {shape_evaluate (shape, t[0]),
	           shape_evaluate (shape, t[1]),
		   shape_evaluate (shape, t[2])};

      cutoff *= 1.25; /* XXX */

      pass = (fabs (v[0]) <= cutoff &&
	      fabs (v[1]) <= cutoff &&
	      fabs (v[2]) <= cutoff);
    }

    if (pass) /* difference of coincident surfaces will produce zero-measure zero-isosets */
    {         /* we try to eliminate those by perturbing the mid-point towards the insidide */
      ADDMUL (d, -cutoff, n0, d);
      v0 = shape_evaluate (shape, d);
      pass = (v0 < 0.0);
    }

    if (pass)
    {
      NORMAL (t[0], t[1], t[2], n1);

      if (DOT (n0, n1) > 0.0)
      {
	COPY (t[0], (*out) [*m][0]);
	COPY (t[1], (*out) [*m][1]);
	COPY (t[2], (*out) [*m][2]);
      }
      else
      {
	COPY (t[2], (*out) [*m][0]);
	COPY (t[1], (*out) [*m][1]);
	COPY (t[0], (*out) [*m][2]);
      }

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

/* create octree down to a cutoff edge length */
struct octree* octree_create (REAL extents [6], REAL cutoff)
{
  REAL  d [3], e [6];
  struct octree *oct;

  ERRMEM (oct = calloc (1, sizeof (struct octree)));

  COPY6 (extents, oct->extents);

  SUB (extents+3, extents, d);

  if (d [0] > cutoff || d [1] > cutoff || d [2] > cutoff)
  {
    SCALE (d, 0.5);

    VECTOR (e, extents[0], extents[1], extents[2]);
    VECTOR (e+3, e[0]+d[0], e[1]+d[1], e[2]+d[2]);
    oct->down [0] = octree_create (e, cutoff);
    oct->down [0]->up = oct;

    VECTOR (e, extents[0]+d[0], extents[1], extents[2]);
    VECTOR (e+3, e[0]+d[0], e[1]+d[1], e[2]+d[2]);
    oct->down [1] = octree_create (e, cutoff);
    oct->down [1]->up = oct;

    VECTOR (e, extents[0]+d[0], extents[1]+d[1], extents[2]);
    VECTOR (e+3, e[0]+d[0], e[1]+d[1], e[2]+d[2]);
    oct->down [2] = octree_create (e, cutoff);
    oct->down [2]->up = oct;

    VECTOR (e, extents[0], extents[1]+d[1], extents[2]);
    VECTOR (e+3, e[0]+d[0], e[1]+d[1], e[2]+d[2]);
    oct->down [3] = octree_create (e, cutoff);
    oct->down [3]->up = oct;

    VECTOR (e, extents[0], extents[1], extents[2]+d[2]);
    VECTOR (e+3, e[0]+d[0], e[1]+d[1], e[2]+d[2]);
    oct->down [4] = octree_create (e, cutoff);
    oct->down [4]->up = oct;

    VECTOR (e, extents[0]+d[0], extents[1], extents[2]+d[2]);
    VECTOR (e+3, e[0]+d[0], e[1]+d[1], e[2]+d[2]);
    oct->down [5] = octree_create (e, cutoff);
    oct->down [5]->up = oct;

    VECTOR (e, extents[0]+d[0], extents[1]+d[1], extents[2]+d[2]);
    VECTOR (e+3, e[0]+d[0], e[1]+d[1], e[2]+d[2]);
    oct->down [6] = octree_create (e, cutoff);
    oct->down [6]->up = oct;

    VECTOR (e, extents[0], extents[1]+d[1], extents[2]+d[2]);
    VECTOR (e+3, e[0]+d[0], e[1]+d[1], e[2]+d[2]);
    oct->down [7] = octree_create (e, cutoff);
    oct->down [7]->up = oct;
  }
      
  return oct;
}

/* insert shape and refine octree down to a cutoff edge length */
void octree_insert_shape (struct octree *oct, struct shape *shape, REAL cutoff)
{
  REAL t [5][3][3], p [8][3], (*d) [8], (*s) [3][3];
  REAL *x = oct->extents, e [6], z [3];
  int i, j, k, l, n, m, o, size;
  struct shape **leaf, **tmp;
  struct triang *triang;
  short allacurate;
  char *flagged;

  VECTOR (p[0], x[0], x[1], x[2]);
  VECTOR (p[1], x[0], x[4], x[2]);
  VECTOR (p[2], x[3], x[4], x[2]);
  VECTOR (p[3], x[3], x[1], x[2]);
  VECTOR (p[4], x[0], x[1], x[5]);
  VECTOR (p[5], x[0], x[4], x[5]);
  VECTOR (p[6], x[3], x[4], x[5]);
  VECTOR (p[7], x[3], x[1], x[5]);

  n = shape_unique_leaves (shape, p, cutoff, &leaf);
  size = 100;

  ERRMEM (flagged = calloc (n, 1))
  ERRMEM (tmp = malloc (n * sizeof (struct shape*)));
  ERRMEM (d = malloc (n * sizeof (REAL [8])));
  ERRMEM (s = malloc (size * sizeof (REAL [3][3])));

  allacurate = 1;
  triang = NULL;

  for (i = 0; i < n; i ++)
  {
    for (j = 0; j < 8; j ++) d [i][j] = shape_evaluate (leaf[i], p[j]);

    if (!accurate (p, d[i], leaf[i], cutoff))  /* but not accurate enough */
    {
      allacurate = 0;
      break;
    }

    for (j = 1; j < 8; j ++)
    {
      if (d [i][0] * d [i][j] <= 0.0) /* contains 0-isosurface */
      {
	flagged [i] = 1;
	break;
      }
    }
  }

  /* if allcaturate == 1 extract triangles */

  if (allacurate)
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
	  split (leaf[i], t[j], tmp, k, shape, cutoff, &s, &m, &size);

	  if (m)
	  {
	    if (!triang) ERRMEM (triang = calloc (1, sizeof (struct triang)));

	    ERRMEM (triang->t = realloc (triang->t, (triang->n+m)* sizeof (REAL [3][3])));

	    for (o = 0; o < m; o ++)
	    {
	      COPY (s [o][0], triang->t [triang->n+o][0]);
	      COPY (s [o][1], triang->t [triang->n+o][1]);
	      COPY (s [o][2], triang->t [triang->n+o][2]);
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
    triang->shape = shape;
    triang->next = oct->triang;
    oct->triang = triang;
  }
  else if (!allacurate) /* not enough accuracy */
  {
    if (!oct->down [0])
    {
      COPY6 (x, e);

      SUB (e+3, e, z);

      SCALE (z, 0.5);

      if (z [0] <= cutoff && z [1] <= cutoff && z [2] <= cutoff) return;

      VECTOR (e, x[0], x[1], x[2]);
      VECTOR (e+3, e[0]+z[0], e[1]+z[1], e[2]+z[2]);
      oct->down [0] = octree_create (e, FLT_MAX);
      oct->down [0]->up = oct;

      VECTOR (e, x[0]+z[0], x[1], x[2]);
      VECTOR (e+3, e[0]+z[0], e[1]+z[1], e[2]+z[2]);
      oct->down [1] = octree_create (e, FLT_MAX);
      oct->down [1]->up = oct;

      VECTOR (e, x[0]+z[0], x[1]+z[1], x[2]);
      VECTOR (e+3, e[0]+z[0], e[1]+z[1], e[2]+z[2]);
      oct->down [2] = octree_create (e, FLT_MAX);
      oct->down [2]->up = oct;

      VECTOR (e, x[0], x[1]+z[1], x[2]);
      VECTOR (e+3, e[0]+z[0], e[1]+z[1], e[2]+z[2]);
      oct->down [3] = octree_create (e, FLT_MAX);
      oct->down [3]->up = oct;

      VECTOR (e, x[0], x[1], x[2]+z[2]);
      VECTOR (e+3, e[0]+z[0], e[1]+z[1], e[2]+z[2]);
      oct->down [4] = octree_create (e, FLT_MAX);
      oct->down [4]->up = oct;

      VECTOR (e, x[0]+z[0], x[1], x[2]+z[2]);
      VECTOR (e+3, e[0]+z[0], e[1]+z[1], e[2]+z[2]);
      oct->down [5] = octree_create (e, FLT_MAX);
      oct->down [5]->up = oct;

      VECTOR (e, x[0]+z[0], x[1]+z[1], x[2]+z[2]);
      VECTOR (e+3, e[0]+z[0], e[1]+z[1], e[2]+z[2]);
      oct->down [6] = octree_create (e, FLT_MAX);
      oct->down [6]->up = oct;

      VECTOR (e, x[0], x[1]+z[1], x[2]+z[2]);
      VECTOR (e+3, e[0]+z[0], e[1]+z[1], e[2]+z[2]);
      oct->down [7] = octree_create (e, FLT_MAX);
      oct->down [7]->up = oct;
    }

    for (i = 0; i < 8; i ++) octree_insert_shape (oct->down [i], shape, cutoff);
  }
}

/* free octree memory */
void octree_destroy (struct octree *oct)
{
}
