/*
 * octree.c
 * --------
 */

#include <stdlib.h>
#include <float.h>
#include <stdio.h>
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

/* recursive insert */
static void recursive_insert (struct octree *oct, struct shape *shape, REAL cutoff)
{
  REAL p [8][3], d [8], e [6], f [6];
  REAL *extents = oct->extents;
  struct octcut *cut;
  int i;

  COPY6 (oct->extents, e);
  COPY6 (shape->extents, f);

  if (f [3] < e [0] || f [0] > e [3] ||
      f [4] < e [1] || f [1] > e [4] ||
      f [5] < e [2] || f [2] > e [5]) return;

  VECTOR (p[0], e[0], e[1], e[2]);
  VECTOR (p[1], e[0], e[4], e[2]);
  VECTOR (p[2], e[3], e[4], e[2]);
  VECTOR (p[3], e[3], e[1], e[2]);
  VECTOR (p[4], e[0], e[1], e[5]);
  VECTOR (p[5], e[0], e[4], e[5]);
  VECTOR (p[6], e[3], e[4], e[5]);
  VECTOR (p[7], e[3], e[1], e[5]);

  for (i = 0; i < 8; i ++) d [i] = shape_evaluate (shape, p[i]);

  if (accurate (p, d, shape, cutoff))
  {
    for (i = 1; i < 8; i ++)
      if (d [0] * d [i] <= 0.0) break;

    if (i < 8)
    {
      ERRMEM (cut = malloc (sizeof (struct octcut)));

      for (i = 0; i < 8; i ++) cut->d [i] = d[i];

      cut->shape = shape;
      cut->next = oct->cut;
      oct->cut = cut;
    }
  }
  else
  {
    if (!oct->down [0])
    {
      SUB (e+3, e, d);

      SCALE (d, 0.5);

      VECTOR (e, extents[0], extents[1], extents[2]);
      VECTOR (e+3, e[0]+d[0], e[1]+d[1], e[2]+d[2]);
      oct->down [0] = octree_create (e, FLT_MAX);
      oct->down [0]->up = oct;

      VECTOR (e, extents[0]+d[0], extents[1], extents[2]);
      VECTOR (e+3, e[0]+d[0], e[1]+d[1], e[2]+d[2]);
      oct->down [1] = octree_create (e, FLT_MAX);
      oct->down [1]->up = oct;

      VECTOR (e, extents[0]+d[0], extents[1]+d[1], extents[2]);
      VECTOR (e+3, e[0]+d[0], e[1]+d[1], e[2]+d[2]);
      oct->down [2] = octree_create (e, FLT_MAX);
      oct->down [2]->up = oct;

      VECTOR (e, extents[0], extents[1]+d[1], extents[2]);
      VECTOR (e+3, e[0]+d[0], e[1]+d[1], e[2]+d[2]);
      oct->down [3] = octree_create (e, FLT_MAX);
      oct->down [3]->up = oct;

      VECTOR (e, extents[0], extents[1], extents[2]+d[2]);
      VECTOR (e+3, e[0]+d[0], e[1]+d[1], e[2]+d[2]);
      oct->down [4] = octree_create (e, FLT_MAX);
      oct->down [4]->up = oct;

      VECTOR (e, extents[0]+d[0], extents[1], extents[2]+d[2]);
      VECTOR (e+3, e[0]+d[0], e[1]+d[1], e[2]+d[2]);
      oct->down [5] = octree_create (e, FLT_MAX);
      oct->down [5]->up = oct;

      VECTOR (e, extents[0]+d[0], extents[1]+d[1], extents[2]+d[2]);
      VECTOR (e+3, e[0]+d[0], e[1]+d[1], e[2]+d[2]);
      oct->down [6] = octree_create (e, FLT_MAX);
      oct->down [6]->up = oct;

      VECTOR (e, extents[0], extents[1]+d[1], extents[2]+d[2]);
      VECTOR (e+3, e[0]+d[0], e[1]+d[1], e[2]+d[2]);
      oct->down [7] = octree_create (e, FLT_MAX);
      oct->down [7]->up = oct;
    }

    for (i = 0; i < 8; i ++) recursive_insert (oct->down [i], shape, cutoff);
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

/* insert list of shape and refine octree down to a cutoff edge length */
void octree_insert_shapes (struct octree *oct, struct shape *list, REAL cutoff)
{
  struct shape *shape;

  for (shape = list; shape; shape = shape->next)
  {
    recursive_insert (oct, shape, cutoff);
  }
}

/* free octree memory */
void octree_destroy (struct octree *oct)
{
}
