/*
 * oaktree.h
 * ---------
 */

#ifndef __oaktree__
#define __oaktree__

struct halfplane
{
  REAL p [3], n [3];

  short scolor;
};

struct shape
{
  REAL extents [6];

  enum {ADD, MUL, SUB, HPL} what;

  void *data;

  short vcolor;

  struct shape *left, *right;

  struct shape *next;
};

/* return distance to shape at given point */
REAL shape_evaluate (struct shape *shape, REAL *point);

struct octcut
{
  REAL d [8];

  short color;

  struct shape *shape;

  struct octcut *next;
};

struct octree
{
  REAL extents [6];

  struct octcut *cut;

  struct octree *up, *down [8];
};

/* create octree down to a cutoff edge length */
struct octree* octree_create (REAL extents [6], REAL cutoff);

/* insert list of shape and refine octree down to a cutoff edge length */
void octree_insert_shapes (struct octree *oct, struct shape *list, REAL cutoff);

/* free octree memory */
void octree_destroy (struct octree *oct);

struct simulation
{
  char *outpath;

  REAL duration, step;

  REAL grid, cutoff;

  REAL extents [6];

  struct shape *shape;
};

#endif
