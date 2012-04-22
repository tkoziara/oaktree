/*
 * oaktree.h
 * ---------
 */

#ifndef __oaktree__
#define __oaktree__

struct halfplane
{
  REAL p [3], n [3];

  short vcolor, scolor;
};

struct sphere
{
  REAL c [3], r;

  short vcolor, scolor;
};

#if 0
struct  warp /* deformation as a space warp */
{
};
#endif

struct shape
{
  REAL *extents;

  enum {ADD, MUL, SUB, HPL, SPH} what;

  void *data;

  char *label;

  struct shape *left, *right;

  struct shape *prev, *next;
};

/* copy and label shape */
struct shape* shape_copy (struct shape *shape, char *label);

/* combine two shapes */
struct shape* shape_combine (struct shape *left, short what, struct shape *right);

/* return distance to shape at given point */
REAL shape_evaluate (struct shape *shape, REAL *point);

/* rotate shape about a point using a rotation matrix */
void shape_rotate (struct shape *shape, REAL *point, REAL *matrix);

/* free shape memory */
void shape_destroy (struct shape *shape);

struct octcut
{
  REAL d [8];

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

  struct shape *solids;

  struct octree *octree;

  struct simulation *prev, *next;
};

/* global simulations list */
extern struct simulation *simulation;

#endif
