/*
 * oaktree.h
 * ---------
 */

#ifndef __oaktree__
#define __oaktree__

struct halfplane
{
  REAL oobb [8][3];

  REAL p [3], n [3];

  short vcolor, scolor;
};

struct sphere
{
  REAL oobb [8][3];

  REAL c [3], r;

  short vcolor, scolor;
};

struct cylinder
{
  REAL oobb [8][3];

  REAL p [3], d [3], r;

  short vcolor, scolor;
};

#if 0
struct  warp /* deformation as a space warp */
{
};
#endif

struct shape
{
  enum {ADD, MUL, SUB, HPL, SPH, CYL} what;

  void *data;

  char *label;

  struct shape *left, *right;

  struct shape *prev, *next;
};

/* copy and label shape */
struct shape* shape_copy (struct shape *shape, char *label);

/* combine two shapes */
struct shape* shape_combine (struct shape *left, short what, struct shape *right);

/* output unique shape leaves inside of a grid cell and return their count */
int shape_unique_leaves (struct shape *shape, REAL p [8][3], REAL cutoff, struct shape ***leaves);

/* rotate shape about a point using a rotation matrix */
void shape_rotate (struct shape *shape, REAL *point, REAL *matrix);

/* return distance to shape at given point */
REAL shape_evaluate (struct shape *shape, REAL *point);

/* return value and compute shape normal at given point */
REAL shape_normal (struct shape *shape, REAL *point, REAL *normal);

/* free shape memory */
void shape_destroy (struct shape *shape);

struct triang
{
  REAL (*t) [3][3];

  short n;

  struct shape *shape;

  struct triang *next;
};

struct octree
{
  REAL extents [6];

  struct triang *triang;

  struct octree *up, *down [8];
};

/* create octree down to a cutoff edge length */
struct octree* octree_create (REAL extents [6], REAL cutoff);

/* insert shape and refine octree down to a cutoff edge length */
void octree_insert_shape (struct octree *oct, struct shape *shape, REAL cutoff);

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
