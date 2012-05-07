/*
 * oaktree.h
 * ---------
 */

#ifndef __oaktree__
#define __oaktree__

struct halfplane
{
  REAL p [3], n [3], r, s;

  short vcolor, scolor;
};

struct sphere
{
  REAL c [3], r, s;

  short vcolor, scolor;
};

struct cylinder
{
  REAL p [3], d [3], r, s;

  short vcolor, scolor;
};

struct shape
{
  enum {ADD, MUL, HPL, SPH, CYL} what;

  void *data;

  struct shape *left, *right;
};

/* copy shape */
struct shape* shape_copy (struct shape *shape);

/* return the same inverted shape */
struct shape* shape_invert (struct shape *shape);

/* combine two shapes */
struct shape* shape_combine (struct shape *left, short what, struct shape *right);

/* move shape */
void shape_move (struct shape *shape, REAL *vector);

/* rotate shape about a point using a rotation matrix */
void shape_rotate (struct shape *shape, REAL *point, REAL *matrix);

/* return distance to shape at given point */
REAL shape_evaluate (struct shape *shape, REAL *point);

/* return value and compute shape normal at given point */
REAL shape_normal (struct shape *shape, REAL *point, REAL *normal);

/* output unique shape leaves overlapping (c,r) sphere and return their count */
int shape_unique_leaves (struct shape *shape, REAL c [3], REAL r, struct shape ***leaves);

/* free shape memory */
void shape_destroy (struct shape *shape);

struct triang
{
  REAL (*t) [4][3];

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

/* create octree */
struct octree* octree_create (REAL extents [6]);

/* insert shape and refine octree down to a cutoff edge length */
void octree_insert_shape (struct octree *octree, struct shape *shape, REAL cutoff);

/* free octree memory */
void octree_destroy (struct octree *octree);

struct solid
{
  struct shape *shape;

  char *label;

  struct solid *prev, *next;
};

struct simulation
{
  char *outpath;

  REAL duration, step;

  REAL cutoff;

  REAL extents [6];

  struct solid *solid;

  struct octree *octree;

  struct simulation *prev, *next;
};

/* global simulations list */
extern struct simulation *simulation;

#endif
