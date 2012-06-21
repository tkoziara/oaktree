/*
 * oaktree.h
 * ---------
 */

#ifndef __oaktree__
#define __oaktree__

struct halfspace
{
  REAL p [3], n [3], r, s;

  short scolor;
};

struct sphere
{
  REAL c [3], r, s;

  short scolor;
};

struct cylinder
{
  REAL p [3], d [3], r, s;

  short scolor;
};

struct fillet
{
  REAL r;

  short scolor;
};

struct shape
{
  enum {ADD, MUL, HSP, SPH, CYL, FLT} what;

  void *data;

  struct shape *up, *left, *right;
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

/* insert fillet between surfaces overlapping (c, r) sphere */
void shape_fillet (struct shape *shape, REAL c [3], REAL r, REAL fillet, short scolor);

/* return distance to shape at given point */
REAL shape_evaluate (struct shape *shape, REAL *point);

/* compute shape extents */
void shape_extents (struct shape *shape, REAL *extents);

/* output unique shape leaves overlapping (c,r) sphere and return their count or inside flag if count is zero */
int shape_unique_leaves (struct shape *shape, REAL c [3], REAL r, struct shape ***leaves, char *inside);

/* test whether the leaf is in a union of shapes */
int shape_leaf_in_union (struct shape *leaf);

/* free shape memory */
void shape_destroy (struct shape *shape);

struct domain
{
  struct shape *shape;

  char *label;

  REAL grid;

  struct domain *prev, *next;
};

struct triang
{
  REAL (*t) [4][3];

  short n;
};

struct cell
{
  struct triang *triang;

  struct domain *domain;

  struct cell **adj;

  short nadj;

  struct cell *next;
};

struct octree
{
  REAL extents [6];

  struct cell *cell;

  struct octree *up, *down [8];
};

/* create octree */
struct octree* octree_create (REAL extents [6]);

/* insert domain and refine octree down to a cutoff edge length */
void octree_insert_domain (struct octree *octree, struct domain *domain, REAL cutoff);

/* free octree memory */
void octree_destroy (struct octree *octree);

struct simulation
{
  char *outpath;

  REAL duration;

  REAL step;

  REAL cutoff;

  REAL extents [6];

  struct domain *domain;

  struct octree *octree;

  struct simulation *prev, *next;
};

/* global simulations list */
extern struct simulation *simulation;

#endif
