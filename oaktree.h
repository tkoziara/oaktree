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

struct shape
{
  enum {ADD, MUL, HSP, SPH, CYL} what;

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

/* output unique shape leaves overlapping (c,r) sphere and return their count or inside flag if count is zero */
int shape_unique_leaves (struct shape *shape, REAL c [3], REAL r, struct shape ***leaves, int *inside);

/* compute shape extents */
void shape_extents (struct shape *shape, REAL *extents);

/* free shape memory */
void shape_destroy (struct shape *shape);

struct solid
{
  struct shape *shape;

  char *label;

  struct solid *prev, *next;
};

struct triang
{
  REAL (*t) [4][3];

  short n;
};

struct node
{
  unsigned char hanging;

  struct node *next;
};

struct element
{
  struct node *node [8];

  struct triang *triang;

  struct solid *solid;

  struct element *next;
};

struct octree
{
  REAL extents [6];

  struct element *element;

  struct octree *up, *down [8];
};

/* create octree */
struct octree* octree_create (REAL extents [6]);

/* insert solid and refine octree down to a cutoff edge length */
struct node* octree_insert_solid (struct octree *octree, struct solid *solid, REAL cutoff);

/* free octree memory */
void octree_destroy (struct octree *octree);

struct simulation
{
  char *outpath;

  REAL duration;

  REAL step;

  REAL cutoff;

  REAL extents [6];

  struct solid *solid;

  struct octree *octree;

  struct node *node;

  struct simulation *prev, *next;
};

/* global simulations list */
extern struct simulation *simulation;

#endif
