/*
 * oaktree.h
 * ---------
 */

#ifndef __oaktree__
#define __oaktree__

struct superellipsoid
{
  REAL c [3], u [3], v [3], w [3], p, q; /* F(x) = (|<c-x,u>|**p + |<c-x,v>|**p)**(q/p) + |<c-x,w>|**q - 1 */

  short vcolor, scolor;
};

#if 0
struct warp /* space warp solely defines physical deformation */
{
};
#endif

struct shape /* general shape */
{
  REAL extents [6];

  enum {ADD, MUL, SUB, ELP} what;

  void *data;

  struct shape *left, *right;

  struct shape *next;
};

/* return distance to shape at given point */
REAL shape_evaluate (struct shape *shape, REAL *point);

struct contact
{
  REAL p [3], n [3];

  short color [2];

  struct shape *shape [2];

  struct contact *next;
};

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

  struct contact *con;

  struct octree *up, *down [8];
};

/* create octree down to a cutoff edge length */
struct octree* octree_create (REAL extents [6], REAL cutoff);

/* insert list of shape and refine octree down to a cutoff edge length */
void octree_insert_shapes (struct octree *oct, struct shape *list, REAL cutoff);

/* extract contact points and coarsen them up to a cutoff distance */
struct contact* octree_extract_contacts (struct octree *oct, REAL cutoff);

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
