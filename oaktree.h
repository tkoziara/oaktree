/*
 * oaktree.h
 * ---------
 */

#ifndef __oaktree__
#define __oaktree__

#define REAL float

struct superellipsoid
{
  REAL extents [6];

  REAL c [3], v [3][3], r, t; /* F(x) = (|<c-x,v[0]>|**r + |<c-x,v[1]>|**r)**(t/r) + |z/v[2]|**t - 1 */

  short color;
};

struct shape /* general shape */
{
  REAL extents [6];

  enum {ADD, MUL, SUB, ELP} what;

  void *data;

  struct shape *left, *right;

  struct shape *next;
};

/* return distance to shape at given point */
REAL shape_evaluate (struct shape *shp, REAL *point);

struct contact
{
  REAL p [3], n [3];

  short color [2];

  struct shape *shp [2];

  struct contact *next;
};

struct octcut
{
  REAL d [8];

  short color;

  struct shape *shp;

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

#endif
