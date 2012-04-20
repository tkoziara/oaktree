/*
 * octree.c
 * --------
 */

#include <stdlib.h>
#include "oaktree.h"
#include "error.h"
#include "alg.h"

/* create octree down to a cutoff edge length */
struct octree* octree_create (REAL extents [6], REAL cutoff)
{
  return NULL;
}

/* insert list of shape and refine octree down to a cutoff edge length */
void octree_insert_shapes (struct octree *oct, struct shape *list, REAL cutoff)
{
}

/* free octree memory */
void octree_destroy (struct octree *oct)
{
}
