/*
 * render.h
 * --------
 */

#include "oaktree.h"

#ifndef __render__
#define __render__

/* render octree itself */
void render_octree (struct octree *octree);

/* render shapes */
void render_shapes (struct octree *octree, REAL cutoff);

#endif
