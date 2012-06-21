/*
 * render.h
 * --------
 */

#include "oaktree.h"

#ifndef __render__
#define __render__

/* render octree itself */
void render_octree (struct octree *octree);

/* render domains */
void render_domains (struct octree *octree);

/* render cells */
void render_cells (struct octree *octree);

#endif
