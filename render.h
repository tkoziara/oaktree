/*
 * render.h
 * --------
 */

#include "oaktree.h"

#ifndef __render__
#define __render__

/* render octree itself */
void render_octree (struct octree *octree);

/* render solids */
void render_solids (struct octree *octree);

/* render elements */
void render_elements (struct octree *octree);

/* render nodes */
void render_nodes (struct octree *octree);

#endif
