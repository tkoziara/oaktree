/*
 * render.h
 * --------
 */

#include "oaktree.h"

#ifndef __render__
#define __render__

/* render octree itself */
void render_octree (struct octree *oct);

/* render shapes */
void render_shapes (struct octree *oct, REAL cutoff);

/* render contacts */
void render_contacts (struct octree *oct);

#endif
