/*
 * render.c
 * --------
 */

#include <stdio.h>
#if __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
  #include <GL/glext.h>
#endif
#include "render.h"
#include "error.h"
#include "alg.h"

/* render octree itself */
void render_octree (struct octree *octree)
{
  REAL *e;
  int i;

  if (octree->down [0])
  {
    for (i = 0; i < 8; i ++) render_octree (octree->down [i]);
  }
  else
  {
    glBegin (GL_LINES);

    glColor4f (0.1, 0.1, 0.1, 0.1);

    e = octree->extents;

    /* lower base */
    glVertex3f (e[0], e[1], e[2]);
    glVertex3f (e[3], e[1], e[2]);

    glVertex3f (e[0], e[1], e[2]);
    glVertex3f (e[0], e[4], e[2]);

    glVertex3f (e[3], e[1], e[2]);
    glVertex3f (e[3], e[4], e[2]);

    glVertex3f (e[0], e[4], e[2]);
    glVertex3f (e[3], e[4], e[2]);

    /* upper base */
    glVertex3f (e[0], e[1], e[5]);
    glVertex3f (e[3], e[1], e[5]);

    glVertex3f (e[0], e[1], e[5]);
    glVertex3f (e[0], e[4], e[5]);

    glVertex3f (e[3], e[1], e[5]);
    glVertex3f (e[3], e[4], e[5]);

    glVertex3f (e[0], e[4], e[5]);
    glVertex3f (e[3], e[4], e[5]);

    /* sides */
    glVertex3f (e[0], e[1], e[2]);
    glVertex3f (e[0], e[1], e[5]);

    glVertex3f (e[3], e[1], e[2]);
    glVertex3f (e[3], e[1], e[5]);

    glVertex3f (e[3], e[4], e[2]);
    glVertex3f (e[3], e[4], e[5]);

    glVertex3f (e[0], e[4], e[2]);
    glVertex3f (e[0], e[4], e[5]);

    glEnd ();
  }
}

/* render domains */
void render_domains (struct octree *octree)
{
  struct cell *cell;
  struct face *face;
  int i;

  if (octree->down [0])
  {
    for (i = 0; i < 8; i ++) render_domains (octree->down [i]);
  }

  glBegin (GL_TRIANGLES);

  glColor3f (0.5, 0.5, 0.5);

  for (cell = octree->cell; cell; cell = cell->next)
  {
    for (face = cell->face; face; face = face->next)
    {
      if (face->leaf == NULL) continue;

      REAL (*t) [3][3] = face->t;

      for (i = 0; i < face->n; i ++)
      {
	glNormal3f (face->normal[0], face->normal[1], face->normal[2]);
	glVertex3f (t[i][0][0], t[i][0][1], t[i][0][2]);
	glVertex3f (t[i][1][0], t[i][1][1], t[i][1][2]);
	glVertex3f (t[i][2][0], t[i][2][1], t[i][2][2]);
      }
    }
  }

  glEnd ();
}

/* render cells */
void render_cells (struct octree *octree)
{
  REAL p [3], q [3], *x;
  struct cell *cell;
  struct face *face;
  int i, j;

  if (octree->down [0])
  {
    for (i = 0; i < 8; i ++) render_cells (octree->down [i]);
  }

  x = octree->extents;

  MID (x, x+3, p);

  for (cell = octree->cell; cell; cell = cell->next)
  {
    glBegin (GL_TRIANGLES);

    glColor4f (0.6, 0.6, 0.6, 0.3);

    for (face = cell->face; face; face = face->next)
    {
      glNormal3f (face->normal[0], face->normal[1], face->normal[2]);

      for (i = 0; i < face->n; i++)
      {
	for (j = 0; j < 3; j ++)
	{
	  SUB (face->t[i][j], p, q);
	  ADDMUL (p, 0.8, q, q);
	  glVertex3f (q[0], q[1], q[2]);
	}
      }
    }

    glEnd ();
  }
}
