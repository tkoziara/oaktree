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
void render_octree (struct octree *oct)
{
  REAL *e;
  int i;

  if (oct->down [0])
  {
    for (i = 0; i < 8; i ++) render_octree (oct->down [i]);
  }
  else
  {
    glBegin (GL_LINES);

    glColor3f (0.2, 0.2, 0.2);

    e = oct->extents;

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

/* render shapes */
void render_shapes (struct octree *oct, REAL cutoff)
{
  struct triang *triang;
  int i;

  if (oct->down [0])
  {
    for (i = 0; i < 8; i ++) render_shapes (oct->down [i], cutoff);
  }

  glBegin (GL_TRIANGLES);

  glColor3f (0.5, 0.5, 0.5);

  for (triang = oct->triang; triang; triang = triang->next)
  {
    REAL (*t) [4][3] = triang->t;

    for (i = 0; i < triang->n; i ++)
    {
      glNormal3f (t[i][3][0], t[i][3][1], t[i][3][2]);
      glVertex3f (t[i][0][0], t[i][0][1], t[i][0][2]);
      glVertex3f (t[i][1][0], t[i][1][1], t[i][1][2]);
      glVertex3f (t[i][2][0], t[i][2][1], t[i][2][2]);
    }
  }

  glEnd ();
}
