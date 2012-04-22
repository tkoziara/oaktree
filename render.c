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
#include "polygon.h"
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
  REAL p [8][3], t [5][3][3], *e = oct->extents;
  struct octcut *cut;
  int i, j;

  if (oct->down [0])
  {
    for (i = 0; i < 8; i ++) render_shapes (oct->down [i], cutoff);
  }

  VECTOR (p[0], e[0], e[1], e[2]);
  VECTOR (p[1], e[0], e[4], e[2]);
  VECTOR (p[2], e[3], e[4], e[2]);
  VECTOR (p[3], e[3], e[1], e[2]);
  VECTOR (p[4], e[0], e[1], e[5]);
  VECTOR (p[5], e[0], e[4], e[5]);
  VECTOR (p[6], e[3], e[4], e[5]);
  VECTOR (p[7], e[3], e[1], e[5]);

  glBegin (GL_TRIANGLES);

  glColor3f (0.5, 0.5, 0.5);

  for (cut = oct->cut; cut; cut = cut->next)
  {
    j = polygonise (p, cut->d, 0.0, 0.01*cutoff, t);

    for (i = 0; i < j; i ++)
    {
      REAL n [3];
      NORMAL (t[i][0], t[i][1], t[i][2], n); /* FIXME: let OpenGL do this */
      glNormal3f (n[0], n[1], n[2]);
      glVertex3f (t[i][0][0], t[i][0][1], t[i][0][2]);
      glVertex3f (t[i][1][0], t[i][1][1], t[i][1][2]);
      glVertex3f (t[i][2][0], t[i][2][1], t[i][2][2]);
    }
  }

  glEnd ();
}
