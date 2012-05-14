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

/* render solids */
void render_solids (struct octree *octree)
{
  struct element *element;
  int i;

  if (octree->down [0])
  {
    for (i = 0; i < 8; i ++) render_solids (octree->down [i]);
  }

  glBegin (GL_TRIANGLES);

  glColor3f (0.5, 0.5, 0.5);

  for (element = octree->element; element; element = element->next)
  {
    struct triang *triang = element->triang;

    if (triang)
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
  }

  glEnd ();
}

/* render elements */
void render_elements (struct octree *octree)
{
  REAL e [6], w, *x;
  int i;

  if (octree->down [0])
  {
    for (i = 0; i < 8; i ++) render_elements (octree->down [i]);
  }

  if (octree->element)
  {
    glBegin (GL_QUADS);

    glColor4f (0.6, 0.6, 0.6, 0.3);

    x = octree->extents;
    w = x[3] - x[0];
    e[0] = x[0] + 0.07*w;
    e[1] = x[1] + 0.07*w;
    e[2] = x[2] + 0.07*w;
    e[3] = x[0] + 0.93*w;
    e[4] = x[1] + 0.93*w;
    e[5] = x[2] + 0.93*w;

    /* lower base */
    glNormal3f (0, 0, -1);
    glVertex3f (e[0], e[1], e[2]);
    glVertex3f (e[0], e[4], e[2]);
    glVertex3f (e[3], e[4], e[2]);
    glVertex3f (e[3], e[1], e[2]);

    /* upper base */
    glNormal3f (0, 0, 1);
    glVertex3f (e[0], e[1], e[5]);
    glVertex3f (e[3], e[1], e[5]);
    glVertex3f (e[3], e[4], e[5]);
    glVertex3f (e[0], e[4], e[5]);

    /* y low */
    glNormal3f (0, -1, 0);
    glVertex3f (e[0], e[1], e[2]);
    glVertex3f (e[3], e[1], e[2]);
    glVertex3f (e[3], e[1], e[5]);
    glVertex3f (e[0], e[1], e[5]);

    /* y high */
    glNormal3f (0, 1, 0);
    glVertex3f (e[0], e[4], e[2]);
    glVertex3f (e[0], e[4], e[5]);
    glVertex3f (e[3], e[4], e[5]);
    glVertex3f (e[3], e[4], e[2]);

    /* x low */
    glNormal3f (-1, 0, 0);
    glVertex3f (e[0], e[1], e[2]);
    glVertex3f (e[0], e[1], e[5]);
    glVertex3f (e[0], e[4], e[5]);
    glVertex3f (e[0], e[4], e[2]);

    /* x high */
    glNormal3f (1, 0, 0);
    glVertex3f (e[3], e[1], e[2]);
    glVertex3f (e[3], e[4], e[2]);
    glVertex3f (e[3], e[4], e[5]);
    glVertex3f (e[3], e[1], e[5]);

    glEnd ();
  }
}

/* render nodes */
void render_nodes (struct octree *octree)
{
  struct node **node;
  REAL p [8][3], *x;
  int i;

  if (octree->down [0])
  {
    for (i = 0; i < 8; i ++) render_nodes (octree->down [i]);
  }

  if (octree->element)
  {
    node = octree->element->node;
    x = octree->extents;

    VECTOR (p[0], x[0], x[1], x[2]);
    VECTOR (p[1], x[0], x[4], x[2]);
    VECTOR (p[2], x[3], x[4], x[2]);
    VECTOR (p[3], x[3], x[1], x[2]);
    VECTOR (p[4], x[0], x[1], x[5]);
    VECTOR (p[5], x[0], x[4], x[5]);
    VECTOR (p[6], x[3], x[4], x[5]);
    VECTOR (p[7], x[3], x[1], x[5]);

    glPointSize (4.0);
    glBegin (GL_POINTS);

    for (i = 0; i < 8; i ++)
    {
      if (node [i]->hanging)
        glColor3f (1, 0, 0);
      else glColor3f (0, 0, 1);

      glVertex3f (p[i][0], p[i][1], p[i][2]);
    }

    glEnd ();
    glPointSize (1.0);
  }
}
