/*
 * oaktree.c
 * ---------
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include "oaktree.h"
#include "viewer.h"
#include "render.h"
#include "input.h"
#include "timer.h"
#include "error.h"
#include "alg.h"

/* global simulations list */
struct simulation *simulation = NULL;

#if OPENGL
/* global viewer data */
static int vieweron = 0,
	   width = 512,
	   height = 512; 

/* menu callback */
static int  menu (char ***names, int **codes)
{
  return 0;
}

/* init callback */
static void init ()
{
  if (simulation)
  {
    viewer_update_extents (simulation->extents);
  }
}

/* idle callback */
static int  idle ()
{
  return 0;
}

/* quit callback */
static void quit ()
{
}

/* render callback */
static void render ()
{
  if (simulation)
  {
    render_shapes (simulation->octree, simulation->cutoff);

    render_octree (simulation->octree);
  }
}

/* key callback */
static void key (int key, int x, int y)
{
}

/* keyspec callback */
static void keyspec (int key, int x, int y)
{
}

/* mouse callback */
static void mouse (int button, int state, int x, int y)
{
}

/* motion callback */
static void motion (int x, int y)
{
}

/* passive callback */
static void passive (int x, int y)
{
}
#endif

/* initialize simulation */
static void initialize (struct simulation *simulation)
{
  struct solid *solid;
  REAL e [6], g [6];
  struct timing t;
  double dt;

  timerstart (&t);

  g [0] =  FLT_MAX;
  g [1] =  FLT_MAX;
  g [2] =  FLT_MAX;
  g [3] = -FLT_MAX;
  g [4] = -FLT_MAX;
  g [5] = -FLT_MAX;

  for (solid = simulation->solid; solid; solid = solid->next)
  {
    shape_extents (solid->shape, e);

    if (e[0] < g[0]) g[0] = e[0];
    if (e[1] < g[1]) g[1] = e[1];
    if (e[2] < g[2]) g[2] = e[2];
    if (e[3] > g[3]) g[3] = e[3];
    if (e[4] > g[4]) g[4] = e[4];
    if (e[5] > g[5]) g[5] = e[5];
  }

  SUB (g+3, g, e);
  MAXABS (e, e[3]);
  MID (g+3, g, e);
  e[3] = 0.5*e[3];
  g [0] = e[0] - e[3];
  g [1] = e[1] - e[3];
  g [2] = e[2] - e[3];
  g [3] = e[0] + e[3];
  g [4] = e[1] + e[3];
  g [5] = e[2] + e[3];
  COPY6 (g, simulation->extents); /* centered cube */

  simulation->octree = octree_create (simulation->extents);

  for (solid = simulation->solid; solid; solid = solid->next)
    octree_insert_shape (simulation->octree, solid->shape, simulation->cutoff);

  dt = timerend (&t);

  printf ("Simulation initialized in %g s.\n", dt);
}

/* run simulation */
static void run (struct simulation *simulation)
{
}

/* finalize simulation */
static void finalize (struct simulation *simulation)
{
}

/* return input file path and parse arguments */
static char* getfile (int argc, char **argv)
{
  char *path;
  FILE *f;
  int n;

  for (n = 1, path = NULL; n < argc; n ++)
  {
    if (path == NULL && (f = fopen (argv [n], "r")))
    {
      path = argv [n];
      fclose (f);
    }
#if OPENGL
    else if (strcmp (argv [n], "-v") == 0) vieweron = 1;
    else if (strcmp (argv [n], "-g") == 0)
    {
      if (++ n < argc)
      {
	sscanf (argv [n], "%dx%d", &width, &height);
      }
    }
#endif
  }

  return path;
}

/* entry point */
int main (int argc, char **argv)
{
  struct simulation *s, *n;
  int inputerror;

#if OPENGL
  char *synopsis = "SYNOPSIS: oaktree [-v] [-g WIDTHxHEIGHT] path\n";
#else
  char *synopsis = "SYNOPSIS: oaktree path\n";
#endif
  char *path = getfile (argc, argv);

  if (!path) printf (synopsis);
  else inputerror = input (path);

  if (!inputerror)
  {
    for (s = simulation; s; s = s->next)
    {
      initialize (s);
    }
  }

#if OPENGL
  if (vieweron && !inputerror)
  {
    REAL extents [6] = {-1, -1, -1, 1, 1, 1};

    viewer (&argc, argv, "oeaktree", width, height, extents, menu,
      init, idle, quit, render, key, keyspec, mouse, motion, passive);
  }
  else
#endif
  {
    for (s = simulation; s; s = s->next)
    {
      run (s);
    }
  }

  for (s = simulation; s; s = n)
  {
    n = s->next;
    finalize (s);
    free (s);
  }

  return 0;
}
