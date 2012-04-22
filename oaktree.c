/*
 * oaktree.c
 * ---------
 */

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "oaktree.h"
#include "viewer.h"
#include "render.h"
#include "input.h"
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
    simulation->octree = octree_create (simulation->extents, simulation->grid);

    octree_insert_shapes (simulation->octree, simulation->solids, simulation->cutoff);

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
    //render_octree (simulation->octree);
   
    render_shapes (simulation->octree, simulation->cutoff);
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

/* run analysis */
static void runanalysis ()
{
  if (simulation)
  {
    simulation->octree = octree_create (simulation->extents, simulation->grid);

    octree_insert_shapes (simulation->octree, simulation->solids, simulation->cutoff);
  }
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
  int inputerror;

#if OPENGL
  char *synopsis = "SYNOPSIS: oaktree [-v] [-g WIDTHxHEIGHT] path\n";
#else
  char *synopsis = "SYNOPSIS: oaktree path\n";
#endif
  char *path = getfile (argc, argv);

  if (!path) printf (synopsis);
  else inputerror = input (path);

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
      runanalysis ();
    }

  return 0;
}
