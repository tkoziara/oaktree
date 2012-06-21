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
#if __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
  #include <GL/glext.h>
#endif

static int vieweron = 0, width = 512, height = 512;  /* viewer flag, initial window width and height */
enum {MENU_SIMULATION = 0, MENU_RENDER, MENU_LAST}; /* menu identifiers */
static char* menu_name [MENU_LAST];  /* menu names */
static int menu_code [MENU_LAST]; /* menu codes */
enum {SIMULATION_NEXT, SIMULATION_PREVIOUS, RENDER_DOMAINS, RENDER_CELLS, RENDER_OCTREE}; /* menu items */
static enum {DOMAINS = 1 << 0, CELLS = 1 << 1, OCTREE = 1 << 2} render_item = DOMAINS; /* render item */

/* simulation menu callback */
static void menu_simulation (int item)
{
  switch (item)
  {
  case SIMULATION_NEXT:
    if (simulation->next) simulation = simulation->next;
    viewer_update_extents (simulation->extents);
    break;
  case SIMULATION_PREVIOUS:
    if (simulation->prev) simulation = simulation->prev;
    viewer_update_extents (simulation->extents);
    break;
  }

  viewer_redraw_all ();
}

/* render menu callback */
static void menu_render (int item)
{
  switch (item)
  {
  case RENDER_DOMAINS:
    if (render_item & DOMAINS) render_item &= ~DOMAINS;
    else render_item |= DOMAINS;
    break;
  case RENDER_CELLS:
    if (render_item & CELLS) render_item &= ~CELLS;
    else render_item |= CELLS;
    break;
  case RENDER_OCTREE:
    if (render_item & OCTREE) render_item &= ~OCTREE;
    else render_item |= OCTREE;
    break;
  }

  viewer_redraw_all ();
}

/* menu callback */
static int  menu (char ***names, int **codes)
{
  ASSERT (simulation, "No simulation defined!");

  menu_name [MENU_SIMULATION] = "simulation";
  menu_code [MENU_SIMULATION] = glutCreateMenu (menu_simulation);
  glutAddMenuEntry ("previous /</", SIMULATION_PREVIOUS);
  glutAddMenuEntry ("next />/", SIMULATION_NEXT);

  menu_name [MENU_RENDER] = "render";
  menu_code [MENU_RENDER] = glutCreateMenu (menu_render);
  glutAddMenuEntry ("domains /d/", RENDER_DOMAINS);
  glutAddMenuEntry ("cells /c/", RENDER_CELLS);
  glutAddMenuEntry ("octree /o/", RENDER_OCTREE);

  *names = menu_name;
  *codes = menu_code;
  return MENU_LAST;
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
    if (render_item & DOMAINS) render_domains (simulation->octree);

    if (render_item & CELLS)
    {
      if (render_item & DOMAINS)
      {
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      }

      render_cells (simulation->octree);

      if (render_item & DOMAINS) glDisable (GL_BLEND);
    }

    if (render_item & OCTREE)
    {
      if (render_item & (DOMAINS|CELLS))
      {
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      }

      render_octree (simulation->octree);

      if (render_item & (DOMAINS|CELLS)) glDisable (GL_BLEND);
    }
  }
}

/* key callback */
static void key (int key, int x, int y)
{
  switch (key)
  {
  case 27:
    break;
  case '<':
    menu_simulation (SIMULATION_PREVIOUS);
    break;
  case '>':
    menu_simulation (SIMULATION_NEXT);
    break;
  case 'd':
    menu_render (RENDER_DOMAINS);
    break;
  case 'c':
    menu_render (RENDER_CELLS);
    break;
  case 'o':
    menu_render (RENDER_OCTREE);
    break;
  }
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
  struct domain *domain;
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

  for (domain = simulation->domain; domain; domain = domain->next)
  {
    shape_extents (domain->shape, e);

    if (e[0] < g[0]) g[0] = e[0];
    if (e[1] < g[1]) g[1] = e[1];
    if (e[2] < g[2]) g[2] = e[2];
    if (e[3] > g[3]) g[3] = e[3];
    if (e[4] > g[4]) g[4] = e[4];
    if (e[5] > g[5]) g[5] = e[5];
  }

  g [0] -= simulation->cutoff;
  g [1] -= simulation->cutoff;
  g [2] -= simulation->cutoff;
  g [3] += simulation->cutoff;
  g [4] += simulation->cutoff;
  g [5] += simulation->cutoff;

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

  for (domain = simulation->domain; domain; domain = domain->next)
  {
    octree_insert_domain (simulation->octree, domain, simulation->cutoff);
  }

  dt = timerend (&t);

  printf ("Simulation [%s] initialized in %g s.\n", simulation->outpath, dt);
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
