#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "error.h"

/* STL read */
REAL* stlread (const char *path, int *count)
{
  FILE *f = fopen (path, "r");
  if (!f) return NULL;
  int size = 1024, i;
  char buff[1024];
  REAL *triang, *t;
  (*count) = 0;

  ERRMEM (triang = malloc (size * sizeof(REAL[9])));
  t = triang;

  while (fscanf (f, "%s", buff) != EOF)
  {
    if (strcmp (buff, "outer") == 0)
    {
      fscanf (f, "%s", buff);
      ASSERT (strcmp (buff, "loop") == 0, "STL fromat error");
      for (i = 0; i < 3; i ++)
      {
	fscanf (f, "%s", buff);
        ASSERT (strcmp (buff, "vertex") == 0, "STL fromat error");
	fscanf (f, "%f", t+3*i+0); /* FIXME: detect REAL type */
	fscanf (f, "%f", t+3*i+1);
	fscanf (f, "%f", t+3*i+2);
      }
      t += 9;
      (*count) ++;

      if ((*count) == size)
      {
	size *= 2;
	ERRMEM (triang = realloc (triang, size * sizeof(REAL[9])));
	t = &triang[(*count)*9];
      }
    }
  }

  return triang;
}
