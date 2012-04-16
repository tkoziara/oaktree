/*
 * err.h
 * -----
 */

#ifndef __error__
#define __error__

/* textual assertion */
#define ASSERT(__test__, ...)\
  do {\
  if (!(__test__)) { fprintf (stderr, "%s: %d => ", __FILE__, __LINE__);\
    fprintf (stderr, __VA_ARGS__);\
    fprintf (stderr, "\n"); exit (1); } } while (0)

/* memory validity assertion */
#define ERRMEM(__pointer__) ASSERT (__pointer__, "Out of memory!");

#endif
