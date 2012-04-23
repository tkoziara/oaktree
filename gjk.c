/*
 * gjk.c
 * ------
 */

#include <float.h>
#include "gjk.h"
#include "alg.h"

/* auxiliary point structure used to describe
 * points in the Minkowski difference set (A - B) */
typedef struct { REAL w[3], *a, *b; } point;

/* simple bit lighting macro for vertex sets */
#define set(i,j,k,l) ((i*8)|(j*4)|(k*2)|(l*1))

/* number of bits lightened in numbers from 0 to 15;
 * they correspond to the cardinality of a vertex set
 * describing the inner simplex approximation of (A - B) */
static const int card [] =
{0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};

/* vertex schemes; the first three sets use only indices from {0, 1},
 * the first seven schemes use only indices {0, 1, 2}, and the complete
 * fifteen of them employs the complete set of indices {0, 1, 2, 3} */
static const int base [] =
{set(0,0,0,1),
 set(0,0,1,0),
 set(0,0,1,1),
 set(0,1,0,0),
 set(0,1,0,1),
 set(0,1,1,0),
 set(0,1,1,1),
 set(1,0,0,0),
 set(1,0,0,1),
 set(1,0,1,0),
 set(1,1,0,0),
 set(1,0,1,1),
 set(1,1,0,1),
 set(1,1,1,0),
 set(1,1,1,1)};

/* last [n] is the number of last base sets
 * to be checked during the projection on
 * a simplex with n vertices */
static const int last [] = 
{0, 0, 3, 7, 15};

/* find minimal point in set (c, n) along the direction of 'v' */
inline static REAL* minimal_support_point (REAL *c, int n, REAL *v)
{
  REAL dot, dotmin = DBL_MAX, *out = NULL;

  for (; n > 0; n --, c += 3)
  {
    dot = DOT (c, v);
    if (dot < dotmin) {dotmin = dot; out = c;}
  }

  return out;
}

/* find maximal point in set (c, n) along the direction of 'v' */
inline static REAL* maximal_support_point (REAL *c, int n, REAL *v)
{
  REAL dot, dotmax = -DBL_MAX, *out = NULL;

  for (; n > 0; n --, c += 3)
  {
    dot = DOT (c, v);
    if (dot > dotmax) {dotmax = dot; out = c;}
  }

  return out;
}

/* projection of the zero point (0, 0, 0) onto the convex
 * hull spanned by points (w, n); on exit set 'w' is overwritten
 * by the smallest simplex in (w, n) such that projection of (0, 0, 0) 
 * onto it is the same as the projection onto the original hull; the
 * new dimension of 'w' is returend; 'l' contains barycentric coordinates
 * of the projection point; 'v' contains the point itslef */
static int project (point *w, int n, REAL *l, REAL *v)
{
  REAL dot [4][4], delta [16][4], sum;
  int i, j, k, m, s, c, o, f;

  if (n == 1)
  {
    l[0] = 1.0;
    COPY (w[0].w, v);
    return 1;
  }

  /* calculate all necessary dot products;
   * this could be optimised by calculating
   * dots related only to the new point */
  for (i = 0; i < n; i ++)
  {
    for (j = i; j < n; j ++)
    {
      dot [i][j] = DOT (w[i].w,w[j].w);
      dot [j][i] = dot [i][j];
    }
  }

  /* loop and calculate components
   * of the determinant expansion
   * according to the recursive formula
   * from the reference paper */
  for (m = 0; m < last [n]; m ++)
  {
    s = base [m];
    c = ~s;
    f = 1;
    if (card [s] == 1)
    {
      for (i = 0; i < n; i ++)
      {
	if (s & (1<<i))
	{
	  delta [s][i] = 1.0;
	  k = i;
	  break;
	}
      }
      for (j = 0; j < n; j ++)
      {
	if (c & (1<<j))
	{
	  delta [s][j] = dot[i][k] - dot[i][j];
	  if (delta [s][j] > 0.0) f = 0; /* cancel the flag => this is not yet the right sub-simplex */
	}
      }
    }
    else
    {
      for (i = 0, k = -1; i < n; i ++)
      {
	if (s & (1<<i))
	{
	  o = s;
	  o &= ~(1<<i);
	  if (k < 0) k = i;
	  delta [s][i] =  delta[o][i];
	  if (delta [s][i] <= 0.0) f = 0; 
	}
      }
      for (j = 0; j < n; j ++)
      {
	if (c & (1<<j))
	{
	  delta [s][j] = 0.0;
	  for (i = 0; i < n; i ++)
	  {
	    if (s & (1<<i))
	    {
	      delta [s][j] += delta [s][i] * (dot [i][k] - dot [i][j]); 
	    }
	  }
	  if (delta [s][j] > 0.0) f = 0; 
	}
      }
    }

    /* the flag was not canceled, that is
     * the projection on the current feature
     * was succesful => refresh vertex table
     * and compute barycentric coordinates
     * of the projection point */
    if (f)
    {
      sum = 0.0;
      for (i = o = 0; i < n; i ++)
      {
	if (s & (1<<i))
	{
	  w[o] = w[i];
	  l[o] = delta [s][i];
          sum += l[o];
	  o ++;
	}
      }
      SET (v, 0.0);
      for(i = 0; i < o; i ++)
      {
	l[i] /= sum; /* current coordinate */
	ADDMUL (v, l[i], w[i].w, v); /* convex combination == projection point */
      }
      return o;
    }
  }

  return 0;
}

/* public driver routine => input two polytopes A = (a, na) and B = (b, nb); outputs
 * p in A and q in B such that d = |p - q| is minimal; the distance d is returned */
REAL gjk (REAL epsilon, REAL *a, int na, REAL *b, int nb, REAL *p, REAL *q)
{
  point w [4];
  REAL v [3],
	 vlen,
	 delta,
	 mi = 0.0,
	 l [4];

  int toofar = 1,
      n = 0,
      j = 0,
      k = (na+nb)*(na+nb);

  SUB (a, b, v);
  vlen = LEN (v);

  while (toofar && vlen > epsilon && n < 4 && j ++ < k) /* (#) see below */
  {
    w[n].a = minimal_support_point (a, na, v);
    w[n].b = maximal_support_point (b, nb, v);
    SUB (w[n].a, w[n].b, w[n].w);
    delta = DOT (v, w[n].w) / vlen;
    mi = MAX (mi, delta);
    toofar = (vlen - mi) > epsilon;
    if (toofar)
    {
      n = project (w, ++n, l, v); /* (#) n = 4 can happen due to roundoff */
      vlen = LEN (v);
    }
  }

  if (n)
  {
    SET (p, 0);
    SET (q, 0);
    for (n --; n >= 0; n --)
    {
      ADDMUL (p, l[n], w[n].a, p);
      ADDMUL (q, l[n], w[n].b, q);
    }
  }
  else /* the while loop was never entered */
  {
    COPY (a, p);
    COPY (b, q);
  }

  return vlen;
}
