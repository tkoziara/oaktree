/*
 * sort.h
 * ------
 */

#include <limits.h>

#ifndef __sort__
#define __sort__

#define DOUBLY_LINKED(PREV, NEXT) for (i = list; i; j = i, i = i->NEXT) i->PREV = j;

#define SINGLE_LINKED(PREV, NEXT)

#define IMPLEMENT_LIST_SORT(KIND, CALL, LIST, PREV, NEXT, LE)\
static LIST* CALL (LIST *list)\
{\
  LIST *i, *j, *k, *h, *t;\
  int l, m, n;\
\
  for (l = 1; l < INT_MAX;l *= 2)\
  {\
    h = t = NULL;\
\
    for (j = list;;)\
    {\
      i = j;\
\
      for (m = 0; m < l && j; j = j->NEXT, m ++);\
      for (n = 0, k = j; n < l && k; k = k->NEXT, n ++);\
\
      if (!j && i == list)\
      {\
	KIND (PREV, NEXT)\
	return list;\
      }\
      else if (!(m+n)) break;\
\
      if (!h) h = (LE (i, j) ? i : j);\
\
      for (; m && n;)\
      {\
	if (LE (i, j))\
	{\
	  if (t) t->NEXT = i;\
	  t = i;\
	  i = i->NEXT;\
	  m --;\
	}\
	else\
	{\
	  if (t) t->NEXT = j;\
	  t = j;\
	  j = j->NEXT;\
	  n --;\
	}\
      }\
\
      while (m)\
      {\
	t->NEXT = i;\
	t = i;\
	i = i->NEXT;\
	m --;\
      }\
\
      while (n)\
      {\
	t->NEXT = j;\
	t = j;\
	j = j->NEXT;\
	n --;\
      }\
    }\
\
    t->NEXT = NULL;\
    list = h;\
  }\
\
  return list;\
}

#endif
