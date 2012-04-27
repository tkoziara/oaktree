/*
 * timer.h
 * ---------
 */

#include <sys/time.h>

#ifndef __timer__
#define __timer__

struct timing
{
  struct timeval time;
  double sec, total;
};

static inline void timerstart (struct timing *t)
{
  gettimeofday (&t->time, NULL);
}

static inline double timerend (struct timing *t)
{
  struct timeval newtime;
  gettimeofday (&newtime, NULL);
  t->sec = (((double)newtime.tv_sec - (double)t->time.tv_sec) * 1000000. +
    ((double)newtime.tv_usec - (double)t->time.tv_usec)) / 1000000.;
  t->total += t->sec;
  return t->sec;
}

#endif
