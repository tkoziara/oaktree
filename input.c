/*
 * input.h
 * ---------
 */

#include <Python.h>
#include <structmember.h>
#include <float.h>
#include "oaktree.h"
#include "input.h"
#include "error.h"
#include "alg.h"

#ifndef Py_RETURN_FALSE
#define Py_RETURN_FALSE return Py_INCREF(Py_False), Py_False
#endif

#ifndef Py_RETURN_TRUE
#define Py_RETURN_TRUE return Py_INCREF(Py_True), Py_True
#endif

#ifndef Py_RETURN_NONE
#define Py_RETURN_NONE return Py_INCREF(Py_None), Py_None
#endif

/* string buffer length */
#define BUFLEN 512

/* minimal type initialization */
#define TYPEINIT(typedesc, type, name, flags, dealloc, new, methods, members, getset)\
memset (&(typedesc), 0, sizeof (PyTypeObject));\
(typedesc).tp_basicsize = sizeof (type);\
(typedesc).tp_name = name;\
(typedesc).tp_flags = flags;\
(typedesc).tp_dealloc = (destructor)dealloc;\
(typedesc).tp_new = new;\
(typedesc).tp_methods = methods;\
(typedesc).tp_members = members;\
(typedesc).tp_getset = getset

/*
 * utilities
 */

/* string test */
static int is_string (PyObject *obj, char *var)
{
  if (obj)
  {
    if (!PyUnicode_Check (obj))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a string", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* positive test */
static int is_positive (double num, char *var)
{
  if (num <= 0)
  {
    char buf [BUFLEN];
    sprintf (buf, "'%s' must be positive", var);
    PyErr_SetString (PyExc_ValueError, buf);
    return 0;
  }

  return 1;
}

/* tuple test */
static int is_tuple (PyObject *obj, char *var, int len)
{
  if (obj)
  {
    if (!PyTuple_Check (obj))
    {
      char buf [BUFLEN];
      snprintf (buf, BUFLEN, "'%s' must be a tuple", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    if (len > 0 && PyTuple_Size (obj) != len)
    {
      char buf [BUFLEN];
      snprintf (buf, BUFLEN, "'%s' must have %d elements", var, len);
      PyErr_SetString (PyExc_ValueError, buf);
      return 0;
    }
  }

  return 1;
}

/* list of tuples test => returns length of the list or zero */
static int is_list_of_tuples (PyObject *obj, char *var, int min_length, int tuple_length)
{
  if (obj)
  {
    if (!PyList_Check (obj))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a list object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    int i, j, n = PyList_Size (obj);

    if (n < min_length)
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must have at least %d items", var, min_length);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    for (i = 0; i < n; i ++)
    {
      PyObject *item = PyList_GetItem (obj, i);

      if (!PyTuple_Check (item))
      {
	char buf [BUFLEN];
	sprintf (buf, "'%s' must be a list of tuples: item %d is not a tuple", var, i);
	PyErr_SetString (PyExc_ValueError, buf);
	return 0;
      }

      j = PyTuple_Size (item);

      if (j != tuple_length)
      {
	char buf [BUFLEN];
	sprintf (buf, "'%s' list items must be tuples of length %d: item %d has length %d", var, tuple_length, i, j);
	PyErr_SetString (PyExc_ValueError, buf);
	return 0;
      }
    }

    return n;
  }

  return 1;
}

/* define keywords */
#define KEYWORDS(...) char *kwl [] = {__VA_ARGS__, NULL}

/* parse arguments with keywords */
#define PARSEKEYS(fmt, ...) if (!PyArg_ParseTupleAndKeywords (args, kwds, fmt, kwl, __VA_ARGS__)) return NULL

/* parse arguments without keywords */
#define PARSE(fmt, ...) if (!PyArg_ParseTuple (args, fmt, __VA_ARGS__)) return NULL

/* object types assertion */
#define TYPETEST(test) if(!(test)) return NULL

/* string argument if block comparison */
#define IFIS(obj, val) if (strcmp (PyUnicode_AsUTF8 (obj), val) == 0)
#define ELIF(obj, val) else if (strcmp (PyUnicode_AsUTF8 (obj), val) == 0)
#define ELSE else

/*
 * SIMULATION
 */

static PyTypeObject SIMULATION_TYPE;

typedef struct {
  PyObject_HEAD
  struct simulation *ptr;
} SIMULATION;

/* SIMULATION test */
static int is_simulation (SIMULATION *obj, char *var)
{
  if (obj)
  {
    if (!PyObject_IsInstance ((PyObject*)obj, (PyObject*)&SIMULATION_TYPE))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a SIMULATION", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* constructor */
static PyObject* SIMULATION_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("outpath", "duration", "step", "cutoff");
  double duration, step, cutoff;
  struct simulation *simu;
  PyObject *outpath;
  SIMULATION *self;

  self = (SIMULATION*)type->tp_alloc (type, 0);

  if (self)
  {
    PARSEKEYS ("Oddd", &outpath, &duration, &step, &cutoff);

    TYPETEST (is_string (outpath, kwl [0]) && is_positive (duration, kwl[1]) &&
	      is_positive (step, kwl[2]) && is_positive (cutoff, kwl[4]));

    ERRMEM (simu = calloc (1, sizeof (struct simulation)));
    ERRMEM (simu->outpath = malloc (strlen (PyUnicode_AsUTF8 (outpath)) + 1));
    strcpy (simu->outpath, PyUnicode_AsUTF8 (outpath));
    simu->duration = duration;
    simu->step = step;
    simu->cutoff = cutoff;

    self->ptr = simu;

    if (simulation) simulation->prev = simu;
    simu->next = simulation;
    simulation = simu; /* insert into global list */
  }

  return (PyObject*)self;
}

/* destructor */
static void SIMULATION_dealloc (SIMULATION *self)
{
  Py_TYPE(self)->tp_free((PyObject*)self);
}

/* SIMULATION methods */
static PyMethodDef SIMULATION_methods [] =
{ {NULL, NULL, 0, NULL} };

/* SIMULATION members */
static PyMemberDef SIMULATION_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* SIMULATION getset */
static PyGetSetDef SIMULATION_getset [] =
{ {NULL, 0, 0, NULL, NULL} };

/*
 * SHAPE
 */

static PyTypeObject SHAPE_TYPE;

typedef struct {
  PyObject_HEAD
  struct shape *ptr;
} SHAPE;

/* SHAPE test */
static int is_shape (SHAPE *obj, char *var)
{
  if (obj)
  {
    if (!PyObject_IsInstance ((PyObject*)obj, (PyObject*)&SHAPE_TYPE))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a SHAPE", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* constructor */
static PyObject* SHAPE_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyErr_SetString (PyExc_RuntimeError, "Cannot create unspecified shape!");
  return NULL;
}

/* destructor */
static void SHAPE_dealloc (SHAPE *self)
{
  shape_destroy (self->ptr);

  Py_TYPE(self)->tp_free((PyObject*)self);
}

/* SHAPE methods */
static PyMethodDef SHAPE_methods [] =
{ {NULL, NULL, 0, NULL} };

/* SHAPE members */
static PyMemberDef SHAPE_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* SHAPE getset */
static PyGetSetDef SHAPE_getset [] =
{ {NULL, 0, 0, NULL, NULL} };

/*
 * DOMAIN
 */

static PyTypeObject DOMAIN_TYPE;

#ifdef DOMAIN
#undef DOMAIN
#endif

typedef struct {
  PyObject_HEAD
  struct domain *ptr;
} DOMAIN;

#if 0
/* DOMAIN test */
static int is_domain (DOMAIN *obj, char *var)
{
  if (obj)
  {
    if (!PyObject_IsInstance ((PyObject*)obj, (PyObject*)&DOMAIN_TYPE))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a DOMAIN", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}
#endif

/* constructor */
static PyObject* DOMAIN_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("simu", "shape", "label", "grid");
  struct domain *domain;
  SIMULATION *simu;
  PyObject *label;
  SHAPE *shape;
  double grid;
  DOMAIN *self;

  self = (DOMAIN*)type->tp_alloc (type, 0);

  if (self)
  {
    grid = FLT_MAX;
    label = NULL;

    PARSEKEYS ("OO|Od", &simu, &shape, &label, &grid);

    TYPETEST (is_simulation (simu, kwl[0]) && is_shape (shape, kwl[1]) &&
	      is_string (label, kwl[2]) && is_positive (grid, kwl[3]));

    if (grid <= simu->ptr->cutoff)
    {
      PyErr_SetString (PyExc_RuntimeError, "Grid size must be larger than simulation cutoff!");
      return NULL;
    }

    ERRMEM (domain = calloc (1, sizeof (struct domain)));
    domain->shape = shape_copy (shape->ptr);
    if (label)
    {
      ERRMEM (domain->label = malloc (strlen (PyUnicode_AsUTF8 (label)) + 1));
      strcpy (domain->label, PyUnicode_AsUTF8 (label));
    }
    domain->grid = grid;

    if (simu->ptr->domain) simu->ptr->domain->prev = domain;
    domain->next = simu->ptr->domain;
    simu->ptr->domain = domain;
    domain->prev = NULL;

    self->ptr = domain;
  }

  return (PyObject*)self;
}

/* destructor */
static void DOMAIN_dealloc (DOMAIN *self)
{
  Py_TYPE(self)->tp_free((PyObject*)self);
}

/* DOMAIN methods */
static PyMethodDef DOMAIN_methods [] =
{ {NULL, NULL, 0, NULL} };

/* DOMAIN members */
static PyMemberDef DOMAIN_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* DOMAIN getset */
static PyGetSetDef DOMAIN_getset [] =
{ {NULL, 0, 0, NULL, NULL} };

/*
 * SUBROUTINES
 */

/* create sphere */
static PyObject* SPHERE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("center", "r", "scolor");
  struct sphere *sphere;
  PyObject *center;
  int scolor;
  double r;
  SHAPE *out;

  out = (SHAPE*)SHAPE_TYPE.tp_alloc (&SHAPE_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("Odi", &center, &r, &scolor);

    TYPETEST (is_tuple (center, kwl[0], 3) && is_positive (r, kwl[1]));

    ERRMEM (out->ptr = calloc (1, sizeof (struct shape)));
    ERRMEM (sphere = malloc (sizeof (struct sphere)));

    sphere->c [0] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (center, 0));
    sphere->c [1] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (center, 1));
    sphere->c [2] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (center, 2));
    sphere->r = r;
    sphere->s = 1.0;
    sphere->scolor = scolor;

    out->ptr->what = SPH;
    out->ptr->data = sphere;
  }

  return (PyObject*)out;
}

/* create cylinder */
static PyObject* CYLINDER (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("base", "h", "r", "scolor");
  struct shape *sa, *sb, *sc;
  struct halfspace *a, *b;
  PyObject *base, *scolor;
  struct cylinder *c;
  double r, h;
  REAL x [3];
  SHAPE *out;

  out = (SHAPE*)SHAPE_TYPE.tp_alloc (&SHAPE_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OddO", &base, &h, &r, &scolor);

    TYPETEST (is_tuple (base, kwl[0], 3) && is_positive (h, kwl[1]) && is_positive (r, kwl[2]) && is_tuple (scolor, kwl[4], 3));

    ERRMEM (sa = calloc (1, sizeof (struct shape)));
    ERRMEM (a = malloc (sizeof (struct halfspace)));

    a->p [0] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (base, 0));
    a->p [1] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (base, 1));
    a->p [2] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (base, 2));
    VECTOR (a->n, 0, 0, -1);
    a->scolor = PyLong_AsLong(PyTuple_GetItem(scolor, 0));
    a->r = r;
    a->s = 1.0;
    sa->what = HSP;
    sa->data = a;

    ERRMEM (sb = calloc (1, sizeof (struct shape)));
    ERRMEM (b = malloc (sizeof (struct halfspace)));

    VECTOR (b->p, a->p[0], a->p[1], a->p[2]+h);
    VECTOR (b->n, 0, 0, 1);
    b->scolor = PyLong_AsLong(PyTuple_GetItem(scolor, 2));
    b->r = r;
    b->s = 1.0;
    sb->what = HSP;
    sb->data = b;

    ERRMEM (sc = calloc (1, sizeof (struct shape)));
    ERRMEM (c = malloc (sizeof (struct cylinder)));

    COPY (a->p, c->p);
    VECTOR (c->d, 0, 0, 1);
    c->r = r;
    c->s = 1.0;
    c->scolor = PyLong_AsLong(PyTuple_GetItem(scolor, 1));
    MID (a->p, b->p, x);
    sc->what = CYL;
    sc->data = c;

    out->ptr = shape_combine (sc, MUL, shape_combine (sa, MUL, sb));
  }

  return (PyObject*)out;
}

/* create cube */
static PyObject* CUBE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("corner", "u", "v", "w", "scolor");
  struct shape *a, *b, *c, *d, *e, *f;
  PyObject *corner, *scolor;
  struct halfspace *h;
  double u, v, w;
  double p [3];
  SHAPE *out;

  out = (SHAPE*)SHAPE_TYPE.tp_alloc (&SHAPE_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OdddO", &corner, &u, &v, &w, &scolor);

    TYPETEST (is_tuple (corner, kwl[0], 3) && is_positive (u, kwl [1]) &&
	      is_positive (v, kwl [2]) && is_positive (w, kwl [3]) &&
	      is_tuple (scolor, kwl[5], 6));

    p [0] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (corner, 0));
    p [1] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (corner, 1)); 
    p [2] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (corner, 2));

    /* -1, 0, 0 */
    ERRMEM (a = calloc (1, sizeof (struct shape)));
    ERRMEM (h = malloc (sizeof (struct halfspace)));
    VECTOR (h->p, p[0], p[1]+0.5*v, p[2]+0.5*w);
    VECTOR (h->n, -1, 0, 0);
    h->r = ALG_SQR2 * MAX (v, w) / 2.;
    h->s = 1.0;
    h->scolor = PyLong_AsLong (PyTuple_GetItem (scolor, 0));
    a->what = HSP;
    a->data = h;

   /* 0, -1, 0 */
    ERRMEM (b = calloc (1, sizeof (struct shape)));
    ERRMEM (h = malloc (sizeof (struct halfspace)));
    VECTOR (h->p, p[0]+0.5*u, p[1], p[2]+0.5*w);
    VECTOR (h->n, 0, -1, 0);
    h->r = ALG_SQR2 * MAX (u, w) / 2.;
    h->s = 1.0;
    h->scolor = PyLong_AsLong (PyTuple_GetItem (scolor, 1));
    b->what = HSP;
    b->data = h;

    /* 0, 0, -1 */
    ERRMEM (c = calloc (1, sizeof (struct shape)));
    ERRMEM (h = malloc (sizeof (struct halfspace)));
    VECTOR (h->p, p[0]+0.5*u, p[1]+0.5*v, p[2]);
    VECTOR (h->n, 0, 0, -1);
    h->r = ALG_SQR2 * MAX (u, v) / 2.;
    h->s = 1.0;
    h->scolor = PyLong_AsLong (PyTuple_GetItem (scolor, 2));
    c->what = HSP;
    c->data = h;

    /* 1, 0, 0 */
    ERRMEM (d = calloc (1, sizeof (struct shape)));
    ERRMEM (h = malloc (sizeof (struct halfspace)));
    VECTOR (h->p, p[0]+u, p[1]+0.5*v, p[2]+0.5*w);
    VECTOR (h->n, 1, 0, 0);
    h->r = ALG_SQR2 * MAX (v, w) / 2.;
    h->s = 1.0;
    h->scolor = PyLong_AsLong (PyTuple_GetItem (scolor, 3));
    d->what = HSP;
    d->data = h;

   /* 0, 1, 0 */
    ERRMEM (e = calloc (1, sizeof (struct shape)));
    ERRMEM (h = malloc (sizeof (struct halfspace)));
    VECTOR (h->p, p[0]+0.5*u, p[1]+v, p[2]+0.5*w);
    VECTOR (h->n, 0, 1, 0);
    h->r = ALG_SQR2 * MAX (u, w) / 2.;
    h->s = 1.0;
    h->scolor = PyLong_AsLong (PyTuple_GetItem (scolor, 4));
    e->what = HSP;
    e->data = h;

    /* 0, 0, 1 */
    ERRMEM (f = calloc (1, sizeof (struct shape)));
    ERRMEM (h = malloc (sizeof (struct halfspace)));
    VECTOR (h->p, p[0], p[1], p[2]+w);
    VECTOR (h->p, p[0]+0.5*u, p[1]+0.5*v, p[2]+w);
    VECTOR (h->n, 0, 0, 1);
    h->r = ALG_SQR2 * MAX (u, v) / 2.;
    h->s = 1.0;
    h->scolor = PyLong_AsLong (PyTuple_GetItem (scolor, 5));
    f->what = HSP;
    f->data = h;
    
    out->ptr = shape_combine (shape_combine (shape_combine (a, MUL, d), MUL, shape_combine (b, MUL, e)), MUL, shape_combine (c, MUL, f));
  }

  return (PyObject*)out;
}

/* create polygon */
static PyObject* POLYGON (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("polygon", "h", "scolor");
  PyObject *polygon, *scolor;
  struct halfspace *a;
  struct shape **s;
  double h;
  int i, n;
  SHAPE *out;

  out = (SHAPE*)SHAPE_TYPE.tp_alloc (&SHAPE_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OdO", &polygon, &h, &scolor);

    n = is_list_of_tuples (polygon, kwl[0], 3, 2);

    TYPETEST (n && is_positive (h, kwl[1]) && is_tuple (scolor, kwl[2], n+2));

    REAL e [6] = {FLT_MAX, FLT_MAX, 0, -FLT_MAX, -FLT_MAX, 0}, up [3] = {0, 0, 1}, v [3], w [3];

    ERRMEM (s = calloc (n+2, sizeof (struct shape*)));

    /* sides */
    for (i = 0; i < n; i ++)
    {
      PyObject *item = PyList_GetItem (polygon, i);
      REAL p [3] = {PyFloat_AsDouble (PyTuple_GetItem (item, 0)), PyFloat_AsDouble (PyTuple_GetItem (item, 1)), 0.0}, q [3];
      if (i+1 < n) item = PyList_GetItem (polygon, i+1);
      else item = PyList_GetItem (polygon, 0);
      q [0] = PyFloat_AsDouble (PyTuple_GetItem (item, 0)),
      q [1] = PyFloat_AsDouble (PyTuple_GetItem (item, 1)),
      q [2] = 0.0;

      if (p[0] < e[0]) e[0] = p[0];
      if (p[1] < e[1]) e[1] = p[1];
      if (p[0] > e[3]) e[3] = p[0];
      if (p[1] > e[4]) e[4] = p[1];

      ERRMEM (s [i+1] = calloc (1, sizeof (struct shape)));
      ERRMEM (a = malloc (sizeof (struct halfspace)));
      s [i+1]->what = HSP;
      s [i+1]->data = a;

      SUB (q, p, v);
      PRODUCT (v, up, a->n);
      NORMALIZE (a->n);
      MID (q, p, a->p);
      a->p [2] = 0.5*h;
      SUB (a->p, p, v);
      a->r = LEN (v);
      a->s = 1.0;
      a->scolor = PyLong_AsLong (PyTuple_GetItem (scolor, i+1));
    }

    /* base */
    ERRMEM (s [0] = calloc (1, sizeof (struct shape)));
    ERRMEM (a = malloc (sizeof (struct halfspace)));
    s [0]->what = HSP;
    s [0]->data = a;

    MID (e, e+3, a->p);
    VECTOR (a->n, 0, 0, -1);
    SUB (a->p, e, v);
    a->r = LEN (v);
    a->s = 1.0;
    a->scolor = PyLong_AsLong (PyTuple_GetItem (scolor, 0));

    /* top */
    ERRMEM (s [n+1] = calloc (1, sizeof (struct shape)));
    ERRMEM (a = malloc (sizeof (struct halfspace)));
    s [n+1]->what = HSP;
    s [n+1]->data = a;

    MID (e, e+3, a->p); a->p [2] += h;
    VECTOR (a->n, 0, 0, 1);
    a->r = LEN (v);
    a->s = 1.0;
    a->scolor = PyLong_AsLong (PyTuple_GetItem (scolor, n+1));

    /* combine */
    out->ptr = shape_combine (s[0], MUL, s[n+1]); /* base and top */

    struct list /* circular list of side planes */
    {
      struct shape *shape;
      struct list *prev, *next;
    };

    struct shape *s_l, *s_m, *s_r, *conc;
    struct list *list, *head, *item;
    struct halfspace *l, *m, *r;

    ERRMEM (list = malloc (n * sizeof (struct list)));

    /* create list */
    for (i = 0; i < n; i ++)
    {
      list [i].shape = s[i+1];

      if (i == 0)
      {
	list [i].prev = &list[n-1];
	list [i].next = &list[1];
      }
      else if (i+1 < n)
      {
	list [i].prev = &list[i-1];
	list [i].next = &list[i+1];
      }
      else
      {
	list [i].prev = &list[i-1];
	list [i].next = &list[0];
      }
    }

    /* create outer hull */
    item = list;
    head = NULL;
    do
    {
      s_l = item->prev->shape;
      s_m = item->shape;
      s_r = item->next->shape;

      l = s_l->data;
      m = s_m->data;
      r = s_r->data;

      PRODUCT (l->n, m->n, v);
      PRODUCT (m->n, r->n, w);

      if (v[2] >= 0.0 && w[2] >= 0.0) /* outer hull */
      {
	out->ptr = shape_combine (out->ptr, MUL, shape_copy (s_m));
	if (head == NULL) head = item; /* list head may be on a convex face */
      }
      else if (v[2] >= 0.0 && w[2] < 0.0) /* concavity starts */
      {
	if (head == NULL) head = item; /* list head may be here as well (e.g. gears don't have convex faces) */
      }

      item = item->next;
    } while (item != list);

    if (!head)
    {
      PyErr_SetString (PyExc_ValueError, "Your polygon definition must be wrong");
      return NULL;
    }

    /* subtract concavities */
    conc = NULL;
    item = head;
    do
    {
      s_l = item->prev->shape;
      s_m = item->shape;
      s_r = item->next->shape;

      l = s_l->data;
      m = s_m->data;
      r = s_r->data;

      PRODUCT (l->n, m->n, v);
      PRODUCT (m->n, r->n, w);

      if (v[2] >= 0.0 && w[2] < 0.0) /* concavity starts */
      {
	ASSERT (!conc, "Algorithmic error!");
	s_m = shape_copy (s_m);
	m = s_m->data;
	SCALE (m->n, -1.0);
	conc = s_m;
      }
      else if (v[2] < 0.0 && w[2] < 0.0) /* continuous concavity */
      {
	ASSERT (conc, "Algorithmic error!");
	s_m = shape_copy (s_m);
	m = s_m->data;
	SCALE (m->n, -1.0);
	conc = shape_combine (conc, MUL, s_m);
      }
      else if (v[2] < 0.0 && w[2] >= 0.0) /* concavity ends */
      {
	ASSERT (conc, "Algorithmic error!");
	s_m = shape_copy (s_m);
	m = s_m->data;
	SCALE (m->n, -1.0);
	conc = shape_combine (conc, MUL, s_m);
	shape_invert (conc);
	out->ptr = shape_combine (out->ptr, MUL, conc);
	conc = NULL;
      }

      item = item->next;
    }
    while (item != head);

    /* clean up */
    for (i = 1; i <= n; i ++) shape_destroy (s [i]);
    free (s);
    free (list);
  }

  return (PyObject*)out;
}

/* create sphere */
static PyObject* MLS__ (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("op", "r", "scolor");
  int scolor, n, i, j;
  PyObject *op, *x;
  struct mls *mls;
  double r;
  SHAPE *out;

  out = (SHAPE*)SHAPE_TYPE.tp_alloc (&SHAPE_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("Odi", &op, &r, &scolor);

    n = is_list_of_tuples (op, kwl[0], 1, 6);

    TYPETEST (n && is_positive (r, kwl[1]));

    ERRMEM (out->ptr = calloc (1, sizeof (struct shape)));
    ERRMEM (mls = malloc (sizeof (struct mls)));
    ERRMEM (mls->op = malloc (n * sizeof (REAL [6])));

    for (i = 0; i < n; i ++)
    {
      x = PyList_GetItem (op, i);

      for (j = 0; j < 6; j ++)
      {
        mls->op [i][j] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (x, j));
      }

      NORMALIZE (mls->op[i]+3);
    }
    mls->nop = n;
    mls->r = r;
    mls->s = 1.0;
    mls->scolor = scolor;

    out->ptr->what = MLS;
    out->ptr->data = mls;
  }

  return (PyObject*)out;
}

/* copy shape */
static PyObject* COPY__ (PyObject *self, PyObject *args, PyObject *kwds) /* COPY__ => alg.h has a macro COPY */
{
  KEYWORDS ("shape");
  SHAPE *shape, *out;

  out = (SHAPE*)SHAPE_TYPE.tp_alloc (&SHAPE_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("O", &shape);

    TYPETEST (is_shape (shape, kwl[0]));

    out->ptr = shape_copy (shape->ptr);
  }

  return (PyObject*)out;
}

/* union of shapes */
static PyObject* UNION (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("shape1", "shape2");
  SHAPE *shape1, *shape2, *out;

  out = (SHAPE*)SHAPE_TYPE.tp_alloc (&SHAPE_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OO", &shape1, &shape2);

    TYPETEST (is_shape (shape1, kwl[0]) && is_shape (shape2, kwl[1]));

    out->ptr = shape_combine (shape_copy (shape1->ptr), ADD, shape_copy (shape2->ptr));
  }

  return (PyObject*)out;
}

/* intersection of shapes */
static PyObject* INTERSECTION (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("shape1", "shape2");
  SHAPE *shape1, *shape2, *out;

  out = (SHAPE*)SHAPE_TYPE.tp_alloc (&SHAPE_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OO", &shape1, &shape2);

    TYPETEST (is_shape (shape1, kwl[0]) && is_shape (shape2, kwl[1]));

    out->ptr = shape_combine (shape_copy (shape1->ptr), MUL, shape_copy (shape2->ptr));
  }

  return (PyObject*)out;
}

/* difference of shapes */
static PyObject* DIFFERENCE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("shape1", "shape2");
  SHAPE *shape1, *shape2, *out;

  out = (SHAPE*)SHAPE_TYPE.tp_alloc (&SHAPE_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OO", &shape1, &shape2);

    TYPETEST (is_shape (shape1, kwl[0]) && is_shape (shape2, kwl[1]));

    out->ptr = shape_combine (shape_copy (shape1->ptr), MUL, shape_invert (shape_copy (shape2->ptr)));
  }

  return (PyObject*)out;
}

/* move shape */
static PyObject* MOVE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("shape", "vector");
  PyObject *vector;
  SHAPE *shape;
  REAL v [3];

  PARSEKEYS ("OO", &shape, &vector);

  TYPETEST (is_shape (shape, kwl[0]) && is_tuple (vector, kwl[2], 3));

  v [0] = PyFloat_AsDouble (PyTuple_GetItem (vector, 0));
  v [1] = PyFloat_AsDouble (PyTuple_GetItem (vector, 1));
  v [2] = PyFloat_AsDouble (PyTuple_GetItem (vector, 2));

  shape_move (shape->ptr, v);

  Py_RETURN_NONE;
}

/* rotate shape */
static PyObject* ROTATE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("shape", "point", "vector", "angle");
  REAL r [9], p [3], v [3];
  PyObject *point, *vector;
  double angle;
  SHAPE *shape;

  PARSEKEYS ("OOOd", &shape, &point, &vector, &angle);

  TYPETEST (is_shape (shape, kwl[0]) && is_tuple (point, kwl[1], 3) && is_tuple (vector, kwl[2], 3));

  p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
  p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
  p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

  v [0] = PyFloat_AsDouble (PyTuple_GetItem (vector, 0));
  v [1] = PyFloat_AsDouble (PyTuple_GetItem (vector, 1));
  v [2] = PyFloat_AsDouble (PyTuple_GetItem (vector, 2));

  if (LEN (v) == 0)
  {
    PyErr_SetString (PyExc_ValueError, "Rotation direction is zero");
    return NULL;
  }

  angle = (ALG_PI * angle / 180.0);

  ROTATION_MATRIX (v, angle, r);

  shape_rotate (shape->ptr, p, r);

  Py_RETURN_NONE;
}

/* create fillet */
static PyObject* FILLET (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("shape", "c", "r", "fillet", "scolor");
  double r, fillet;
  PyObject *cobj;
  SHAPE *shape;
  int scolor;
  REAL c[3];

  PARSEKEYS ("OOddi", &shape, &cobj, &r, &fillet, &scolor);

  TYPETEST (is_shape (shape, kwl[0]) && is_tuple (cobj, kwl[1], 3) && is_positive (r, kwl[2]) && is_positive (fillet, kwl[3]));

  c [0] = PyFloat_AsDouble (PyTuple_GetItem (cobj, 0));
  c [1] = PyFloat_AsDouble (PyTuple_GetItem (cobj, 1));
  c [2] = PyFloat_AsDouble (PyTuple_GetItem (cobj, 2));

  shape_fillet (shape->ptr, c, r, fillet, scolor);

  Py_RETURN_NONE;
}

static PyMethodDef methods [] =
{
  {"SPHERE", (PyCFunction)SPHERE, METH_VARARGS|METH_KEYWORDS, "Create sphere"},
  {"CYLINDER", (PyCFunction)CYLINDER, METH_VARARGS|METH_KEYWORDS, "Create cylinder"},
  {"CUBE", (PyCFunction)CUBE, METH_VARARGS|METH_KEYWORDS, "Create cube"},
  {"POLYGON", (PyCFunction)POLYGON, METH_VARARGS|METH_KEYWORDS, "Create polygon"},
  {"MLS", (PyCFunction)MLS__, METH_VARARGS|METH_KEYWORDS, "Create moving least square fit"},
  {"COPY", (PyCFunction)COPY__, METH_VARARGS|METH_KEYWORDS, "Copy shape"},
  {"UNION", (PyCFunction)UNION, METH_VARARGS|METH_KEYWORDS, "Union of shapes"},
  {"INTERSECTION", (PyCFunction)INTERSECTION, METH_VARARGS|METH_KEYWORDS, "Intersection of shapes"},
  {"DIFFERENCE", (PyCFunction)DIFFERENCE, METH_VARARGS|METH_KEYWORDS, "Difference of shapes"},
  {"MOVE", (PyCFunction)MOVE, METH_VARARGS|METH_KEYWORDS, "Move shape"},
  {"ROTATE", (PyCFunction)ROTATE, METH_VARARGS|METH_KEYWORDS, "Rotate shape"},
  {"FILLET", (PyCFunction)FILLET, METH_VARARGS|METH_KEYWORDS, "Create fillet"},
  {NULL, 0, 0, NULL}
};

/* 
 * initialization
 */

static struct PyModuleDef inputmodule = {
  PyModuleDef_HEAD_INIT,
  "oaktree",
  NULL,
  -1,
  methods
};

PyMODINIT_FUNC PyInit_input(void) {
  PyObject* m;

  /* Initialize SIMULATION_TYPE */
  TYPEINIT(SIMULATION_TYPE,
           SIMULATION,
           "oaktree.SIMULATION",
           Py_TPFLAGS_DEFAULT,
           (destructor)SIMULATION_dealloc,
           (newfunc)SIMULATION_new,
           SIMULATION_methods,
           SIMULATION_members,
           SIMULATION_getset);

  /* Initialize SHAPE_TYPE */
  TYPEINIT(SHAPE_TYPE,
           SHAPE,
           "oaktree.SHAPE",
           Py_TPFLAGS_DEFAULT,
           (destructor)SHAPE_dealloc,
           (newfunc)SHAPE_new,
           SHAPE_methods,
           SHAPE_members,
           SHAPE_getset);

  /* Initialize DOMAIN_TYPE */
  TYPEINIT(DOMAIN_TYPE,
           DOMAIN,
           "oaktree.DOMAIN",
           Py_TPFLAGS_DEFAULT,
           (destructor)DOMAIN_dealloc,
           (newfunc)DOMAIN_new,
           DOMAIN_methods,
           DOMAIN_members,
           DOMAIN_getset);

  if (PyType_Ready(&SIMULATION_TYPE) < 0)
    return NULL;

  if (PyType_Ready(&SHAPE_TYPE) < 0)
    return NULL;

  if (PyType_Ready(&DOMAIN_TYPE) < 0)
    return NULL;

  m = PyModule_Create(&inputmodule);
  if (m == NULL)
    return NULL;

  Py_INCREF(&SIMULATION_TYPE);
  PyModule_AddObject(m, "SIMULATION", (PyObject*)&SIMULATION_TYPE);

  Py_INCREF(&SHAPE_TYPE);
  PyModule_AddObject(m, "SHAPE", (PyObject*)&SHAPE_TYPE);

  Py_INCREF(&DOMAIN_TYPE);
  PyModule_AddObject(m, "DOMAIN", (PyObject*)&DOMAIN_TYPE);

  return m;
}

/* 
 * interface
 */

/* interpret an input file (return 0 on success)*/
int input (const char *path)
{
  FILE *file;
  char *line;
  int error;

  ASSERT (file = fopen (path, "r"), "File open failed!");

  Py_Initialize();

  PyObject *module = PyInit_input();
  if (!module) return -1;

  /* Add the module to sys.modules */
  PyObject *sys_modules = PyImport_GetModuleDict();
  PyDict_SetItemString(sys_modules, "oaktree", module);

  /* Now we can safely run the initialization code */
  PyRun_SimpleString("import sys\n"
                    "from oaktree import SIMULATION\n"
                    "from oaktree import SHAPE\n"
                    "from oaktree import SPHERE\n"
                    "from oaktree import CYLINDER\n"
                    "from oaktree import CUBE\n"
                    "from oaktree import POLYGON\n"
                    "from oaktree import MLS\n"
                    "from oaktree import COPY\n"
                    "from oaktree import UNION\n"
                    "from oaktree import INTERSECTION\n"
                    "from oaktree import DIFFERENCE\n"
                    "from oaktree import MOVE\n"
                    "from oaktree import ROTATE\n"
                    "from oaktree import FILLET\n"
                    "from oaktree import DOMAIN\n");

  ERRMEM (line = malloc (128 + strlen (path)));
  sprintf (line, "exec(open('%s').read())", path);

  error = PyRun_SimpleString (line);
  fclose (file);
  free (line);

  return error;
}
