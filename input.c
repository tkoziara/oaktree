/*
 * input.h
 * ---------
 */

#include <Python.h>
#include <structmember.h>
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
    if (!PyString_Check (obj))
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

/* define keywords */
#define KEYWORDS(...) char *kwl [] = {__VA_ARGS__, NULL}

/* parse arguments with keywords */
#define PARSEKEYS(fmt, ...) if (!PyArg_ParseTupleAndKeywords (args, kwds, fmt, kwl, __VA_ARGS__)) return NULL

/* parse arguments without keywords */
#define PARSE(fmt, ...) if (!PyArg_ParseTuple (args, fmt, __VA_ARGS__)) return NULL

/* object types assertion */
#define TYPETEST(test) if(!(test)) return NULL

/* string argument if block comparison */
#define IFIS(obj, val) if (strcmp (PyString_AsString (obj), val) == 0)
#define ELIF(obj, val) else if (strcmp (PyString_AsString (obj), val) == 0)
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
  KEYWORDS ("outpath", "duration", "step", "grid", "cutoff", "extents");
  double duration, step, grid, cutoff;
  PyObject *outpath, *extents;
  struct simulation *simu;
  SIMULATION *self;

  self = (SIMULATION*)type->tp_alloc (type, 0);

  if (self)
  {
    PARSEKEYS ("OddddO", &outpath, &duration, &step, &grid, &cutoff, &extents);

    TYPETEST (is_string (outpath, kwl [0]) && is_positive (duration, kwl[1]) &&
	      is_positive (step, kwl[2]) && is_positive (grid, kwl[3]) &&
	      is_positive (cutoff, kwl[4]) && is_tuple (extents, kwl[5], 6));

    ERRMEM (simu = calloc (1, sizeof (struct simulation)));
    ERRMEM (simu->outpath = malloc (strlen (PyString_AsString (outpath)) + 1));
    strcpy (simu->outpath, PyString_AsString (outpath));
    simu->duration = duration;
    simu->step = step;
    simu->grid = grid;
    simu->cutoff = cutoff;
    simu->extents [0] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (extents, 0));
    simu->extents [1] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (extents, 1));
    simu->extents [2] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (extents, 2));
    simu->extents [3] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (extents, 3));
    simu->extents [4] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (extents, 4));
    simu->extents [5] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (extents, 5));

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
  self->ob_type->tp_free ((PyObject*)self);
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

  self->ob_type->tp_free ((PyObject*)self);
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
 * SUBROUTINES
 */

/* create sphere */
static PyObject* SPHERE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("center", "r", "vcolor", "scolor");
  struct sphere *sphere;
  int vcolor, scolor;
  PyObject *center;
  double r;
  SHAPE *out;

  out = (SHAPE*)SHAPE_TYPE.tp_alloc (&SHAPE_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("Odii", &center, &r, &vcolor, &scolor);

    TYPETEST (is_tuple (center, kwl[0], 3) && is_positive (r, kwl[1]));

    ERRMEM (out->ptr = calloc (1, sizeof (struct shape)));
    ERRMEM (sphere = malloc (sizeof (struct sphere)));

    sphere->c [0] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (center, 0));
    sphere->c [1] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (center, 1));
    sphere->c [2] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (center, 2));
    sphere->r = r;
    sphere->vcolor = vcolor;
    sphere->scolor = scolor;

    out->ptr->what = SPH;
    out->ptr->data = sphere;
  }

  return (PyObject*)out;
}

/* create cylinder */
static PyObject* CYLINDER (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("base", "h", "r", "vcolor", "scolor");
  struct shape *sa, *sb, *sc;
  struct halfplane *a, *b;
  PyObject *base, *scolor;
  struct cylinder *c;
  double r, h;
  REAL x [3];
  int vcolor;
  SHAPE *out;

  out = (SHAPE*)SHAPE_TYPE.tp_alloc (&SHAPE_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OddiO", &base, &h, &r, &vcolor, &scolor);

    TYPETEST (is_tuple (base, kwl[0], 3) && is_positive (h, kwl[1]) && is_positive (r, kwl[2]) && is_tuple (scolor, kwl[4], 3));

    ERRMEM (sa = calloc (1, sizeof (struct shape)));
    ERRMEM (a = malloc (sizeof (struct halfplane)));

    a->p [0] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (base, 0));
    a->p [1] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (base, 1));
    a->p [2] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (base, 2));
    VECTOR (a->n, 0, 0, -1);
    a->vcolor = vcolor;
    a->scolor = PyInt_AsLong (PyTuple_GetItem (scolor, 0));
    a->r = r;
    sa->what = HPL;
    sa->data = a;

    ERRMEM (sb = calloc (1, sizeof (struct shape)));
    ERRMEM (b = malloc (sizeof (struct halfplane)));

    VECTOR (b->p, a->p[0], a->p[1], a->p[2]+h);
    VECTOR (b->n, 0, 0, 1);
    b->vcolor = vcolor;
    b->scolor = PyInt_AsLong (PyTuple_GetItem (scolor, 2));
    b->r = r;
    sb->what = HPL;
    sb->data = b;

    ERRMEM (sc = calloc (1, sizeof (struct shape)));
    ERRMEM (c = malloc (sizeof (struct cylinder)));

    COPY (a->p, c->p);
    VECTOR (c->d, 0, 0, 1);
    c->r = r;
    c->vcolor = vcolor;
    c->scolor = PyInt_AsLong (PyTuple_GetItem (scolor, 1));
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
  KEYWORDS ("corner", "u", "v", "w", "vcolor", "scolor");
  struct shape *shape, *a, *b, *c, *d, *e, *f;
  PyObject *corner, *scolor;
  struct halfplane *h;
  double u, v, w;
  double p [3];
  int vcolor;
  SHAPE *out;

  out = (SHAPE*)SHAPE_TYPE.tp_alloc (&SHAPE_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OdddiO", &corner, &u, &v, &w, &vcolor, &scolor);

    TYPETEST (is_tuple (corner, kwl[0], 3) && is_positive (u, kwl [1]) &&
	      is_positive (v, kwl [2]) && is_positive (w, kwl [3]) &&
	      is_tuple (corner, kwl[5], 3));

    p [0] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (corner, 0));
    p [1] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (corner, 1)); 
    p [2] = (REAL) PyFloat_AsDouble (PyTuple_GetItem (corner, 2));

    /* -1, 0, 0 */
    ERRMEM (a = calloc (1, sizeof (struct shape)));
    ERRMEM (h = malloc (sizeof (struct halfplane)));
    VECTOR (h->p, p[0], p[1]+0.5*v, p[2]+0.5*w);
    VECTOR (h->n, -1, 0, 0);
    h->r = ALG_SQR2 * MAX (v, w) / 2.;
    h->vcolor = vcolor;
    h->scolor = PyInt_AsLong (PyTuple_GetItem (scolor, 0));
    a->what = HPL;
    a->data = h;

   /* 0, -1, 0 */
    ERRMEM (b = calloc (1, sizeof (struct shape)));
    ERRMEM (h = malloc (sizeof (struct halfplane)));
    VECTOR (h->p, p[0]+0.5*u, p[1], p[2]+0.5*w);
    VECTOR (h->n, 0, -1, 0);
    h->r = ALG_SQR2 * MAX (u, w) / 2.;
    h->vcolor = vcolor;
    h->scolor = PyInt_AsLong (PyTuple_GetItem (scolor, 1));
    b->what = HPL;
    b->data = h;

    /* 0, 0, -1 */
    ERRMEM (c = calloc (1, sizeof (struct shape)));
    ERRMEM (h = malloc (sizeof (struct halfplane)));
    VECTOR (h->p, p[0]+0.5*u, p[1]+0.5*v, p[2]);
    VECTOR (h->n, 0, 0, -1);
    h->r = ALG_SQR2 * MAX (u, v) / 2.;
    h->vcolor = vcolor;
    h->scolor = PyInt_AsLong (PyTuple_GetItem (scolor, 2));
    c->what = HPL;
    c->data = h;

    /* 1, 0, 0 */
    ERRMEM (d = calloc (1, sizeof (struct shape)));
    ERRMEM (h = malloc (sizeof (struct halfplane)));
    VECTOR (h->p, p[0]+u, p[1]+0.5*v, p[2]+0.5*w);
    VECTOR (h->n, 1, 0, 0);
    h->r = ALG_SQR2 * MAX (v, w) / 2.;
    h->vcolor = vcolor;
    h->scolor = PyInt_AsLong (PyTuple_GetItem (scolor, 3));
    d->what = HPL;
    d->data = h;

   /* 0, 1, 0 */
    ERRMEM (e = calloc (1, sizeof (struct shape)));
    ERRMEM (h = malloc (sizeof (struct halfplane)));
    VECTOR (h->p, p[0]+0.5*u, p[1]+v, p[2]+0.5*w);
    VECTOR (h->n, 0, 1, 0);
    h->r = ALG_SQR2 * MAX (u, w) / 2.;
    h->vcolor = vcolor;
    h->scolor = PyInt_AsLong (PyTuple_GetItem (scolor, 4));
    e->what = HPL;
    e->data = h;

    /* 0, 0, 1 */
    ERRMEM (f = calloc (1, sizeof (struct shape)));
    ERRMEM (h = malloc (sizeof (struct halfplane)));
    VECTOR (h->p, p[0], p[1], p[2]+w);
    VECTOR (h->p, p[0]+0.5*u, p[1]+0.5*v, p[2]+w);
    VECTOR (h->n, 0, 0, 1);
    h->r = ALG_SQR2 * MAX (u, v) / 2.;
    h->vcolor = vcolor;
    h->scolor = PyInt_AsLong (PyTuple_GetItem (scolor, 5));
    f->what = HPL;
    f->data = h;

    shape = shape_combine (shape_combine (shape_combine (a, MUL, d), MUL, shape_combine (b, MUL, e)), MUL, shape_combine (c, MUL, f));

    out->ptr = shape;
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

    out->ptr = shape_combine (shape_copy (shape1->ptr, NULL), ADD, shape_copy (shape2->ptr, NULL));
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

    out->ptr = shape_combine (shape_copy (shape1->ptr, NULL), MUL, shape_copy (shape2->ptr, NULL));
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

    out->ptr = shape_combine (shape_copy (shape1->ptr, NULL), SUB, shape_copy (shape2->ptr, NULL));
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

/* create solid */
static PyObject* SOLID (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("simu", "shape", "label");
  struct shape *copy;
  SIMULATION *simu;
  PyObject *label;
  SHAPE *shape;

  PARSEKEYS ("OOO", &simu, &shape, &label);

  TYPETEST (is_simulation (simu, kwl[0]) && is_shape (shape, kwl[1]) && is_string (label, kwl[2]));

  copy = shape_copy (shape->ptr, PyString_AsString (label));

  if (simu->ptr->solids)
    simu->ptr->solids->prev = copy;
  copy->next = simu->ptr->solids;
  simu->ptr->solids = copy; /* insert copy into list of solids */

  Py_RETURN_NONE;
}

static PyMethodDef methods [] =
{
  {"SPHERE", (PyCFunction)SPHERE, METH_VARARGS|METH_KEYWORDS, "Create sphere"},
  {"CYLINDER", (PyCFunction)CYLINDER, METH_VARARGS|METH_KEYWORDS, "Create cylinder"},
  {"CUBE", (PyCFunction)CUBE, METH_VARARGS|METH_KEYWORDS, "Create cube"},
  {"UNION", (PyCFunction)UNION, METH_VARARGS|METH_KEYWORDS, "Union of shapes"},
  {"INTERSECTION", (PyCFunction)INTERSECTION, METH_VARARGS|METH_KEYWORDS, "Intersection of shapes"},
  {"DIFFERENCE", (PyCFunction)DIFFERENCE, METH_VARARGS|METH_KEYWORDS, "Difference of shapes"},
  {"MOVE", (PyCFunction)MOVE, METH_VARARGS|METH_KEYWORDS, "Move shape"},
  {"ROTATE", (PyCFunction)ROTATE, METH_VARARGS|METH_KEYWORDS, "Rotate shape"},
  {"SOLID", (PyCFunction)SOLID, METH_VARARGS|METH_KEYWORDS, "Create solid"},
  {NULL, 0, 0, NULL}
};

/* 
 * initialization
 */

static void initinput (void)
{
  PyObject *m;

  if (!(m =  Py_InitModule3 ("oaktree", methods, "oaktree module"))) return;

  TYPEINIT (SHAPE_TYPE, SHAPE, "solfec.SHAPE",
    Py_TPFLAGS_DEFAULT, SHAPE_dealloc, SHAPE_new,
    SHAPE_methods, SHAPE_members, SHAPE_getset);
  if (PyType_Ready (&SHAPE_TYPE) < 0) return;
  Py_INCREF (&SHAPE_TYPE);
  PyModule_AddObject (m, "SHAPE", (PyObject*)&SHAPE_TYPE);

  TYPEINIT (SIMULATION_TYPE, SIMULATION, "solfec.SIMULATION",
    Py_TPFLAGS_DEFAULT, SIMULATION_dealloc, SIMULATION_new,
    SIMULATION_methods, SIMULATION_members, SIMULATION_getset);
  if (PyType_Ready (&SIMULATION_TYPE) < 0) return;
  Py_INCREF (&SIMULATION_TYPE);
  PyModule_AddObject (m, "SIMULATION", (PyObject*)&SIMULATION_TYPE);
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

  initinput ();

  PyRun_SimpleString("from oaktree import SIMULATION\n"
                     "from oaktree import SHAPE\n"
                     "from oaktree import SPHERE\n"
                     "from oaktree import CYLINDER\n"
                     "from oaktree import CUBE\n"
                     "from oaktree import UNION\n"
                     "from oaktree import INTERSECTION\n"
                     "from oaktree import DIFFERENCE\n"
                     "from oaktree import MOVE\n"
                     "from oaktree import ROTATE\n"
                     "from oaktree import SOLID\n");

  ERRMEM (line = malloc (128 + strlen (path)));
  sprintf (line, "execfile ('%s')", path);

  error = PyRun_SimpleString (line); /* we do not run a file directly because FILE destriptors differe
					between WIN32 and UNIX while Python is often provided in binary form */
  fclose (file);
  free (line);

  return error;
}
