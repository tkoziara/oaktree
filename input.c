/*
 * input.h
 * ---------
 */

#include <Python.h>
#include <structmember.h>
#include "oaktree.h"
#include "input.h"
#include "alg.h"
#include "err.h"

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

/* test whether an object is a tuple of length len */
static int is_tuple (PyObject *obj, char *var, int len)
{
  if (obj)
  {
    if (!PyTuple_Check (obj))
    {
      char buf [BUFLEN];
      snprintf (buf, BUFLEN, "'%s' must be a tuple object", var);
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
 * SHAPE
 */

static PyTypeObject lng_SHAPE_TYPE;

typedef struct {
  PyObject_HEAD
  struct shape *shp;
} lng_SHAPE;

/* constructor */
static PyObject* lng_SHAPE_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyErr_SetString (PyExc_RuntimeError, "Cannot create unspecified shape!");
  return NULL;
}

/* destructor */
static void lng_SHAPE_dealloc (lng_SHAPE *self)
{
  self->ob_type->tp_free ((PyObject*)self);
}

/* SHAPE methods */
static PyMethodDef lng_SHAPE_methods [] =
{ {NULL, NULL, 0, NULL} };

/* SHAPE members */
static PyMemberDef lng_SHAPE_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* SHAPE getset */
static PyGetSetDef lng_SHAPE_getset [] =
{ {NULL, 0, 0, NULL, NULL} };

/*
 * SUBROUTINES
 */

/* create superellipsoid */
static PyObject* lng_SUPERELLIPSOID (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("center", "radii", "r", "t", "color");
  PyObject *center, *radii;
  lng_SHAPE *out;
  double r, t;
  int color;

  out = (lng_SHAPE*)lng_SHAPE_TYPE.tp_alloc (&lng_SHAPE_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OOddi", &center, &radii, &r, &t, &color);

    TYPETEST (is_tuple (center, kwl[0], 3) && is_tuple (radii, kwl[1], 3)); /* TODO */
  }

  return (PyObject*)out;
}

static PyMethodDef lng_methods [] =
{
  {"SUPERELLIPSOID", (PyCFunction)lng_SUPERELLIPSOID, METH_VARARGS|METH_KEYWORDS, "Create superellipsoid"},
  {NULL, 0, 0, NULL}
};

/* 
 * initialization
 */

static void initinput (void)
{
  PyObject *m;

  TYPEINIT (lng_SHAPE_TYPE, lng_SHAPE, "solfec.SHAPE",
    Py_TPFLAGS_DEFAULT, lng_SHAPE_dealloc, lng_SHAPE_new,
    lng_SHAPE_methods, lng_SHAPE_members, lng_SHAPE_getset);

  if (PyType_Ready (&lng_SHAPE_TYPE) < 0) return;

  if (!(m =  Py_InitModule3 ("oaktree", lng_methods, "oaktree module"))) return;

  Py_INCREF (&lng_SHAPE_TYPE);

  PyModule_AddObject (m, "SHAPE", (PyObject*)&lng_SHAPE_TYPE);
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

  PyRun_SimpleString("from oaktree import SHAPE\n"
                     "from oaktree import SUPERELLIPSOID\n");

  ERRMEM (line = malloc (128 + strlen (path)));
  sprintf (line, "execfile ('%s')", path);

  error = PyRun_SimpleString (line); /* we do not run a file directly because FILE destriptors differe
					between WIN32 and UNIX while Python is often provided in binary form */
  fclose (file);
  free (line);

  return error;
}
