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

static PyTypeObject SHAPE_TYPE;

typedef struct {
  PyObject_HEAD
  struct shape *shp;
} SHAPE;

/* constructor */
static PyObject* SHAPE_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyErr_SetString (PyExc_RuntimeError, "Cannot create unspecified shape!");
  return NULL;
}

/* destructor */
static void SHAPE_dealloc (SHAPE *self)
{
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

/* create superellipsoid */
static PyObject* SUPERELLIPSOID (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("center", "radii", "p", "q", "vcolor", "scolor");
  PyObject *center, *radii;
  int vcolor, scolor;
  double p, q;
  SHAPE *out;

  out = (SHAPE*)SHAPE_TYPE.tp_alloc (&SHAPE_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OOddii", &center, &radii, &p, &q, &vcolor, &scolor);

    TYPETEST (is_tuple (center, kwl[0], 3) && is_tuple (radii, kwl[1], 3)); /* TODO */
  }

  return (PyObject*)out;
}

static PyMethodDef methods [] =
{
  {"SUPERELLIPSOID", (PyCFunction)SUPERELLIPSOID, METH_VARARGS|METH_KEYWORDS, "Create superellipsoid"},
  {NULL, 0, 0, NULL}
};

/* 
 * initialization
 */

static void initinput (void)
{
  PyObject *m;

  TYPEINIT (SHAPE_TYPE, SHAPE, "solfec.SHAPE",
    Py_TPFLAGS_DEFAULT, SHAPE_dealloc, SHAPE_new,
    SHAPE_methods, SHAPE_members, SHAPE_getset);

  if (PyType_Ready (&SHAPE_TYPE) < 0) return;

  if (!(m =  Py_InitModule3 ("oaktree", methods, "oaktree module"))) return;

  Py_INCREF (&SHAPE_TYPE);

  PyModule_AddObject (m, "SHAPE", (PyObject*)&SHAPE_TYPE);
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
