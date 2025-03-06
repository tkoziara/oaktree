#
# double / float
#

REAL = float

#
# Specify C compiler
#

CC = cc

# 
# Debug or optimized version switch (yes/no)
#

DEBUG = yes
PROFILE = no

#
# BLAS
#

BLAS = -L/usr/lib -lblas

#
# LAPACK
#

LAPACK = -L/usr/lib -llapack

#
# Python
#

PYTHON = -I/usr/include/python3.10
PYTHONLIB = -L/usr/lib -lpython3.10

#
# OpenGL (yes/no)
#

OPENGL = yes
GLINC =
GLLIB = -framework GLUT -framework OpenGL

#
# VBO  (OPENGL == yes)
#

VBO = yes

#
# MPI (yes/no)
#

MPI = no
MPICC = mpicc
