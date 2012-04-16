#
# Compilation flags setup
#

ifeq ($(OPENGL),yes)
  ifeq ($(VBO),yes)
    OPENGL = -DOPENGL -DVBO $(GLINC)
  else
    OPENGL = -DOPENGL $(GLINC)
  endif
else
  OPENGL =
  GLLIB = 
endif

ifeq ($(DEBUG),yes)
  DEBUG =  -W -Wall -Wno-unused-parameter -pedantic -g -DDEBUG
  ifeq ($(PROFILE),yes)
    PROFILE = -p
  else
    PROFILE =
  endif
else
  DEBUG =  -w -pedantic -O3 -funroll-loops
  PROFILE =
endif

ifeq ($(MPI),yes)
  MPIFLG = -DMPI $(ZOLTANINC)
  MPILIBS = $(ZOLTANLIB)
endif
