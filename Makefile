include Config.mak

OB0 =   obj/input.o \
	obj/polygon.o \
	obj/octree.o \
	obj/shape.o \

ifeq ($(OPENGL),yes)

OBJ =   obj/viewer.o \
	obj/render.o \
	$(OB0) \

else

OBJ =  $(OB0)

endif

include Flags.mak

CFLAGS = -std=c99 $(DEBUG) $(PROFILE) $(REAL)

LIB = -lm $(LAPACK) $(BLAS) $(GLLIB) $(PYTHONLIB)

ifeq ($(MPI),yes)
  LIBMPI = -lm $(LAPACK) $(BLAS) $(PYTHONLIB) $(ZOLTANLIB)
endif

ifeq ($(MPI),yes)

all: oaktree oaktree-mpi

oaktree-mpi: obj/oaktree-mpi.o $(OB0)
	$(MPICC) $(PROFILE) -o $@ $< $(OB0) $(LIBMPI)

else

all: oaktree

endif

oaktree: obj/oaktree.o $(OBJ)
	$(CC) $(PROFILE) -o $@ $< $(OBJ) $(LIB)

del:
	rm -fr out/*
	rm -fr *cubin
	rm -fr *dSYM

clean:
	rm -f oaktree
	rm -f oaktree-mpi
	rm -fr out/*
	rm -f core obj/*.o
	rm -f obj/*.a
	rm -fr *dSYM
	rm -fr *cubin

obj/viewer.o: viewer.c viewer.h error.h alg.h
	$(CC) $(CFLAGS) $(OPENGL) -c -o $@ $<

obj/render.o: render.c render.h oaktree.h error.h alg.h
	$(CC) $(CFLAGS) $(OPENGL) -c -o $@ $<

obj/input.o: input.c input.h oaktree.h error.h alg.h
	$(CC) $(CFLAGS) $(PYTHON) -c -o $@ $<

obj/polygon.o: polygon.c polygon.h error.h alg.h
	$(CC) $(CFLAGS) -c -o $@ $<

obj/octree.o: octree.c oaktree.h polygon.h error.h alg.h
	$(CC) $(CFLAGS) -c -o $@ $<

obj/shape.o: shape.c oaktree.h error.h alg.h
	$(CC) $(CFLAGS) -c -o $@ $<

obj/oaktree.o: oaktree.c oaktree.h viewer.h render.h input.h timer.h error.h alg.h
	$(CC) $(CFLAGS) $(OPENGL) -c -o $@ $<

# MPI

obj/oaktree-mpi.o: oaktree.c oaktree.h input.h error.h alg.h
	$(CC) $(CFLAGS) $(MPIFLAGS) -c -o $@ $<
