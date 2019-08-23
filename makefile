CC=g++
CFLAGS= -I$(IDIR) -I.

LIBS = -lm

DEPS = all.h iterative.h globals.h point.h setup.h alglibinternal.h alglibmisc.h integration.h stdafx.h ap.h linalg.h dataanalysis.h diffequations.h fasttransforms.h interpolation.h optimization.h solvers.h specialfunctions.h statistics.h misc_calc.h iterative_exp.h

OBJ =  main.o iterative.o globals.o point.o setup.o alglibinternal.o alglibmisc.o integration.o ap.o linalg.o dataanalysis.o diffequations.o fasttransforms.o interpolation.o optimization.o solvers.o specialfunctions.o statistics.o misc_calc.o iterative_exp.o

%.o: %.c -o $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

FBILI2D: $(OBJ)
	g++ -o $@ $^ $(CFLAGS)

