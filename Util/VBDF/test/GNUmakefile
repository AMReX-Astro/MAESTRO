
F90    ?= gfortran
FFLAGS  = -Ibuild -Jbuild -g -Wall -Wno-unused-dummy-argument
CC      = gcc
CFLAGS  = -std=c99 -Wall -g -DBL_Linux -DBL_FORT_USE_UNDERSCORE 

BOXLIB = build/bl_types.o build/timer_c.o build/parallel_stubs.o build/bl_error.o build/bl_string.o build/bl_IO.o build/bl_timer.o
LINALG = build/idamax.o build/dscal.o build/daxpy.o build/vddot.o build/dgefa.o build/dgesl.o

vpath %.f90 ../ /home/ajacobs/Research/Codebase/BoxLib/Src/F_BaseLib/
vpath %.f   ../
vpath %.f   ../../BLAS
vpath %.f   ../../VODE
vpath %.c /home/ajacobs/Research/Codebase/BoxLib/Src/F_BaseLib/

all: t1.exe coeffs.exe

#
# rules
#

%.exe: %.f90 $(BOXLIB) $(LINALG) build/bdf.o 
	$(F90) $(FFLAGS) $^ -o $@

build/%.o: %.c
	@mkdir -p build
	$(CC) $(CFLAGS) -c $^ -o $@

build/%.o: %.f
	@mkdir -p build
	$(F90) -c $(FFLAGS) $^ -o $@

build/%.o: %.f90
	@mkdir -p build
	$(F90) -c $(FFLAGS) $^ -o $@

clean::
	rm -rf build/ *.exe *~

