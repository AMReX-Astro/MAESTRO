
F90    ?= gfortran
FFLAGS  = -Ibuild -Jbuild -g -Wall -Wno-unused-dummy-argument
CC      = gcc
CFLAGS  = -std=c99 -Wall -g -DBL_Linux -DBL_FORT_USE_UNDERSCORE 

vpath %.f90 ../ /home/ajacobs/Research/Codebase/BoxLib/Src/F_BaseLib/
vpath %.f   ../
vpath %.c /home/ajacobs/Research/Codebase/BoxLib/Src/F_BaseLib/

all: t1.exe

#
# rules
#

%.exe: %.f90 build/bl_types.o build/timer_c.o build/parallel_stubs.o build/bl_error.o build/bl_string.o build/bl_IO.o build/bl_timer.o build/bdf.o build/LinAlg.o 
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
	rm -r build/ *.exe *~

