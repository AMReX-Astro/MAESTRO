NDEBUG := t
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

# location of root directory
MAESTRO_TOP_DIR := $(MAESTRO_HOME)


# programs to be compiled
programs += t1

#EOS_DIR := helmeos
#CONDUCTIVITY_DIR := constant

#NETWORK_TOP_DIR := $(ASTRODEV_DIR)/networks
#NETWORK_DIR := triple_alpha_plus_cago

BOXLIB_ONLY := t

EXTRA_TOP_DIR := $(ASTRODEV_DIR)
EXTRA_DIR =  bdf/tests

NEED_VBDF := t

testAll:

include $(MAESTRO_TOP_DIR)/GMaestro.mak

testAll: $(pnames)


%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif



#F90    ?= gfortran
#FFLAGS  = -Ibuild -Jbuild -g -Wall -Wno-unused-dummy-argument
#
#vpath %.f90 ../ /home/ajacobs/Research/Codebase/BoxLib/Src/F_BaseLib/
#vpath %.f   ../
#
#all: t1.exe
#
##
## rules
##
#
#%.exe: %.f90 build/bl_types.o build/parallel.o build/bl_error.o build/bdf.o build/LinAlg.o 
#	$(F90) $(FFLAGS) $^ -o $@
#
#build/%.o: %.f
#	@mkdir -p build
#	$(F90) -c $(FFLAGS) $^ -o $@
#
#build/%.o: %.f90
#	@mkdir -p build
#	$(F90) -c $(FFLAGS) $^ -o $@
#
#build/%.mod: %.f90 
#	@mkdir -p build
#	$(F90) -c $(FFLAGS) $^ -o $*.o
#
#clean::
#	rm -r build/ *.exe 
#
