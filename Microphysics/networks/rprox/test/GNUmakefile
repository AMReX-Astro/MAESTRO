NDEBUG := t
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

MAESTRO_TOP_DIR :=  ../../../..

# programs to be compiled
programs += testburn
programs += testjacobian
#programs += compare_new_rates.f90
#programs += rates.f90
#programs += evolve.f90

EOS_DIR := helmeos
CONDUCTIVITY_DIR := constant
NETWORK_DIR := rprox

BOXLIB_ONLY := t

EXTRA_DIR = Microphysics/networks/rprox/test

testAll: table

include $(MAESTRO_TOP_DIR)/GMaestro.mak

testAll: $(pnames)

%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif

