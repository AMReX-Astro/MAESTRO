NDEBUG := t
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

# location of root directory
MAESTRO_HOME ?= ../../../../MAESTRO
MAESTRO_TOP_DIR := $(MAESTRO_HOME)


# programs to be compiled
programs += testburn
programs += testjacobian
#programs += compare_new_rates.f90
#programs += rates.f90
#programs += evolve.f90

EOS_DIR := helmeos
CONDUCTIVITY_DIR := constant

NETWORK_TOP_DIR := $(ASTRODEV_DIR)/networks
NETWORK_DIR := rprox

BOXLIB_ONLY := t

EXTRA_TOP_DIR := $(ASTRODEV_DIR)
EXTRA_DIR = networks/rprox/test

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

