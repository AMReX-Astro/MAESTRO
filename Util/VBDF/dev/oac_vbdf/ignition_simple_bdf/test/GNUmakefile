NDEBUG := 
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

MAESTRO_TOP_DIR := $(MAESTRO_HOME)

programs += testburn
programs += eval             
programs += testjacobian

EOS_DIR := helmeos
CONDUCTIVITY_DIR := constant
NETWORK_TOP_DIR := ../..
NETWORK_DIR := ignition_simple_bdf

BOXLIB_ONLY := t
EXTRA_DIR = ../AstroDev/networks/ignition_simple_bdf/test

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



