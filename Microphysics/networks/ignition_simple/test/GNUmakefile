NDEBUG := 
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

MAESTRO_TOP_DIR := ../../../..

programs += testburn
programs += eval
programs += testjacobian

EOS_DIR := gamma_law_general
CONDUCTIVITY_DIR := constant
NETWORK_DIR := ignition_simple

BOXLIB_ONLY := t
EXTRA_DIR = Microphysics/networks/ignition_simple/test

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



