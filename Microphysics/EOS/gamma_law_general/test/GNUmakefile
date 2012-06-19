NDEBUG := t
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

MAESTRO_TOP_DIR := ../../../..

EOS_DIR := gamma_law_general
CONDUCTIVITY_DIR := constant

NETWORK_DIR := general_null
GENERAL_NET_INPUTS := $(MAESTRO_TOP_DIR)/Microphysics/networks/$(NETWORK_DIR)/ignition.net

BOXLIB_ONLY := t
EXTRA_DIR = Microphysics/EOS/helmeos/test

testALL: testeos

# define the packages to build this problem
#   note: we only need to include the problem's directory if there are 
#   unique files there (as specified in a GPackage.mak).  The problem 
#   directory is always placed at the start of the vpath by the 
#   GMakerules.mak.

include $(MAESTRO_TOP_DIR)/GMaestro.mak

# define the build instructions for the executable
testeos: $(objects)
	$(LINK.f90) -o testeos.$(suf).exe $(objects) $(libraries)
	@echo SUCCESS



