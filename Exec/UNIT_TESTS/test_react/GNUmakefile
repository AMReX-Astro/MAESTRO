NDEBUG := t
MPI    := 
OMP    := t
ACC    :=


COMP := gfortran

PARTICLES := t

MKVERBOSE := t

# define the location of the MAESTRO top directory
MAESTRO_TOP_DIR := ../../..

# this unit test requires the Microphysics repo -- make sure we have it
ifndef MICROPHYSICS_HOME
   $(error test_react requires the Microphysics repo. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)
endif

# define the physics packages to build this problem
EOS_DIR := helmholtz
CONDUCTIVITY_DIR := stellar

#NETWORK_TOP_DIR := $(ASTRODEV_DIR)/networks
#NETWORK_DIR := rprox
#NETWORK_DIR := triple_alpha
#NETWORK_DIR := ignition_chamulak
#NETWORK_DIR := aprox13
NETWORK_DIR := ignition_simple

# define the integrator you would like to use, default is VODE
#INTEGRATOR_DIR := VBDF

# a special case for the test suite
ifeq ($(NETWORK_DIR), ignition_chamulak)
  NETWORK_TOP_DIR := $(MAESTRO_TOP_DIR)/Microphysics/networks
endif

# define the special directories needed to build this problem.  Note:
# we only need to include the problem's directory if there are unique
# files there (as specified in a GPackage.mak).  The problem directory
# is always placed at the start of the vpath by the GMakerules.mak.
EXTRA_DIR := Exec/UNIT_TESTS/test_react

UNIT_TEST := t

# include the MAESTRO build stuff
include $(MAESTRO_TOP_DIR)/GMaestro.mak
