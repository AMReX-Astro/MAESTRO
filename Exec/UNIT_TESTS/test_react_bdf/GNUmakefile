NDEBUG :=
MPI    :=
OMP    :=
ACC    := t

COMP := gfortran

MKVERBOSE := t

# define the location of the MAESTRO top directory
MAESTRO_TOP_DIR := $(MAESTRO_HOME)

# this unit test requires the Microphysics repo -- make sure we have it
ifndef MICROPHYSICS_HOME
   $(error test_react requires the Microphysics repo. Please ensure that you have downloaded it and set $$MICROPHYSICS_HOME appropriately)
endif

# define the physics packages to build this problem
EOS_DIR := helmholtz
CONDUCTIVITY_DIR := timmes_stellar

#NETWORK_DIR := ignition_simple_bdf
NETWORK_DIR := rprox_bdf

# define the special directories needed to build this problem.  Note:
# we only need to include the problem's directory if there are unique
# files there (as specified in a GPackage.mak).  The problem directory
# is always placed at the start of the vpath by the GMakerules.mak.
EXTRA_TOP_DIR := $(MAESTRO_HOME)/Exec/UNIT_TESTS
EXTRA_DIR := test_react_bdf

UNIT_TEST := t

# include the MAESTRO build stuff
include $(MAESTRO_TOP_DIR)/GMaestro.mak
