NDEBUG :=
MPI    :=
OMP    :=
ACC    := t

#COMP := Cray
#COMP := PGI
COMP := gfortran

MKVERBOSE := t

# define the location of the MAESTRO top directory
MAESTRO_TOP_DIR := $(MAESTRO_HOME)

# define the physics packages to build this problem
MICROPHYSICS_HOME := /home/ajacobs/Codebase/MAESTRO/Microphysics
EOS_DIR := helmeos
CONDUCTIVITY_DIR := timmes_stellar

NETWORK_TOP_DIR := $(MAESTRO_HOME)/Util/VBDF/dev/oac_vbdf
NETWORK_DIR := ignition_simple_bdf
#NETWORK_DIR := rprox_vbdf

# define the special directories needed to build this problem.  Note:
# we only need to include the problem's directory if there are unique
# files there (as specified in a GPackage.mak).  The problem directory
# is always placed at the start of the vpath by the GMakerules.mak.
EXTRA_TOP_DIR := $(MAESTRO_HOME)/Util/VBDF/dev
EXTRA_DIR := oac_vbdf/test_react

UNIT_TEST := t

# include the MAESTRO build stuff
include $(MAESTRO_TOP_DIR)/GMaestro.mak
