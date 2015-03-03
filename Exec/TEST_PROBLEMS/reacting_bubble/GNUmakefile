NDEBUG := 
MPI    :=
OMP    :=

SDC := 

COMP := gfortran

MKVERBOSE := t

# define the location of the MAESTRO root directory
MAESTRO_TOP_DIR := ../../..


# define the physics packages to build this problem
EOS_DIR := helmeos 
CONDUCTIVITY_DIR := timmes_stellar

ifdef SDC
  NETWORK_DIR := ignition_simple_SDC
else
  NETWORK_DIR := ignition_simple
endif

# define the special directories needed to build this problem.  Note:
# we only need to include the problem's directory if there are unique
# files there (as specified in a GPackage.mak).  The problem directory
# is always placed at the start of the vpath by the GMakerules.mak.
EXTRA_DIR := 


# include the MAESTRO build stuff
include $(MAESTRO_TOP_DIR)/GMaestro.mak


