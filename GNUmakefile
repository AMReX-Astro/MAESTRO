NDEBUG := t
MPI    := t
OMP    :=

SDC :=

COMP := gfortran

MKVERBOSE := t

# define the location of the MAESTRO top directory
MAESTRO_TOP_DIR := $(MAESTRO_HOME)


# define the physics packages to build this problem 
EOS_DIR := helmeos
CONDUCTIVITY_DIR := timmes_stellar

NETWORK_TOP_DIR := $(ASTRODEV_DIR)/networks

ifdef SDC
  NETWORK_DIR := approx8_SDC
else
  NETWORK_DIR := approx8
endif

# define the special directories needed to build this problem.  Note:
# we only need to include the problem's directory if there are unique
# files there (as specified in a GPackage.mak).  The problem directory
# is always placed at the start of the vpath by the GMakerules.mak.
EXTRA_DIR :=


# include the MAESTRO build stuff
include $(MAESTRO_TOP_DIR)/GMaestro.mak


