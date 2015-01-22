NDEBUG := t
MPI    := t
OMP    :=

COMP := gfortran

MKVERBOSE := t

# define the location of the MAESTRO top directory
MAESTRO_TOP_DIR := $(MAESTRO_HOME)


# define the physics packages to build this problem
EOS_DIR := multigamma
CONDUCTIVITY_DIR := constant
NETWORK_DIR := general_null

GENERAL_NET_INPUTS := $(MAESTRO_TOP_DIR)/Microphysics/networks/$(NETWORK_DIR)/simple3.net


# define the special directories needed to build this problem.  Note:
# we only need to include the problem's directory if there are unique
# files there (as specified in a GPackage.mak).  The problem directory
# is always placed at the start of the vpath by the GMakerules.mak.
EXTRA_DIR :=


# include the MAESTRO build stuff
include $(MAESTRO_TOP_DIR)/GMaestro.mak


