NDEBUG := t
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

# define the location of the fParallel root directory
FPARALLEL := ../../..


# define the physics packages to build this problem
EOS_DIR := extern/EOS/gamma_law_general 
CONDUCTIVITY_DIR := extern/conductivity/constant

NETWORK_DIR := extern/networks/general_null
GENERAL_NET_INPUTS := simple.net

# define the special directories needed to build this problem.  Note:
# we only need to include the problem's directory if there are unique
# files there (as specified in a GPackage.mak).  The problem directory
# is always placed at the start of the vpath by the GMakerules.mak.
EXTRA_DIR := extern/random


# include the MAESTRO build stuff
include $(FPARALLEL)/MAESTRO/GMaestro.mak


