NDEBUG := t
MPI    := t
OMP    := t

COMP := PathScale

MKVERBOSE := t

# define the location of the fParallel root directory
FPARALLEL := ../../../..

# include the main Makefile stuff
include $(FPARALLEL)/mk/GMakedefs.mak

# define the physics packages to build this problem
EOS_DIR := extern/EOS/helmeos 
CONDUCTIVITY_DIR := extern/conductivity/constant
NETWORK_DIR := extern/networks/H_core_1spec

# define the special directories needed to build this problem.  Note:
# we only need to include the problem's directory if there are unique
# files there (as specified in a GPackage.mak).  The problem directory
# is always placed at the start of the vpath by the GMakerules.mak.
#
# Include MAESTRO/Hcore to use my version of the model parser that 
# contains arrays for Abar and Zbar
EXTRA_DIR := extern/random MAESTRO/H_core


# define the runtime parameter files
PROBIN_PARAMETERS := ../../_parameters ../_parameters.hcore


# include the MAESTRO build stuff
include $(FPARALLEL)/MAESTRO/GMaestro.mak

