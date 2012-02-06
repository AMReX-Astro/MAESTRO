NDEBUG := t
MPI    :=
OMP    :=

COMP := Intel

MKVERBOSE := t

# define the location of the fParallel root directory
FPARALLEL := ../../..


# define the physics packages to build this problem 
EOS_DIR := extern/EOS/helmeos 
CONDUCTIVITY_DIR := extern/conductivity/timmes_stellar
NETWORK_DIR := extern/networks/null

# define the special directories needed to build this problem.  Note:
# we only need to include the problem's directory if there are unique
# files there (as specified in a GPackage.mak).  The problem directory
# is always placed at the start of the vpath by the GMakerules.mak.
EXTRA_DIR := 


# define the runtime parameter files
PROBIN_PARAMETERS := ../_parameters


# include the MAESTRO build stuff
include $(FPARALLEL)/MAESTRO/GMaestro.mak


