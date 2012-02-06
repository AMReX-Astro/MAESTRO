NDEBUG := t
MPI    := t
OMP    :=

COMP := gfortran

MKVERBOSE := t

# define the location of the fParallel root directory
FPARALLEL := ../../..


# define the physics packages to build this problem
EOS_DIR := extern/EOS/gamma_law_general
CONDUCTIVITY_DIR := extern/conductivity/timmes_stellar
NETWORK_DIR := extern/networks/gammalaw

# define the special directories needed to build this problem.  Note:
# we only need to include the problem's directory if there are unique
# files there (as specified in a GPackage.mak).  The problem directory
# is always placed at the start of the vpath by the GMakerules.mak.
EXTRA_DIR :=


# include the MAESTRO build stuff
include $(FPARALLEL)/MAESTRO/GMaestro.mak


