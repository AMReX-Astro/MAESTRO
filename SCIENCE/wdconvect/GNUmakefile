NDEBUG :=  
MPI    := t
OMP    :=

HYPRE :=

COMP := Intel

MKVERBOSE := t

# define the location of the fParallel root directory
FPARALLEL := ../..

# define the physics packages to build this problem
EOS_DIR := extern/EOS/helmeos   
CONDUCTIVITY_DIR := extern/conductivity/timmes_stellar
NETWORK_DIR := extern/networks/ignition_chamulak 

# define the special directories needed to build this problem.  Note:
# we only need to include the problem's directory if there are unique
# files there (as specified in a GPackage.mak).  The problem directory
# is always placed at the start of the vpath by the GMakerules.mak.
EXTRA_DIR := extern/random

ifdef HYPRE
 
# on manda
  HYPRE_DIR = ../../extern/hypre
 
# on hopper
# HYPRE_DIR =  /usr/common/acts/Hypre/2.6.0b/src/cray-xt5_O
 
  LDFLAGS += -L $(HYPRE_DIR)/lib
  FINCLUDE_LOCATIONS += $(HYPRE_DIR)/include
  hypre_libraries = -lHYPRE
 
endif

# define the runtime parameter files
PROBIN_PARAMETERS := ../_parameters ./_parameters.wdconvect


# include the MAESTRO build stuff
include $(FPARALLEL)/MAESTRO/GMaestro.mak


