NDEBUG :=  
MPI    := t
OMP    :=

HYPRE :=

COMP := gfortran

MKVERBOSE := t

# define the location of the MAESTRO top directory
MAESTRO_TOP_DIR := ../../..

# define the physics packages to build this problem
EOS_DIR := helmeos   
CONDUCTIVITY_DIR := timmes_stellar
NETWORK_DIR := ignition_chamulak 

# define the special directories needed to build this problem.  Note:
# we only need to include the problem's directory if there are unique
# files there (as specified in a GPackage.mak).  The problem directory
# is always placed at the start of the vpath by the GMakerules.mak.
EXTRA_DIR := 

ifdef HYPRE

# on manda
  HYPRE_DIR = ../../extern/hypre

# on hopper
# HYPRE_DIR =  /usr/common/acts/Hypre/2.6.0b/src/cray-xt5_O

  LDFLAGS += -L $(HYPRE_DIR)/lib
  FINCLUDE_LOCATIONS += $(HYPRE_DIR)/include
  hypre_libraries = -lHYPRE

endif


# include the MAESTRO build stuff
include $(MAESTRO_TOP_DIR)/GMaestro.mak


