NDEBUG := f
MPI    := t
OMP    := 

HYPRE :=

COMP := gfortran

MKVERBOSE := t

# define the location of the MAESTRO top directory
MAESTRO_TOP_DIR := $(MAESTRO_HOME)

# define the physics packages to build this problem
EOS_DIR := helmeos   
CONDUCTIVITY_DIR := timmes_stellar
NETWORK_DIR := MESA

GENERAL_NET_INPUTS := $(MAESTRO_TOP_DIR)/Microphysics/networks/$(NETWORK_DIR)/Net-Files/test_reduced4.net

# define the special directories needed to build this problem.  Note:
# we only need to include the problem's directory if there are unique
# files there (as specified in a GPackage.mak).  The problem directory
# is always placed at the start of the vpath by the GMakerules.mak.
EXTRA_DIR := 



# include the MAESTRO build stuff
include $(MAESTRO_TOP_DIR)/GMaestro.mak


