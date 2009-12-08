NDEBUG := t
MPI    := t
OMP    :=

COMP := Intel

MKVERBOSE := t

# define the location of the fParallel root directory
FPARALLEL := ../..


# include the main Makefile stuff
include $(FPARALLEL)/mk/GMakedefs.mak


# define the packages to build this problem
#   note: we only need to include the problem's directory if there are 
#   unique files there (as specified in a GPackage.mak).  The problem 
#   directory is always placed at the start of the vpath by the 
#   GMakerules.mak.
Fmdirs := boxlib \
          mg \
	  extern/EOS/gamma_law_general \
	  MAESTRO/test_diffusion \
          extern/constants \
          extern/networks/triple_alpha \
          extern/VODE

# we dont' want to include EVERYTHING in the MAESTRO directory, but
# we should add it to the VPATH for source file searching
VPATH_LOCATIONS += $(FPARALLEL)/MAESTRO/

# define the runtime parameter files
PROBIN_PARAMETERS := ../_parameters ./_parameters

# include the MAESTRO build stuff
include $(FPARALLEL)/MAESTRO/GMaestro.mak

# include the Makefile rules
include $(FPARALLEL)/mk/GMakerules.mak

