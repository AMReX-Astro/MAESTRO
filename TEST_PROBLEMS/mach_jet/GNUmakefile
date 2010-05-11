NDEBUG := t
MPI    :=
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
          MAESTRO \
          extern/EOS/gamma_law_general \
          extern/EOS/conducteos \
          extern/constants \
          extern/networks/gammalaw \
          extern/VODE

# directories containing files that are 'include'-d via Fortran
Fmincludes := extern/EOS/helmeos

# define the runtime parameter files
PROBIN_PARAMETERS := ../_parameters


# include the MAESTRO build stuff
include $(FPARALLEL)/MAESTRO/GMaestro.mak

# include the Makefile rules
include $(FPARALLEL)/mk/GMakerules.mak

