NDEBUG := 
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
Fmdirs := MAESTRO/wdconvect \
          boxlib \
          mg \
          MAESTRO \
          extern/helmeos \
          extern/conducteos \
          extern/constants \
          extern/networks/ignition \
          extern/VODE

# directories containing files that are 'include'-d via Fortran
Fmincludes := extern/helmeos

# define the runtime parameter files
PROBIN_PARAMETERS := ../_parameters ./_parameters.wdconvect
PROBIN_TEMPLATE   := ../probin.template

Fmpack := $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir))
Fmincs := $(foreach dir, $(Fmincludes), $(FPARALLEL)/$(dir))


# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs)

# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)


# define the build instructions for the executable
main.$(suf).exe: $(objects)
	$(HPCLINK) $(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)
	@echo SUCCESS


probin.f90: $(PROBIN_PARAMETERS) $(PROBIN_TEMPLATE)
	../write_probin.py -t $(PROBIN_TEMPLATE) $(PROBIN_PARAMETERS)



# include the Makefile rules
include $(FPARALLEL)/mk/GMakerules.mak

