NDEBUG :=  
NDEBUG := 
MPI    :=
OMP    :=

COMP := Intel

MKVERBOSE := t

# define the location of the fParallel root directory
FPARALLEL := ../../../..


# include the main Makefile stuff
include $(FPARALLEL)/mk/GMakedefs.mak


# define the programs we want to build
programs += testburn
programs += eval
programs += testjacobian

# define the packages to build this problem
Fmdirs := boxlib \
          extern/helmeos \
          extern/networks/ignition_chamulak \
          extern/constants \
          extern/VODE

Fmpack := $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir))

# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs)

FINCLUDE_LOCATIONS += $(Fmlocs)


all: $(pnames)


# include the Makefile rules
include $(FPARALLEL)/mk/GMakerules.mak


# define the build instructions for the executable
%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif



