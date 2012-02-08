# A set of useful macros for putting together a MAESTRO application.

# include the main Makefile stuff
include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

# default target (make just takes the one that appears first)
ALL: main.$(suf).exe


#-----------------------------------------------------------------------------
# core BoxLib directories
BOXLIB_CORE := Src/F_BaseLib \
               Src/LinearSolvers/F_MG


#-----------------------------------------------------------------------------
# core MAESTRO directories
MAESTRO_CORE := 

# if we are doing the SDC algorithm, then look first for source
# files living in Source_SDC/ then look in Source/
ifdef SDC
  MAESTRO_CORE += MAESTRO/Source_SDC
endif

# next look for the files in Source/ itself
#   Note: a unit test tests only a single component of the MAESTRO
#   algorithm, so we don't, in general, want to build all of the
#   source in the MAESTRO/Source directory.  So, for unit tests, we
#   leave it off the list of core directories 
ifndef UNIT_TEST
  MAESTRO_CORE += MAESTRO/Source 
endif


#-----------------------------------------------------------------------------
# core extern directories needed by every MAESTRO build
EXTERN_CORE := extern/constants \
	       extern/model_parser \
               extern/LAPACK \
               extern/random \
               extern/BLAS 

# add in the network, EOS, and conductivity
EXTERN_CORE += $(EOS_DIR) \
               $(NETWORK_DIR) \
               $(CONDUCTIVITY_DIR) 

# networks in general need the VODE 
ifneq ($(findstring null, $(NETWORK_DIR)), null)
  EXTERN_CORE += extern/VODE 
endif


#-----------------------------------------------------------------------------
# compile in support for particles
PARTICLES := t


#-----------------------------------------------------------------------------
# The directories listed in Fmdirs should contain a GPackage.mak,
# which specifies to the MAESTRO build system the list of files
# to build.  These directories are also put in the vpath to define
# the locations that Make will look for source files.

# The directories listed in Fmincludes contain files that are included
# in source files, and thus specified using -I in the compiler flags.

Fmdirs += $(EXTRA_DIR) \
          $(EXTERN_CORE) \
          $(MAESTRO_CORE)


# the helmeos has an include file
ifeq ($(findstring helmeos, $(EOS_DIR)), helmeos)
  Fmincludes := extern/EOS/helmeos
endif


Fmpack := $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir))
Fmincs := $(foreach dir, $(Fmincludes), $(FPARALLEL)/$(dir))

Fmpack += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir))




# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files

#  Note: GMakerules.mak will include '.' at the start of the
#  VPATH_LOCATIONS to first search in the problem directory
VPATH_LOCATIONS += $(Fmlocs)

# we always want to search the MAESTRO/Source directory, even for
# unit tests, since they may build individual files there.
ifdef UNIT_TEST
  VPATH_LOCATIONS += $(FPARALLEL)/MAESTRO/Source
endif



# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)


#-----------------------------------------------------------------------------
# define the build instructions for the executable
main.$(suf).exe: $(objects)
	$(HPCLINK) $(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)
	@echo SUCCESS


#-----------------------------------------------------------------------------
# runtime parameter stuff (probin.f90)

# template used by write_probin.py to build probin.f90
PROBIN_TEMPLATE := $(FPARALLEL)/MAESTRO/probin.template

# list of the directories to search for _parameters files
PROBIN_PARAMETER_DIRS = ./ ../ ../../

# list of all valid _parameters files for probin
PROBIN_PARAMETERS := $(shell $(BOXLIB_HOME)/Tools/F_scripts/findparams.py $(PROBIN_PARAMETER_DIRS))

# list of all valid _parameters files for extern
EXTERN_PARAMETER_DIRS += $(foreach dir, $(EXTERN_CORE), $(FPARALLEL)/$(dir))
EXTERN_PARAMETERS := $(shell $(BOXLIB_HOME)/Tools/F_scripts/findparams.py $(EXTERN_PARAMETER_DIRS))

probin.f90: $(PROBIN_PARAMETERS) $(EXTERN_PARAMETERS) $(PROBIN_TEMPLATE)
	@echo "---------------------------------------------------------------------------"
	@echo "WRITING probin.f90"
	$(BOXLIB_HOME)/Tools/F_scripts/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
           --pa "$(PROBIN_PARAMETERS)" --pb "$(EXTERN_PARAMETERS)"
	@echo "---------------------------------------------------------------------------"
	@echo " "


#-----------------------------------------------------------------------------
# build_info stuff
deppairs: build_info.f90

build_info.f90: 
	@echo "---------------------------------------------------------------------------"
	@echo "WRITING build_info.f90"
	$(BOXLIB_HOME)/Tools/F_scripts/make_build_info \
            "$(Fmdirs)" "$(COMPILE.f90)" "$(COMPILE.f)" \
            "$(COMPILE.c)" "$(LINK.f90)" "$(BOXLIB_HOME)"
	@echo "---------------------------------------------------------------------------"
	@echo " "

$(odir)/build_info.o: build_info.f90
	$(COMPILE.f90) $(OUTPUT_OPTION) build_info.f90
	rm -f build_info.f90


#-----------------------------------------------------------------------------
# include the fParallel Makefile rules
include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak


#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. Fmlocs, simply do "make print-Fmlocs".  This will
# print out the value.
print-%: ; @echo $* is $($*)


#-----------------------------------------------------------------------------
# cleaning.  Add more actions to 'clean' and 'realclean' to remove 
# probin.f90 and build_info.f90 -- this is where the '::' in make comes
# in handy
clean:: 
	$(RM) probin.f90 
	$(RM) build_info.f90




