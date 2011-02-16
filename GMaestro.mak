# A set of useful macros for putting together a MAESTRO application.

# include the main Makefile stuff
include $(FPARALLEL)/mk/GMakedefs.mak


# core MAESTRO directories
MAESTRO_CORE := boxlib \
                mg \
                extern/constants \
		extern/model_parser \
                extern/LAPACK \
                extern/random

# a unit test tests only a single component of the MAESTRO algorithm,
# so we don't, in general, want to build all of the source in the
# MAESTRO/Source directory.  So, for unit tests, we leave it off the
# list of core directories
ifndef UNIT_TEST
       MAESTRO_CORE += MAESTRO/Source
endif


# compile in support for particles
PARTICLES := t

# The directories listed in Fmdirs should contain a GPackage.mak,
# which specifies to the MAESTRO build system the list of files
# to build.  These directories are also put in the vpath to define
# the locations that Make will look for source files.

# The directories listed in Fmincludes contain files that are included
# in source files, and thus specified using -I in the compiler flags.

# add in the problem specific stuff ("extras"), EOS, network, and
# conductivity
Fmdirs += $(EXTRA_DIR) \
          $(EOS_DIR) \
          $(NETWORK_DIR) \
          $(CONDUCTIVITY_DIR) \
          $(MAESTRO_CORE)


# networks in general need the VODE + other packages
ifneq ($(findstring null, $(NETWORK_DIR)), null)
	Fmdirs += extern/VODE 
endif


# the helmeos has an include file
ifeq ($(findstring helmeos, $(EOS_DIR)), helmeos)
	Fmincludes := extern/EOS/helmeos
endif


Fmpack := $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(Fmdirs), $(FPARALLEL)/$(dir))
Fmincs := $(foreach dir, $(Fmincludes), $(FPARALLEL)/$(dir))


# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs)

# we always want to search the MAESTRO/Source directory, even for
# unit tests, since they may build individual files there.
ifdef UNIT_TEST
      VPATH_LOCATIONS += $(FPARALLEL)/MAESTRO/Source
endif



# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)


# define the build instructions for the executable
main.$(suf).exe: $(objects)
	$(HPCLINK) $(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)
	@echo SUCCESS

# runtime parameter stuff
PROBIN_TEMPLATE := $(FPARALLEL)/MAESTRO/probin.template

probin.f90: $(PROBIN_PARAMETERS) $(PROBIN_TEMPLATE)
	$(FPARALLEL)/MAESTRO/write_probin.py -t $(PROBIN_TEMPLATE) $(PROBIN_PARAMETERS)


# build_info stuff
deppairs: build_info.f90

build_info.f90: 
	$(FPARALLEL)/scripts/make_build_info \
            "$(Fmdirs)" "$(COMPILE.f90)" "$(COMPILE.f)" \
            "$(COMPILE.c)" "$(LINK.f90)"


$(odir)/build_info.o: build_info.f90
	$(COMPILE.f90) $(OUTPUT_OPTION) build_info.f90
	rm -f build_info.f90


# include the fParallel Makefile rules
include $(FPARALLEL)/mk/GMakerules.mak


# for debugging, comment out this line, and do, eg, "make
# print-Fmlocs" to have make print out the contents of 
# the Fmlocs variable
print-%: ; @echo $* is $($*)
