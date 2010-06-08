# A set of useful macros for putting together a MAESTRO application.

# The directories listed in Fmdirs should contain a GPackage.mak,
# which specifies to the MAESTRO build system the list of files
# to build.  These directories are also put in the vpath to define
# the locations that Make will look for source files.

# The directories listed in Fmincludes contain files that are included
# in source files, and thus specified using -I in the compiler flags.

# core MAESTRO directories
MAESTRO_CORE := boxlib \
                mg \
                MAESTRO/Source \
                extern/constants \
                extern/LAPACK



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

deppairs: build_info.f90

build_info.f90: 
	$(FPARALLEL)/scripts/make_build_info "$(Fmdirs)" "$(COMPILE.f90)" "$(COMPILE.f)" "$(COMPILE.c)" "$(LINK.f90)"


$(odir)/build_info.o: build_info.f90
	$(COMPILE.f90) $(OUTPUT_OPTION) build_info.f90
	rm -f build_info.f90


# unsure why this method does not work right -- it gives circular dependencies
# but this would force build_info.f90 to only be rebuilt if another source
# file changed.

#$(odir)/build_info.o: $(fsources) $(f90sources) $(csources)
#       $(FPARALLEL)/scripts/make_build_info
#       $(COMPILE.f90) $(OUTPUT_OPTION) build_info.f90

print-%: ; @echo $* is $($*).
