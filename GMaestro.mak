# A set of useful macros for putting together a MAESTRO application.

# The directories listed in Fmdirs should contain a GPackage.mak,
# which specifies to the MAESTRO build system the list of files
# to build.  These directories are also put in the vpath to define
# the locations that Make will look for source files.

# The directories listed in Fmincludes contain files that are included
# in source files, and thus specified using -I in the compiler flags.

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


