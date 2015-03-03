# A set of useful macros for putting together a MAESTRO application.

# check the version number -- this comes from the GNU Make cookbook
NEED := 3.81
OK := $(filter $(NEED),$(firstword $(sort $(MAKE_VERSION) $(NEED))))

ifndef OK
  $(error your version of GNU make is too old.  You need atleast version $(NEED))
endif

ifeq ($(findstring ~, $(BOXLIB_HOME)), ~)
   $(error you cannot include the ~ character in your BOXLIB_HOME variable)
endif

# include the main Makefile stuff
include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

# default target (make just takes the one that appears first)
ALL: main.$(suf).exe


#-----------------------------------------------------------------------------
# core BoxLib directories
BOXLIB_CORE := Src/F_BaseLib \
               Src/LinearSolvers/F_MG

# include the random number generator stuff
RANDOM := t

#-----------------------------------------------------------------------------
# core MAESTRO directories
MAESTRO_CORE := 

# path to SDC files -- note this must come before Source/ in the vpath
SDC_CORE := 

ifdef SDC
  SDC_CORE += $(MAESTRO_TOP_DIR)/Source_SDC
endif

# next look for the files in Source/ itself 
#
#   Note: a unit test (UNIT_TEST := t) tests only a single component
#   of the MAESTRO algorithm, so we don't, in general, want to build
#   all of the source in the MAESTRO/Source directory.  So, for unit
#   tests, we leave it off the list of core directories, but do
#   include it in the VPATH 
#
#   Setting BOXLIB_ONLY := t means that we don't even want the
#   MAESTRO/Source directory in our VPATH

ifndef UNIT_TEST 
  ifndef BOXLIB_ONLY 
    MAESTRO_CORE += Source 
  endif 
endif


#-----------------------------------------------------------------------------
# core extern directories needed by every MAESTRO build
UTIL_CORE := 
ifndef BOXLIB_ONLY
  UTIL_CORE := Util/model_parser 
endif


#-----------------------------------------------------------------------------
# microphysics
MICROPHYS_CORE := $(MAESTRO_TOP_DIR)/Microphysics/EOS

# locations of the microphysics 
ifndef EOS_TOP_DIR 
  EOS_TOP_DIR := $(MAESTRO_TOP_DIR)/Microphysics/EOS
endif

ifndef NETWORK_TOP_DIR 
  NETWORK_TOP_DIR := $(MAESTRO_TOP_DIR)/Microphysics/networks
endif

ifndef CONDUCTIVITY_TOP_DIR
  CONDUCTIVITY_TOP_DIR := $(MAESTRO_TOP_DIR)/Microphysics/conductivity
endif

# add in the network, EOS, and conductivity
MICROPHYS_CORE += $(EOS_TOP_DIR)/$(EOS_DIR) \
                  $(NETWORK_TOP_DIR)/$(NETWORK_DIR) \
                  $(CONDUCTIVITY_TOP_DIR)/$(CONDUCTIVITY_DIR) 

# get any additional network dependencies
include $(NETWORK_TOP_DIR)/$(strip $(NETWORK_DIR))/NETWORK_REQUIRES

ifdef NEED_VODE
  UTIL_CORE += Util/VODE
endif

ifdef NEED_BLAS
  UTIL_CORE += Util/BLAS
endif

ifdef NEED_VBDF
  UTIL_CORE += Util/VBDF
endif

# the helmeos has an include file -- also add a target to link the table
# into the problem directory.
ifeq ($(findstring helmeos, $(EOS_DIR)), helmeos)
  Fmincludes := Microphysics/EOS/helmeos
  EOS_PATH := $(MAESTRO_TOP_DIR)/Microphysics/EOS/$(strip $(EOS_DIR))
  ALL: table
endif


table:
	@if [ ! -f helm_table.dat ]; then echo ${bold}Linking helm_table.dat${normal}; ln -s $(EOS_PATH)/helm_table.dat .;  fi


#-----------------------------------------------------------------------------
# extra directory
ifndef EXTRA_TOP_DIR 
  EXTRA_TOP_DIR := $(MAESTRO_TOP_DIR)/
endif

EXTRAS := $(addprefix $(EXTRA_TOP_DIR)/, $(EXTRA_DIR))

ifndef EXTRA_TOP_DIR2 
  EXTRA_TOP_DIR2 := $(MAESTRO_TOP_DIR)/
endif

EXTRAS += $(addprefix $(EXTRA_TOP_DIR2)/, $(EXTRA_DIR2))


#-----------------------------------------------------------------------------
# compile in support for particles
PARTICLES := t


#-----------------------------------------------------------------------------
# Fmpack is the list of all the GPackage.mak files that we need to
# include into the build system to define the list of source files.
#
# Fmlocs is the list of all the directories that we want to search
# for the source files in -- this is usually going to be the
# same as the list of directories containing GPackage.mak defined
# above.
#
# Fincs is the list of directories that have include files that
# we need to tell the compiler about.


# SDC
Fmpack := $(foreach dir, $(SDC_CORE), $(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(SDC_CORE), $(dir))

# Maestro and Util modules
Fmdirs += $(UTIL_CORE) \
          $(MAESTRO_CORE)

Fmpack += $(foreach dir, $(Fmdirs), $(MAESTRO_TOP_DIR)/$(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(Fmdirs), $(MAESTRO_TOP_DIR)/$(dir))

# Microphysics
Fmpack += $(foreach dir, $(MICROPHYS_CORE), $(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(MICROPHYS_CORE), $(dir))

# Extras
Fmpack += $(foreach dir, $(EXTRAS), $(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(EXTRAS), $(dir))

# BoxLib
Fmpack += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(BOXLIB_CORE), $(BOXLIB_HOME)/$(dir))


# any include directories
Fmincs := $(foreach dir, $(Fmincludes), $(MAESTRO_TOP_DIR)/$(dir))


# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

# vpath defines the directories to search for the source files

# we always want to search the MAESTRO/Source directory, even for
# unit tests, since they may build individual files there.
ifdef UNIT_TEST
  VPATH_LOCATIONS += $(MAESTRO_TOP_DIR)/Source
endif

#  Note: GMakerules.mak will include '.' at the start of the
#  VPATH_LOCATIONS to first search in the problem directory
VPATH_LOCATIONS += $(Fmlocs)  $(EXTRA_LOCATIONS)


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
ifndef BOXLIB_ONLY
  PROBIN_TEMPLATE := $(MAESTRO_TOP_DIR)/Source/probin.template
else
  PROBIN_TEMPLATE := $(MAESTRO_TOP_DIR)/Util/parameters/dummy.probin.template
endif

# list of the directories to search for _parameters files
PROBIN_PARAMETER_DIRS = ./ 

ifndef BOXLIB_ONLY 
  PROBIN_PARAMETER_DIRS += $(MAESTRO_TOP_DIR)/Source
endif

# list of all valid _parameters files for probin
PROBIN_PARAMETERS := $(shell $(BOXLIB_HOME)/Tools/F_scripts/findparams.py $(PROBIN_PARAMETER_DIRS))

# list of all valid _parameters files for extern
EXTERN_PARAMETER_DIRS += $(MICROPHYS_CORE)
EXTERN_PARAMETERS := $(shell $(BOXLIB_HOME)/Tools/F_scripts/findparams.py $(EXTERN_PARAMETER_DIRS))

probin.f90: $(PROBIN_PARAMETERS) $(EXTERN_PARAMETERS) $(PROBIN_TEMPLATE)
	@echo " "
	@echo "${bold}WRITING probin.f90${normal}"
	$(BOXLIB_HOME)/Tools/F_scripts/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
           --pa "$(PROBIN_PARAMETERS)" --pb "$(EXTERN_PARAMETERS)"
	@echo " "


#-----------------------------------------------------------------------------
# build_info stuff
deppairs: build_info.f90

build_info.f90: 
	@echo " "
	@echo "${bold}WRITING build_info.f90${normal}"
	$(BOXLIB_HOME)/Tools/F_scripts/makebuildinfo.py \
           --modules "$(Fmdirs) $(MICROPHYS_CORE)" \
           --FCOMP "$(COMP)" \
           --FCOMP_version "$(FCOMP_VERSION)" \
           --f90_compile_line "$(COMPILE.f90)" \
           --f_compile_line "$(COMPILE.f)" \
           --C_compile_line "$(COMPILE.c)" \
           --link_line "$(LINK.f90)" \
           --boxlib_home "$(BOXLIB_HOME)" \
           --source_home "$(MAESTRO_TOP_DIR)" \
           --extra_home "$(ASTRODEV_DIR)"
	@echo " "

$(odir)/build_info.o: build_info.f90
	$(COMPILE.f90) $(OUTPUT_OPTION) build_info.f90
	rm -f build_info.f90


#-----------------------------------------------------------------------------
# include the BoxLib Fortran Makefile rules
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




