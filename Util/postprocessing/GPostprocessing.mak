# A set of useful macros for putting together postprocessing routines

# include the main Makefile stuff
include $(AMREX_HOME)/Tools/F_mk/GMakedefs.mak


#-----------------------------------------------------------------------------
# core AMReX directories
AMREX_CORE := Src/F_BaseLib


#-----------------------------------------------------------------------------
# MAESTRO directories needed
Fmdirs := Microphysics/EOS 

# the helmeos has an include file -- also add a target to link the table
# into the problem directory.
ifeq ($(findstring helmeos, $(EOS_DIR)), helmeos)
  EOS_DIR := helmholtz
  EOS_TOP_DIR := $(MICROPHYSICS_HOME)/EOS
  EOS_PATH := $(EOS_TOP_DIR)/helmholtz
  ALL: table
endif

table:
	@if [ ! -f helm_table.dat ]; then echo ${bold}Linking helm_table.dat${normal}; ln -s $(EOS_PATH)/helm_table.dat .;  fi

ifeq ($(findstring multigamma, $(EOS_DIR)), multigamma)
  EOS_TOP_DIR := $(MICROPHYSICS_HOME)/EOS
endif

MICROPHYS_CORE := $(MAESTRO_TOP_DIR)/Microphysics/EOS $(MAESTRO_TOP_DIR)/Microphysics/screening

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

ifndef CONDUCTIVITY_DIR
  CONDUCTIVITY_DIR = constant
endif

# add in the network, EOS, and conductivity
MICROPHYS_CORE += $(EOS_TOP_DIR)/$(EOS_DIR) \
                  $(NETWORK_TOP_DIR)/$(NETWORK_DIR) \
                  $(CONDUCTIVITY_TOP_DIR)/$(CONDUCTIVITY_DIR) 

# get any additional network dependencies
include $(NETWORK_TOP_DIR)/$(strip $(NETWORK_DIR))/NETWORK_REQUIRES

ifdef NEED_VODE
  Fmdirs += Util/VODE Util/LINPACK Util/BLAS
endif

ifdef NEED_BLAS
  Fmdirs += Util/BLAS
endif

ifdef NEED_LINPACK
  UTIL_CORE += Util/LINPACK
endif

ifdef NEED_VBDF
  UTIL_CORE += Util/VBDF
endif

Fmdirs += Util/model_parser Util/simple_log


# explicitly add in any source defined in the build directory
f90sources += $(POST_SOURCES)


#-----------------------------------------------------------------------------
# core AMReX directories
Fmpack := $(foreach dir, $(AMREX_CORE), $(AMREX_HOME)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(AMREX_CORE), $(AMREX_HOME)/$(dir))
Fmincs :=

# auxillary directories
Fmpack += $(foreach dir, $(Fmdirs), $(MAESTRO_TOP_DIR)/$(dir)/GPackage.mak)
Fmpack += $(foreach dir, $(MICROPHYS_CORE), $(dir)/GPackage.mak)

Fmlocs += $(foreach dir, $(Fmdirs), $(MAESTRO_TOP_DIR)/$(dir))
Fmlocs += $(foreach dir, $(MICROPHYS_CORE), $(dir))

Fmincs += $(foreach dir, $(Fmincludes), $(MAESTRO_TOP_DIR)/$(dir))


# include the necessary GPackage.mak files that define this setup
include $(Fmpack)



# we need a probin.f90, since the various microphysics routines can
# have runtime parameters
f90sources += probin.f90

PROBIN_TEMPLATE := $(MAESTRO_TOP_DIR)/Util/parameters/dummy.probin.template
PROBIN_PARAMETER_DIRS = $(MAESTRO_TOP_DIR)/Util/initial_models/
EXTERN_PARAMETER_DIRS += $(MICROPHYS_CORE) $(NETWORK_TOP_DIR)


PROBIN_PARAMETERS := $(shell $(AMREX_HOME)/Tools/F_scripts/findparams.py $(PROBIN_PARAMETER_DIRS))
EXTERN_PARAMETERS := $(shell $(AMREX_HOME)/Tools/F_scripts/findparams.py $(EXTERN_PARAMETER_DIRS))

probin.f90: $(PROBIN_PARAMETERS) $(EXTERN_PARAMETERS) $(PROBIN_TEMPLATE)
	@echo " "
	@echo "${bold}WRITING probin.f90${normal}"
	$(AMREX_HOME)/Tools/F_scripts/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
           --pa "$(PROBIN_PARAMETERS)" --pb "$(EXTERN_PARAMETERS)"
	@echo " "





# vpath defines the directories to search for the source files

#  VPATH_LOCATIONS to first search in the problem directory      
#  Note: GMakerules.mak will include '.' at the start of the
VPATH_LOCATIONS += $(Fmlocs)


# we need the MAESTRO constants
f90sources += constants_cgs.f90
VPATH_LOCATIONS += $(MAESTRO_TOP_DIR)/Source


# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)


init_1d.$(suf).exe: $(objects)
	$(LINK.f90) -o init_1d.$(suf).exe $(objects) $(libraries)
	@echo SUCCESS


# include the fParallel Makefile rules
include $(AMREX_HOME)/Tools/F_mk/GMakerules.mak


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
	$(RM) helm_table.dat
