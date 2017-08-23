NDEBUG := t
MPI    :=
OMP    :=

MKVERBOSE := t

COMP := gfortran

# some routines need an eos/network (i.e. to compute thermodynamic
# quantities.  If that is the case, set NEED_EOS_NETWORK := t
NEED_EOS_NETWORK := t

# include linear algebra packages
NEED_BLAS := t
NEED_LINPACK := t

# define the location of the MAESTRO root directory
MAESTRO_TOP_DIR := $(MAESTRO_HOME)

# include the main Makefile stuff
include $(AMREX_HOME)/Tools/F_mk/GMakedefs.mak

# core AMReX directories
AMREX_CORE := Src/F_BaseLib

# other packages needed for data_processing
Fmdirs :=

# directories containing files that are 'include'-d via Fortran
Fmincludes :=

ifdef NEED_EOS_NETWORK
  EOS_TOP_DIR := $(MICROPHYSICS_HOME)/EOS
  EOS_DIR := helmholtz

  NETWORK_TOP_DIR := $(MAESTRO_TOP_DIR)/Microphysics/networks
  NETWORK_DIR := triple_alpha_plus_cago

  MICROPHYS_CORE += $(MAESTRO_HOME)/Microphysics/EOS \
                    $(EOS_TOP_DIR)/$(EOS_DIR) \
                    $(NETWORK_TOP_DIR)/$(NETWORK_DIR) \
                    $(MAESTRO_TOP_DIR)/Microphysics/screening 

  ifeq ($(findstring helmholtz, $(EOS_DIR)), helmholtz)
    EOS_PATH := $(MICROPHYSICS_HOME)/EOS/$(strip $(EOS_DIR))
    all: table
  endif

  f90sources += probin.f90

  include $(NETWORK_TOP_DIR)/$(strip $(NETWORK_DIR))/NETWORK_REQUIRES

  ifdef NEED_VODE
    Fmdirs += Util/VODE
  endif

  ifdef NEED_BLAS
    Fmdirs += Util/BLAS
  endif

  ifdef NEED_LINPACK
    Fmdirs += Util/LINPACK
  endif
endif

# Maestro stuff
Fmpack := $(foreach dir, $(Fmdirs), $(MAESTRO_TOP_DIR)/$(dir)/GPackage.mak)
Fmlocs := $(foreach dir, $(Fmdirs), $(MAESTRO_TOP_DIR)/$(dir))

# AMReX stuff
Fmpack += $(foreach dir, $(AMREX_CORE), $(AMREX_HOME)/$(dir)/GPackage.mak)
Fmlocs += $(foreach dir, $(AMREX_CORE), $(AMREX_HOME)/$(dir))

# Microphysics
Fmpack += $(foreach dir, $(MICROPHYS_CORE), $(dir)/GPackage.mak)
Fmpack += $(MICROPHYS_CORE/GPackage.mak)
Fmlocs += $(foreach dir, $(MICROPHYS_CORE), $(dir))

# any include files
Fmincs := $(foreach dir, $(Fmincludes), $(MAESTRO_TOP_DIR)/$(dir))

# include the necessary GPackage.mak files that define this setup
include $(Fmpack)

#Add module
f90sources += fsubchandra_mod.f90

# vpath defines the directories to search for the source files
VPATH_LOCATIONS += $(Fmlocs)

# list of directories to put in the Fortran include path
FINCLUDE_LOCATIONS += $(Fmincs)

programs += fsubchandra

all: $(pnames) table

# probin stuff
PROBIN_PARAMETER_DIRS := 
PROBIN_PARAMETERS := $(shell $(AMREX_HOME)/Tools/F_scripts/findparams.py $(PROBIN_PARAMETER_DIRS))

EXTERN_PARAMETER_DIRS += $(MICROPHYS_CORE) $(MAESTRO_TOP_DIR)/Microphysics/networks
EXTERN_PARAMETERS := $(shell $(AMREX_HOME)/Tools/F_scripts/findparams.py $(EXTERN_PARAMETER_DIRS))

PROBIN_TEMPLATE := $(MAESTRO_TOP_DIR)/Util/parameters/dummy.probin.template

probin.f90: $(PROBIN_PARAMETERS) $(EXTERN_PARAMETERS) $(PROBIN_TEMPLATE)
	@echo " "
	$(AMREX_HOME)/Tools/F_scripts/write_probin.py \
           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
           --pa "$(PROBIN_PARAMETERS)" --pb "$(EXTERN_PARAMETERS)"

include $(AMREX_HOME)/Tools/F_mk/GMakerules.mak

%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif

table:
	@if [ ! -f helm_table.dat ]; then echo ${bold}Linking helm_table.dat${normal}; ln -s $(EOS_PATH)/helm_table.dat .;  fi



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

realclean ::
	$(RM) helm_table.dat
