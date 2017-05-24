NDEBUG := t
MPI    :=
OMP    := t

COMP := gfortran

MKVERBOSE := t 

NEED_VODE := t
NEED_BLAS := t

EOS_DIR := helmeos
NETWORK_DIR := triple_alpha_plus_cago

# define the location of the MAESTRO root directory
MAESTRO_TOP_DIR := $(MAESTRO_HOME)

NETWORK_DIR := triple_alpha_plus_cago

# other packages needed for data_processing
#Fmdirs := 

# directories containing files that are 'include'-d via Fortran
#Fmincludes :=

#ifdef NEED_EOS_NETWORK
#  EOS_TOP_DIR := $(MICROPHYSICS_HOME)/EOS
#  EOS_DIR := helmholtz
#
#  NETWORK_TOP_DIR := $(MAESTRO_TOP_DIR)/Microphysics/networks
#  NETWORK_DIR := triple_alpha_plus_cago
#
#
#  #MICROPHYS_CORE := $(MAESTRO_TOP_DIR)/Microphysics/EOS
#  MICROPHYS_CORE := $(EOS_TOP_DIR)
#
#  MICROPHYS_CORE += $(EOS_TOP_DIR)/$(EOS_DIR) \
#                    $(NETWORK_TOP_DIR)/$(NETWORK_DIR)
#
#  ifeq ($(findstring helmholtz, $(EOS_DIR)), helmholtz)
#    Fmincludes += Microphysics/EOS/helmholtz
#  endif
#
#  f90sources += probin.f90
#
#  include $(NETWORK_TOP_DIR)/$(strip $(NETWORK_DIR))/NETWORK_REQUIRES
#
#  ifdef NEED_VODE
#    Fmdirs += Util/VODE
#  endif
#
#  ifdef NEED_BLAS
#    Fmdirs += Util/BLAS
#  endif
#endif

#Maestro stuff
#Fmpack := $(foreach dir, $(Fmdirs), $(MAESTRO_TOP_DIR)/$(dir)/GPackage.mak)
#Fmlocs := $(foreach dir, $(Fmdirs), $(MAESTRO_TOP_DIR)/$(dir))
#
##AMReX stuff
#Fmpack += $(foreach dir, $(AMREX_CORE), $(AMREX_HOME)/$(dir)/GPackage.mak)
#Fmlocs += $(foreach dir, $(AMREX_CORE), $(AMREX_HOME)/$(dir))
#
##Microphysics
#Fmpack += $(foreach dir, $(MICROPHYS_CORE), $(dir)/GPackage.mak)
#Fmlocs += $(foreach dir, $(MICROPHYS_CORE), $(dir))

# any include files
#Fmincs := $(foreach dir, $(Fmincludes), $(MAESTRO_TOP_DIR)/$(dir))

# include the necessary GPackage.mak files that define this setup
#include $(Fmpack)

#Add module
#f90sources += fsubchandra_mod.f90

# vpath defines the directories to search for the source files
#VPATH_LOCATIONS += $(Fmlocs)

# list of directories to put in the Fortran include path
#FINCLUDE_LOCATIONS += $(Fmincs)

#programs += fsubchandra
POST_SOURCES := fsubchandra_mod.f90
POST_SOURCES += fsubchandra.f90

#all: $(pnames)

include ../GPostprocessing.mak

ALL: fsubchandra.$(suf).exe

# probin stuff
#PROBIN_PARAMETER_DIRS :=
#EXTERN_PARAMETER_DIRS += $(MICROPHYS_CORE)
#EXTERN_PARAMETERS := $(shell $(AMREX_HOME)/Tools/F_scripts/findparams.py $(EXTERN_PARAMETER_DIRS))
#
#PROBIN_TEMPLATE := $(MAESTRO_TOP_DIR)/Util/parameters/dummy.probin.template
#
#probin.f90: $(PROBIN_PARAMETERS) $(EXTERN_PARAMETERS) $(PROBIN_TEMPLATE)
#	@echo " "
#	$(AMREX_HOME)/Tools/F_scripts/write_probin.py \
#           -t $(PROBIN_TEMPLATE) -o probin.f90 -n probin \
#           --pa "$(PROBIN_PARAMETERS)" --pb "$(EXTERN_PARAMETERS)"
#
#include $(AMREX_HOME)/Tools/F_mk/GMakerules.mak
#
#
#%.$(suf).exe:%.f90 $(objects)
#ifdef MKVERBOSE
#	$(LINK.f90) -o $@ $< $(objects) $(libraries)
#else	
#	@echo "Linking $@ ... "
#	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
#endif
