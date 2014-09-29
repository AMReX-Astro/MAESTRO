NDEBUG := t
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

# location of root directory
MAESTRO_HOME ?= ../../../../MAESTRO
MAESTRO_TOP_DIR := $(MAESTRO_HOME)


# programs to be compiled
programs += evolve.f90

EOS_DIR := helmeos
#EOS_TOP_DIR := $(ASTRODEV_DIR)/EOS
#EOS_DIR := MESA
CONDUCTIVITY_DIR := constant

NETWORK_TOP_DIR := $(ASTRODEV_DIR)/networks
NETWORK_DIR := rprox
#NETWORK_DIR := MESA
#GENERAL_NET_INPUTS := $(NETWORK_TOP_DIR)/MESA/Net-Files/rprox_extended.net

BOXLIB_ONLY := t

EXTRA_TOP_DIR := $(ASTRODEV_DIR)
EXTRA_DIR = networks/rprox/test

testAll:

include $(MAESTRO_TOP_DIR)/GMaestro.mak

testAll: $(pnames) table

table:
	@if [ ! -f helm_table.dat ]; then echo ${bold}Linking helm_table.dat${normal}; ln -s $(EOS_PATH)/helm_table.dat .;  fi


%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif

