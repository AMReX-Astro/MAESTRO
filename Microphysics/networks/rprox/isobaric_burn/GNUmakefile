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
CONDUCTIVITY_DIR := constant
NETWORK_DIR := rprox

BOXLIB_ONLY := t

EXTRA_DIR = Microphysics/networks/rprox/test

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

