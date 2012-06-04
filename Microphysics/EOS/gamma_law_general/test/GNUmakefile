NDEBUG :=  
NDEBUG := t
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

# define the location of the MAESTRO root directory
MAESTRO_TOP_DIR := ../../../..


# define the physics packages to build this problem
EOS_DIR := gamma_law_general
CONDUCTIVITY_DIR := constant

NETWORK_DIR := general_null
GENERAL_NET_INPUTS := $(MAESTRO_TOP_DIR)/Microphysics/networks/$(NETWORK_DIR)/ignition.net

EXTRA_DIR := Microphysics/EOS/gamma_law_general/test

UNIT_TEST := t


include $(MAESTRO_TOP_DIR)/GMaestro.mak


