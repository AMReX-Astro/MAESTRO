NDEBUG := 
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

SDC    :=

# location of top-level MAESTRO directory
MAESTRO_TOP_DIR = ../../../..


# EOS and network
EOS_DIR := helmeos

NETWORK_DIR := general_null
GENERAL_NET_INPUTS := $(MAESTRO_TOP_DIR)/Microphysics/networks/$(NETWORK_DIR)/ignition_wdconvect_urca.net


# the local sources needed to build the initial model generator
MODEL_SOURCES := init_1d.f90
MODEL_SOURCES += urca_composition.f90


include ../../GInitModels.mak


