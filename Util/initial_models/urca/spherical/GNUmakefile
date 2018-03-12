NDEBUG := 
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t

SDC    :=

# location of top-level MAESTRO directory
MAESTRO_TOP_DIR = ../../../..

# EOS and network
EOS_TOP_DIR := $(MICROPHYSICS_HOME)/EOS
EOS_DIR := helmholtz
CONDUCTIVITY_DIR := stellar
NETWORK_TOP_DIR := $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := ignition_reaclib/URCA-simple

# the local sources needed to build the initial model generator
MODEL_SOURCES := init_1d.f90
MODEL_SOURCES += urca_composition.f90

include GInitModels.mak
