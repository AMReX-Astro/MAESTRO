NDEBUG := 
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t


# location of top-level MAESTRO directory
MAESTRO_TOP_DIR = ../../..


# EOS and network
EOS_TOP_DIR := $(MICROPHYSICS_HOME)/EOS
EOS_DIR := helmholtz

CONDUCTIVITY_DIR := stellar
NETWORK_TOP_DIR := $(MICROPHYSICS_HOME)/networks
NETWORK_DIR := ignition_reaclib/URCA-simple

# the local sources needed to build
POST_SOURCES := central_angle_average.f90

include GPostprocessing.mak

ALL: central_angle_average.$(suf).exe
