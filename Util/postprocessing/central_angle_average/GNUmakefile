NDEBUG := 
MPI    :=
OMP    :=

COMP := gfortran

MKVERBOSE := t


# location of top-level MAESTRO directory
MAESTRO_TOP_DIR = ../../..


# EOS and network
EOS_DIR := helmeos

NETWORK_DIR := general_null
GENERAL_NET_INPUTS := $(MAESTRO_TOP_DIR)/Microphysics/networks/$(NETWORK_DIR)/ignition_wdconvect_urca.net

# the local sources needed to build
POST_SOURCES := central_angle_average.f90

include ../GPostprocessing.mak

ALL: central_angle_average.$(suf).exe
