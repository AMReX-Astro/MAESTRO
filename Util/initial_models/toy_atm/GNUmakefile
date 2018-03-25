# Makefile for the init_1d initial model generator

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
GENERAL_NET_INPUTS := $(MAESTRO_TOP_DIR)/Microphysics/networks/$(NETWORK_DIR)/triple_alpha_plus_o.net



# the local sources needed to build the initial model generator
MODEL_SOURCES := init_1d_tanh2.f90


include ../GInitModels.mak




