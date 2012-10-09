NDEBUG := t
MPI    := t
OMP    := 

#COMP   := blrts_xlf90
COMP   := gfortran

MKVERBOSE := t

MAESTRO_TOP_DIR := ../..

EOS_TOP_DIR      := $(MAESTRO_TOP_DIR)/TEST_PROBLEMS/polytrope_test
EOS_DIR          := idealgas

NETWORK_DIR      := null
CONDUCTIVITY_DIR := constant

EXTRA_DIR := Util/random \
             TEST_PROBLEMS/polytrope_test

include $(MAESTRO_TOP_DIR)/GMaestro.mak
