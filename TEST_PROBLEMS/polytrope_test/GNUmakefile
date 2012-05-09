NDEBUG := t
MPI    := t
OMP    := 

#COMP   := blrts_xlf90
COMP   := gfortran

MKVERBOSE := t

FPARALLEL := ../../..

Fmdirs := MAESTRO/TEST_PROBLEMS/polytrope_test

EOS_DIR          := MAESTRO/TEST_PROBLEMS/polytrope_test/idealgas
NETWORK_DIR      := extern/networks/null
CONDUCTIVITY_DIR := extern/conductivity/constant

EXTRA_DIR := extern/random

include $(FPARALLEL)/MAESTRO/GMaestro.mak
