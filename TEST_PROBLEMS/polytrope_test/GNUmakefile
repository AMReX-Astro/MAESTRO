NDEBUG := t
MPI    := t
OMP    := 

#COMP   := blrts_xlf90
COMP   := gfortran

MKVERBOSE := t

FPARALLEL := ../..

Fmdirs := MAESTRO/polytrope_test_new

EOS_DIR := MAESTRO/polytrope_test/idealgas
CONDUCTIVITY_DIR := extern/conductivity/constant
NETWORK_DIR := extern/networks/null

EXTRA_DIR := extern/random

include $(FPARALLEL)/MAESTRO/GMaestro.mak
