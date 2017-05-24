NDEBUG := t
MPI    :=
OMP    := t

COMP := gfortran

MKVERBOSE := t 

# Needed utility code
NEED_VODE := t
NEED_BLAS := t

# Microphysics configuration
EOS_DIR := helmeos
NETWORK_DIR := triple_alpha_plus_cago

# define the location of the MAESTRO root directory
MAESTRO_TOP_DIR := $(MAESTRO_HOME)

# Local sources
POST_SOURCES := fsubchandra_mod.f90
POST_SOURCES += fsubchandra.f90

include ../GPostprocessing.mak

ALL: fsubchandra.$(suf).exe
