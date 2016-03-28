FCOMP ?= PGI-titan
ACC ?= t

ifeq ($(FCOMP),PGI-local)
  #For local workstation installations of PGI compilers
  F90     = pgf95

  ifdef ACC
    FFLAGS  = -module build -Ibuild -acc -Minfo=acc -Mcuda=cuda7.0 -ta=nvidia:maxwell
   
    #You can use this set of flags to play with PGI's beta implementation of managed memory on NVIDIA cards	
	 #FFLAGS  = -module build -Ibuild -acc -Minfo=acc -Mcuda=cuda7.0 -ta=nvidia:maxwell,managed
  else
    FFLAGS  = -module build -Ibuild -g -O0
  endif
else ifeq ($(FCOMP),PGI-titan)
  #For compiling with PGI on OLCF's Titan
  #Note, this has been tested with the following modules:
  #   PrgEnv-pgi (the default version loaded)
  #   pgi/16.1.lustre  (from Dave Norton's set of PGI modules, add `/autofs/nccs-svm1_home1/norton/.modules` to $MODULEPATH
  #   cudatoolkit
  F90     = ftn
  
  ifdef ACC
	 #Note that the appropriate -ta should be implied in the `ftn` wrapper on Titan
    FFLAGS  = -module build -Ibuild -acc -Minfo=acc
  else
    FFLAGS  = -module build -Ibuild -g -O0 -noacc
  endif
else ifeq ($(FCOMP),GNU)
  F90     = gfortran
  FFLAGS  = -Ibuild -Jbuild -g -Wall -Wno-unused-dummy-argument -O0 -fbounds-check -fbacktrace -Wuninitialized -Wunused -ffpe-trap=invalid,zero,overflow,underflow -finit-real=snan

else ifeq ($(FCOMP),Cray)
  F90     = ftn
  ifdef ACC
    FFLAGS  = -Ibuild -Jbuild -h msgs -h acc -lcudart
  else
    FFLAGS  = -Ibuild -Jbuild -h msgs -h noacc
  endif
else
  $(error ERROR: compiler $(FCOMP) invalid)
endif

all: t1.exe

#
# rules
#

%.exe: %.f90 build/bl_types.o build/bl_constants.o build/vddot.o build/idamax.o build/dscal.o build/daxpy.o build/dgefa.o build/dgesl.o build/bdf.o 
	$(F90) $(FFLAGS) $^ -o $@

build/%.o: %.f
	@mkdir -p build
	$(F90) -c $(FFLAGS) $^ -o $@

build/%.o: %.f90
	@mkdir -p build
	$(F90) -c $(FFLAGS) $^ -o $@

clean: 
	rm -rf build
	rm -f *.o
	rm -f t1.exe
