NDEBUG :=  
NDEBUG := t
MPI    :=
OMP    :=

COMP := PathScale

MKVERBOSE := t

FPARALLEL := ../../../..

include $(FPARALLEL)/mk/GMakedefs.mak

programs += testburn

include $(FPARALLEL)/boxlib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib

include $(FPARALLEL)/extern/networks/dvode_test/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/networks/dvode_test

include $(FPARALLEL)/extern/VODE/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/VODE


all: $(pnames)

include $(FPARALLEL)/mk/GMakerules.mak

%.$(suf).exe:%.f90 $(objects)
ifdef MKVERBOSE
	$(LINK.f90) -o $@ $< $(objects) $(libraries)
else
	@echo "Linking $@ ... "
	@$(LINK.f90) -o $@ $< $(objects) $(libraries)
endif



