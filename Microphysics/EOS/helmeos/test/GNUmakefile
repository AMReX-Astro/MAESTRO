NDEBUG :=  
NDEBUG := t
MPI    :=
OMP    :=

COMP := Intel

MKVERBOSE := t

FPARALLEL := ../../..

include $(FPARALLEL)/mk/GMakedefs.mak

f90sources += testeos.f90

include $(FPARALLEL)/boxlib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib

include $(FPARALLEL)/extern/helmeos/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/helmeos
FINCLUDE_LOCATIONS += $(FPARALLEL)/extern/helmeos

include $(FPARALLEL)/extern/networks/ignition/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/networks/ignition

include $(FPARALLEL)/extern/VODE/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/VODE

testeos.$(suf).exe: $(objects)
	$(LINK.f90) -o testeos.$(suf).exe $(objects) $(libraries)

include $(FPARALLEL)/mk/GMakerules.mak

