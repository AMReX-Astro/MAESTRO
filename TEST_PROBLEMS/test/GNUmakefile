NDEBUG :=  
NDEBUG := t
MPI    :=
OMP    :=

COMP := Intel

MKVERBOSE := t

FPARALLEL := ../..

include $(FPARALLEL)/mk/GMakedefs.mak

f90sources += main.f90

include $(FPARALLEL)/boxlib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib

include $(FPARALLEL)/mg/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/mg

include $(FPARALLEL)/nova_varden_heating/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/nova_varden_heating

include ../../extern/helmeos/GPackage.mak
VPATH_LOCATIONS += ../../extern/helmeos
FINCLUDE_LOCATIONS += ../../extern/helmeos

include $(FPARALLEL)/varden/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/varden

include $(FPARALLEL)/extern/networks/ignition/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/networks/ignition

main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(FPARALLEL)/mk/GMakerules.mak

