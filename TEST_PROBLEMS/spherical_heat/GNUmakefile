NDEBUG :=  
MPI    :=  

NDEBUG := t
MPI    := t

OMP    :=

COMP := Intel

MKVERBOSE := t

FPARALLEL := ../..

include $(FPARALLEL)/mk/GMakedefs.mak

include $(FPARALLEL)/spherical_heat/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/spherical_heat

include $(FPARALLEL)/boxlib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib

include $(FPARALLEL)/mg/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/mg

include $(FPARALLEL)/MAESTRO/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/MAESTRO

include $(FPARALLEL)/extern/helmeos/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/helmeos
FINCLUDE_LOCATIONS += $(FPARALLEL)/extern/helmeos

include $(FPARALLEL)/extern/conducteos/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/conducteos
FINCLUDE_LOCATIONS += $(FPARALLEL)/extern/conducteos

include $(FPARALLEL)/extern/constants/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/constants

include $(FPARALLEL)/extern/networks/null/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/networks/null

include $(FPARALLEL)/extern/VODE/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/VODE

main.$(suf).exe: $(objects)
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)
	@echo SUCCESS

include $(FPARALLEL)/mk/GMakerules.mak

