NDEBUG := 
NDEBUG := t
MPI    :=
OMP    :=

COMP := Intel

MKVERBOSE := t

FPARALLEL := ../..

include $(FPARALLEL)/mk/GMakedefs.mak

VPATH_LOCATIONS += $(FPARALLEL)/MAESTRO/test_eta_stuff

include $(FPARALLEL)/boxlib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib

include $(FPARALLEL)/mg/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/mg

include $(FPARALLEL)/MAESTRO/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/MAESTRO

include $(FPARALLEL)/extern/EOS/helmeos/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/EOS/helmeos
FINCLUDE_LOCATIONS += $(FPARALLEL)/extern/EOS/helmeos

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


PROBIN_PARAMETERS := ../_parameters
PROBIN_TEMPLATE   := ../probin.template

probin.f90: $(PROBIN_PARAMETERS) $(PROBIN_TEMPLATE)
	../write_probin.py -t $(PROBIN_TEMPLATE) $(PROBIN_PARAMETERS)


include $(FPARALLEL)/mk/GMakerules.mak

