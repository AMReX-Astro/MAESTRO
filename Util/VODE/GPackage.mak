fsources += dvode.f

fsources += dvhin.f
fsources += dvindy.f
fsources += dvjac.f
fsources += dvjust.f
fsources += dvnlsd.f
fsources += dvnorm.f
fsources += dvset.f
fsources += dvsol.f
# dvsrco.f
fsources += dvstep.f
fsources += xerrwd.f
fsources += dewset.f
fsources += ixsav.f
fsources += dumach.f
fsources += iumach.f

fsources += dacopy.f
fsources += dgbfa.f
fsources += dgbsl.f
fsources += dgefa.f
fsources += dgesl.f
# xsetf.f
# xsetun.f

include $(FPARALLEL)/extern/BLAS/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/BLAS

