# For reference, below is the original GPackage.mak shipped with
# BDF in the Combustion codebase.

# For F90 BoxLib based codes

fsources += vode.f LinAlg.f math_d.f tranlib_d.f

ifdef USE_EGZ
  f90sources += egz_module.f90
else
  f90sources += eglib_module.f90
  fsources += EGini.f
  ifdef USE_EGM
     fsources += EGMlib.f
  else
     ifdef USE_EGF
        fsources += EGFlib.f
     else
        fsources += EGSlib.f
     endif
  endif  
endif

ifdef USE_WORK_SPACE_MODULE
  # vode, eglib & tranlib all need work space.  
  # Here, the work spaces are in F90 modules. 
  f90sources += vode_module.f90 tranlib_module.f90
  ifndef USE_EGZ
     f90sources += eglib_module.f90
  endif
endif
