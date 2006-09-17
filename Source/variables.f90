! a module to provide integer indices into the various storage arrays
! for accessing the different variables by name

module variables

  implicit none

  integer, save :: rho_comp, rhoh_comp, spec_comp, trac_comp, press_comp, derive_comp

contains

  subroutine init_variables(dm, nscal, nspec)

    integer, intent(in) :: dm, nscal, nspec

    rho_comp    = 1
    rhoh_comp   = 2
    spec_comp   = rhoh_comp + 1
    trac_comp = spec_comp + nspec
    press_comp  = dm + nscal + 1
    derive_comp = dm + nscal + 1

  end subroutine init_variables

end module variables
