! a module to provide integer indices into the various storage arrays
! for accessing the different variables by name

module variables

  implicit none

  integer, save :: rho_comp, rhoh_comp, press_comp, derive_comp

contains

  subroutine init_variables(dm, nscal)

    integer, intent(in) :: dm, nscal

    rho_comp = dm + 1
    rhoh_comp = rho_comp + 1
    press_comp = dm + nscal + 1
    derive_comp = rho_comp + nscal

  end subroutine init_variables
end module variables
