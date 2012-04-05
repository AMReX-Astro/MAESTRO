! inlet_bc_module is a simple container module that holds the parameters
! that are used by physbc to implement the inlet boundary conditions.
! As these are problem-specific, any problem needing inlet boundary 
! conditions should create its own version of this module, using this
! outline.

module inlet_bc_module

  use bl_types
  use bl_constants_module
  use bl_space
  use network
  use eos_module

  implicit none

  ! parameters that would be used by physbc in the EXT_DIR sections
  ! would be stored here with the 'save' attribute
  real(dp_t), save :: INLET_RHO
  real(dp_t), save :: INLET_RHOH
  real(dp_t), save :: INLET_TEMP
  real(dp_t), save :: INLET_CS

  logical, save :: inlet_bc_initialized = .false.

contains

  ! here we would initialize the parameters that are module variables.
  ! this routine is called when the base state is defined initially,
  ! and upon restart, just after the base state is read in.
  subroutine set_inlet_bcs()

    temp_eos = 10.d0
    den_eos = 1.d-3
    p_eos = 1.d6
    xn_eos(:) = 1.d0

    call eos(eos_input_rp, den_eos, temp_eos, &
             xn_eos, &
             p_eos, h_eos, e_eos, & 
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             .false.)

    INLET_CS   = cs_eos

    temp_eos = 10.d0
    den_eos = 5.d-4
    p_eos = 1.d6
    xn_eos(:) = 1.d0

    call eos(eos_input_rp, den_eos, temp_eos, &
             xn_eos, &
             p_eos, h_eos, e_eos, & 
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             .false.)

    INLET_RHO  = den_eos
    INLET_RHOH = den_eos*h_eos
    INLET_TEMP = temp_eos

    inlet_bc_initialized = .true.

  end subroutine set_inlet_bcs

end module inlet_bc_module
