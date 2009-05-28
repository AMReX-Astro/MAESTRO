! inlet_bc_module serves as a container to hold the inflow boundary 
! condition information.
!
! These quantities are initialized through a call to set_inlet_bcs(),
! which should be done on initialization and restart.

module inlet_bc_module

  use bl_types
  use bl_constants_module
  use bl_space
  use network

  implicit none

  real(dp_t), save    :: INLET_VEL         ! normal velocity through boundary
  real(dp_t), save    :: INLET_RHO
  real(dp_t), save    :: INLET_RHOH
  real(dp_t), save    :: INLET_TEMP
  real(dp_t), save    :: INLET_RHOX(nspec)
  real(dp_t), save    :: INLET_TRA

  logical, save :: inlet_bc_initialized = .false.

contains

  subroutine set_inlet_bcs(rho_in, T_in, xn_in, vel_in)

    use eos_module
    use network

    real(dp_t) :: rho_in, T_in
    real(dp_t) :: xn_in(nspec)
    real(dp_t) :: vel_in

    den_eos(1)  = rho_in
    temp_eos(1) = T_in
    xn_eos(1,:) = xn_in(:)
    
    call eos(eos_input_rt, den_eos, temp_eos, &
             npts, nspec, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             do_diag)

    INLET_RHO     = rho_in
    INLET_RHOH    = rho_in*h_eos(1)
    INLET_TEMP    = T_in
    INLET_RHOX(:) = rho_in*xn_in(:)
    INLET_VEL     = vel_in
    INLET_TRA     = ZERO

    inlet_bc_initialized = .true.

  end subroutine set_inlet_bcs

end module inlet_bc_module
