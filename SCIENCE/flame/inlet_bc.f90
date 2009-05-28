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

  subroutine set_inlet_bcs()

    ! initialize the inflow boundary condition variables
    use eos_module
    use network
    use probin_module, ONLY: dens_fuel, temp_fuel, xc12_fuel, vel_fuel

    integer :: ic12, io16

    ! figure out the indices for different species
    ic12  = network_species_index("carbon-12")
    io16  = network_species_index("oxygen-16")

    den_eos(1)  = dens_fuel
    temp_eos(1) = temp_fuel

    xn_eos(1,:) = ZERO
    xn_eos(1,ic12) = xc12_fuel
    xn_eos(1,io16) = 1.d0 - xc12_fuel

    
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

    INLET_RHO     = dens_fuel
    INLET_RHOH    = dens_fuel*h_eos(1)
    INLET_TEMP    = temp_fuel
    INLET_RHOX(:) = dens_fuel*xn_eos(1,:)
    INLET_VEL     = vel_fuel
    INLET_TRA     = ZERO

    inlet_bc_initialized = .true.

  end subroutine set_inlet_bcs

end module inlet_bc_module
