! this defines the "base state" for a small-scale laminar flame problem.
! basically, we set everything but the pressure to 0.  This will mimic
! the way that the original smallscale low Mach number code dealt with 
! flames.
!
! This routine will also define the inflow boundary condition state

module base_state_module

  use bl_types

  implicit none

  private

  public :: init_base_state

contains

  subroutine init_base_state(n,model_file,s0_init,p0_init,dx)

    use bl_constants_module
    use network
    use eos_module
    use probin_module, ONLY: dens_fuel, temp_fuel, xc12_fuel, vel_fuel
    use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp, trac_comp, ntrac
    use geometry, only: dr, spherical, nr, dm
    use inlet_bc_module, only: set_inlet_bcs

    integer,             intent(in   ) :: n
    character(len=256),  intent(in   ) :: model_file ! Not used
    real(kind=dp_t),     intent(inout) :: s0_init(0:,:)
    real(kind=dp_t),     intent(inout) :: p0_init(0:)
    real(kind=dp_t),     intent(in   ) :: dx(:)

    ! local
    integer :: ic12, io16
    real(kind=dp_t) :: p_ambient
    

    ! figure out the indices for different species
    ic12  = network_species_index("carbon-12")
    io16  = network_species_index("oxygen-16")
    

    ! determine the ambient pressure from the inflow conditions
    ! (which define the fuel state)
    den_eos(1)  = dens_fuel
    temp_eos(1) = temp_fuel

    xn_eos(1,:) = ZERO
    xn_eos(1,ic12) = xc12_fuel
    xn_eos(1,io16) = 1.d0 - xc12_fuel
    
    ! given rho, T, and X, compute h
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

    p_ambient = p_eos(1)


    ! set the base state quantities (except pressure) to be 0, and
    ! set the base state pressure to be the ambient pressure.
    s0_init(:,rho_comp)  = ZERO
    s0_init(:,rhoh_comp) = ZERO
    s0_init(:,spec_comp:spec_comp+nspec-1) = ZERO
    s0_init(:,temp_comp) = ZERO
    p0_init(:) = p_ambient
    if (ntrac > 0) then
       s0_init(:,trac_comp:trac_comp+ntrac-1) = ZERO
    endif


    ! define the inflow boundary condition parameters
    call set_inlet_bcs()

  end subroutine init_base_state

end module base_state_module
