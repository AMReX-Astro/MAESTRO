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
    use bl_error_module
    use network
    use eos_module
    use probin_module, ONLY: dens_fuel, temp_fuel, xc12_fuel, vel_fuel, frac, &
         anelastic_cutoff, base_cutoff_density, prob_lo, prob_hi, temp_ash, temp_fuel
    use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp, trac_comp, ntrac
    use geometry, only: dr, spherical, nr, dm
    use inlet_bc_module, only: set_initial_inlet_bcs

    integer,             intent(in   ) :: n
    character(len=256),  intent(in   ) :: model_file ! Not used
    real(kind=dp_t),     intent(inout) :: s0_init(0:,:)
    real(kind=dp_t),     intent(inout) :: p0_init(0:)
    real(kind=dp_t),     intent(in   ) :: dx(:)

    ! local
    real(kind=dp_t) :: rlen, rloc
    integer :: r, ic12, io16, img24
    real(kind=dp_t) :: p_ambient, dens_ash, rhoh_fuel, rhoh_ash
    real(kind=dp_t) :: xn_fuel(nspec), xn_ash(nspec)
    

    ! figure out the indices for different species
    ic12  = network_species_index("carbon-12")
    io16  = network_species_index("oxygen-16")
    img24  = network_species_index("magnesium-24")    

    ! length of the domain
    rlen = (prob_hi(dm) - prob_lo(dm))

    ! figure out the thermodynamics of the fuel and ash state

    ! fuel
    xn_fuel(:)    = ZERO
    xn_fuel(ic12) = xc12_fuel
    xn_fuel(io16) = 1.d0 - xc12_fuel

    den_eos(1)  = dens_fuel
    temp_eos(1) = temp_fuel
    xn_eos(1,:) = xn_fuel(:)

    call eos(eos_input_rt, den_eos, temp_eos, &
             npts, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             .false.)

    ! note: p_ambient should be = p0_init
    p_ambient = p_eos(1)
    rhoh_fuel = dens_fuel*h_eos(1)

    ! ash
    xn_ash(:)     = ZERO
    xn_ash(ic12)  = ZERO
    xn_ash(io16)  = 1.d0 - xc12_fuel    
    xn_ash(img24) = xc12_fuel

    den_eos(1)  = dens_fuel    ! initial guess
    temp_eos(1) = temp_ash
    xn_eos(1,:) = xn_ash(:)
    p_eos(1) = p_ambient

    call eos(eos_input_tp, den_eos, temp_eos, &
             npts, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             .false.)

    dens_ash = den_eos(1)
    rhoh_ash = dens_ash*h_eos(1)

    do r=0,nr(n)-1

       rloc = prob_lo(dm) + (dble(r)+0.5d0)*dr(n)

       ! the flame propagates in the -y direction.  If we are more than
       ! fuel/ash division is frac through the domain
       if (rloc < prob_lo(dm) + frac*rlen) then
          
          ! fuel
          s0_init(r,rho_comp)  = dens_fuel
          s0_init(r,rhoh_comp) = rhoh_fuel
          s0_init(r,temp_comp) = temp_fuel
          s0_init(r,spec_comp:spec_comp+nspec-1) = dens_fuel*xn_fuel(:)
          s0_init(r,trac_comp:trac_comp+ntrac-1) = ZERO
          
       else
          
          ! ash
          s0_init(r,rho_comp)  = dens_ash
          s0_init(r,rhoh_comp) = rhoh_ash
          s0_init(r,temp_comp) = temp_ash
          s0_init(r,spec_comp:spec_comp+nspec-1) = dens_ash*xn_ash(:)
          s0_init(r,trac_comp:trac_comp+ntrac-1) = ZERO
          
       endif

    enddo

    ! set the base state pressure to be the ambient pressure.
    p0_init(:) = p_ambient

    ! sanity check
    if (dens_fuel < base_cutoff_density .or. &
        dens_fuel < anelastic_cutoff) then
       call bl_error('ERROR: fuel density < (base_cutoff_density or anelastic_cutoff)')
    endif

    call set_initial_inlet_bcs()

  end subroutine init_base_state

end module base_state_module
