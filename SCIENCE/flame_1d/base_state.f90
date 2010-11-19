! This defines the "base state" for a small-scale laminar flame problem.
!
! We start by initializing the base state to represent the fuel and ash
! separated by a smoothed interface.  In initdata, we use this to map into
! the full state.  The base state is then overwritten, resetting everything
! but the pressure to 0.
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
    use probin_module, ONLY: dens_fuel, temp_fuel, xc12_fuel, vel_fuel, &
         interface_pos_frac, smooth_len_frac, &
         anelastic_cutoff, base_cutoff_density, prob_lo, prob_hi, &
         temp_ash, temp_fuel
    use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp, trac_comp, ntrac
    use geometry, only: dr, spherical, nr, dr_fine
    use inlet_bc_module, only: set_inlet_bcs

    integer,             intent(in   ) :: n
    character(len=256),  intent(in   ) :: model_file ! Not used
    real(kind=dp_t),     intent(inout) :: s0_init(0:,:)
    real(kind=dp_t),     intent(inout) :: p0_init(0:)
    real(kind=dp_t),     intent(in   ) :: dx(:)

    ! local
    real(kind=dp_t) :: rlen, rloc
    integer :: r, ic12, io16, img24
    real(kind=dp_t) :: p_ambient, dens_ash, rhoh_fuel, rhoh_ash
    real(kind=dp_t) :: xn_fuel(nspec), xn_ash(nspec), xn_smooth(nspec)

    ! figure out the indices for different species
    ic12  = network_species_index("carbon-12")
    io16  = network_species_index("oxygen-16")
    img24  = network_species_index("magnesium-24")    

    ! length of the domain
    rlen = (prob_hi(size(dx)) - prob_lo(size(dx)))

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

    ! initialize the fuel and ash, but put in a smooth temperature
    ! profile -- this means that we'll need to go back through an
    ! adjust the thermodynamics
    do r=0,nr(n)-1

       rloc = prob_lo(size(dx)) + (dble(r)+0.5d0)*dr(n)

       ! the flame propagates in the -y direction.  The fuel/ash division 
       ! is interface_pos_frac through the domain
       if (rloc < prob_lo(size(dx)) + interface_pos_frac*rlen) then
          
          ! fuel
          s0_init(r,rho_comp)  = dens_fuel
          s0_init(r,rhoh_comp) = rhoh_fuel
          !s0_init(r,temp_comp) = temp_fuel
          s0_init(r,spec_comp:spec_comp+nspec-1) = dens_fuel*xn_fuel(:)
          s0_init(r,trac_comp:trac_comp+ntrac-1) = ZERO
          
       else
          
          ! ash
          s0_init(r,rho_comp)  = dens_ash
          s0_init(r,rhoh_comp) = rhoh_ash
          !s0_init(r,temp_comp) = temp_ash
          s0_init(r,spec_comp:spec_comp+nspec-1) = dens_ash*xn_ash(:)
          s0_init(r,trac_comp:trac_comp+ntrac-1) = ZERO
          
       endif

       ! give the temperature a smooth profile
       s0_init(r,temp_comp) = temp_fuel + (temp_ash - temp_fuel) * &
            HALF * (ONE + &
            tanh( (rloc - (prob_lo(size(dx)) + interface_pos_frac*rlen)) / &
            (smooth_len_frac*rlen) ) )

       ! give the carbon mass fraction a smooth profile too
       xn_smooth(ic12) = xn_fuel(ic12) + (xn_ash(ic12) - xn_fuel(ic12)) * &
            HALF * (ONE + &
            tanh( (rloc - (prob_lo(size(dx)) + interface_pos_frac*rlen)) / &
            (smooth_len_frac*rlen) ) )

       xn_smooth(io16) = xn_fuel(io16)
       xn_smooth(img24) = 1.d0 - xn_smooth(ic12) - xn_smooth(io16)

       ! get the new density and enthalpy 
       den_eos(1)  = s0_init(r,rho_comp)
       temp_eos(1) = s0_init(r,temp_comp)
       xn_eos(1,:) = xn_smooth(:)
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

       s0_init(r,rho_comp)  = den_eos(1)
       s0_init(r,spec_comp:spec_comp+nspec-1) = den_eos(1)*xn_smooth(:)
       s0_init(r,rhoh_comp) = den_eos(1)*h_eos(1)

    enddo

    ! set the base state pressure to be the ambient pressure.
    p0_init(:) = p_ambient

    ! sanity check
    if (dens_fuel < base_cutoff_density .or. &
        dens_fuel < anelastic_cutoff) then
       call bl_error('ERROR: fuel density < (base_cutoff_density or anelastic_cutoff)')
    endif

    ! define the inflow boundary condition parameters
    call set_inlet_bcs()

  end subroutine init_base_state

end module base_state_module
