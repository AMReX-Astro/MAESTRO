! initialize an isentropic layer below an isothermal layer.  The
! transition happens at a height r_mid.

module base_state_module

  use bl_types
  use network, only: nspec
  
  implicit none

  private

  public :: init_base_state

contains

  subroutine init_base_state(n,model_file,s0_init,p0_init,dx)

    use bl_prof_module
    use parallel
    use bl_error_module
    use bl_constants_module
    use eos_module, only: eos, eos_input_rp, gamma_const
    use eos_type_module
    use network, only: spec_names, network_species_index
    use probin_module, only: base_cutoff_density, anelastic_cutoff, &
                             buoyancy_cutoff_factor, grav_const, dens_base, pres_base, &
                             prob_lo, r_mid
    use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp, trac_comp, ntrac
    use geometry, only: dr, spherical, nr
    
    integer           , intent(in   ) :: n
    character(len=256), intent(in   ) :: model_file
    real(kind=dp_t)   , intent(inout) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(inout) :: p0_init(0:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    ! local
    integer         :: j
    real(kind=dp_t) :: H,r,r_base
    real(kind=dp_t) :: dens_mid, dens_zone, pres_zone

    logical :: do_isentropic, first_isothermal

    type(bl_prof_timer), save :: bpt

    type(eos_t) :: eos_state

    call build(bpt, "init_base_state")

    if (spherical .eq. 1) then
       call bl_error("ERROR: base_state.f90 is not valid for spherical")
    endif

    r_base = prob_lo(size(dx))

    

    ! set the pressure at the base -- we integrate HSE from this
    p0_init(0) = pres_base


    ! set an initial guess for the temperature -- this will be reset
    ! by the EOS
    eos_state%T = 1000.d0

    first_isothermal = .true.

    do j = 0, nr(n)-1

       r = dble(j) * dr(1) + r_base

       if (r < r_mid) then
          do_isentropic = .true.
       else
          do_isentropic = .false.
       endif


       if (do_isentropic) then

          ! we can integrate HSE with p = K rho^gamma analytically
          dens_zone = dens_base*(grav_const*dens_base*(gamma_const - 1.0)*(r-r_base)/ &
               (gamma_const*pres_base) + 1.d0)**(1.d0/(gamma_const - 1.d0))

       else

          if (first_isothermal) then
             dens_mid = s0_init(j-1, rho_comp)

             ! compute the pressure scale height (for an isothermal,
             ! ideal-gas atmosphere)
             H = p0_init(j-1) / dens_mid / abs(grav_const)

             first_isothermal = .false.
          endif
      
          ! the density of an isothermal gamma-law atm is exponential
          dens_zone = dens_mid * exp(-(r-r_mid)/H)

       end if

       s0_init(j, rho_comp) = dens_zone

       ! compute the pressure by discretizing HSE
       if (j > 0) then
          p0_init(j) = p0_init(j-1) - &
               dr(1) * HALF * (s0_init(j,rho_comp) + s0_init(j-1,rho_comp)) * &
               abs(grav_const)          
       end if
       
       ! use the EOS to make the state consistent

       ! only initialize the first species
       eos_state%xn(:) = ZERO
       eos_state%xn(1) = ONE

       eos_state%rho = s0_init(j, rho_comp)
       eos_state%p = p0_init(j)

       ! (rho,p) --> T, h
       call eos(eos_input_rp, eos_state, .false.)

       s0_init(j,rhoh_comp) = s0_init(j, rho_comp) * eos_state%h

       s0_init(j,spec_comp:spec_comp-1+nspec) = ZERO
       s0_init(j,spec_comp) = s0_init(j, rho_comp)

       s0_init(j,temp_comp) = eos_state%T
       s0_init(j,trac_comp) = ZERO


    end do

    ! HSE check
    !do j = 1, nr(n)-1
    !   print *,'DPDR RHO*G ',j, (p0_init(j)-p0_init(j-1)) / dr(1), &
    !            .5d0*(s0_init(j,rho_comp)+s0_init(j-1,rho_comp))*grav_const
    !end do

    call destroy(bpt)

  end subroutine init_base_state

end module base_state_module
