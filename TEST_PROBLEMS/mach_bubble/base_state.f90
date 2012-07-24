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
    use eos_module, only: eos_input_rp, eos, gamma_const
    use eos_type_module
    use network, only: spec_names, network_species_index
    use probin_module, only: base_cutoff_density, anelastic_cutoff, &
                             buoyancy_cutoff_factor, grav_const, dens_base, pres_base, do_isentropic

    use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp, trac_comp, ntrac
    use geometry, only: dr, spherical, nr
    
    integer           , intent(in   ) :: n
    character(len=256), intent(in   ) :: model_file
    real(kind=dp_t)   , intent(inout) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(inout) :: p0_init(0:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    ! local
    integer         :: r,j
    real(kind=dp_t) :: H,z,const
    real(kind=dp_t) :: dens_zone, temp_zone, pres_zone
    real(kind=dp_t) :: xn_zone(nspec)

    type (eos_t) :: eos_state

    type(bl_prof_timer), save :: bpt

    call build(bpt, "init_base_state")

    if (spherical .eq. 1) then
       call bl_error("ERROR: base_state.f90 is not valid for spherical")
    endif

    ! for isentropic
    ! an isentropic gas satisfies p ~ rho^gamma
    const = pres_base/dens_base**gamma_const

    ! use the EOS to make the state consistent
    dens_zone  = dens_base
    pres_zone    = pres_base

    ! only initialize the first species
    xn_zone(:) = ZERO
    xn_zone(1) = 1.d0

    p0_init(0) = pres_zone

    ! compute the pressure scale height (for an isothermal, ideal-gas
    ! atmosphere)
    H = pres_base / dens_base / abs(grav_const)

    ! set an initial guess for the temperature -- this will be reset
    ! by the EOS
    temp_zone = 1000.d0

    do j = 0, nr(n)-1

       if (do_isentropic) then

          z = dble(j) * dr(1)

          ! we can integrate HSE with p = K rho^gamma analytically
          dens_zone = dens_base*(grav_const*dens_base*(gamma_const - 1.0)*z/ &
               (gamma_const*pres_base) + 1.d0)**(1.d0/(gamma_const - 1.d0))

       else
      
          z = (dble(j)+HALF) * dr(1)

          ! the density of an isothermal gamma-law atm is exponential
          dens_zone = dens_base * exp(-z/H)

       end if

       s0_init(j, rho_comp) = dens_zone

       ! compute the pressure by discretizing HSE
       if (j.gt.0) then
          p0_init(j) = p0_init(j-1) - &
               dr(1) * HALF * (s0_init(j,rho_comp) + s0_init(j-1,rho_comp)) * &
               abs(grav_const)
          
          pres_zone = p0_init(j)
       end if
       
       ! use the EOS to make the state consistent
       eos_state%rho   = dens_zone
       eos_state%p     = pres_zone
       eos_state%T     = temp_zone
       eos_state%xn(:) = xn_zone(:)

       ! (rho,p) --> T, h
       call eos(eos_input_rp, eos_state, .false.)

       s0_init(j, rho_comp) = dens_zone
       s0_init(j,rhoh_comp) = dens_zone * eos_state%h

       s0_init(j,spec_comp:spec_comp-1+nspec) = ZERO
       s0_init(j,spec_comp) = dens_zone

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
