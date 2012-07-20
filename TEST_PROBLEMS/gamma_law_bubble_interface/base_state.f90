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
    use eos_module, only: eos, eos_input_rp, eos_input_ps, eos_input_tp, gamma_const
    use eos_type_module
    use network, only: spec_names, network_species_index
    use probin_module, only: base_cutoff_density, anelastic_cutoff, &
                             buoyancy_cutoff_factor, grav_const, dens_base, pres_base, &
                             prob_lo, r_mid, s_jump, do_isentropic_above
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
    real(kind=dp_t) :: dens_mid, pres_mid, s_mid, dens_zone, pres_zone, dens_zone_old
    real(kind=dp_t) :: drho, s_want, T_want

    integer :: iter
    real(kind=dp_t) :: tol = 1.e-8_dp_t

    logical :: do_isentropic, first_isothermal, mid_isentropic, reset_isentropic

    type(bl_prof_timer), save :: bpt

    type(eos_t) :: eos_state

    call build(bpt, "init_base_state")

    if (spherical .eq. 1) then
       call bl_error("ERROR: base_state.f90 is not valid for spherical")
    endif

    r_base = prob_lo(size(dx))

    

    ! set the pressure and density at the base -- we integrate HSE from this
    p0_init(0) = pres_base
    s0_init(0, rho_comp) = dens_base


    ! only initialize the first species
    eos_state%xn(:) = ZERO
    eos_state%xn(1) = ONE


    ! set an initial guess for the temperature -- this will be reset
    ! by the EOS
    eos_state%T = 1000.d0


    ! get initial entropy
    eos_state%rho = dens_base
    eos_state%p = pres_base

    call eos(eos_input_rp, eos_state, .false.)


    ! store the thermodynamics for the base 
    s0_init(j,rhoh_comp) = s0_init(j, rho_comp) * eos_state%h

    s0_init(j,spec_comp:spec_comp-1+nspec) = ZERO
    s0_init(j,spec_comp) = s0_init(j, rho_comp)

    s0_init(j,temp_comp) = eos_state%T
    s0_init(j,trac_comp) = ZERO


    ! set the initial entropy to constrain to
    s_want = eos_state%s


    reset_isentropic = .true.
    first_isothermal = .true.

    ! now do the integration
    do j = 1, nr(n)-1

       r = dble(j) * dr(1) + r_base

       ! initial guess for rho
       dens_zone = s0_init(j-1, rho_comp)

       drho = 1.d33
       dens_zone_old = 1.d33

       iter = 1


       ! loop until density converges
       do while (abs(drho) > tol*dens_zone)

          ! what pressure does HSE want?
          pres_zone = p0_init(j-1) - &
               dr(1) * HALF * (dens_zone + s0_init(j-1,rho_comp)) * &
               abs(grav_const)                    


          ! what are we constraining to?
          if (r < r_mid) then

             ! constraining to entropy
             eos_state%p = pres_zone
             eos_state%s = s_want

             call eos(eos_input_ps, eos_state, .false.)
             
             dens_zone = eos_state%rho

          else
             if (do_isentropic_above) then

                ! constraining to 'jumped' entropy

                ! first find the new entropy
                if (reset_isentropic) then

                   eos_state%rho = s0_init(j-1,rho_comp)
                   eos_state%p = p0_init(j-1)

                   ! (rho,p) --> s
                   call eos(eos_input_rp, eos_state, .false.)
             
                   ! jump the entropy
                   s_want = s_jump*eos_state%s

                   reset_isentropic = .false.
                endif

                ! now find the density corresponding to these conditions
                eos_state%p = pres_zone
                eos_state%s = s_want

                call eos(eos_input_ps, eos_state, .false.)

                dens_zone = eos_state%rho

             else
                
                ! constraining to isothermal
                if (first_isothermal) then

                   eos_state%rho = s0_init(j-1,rho_comp)
                   eos_state%p = p0_init(j-1)

                   ! (rho,p) --> T
                   call eos(eos_input_rp, eos_state, .false.)
             
                   ! jump the entropy
                   eos_state%s = s_jump*eos_state%s

                   ! now find the temperature (p,s) -> T
                   call eos(eos_input_ps, eos_state, .false.)

                   T_want = eos_state%T

                   first_isothermal = .false.
                endif

                ! now find the density corresponding to these conditions
                eos_state%p = pres_zone
                eos_state%T = T_want

                call eos(eos_input_tp, eos_state, .false.)

                dens_zone = eos_state%rho

             endif
          endif

          drho = dens_zone - dens_zone_old
          dens_zone_old = dens_zone

          iter = iter + 1
          if (iter > 100) call bl_error("too many iterations in base_state")

       enddo


       s0_init(j, rho_comp) = dens_zone
       p0_init(j) = pres_zone
       
       ! use the EOS to make the state consistent
       eos_state%rho = s0_init(j, rho_comp)
       eos_state%p = p0_init(j)

       ! (rho,p) --> T, h
       call eos(eos_input_rp, eos_state, .false.)

       s0_init(j,rhoh_comp) = s0_init(j, rho_comp) * eos_state%h

       s0_init(j,spec_comp:spec_comp-1+nspec) = ZERO
       s0_init(j,spec_comp) = s0_init(j, rho_comp)

       s0_init(j,temp_comp) = eos_state%T
       s0_init(j,trac_comp) = ZERO

       print *, j, s0_init(j, rho_comp), p0_init(j), s0_init(j,temp_comp), s0_init(j,spec_comp:spec_comp-1+nspec)

    end do

    ! HSE check
    !do j = 1, nr(n)-1
    !   print *,'DPDR RHO*G ',j, (p0_init(j)-p0_init(j-1)) / dr(1), &
    !            .5d0*(s0_init(j,rho_comp)+s0_init(j-1,rho_comp))*grav_const
    !end do

    call destroy(bpt)

  end subroutine init_base_state

end module base_state_module
