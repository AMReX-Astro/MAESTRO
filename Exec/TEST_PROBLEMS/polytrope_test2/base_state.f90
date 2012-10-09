module base_state_module

  use bl_types
  use network

  implicit none

  public :: init_base_state

contains

  subroutine init_base_state(n,model_file,s0_init,p0_init,dx)

    use bl_prof_module
    use eos_module
    use probin_module, only: rho0, grav_const, H, &
         prob_lo, prob_hi, nn, K, dlow, dhigh, Tlow, Thigh
    use bl_constants_module
    use geometry, only: dr, nr
    use inlet_bc_module, only: set_inlet_bcs
    use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp, &
         trac_comp, ntrac

    integer           , intent(in   ) :: n
    character(len=256), intent(in   ) :: model_file
    real(kind=dp_t)   , intent(inout) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(inout) :: p0_init(0:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    type(bl_prof_timer), save :: bpt

    real(kind=dp_t) :: rloc, rmax
    integer :: r

    real(kind=dp_t) :: p0_lower,p0_upper
    real(kind=dp_t) :: xn_lower(nspec), xn_upper(nspec)
    integer :: ia, ib

    real(kind=dp_t) :: d_ambient, p_ambient, xn_ambient(nspec)

    real(kind=dp_t) :: max_hse_error, dpdr, rhog

    real(kind=dp_t), parameter :: SMALL = 1.d-12

    call build(bpt, "init_base_state")

    rmax = (dble(nr(n)-1)+HALF)*dr(n)

    ! fill the base state arrays
    do r=-1,nr(n)

       !height above the bottom of the domain
       if (r .eq. -1) then
          rloc = 0.0d0
       else if (r .eq. nr(n)) then
          rloc = dble(r)*dr(n)
       else
          rloc = (dble(r)+HALF)*dr(n)
       end if

       d_ambient = (rho0**(1.0d0/nn) + grav_const / (nn+1.0d0) / K * (rloc-rmax))**(nn)
       p_ambient = K * (rho0**(1.0d0/nn) + grav_const / (nn+1.0d0) / K * (rloc-rmax))**(nn+1.0d0)
!       d_ambient = rho0
!       p_ambient = 10.d0 + rho0 * grav_const * rloc
       xn_ambient(:) = SMALL
       xn_ambient(1) = ONE - (nspec - 1)*SMALL

       den_eos   = d_ambient
       p_eos     = p_ambient
       xn_eos(:) = xn_ambient(:)

       ! (rho,p) --> T, h
       call eos(eos_input_rp, den_eos, temp_eos, xn_eos, p_eos, h_eos,&
            e_eos, cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
            dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
            dpdX_eos, dhdX_eos, &
            gam1_eos, cs_eos, s_eos, &
            dsdt_eos, dsdr_eos, &
            .false., &
            pt_index_eos)

       if (r .eq. -1) then
          dlow = den_eos
          Tlow = temp_eos
       else if (r .eq. nr(n)) then
          dhigh = den_eos
          Thigh = temp_eos
       else
          s0_init(r, rho_comp) = den_eos
          s0_init(r,rhoh_comp) = den_eos * h_eos
          s0_init(r,spec_comp:spec_comp+nspec-1) = den_eos * xn_eos(:)
          p0_init(r) = p_eos
          s0_init(r,temp_comp) = temp_eos
       end if

    end do

    ! find size of HSE errors
    max_hse_error = -1.d30
    do r=1,nr(n)-1

       rloc = prob_lo(size(dx)) + (dble(r) + HALF)*dr(n)

       dpdr = (p0_init(r)-p0_init(r-1))/dr(n)
       rhog = HALF*(s0_init(r,rho_comp) + s0_init(r-1,rho_comp))*grav_const

       max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(dpdr))

    enddo

    call set_inlet_bcs()

    call destroy(bpt)

  end subroutine init_base_state

end module base_state_module
