module base_state_module

  use bl_types

  implicit none

  private

  public :: init_base_state

contains

  subroutine init_base_state(n,model_file,s0_init,p0_init,dx)

    use parallel
    use bl_error_module
    use bl_constants_module
    use eos_module
    use network, only: spec_names, nspec
    use probin_module, ONLY: prob_lo, dm_in
    use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp, trac_comp
    use geometry, only: dr, spherical, nr
    use inlet_bc_module

    integer,             intent(in   ) :: n
    character(len=256),  intent(in   ) :: model_file ! Not used
    real(kind=dp_t),     intent(inout) :: s0_init(0:,:)
    real(kind=dp_t),     intent(inout) :: p0_init(0:)
    real(kind=dp_t),     intent(in   ) :: dx(:)

    ! local
    integer :: ndum,r,comp,dm
    parameter (ndum = 30)

    real(kind=dp_t) :: state1d(ndum),Pamb,starting_rad
    real(kind=dp_t) :: loloc,hiloc,flameloc

    dm = dm_in

    ! set bottom of domain
    starting_rad = prob_lo(dm)

    ! set flame location
    flameloc = ONE

    ! first set the inflow boundary condition
    call set_inlet_bcs()

    ! load in the background pressure
    call asin1d('flame_4.e7_screen_left.out', -.00125d0, 0.d0, state1d, ndum, .false.)
    Pamb = state1d(18)

    do r=0,nr(n)-1

       loloc = starting_rad +  dble(r)     *dx(dm) - flameloc
       hiloc = starting_rad + (dble(r)+ONE)*dx(dm) - flameloc

       call asin1d('flame_4.e7_screen_left.out', loloc, hiloc, state1d, ndum, .false.)

       p_eos(1) = Pamb
       den_eos(1) = state1d(3)
       temp_eos(1) = state1d(9)
       do comp=1,nspec
          if(spec_names(comp) .eq. "carbon-12") then
             xn_eos(1,comp) = state1d(21)
          else if(spec_names(comp) .eq. "magnesium-24") then
             xn_eos(1,comp) = state1d(22)
          else if(spec_names(comp) .eq. "oxygen-16") then
             xn_eos(1,comp) = state1d(23)
          endif
       enddo

       ! given P, T, and X, compute rho
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

       ! given rho, T, and X, compute h.
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

       s0_init(r,rho_comp) = den_eos(1)
       s0_init(r,rhoh_comp) = den_eos(1)*h_eos(1)

       do comp=1,nspec
          s0_init(r,spec_comp+comp-1) = den_eos(1)*xn_eos(1,comp)
       enddo
       s0_init(r,trac_comp) = 0.0d0
       s0_init(r,temp_comp) = temp_eos(1)
       p0_init(r) = pamb

    enddo

  end subroutine init_base_state

end module base_state_module
