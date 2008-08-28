module base_state_module

  use bl_types

  implicit none

  private

  public :: init_base_state

contains

  subroutine init_base_state(n,model_file,s0_init,p0_init,dx)

    use bc_module
    use setbc_module
    use multifab_module
    use define_bc_module
    use bl_constants_module
    use eos_module
    use probin_module, ONLY: prob_lo_y, prob_lo_z
    use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp, trac_comp
    use geometry, only: dr, spherical, r_start_coord, r_end_coord, numdisjointchunks, dm
    use inlet_bc_module

    integer,             intent(in   ) :: n
    character(len=256),  intent(in   ) :: model_file ! I'm not using this anymore
    real(kind=dp_t),     intent(inout) :: s0_init(0:,:)
    real(kind=dp_t),     intent(inout) :: p0_init(0:)
    real(kind=dp_t),     intent(in   ) :: dx(:)

    ! local
    integer :: ndum,r,comp,i
    real(dp_t) :: starting_rad

    parameter (ndum = 30)

    character(len=128) :: lamsolfile
    real(kind=dp_t) :: state1d(ndum),Pamb,temporary
    real(kind=dp_t) :: loloc,hiloc,flameloc

    if (dm .eq. 2) then
       starting_rad = prob_lo_y
    else if(dm .eq. 3) then
       starting_rad = prob_lo_z
    endif

    lamsolfile = 'flame_4.e7_screen_left.out'

    ! first set the inflow boundary condition
    call asin1d(lamsolfile, -.00125d0, 0.d0, state1d, ndum, .false.)

    Pamb = state1d(18)
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
       else
          print*,"In initdata, spec_names(",comp,") invalid"
       endif
    enddo

    ! given P, T, and X, compute rho
    call eos(eos_input_tp, den_eos, temp_eos, &
             npts, nspec, &
             xn_eos, &
             p_eos, h_eos, e_eos, & 
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             do_diag)

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

    INLET_VN = 0.0d0
    INLET_VT = 0.0d0
    INLET_RHO = den_eos(1)
    INLET_RHOH = den_eos(1)*h_eos(1)

    do comp=1,nspec
       if(spec_names(comp) .eq. "carbon-12") then
          INLET_RHOC12 = den_eos(1)*xn_eos(1,comp)
       else if(spec_names(comp) .eq. "magnesium-24") then
          INLET_RHOMG24 = den_eos(1)*xn_eos(1,comp)
       else if(spec_names(comp) .eq. "oxygen-16") then
          INLET_RHOO16 = den_eos(1)*xn_eos(1,comp)
       endif
    enddo
    INLET_TEMP = temp_eos(1)
    INLET_TRA = 0.0d0

    ! Now do the interior cells
    flameloc = ONE

    do i=1,numdisjointchunks(n)
       do r=r_start_coord(n,i),r_end_coord(n,i)

          loloc = starting_rad +  dble(r)     *dx(dm) - flameloc
          hiloc = starting_rad + (dble(r)+ONE)*dx(dm) - flameloc

          call asin1d(lamsolfile, loloc, hiloc, state1d, ndum, .false.)

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
               npts, nspec, &
               xn_eos, &
               p_eos, h_eos, e_eos, & 
               cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
               dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
               dpdX_eos, dhdX_eos, &
               gam1_eos, cs_eos, s_eos, &
               dsdt_eos, dsdr_eos, &
               do_diag)

          ! given rho, T, and X, compute h.
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

          s0_init(r,rho_comp) = den_eos(1)
          s0_init(r,rhoh_comp) = den_eos(1)*h_eos(1)

          do comp=1,nspec
             s0_init(r,spec_comp+comp-1) = den_eos(1)*xn_eos(1,comp)
          enddo
          s0_init(r,trac_comp) = 0.0d0
          s0_init(r,temp_comp) = temp_eos(1)
          p0_init(r) = pamb

       enddo
    end do

  end subroutine init_base_state

end module base_state_module
