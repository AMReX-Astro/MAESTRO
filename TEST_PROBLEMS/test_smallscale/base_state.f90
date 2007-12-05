module base_state_module

  use bl_types
  use bl_constants_module
  use bc_module
  use setbc_module
  use define_bc_module
  use multifab_module
  use eos_module
  use variables
  use network
  use geometry

  implicit none

  private
  public :: init_base_state

contains

  subroutine init_base_state (model_file,n_base,s0,p0,gam1,dx,prob_lo,prob_hi)

    character (len=256), intent(in   ) :: model_file ! I'm not using this anymore
    integer        ,     intent(in   ) :: n_base
    real(kind=dp_t),     intent(inout) ::    s0(0:,:)
    real(kind=dp_t),     intent(inout) ::    p0(0:)
    real(kind=dp_t),     intent(inout) ::  gam1(0:)
    real(kind=dp_t),     intent(in   ) :: prob_lo(:)
    real(kind=dp_t),     intent(in   ) :: prob_hi(:)
    real(kind=dp_t),     intent(in   ) :: dx(:)

    ! local
    integer ndum,i,dm,nspec,comp

    parameter (ndum = 30)
    parameter (nspec = 3)

    character(len=128) :: lamsolfile
    real(kind=dp_t) :: state1d(ndum),Pamb,temporary
    real(kind=dp_t) :: loloc,hiloc,flameloc,qreact
    
    dm = size(dx)

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
    if(use_big_h) then
       qreact = 0.0d0
       do comp=1,nspec
          qreact = qreact + ebin(comp)*xn_eos(1,comp)
       enddo
       INLET_RHOH = den_eos(1)*(h_eos(1) + qreact)
    else
       INLET_RHOH = den_eos(1)*h_eos(1)
    endif
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

    do i=0,n_base-1

       loloc = dble(i)*dx(dm) - flameloc
       hiloc = (dble(i) + ONE)*dx(dm) - flameloc

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

       s0(i,rho_comp) = den_eos(1)
       if(use_big_h) then
          qreact = ZERO
          do comp=1,nspec
             qreact = qreact + ebin(comp)*xn_eos(1,comp)
          enddo
          temporary = h_eos(1) + qreact
          s0(i,rhoh_comp) = den_eos(1)*temporary
       else
          s0(i,rhoh_comp) = den_eos(1)*h_eos(1)
       endif
       do comp=1,nspec
          s0(i,spec_comp+comp-1) = den_eos(1)*xn_eos(1,comp)
       enddo
       s0(i,trac_comp) = 0.0d0
       s0(i,temp_comp) = temp_eos(1)
       p0(i) = pamb
       gam1(i) = gam1_eos(1)

    enddo

  end subroutine init_base_state

end module base_state_module
