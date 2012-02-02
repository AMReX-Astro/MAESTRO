! inlet_bc_module serves as a container to hold the inflow boundary 
! condition information.
!
! These quantities are initialized through a call to set_inlet_bcs(),
! which should be done on initialization and restart.

module inlet_bc_module

  use bl_types
  use bl_constants_module
  use bl_space
  use network

  implicit none

  real(dp_t), save :: INLET_VN
  real(dp_t), save :: INLET_VT
  real(dp_t), save :: INLET_RHO
  real(dp_t), save :: INLET_RHOH
  real(dp_t), save :: INLET_RHOC12
  real(dp_t), save :: INLET_RHOO16
  real(dp_t), save :: INLET_RHOMG24
  real(dp_t), save :: INLET_TEMP
  real(dp_t), save :: INLET_TRA

  logical, save :: inlet_bc_initialized = .false.

contains

  subroutine set_inlet_bcs()

    use eos_module

    ! local variables
    integer :: ndum,comp
    parameter (ndum=30)

    real(kind=dp_t) :: state1d(ndum)

    ! now reset inflow boundary conditions
    call asin1d('flame_4.e7_screen_left.out', -.00125d0, 0.d0, state1d, ndum, .false.)

       p_eos = state1d(18)
     den_eos = state1d(3)
    temp_eos = state1d(9)

    do comp=1,nspec
       if(spec_names(comp) .eq. "carbon-12") then
          xn_eos(comp) = state1d(21)
       else if(spec_names(comp) .eq. "magnesium-24") then
          xn_eos(comp) = state1d(22)
       else if(spec_names(comp) .eq. "oxygen-16") then
          xn_eos(comp) = state1d(23)
       else
          print*,"In initdata, spec_names(",comp,") invalid"
       endif
    enddo

    ! given P, T, and X, compute rho
    call eos(eos_input_tp, den_eos, temp_eos, &
             xn_eos, p_eos, h_eos, e_eos, & 
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             .false.)

    ! given rho, T, and X, compute h
    call eos(eos_input_rt, den_eos, temp_eos, &
             xn_eos, p_eos, h_eos, e_eos, & 
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             .false.)

    INLET_VN = 0.0d0
    INLET_VT = 0.0d0
    INLET_RHO = den_eos
    INLET_RHOH = den_eos*h_eos

    do comp=1,nspec
       if(spec_names(comp) .eq. "carbon-12") then
          INLET_RHOC12 = den_eos*xn_eos(comp)
       else if(spec_names(comp) .eq. "magnesium-24") then
          INLET_RHOMG24 = den_eos*xn_eos(comp)
       else if(spec_names(comp) .eq. "oxygen-16") then
          INLET_RHOO16 = den_eos*xn_eos(comp)
       endif
    enddo
    INLET_TEMP = temp_eos
    INLET_TRA = 0.0d0

    inlet_bc_initialized = .true.

  end subroutine set_inlet_bcs

end module inlet_bc_module
