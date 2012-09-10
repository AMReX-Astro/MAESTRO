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

    use eos_module, only: eos_input_rt, eos_input_rp, eos
    use eos_type_module

    ! local variables
    integer :: ndum,comp
    parameter (ndum=30)

    real(kind=dp_t) :: state1d(ndum)

    type (eos_t) :: eos_state

    ! now reset inflow boundary conditions
    call asin1d('flame_4.e7_screen_left.out', -.00125d0, 0.d0, state1d, ndum, .false.)

    eos_state%p   = state1d(18)
    eos_state%rho = state1d(3)
    eos_state%T   = state1d(9)

    eos_state%xn(:) = ZERO

    do comp=1,nspec
       if(spec_names(comp) .eq. "carbon-12") then
          eos_state%xn(comp) = state1d(21)

       else if(spec_names(comp) .eq. "magnesium-24") then
          eos_state%xn(comp) = state1d(22)

       else if(spec_names(comp) .eq. "oxygen-16") then
          eos_state%xn(comp) = state1d(23)

       else
          print*,"In initdata, spec_names(",comp,") invalid"
       endif
    enddo

    ! given P, T, and X, compute rho
    call eos(eos_input_tp, eos_state, .false.)

    ! given rho, T, and X, compute h
    call eos(eos_input_rt, eos_state, .false.)  ! not sure why this is needed

    INLET_VN   = 0.0d0
    INLET_VT   = 0.0d0
    INLET_RHO  = eos_state%rho
    INLET_RHOH = eos_state%rho * eos_state%h

    do comp=1,nspec
       if(spec_names(comp) .eq. "carbon-12") then
          INLET_RHOC12  = eos_state%rho * eos_state%xn(comp)

       else if(spec_names(comp) .eq. "magnesium-24") then
          INLET_RHOMG24 = eos_state%rho * eos_state%xn(comp)

       else if(spec_names(comp) .eq. "oxygen-16") then
          INLET_RHOO16  = eos_state%rho * eos_state%xn(comp)
       endif
    enddo
    INLET_TEMP = eos_state%T
    INLET_TRA = 0.0d0

    inlet_bc_initialized = .true.

  end subroutine set_inlet_bcs

end module inlet_bc_module
