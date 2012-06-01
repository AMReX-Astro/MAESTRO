module eos_module

!  This implements the ideal gas law
!  P=d*T, where P is the pressure, d is the density
!  and T is the temperature.  The internal energy
!  e is given by e=1/(gamma-1) T.

  use bl_types
  use bl_error_module
  use bl_space
  use network, only: nspec
  use eos_type_module

  implicit none

  real(kind=dp_t), save :: gamma_const

  integer, parameter, public :: eos_input_rt = 1
  ! density, temperature are inputs

  integer, parameter, public :: eos_input_rh = 2
  ! density, enthalpy are inputs

  integer, parameter, public :: eos_input_tp = 3
  ! temperature, pressure are inputs

  integer, parameter, public :: eos_input_rp = 4
  ! density, pressure are inputs

  integer, parameter, public :: eos_input_re = 5
  ! density, internal energy are inputs

  integer, parameter, public :: eos_input_ps = 6
  ! pressure, entropy are inputs

  real(kind=dp_t), public :: xn_eos(nspec)
  real(kind=dp_t), public :: temp_eos
  real(kind=dp_t), public :: den_eos
  real(kind=dp_t), public :: p_eos
  real(kind=dp_t), public :: h_eos
! useless variables I need to include...
  real(kind=dp_t), public :: e_eos
  real(kind=dp_t), public :: abar_eos
  real(kind=dp_t), public :: zbar_eos
  real(kind=dp_t), public :: cv_eos
  real(kind=dp_t), public :: cp_eos
  real(kind=dp_t), public :: xne_eos
  real(kind=dp_t), public :: eta_eos
  real(kind=dp_t), public :: pele_eos
  real(kind=dp_t), public :: dpdt_eos
  real(kind=dp_t), public :: dpdr_eos
  real(kind=dp_t), public :: dedr_eos
  real(kind=dp_t), public :: dedt_eos
  real(kind=dp_t), public :: gam1_eos
  real(kind=dp_t), public ::   cs_eos
  real(kind=dp_t), public ::    s_eos
  real(kind=dp_t), public :: dsdt_eos
  real(kind=dp_t), public :: dsdr_eos
  real(kind=dp_t), public :: dpdX_eos(nspec)
  real(kind=dp_t), public :: dhdX_eos(nspec)
  real(kind=dp_t), public :: conduct_eos

  integer, public         :: pt_index_eos(MAX_SPACEDIM)

  interface eos
     module procedure eos_old
     module procedure eos_new
  end interface eos

contains

  subroutine eos_init(use_eos_coulomb, small_temp, small_dens, gamma_in)

    implicit none

    logical        , intent(in), optional :: use_eos_coulomb
    real(kind=dp_t), intent(in), optional :: small_temp, small_dens
    real(kind=dp_t), intent(in), optional :: gamma_in

    if (present(gamma_in)) then
       gamma_const = gamma_in
    else
       gamma_const = 1.666666d0
    end if

  end subroutine eos_init


  !---------------------------------------------------------------------------
  ! new interface
  !---------------------------------------------------------------------------
  subroutine eos_new(input, eos_state, do_eos_diag, pt_index)

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: eos_state
    logical,           intent(in   ) :: do_eos_diag
    integer, optional, intent(in   ) :: pt_index(:)

    call eos_old(input, eos_state%rho, eos_state%T, &
                 eos_state%xn, &
                 eos_state%p, eos_state%h, eos_state%e, &
                 eos_state%cv, eos_state%cp, eos_state%xne, &
                 eos_state%eta, eos_state%pele, &
                 eos_state%dpdT, eos_state%dpdr, &
                 eos_state%dedT, eos_state%dedr, &
                 eos_state%dpdX, eos_state%dhdX, &
                 eos_state%gam1, eos_state%cs, eos_state%s, &
                 eos_state%dsdT, eos_state%dsdr, &
                 do_eos_diag, pt_index)

  end subroutine eos_new


  !---------------------------------------------------------------------------
  ! The main interface -- this is used directly by MAESTRO
  !---------------------------------------------------------------------------
  subroutine eos_old(input, rho, temp, xmass, pres, h, &
                     eint, c_v, c_p, ne, eta, pele, &
                     dPdT, dPdR, dEdT, dEdR, dPdX, dhdX, &
                     gam1, cs, entropy, dsdT, dsdR, do_eos_diag, pt_index)

!  input = 1 means rho, temp are inputs
!        = 2 means temp, pres are inputs
!        = 3 means rho, pres are inputs

    implicit none

    integer, intent(in) :: input

    integer :: k,n

    real(kind=dp_t), intent(inout) :: rho, temp, pres, eint, h

! all the silly variables that I will not use
    logical do_eos_diag
    real(kind=dp_t) :: xmass(nspec), c_v, c_p
    real(kind=dp_t) :: ne, eta, pele
    real(kind=dp_t) :: dPdT, dPdR, dEdT
    real(kind=dp_t) :: dEdR, gam1, entropy
    real(kind=dp_t) :: cs, dPdX(nspec)
    real(kind=dp_t) :: dhdX(nspec), dsdT, dsdR
    integer, optional, intent(in   ) :: pt_index(:)

    ! first find rho, pres, and temp

    if (input .EQ. eos_input_rt) then

       ! input = 1: rho, temp are inputs

       pres = rho*temp

    else if (input .EQ. eos_input_rh) then

       ! input = 2: rho, h are inputs

       temp = h*(gamma_const - 1.0)/gamma_const
       pres = temp*rho

    else if (input .EQ. eos_input_tp) then

       ! input = 3: temp, pres are inputs

       rho = pres/temp

    else if (input .EQ. eos_input_rp) then

       !input = 4: density, pressure are inputs

       temp = pres/rho

    else if (input .EQ. eos_input_re) then

       !input = 5: density, internal energy are inputs

       temp = (gamma_const-1.0)*eint
       pres = temp*rho

    else

       call bl_error('invalid input')

    end if

    ! then calculate everything else

    h = gamma_const/(gamma_const - 1.0) * temp
    eint = temp/(gamma_const - 1.0)
    gam1 = gamma_const
    cs = sqrt(gamma_const*pres/rho)
    dPdT = pres/temp
    dPdR = pres/rho
    dedT = eint/temp
    dedR = 0.d0
    dsdT = 0.d0
    dsdR = 0.d0
    entropy = pres/rho**gam1
    pele = 0.d0
    eta  = 0.d0
    ne   = 0.d0

    do n = 1, nspec
       dPdX(n) = 0.d0
       dhdX(n) = 0.d0
    enddo

    c_v = dedT
    c_p = gamma_const*c_v


  end subroutine eos_old

end module eos_module
