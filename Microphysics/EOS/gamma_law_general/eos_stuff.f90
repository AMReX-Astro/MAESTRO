! This is a constant gamma equation of state, using and ideal gas.
!
! We compute the mean molecular weight, mu, based on the mass fractions
! of the different species.
!
!   NOTE: in the helmholtz EOS, we use Abar, which is the weighted
!   ion mass.  Abar does not include any free electron contribution.
!
! There are 2 ways to compute mu, depending on whether the gas is
! completely ionized or not.
!
! For a neutral gas, the mean molecular weight is:
!
!   1/mu = sum_k { X_k / A_k }
!
! For a completely ionized gas, the mean molecular weight is:
!
!   1/mu = sum_k { (1 + Z_k) X_k / A_k }
!
! At the moment, we will select between these 2 ionization extremes
! (completely neutral vs. completely ionized) by a hard-coded
! parameter, eos_assume_neutral
!
! The ratio of specific heats (gamma) is allowed to vary.  NOTE: the
! expression for entropy is only valid for an ideal MONATOMIC gas
! (gamma = 5/3).

module eos_module

  use bl_types
  use bl_space
  use bl_constants_module, only: M_PI, ONE, TWO
  use network, only: nspec, aion, zion
  use eos_type_module
  use eos_data_module

  implicit none

  private

  real(kind=dp_t), save, private :: smallt
  real(kind=dp_t), save, private :: smalld
  real(kind=dp_t), save, public  :: gamma_const

  logical, save, private :: initialized = .false.

  logical, parameter, private :: eos_assume_neutral = .true.

  private nspec, aion, zion

  public eos_init, eos_finalize, eos
  public eos_input_rt, eos_input_rh, eos_input_tp, eos_input_rp, &
         eos_input_re, eos_input_ps, eos_input_ph

  interface eos
     module procedure eos_new
  end interface eos

contains

  ! EOS initialization routine
  subroutine eos_init(small_temp, small_dens, gamma_in)

    implicit none

    real(kind=dp_t), intent(in), optional :: small_temp
    real(kind=dp_t), intent(in), optional :: small_dens
    real(kind=dp_t), intent(in), optional :: gamma_in

    ! constant ratio of specific heats
    if (present(gamma_in)) then
       gamma_const = gamma_in
    else
       gamma_const = 5.d0/3.d0
    end if

    ! small temperature and density parameters
    if (present(small_temp)) then
       smallt = small_temp
    else
       smallt = 0.d0
    endif

    if (present(small_dens)) then
       smalld = small_dens
    else
       smalld = 0.d0
    endif

    initialized = .true.

  end subroutine eos_init


  subroutine eos_finalize()

  end subroutine eos_finalize


  !---------------------------------------------------------------------------
  ! new interface
  !---------------------------------------------------------------------------
  subroutine eos_new(input, eos_state, do_eos_diag, pt_index)

    use bl_error_module
    use fundamental_constants_module, only: k_B, n_A, hbar

    implicit none

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: eos_state
    logical, optional, intent(in   ) :: do_eos_diag
    integer, optional, intent(in   ) :: pt_index(:)

    double precision, parameter :: fac = ONE / (TWO*M_PI*hbar*hbar)**1.5d0


! All state information comes in through the eos_t derived type

! rho      -- mass density (g/cc)
! T        -- temperature (K)
! xn       -- the mass fractions of the individual isotopes
! p        -- the pressure (dyn/cm**2)
! h        -- the enthalpy (erg/g)
! e        -- the internal energy (erg/g)
! cv       -- specific heat at constant volume
! cp       -- specific heat at constant pressure
! xne      -- number density of electrons + positrons
! eta      -- degeneracy parameter
! pele     -- electron pressure + positron pressure
! dpdT     -- d pressure/ d temperature
! dpdr     -- d pressure/ d density
! dedT     -- d energy/ d temperature
! dedr     -- d energy/ d density
! dpdX     -- d pressure / d xmass(k)
! dhdX     -- d enthalpy / d xmass(k)  -- AT CONSTANT PRESSURE!!!
! gam1     -- first adiabatic index (d log P/ d log rho) |_s
! cs       -- sound speed -- note that this is the non-relativistic one
!             (we compute it in this wrapper as sqrt(gam1 p /rho) instead
!             of taking the relativistic version from helmeos.
! s        -- entropy (erg/g/K)  NOTE: presently the entropy expression is
!             valid only for an ideal MONATOMIC gas (gamma = 5/3).
!
! input = 1 means dens, temp    , and xmass are inputs, return enthalpy, eint
!       = 2 means dens, enthalpy, and xmass are inputs, return temp    , eint
!                (note, temp should be filled with an initial guess)
!       = 3 means temp, pres    , and xmass are inputs, return dens    , etc
!       = 4 means dens, pres    , and xmass are inputs, return temp    , etc
!       = 5 means dens, eint    , and xmass are inputs, return temp    , etc
!       = 6 means pres, entr    , and xmass are inputs, return temp    , etc
!
!
! derivatives wrt X_k:
!
!   For an ideal gas, the thermodynamic quantities only depend on composition
!   through the mean molecular weight, mu.
!
!   Using the chain rule:
!
!   dp/dX_k = dp/d(mu) d(mu)/dX_k
!


    ! local variables
    real(kind=dp_t) :: ymass(nspec)
    real(kind=dp_t) :: mu
    real(kind=dp_t) :: dmudX, sum_y

    real(kind=dp_t) :: dedX(nspec)

    ! get the mass of a nucleon from Avogadro's number.
    real(kind=dp_t), parameter :: m_nucleon = 1.d0/n_A

    integer :: k, n

    ! general sanity checks
    if (.not. initialized) call bl_error('EOS: not initialized')


    !-------------------------------------------------------------------------
    ! compute mu -- the mean molecular weight
    !-------------------------------------------------------------------------
    if (eos_assume_neutral) then
       ! assume completely neutral atoms

       sum_y  = 0.d0

       do n = 1, nspec
          ymass(n) = eos_state%xn(n)/aion(n)
          sum_y = sum_y + ymass(n)
       enddo

       mu = 1.d0/sum_y

    else
       ! assume completely ionized species

       sum_y  = 0.d0

       do n = 1, nspec
          ymass(n) = eos_state%xn(n)*(1.d0 + zion(n))/aion(n)
          sum_y = sum_y + ymass(n)
       enddo

       mu = 1.d0/sum_y

    endif

    !-------------------------------------------------------------------------
    ! for all EOS input modes EXCEPT eos_input_rt, first compute dens
    ! and temp as needed from the inputs
    !-------------------------------------------------------------------------
    if (input .EQ. eos_input_rh) then

       ! dens, enthalpy, and xmass are inputs

       ! Solve for the temperature:
       ! h = e + p/rho = (p/rho)*[1 + 1/(gamma-1)] = (p/rho)*gamma/(gamma-1)
       eos_state%T = (eos_state%h*mu*m_nucleon/k_B)*(gamma_const - 1.0_dp_t)/gamma_const


    else if (input .EQ. eos_input_tp ) then

       ! temp, pres, and xmass are inputs

       ! Solve for the density:
       ! p = rho k T / (mu m_nucleon)
       eos_state%rho = eos_state%p*mu*m_nucleon/(k_B*eos_state%T)


    else if (input .EQ. eos_input_rp ) then

       ! dens, pres, and xmass are inputs

       ! Solve for the temperature:
       ! p = rho k T / (mu m_nucleon)
       eos_state%T = eos_state%p*mu*m_nucleon/(k_B*eos_state%rho)


    else if (input .EQ. eos_input_re) then

       ! dens, energy, and xmass are inputs

       ! Solve for the temperature
       ! e = k T / [(mu m_nucleon)*(gamma-1)]
       eos_state%T = eos_state%e*mu*m_nucleon*(gamma_const-1.0_dp_t)/k_B


    else if (input .EQ. eos_input_ps) then

       ! pressure and entropy are inputs

       ! Solve for the temperature
       ! Invert Sackur-Tetrode eqn (below) using
       ! rho = p mu m_nucleon / (k T)
       eos_state%T = eos_state%p**(2.0_dp_t/5.0_dp_t) * &
            ( 2.0_dp_t*M_PI*hbar*hbar/(mu*m_nucleon) )**(3.0_dp_t/5.0_dp_t) * &
            dexp(2.0_dp_t*mu*m_nucleon*eos_state%s/(5.0_dp_t*k_B) - 1.0_dp_t) / &
            k_B

       ! Solve for the density
       ! rho = p mu m_nucleon / (k T)
       eos_state%rho = eos_state%p*mu*m_nucleon/(k_B*eos_state%T)

    else if (input .EQ. eos_input_ph) then
       call bl_error("ERROR: eos_input_ph not implemented")

    endif

    !-------------------------------------------------------------------------
    ! now we have the density and temperature (and mass fractions /
    ! mu), regardless of the inputs.
    !-------------------------------------------------------------------------

    ! compute the pressure simply from the ideal gas law, and the
    ! specific internal energy using the gamma-law EOS relation
    eos_state%p = eos_state%rho*k_B*eos_state%T/(mu*m_nucleon)
    eos_state%e = eos_state%p/(gamma_const - 1.0_dp_t)/eos_state%rho

    ! enthalpy is h = e + p/rho
    eos_state%h = eos_state%e + eos_state%p/eos_state%rho

    ! entropy (per gram) of an ideal monoatomic gas (the Sackur-Tetrode equation)
    ! NOTE: this expression is only valid for gamma = 5/3.
    eos_state%s = (k_B/(mu*m_nucleon))*(2.5_dp_t + &
         log( ( (mu*m_nucleon)**2.5/eos_state%rho )*(k_B*eos_state%T)**1.5_dp_t * fac) )

    ! compute the thermodynamic derivatives and specific heats
    eos_state%dpdT = eos_state%p/eos_state%T
    eos_state%dpdr = eos_state%p/eos_state%rho
    eos_state%dedT = eos_state%e/eos_state%T
    eos_state%dedr = 0.d0
    eos_state%dsdT = 1.5_dp_t * (k_B / (mu*m_nucleon*eos_state%T))
    eos_state%dsdR = - (k_B / (mu* m_nucleon*eos_State%rho) )

    eos_state%cv = eos_state%dedT
    eos_state%cp = gamma_const*eos_state%cv

    eos_state%gam1 = gamma_const

    do n = 1, nspec

       ! the species only come into p and e (and therefore h)
       ! through mu, so first compute dmu/dX
       !
       ! NOTE: an extra, constant term appears in dmudx, which
       ! results from writing mu = sum { X_k} / sum {X_k / A_k}
       ! (for the neutral, analogous for the ionized).  The
       ! numerator is simply 1, but we can differentiate
       ! wrt it, giving the constant mu(k) term in dmudx.  Since
       ! dPdX only appears in a sum over species creation rate
       ! (omegadot) and sum{omegadot} = 0, this term has no effect.
       ! If is added simply for completeness.

       if (eos_assume_neutral) then
          dmudX =  (mu/aion(n))*(aion(n) - mu)
       else
          dmudX =  (mu/aion(n))*(aion(n) - mu*(1.0_dp_t + zion(n)))
       endif

       eos_state%dpdX(n) = -(eos_state%p/mu)*dmudX
       dedX(n) = -(eos_state%e/mu)*dmudX

       ! dhdX is at constant pressure -- see paper III for details
       eos_state%dhdX(n) = dedX(n) + &
            (eos_state%p/eos_state%rho**2 - eos_state%dedr)*eos_state%dpdX(n)/eos_state%dpdr
    enddo

    ! electron-specific stuff (really for the degenerate EOS)
    eos_state%xne = 0.d0
    eos_state%eta = 0.d0
    eos_state%pele = 0.d0

    ! sound speed
    eos_state%cs = sqrt(gamma_const*eos_state%p/eos_state%rho)

    return
  end subroutine eos_new

end module eos_module
