! This is a constant gamma equation of state, using and ideal gas.
!
! We compute the mean molecular weight, mu, based on the mass fractions
! of the different species.
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
!   1/mu = sum_k { (1 + Z_k) X_k / A_k }  = 1/mu_e + 1/mu
!
! It is expected that composition() is called the input eos_t 
! state before calling this eos routine.
!
! The ratio of specific heats (gamma) is allowed to vary.  NOTE: the
! expression for entropy is only valid for an ideal MONATOMIC gas
! (gamma = 5/3).

module actual_eos_module

  use bl_types
  use bl_space
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module

  implicit none

  private

  character (len=64) :: eos_name = "gamma_law_general"

  real (kind=dp_t), save :: gamma_const

  logical, save :: assume_neutral

  public actual_eos, actual_eos_init, eos_name

contains

  subroutine actual_eos_init

    use extern_probin_module, only: eos_gamma, eos_assume_neutral

    if (eos_gamma > ZERO) then
       gamma_const = eos_gamma
    else
       call bl_error("eos_gamma cannot be < 0")
    endif

    assume_neutral = eos_assume_neutral

  end subroutine actual_eos_init


  subroutine actual_eos(input, state)

    use fundamental_constants_module, only: k_B, n_A, hbar

    implicit none

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: state

    ! get the mass of a nucleon from Avogadro's number.
    real (kind=dp_t), parameter :: m_nucleon = ONE / n_A
    real (kind=dp_t), parameter :: fac = ONE / (TWO*M_PI*hbar*hbar)**1.5d0

    real (kind=dp_t) :: Tinv, rhoinv

    ! calculate mu

    if (assume_neutral) then
       state % mu = state % abar
    else
       state % mu = ONE / sum( (ONE + zion(:)) * state % xn(:) / aion(:) )
    endif

    !-------------------------------------------------------------------------
    ! for all EOS input modes EXCEPT eos_input_rt, first compute dens
    ! and temp as needed from the inputs
    !-------------------------------------------------------------------------

    select case (input)
       
    case (eos_input_rt)

       ! dens, temp, and xmass are inputs

       ! we don't need to do anything here
       continue

    case (eos_input_rh)

       ! dens, enthalpy, and xmass are inputs

       ! Solve for the temperature:
       ! h = e + p/rho = (p/rho)*[1 + 1/(gamma-1)] = (p/rho)*gamma/(gamma-1)

       state % T = (state % h * state % mu * m_nucleon / k_B)*(gamma_const - ONE)/gamma_const

    case (eos_input_tp ) 

       ! temp, pres, and xmass are inputs

       ! Solve for the density:
       ! p = rho k T / (mu m_nucleon)

       state%rho = state%p * state%mu * m_nucleon / (k_B * state%T)

    case (eos_input_rp )

       ! dens, pres, and xmass are inputs

       ! Solve for the temperature:
       ! p = rho k T / (mu m_nucleon)

       state%T = state%p * state%mu * m_nucleon / (k_B * state%rho)

    case (eos_input_re) 

       ! dens, energy, and xmass are inputs

       ! Solve for the temperature
       ! e = k T / [(mu m_nucleon)*(gamma-1)]

       state%T = state%e * state%mu * m_nucleon * (gamma_const - ONE) / k_B

    case (eos_input_ps)

       ! pressure, entropy, and xmass are inputs

       ! Solve for the temperature
       ! Invert Sackur-Tetrode eqn (below) using
       ! rho = p mu m_nucleon / (k T)

       state%T = state%p**(TWO/FIVE) * &
            ( TWO*M_PI*hbar*hbar/(state%mu*m_nucleon) )**(THREE/FIVE) * &
            dexp(TWO*state%mu*m_nucleon*state%s/(FIVE*k_B) - ONE) / k_B

       ! Solve for the density
       ! rho = p mu m_nucleon / (k T)

       state%rho = state%p*state%mu*m_nucleon/(k_B*state%T)

    case (eos_input_ph)

       call bl_error("EOS: eos_input_ph not implemented")


    case (eos_input_th)

       ! temperature, enthalpy, and xmass are inputs
       
       ! this system is underconstrained
       
       call bl_error("EOS: eos_input_th is not a valid input for the gamma law EOS")


    case default

       call bl_error("EOS: invalid input")

    end select

    !-------------------------------------------------------------------------
    ! now we have the density and temperature (and mass fractions /
    ! mu), regardless of the inputs.
    !-------------------------------------------------------------------------

    Tinv = ONE / state % T
    rhoinv = ONE / state % rho

    ! compute the pressure simply from the ideal gas law, and the
    ! specific internal energy using the gamma-law EOS relation
    state%p = state%rho * k_B * state%T / (state%mu * m_nucleon)
    state%e = state%p / (gamma_const - ONE) * rhoinv

    ! enthalpy is h = e + p/rho
    state%h = state%e + state%p * rhoinv

    ! entropy (per gram) of an ideal monoatomic gas (the Sackur-Tetrode equation)
    ! NOTE: this expression is only valid for gamma = 5/3.
    state%s = (k_B/(state%mu*m_nucleon))*(2.5_dp_t + &
         log( ( (state%mu*m_nucleon)**2.5_dp_t*rhoinv)*(k_B*state%T)**1.5_dp_t * fac) )

    ! compute the thermodynamic derivatives and specific heats
    state % dpdT = state%p * Tinv
    state % dpdr = state%p * rhoinv
    state % dedT = state%e * Tinv
    state % dedr = ZERO
    state % dsdT = 1.5_dp_t * (k_B / (state%mu * m_nucleon*state%T))
    state % dsdr = - (k_B / (state%mu* m_nucleon*state%rho) )
    state % dhdT = state % dedT + state % dpdT * rhoinv
    state % dhdr = ZERO

    state%cv = state%dedT
    state%cp = gamma_const*state%cv

    state%gam1 = gamma_const

    ! electron-specific stuff (really for the degenerate EOS)
    state%xne = 0.d0
    state%eta = 0.d0
    state%pele = 0.d0

    ! sound speed
    state%cs = sqrt(gamma_const * state%p * rhoinv)

    state % dpdA = - state % p * (ONE/state % abar)
    state % dedA = - state % e * (ONE/state % abar)

    if (assume_neutral) then
       state % dpdZ = ZERO
       state % dedZ = ZERO
    else
       state % dpdZ = state % p * (ONE/(ONE + state % zbar))
       state % dedZ = state % e * (ONE/(ONE + state % zbar))
    endif

  end subroutine actual_eos

end module actual_eos_module
