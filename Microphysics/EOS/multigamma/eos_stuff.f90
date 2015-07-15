! This is a multi-gamma EOS.  Each species can have its own gamma, but
! otherwise, they all act as an ideal gas

! Note: at the moment, it is not clear what the proper expression for
! a multi-gamma entropy should be, so do not rely on the entropy.

module eos_module

  use bl_types
  use bl_space
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion, network_species_index
  use eos_type_module
  use eos_data_module

  implicit none

  private

  real(kind=dp_t), save, private :: smallt
  real(kind=dp_t), save, private :: smalld
  real(kind=dp_t), save :: gammas(nspec)

  logical, save, private :: initialized = .false.

  public eos_init, eos_finalize, eos
  public eos_input_rt, eos_input_rh, eos_input_tp, eos_input_rp, &
         eos_input_re, eos_input_ps, eos_input_ph
  public gammas

contains

  ! EOS initialization routine
  subroutine eos_init(small_temp, small_dens, gamma_in)

    use extern_probin_module, only: eos_gamma_default, &
                             species_a_name, species_a_gamma, &
                             species_b_name, species_b_gamma, &
                             species_c_name, species_c_gamma

    implicit none
 
    real(kind=dp_t), intent(in), optional :: small_temp
    real(kind=dp_t), intent(in), optional :: small_dens
    real(kind=dp_t), intent(in), optional :: gamma_in  ! this is ignored

    integer :: idx

    ! set the gammas for the species -- we have some runtime parameters
    ! that can override the default gammas for a few named species
    gammas(:) = eos_gamma_default

    print *, "in eos_init, species, gamma: "
    print *, "species a: ", trim(species_a_name), species_a_gamma
    print *, "species b: ", trim(species_b_name), species_b_gamma
    print *, "species c: ", trim(species_c_name), species_c_gamma

    idx = network_species_index(species_a_name)
    if (idx > 0) gammas(idx) = species_a_gamma

    idx = network_species_index(species_b_name)
    if (idx > 0) gammas(idx) = species_b_gamma

    idx = network_species_index(species_c_name)
    if (idx > 0) gammas(idx) = species_c_gamma


    ! Small temperature and density parameters
    smallt = ZERO
    if (present(small_temp)) smallt = small_temp

    smalld = ZERO
    if (present(small_dens)) smalld = small_dens

    initialized = .true.
 
  end subroutine eos_init



  !---------------------------------------------------------------------------
  ! The main interface
  !---------------------------------------------------------------------------
  subroutine eos(input, state, do_eos_diag, pt_index)

    use fundamental_constants_module, only: k_B, n_A, hbar

    implicit none

    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: state
    logical, optional, intent(in   ) :: do_eos_diag
    integer, optional, intent(in   ) :: pt_index(:)

    ! Local variables
    real(kind=dp_t) :: ymass
    real(kind=dp_t) :: sumY_gm1, sumYg_gm1
    real(kind=dp_t) :: sum_y, abar
    real(kind=dp_t) :: dens, temp

    logical :: eos_diag

    ! Get the mass of a nucleon from Avogadro's number.
    real(kind=dp_t), parameter :: m_nucleon = ONE / n_A

    real(kind=dp_t) :: dedX(nspec)

    integer :: k, n

    ! general sanity checks
    if (.not. initialized) call bl_error('EOS: not initialized')


    ! Get abar (= mu)
    sum_y = 0.d0

    do n = 1, nspec
       ymass = state % xn(n) / aion(n)
       sum_y = sum_y + ymass
    enddo

    abar = 1.d0/sum_y

    ! special gamma factors
    sumY_gm1 = ZERO
    sumYg_gm1 = ZERO

    do n = 1, nspec
       sumY_gm1 = sumY_gm1 + state % xn(n)/(aion(n)*(gammas(n)-ONE))
       sumYg_gm1 = sumYg_gm1 + state % xn(n)*gammas(n)/(aion(n)*(gammas(n)-ONE))
    enddo
       
    

    !-------------------------------------------------------------------------
    ! For all EOS input modes EXCEPT eos_input_rt, first compute dens
    ! and temp as needed from the inputs.
    !-------------------------------------------------------------------------

    select case (input)

    case (eos_input_rt)

       ! dens, temp and xmass are inputs

       ! We don't need to do anything here
       temp = state % T
       dens = state % rho


    case (eos_input_rh)

       ! dens, enthalpy, and xmass are inputs

       ! Solve for the temperature:
       ! h = e + p/rho = (p/rho)*[1 + 1/(gamma-1)] = (p/rho)*gamma/(gamma-1)
       dens = state % rho
       temp = (state % h * m_nucleon / k_B)/sumYg_gm1


    case (eos_input_tp)

       ! temp, pres, and xmass are inputs
          
       ! Solve for the density:
       ! p = rho k T / (abar m_nucleon)
       dens = state % p * abar * m_nucleon / (k_B * state % T)
       temp = state % T


    case (eos_input_rp)

       ! dens, pres, and xmass are inputs

       ! Solve for the temperature:
       ! p = rho k T / (mu m_nucleon)
       dens = state % rho
       temp = state % p * abar * m_nucleon / (k_B * state % rho)


    case (eos_input_re)

       ! dens, energy, and xmass are inputs

       ! Solve for the temperature
       ! e = k T / [(mu m_nucleon)*(gamma-1)]
       dens = state % rho
       temp = state % e * m_nucleon / (k_B * sumY_gm1)


    case (eos_input_ps)
       
       ! pressure entropy, and xmass are inputs
       call bl_error("eos_input_ps is not supported")


    case (eos_input_ph)

       ! pressure, enthalpy and xmass are inputs

       ! Solve for temperature and density
       dens = state % p * abar / state % h * sumYg_gm1
       temp = state % p * abar * m_nucleon / (k_B * dens)

    case default

       call bl_error('EOS: invalid input.')

    end select

    !-------------------------------------------------------------------------
    ! Now we have the density and temperature (and mass fractions /
    ! mu), regardless of the inputs.
    !-------------------------------------------------------------------------

    state % T   = temp
    state % rho = dens

    ! Compute the pressure simply from the ideal gas law, and the
    ! specific internal energy using the gamma-law EOS relation.
    state % p = dens*k_B*temp/(abar*m_nucleon)
    state % e = k_B*temp*sumY_gm1/m_nucleon

    
    ! enthalpy is h = e + p/rho
    state % h = state % e + state % p / dens
    
    ! entropy (per gram) -- this is wrong. Not sure what the expression
    ! is for a multigamma gas
    state % s = (k_B/(abar*m_nucleon))*(2.5_dp_t + &
         log( ( (abar*m_nucleon)**2.5/dens )*(k_B*temp)**1.5_dp_t / (2.0_dp_t*M_PI*hbar*hbar)**1.5_dp_t ) )


    ! Compute the thermodynamic derivatives and specific heats 
    state % dpdT = state % p / temp
    state % dpdr = state % p / dens

    state % dedT = state % e / temp
    state % dedr = ZERO
    state % dsdT = ZERO
    state % dsdr = ZERO

    state % cv = state % dedT
    state % cp = state % h / state % T

    state % gam1 = state % cp / state % cv

 
    ! Get dpdX, dedX, dhdX.
    state % dpdX(:) = state % rho * k_B * state % T / (m_nucleon * aion(:))
    dedX(:) = k_B * state % T / (m_nucleon * aion(:) * (gammas(:) - ONE))
  
    state % dhdX(:) = dedX(:) + (state % p / state % rho**2 - state % dedr) &
         *  state % dpdX(:) / state % dpdr
  

    ! sound speed
    state % cs = sqrt(state % gam1 * state % p / dens)

    return
  end subroutine eos

  subroutine eos_finalize()

  end subroutine eos_finalize

end module eos_module
