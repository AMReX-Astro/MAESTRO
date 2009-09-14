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


module eos_module

  use bl_types
  use bl_constants_module
  use network, only: nspec, aion, zion

  implicit none

  integer, parameter :: NP = 1
  integer, parameter :: npts = 1

  real(kind=dp_t) :: xn_eos(NP,nspec)

  real(kind=dp_t) :: temp_eos(NP)
  real(kind=dp_t) :: den_eos(NP)
  real(kind=dp_t) :: abar_eos(NP)
  real(kind=dp_t) :: zbar_eos(NP)
  real(kind=dp_t) :: e_eos(NP)
  real(kind=dp_t) :: p_eos(NP)
  real(kind=dp_t) :: h_eos(NP)
  real(kind=dp_t) :: cv_eos(NP)
  real(kind=dp_t) :: cp_eos(NP)
  real(kind=dp_t) :: xne_eos(NP)
  real(kind=dp_t) :: eta_eos(NP)
  real(kind=dp_t) :: pele_eos(NP)
  real(kind=dp_t) :: dpdt_eos(NP)
  real(kind=dp_t) :: dpdr_eos(NP)
  real(kind=dp_t) :: dedr_eos(NP)
  real(kind=dp_t) :: dedt_eos(NP)
  real(kind=dp_t) :: gam1_eos(NP)
  real(kind=dp_t) ::   cs_eos(NP)
  real(kind=dp_t) ::    s_eos(NP)
  real(kind=dp_t) :: dsdt_eos(NP)
  real(kind=dp_t) :: dsdr_eos(NP)
  real(kind=dp_t) :: dpdX_eos(NP,nspec)
  real(kind=dp_t) :: dhdX_eos(NP,nspec)
  real(kind=dp_t) :: conduct_eos(NP)

  integer, parameter :: eos_input_rt = 1   ! density, temperature are inputs
  integer, parameter :: eos_input_rh = 2   ! density, enthalpy are inputs
  integer, parameter :: eos_input_tp = 3   ! temperature, pressure are inputs
  integer, parameter :: eos_input_rp = 4   ! density, pressure are inputs
  integer, parameter :: eos_input_re = 5   ! density, internal energy are inputs
 
  logical :: do_diag

  real(kind=dp_t), save, private :: smallt
  real(kind=dp_t), save, private :: smalld
  real(kind=dp_t), save          :: gamma

  logical, save, private :: initialized = .false.

  logical, parameter :: eos_assume_neutral = .true.

  private nspec, aion, zion

contains

  ! EOS initialization routine -- this is used by both MAESTRO and Castro
  ! For this general EOS, this calls helmeos_init() which reads in the 
  ! table with the electron component's properties.
  subroutine eos_init(use_eos_coulomb, small_temp, small_dens, gamma_in)

    implicit none
 
    logical        , intent(in), optional :: use_eos_coulomb
    real(kind=dp_t), intent(in), optional :: small_temp
    real(kind=dp_t), intent(in), optional :: small_dens
    real(kind=dp_t), intent(in), optional :: gamma_in
 
    ! constant ratio of specific heats
    if (present(gamma_in)) then
       gamma = gamma_in
    else
       gamma = 1.4d0
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

  !---------------------------------------------------------------------------
  ! The main interface -- this is used directly by MAESTRO
  !---------------------------------------------------------------------------
  subroutine eos(input, dens, temp, &
                 npoints, nspecies, &
                 xmass, &
                 pres, enthalpy, eint, &
                 c_v, c_p, ne, eta, pele, &
                 dPdT, dPdR, dEdT, dEdR, &
                 dPdX, dhdX, &
                 gam1, cs, entropy, &
                 dsdT, dsdR, &
                 do_eos_diag)

    use bl_error_module
    use fundamental_constants_module, only: k_B, n_A

! dens     -- mass density (g/cc)
! temp     -- temperature (K)
! npoints     -- the number of elements in input/output arrays
! nspecies -- the number of isotopes
! xmass    -- the mass fractions of the individual isotopes
! pres     -- the pressure (dyn/cm**2)
! enthalpy -- the enthalpy (erg/g)
! eint     -- the internal energy (erg/g)
! c_v      -- specific heat at constant volume
! c_p      -- specific heat at constant pressure
! ne       -- number density of electrons + positrons
! eta      -- degeneracy parameter
! pele     -- electron pressure + positron pressure
! dPdT     -- d pressure/ d temperature
! dPdR     -- d pressure/ d density
! dEdT     -- d energy/ d temperature
! dEdR     -- d energy/ d density
! dPdX     -- d pressure / d xmass(k)
! dhdX     -- d enthalpy / d xmass(k)  -- AT CONSTANT PRESSURE!!!
! gam1     -- first adiabatic index (d log P/ d log rho) |_s
! cs       -- sound speed -- note that this is the non-relativistic one
!             (we compute it in this wrapper as sqrt(gam1 p /rho) instead
!             of taking the relativistic version from helmeos.
! entropy  -- entropy (erg/g/K)
!
! input = 1 means dens, temp    , and xmass are inputs, return enthalpy, eint
!       = 2 means dens, enthalpy, and xmass are inputs, return temp    , eint
!                (note, temp should be filled with an initial guess)
!       = 3 means temp, pres    , and xmass are inputs, return dens    , etc
!       = 4 means dens, pres    , and xmass are inputs, return temp    , etc
!       = 5 means dens, eint    , and xmass are inputs, return temp    , etc
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

    implicit none

    logical do_eos_diag
    integer, intent(in) :: input
    integer, intent(in) :: npoints
    integer, intent(in) :: nspecies

    real(kind=dp_t) :: dens(npoints), temp(npoints)
    real(kind=dp_t) :: xmass(npoints,nspecies)
    real(kind=dp_t) :: pres(npoints), enthalpy(npoints), eint(npoints)
    real(kind=dp_t) :: c_v(npoints), c_p(npoints)
    real(kind=dp_t) :: ne(npoints), eta(npoints), pele(npoints)
    real(kind=dp_t) :: dPdT(npoints), dPdR(npoints), dedT(npoints), dedR(npoints)
    real(kind=dp_t) :: gam1(npoints), entropy(npoints), cs(npoints)
    real(kind=dp_t) :: dPdX(npoints,nspecies), dedX(npoints,nspecies), dhdX(npoints,nspecies)
    real(kind=dp_t) :: dsdT(npoints), dsdR(npoints)

    ! local variables
    real(kind=dp_t) :: ymass(npoints,nspecies)    
    real(kind=dp_t) :: mu(npoints)
    real(kind=dp_t) :: dmudX, sum_y

    ! get the mass of a nucleon from Avogadro's number.
    real(kind=dp_t), parameter :: m_nucleon = 1.d0/n_A

    integer :: k, n

    ! general sanity checks
    if (.not. initialized) call bl_error('EOS: not initialized')
      
    if (nspecies /= nspec) then
       call bl_error('EOS: too many species')
    endif

    if (npoints > NP) then
       call bl_error('EOS: eos called with too large of a vector size')
    endif

    !-------------------------------------------------------------------------
    ! compute mu -- the mean molecular weight
    !-------------------------------------------------------------------------
    if (eos_assume_neutral) then
       ! assume completely neutral atoms

       do k = 1, npoints
          sum_y  = 0.d0
          
          do n = 1, nspecies
             ymass(k,n) = xmass(k,n)/aion(n)
             sum_y = sum_y + ymass(k,n)
          enddo
          
          mu(k) = 1.d0/sum_y
       enddo

    else
       ! assume completely ionized species

       do k = 1, npoints
          sum_y  = 0.d0
          
          do n = 1, nspecies
             ymass(k,n) = xmass(k,n)*(1.d0 + zion(n))/aion(n)
             sum_y = sum_y + ymass(k,n)
          enddo
          
          mu(k) = 1.d0/sum_y
       enddo

    endif

    !-------------------------------------------------------------------------
    ! for all EOS input modes EXCEPT eos_input_rt, first compute dens
    ! and temp as needed from the inputs
    !-------------------------------------------------------------------------
    if (input .EQ. eos_input_rh) then

       ! dens, enthalpy, and xmass are inputs
       do k = 1, npoints

          ! Solve for the temperature:
          ! h = e + p/rho = (p/rho)*[1 + 1/(gamma-1)] = (p/rho)*gamma/(gamma-1)
          temp(k) = (enthalpy(k)*mu(k)*m_nucleon/k_B)*(gamma - 1.0)/gamma
       enddo


    else if (input .EQ. eos_input_tp ) then

       ! temp, pres, and xmass are inputs
       do k = 1, npoints
          
          ! Solve for the density:
          ! p = rho k T / (mu m_nucleon)
          dens(k) = pres(k)*mu(k)*m_nucleon/(k_B*temp(k))
       enddo


    else if (input .EQ. eos_input_rp ) then

       ! dens, pres, and xmass are inputs
       do k = 1, npoints

          ! Solve for the temperature:
          ! p = rho k T / (mu m_nucleon)
          temp(k) = pres(k)*mu(k)*m_nucleon/(k_B*dens(k))
       enddo


    else if (input .EQ. eos_input_re) then

       ! dens, energy, and xmass are inputs
       do k = 1, npoints

          ! Solve for the temperature
          ! e = k T / [(mu m_nucleon)*(gamma-1)]
          temp(k) = eint(k)*mu(k)*m_nucleon*(gamma-1.0)/k_B

       enddo
    
    endif

    !-------------------------------------------------------------------------
    ! now we have the density and temperature (and mass fractions /
    ! mu), regardless of the inputs.
    !-------------------------------------------------------------------------
    do k = 1, npoints

       ! compute the pressure simply from the ideal gas law, and the
       ! specific internal energy using the gamma-law EOS relation
       pres(k) = dens(k)*k_B*temp(k)/(mu(k)*m_nucleon)
       eint(k) = pres(k)/(gamma - 1.0)/dens(k)

       ! enthalpy is h = e + p/rho
       enthalpy(k) = eint(k) + pres(k)/dens(k)

       ! entropy ?
       entropy(k) = 0.d0

       ! compute the thermodynamic derivatives and specific heats 
       dPdT(k) = pres(k)/temp(k)
       dPdR(k) = pres(k)/dens(k)
       dedT(k) = eint(k)/temp(k)
       dedR(k) = 0.d0
       dsdT(k) = 0.d0
       dsdR(k) = 0.d0

       c_v(k) = dedT(k)
       c_p(k) = gamma*c_v(k)

       gam1(k) = gamma

       do n = 1, nspecies

          ! the species only come into p and e (and therefore h)
          ! through mu, so first compute dmu/dX
          if (eos_assume_neutral) then
             dmudX =  (mu(k)/aion(n))*(aion(n) - mu(k))
          else
             dmudX =  (mu(k)/aion(n))*(aion(n) - mu(k)*(1.0d0 + zion(n)))
          endif

          dPdX(k,n) = -(pres(k)/mu(k))*dmudX
          dedX(k,n) = -(eint(k)/mu(k))*dmudX
          
          ! dhdX is at constant pressure -- see paper III for details
          dhdX(k,n) = dedX(k,n) + &
               (pres(k)/dens(k)**2 - dedR(k))*dPdX(k,n)/dPdr(k)
       enddo

       ! electron-specific stuff (really for the degenerate EOS)
       ne(k) = 0.d0
       eta(k) = 0.d0
       pele(k) = 0.d0

       ! sound speed
       cs(k) = sqrt(gamma*pres(k)/dens(k))

    enddo

    return
  end subroutine eos

end module eos_module
