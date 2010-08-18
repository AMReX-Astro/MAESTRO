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
  use bl_constants_module, only: M_PI, ONE
  use network, only: nspec, aion, zion

  implicit none

  private

  integer, parameter, public :: NP = 1
  integer, parameter, public :: npts = 1

  real(kind=dp_t), public :: xn_eos(NP,nspec)
  real(kind=dp_t), public :: temp_eos(NP)
  real(kind=dp_t), public :: den_eos(NP)
  real(kind=dp_t), public :: abar_eos(NP)
  real(kind=dp_t), public :: zbar_eos(NP)
  real(kind=dp_t), public :: e_eos(NP)
  real(kind=dp_t), public :: p_eos(NP)
  real(kind=dp_t), public :: h_eos(NP)
  real(kind=dp_t), public :: cv_eos(NP)
  real(kind=dp_t), public :: cp_eos(NP)
  real(kind=dp_t), public :: xne_eos(NP)
  real(kind=dp_t), public :: eta_eos(NP)
  real(kind=dp_t), public :: pele_eos(NP)
  real(kind=dp_t), public :: dpdt_eos(NP)
  real(kind=dp_t), public :: dpdr_eos(NP)
  real(kind=dp_t), public :: dedr_eos(NP)
  real(kind=dp_t), public :: dedt_eos(NP)
  real(kind=dp_t), public :: gam1_eos(NP)
  real(kind=dp_t), public ::   cs_eos(NP)
  real(kind=dp_t), public ::    s_eos(NP)
  real(kind=dp_t), public :: dsdt_eos(NP)
  real(kind=dp_t), public :: dsdr_eos(NP)
  real(kind=dp_t), public :: dpdX_eos(NP,nspec)
  real(kind=dp_t), public :: dhdX_eos(NP,nspec)
  real(kind=dp_t), public :: conduct_eos(NP)

  integer, public         :: pt_index_eos(MAX_SPACEDIM)

  integer, parameter, public :: eos_input_rt = 1   ! density, temperature are inputs
  integer, parameter, public :: eos_input_rh = 2   ! density, enthalpy are inputs
  integer, parameter, public :: eos_input_tp = 3   ! temperature, pressure are inputs
  integer, parameter, public :: eos_input_rp = 4   ! density, pressure are inputs
  integer, parameter, public :: eos_input_re = 5   ! density, internal energy are inputs
  integer, parameter, public :: eos_input_ps = 6   ! pressure, entropy are inputs
 

  real(kind=dp_t), save, private :: smallt
  real(kind=dp_t), save, private :: smalld
  real(kind=dp_t), save, public  :: gamma_const

  logical, save, private :: initialized = .false.

  logical, parameter, private :: eos_assume_neutral = .true.

  private nspec, aion, zion

  public eos_init, eos_get_small_temp, eos_get_small_dens, eos_given_ReX, &
       eos_e_given_RPX, eos_S_given_ReX, eos_given_RTX, eos_dpdr_given_RTX, &
       eos_given_TPX, eos_given_PSX, eos

contains

  ! EOS initialization routine -- this is used by both MAESTRO and Castro
  subroutine eos_init(use_eos_coulomb, small_temp, small_dens, gamma_in)

    implicit none
 
    logical        , intent(in), optional :: use_eos_coulomb
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


  !---------------------------------------------------------------------------
  ! Castro interfaces 
  !---------------------------------------------------------------------------
  subroutine eos_get_small_temp(small_temp_out)
 
    real(kind=dp_t), intent(out) :: small_temp_out
 
    small_temp_out = smallt
 
  end subroutine eos_get_small_temp
 
  subroutine eos_get_small_dens(small_dens_out)
 
    real(kind=dp_t), intent(out) :: small_dens_out
 
    small_dens_out = smalld
 
  end subroutine eos_get_small_dens

  subroutine eos_given_ReX(G, P, C, T, dpdr_e, dpde, R, e, X, pt_index)

    ! note: here, dpdr_e is partial p / partial rho at constant e   
    !       and   dpde is partial p / partial e   at constant rho 


     ! In/out variables
     real(kind=dp_t), intent(  out) :: G, P, C, dpdr_e, dpde
     real(kind=dp_t), intent(inout) :: T
     real(kind=dp_t), intent(in   ) :: R, e, X(:)
     integer, optional, intent(in   ) :: pt_index(:)

     ! Local variables
     logical :: do_diag

     do_diag = .false.

     temp_eos(1) = T
      den_eos(1) = R
        e_eos(1) = e
      xn_eos(1,1:nspec) = X(1:nspec)

     call eos(eos_input_re, den_eos, temp_eos, &
              npts, &
              xn_eos, &
              p_eos, h_eos, e_eos, &
              cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
              dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
              dpdX_eos, dhdX_eos, &
              gam1_eos, cs_eos, s_eos, &
              dsdt_eos, dsdr_eos, &
              do_diag)

    G  = gam1_eos(1)
    P  =    p_eos(1)
    C  =   cs_eos(1)
    T  = temp_eos(1)
    dpdr_e = dpdr_eos(1) - dpdt_eos(1)*dedr_eos(1)/dedt_eos(1)
    dpde = dpdt_eos(1) / dedt_eos(1)

  end subroutine eos_given_ReX

  subroutine eos_e_given_RPX(e, T, R, P, X, pt_index)

     ! In/out variables
     real(kind=dp_t), intent(  out) :: e
     real(kind=dp_t), intent(in   ) :: R, p, X(:)
     real(kind=dp_t), intent(inout) :: T
     integer, optional, intent(in   ) :: pt_index(:)

     ! Local variables
     logical :: do_diag

     do_diag = .false.

     temp_eos(1) = T
      den_eos(1) = R
        p_eos(1) = P
      xn_eos(1,1:nspec) = X(1:nspec)

     call eos(eos_input_rp, den_eos, temp_eos, &
              npts, &
              xn_eos, &
              p_eos, h_eos, e_eos, &
              cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
              dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
              dpdX_eos, dhdX_eos, &
              gam1_eos, cs_eos, s_eos, &
              dsdt_eos, dsdr_eos, &
              do_diag)

    e  =    e_eos(1)
    T  = temp_eos(1)

  end subroutine eos_e_given_RPX

  subroutine eos_S_given_ReX(S, R, e, T, X, pt_index)

     implicit none

     ! In/out variables
     real(kind=dp_t), intent(  out) :: S
     real(kind=dp_t), intent(in   ) :: R, e, T, X(:)
     integer, optional, intent(in   ) :: pt_index(:)

     ! Local variables
     logical :: do_diag

     do_diag = .false.

     temp_eos(1) = T
      den_eos(1) = R
        e_eos(1) = e
      xn_eos(1,1:nspec) = X(1:nspec)

     call eos(eos_input_re, den_eos, temp_eos, &
              npts, &
              xn_eos, &
              p_eos, h_eos, e_eos, &
              cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
              dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
              dpdX_eos, dhdX_eos, &
              gam1_eos, cs_eos, s_eos, &
              dsdt_eos, dsdr_eos, &
              do_diag)

    S  = s_eos(1)

  end subroutine eos_S_given_ReX

  subroutine eos_given_RTX(e, P, R, T, X, pt_index)

     ! In/out variables
     real(kind=dp_t), intent(  out) :: e, P
     real(kind=dp_t), intent(in   ) :: R, T, X(:)
     integer, optional, intent(in   ) :: pt_index(:)

     ! Local variables
     logical :: do_diag

     do_diag = .false.

      den_eos(1) = R
     temp_eos(1) = T
      xn_eos(1,1:nspec) = X(1:nspec)

     call eos(eos_input_rt, den_eos, temp_eos, &
              npts, &
              xn_eos, &
              p_eos, h_eos, e_eos, &
              cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
              dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
              dpdX_eos, dhdX_eos, &
              gam1_eos, cs_eos, s_eos, &
              dsdt_eos, dsdr_eos, &
              do_diag)

    P  =    p_eos(1)
    e  =    e_eos(1)

  end subroutine eos_given_RTX

  subroutine eos_dpdr_given_RTX(e, P, R, T, X, dpdr, pt_index)

    ! note: here, dpdr is partial p / partial rho at constant T
    ! this is different than the dpdr_e that Castro uses for source
    ! terms in the primitive variable formulation.

    ! In/out variables
    real(kind=dp_t), intent(  out) :: e, P, dpdr
    real(kind=dp_t), intent(in   ) :: R, T, X(:)
    integer, optional, intent(in   ) :: pt_index(:)

    ! Local variables
    logical :: do_diag

    do_diag = .false.

    den_eos(1) = R
    temp_eos(1) = T
    xn_eos(1,1:nspec) = X(1:nspec)

    call eos(eos_input_rt, den_eos, temp_eos, &
             npts, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             do_diag)

    P  =    p_eos(1)
    e  =    e_eos(1)
    dpdr =  dpdr_eos(1)

  end subroutine eos_dpdr_given_RTX

  subroutine eos_given_TPX(e, P, R, T, X, pt_index)

     ! In/out variables
     real(kind=dp_t), intent(  out) :: e, R
     real(kind=dp_t), intent(in   ) :: P, T, X(:)
     integer, optional, intent(in   ) :: pt_index(:)

     ! Local variables
     logical :: do_diag

     do_diag = .false.

     ! An initial guess of density needs to be given
     den_eos(1) = R 
     p_eos(1) = P
     temp_eos(1) = T
     xn_eos(1,1:nspec) = X(1:nspec)

     call eos(eos_input_tp, den_eos, temp_eos, &
              npts, &
              xn_eos, &
              p_eos, h_eos, e_eos, &
              cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
              dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
              dpdX_eos, dhdX_eos, &
              gam1_eos, cs_eos, s_eos, &
              dsdt_eos, dsdr_eos, &
              do_diag)

    R  =    den_eos(1)
    e  =    e_eos(1)

  end subroutine eos_given_TPX

  subroutine eos_given_PSX(P, S, X, R, T, e, pt_index)

    ! In/out variables
    real(kind=dp_t), intent(  out) :: e, R, T
    real(kind=dp_t), intent(in   ) :: P, S, X(:)
    integer, optional, intent(in   ) :: pt_index(:)

    ! Local variables
    logical :: do_diag

    do_diag = .false.

    p_eos(1) = P
    s_eos(1) = S
    xn_eos(1,1:nspec) = X(1:nspec)

    call eos(eos_input_ps, den_eos, temp_eos, &
             npts, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             do_diag)

    R = den_eos(1)
    T = temp_eos(1)
    e = e_eos(1)

  end subroutine eos_given_PSX


  !---------------------------------------------------------------------------
  ! The main interface -- this is used directly by MAESTRO
  !---------------------------------------------------------------------------
  subroutine eos(input, dens, temp, &
                 npoints, &
                 xmass, &
                 pres, enthalpy, eint, &
                 c_v, c_p, ne, eta, pele, &
                 dPdT, dPdR, dEdT, dEdR, &
                 dPdX, dhdX, &
                 gam1, cs, entropy, &
                 dsdT, dsdR, &
                 do_eos_diag, &
                 pt_index)

    use bl_error_module
    use fundamental_constants_module, only: k_B, n_A, hbar

! dens     -- mass density (g/cc)
! temp     -- temperature (K)
! npoints     -- the number of elements in input/output arrays
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
! entropy  -- entropy (erg/g/K)  NOTE: presently the entropy expression is 
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

    implicit none

    logical do_eos_diag
    integer, intent(in) :: input
    integer, intent(in) :: npoints

    real(kind=dp_t) :: dens(npoints), temp(npoints)
    real(kind=dp_t) :: xmass(npoints,nspec)
    real(kind=dp_t) :: pres(npoints), enthalpy(npoints), eint(npoints)
    real(kind=dp_t) :: c_v(npoints), c_p(npoints)
    real(kind=dp_t) :: ne(npoints), eta(npoints), pele(npoints)
    real(kind=dp_t) :: dPdT(npoints), dPdR(npoints), dedT(npoints), dedR(npoints)
    real(kind=dp_t) :: gam1(npoints), entropy(npoints), cs(npoints)
    real(kind=dp_t) :: dPdX(npoints,nspec), dedX(npoints,nspec), dhdX(npoints,nspec)
    real(kind=dp_t) :: dsdT(npoints), dsdR(npoints)

    integer, optional, intent(in   ) :: pt_index(:)


    ! local variables
    real(kind=dp_t) :: ymass(npoints,nspec)    
    real(kind=dp_t) :: mu(npoints)
    real(kind=dp_t) :: dmudX, sum_y

    ! get the mass of a nucleon from Avogadro's number.
    real(kind=dp_t), parameter :: m_nucleon = 1.d0/n_A

    integer :: k, n

    ! general sanity checks
    if (.not. initialized) call bl_error('EOS: not initialized')
      
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
          
          do n = 1, nspec
             ymass(k,n) = xmass(k,n)/aion(n)
             sum_y = sum_y + ymass(k,n)
          enddo
          
          mu(k) = 1.d0/sum_y
       enddo

    else
       ! assume completely ionized species

       do k = 1, npoints
          sum_y  = 0.d0
          
          do n = 1, nspec
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
          temp(k) = (enthalpy(k)*mu(k)*m_nucleon/k_B)*(gamma_const - 1.0)/gamma_const
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
          temp(k) = eint(k)*mu(k)*m_nucleon*(gamma_const-1.0)/k_B

       enddo

    else if (input .EQ. eos_input_ps) then
       
       ! pressure and entropy are inputs
       do k = 1, npoints

          ! Solve for the temperature
          ! Invert Sackur-Tetrode eqn (below) using 
          ! rho = p mu m_nucleon / (k T)
          temp(k) = pres(k)**(2.0_dp_t/5.0_dp_t) * &
                    ( 2.0_dp_t*M_PI*hbar*hbar/(mu(k)*m_nucleon) )**(3.0_dp_t/5.0_dp_t) * &
                    dexp(2.0_dp_t*mu(k)*m_nucleon*entropy(k)/(5.0_dp_t*k_B) - 1.0_dp_t) / &
                    k_B

          ! Solve for the density
          ! rho = p mu m_nucleon / (k T)
          dens(k) = pres(k)*mu(k)*m_nucleon/(k_B*temp(k))

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
       eint(k) = pres(k)/(gamma_const - 1.0)/dens(k)

       ! enthalpy is h = e + p/rho
       enthalpy(k) = eint(k) + pres(k)/dens(k)

       ! entropy (per gram) of an ideal monoatomic gas (the Sactur-Tetrode equation)
       ! NOTE: this expression is only valid for gamma = 5/3.
       entropy(k) = (k_B/(mu(k)*m_nucleon))*(2.5_dp_t + &
            log( ( (mu(k)*m_nucleon)**2.5/dens(k) )*(k_B*temp(k))**1.5_dp_t / (2.0_dp_t*M_PI*hbar*hbar)**1.5_dp_t ) )

       ! compute the thermodynamic derivatives and specific heats 
       dPdT(k) = pres(k)/temp(k)
       dPdR(k) = pres(k)/dens(k)
       dedT(k) = eint(k)/temp(k)
       dedR(k) = 0.d0
       dsdT(k) = 0.d0
       dsdR(k) = 0.d0

       c_v(k) = dedT(k)
       c_p(k) = gamma_const*c_v(k)

       gam1(k) = gamma_const

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
       cs(k) = sqrt(gamma_const*pres(k)/dens(k))

    enddo

    return
  end subroutine eos

end module eos_module
