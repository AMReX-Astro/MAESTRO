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

  logical, save, private :: do_coulomb
  real(kind=dp_t), save, private :: smallt
  real(kind=dp_t), save, private :: smalld

  logical, save, private :: initialized = .false.

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

    ! gamma_in is a dummy variable -- it is needed in a generic interface
    ! for an EOS, but only used in a gamma-law EOS, not this general EOS
    real(kind=dp_t), intent(in), optional :: gamma_in
 
    if (present(use_eos_coulomb)) then
       do_coulomb = use_eos_coulomb
    else
       do_coulomb = .true.
    endif
 
    if (present(small_temp)) then
      if (small_temp > 0.d0) then
       smallt = small_temp
      else
       smallt = 5.d6
      end if
    else
       smallt = 5.d6
    endif
 
    if (present(small_dens)) then
       if (small_dens > 0.d0) then
         smalld = small_dens
       else
         smalld = 1.d-5
       end if
    else
       smalld = 1.d-5
    endif

    ! call the helmeos initialization routine and read in the table 
    ! containing the electron contribution.
    call helmeos_init()
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

  subroutine eos_given_ReX(G, P, C, T, dpdr, dpde, R, e, X)

     ! In/out variables
     real(kind=dp_t), intent(  out) :: G, P, C, dpdr, dpde
     real(kind=dp_t), intent(inout) :: T
     real(kind=dp_t), intent(in   ) :: R, e, X(:)

     ! Local variables
     logical :: do_diag

     do_diag = .false.

     temp_eos(1) = T
      den_eos(1) = R
        e_eos(1) = e
      xn_eos(1,1:nspec) = X(1:nspec)

     call eos(eos_input_re, den_eos, temp_eos, &
              npts, nspec, &
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
    dpdr = dpdr_eos(1) - dpdt_eos(1)*dedr_eos(1)/dedt_eos(1)
    dpde = dpdt_eos(1) / dedt_eos(1)

  end subroutine eos_given_ReX

  subroutine eos_e_given_RPX(e, T, R, P, X)

     ! In/out variables
     real(kind=dp_t), intent(  out) :: e
     real(kind=dp_t), intent(in   ) :: R, p, X(:)
     real(kind=dp_t), intent(inout) :: T

     ! Local variables
     logical :: do_diag

     do_diag = .false.

     temp_eos(1) = T
      den_eos(1) = R
        p_eos(1) = P
      xn_eos(1,1:nspec) = X(1:nspec)

      print *, 'calling EOS: ', temp_eos(1)

     call eos(eos_input_rp, den_eos, temp_eos, &
              npts, nspec, &
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

  subroutine eos_S_given_ReX(S, R, e, T, X)

     implicit none

     ! In/out variables
     real(kind=dp_t), intent(  out) :: S
     real(kind=dp_t), intent(in   ) :: R, e, T, X(:)

     ! Local variables
     logical :: do_diag

     do_diag = .false.

     temp_eos(1) = T
      den_eos(1) = R
        e_eos(1) = e
      xn_eos(1,1:nspec) = X(1:nspec)

     call eos(eos_input_re, den_eos, temp_eos, &
              npts, nspec, &
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

  subroutine eos_given_RTX(e, P, R, T, X)

     ! In/out variables
     real(kind=dp_t), intent(  out) :: e, P
     real(kind=dp_t), intent(in   ) :: R, T, X(:)

     ! Local variables
     logical :: do_diag

     do_diag = .false.

      den_eos(1) = R
     temp_eos(1) = T
      xn_eos(1,1:nspec) = X(1:nspec)

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

    P  =    p_eos(1)
    e  =    e_eos(1)

  end subroutine eos_given_RTX

  subroutine eos_given_TPX(e, P, R, T, X)

     ! In/out variables
     real(kind=dp_t), intent(  out) :: e, R
     real(kind=dp_t), intent(in   ) :: P, T, X(:)

     ! Local variables
     logical :: do_diag

     do_diag = .false.

     ! An initial guess of density needs to be given
     den_eos(1) = R 
     p_eos(1) = P
     temp_eos(1) = T
     xn_eos(1,1:nspec) = X(1:nspec)

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

    R  =    den_eos(1)
    e  =    e_eos(1)

  end subroutine eos_given_TPX


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

! a generic wrapper for the Helmholtz based electron/positron degenerate
! EOS.  
!
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
!   The EOS does not return the thermodynamic derivatives with respect
!   to the mass fractions, but rather, only due to the average atomic
!   mass (abar) and average proton number (zbar):
!
!     abar = ( sum_k {X_k} ) / ( sum_k {X_k/A_k} )
!
!     zbar = ( sum_k {Z_k X_k/ A_k} ) / ( sum_k {X_k/A_k} )
!
!   using the chain rule:
!
!   dp/dX_k = dp/d(abar) d(abar)/dX_k  +  dp/d(zbar) d(zbar)/dX_k
!
!   and using the above definitions of abar and zbar and sum_k {X_k} = 1
!
!   d(abar)/dX_k = abar * (1 - abar/A_k)
!   d(zbar)/dX_k = (Z_k - zbar) / ( A_k * sum_i {X_i/A_i} )
!
    implicit none

!    include 'vector_eos.dek'


!     ::::: Arguments
    logical             :: do_eos_diag
    integer, intent(in) :: input, npoints, nspecies

    ! some of these quantites can be inputs or outputs
    real(kind=dp_t), intent(inout) :: dens(npoints), temp(npoints)
    real(kind=dp_t), intent(in)    :: xmass(npoints,nspecies)
    real(kind=dp_t), intent(inout) :: pres(npoints), enthalpy(npoints), &
                                      eint(npoints)

    ! these quantities are always outputs
    real(kind=dp_t), intent(out) :: c_v(npoints), c_p(npoints)
    real(kind=dp_t), intent(out) :: ne(npoints), eta(npoints), pele(npoints)
    real(kind=dp_t), intent(out) :: dPdT(npoints), dPdR(npoints), &
                                    dedT(npoints), dedR(npoints)
    real(kind=dp_t), intent(out) :: gam1(npoints)
    real(kind=dp_t), intent(out) :: entropy(npoints), cs(npoints)
    real(kind=dp_t), intent(out) :: dPdX(npoints,nspecies), &
                                    dhdX(npoints,nspecies)
    real(kind=dp_t), intent(out) :: dsdT(npoints), dsdR(npoints)
    
!     ::::: Local variables and arrays

    integer :: i, k, n, iter, niter, max_newton
    parameter (max_newton = 100)
    
    real(kind=dp_t) :: error
    real(kind=dp_t) :: ymass(npoints,nspecies)
    real(kind=dp_t) :: abar(npoints), zbar(npoints)
    real(kind=dp_t) :: energy_want(npoints)
    real(kind=dp_t) :: enthalpy_want(npoints)
    real(kind=dp_t) :: pres_want(npoints)
    real(kind=dp_t) :: dhdt(npoints)
    real(kind=dp_t) :: tnew(npoints)
    real(kind=dp_t) :: dnew(npoints)
    real(kind=dp_t) :: enth1(npoints)
    real(kind=dp_t) :: ener1(npoints)
    real(kind=dp_t) :: dedX(npoints,nspecies)

    real(kind=dp_t) :: dpdd, pres1

    real(kind=dp_t) :: ttol
    parameter (ttol = 1.0d-8)
    real(kind=dp_t) :: dtol
    parameter (dtol = 1.0d-8)
    real(kind=dp_t) :: stol
    parameter (stol = 1.0d-8)

    logical eosfail

!     ::::: Input/Output arrays for call to helmeos
    real(kind=dp_t) :: temp_row(npoints), den_row(npoints), abar_row(npoints), &
                     zbar_row(npoints), etot_row(npoints), ptot_row(npoints), &
                     cv_row(npoints), cp_row(npoints), &
                     xne_row(npoints), xnp_row(npoints), etaele_row(npoints), &
                     pele_row(npoints), ppos_row(npoints), dpd_row(npoints), &
                     dpt_row(npoints), dpa_row(npoints), dpz_row(npoints), &
                     ded_row(npoints), det_row(npoints), dea_row(npoints), &
                     dez_row(npoints), &
                     stot_row(npoints), dsd_row(npoints), dst_row(npoints)
    real(kind=dp_t) :: gam1_row(npoints), cs_row(npoints)
      

! Dont allow multi-streaming of this function on CRAY X1
!DIR$  NOSTREAM
    if (.not. initialized) call bl_error('EOS: not initialized')
      
    if (nspecies /= nspec) then
       call bl_error('EOS: too many species')
    endif

    if (npoints > NP) then
       call bl_error('EOS: eos called with too large of a vector size')
    endif

!      print*,'EOS with npoints, nspecies = ',npoints,nspecies

! get the average atomic mass and the average proton number for the current
! composition
! MLW: use tnew and dnew as scratch arrays for azbar
!      call azbar(xmass,nxpts,aion,zion,npoints,nspecies,ymass,
!     +           abar,zbar,tnew,dnew)
!
!
! MLW: use tnew and dnew as scratch arrays 
!..mass fractions     = xmass
!..number of nucleons = aion
!..charge of nucleus  = zion
!..output:
!..molar abundances        = ymass
!..mean number of nucleons = abar
!..mean nucleon charge     = zbar
    do k = 1, npoints
       tnew(k)  = 0.0d0
       dnew(k)   = 0.0d0
    enddo

    do i=1,nspecies
       do k = 1, npoints
          ymass(k,i) = xmass(k,i)/aion(i)
          dnew(k)    = dnew(k) + ymass(k,i)
          tnew(k)    = tnew(k) + zion(i) * ymass(k,i)
       enddo
    enddo

    do k = 1, npoints
       abar(k) = 1.0d0/dnew(k)
       zbar(k) = tnew(k) * abar(k)
    enddo

    if (input .EQ. eos_input_rt) then

!---------------------------------------------------------------------------
! input = 1: dens, temp, and xmass are inputs
!---------------------------------------------------------------------------

! we are taking density, temperature, and composition as given
       do k = 1, npoints
          temp_row(k) = temp(k)
          den_row(k)  = dens(k)
          abar_row(k) = abar(k)
          zbar_row(k) = zbar(k)
       enddo
         
! call the eos
       call helmeos(do_coulomb,npoints,eosfail, &
                   temp_row, den_row, abar_row, zbar_row, &
                   etot_row, ptot_row, &
                   cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                   pele_row, ppos_row, &
                   dpd_row, dpt_row, dpa_row, dpz_row, &
                   ded_row, det_row, dea_row, dez_row, & 
                   gam1_row, cs_row, stot_row, &
                   dsd_row, dst_row)

       if (eosfail) call bl_error('EOS: error in the EOS')

! fill the outputs
       do k = 1, npoints
          pres(k) = ptot_row(k)
          eint(k) = etot_row(k)
         
          enthalpy(k) = eint(k) + pres(k)/dens(k)

          c_v(k) = cv_row(k)
          c_p(k) = cp_row(k)

! store the number density of electrons and positrons, the degeneracy
! parameter, and the total electron/positron pressure
          ne(k)   = xne_row(k) + xnp_row(k)
          eta(k)  = etaele_row(k)
          pele(k) = pele_row(k) + ppos_row(k)

          dPdR(k) = dpd_row(k)
          dPdT(k) = dpt_row(k)
          dEdR(k) = ded_row(k)
          dEdT(k) = det_row(k)
          gam1(k) = gam1_row(k)
!          cs(k) =   cs_row(k)
          cs(k) =   sqrt(gam1(k)*pres(k)/dens(k))
          entropy(k) = stot_row(k)
          dsdT(k) = dst_row(k)
          dsdR(k) = dsd_row(k)

          do n = 1, nspecies
             dpdX(k,n) = dpa_row(k) * (abar(k)/aion(n))* &
                                      (aion(n) - abar(k)) + &
                         dpz_row(k) * (abar(k)/aion(n))* &
                                      (zion(n) - zbar(k))

             dEdX(k,n) = dea_row(k) * (abar(k)/aion(n))* &
                                      (aion(n) - abar(k)) + &
                         dez_row(k) * (abar(k)/aion(n))* &
                                      (zion(n) - zbar(k))

             ! create the enthalpy derivatives wrt average composition --
             ! hold pressure constant!!!
             dhdX(k,n) = dEdX(k,n) + &
                  (pres(k)/dens(k)**2 - dEdR(k))*dpdX(k,n)/dPdr(k)
          enddo

       enddo

    else if (input .EQ. eos_input_rh) then

!---------------------------------------------------------------------------
! input = 2: dens, enthalpy, and xmass are inputs
!---------------------------------------------------------------------------

       ! load the initial guess
       do k = 1, npoints
          temp_row(k) = temp(k)
          den_row(k)  = dens(k)
          abar_row(k) = abar(k)
          zbar_row(k) = zbar(k)

          if (do_eos_diag) print*,'T/D INIT ',temp(k),dens(k)

          ! we want to converge to the given enthalpy
          enthalpy_want(k) = enthalpy(k)

          if (do_eos_diag) print*,'WANT H ',enthalpy(k)
       enddo

       call helmeos(do_coulomb,npoints, eosfail, &
                    temp_row, den_row, abar_row, zbar_row, &
                    etot_row, ptot_row, &
                    cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                    pele_row, ppos_row, &
                    dpd_row, dpt_row, dpa_row, dpz_row, &
                    ded_row, det_row, dea_row, dez_row, &
                    gam1_row, cs_row, stot_row, &
                    dsd_row, dst_row)

       if (eosfail) call bl_error('EOS: error in the EOS')

       do iter = 1, max_newton

          niter = iter

          ! recompute the enthalpy and it's temperature derivative
          do k = 1, npoints

             enth1(k) = etot_row(k) + ptot_row(k)/dens(k)
             if (do_eos_diag) print*,'ENTH1 ',iter,enth1(k)

             dhdt(k) = det_row(k) + dpt_row(k)/dens(k)
             if (do_eos_diag) print*,'DHDT ',iter,dhdt(k)

             tnew(k) = temp_row(k) - &
                  (enth1(k) - enthalpy_want(k))/dhdt(k)

             if (do_eos_diag) then
                print *, 'TNEW FIRST ', temp_row(k), ' - ', &
                     enth1(k) - enthalpy_want(k), ' / ', dhdt(k)
             endif

             ! don't let the temperature change by more than a factor of two
             tnew(k) = max(.5d0*temp_row(k), &
                           min(tnew(k), 2.d0*temp_row(k)))

             ! don't let us freeze
             tnew(k) = max(smallt, tnew(k))

             if (do_eos_diag) print*,'TNEW AFTER ',iter,tnew(1)
          enddo

          ! compute the error
          error = 0.0d0
          do k = 1, npoints
             error = max(error,abs(tnew(k) - temp_row(k))/temp_row(k))

             ! store the new temperature
             temp_row(k) = tnew(k)
          enddo

          if (error .LT. ttol) goto 70
        
          call helmeos(do_coulomb,npoints, eosfail, &
                       temp_row, den_row, abar_row, zbar_row, &
                       etot_row, ptot_row, &
                       cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                       pele_row, ppos_row, &
                       dpd_row, dpt_row, dpa_row, dpz_row, &
                       ded_row, det_row, dea_row, dez_row, &
                       gam1_row, cs_row, stot_row, &
                       dsd_row, dst_row)

          if (eosfail) call bl_error('EOS: error in the EOS')
        
       enddo

       ! Land here if too many iterations are needed

       continue

       call bl_error('EOS: Newton-Raphson failed:2: too many iterations')

70     continue

       ! store the end result
       do k = 1, npoints
          temp(k) = tnew(k)
          pres(k) = ptot_row(k)
          eint(k) = etot_row(k)
          
          c_v(k) = cv_row(k)
          c_p(k) = cp_row(k)

          ! store the number density of electrons and positrons, the degeneracy
          ! parameter, and the total electron/positron pressure
          ne(k)   = xne_row(k) + xnp_row(k)
          eta(k)  = etaele_row(k)
          pele(k) = pele_row(k) + ppos_row(k)
          
          dPdR(k) = dpd_row(k)
          dPdT(k) = dpt_row(k)
          dEdR(k) = ded_row(k)
          dEdT(k) = det_row(k)   ! c_v
          gam1(k) = gam1_row(k)
!          cs(k) =   cs_row(k)
          cs(k) =   sqrt(gam1(k)*pres(k)/dens(k))
          entropy(k) = stot_row(k)
          dsdT(k) = dst_row(k)
          dsdR(k) = dsd_row(k)

          do n = 1, nspecies
             dpdX(k,n) = dpa_row(k) * (abar(k)/aion(n))* &
                                      (aion(n) - abar(k)) + &
                         dpz_row(k) * (abar(k)/aion(n))* &
                                      (zion(n) - zbar(k))

             dEdX(k,n) = dea_row(k) * (abar(k)/aion(n))* &
                                      (aion(n) - abar(k)) + &
                         dez_row(k) * (abar(k)/aion(n))* &
                                      (zion(n) - zbar(k))

             ! create the enthalpy derivatives wrt average composition --
             ! hold pressure constant!!!
             dhdX(k,n) = dEdX(k,n) + &
                  (pres(k)/dens(k)**2 - dEdR(k))*dpdX(k,n)/dPdr(k)

          enddo

       enddo

    else if (input .EQ. eos_input_tp ) then

!---------------------------------------------------------------------------
! input = 3: temp, pres, and xmass are inputs
!---------------------------------------------------------------------------

       ! load the initial guess
       do k = 1, npoints
          temp_row(k) = temp(k)
          den_row(k)  = dens(k)
          abar_row(k) = abar(k)
          zbar_row(k) = zbar(k)

          ! we want to converge to the given pressure
          pres_want(k) = pres(k)

          if (pres_want(k) < ZERO) &
               call bl_error('ERROR: pressure < 0 in the EOS')
       enddo
         
       call helmeos(do_coulomb,npoints, eosfail, &
                    temp_row, den_row, abar_row, zbar_row, &
                    etot_row, ptot_row, &
                    cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                    pele_row, ppos_row, &
                    dpd_row, dpt_row, dpa_row, dpz_row, &
                    ded_row, det_row, dea_row, dez_row, &
                    gam1_row, cs_row, stot_row, &
                    dsd_row, dst_row)
       
       if (eosfail) call bl_error('EOS: error in the EOS')

       do iter = 1, max_newton

          niter = iter

          ! recompute the density and it's temperature derivative
          do k = 1, npoints
             pres1 = ptot_row(k)
             dpdd  = dpd_row(k)
             
             dnew(k) = den_row(k) - &
                  (pres1 - pres_want(k))/dpdd

             ! don't let the density change by more than an order of magnitude
             dnew(k) = max(.5d0*den_row(k), &
                           min(dnew(k), 2.d0*den_row(k)))

          enddo

          error = 0.0d0
          do k = 1, npoints

             ! compute the error
             error = max(error,abs(dnew(k) - den_row(k))/den_row(k))

             ! store the new density
             den_row(k) = dnew(k)
          enddo

          ! check if we are evacuating, if so, set the density to smalld, and adjust
          ! the error so we iterate on this one
          do k = 1, npoints
             if (den_row(k) .LT. smalld) then
                den_row(k) = smalld
                error = 1.1d0*dtol
             endif
          enddo

          if (error .LT. dtol) goto 170
        
          call helmeos(do_coulomb,npoints, eosfail, &
                       temp_row, den_row, abar_row, zbar_row, &
                       etot_row, ptot_row, &
                       cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                       pele_row, ppos_row, &
                       dpd_row, dpt_row, dpa_row, dpz_row, &
                       ded_row, det_row, dea_row, dez_row, &
                       gam1_row, cs_row, stot_row, &
                       dsd_row, dst_row)
        
          if (eosfail) call bl_error('EOS: error in the EOS')

       end do

       ! Land here if too many iterations are needed

       call bl_error('EOS: Newton-Raphson failed:3: too many iterations')         

170    continue

       ! store the end result
       do k = 1, npoints
          dens(k) = dnew(k)
          temp(k) = temp_row(k)
          eint(k) = etot_row(k)
          enthalpy(k) = eint(k) + ptot_row(k)/dens(k)

          c_v(k) = cv_row(k)
          c_p(k) = cp_row(k)

          ! store the number density of electrons and positrons, the degeneracy
          ! parameter, and the total electron/positron pressure
          ne(k)   = xne_row(k) + xnp_row(k)
          eta(k)  = etaele_row(k)
          pele(k) = pele_row(k) + ppos_row(k)
          
          dPdR(k) = dpd_row(k)
          dPdT(k) = dpt_row(k)
          dEdR(k) = ded_row(k)
          dEdT(k) = det_row(k)   ! c_v
          gam1(k) = gam1_row(k)
!          cs(k) =   cs_row(k)
          cs(k) =   sqrt(gam1(k)*pres(k)/dens(k))
          entropy(k) = stot_row(k)
          dsdT(k) = dst_row(k)
          dsdR(k) = dsd_row(k)

          do n = 1, nspecies
             dpdX(k,n) = dpa_row(k) * (abar(k)/aion(n))* &
                                      (aion(n) - abar(k)) + &
                         dpz_row(k) * (abar(k)/aion(n))* &
                                      (zion(n) - zbar(k))

             dEdX(k,n) = dea_row(k) * (abar(k)/aion(n))* &
                                      (aion(n) - abar(k)) + &
                         dez_row(k) * (abar(k)/aion(n))* &
                                      (zion(n) - zbar(k))

             ! create the enthalpy derivatives wrt average composition --
             ! hold pressure constant!!!
             dhdX(k,n) = dEdX(k,n) + &
                  (pres(k)/dens(k)**2 - dEdR(k))*dpdX(k,n)/dPdr(k)

          enddo

       enddo

    else if (input .EQ. eos_input_rp ) then

!---------------------------------------------------------------------------
! input = 4: dens, pres, and xmass are inputs
!---------------------------------------------------------------------------

       ! Load the initial guess
       do k = 1, npoints
          temp_row(k) = temp(k)
          den_row(k)  = dens(k)
          abar_row(k) = abar(k)
          zbar_row(k) = zbar(k)
          if (do_eos_diag) print*,'T/D INIT ',temp(k),dens(k)
       enddo

       ! We want to converge to the given pressure
       do k = 1, npoints
          pres_want(k) = pres(k)
          if (do_eos_diag) print*,'P WANT ',pres(k)

          if (pres_want(k) < ZERO) &
               call bl_error('ERROR: pressure < 0 in the EOS')
       enddo

       ! First pass
       call helmeos(do_coulomb,npoints, eosfail, &
                    temp_row, den_row, abar_row, zbar_row, &
                    etot_row, ptot_row, &
                    cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                    pele_row, ppos_row, &
                    dpd_row, dpt_row, dpa_row, dpz_row, &
                    ded_row, det_row, dea_row, dez_row, &
                    gam1_row, cs_row, stot_row, &
                    dsd_row, dst_row)

       if (eosfail) call bl_error('EOS: error in the EOS')
       
       do iter = 1, max_newton

          niter = iter

          ! recompute the temperature but dont allow it to change too much
          do k = 1, npoints
             tnew(k) = temp_row(k) - &
                  (ptot_row(k) - pres_want(k))/dpt_row(k)

             if (do_eos_diag) print*,'PRES ',ptot_row(k),pres_want(k)
             if (do_eos_diag) print*,'PRES DIFF ',ptot_row(k)-pres_want(k)
             
             if (do_eos_diag) print*,'DPDT FAC ', 1.0/dpt_row(k)

             if (do_eos_diag) print*,'TNEW BEFORE MAX ',iter,tnew(k)

             ! don't let the temperature change by more than a factor of 2
             tnew(k) = max(.5d0*temp_row(k), &
                           min(tnew(k), 2.d0*temp_row(k)))

             ! don't let us freeze
             tnew(k) = max(smallt, tnew(k))

             if (do_eos_diag) print*,'TNEW AFTER MAX ',iter,tnew(k)
             if (do_eos_diag) print*,' '
          enddo

          error = 0.0d0

          ! compute the error and store the new temperature
          do k = 1, npoints
             error = max(error,abs(tnew(k) - temp_row(k))/temp_row(k))
             if (do_eos_diag) print *,'ERROR  ',iter,error
             if (do_eos_diag) print*,' '
             temp_row(k) = tnew(k)
          enddo

          if (error .LT. ttol) goto 870
        
          call helmeos(do_coulomb,npoints, eosfail, &
                       temp_row, den_row, abar_row, zbar_row, &
                       etot_row, ptot_row, &
                       cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                       pele_row, ppos_row, &
                       dpd_row, dpt_row, dpa_row, dpz_row, &
                       ded_row, det_row, dea_row, dez_row, &
                       gam1_row, cs_row, stot_row, &
                       dsd_row, dst_row)

          if (eosfail) call bl_error('EOS: error in the EOS')
             
        
       end do

       ! Land here if too many iterations are needed

       print *, 'helmeos input==4 failed to converge, iter = ',niter
       do k = 1, npoints
          print *, 'k, error, temp_row(k), den_row(k) = ', &
                k,error,temp_row(k),den_row(k)
       enddo

       call bl_error('EOS: Newton-Raphson failed:4: too many iterations')

870    continue

       ! store the end result
       do k = 1, npoints

          ! jbb
          ! temp(k) = tnew(k)
          temp(k) = temp_row(k)
          dens(k) = den_row(k)
          eint(k) = etot_row(k)
          c_v(k) = cv_row(k)
          c_p(k) = cp_row(k)

          enthalpy(k) = eint(k) + ptot_row(k)/dens(k)

          ! store the number density of electrons and positrons, the degeneracy
          ! parameter, and the total electron/positron pressure
          ne(k)   = xne_row(k) + xnp_row(k)
          eta(k)  = etaele_row(k)
          pele(k) = pele_row(k) + ppos_row(k)
          
          dPdR(k) = dpd_row(k)
          dPdT(k) = dpt_row(k)
          dEdR(k) = ded_row(k)
          dEdT(k) = det_row(k)   ! c_v
          gam1(k) = gam1_row(k)
!          cs(k) =   cs_row(k)
          cs(k) =   sqrt(gam1(k)*pres(k)/dens(k))
          entropy(k) = stot_row(k)
          dsdT(k) = dst_row(k)
          dsdR(k) = dsd_row(k)

          do n = 1, nspecies
             dpdX(k,n) = dpa_row(k) * (abar(k)/aion(n))* &
                                      (aion(n) - abar(k)) + &
                         dpz_row(k) * (abar(k)/aion(n))* &
                                      (zion(n) - zbar(k))

             dEdX(k,n) = dea_row(k) * (abar(k)/aion(n))* &
                                      (aion(n) - abar(k)) + &
                         dez_row(k) * (abar(k)/aion(n))* &
                                      (zion(n) - zbar(k))

             ! create the enthalpy derivatives wrt average composition --
             ! hold pressure constant!!!
             dhdX(k,n) = dEdX(k,n) + &
                  (pres(k)/dens(k)**2 - dEdR(k))*dpdX(k,n)/dPdr(k)

          enddo

       enddo

    else if (input .EQ. eos_input_re) then

!---------------------------------------------------------------------------
! input = 5: dens, energy, and xmass are inputs
!---------------------------------------------------------------------------

!      do_eos_diag = .true.
       ! load the initial guess
       do k = 1, npoints
          temp_row(k) = temp(k)
          den_row(k)  = dens(k)
          abar_row(k) = abar(k)
          zbar_row(k) = zbar(k)

          if (do_eos_diag) print*,'T/D INIT ',temp(k),dens(k)

          ! we want to converge to the given enthalpy
          energy_want(k) = eint(k)

          if (energy_want(k) < ZERO) &
               call bl_error('ERROR: energy < 0 in the EOS')

          if (do_eos_diag) print*,'WANT e ',energy_want(k)
       enddo

       call helmeos(do_coulomb,npoints, eosfail, &
                    temp_row, den_row, abar_row, zbar_row, &
                    etot_row, ptot_row, &
                    cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                    pele_row, ppos_row, &
                    dpd_row, dpt_row, dpa_row, dpz_row, &
                    ded_row, det_row, dea_row, dez_row, &
                    gam1_row, cs_row, stot_row, &
                    dsd_row, dst_row)

       if (eosfail) call bl_error('EOS: error in the EOS')

       do iter = 1, max_newton

          niter = iter

          ! recompute the energy and its temperature derivative
          do k = 1, npoints

             ener1(k) = etot_row(k)

             if (do_eos_diag) then
                print*,'ENER1 ',iter,ener1(k)
                print*,'DEDT ',iter,det_row(k)
             end if

             tnew(k) = temp_row(k) - (ener1(k)-energy_want(k))/det_row(k)

             if (do_eos_diag) then
                print *, 'TNEW FIRST ', temp_row(k), ' - ', &
                     ener1(k) - energy_want(k), ' / ', det_row(k)
             endif

             ! don't let the temperature change by more than a factor of two
             tnew(k) = max(.5d0*temp_row(k), &
                           min(tnew(k), 2.d0*temp_row(k)))

             ! don't let us freeze
             tnew(k) = max(smallt, tnew(k))

             if (do_eos_diag) print*,'TNEW AFTER ',iter,tnew(1)
          enddo

          ! compute the error
          error = 0.0d0
          do k = 1, npoints
             error = max(error,abs(tnew(k) - temp_row(k))/temp_row(k))

             ! store the new temperature
             temp_row(k) = tnew(k)
          enddo

          if (error .LT. ttol) goto 270
        
          call helmeos(do_coulomb,npoints, eosfail, &
                       temp_row, den_row, abar_row, zbar_row, &
                       etot_row, ptot_row, &
                       cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                       pele_row, ppos_row, &
                       dpd_row, dpt_row, dpa_row, dpz_row, &
                       ded_row, det_row, dea_row, dez_row, &
                       gam1_row, cs_row, stot_row, &
                       dsd_row, dst_row)

          if (eosfail) call bl_error('EOS: error in the EOS')
        
       enddo

       ! Land here if too many iterations are needed

       continue

       call bl_error('EOS: Newton-Raphson failed:2: too many iterations')

270     continue

       ! store the end result
       do k = 1, npoints
          temp(k) = tnew(k)
          pres(k) = ptot_row(k)
          
          c_v(k) = cv_row(k)
          c_p(k) = cp_row(k)

          ! store the number density of electrons and positrons, the degeneracy
          ! parameter, and the total electron/positron pressure
          ne(k)   = xne_row(k) + xnp_row(k)
          eta(k)  = etaele_row(k)
          pele(k) = pele_row(k) + ppos_row(k)
          
          dPdR(k) = dpd_row(k)
          dPdT(k) = dpt_row(k)
          dEdR(k) = ded_row(k)
          dEdT(k) = det_row(k)   ! c_v
          gam1(k) = gam1_row(k)
!          cs(k) =   cs_row(k)
          cs(k) =   sqrt(gam1(k)*pres(k)/dens(k))
          entropy(k) = stot_row(k)
          dsdT(k) = dst_row(k)
          dsdR(k) = dsd_row(k)

          do n = 1, nspecies
             dpdX(k,n) = dpa_row(k) * (abar(k)/aion(n))* &
                                      (aion(n) - abar(k)) + &
                         dpz_row(k) * (abar(k)/aion(n))* &
                                      (zion(n) - zbar(k))

             dEdX(k,n) = dea_row(k) * (abar(k)/aion(n))* &
                                      (aion(n) - abar(k)) + &
                         dez_row(k) * (abar(k)/aion(n))* &
                                      (zion(n) - zbar(k))

             ! create the enthalpy derivatives wrt average composition --
             ! hold pressure constant!!!
             dhdX(k,n) = dEdX(k,n) + &
                  (pres(k)/dens(k)**2 - dEdR(k))*dpdX(k,n)/dPdr(k)

          enddo

       enddo

    else 

       call bl_error('EOS: invalid input')

    endif


    return
  end subroutine eos

end module eos_module
