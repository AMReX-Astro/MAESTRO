!
! This eos is hacked to use a time-INdependent Abar and Zbar
!  both of which are functions of radius
!

module eos_module

  use bl_space, only: MAX_SPACEDIM
  use bl_types
  use bl_constants_module
  use network, only: nspec, aion, zion

  implicit none

  private 

  integer, parameter :: NP = 1
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

  common /eos_common/ xn_eos,temp_eos,den_eos,abar_eos,zbar_eos,e_eos,p_eos,h_eos
  common /eos_common/ cv_eos,cp_eos,xne_eos,eta_eos,pele_eos,dpdt_eos,dpdr_eos,dedr_eos
  common /eos_common/ dedt_eos,gam1_eos,cs_eos,s_eos,dsdt_eos,dsdr_eos,dpdX_eos,dhdX_eos
  common /eos_common/ conduct_eos,pt_index_eos
  SAVE /eos_common/
!$omp threadprivate(/eos_common/)

  integer, parameter, public :: eos_input_rt = 1  ! rho, T are inputs
  integer, parameter, public :: eos_input_rh = 2  ! rho, h are inputs
  integer, parameter, public :: eos_input_tp = 3  ! T, p are inputs
  integer, parameter, public :: eos_input_rp = 4  ! rho, p are inputs
  integer, parameter, public :: eos_input_re = 5  ! rho, e are inputs
  integer, parameter, public :: eos_input_ps = 6  ! p, s are inputs
 

  logical, save, private :: do_coulomb
  real(kind=dp_t), save, private :: smallt
  real(kind=dp_t), save, private :: smalld

  logical, save, private :: initialized = .false.

  private nspec, aion, zion

  public eos_init, eos_get_small_temp, eos_get_small_dens, eos


contains

  ! EOS initialization routine -- this is used by both MAESTRO and Castro
  ! For this general EOS, this calls helmeos_init() which reads in the 
  ! table with the electron component's properties.
  subroutine eos_init(small_temp, small_dens, gamma_in)

    use extern_probin_module, only: use_eos_coulomb

    implicit none
 
    real(kind=dp_t), intent(in), optional :: small_temp
    real(kind=dp_t), intent(in), optional :: small_dens

    ! gamma_in is a dummy variable -- it is needed in a generic interface
    ! for an EOS, but only used in a gamma-law EOS, not this general EOS
    real(kind=dp_t), intent(in), optional :: gamma_in
 
    do_coulomb = use_eos_coulomb
 
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


  !---------------------------------------------------------------------------
  ! The main interface -- this is used directly by MAESTRO
  !---------------------------------------------------------------------------
  subroutine eos(radius, input, dens, temp, &
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

! a generic wrapper for the Helmholtz based electron/positron degenerate
! EOS.  
!
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
! entropy  -- entropy (erg/g/K)
!
! input = 1 means dens, temp    , and xmass are inputs, return enthalpy, eint
!       = 2 means dens, enthalpy, and xmass are inputs, return temp    , eint
!                (note, temp should be filled with an initial guess)
!       = 3 means temp, pres    , and xmass are inputs, return dens    , etc
!       = 4 means dens, pres    , and xmass are inputs, return temp    , etc
!       = 5 means dens, eint    , and xmass are inputs, return temp    , etc
!       = 6 means pres, entropy , and xmass are inputs, return dens    , etc
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
    integer, intent(in) :: input, npoints

    ! some of these quantites can be inputs or outputs
    real(kind=dp_t), intent(in)    :: radius
    real(kind=dp_t), intent(inout) :: dens(npoints), temp(npoints)
    real(kind=dp_t), intent(in)    :: xmass(npoints,nspec)
    real(kind=dp_t), intent(inout) :: pres(npoints), enthalpy(npoints), &
                                      eint(npoints), entropy(npoints)

    ! these quantities are always outputs
    real(kind=dp_t), intent(out) :: c_v(npoints), c_p(npoints)
    real(kind=dp_t), intent(out) :: ne(npoints), eta(npoints), pele(npoints)
    real(kind=dp_t), intent(out) :: dPdT(npoints), dPdR(npoints), &
                                    dedT(npoints), dedR(npoints)
    real(kind=dp_t), intent(out) :: gam1(npoints)
    real(kind=dp_t), intent(out) :: cs(npoints)
    real(kind=dp_t), intent(out) :: dPdX(npoints,nspec), &
                                    dhdX(npoints,nspec)
    real(kind=dp_t), intent(out) :: dsdT(npoints), dsdR(npoints)

    integer, optional, intent(in   ) :: pt_index(:)

    
!     ::::: Local variables and arrays

    integer :: i, k, n, iter, niter, max_newton
    parameter (max_newton = 100)
    
    real(kind=dp_t) :: error, error2
    real(kind=dp_t) :: ymass(npoints,nspec)
    real(kind=dp_t) :: abar(npoints), zbar(npoints)
    real(kind=dp_t) :: energy_want(npoints)
    real(kind=dp_t) :: enthalpy_want(npoints)
    real(kind=dp_t) :: pres_want(npoints)
    real(kind=dp_t) :: entropy_want(npoints)
    real(kind=dp_t) :: dhdt(npoints)
    real(kind=dp_t) :: tnew(npoints)
    real(kind=dp_t) :: dnew(npoints)
    real(kind=dp_t) :: enth1(npoints)
    real(kind=dp_t) :: ener1(npoints)
    real(kind=dp_t) :: dedX(npoints,nspec)

    real(kind=dp_t) :: dpdd, pres1, entr1
    real(kind=dp_t) :: f, g, dfdd, dfdt, dgdd, dgdt, deld

    real(kind=dp_t) :: ttol
    parameter (ttol = 1.0d-8)
    real(kind=dp_t) :: dtol
    parameter (dtol = 1.0d-8)
    real(kind=dp_t) :: stol
    parameter (stol = 1.0d-8)

    logical eosfail
    integer dim_ptindex

    ! err_string is used to convert the pt_index information into a string
    character (len=64) :: err_string  

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


    if (.not. initialized) call bl_error('EOS: not initialized')
      
    if (npoints > NP) then
       call bl_error('EOS: eos called with too large of a vector size')
    endif

    ! this format statement is for writing into err_string -- make sure that
    ! the len of err_string can accomodate this format specifier
1001 format(1x,"zone index info: i = ", i5)
1002 format(1x,"zone index info: i = ", i5, '  j = ', i5)
1003 format(1x,"zone index info: i = ", i5, '  j = ', i5, '  k = ', i5)


!      print*,'EOS with npoints, nspec = ',npoints,nspec

! get the average atomic mass and the average proton number for the current
! composition
! MLW: use tnew and dnew as scratch arrays for azbar
!      call azbar(xmass,nxpts,aion,zion,npoints,nspec,ymass,
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

! npoints = 1 as set at top of module
    call interpolate(radius, abar(1), zbar(1))

!    write(*,*)abar(1), zbar(1)

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

       if (eosfail) then
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('EOS: error in the EOS', err_string)
          else
             call bl_error('EOS: error in the EOS')
          endif
       endif


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

          do n = 1, nspec
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

       if (eosfail) then
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('EOS: error in the EOS', err_string)
          else
             call bl_error('EOS: error in the EOS')
          endif
       endif

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

          if (eosfail) then
             if (present(pt_index)) then
                if (dim_ptindex .eq. 1) then 
                   write (err_string,1001) pt_index(1)
                else if (dim_ptindex .eq. 2) then 
                   write (err_string,1002) pt_index(1), pt_index(2)
                else if (dim_ptindex .eq. 3) then 
                   write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
                end if
                call bl_error('EOS: error in the EOS', err_string)
             else
                call bl_error('EOS: error in the EOS')
             endif
          endif
        
       enddo

       ! Land here if too many iterations are needed

       continue

       if (present(pt_index)) then
          if (dim_ptindex .eq. 1) then 
             write (err_string,1001) pt_index(1)
          else if (dim_ptindex .eq. 2) then 
             write (err_string,1002) pt_index(1), pt_index(2)
          else if (dim_ptindex .eq. 3) then 
             write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
          end if
          call bl_error('EOS: Newton-Raphson failed:2: too many iterations', err_string)
       else
          call bl_error('EOS: Newton-Raphson failed:2: too many iterations')
       endif

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

          do n = 1, nspec
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

          if (pres_want(k) < ZERO) then
             if (present(pt_index)) then
                if (dim_ptindex .eq. 1) then 
                   write (err_string,1001) pt_index(1)
                else if (dim_ptindex .eq. 2) then 
                   write (err_string,1002) pt_index(1), pt_index(2)
                else if (dim_ptindex .eq. 3) then 
                   write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
                end if
                call bl_error('ERROR: pressure < 0 in the EOS', err_string)
             else
                call bl_error('ERROR: pressure < 0 in the EOS')
             endif
          endif
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
       
       if (eosfail) then
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('EOS: error in the EOS', err_string)
          else
             call bl_error('EOS: error in the EOS')
          endif
       endif

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
        
          if (eosfail) then
             if (present(pt_index)) then
                if (dim_ptindex .eq. 1) then 
                   write (err_string,1001) pt_index(1)
                else if (dim_ptindex .eq. 2) then 
                   write (err_string,1002) pt_index(1), pt_index(2)
                else if (dim_ptindex .eq. 3) then 
                   write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
                end if
                call bl_error('EOS: error in the EOS', err_string)
             else
                call bl_error('EOS: error in the EOS')
             endif
          endif

       end do

       ! Land here if too many iterations are needed

       if (present(pt_index)) then
          if (dim_ptindex .eq. 1) then 
             write (err_string,1001) pt_index(1)
          else if (dim_ptindex .eq. 2) then 
             write (err_string,1002) pt_index(1), pt_index(2)
          else if (dim_ptindex .eq. 3) then 
             write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
          end if
          call bl_error('EOS: Newton-Raphson failed:3: too many iterations', err_string)         
       else
          call bl_error('EOS: Newton-Raphson failed:3: too many iterations')         
       endif

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

          do n = 1, nspec
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

          if (pres_want(k) < ZERO) then
             if (present(pt_index)) then
                if (dim_ptindex .eq. 1) then 
                   write (err_string,1001) pt_index(1)
                else if (dim_ptindex .eq. 2) then 
                   write (err_string,1002) pt_index(1), pt_index(2)
                else if (dim_ptindex .eq. 3) then 
                   write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
                end if
                call bl_error('ERROR: pressure < 0 in the EOS', err_string)
             else
                call bl_error('ERROR: pressure < 0 in the EOS')
             endif
          endif
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

       if (eosfail) then
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('EOS: error in the EOS', err_string)
          else
             call bl_error('EOS: error in the EOS')
          endif
       endif
       
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

          if (eosfail) then
             if (present(pt_index)) then
                if (dim_ptindex .eq. 1) then 
                   write (err_string,1001) pt_index(1)
                else if (dim_ptindex .eq. 2) then 
                   write (err_string,1002) pt_index(1), pt_index(2)
                else if (dim_ptindex .eq. 3) then 
                   write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
                end if
                call bl_error('EOS: error in the EOS', err_string)
             else
                call bl_error('EOS: error in the EOS')
             endif
          endif
             
        
       end do

       ! Land here if too many iterations are needed

       print *, 'helmeos input==4 failed to converge, iter = ',niter
       do k = 1, npoints
          print *, 'k, error, temp_row(k), den_row(k) = ', &
                k,error,temp_row(k),den_row(k)
       enddo

       if (present(pt_index)) then
          if (dim_ptindex .eq. 1) then 
             write (err_string,1001) pt_index(1)
          else if (dim_ptindex .eq. 2) then 
             write (err_string,1002) pt_index(1), pt_index(2)
          else if (dim_ptindex .eq. 3) then 
             write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
          end if
          call bl_error('EOS: Newton-Raphson failed:4: too many iterations', err_string)          
       else
          call bl_error('EOS: Newton-Raphson failed:4: too many iterations')          
       endif

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
          cs(k) =   cs_row(k)
!          cs(k) =   sqrt(gam1(k)*pres(k)/dens(k))
          entropy(k) = stot_row(k)
          dsdT(k) = dst_row(k)
          dsdR(k) = dsd_row(k)

          do n = 1, nspec
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

          if (energy_want(k) < ZERO) then
             print *,'BAD HERE ',pt_index(1)
             if (present(pt_index)) then
                if (dim_ptindex .eq. 1) then 
                   write (err_string,1001) pt_index(1)
                else if (dim_ptindex .eq. 2) then 
                   write (err_string,1002) pt_index(1), pt_index(2)
                else if (dim_ptindex .eq. 3) then 
                   write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
                end if
                call bl_error('ERROR: energy < 0 in the EOS', err_string)
             else
                call bl_error('ERROR: energy < 0 in the EOS')
             endif
          endif

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

       if (eosfail) then
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('EOS: error in the EOS', err_string)
          else
             call bl_error('EOS: error in the EOS')
          endif
       endif

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

          if (eosfail) then
             if (present(pt_index)) then
                if (dim_ptindex .eq. 1) then 
                   write (err_string,1001) pt_index(1)
                else if (dim_ptindex .eq. 2) then 
                   write (err_string,1002) pt_index(1), pt_index(2)
                else if (dim_ptindex .eq. 3) then 
                   write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
                end if
                call bl_error('EOS: error in the EOS', err_string)
             else
                call bl_error('EOS: error in the EOS')
             endif
          endif
        
       enddo

       ! Land here if too many iterations are needed

       continue

       if (present(pt_index)) then
          if (dim_ptindex .eq. 1) then 
             write (err_string,1001) pt_index(1)
          else if (dim_ptindex .eq. 2) then 
             write (err_string,1002) pt_index(1), pt_index(2)
          else if (dim_ptindex .eq. 3) then 
             write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
          end if
          call bl_error('EOS: Newton-Raphson failed:2: too many iterations', err_string)
       else
          call bl_error('EOS: Newton-Raphson failed:2: too many iterations')
       endif

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

          do n = 1, nspec
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

    else if (input .EQ. eos_input_ps) then
!---------------------------------------------------------------------------
! input = 6: pres, entropy, and xmass are inputs
!---------------------------------------------------------------------------

       ! load the initial guess
       do k = 1, npoints
          temp_row(k) = temp(k)
          den_row(k)  = dens(k)
          abar_row(k) = abar(k)
          zbar_row(k) = zbar(k)

          if (do_eos_diag) print*,'T/D INIT ',temp(k),dens(k)

          ! we want to converge to the given entropy and pressure
          entropy_want(k) = entropy(k)
          pres_want(k)    = pres(k)

          if (entropy_want(k) < ZERO) then
             if (present(pt_index)) then
                if (dim_ptindex .eq. 1) then 
                   write (err_string,1001) pt_index(1)
                else if (dim_ptindex .eq. 2) then 
                   write (err_string,1002) pt_index(1), pt_index(2)
                else if (dim_ptindex .eq. 3) then 
                   write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
                end if
                call bl_error('ERROR: entropy < 0 in the EOS', err_string)
             else
                call bl_error('ERROR: entropy < 0 in the EOS')
             endif
          endif

          if (pres_want(k) < ZERO) then
             if (present(pt_index)) then
                if (dim_ptindex .eq. 1) then 
                   write (err_string,1001) pt_index(1)
                else if (dim_ptindex .eq. 2) then 
                   write (err_string,1002) pt_index(1), pt_index(2)
                else if (dim_ptindex .eq. 3) then 
                   write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
                end if
                call bl_error('ERROR: pressure < 0 in the EOS', err_string)
             else
                call bl_error('ERROR: pressure < 0 in the EOS')
             endif
          endif

          if (do_eos_diag) then
             print*,'WANT s ',entropy_want(k)
             print*,'WANT pres', pres_want(k)
          endif
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

       if (eosfail) then
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('EOS: error in the EOS', err_string)
          else
             call bl_error('EOS: error in the EOS')
          endif
       endif

       do iter = 1, max_newton

          niter = iter

          ! correct density and temperature
          do k = 1, npoints

             pres1 = ptot_row(k)
             entr1 = stot_row(k)

             if (do_eos_diag) then
                print*,'PRES1 ',iter,pres1
                print*,'ENTR1 ',iter,entr1
             end if

             ! two functions, f and g, to iterate over
             f = pres_want(k) - pres1
             dfdd = -dpd_row(k)
             dfdt = -dpt_row(k)
             
             g = entropy_want(k) - entr1
             dgdd = -dsd_row(k)
             dgdt = -dst_row(k)
             !
             ! 0 = f + dfdd * deld + dfdt * delt
             ! 0 = g + dgdd * deld + dgdt * delt
             !
             deld = (f*dgdt - g*dfdt) / (dgdd*dfdt - dgdt*dfdd)

             dnew(k) = den_row(k) + deld

             tnew(k) = temp_row(k) - (f + dfdd*deld) / dfdt

             if (do_eos_diag) then
                print *, 'DNEW FIRST ', den_row(k), ' + ', &
                     f*dgdt - g*dfdt, ' / ', dgdd*dfdt - dgdt*dfdd
                print *, 'TNEW FIRST ', temp_row(k), ' - ', &
                     f + dfdd*deld, ' / ', dfdt
             endif

             ! don't let the temperature or density change by more
             ! than a factor of two
             tnew(k) = max(HALF*temp_row(k), &
                           min(tnew(k), TWO*temp_row(k)))
             dnew(k) = max(HALF*den_row(k), &
                           min(dnew(k), TWO*den_row(k)))

             ! don't let us freeze or evacuate
             tnew(k) = max(smallt, tnew(k))
             dnew(k) = max(smalld, dnew(k))

             if (do_eos_diag) then
                print*,'DNEW AFTER ',iter,dnew(1)
                print*,'TNEW AFTER ',iter,tnew(1)
             endif
          enddo

          ! compute the errors
          error = ZERO
          error2 = ZERO
          do k = 1, npoints
             error  = max(error ,abs(dnew(k) - den_row(k))/den_row(k))
             error2 = max(error2,abs(tnew(k) - temp_row(k))/temp_row(k))

             ! store the new temperature and density
             den_row(k) = dnew(k)
             temp_row(k) = tnew(k)
          enddo

          if (error .LT. dtol .and. error2 .LT. ttol) goto 370
        
          call helmeos(do_coulomb,npoints, eosfail, &
                       temp_row, den_row, abar_row, zbar_row, &
                       etot_row, ptot_row, &
                       cv_row, cp_row, xne_row, xnp_row, etaele_row, &
                       pele_row, ppos_row, &
                       dpd_row, dpt_row, dpa_row, dpz_row, &
                       ded_row, det_row, dea_row, dez_row, &
                       gam1_row, cs_row, stot_row, &
                       dsd_row, dst_row)

          if (eosfail) then
             if (present(pt_index)) then
                if (dim_ptindex .eq. 1) then 
                   write (err_string,1001) pt_index(1)
                else if (dim_ptindex .eq. 2) then 
                   write (err_string,1002) pt_index(1), pt_index(2)
                else if (dim_ptindex .eq. 3) then 
                   write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
                end if
                call bl_error('EOS: error in the EOS', err_string)
             else
                call bl_error('EOS: error in the EOS')
             endif
          endif
        
       enddo

       ! Land here if too many iterations are needed

       continue

       if (present(pt_index)) then
          if (dim_ptindex .eq. 1) then 
             write (err_string,1001) pt_index(1)
          else if (dim_ptindex .eq. 2) then 
             write (err_string,1002) pt_index(1), pt_index(2)
          else if (dim_ptindex .eq. 3) then 
             write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
          end if
          call bl_error('EOS: Newton-Raphson failed:2: too many iterations', err_string)
       else
          call bl_error('EOS: Newton-Raphson failed:2: too many iterations')
       endif

370     continue

       ! store the end result
       do k = 1, npoints
          dens(k) = dnew(k)
          temp(k) = tnew(k)
          
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
          dsdT(k) = dst_row(k)
          dsdR(k) = dsd_row(k)

          do n = 1, nspec
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

       if (present(pt_index)) then
          if (dim_ptindex .eq. 1) then 
             write (err_string,1001) pt_index(1)
          else if (dim_ptindex .eq. 2) then 
             write (err_string,1002) pt_index(1), pt_index(2)
          else if (dim_ptindex .eq. 3) then 
             write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
          end if
          call bl_error('EOS: invalid input', err_string)
       else
          call bl_error('EOS: invalid input')
       endif

    endif


    return

    contains
      
      subroutine interpolate(r, a, z)
        use model_parser_module, only: meanA, meanZ
        use geometry, only: dr_fine, nr_fine, r_cc_loc, spherical
        use probin_module, only: max_levs

        real(kind=dp_t), intent(in   ) :: r
        real(kind=dp_t), intent(inout) :: a, z
        
        !local
        integer :: index,n
        real(kind=dp_t) :: x0,x1,x2, y0,y1,y2 !, rfac
   
        index  = int(r / dr_fine)

!        rfac = (radius - dble(index)*dr_fine) / dr_fine

        ! FIXME: need to think about this more in the event of restart
        if (spherical .eq. 1) then; n = 1
           ! this really ought to be nlevs
        else; n = max_levs
        endif

        ! piecewise linear
!         if (r .ge. r_cc_loc(n,index)) then
!            if (index .ge. nr_fine-1) then
!               a = meanA(nr_fine-1)
!               z = meanZ(nr_fine-1)
!            else
!               a = meanA(index+1)*(r-r_cc_loc(n,index))/dr_fine &
!                    + meanA(index)*(r_cc_loc(n,index+1)-r)/dr_fine
!               z = meanZ(index+1)*(r-r_cc_loc(n,index))/dr_fine &
!                    + meanZ(index)*(r_cc_loc(n,index+1)-r)/dr_fine
!            endif
!         else
!            if (index .eq. 0) then
!               a = meanA(index)
!               z = meanZ(index)
!            else if (index .gt. nr_fine-1) then
!               a = meanA(nr_fine-1)
!               z = meanZ(nr_fine-1)
!            else
!               a = meanA(index)*(r-r_cc_loc(n,index-1))/dr_fine &
!                    + meanA(index-1)*(r_cc_loc(n,index)-radius)/dr_fine
!               z = meanZ(index)*(r-r_cc_loc(n,index-1))/dr_fine &
!                    + meanZ(index-1)*(r_cc_loc(n,index)-radius)/dr_fine
!            end if
!         end if
                   
        ! quadratic interpolation

        ! index refers to the center point in the quadratic stencil.
        ! we need to modify this if we're too close to the edge
        if (index .eq. 0) then
           index = 1
        else if (index .ge. nr_fine-1) then
           index = nr_fine-2
        end if

        x0 = r_cc_loc(1,index-1)
        x1 = r_cc_loc(1,index)
        x2 = r_cc_loc(1,index+1)

        ! get Abar
        y0 = meanA(index-1) 
        y1 = meanA(index)
        y2 = meanA(index+1)

        a = y0 + (y1-y0)/(x1-x0)*(r-x0) &
           + ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0)*(r-x0)*(r-x1)

        if (a .gt. max(y0,y1,y2)) a = max(y0,y1,y2)
        if (a .lt. min(y0,y1,y2)) a = min(y0,y1,y2)

        ! get Zbar
        y0 = meanZ(index-1) 
        y1 = meanZ(index)
        y2 = meanZ(index+1)

        z = y0 + (y1-y0)/(x1-x0)*(r-x0) &
             + ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0)*(r-x0)*(r-x1)

        if (z .gt. max(y0,y1,y2)) z = max(y0,y1,y2)
        if (z .lt. min(y0,y1,y2)) z = min(y0,y1,y2)

      end subroutine interpolate

      subroutine quad_interp(r, a, z)
        use model_parser_module, only: meanA, meanZ
        use geometry, only: dr_fine, nr_fine, r_cc_loc, spherical
        use probin_module, only: max_levs

        implicit none

        real(kind=dp_t), intent(in   ) :: r
        real(kind=dp_t), intent(inout) :: a, z

        ! local
        real(kind=dp_t) :: x0,x1,x2
        real(kind=dp_t) :: y0,y1,y2
        real(kind=dp_t) :: z0,z1,z2
        integer :: i, id


        ! FIXME: need to think about this more in the event of restart
        if (spherical .eq. 1) then; n = 1
        else; n = max_levs
        endif

        id  = int(r / dr_fine)

        if ( id.eq.0 ) then 
           x0 = r_cc_loc(n,id  )
           x1 = r_cc_loc(n,id+1)
           x2 = r_cc_loc(n,id+2)
           y0 = meanA(id  )
           y1 = meanA(id+1)
           y2 = meanA(id+2)
           z0 = meanZ(id  )
           z1 = meanZ(id+1)
           z2 = meanZ(id+2)
        else if ( id.ge.nr_fine-1 ) then
           id = nr_fine-1
           x0 = r_cc_loc(n,id-2)
           x1 = r_cc_loc(n,id-1)
           x2 = r_cc_loc(n,id  )
           y0 = meanA(id-2)
           y1 = meanA(id-1)
           y2 = meanA(id  )
           z0 = meanZ(id-2)
           z1 = meanZ(id-1)
           z2 = meanZ(id  )
        else
           x0 = r_cc_loc(n,id-1)
           x1 = r_cc_loc(n,id  )
           x2 = r_cc_loc(n,id+1)
           y0 = meanA(id-1)
           y1 = meanA(id  )
           y2 = meanA(id+1)
           z0 = meanZ(id-1)
           z1 = meanZ(id  )
           z2 = meanZ(id+1)
        endif

        a = y0 + (y1-y0)/(x1-x0)*(r-x0) &
             + ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0)*(r-x0)*(r-x1)
        z = z0 + (z1-z0)/(x1-x0)*(r-x0) &
             + ((z2-z1)/(x2-x1)-(z1-z0)/(x1-x0))/(x2-x0)*(r-x0)*(r-x1)

        ! safety check
        if (a .gt. max(y0,y1,y2) .and. r.gt.r_cc_loc(n,1) ) a = max(y0,y1,y2)
        if (a .lt. min(y0,y1,y2)) a = min(y0,y1,y2)
        if (z .gt. max(z0,z1,z2) .and. r.gt.r_cc_loc(n,1) ) z = max(z0,z1,z2)
        if (z .lt. min(z0,z1,z2)) z = min(z0,z1,z2)

        if ( a .le. 0.d0 .OR. z .le. 0.d0 ) then

           write(*,*)'a, z, ', a, z
           write(*,*) 'index ', id
           stop

        end if

        return

      end subroutine quad_interp

  end subroutine eos

end module eos_module
