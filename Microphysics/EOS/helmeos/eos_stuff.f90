module eos_module

  use bl_space, only: MAX_SPACEDIM
  use bl_types
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module

  implicit none

  private 


  integer, parameter, public :: eos_input_rt = 1  ! rho, T are inputs
  integer, parameter, public :: eos_input_rh = 2  ! rho, h are inputs
  integer, parameter, public :: eos_input_tp = 3  ! T, p are inputs
  integer, parameter, public :: eos_input_rp = 4  ! rho, p are inputs
  integer, parameter, public :: eos_input_re = 5  ! rho, e are inputs
  integer, parameter, public :: eos_input_ps = 6  ! p, s are inputs
  integer, parameter, public :: eos_input_ph = 7  ! p, h are inputs
 

  logical, save, private :: do_coulomb
  real(kind=dp_t), save, private :: smallt
  real(kind=dp_t), save, private :: smalld

  logical, save, private :: initialized = .false.

  private nspec, aion, zion

  public eos_init, eos_finalize, eos

  interface eos
     module procedure eos_new
  end interface eos

contains

  ! EOS initialization routine -- this is used by both MAESTRO and Castro
  ! For this general EOS, this calls helmeos_init() which reads in the 
  ! table with the electron component's properties.
  subroutine eos_init(small_temp, small_dens, gamma_in)

    use parallel
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
       smallt = 1.d4
      end if
    else
       smallt = 1.d4
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

    if (parallel_IOProcessor()) print *, 'Initializing helmeos...'
    ! call the helmeos initialization routine and read in the table 
    ! containing the electron contribution.
    call helmeos_init()
    initialized = .true.
 
  end subroutine eos_init


  subroutine eos_finalize()

  end subroutine eos_finalize


  !---------------------------------------------------------------------------
  ! main interface
  !---------------------------------------------------------------------------
  subroutine eos_new(input, eos_state, do_eos_diag, pt_index)

    use bl_error_module
    
    integer,           intent(in   ) :: input
    type (eos_t),      intent(inout) :: eos_state
    logical, optional, intent(in   ) :: do_eos_diag
    integer, optional, intent(in   ) :: pt_index(:)


! a generic wrapper for the Helmholtz based electron/positron degenerate
! EOS.  All the fields are through the eos_t type.  
!
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
! dedR     -- d energy/ d density
! dpdX     -- d pressure / d xmass
! dhdX     -- d enthalpy / d xmass  -- AT CONSTANT PRESSURE!!!
! gam1     -- first adiabatic index (d log P/ d log rho) |_s
! cs       -- sound speed -- note that this is the non-relativistic one
!             (we compute it in this wrapper as sqrt(gam1 p /rho) instead
!             of taking the relativistic version from helmeos.
! s        -- entropy (erg/g/K)
! dsdT     -- d entropy / d temperature
! dsdR     -- d entropy / d density
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

    
!     ::::: Local variables and arrays

    integer :: i, n, iter, niter, max_newton
    parameter (max_newton = 100)
    
    real(kind=dp_t) :: error, error2
    real(kind=dp_t) :: ymass(nspec)
    real(kind=dp_t) :: abar, zbar
    real(kind=dp_t) :: energy_want
    real(kind=dp_t) :: enthalpy_want
    real(kind=dp_t) :: pres_want
    real(kind=dp_t) :: entropy_want
    real(kind=dp_t) :: dhdt
    real(kind=dp_t) :: tnew
    real(kind=dp_t) :: dnew
    real(kind=dp_t) :: enth1
    real(kind=dp_t) :: ener1
    real(kind=dp_t) :: dedX(nspec)

    real(kind=dp_t) :: dpdd, pres1, entr1
    real(kind=dp_t) :: f, g, dfdd, dfdt, dgdd, dgdt, deld

    real(kind=dp_t), parameter :: ttol = 1.0d-8
    real(kind=dp_t), parameter :: dtol = 1.0d-8
    real(kind=dp_t), parameter :: stol = 1.0d-8

    logical eosfail
    integer dim_ptindex

    ! err_string is used to convert the pt_index information into a string
    character (len=64) :: err_string  

!     ::::: Input/Output arrays for call to helmeos
    real(kind=dp_t) :: temp_row, den_row, abar_row, &
                     zbar_row, etot_row, ptot_row, &
                     cv_row, cp_row, &
                     xne_row, xnp_row, etaele_row, &
                     pele_row, ppos_row, dpd_row, &
                     dpt_row, dpa_row, dpz_row, &
                     ded_row, det_row, dea_row, &
                     dez_row, &
                     stot_row, dsd_row, dst_row
    real(kind=dp_t) :: gam1_row, cs_row

    
    if (present(pt_index)) dim_ptindex = size(pt_index,dim=1)
      
    if (.not. initialized) call bl_error('EOS: not initialized')

    ! this format statement is for writing into err_string -- make sure that
    ! the len of err_string can accomodate this format specifier
1001 format(1x,"zone index info: i = ", i5)
1002 format(1x,"zone index info: i = ", i5, '  j = ', i5)
1003 format(1x,"zone index info: i = ", i5, '  j = ', i5, '  k = ', i5)

    tnew  = 0.0d0
    dnew   = 0.0d0

    do i=1,nspec
       ymass(i) = eos_state%xn(i)/aion(i)
       dnew    = dnew + ymass(i)
       tnew    = tnew + zion(i) * ymass(i)
    enddo

    abar = 1.0d0/dnew
    zbar = tnew * abar

    if (input .EQ. eos_input_rt) then

!---------------------------------------------------------------------------
! input = 1: dens, temp, and xmass are inputs
!---------------------------------------------------------------------------

! we are taking density, temperature, and composition as given
       temp_row = eos_state%T
       den_row  = eos_state%rho
       abar_row = abar
       zbar_row = zbar
         
! call the eos
       call helmeos(do_coulomb,eosfail, &
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
       eos_state%p = ptot_row
       eos_state%e = etot_row
         
       eos_state%h = eos_state%e + eos_state%p/eos_state%rho

       eos_state%cv = cv_row
       eos_state%cp = cp_row

       ! store the number density of electrons and positrons, the degeneracy
       ! parameter, and the total electron/positron pressure
       eos_state%xne   = xne_row + xnp_row
       eos_state%eta  = etaele_row
       eos_state%pele = pele_row + ppos_row

       eos_state%dpdr = dpd_row
       eos_state%dpdT = dpt_row
       eos_state%dedr = ded_row
       eos_state%dedT = det_row
       eos_state%gam1 = gam1_row
!       eos_state%cs =   cs_row
       eos_state%cs =   sqrt(eos_state%gam1*eos_state%p/eos_state%rho)
       eos_state%s = stot_row
       eos_state%dsdT = dst_row
       eos_state%dsdr = dsd_row

       do n = 1, nspec
          eos_state%dpdX(n) = dpa_row * (abar/aion(n))* &
                              (aion(n) - abar) + &
                    dpz_row * (abar/aion(n))* &
                              (zion(n) - zbar)

          dEdX(n) = dea_row * (abar/aion(n))* &
                              (aion(n) - abar) + &
                    dez_row * (abar/aion(n))* &
                              (zion(n) - zbar)

          ! create the enthalpy derivatives wrt average composition --
          ! hold pressure constant!!!
          eos_state%dhdX(n) = dEdX(n) + &
               (eos_state%p/eos_state%rho**2 - eos_state%dedr)*eos_state%dpdX(n)/eos_state%dpdr

       enddo

    else if (input .EQ. eos_input_rh) then

!---------------------------------------------------------------------------
! input = 2: dens, enthalpy, and xmass are inputs
!---------------------------------------------------------------------------

       ! load the initial guess
       temp_row = eos_state%T
       den_row  = eos_state%rho
       abar_row = abar
       zbar_row = zbar

       if (do_eos_diag) print*,'T/D INIT ',eos_state%T,eos_state%rho

       ! we want to converge to the given enthalpy
       enthalpy_want = eos_state%h

       if (do_eos_diag) print*,'WANT H ',eos_state%h

       call helmeos(do_coulomb,eosfail, &
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
          enth1 = etot_row + ptot_row/eos_state%rho
          if (do_eos_diag) print*,'ENTH1 ',iter,enth1

          dhdt = det_row + dpt_row/eos_state%rho
          if (do_eos_diag) print*,'DHDT ',iter,dhdt

          tnew = temp_row - &
               (enth1 - enthalpy_want)/dhdt

          if (do_eos_diag) then
             print *, 'TNEW FIRST ', temp_row, ' - ', &
                  enth1 - enthalpy_want, ' / ', dhdt
          endif

          ! don't let the temperature change by more than a factor of two
          tnew = max(.5d0*temp_row, &
                        min(tnew, 2.d0*temp_row))

          ! don't let us freeze
          tnew = max(smallt, tnew)

          if (do_eos_diag) print*,'TNEW AFTER ',iter,tnew

          ! compute the error
          error = 0.0d0
          error = max(error,abs(tnew - temp_row)/temp_row)

          ! store the new temperature
          temp_row = tnew

          if (error .LT. ttol) goto 70
        
          call helmeos(do_coulomb,eosfail, &
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
       eos_state%T = tnew
       eos_state%p = ptot_row
       eos_state%e = etot_row
       
       eos_state%cv = cv_row
       eos_state%cp = cp_row

       ! store the number density of electrons and positrons, the degeneracy
       ! parameter, and the total electron/positron pressure
       eos_state%xne   = xne_row + xnp_row
       eos_state%eta  = etaele_row
       eos_state%pele = pele_row + ppos_row
       
       eos_state%dpdr = dpd_row
       eos_state%dpdT = dpt_row
       eos_state%dedr = ded_row
       eos_state%dedT = det_row   ! c_v
       eos_state%gam1 = gam1_row
!      eos_state%cs =   cs_row
       eos_state%cs =   sqrt(eos_state%gam1*eos_state%p/eos_state%rho)
       eos_state%s = stot_row
       eos_state%dsdT = dst_row
       eos_state%dsdr = dsd_row

       do n = 1, nspec
          eos_state%dpdX(n) = dpa_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dpz_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          dEdX(n) = dea_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dez_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          ! create the enthalpy derivatives wrt average composition --
          ! hold pressure constant!!!
          eos_state%dhdX(n) = dEdX(n) + &
               (eos_state%p/eos_state%rho**2 - eos_state%dedr)*eos_state%dpdX(n)/eos_state%dpdr

       enddo

    else if (input .EQ. eos_input_tp ) then

!---------------------------------------------------------------------------
! input = 3: temp, pres, and xmass are inputs
!---------------------------------------------------------------------------

       ! load the initial guess
       temp_row = eos_state%T
       den_row  = eos_state%rho
       abar_row = abar
       zbar_row = zbar

       ! we want to converge to the given pressure
       pres_want = eos_state%p

       if (pres_want < ZERO) then
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
         
       call helmeos(do_coulomb,eosfail, &
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
          pres1 = ptot_row
          dpdd  = dpd_row
          
          dnew = den_row - &
               (pres1 - pres_want)/dpdd

          ! don't let the density change by more than an order of magnitude
          dnew = max(.5d0*den_row, &
                        min(dnew, 2.d0*den_row))

          ! compute the error
          error = 0.0d0
          error = max(error,abs(dnew - den_row)/den_row)

          ! store the new density
          den_row = dnew

          ! check if we are evacuating, if so, set the density to smalld, and adjust
          ! the error so we iterate on this one
          if (den_row .LT. smalld) then
             den_row = smalld
             error = 1.1d0*dtol
          endif

          if (error .LT. dtol) goto 170
        
          call helmeos(do_coulomb,eosfail, &
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
       eos_state%rho = dnew
       eos_state%T = temp_row
       eos_state%e = etot_row
       eos_state%h = eos_state%e + ptot_row/eos_state%rho

       eos_state%cv = cv_row
       eos_state%cp = cp_row

       ! store the number density of electrons and positrons, the degeneracy
       ! parameter, and the total electron/positron pressure
       eos_state%xne   = xne_row + xnp_row
       eos_state%eta  = etaele_row
       eos_state%pele = pele_row + ppos_row
       
       eos_state%dpdr = dpd_row
       eos_state%dpdT = dpt_row
       eos_state%dedr = ded_row
       eos_state%dedT = det_row   ! c_v
       eos_state%gam1 = gam1_row
!      eos_state%cs =   cs_row
       eos_state%cs =   sqrt(eos_state%gam1*eos_state%p/eos_state%rho)
       eos_state%s = stot_row
       eos_state%dsdT = dst_row
       eos_state%dsdr = dsd_row

       do n = 1, nspec
          eos_state%dpdX(n) = dpa_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dpz_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          dEdX(n) = dea_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dez_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          ! create the enthalpy derivatives wrt average composition --
          ! hold pressure constant!!!
          eos_state%dhdX(n) = dEdX(n) + &
               (eos_state%p/eos_state%rho**2 - eos_state%dedr)*eos_state%dpdX(n)/eos_state%dpdr

       enddo

    else if (input .EQ. eos_input_rp ) then

!---------------------------------------------------------------------------
! input = 4: dens, pres, and xmass are inputs
!---------------------------------------------------------------------------

       ! Load the initial guess
       temp_row = eos_state%T
       den_row  = eos_state%rho
       abar_row = abar
       zbar_row = zbar
       if (do_eos_diag) print*,'T/D INIT ',eos_state%T,eos_state%rho

       ! We want to converge to the given pressure
       pres_want = eos_state%p
       if (do_eos_diag) print*,'P WANT ',eos_state%p

       if (pres_want < ZERO) then
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

       ! First pass
       call helmeos(do_coulomb,eosfail, &
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
          tnew = temp_row - &
               (ptot_row - pres_want)/dpt_row

          if (do_eos_diag) print*,'PRES ',ptot_row,pres_want
          if (do_eos_diag) print*,'PRES DIFF ',ptot_row-pres_want
          
          if (do_eos_diag) print*,'DPDT FAC ', 1.0/dpt_row

          if (do_eos_diag) print*,'TNEW BEFORE MAX ',iter,tnew

          ! don't let the temperature change by more than a factor of 2
          tnew = max(.5d0*temp_row, &
                        min(tnew, 2.d0*temp_row))

          ! don't let us freeze
          tnew = max(smallt, tnew)

          if (do_eos_diag) print*,'TNEW AFTER MAX ',iter,tnew
          if (do_eos_diag) print*,' '

          ! compute the error and store the new temperature
          error = 0.0d0
          error = max(error,abs(tnew - temp_row)/temp_row)
          if (do_eos_diag) print *,'ERROR  ',iter,error
          if (do_eos_diag) print*,' '
          temp_row = tnew

          if (error .LT. ttol) goto 870
        
          call helmeos(do_coulomb,eosfail, &
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
       print *, 'error, temp_row, den_row = ', error, temp_row, den_row

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
       ! jbb
       ! temp = tnew
       eos_state%T = temp_row
       eos_state%rho = den_row
       eos_state%e = etot_row
       eos_state%cv = cv_row
       eos_state%cp = cp_row

       eos_state%h = eos_state%e + ptot_row/eos_state%rho

       ! store the number density of electrons and positrons, the degeneracy
       ! parameter, and the total electron/positron pressure
       eos_state%xne   = xne_row + xnp_row
       eos_state%eta  = etaele_row
       eos_state%pele = pele_row + ppos_row
       
       eos_state%dpdr = dpd_row
       eos_state%dpdT = dpt_row
       eos_state%dedr = ded_row
       eos_state%dedT = det_row   ! c_v
       eos_state%gam1 = gam1_row
!      eos_state%cs =   cs_row
       eos_state%cs =   sqrt(eos_state%gam1*eos_state%p/eos_state%rho)
       eos_state%s = stot_row
       eos_state%dsdT = dst_row
       eos_state%dsdr = dsd_row

       do n = 1, nspec
          eos_state%dpdX(n) = dpa_row * (abar/aion(n))* &
                              (aion(n) - abar) + &
                    dpz_row * (abar/aion(n))* &
                              (zion(n) - zbar)

          dEdX(n) = dea_row * (abar/aion(n))* &
                              (aion(n) - abar) + &
                    dez_row * (abar/aion(n))* &
                              (zion(n) - zbar)

          ! create the enthalpy derivatives wrt average composition --
          ! hold pressure constant!!!
          eos_state%dhdX(n) = dEdX(n) + &
               (eos_state%p/eos_state%rho**2 - eos_state%dedr)*eos_state%dpdX(n)/eos_state%dpdr

       enddo

    else if (input .EQ. eos_input_re) then

!---------------------------------------------------------------------------
! input = 5: dens, energy, and xmass are inputs
!---------------------------------------------------------------------------

!      do_eos_diag = .true.
       ! load the initial guess
       temp_row = eos_state%T
       den_row  = eos_state%rho
       abar_row = abar
       zbar_row = zbar

       if (do_eos_diag) print*,'T/D INIT ',eos_state%T,eos_state%rho

       ! we want to converge to the given energy
       energy_want = eos_state%e

       if (energy_want < ZERO) then
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

       if (do_eos_diag) print*,'WANT e ',energy_want

       call helmeos(do_coulomb,eosfail, &
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
          ener1 = etot_row

          if (do_eos_diag) then
             print*,'ENER1 ',iter,ener1
             print*,'DEDT ',iter,det_row
          end if

          tnew = temp_row - (ener1-energy_want)/det_row

          if (do_eos_diag) then
             print *, 'TNEW FIRST ', temp_row, ' - ', &
                  ener1 - energy_want, ' / ', det_row
          endif

          ! don't let the temperature change by more than a factor of two
          tnew = max(.5d0*temp_row, &
                        min(tnew, 2.d0*temp_row))

          ! don't let us freeze
          tnew = max(smallt, tnew)

          if (do_eos_diag) print*,'TNEW AFTER ',iter,tnew

          ! compute the error
          error = 0.0d0
          error = max(error,abs(tnew - temp_row)/temp_row)

          ! store the new temperature
          temp_row = tnew

          if (error .LT. ttol) goto 270
        
          call helmeos(do_coulomb,eosfail, &
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
          call bl_error('EOS: Newton-Raphson failed:5: too many iterations', err_string)
       else
          call bl_error('EOS: Newton-Raphson failed:5: too many iterations')
       endif

270     continue

       ! store the end result
       eos_state%T = tnew
       eos_state%p = ptot_row
       eos_state%h = eos_state%e + ptot_row/eos_state%rho
       
       eos_state%cv = cv_row
       eos_state%cp = cp_row

       ! store the number density of electrons and positrons, the degeneracy
       ! parameter, and the total electron/positron pressure
       eos_state%xne   = xne_row + xnp_row
       eos_state%eta  = etaele_row
       eos_state%pele = pele_row + ppos_row
       
       eos_state%dpdr = dpd_row
       eos_state%dpdT = dpt_row
       eos_state%dedr = ded_row
       eos_state%dedT = det_row   ! c_v
       eos_state%gam1 = gam1_row
!      eos_state%cs =   cs_row
       eos_state%cs =   sqrt(eos_state%gam1*eos_state%p/eos_state%rho)
       eos_state%s = stot_row
       eos_state%dsdT = dst_row
       eos_state%dsdr = dsd_row

       do n = 1, nspec
          eos_state%dpdX(n) = dpa_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dpz_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          dEdX(n) = dea_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dez_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          ! create the enthalpy derivatives wrt average composition --
          ! hold pressure constant!!!
          eos_state%dhdX(n) = dEdX(n) + &
               (eos_state%p/eos_state%rho**2 - eos_state%dedr)*eos_state%dpdX(n)/eos_state%dpdr

       enddo

    else if (input .EQ. eos_input_ps) then
!---------------------------------------------------------------------------
! input = 6: pres, entropy, and xmass are inputs
!---------------------------------------------------------------------------

       ! load the initial guess
       temp_row = eos_state%T
       den_row  = eos_state%rho
       abar_row = abar
       zbar_row = zbar

       if (do_eos_diag) print*,'T/D INIT ',eos_state%T,eos_state%rho

       ! we want to converge to the given entropy and pressure
       entropy_want = eos_state%s
       pres_want    = eos_state%p

       if (entropy_want < ZERO) then
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

       if (pres_want < ZERO) then
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
          print*,'WANT s ',entropy_want
          print*,'WANT pres', pres_want
       endif

       call helmeos(do_coulomb,eosfail, &
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
          pres1 = ptot_row
          entr1 = stot_row

          if (do_eos_diag) then
             print*,'PRES1 ',iter,pres1
             print*,'ENTR1 ',iter,entr1
          end if

          ! two functions, f and g, to iterate over
          f = pres_want - pres1
          dfdd = -dpd_row
          dfdt = -dpt_row
          
          g = entropy_want - entr1
          dgdd = -dsd_row
          dgdt = -dst_row
          !
          ! 0 = f + dfdd * deld + dfdt * delt
          ! 0 = g + dgdd * deld + dgdt * delt
          !
          deld = (f*dgdt - g*dfdt) / (dgdd*dfdt - dgdt*dfdd)

          dnew = den_row + deld

          tnew = temp_row - (f + dfdd*deld) / dfdt

          if (do_eos_diag) then
             print *, 'DNEW FIRST ', den_row, ' + ', &
                  f*dgdt - g*dfdt, ' / ', dgdd*dfdt - dgdt*dfdd
             print *, 'TNEW FIRST ', temp_row, ' - ', &
                  f + dfdd*deld, ' / ', dfdt
          endif

          ! don't let the temperature or density change by more
          ! than a factor of two
          tnew = max(HALF*temp_row, &
                        min(tnew, TWO*temp_row))
          dnew = max(HALF*den_row, &
                        min(dnew, TWO*den_row))

          ! don't let us freeze or evacuate
          tnew = max(smallt, tnew)
          dnew = max(smalld, dnew)

          if (do_eos_diag) then
             print*,'DNEW AFTER ',iter,dnew
             print*,'TNEW AFTER ',iter,tnew
          endif

          ! compute the errors
          error = ZERO
          error2 = ZERO
          error  = max(error ,abs(dnew - den_row)/den_row)
          error2 = max(error2,abs(tnew - temp_row)/temp_row)

          ! store the new temperature and density
          den_row = dnew
          temp_row = tnew

          if (error .LT. dtol .and. error2 .LT. ttol) goto 370
     
          call helmeos(do_coulomb,eosfail, &
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
          call bl_error('EOS: Newton-Raphson failed:6: too many iterations', err_string)
       else
          call bl_error('EOS: Newton-Raphson failed:6: too many iterations')
       endif

370     continue

       ! store the end result
       eos_state%rho = dnew
       eos_state%T = tnew
       eos_state%e = etot_row
       eos_state%h = eos_state%e + ptot_row/eos_state%rho
          
       eos_state%cv = cv_row
       eos_state%cp = cp_row

       ! store the number density of electrons and positrons, the degeneracy
       ! parameter, and the total electron/positron pressure
       eos_state%xne   = xne_row + xnp_row
       eos_state%eta  = etaele_row
       eos_state%pele = pele_row + ppos_row
       
       eos_state%dpdr = dpd_row
       eos_state%dpdT = dpt_row
       eos_state%dedr = ded_row
       eos_state%dedT = det_row   ! c_v
       eos_state%gam1 = gam1_row
!      eos_state%cs =   cs_row
       eos_state%cs =   sqrt(eos_state%gam1*eos_state%p/eos_state%rho)
       eos_state%dsdT = dst_row
       eos_state%dsdr = dsd_row

       do n = 1, nspec
          eos_state%dpdX(n) = dpa_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dpz_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          dEdX(n) = dea_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dez_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          ! create the enthalpy derivatives wrt average composition --
          ! hold pressure constant!!!
          eos_state%dhdX(n) = dEdX(n) + &
               (eos_state%p/eos_state%rho**2 - eos_state%dedr)*eos_state%dpdX(n)/eos_state%dpdr

       enddo


    else if (input .EQ. eos_input_ph) then
!---------------------------------------------------------------------------
! input = 7: pres, enthalpy, and xmass are inputs
!---------------------------------------------------------------------------

       ! load the initial guess
       temp_row = eos_state%T
       den_row  = eos_state%rho
       abar_row = abar
       zbar_row = zbar

       if (do_eos_diag) print*,'T/D INIT ',eos_state%T,eos_state%rho

       ! we want to converge to the given entropy and pressure
       enthalpy_want = eos_state%h
       pres_want    = eos_state%p

       if (enthalpy_want < ZERO) then
          if (present(pt_index)) then
             if (dim_ptindex .eq. 1) then 
                write (err_string,1001) pt_index(1)
             else if (dim_ptindex .eq. 2) then 
                write (err_string,1002) pt_index(1), pt_index(2)
             else if (dim_ptindex .eq. 3) then 
                write (err_string,1003) pt_index(1), pt_index(2), pt_index(3)
             end if
             call bl_error('ERROR: enthalpy < 0 in the EOS', err_string)
          else
             call bl_error('ERROR: enthalpy < 0 in the EOS')
          endif
       endif

       if (pres_want < ZERO) then
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
          print*,'WANT h ',enthalpy_want
          print*,'WANT pres', pres_want
       endif

       call helmeos(do_coulomb,eosfail, &
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
          pres1 = ptot_row
          enth1 = etot_row + ptot_row/den_row

          if (do_eos_diag) then
             print*,'PRES1 ',iter,pres1
             print*,'ENTH1 ',iter,enth1
          end if

          ! two functions, f and g, to iterate over
          f = pres_want - pres1
          dfdd = -dpd_row
          dfdt = -dpt_row
          
          g = enthalpy_want - enth1
          dgdd = -ded_row - dpd_row/den_row + ptot_row/den_row**2
          dgdt = -det_row - dpt_row/den_row
          !
          ! 0 = f + dfdd * deld + dfdt * delt
          ! 0 = g + dgdd * deld + dgdt * delt
          !
          deld = (f*dgdt - g*dfdt) / (dgdd*dfdt - dgdt*dfdd)

          dnew = den_row + deld

          tnew = temp_row - (f + dfdd*deld) / dfdt

          if (do_eos_diag) then
             print *, 'DNEW FIRST ', den_row, ' + ', &
                  f*dgdt - g*dfdt, ' / ', dgdd*dfdt - dgdt*dfdd
             print *, 'TNEW FIRST ', temp_row, ' - ', &
                  f + dfdd*deld, ' / ', dfdt
          endif

          ! don't let the temperature or density change by more
          ! than a factor of two
          tnew = max(HALF*temp_row, &
                        min(tnew, TWO*temp_row))
          dnew = max(HALF*den_row, &
                        min(dnew, TWO*den_row))

          ! don't let us freeze or evacuate
          tnew = max(smallt, tnew)
          dnew = max(smalld, dnew)

          if (do_eos_diag) then
             print*,'DNEW AFTER ',iter,dnew
             print*,'TNEW AFTER ',iter,tnew
          endif

          ! compute the errors
          error = ZERO
          error2 = ZERO
          error  = max(error ,abs(dnew - den_row)/den_row)
          error2 = max(error2,abs(tnew - temp_row)/temp_row)

          ! store the new temperature and density
          den_row = dnew
          temp_row = tnew

          if (error .LT. dtol .and. error2 .LT. ttol) goto 470
     
          call helmeos(do_coulomb,eosfail, &
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
          call bl_error('EOS: Newton-Raphson failed:7: too many iterations', err_string)
       else
          call bl_error('EOS: Newton-Raphson failed:7: too many iterations')
       endif

470     continue

       ! store the end result
       eos_state%rho = dnew
       eos_state%T = tnew
       eos_state%e = etot_row
       eos_state%s = stot_row
          
       eos_state%cv = cv_row
       eos_state%cp = cp_row

       ! store the number density of electrons and positrons, the degeneracy
       ! parameter, and the total electron/positron pressure
       eos_state%xne   = xne_row + xnp_row
       eos_state%eta  = etaele_row
       eos_state%pele = pele_row + ppos_row
       
       eos_state%dpdr = dpd_row
       eos_state%dpdT = dpt_row
       eos_state%dedr = ded_row
       eos_state%dedT = det_row   ! c_v
       eos_state%gam1 = gam1_row
!      eos_state%cs =   cs_row
       eos_state%cs =   sqrt(eos_state%gam1*eos_state%p/eos_state%rho)
       eos_state%dsdT = dst_row
       eos_state%dsdr = dsd_row

       do n = 1, nspec
          eos_state%dpdX(n) = dpa_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dpz_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          dEdX(n) = dea_row * (abar/aion(n))* &
                                 (aion(n) - abar) + &
                    dez_row * (abar/aion(n))* &
                                 (zion(n) - zbar)

          ! create the enthalpy derivatives wrt average composition --
          ! hold pressure constant!!!
          eos_state%dhdX(n) = dEdX(n) + &
               (eos_state%p/eos_state%rho**2 - eos_state%dedr)*eos_state%dpdX(n)/eos_state%dpdr

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
  end subroutine eos_new

end module eos_module
