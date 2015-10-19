! setup a grid of rho, T, X.  Call the EOS to get the thermodynamic
! quantities.  Store a plotfile with all the information for checking.
! Then invert the results with the EOS to see if we recover the original
! T (i.e. use the different eos_input_* types).

subroutine varden()

  use BoxLib
  use f2kcli
  use ml_boxarray_module
  use layout_module
  use multifab_module
  use ml_cc_restriction_module
  use bl_mem_stat_module
  use bl_timer_module
  use box_util_module
  use bl_IO_module
  use fabio_module
  use variables
  use probin_module, only: prob_lo, prob_hi, &
                           test_set, pmask, max_levs, &
                           small_temp, &
                           temp_min, temp_max, &
                           dens_min, dens_max, &
                           metalicity_max, &
                           run_prefix
  use runtime_init_module
  use geometry, only: initialize_dx
  use bl_constants_module
  use network
  use eos_module, only: eos_input_rt, eos_input_rh, eos_input_tp, &
                        eos_input_rp, eos_input_re, eos_input_ps, eos, eos_init
  use eos_type_module

  implicit none

  integer :: i, n
  integer :: ii, jj, kk

  integer :: ng_s
  integer :: nlevs, dm

  type(ml_layout) :: mla

  type(multifab), allocatable :: s(:)

  real(kind=dp_t), pointer :: sp(:,:,:,:)

  real(kind=dp_t), pointer :: dx(:,:)

  integer, allocatable :: lo(:),hi(:)

  type(ml_boxarray) :: mba

  real(kind=dp_t) :: temp_zone, dens_zone, metalicity
  real(kind=dp_t) :: dlogrho, dlogT, dmetal
  real(kind=dp_t) :: xn_zone(nspec)

  integer :: ih1, ihe4

  type (eos_t) :: eos_state

  ! general Maestro initializations
  call runtime_init()

  ! microphysics
  call network_init()
  call eos_init(small_temp=small_temp)

  ! note: custom variables for this problem
  call init_variables()


  ! setup the grid
  call read_a_hgproj_grid(mba, test_set)

  call ml_layout_build(mla,mba,pmask)

  ! check for proper nesting
  if (.not. ml_boxarray_properly_nested(mla%mba, 3, pmask)) then
     call bl_error('ERROR: fixed_grids not properly nested')
  end if

  ! initialize nlevs
  nlevs = mla%nlevel

  if (nlevs .ne. 1) then
     call bl_error('ERROR: only 1 level of refinement supported')
  end if

  ! initialize dm
  dm = mla%dim

  if (dm /= 3) then
     call bl_error('ERROR: grid must be three-dimensional')
  endif


  ! allocate states
  allocate(s(nlevs))

  do n = 1,nlevs
     call multifab_build(s(n), mla%la(n), nscal, 1)
  end do


  ! initialize_dx
  call initialize_dx(dx,mba,nlevs)


  ! other allocations
  allocate(lo(dm))
  allocate(hi(dm))


  ! density, temperature, and metalicity increments
  dlogrho   = (log10(dens_max) - log10(dens_min))/(extent(mla%mba%pd(1),1) - 1)
  dlogT     = (log10(temp_max) - log10(temp_min))/(extent(mla%mba%pd(1),2) - 1)
  dmetal    = (metalicity_max  - ZERO           )/(extent(mla%mba%pd(1),3) - 1)


  ! initialize the thermodynamic cube and do the initial EOS call
  ng_s  = nghost(s(1))

  ih1 = network_species_index('hydrogen-1')
  ihe4 = network_species_index('helium-4')


  do n = 1, nlevs
     do i = 1, nfabs(s(n))

        sp => dataptr(s(n), i)

        lo = lwb(get_box(s(n), i))
        hi = upb(get_box(s(n), i))

        !$OMP PARALLEL DO PRIVATE(ii,jj,kk, metalicity, xn_zone, temp_zone, dens_zone, eos_state)
        do kk = lo(3), hi(3)

           ! set the composition -- approximately solar
           metalicity = ZERO + dble(kk)*dmetal
           xn_zone(:) = metalicity/(nspec - 2)   ! all but H, He
           xn_zone(ih1)  = 0.75_dp_t - HALF*metalicity
           xn_zone(ihe4) = 0.25_dp_t - HALF*metalicity

           do jj = lo(2), hi(2)

              ! set the temperature
              temp_zone = 10.0**(log10(temp_min) + dble(jj)*dlogT)

              do ii = lo(1), hi(1)
                 
                 ! set the density
                 dens_zone = 10.0**(log10(dens_min) + dble(ii)*dlogrho)

                 !------------------------------------------------------------
                 ! call the EOS -- rho, T, X directly
                 ! input: eos_input_rt
                 !------------------------------------------------------------
                 eos_state%T     = temp_zone
                 eos_state%rho   = dens_zone
                 eos_state%xn(:) = xn_zone(:)

                 call eos(eos_input_rt, eos_state)

                 ! store the thermodynamic state
                 sp(ii,jj,kk,rho_comp) = dens_zone
                 sp(ii,jj,kk,rhoh_comp) = dens_zone*eos_state%h
                 sp(ii,jj,kk,temp_comp) = temp_zone
                 sp(ii,jj,kk,spec_comp:spec_comp-1+nspec) = xn_zone(:)

                 sp(ii,jj,kk,h_comp) = eos_state%h
                 sp(ii,jj,kk,p_comp) = eos_state%p
                 sp(ii,jj,kk,e_comp) = eos_state%e
                 sp(ii,jj,kk,s_comp) = eos_state%s


                 !------------------------------------------------------------
                 ! call the EOS with rho, h
                 ! input: eos_input_rh
                 !------------------------------------------------------------

                 ! change initial T guess to make the root find do some work
                 eos_state%T     = HALF*temp_zone   
                 eos_state%rho   = dens_zone
                 eos_state%h     = sp(ii,jj,kk,h_comp)
                 eos_state%xn(:) = xn_zone(:)

                 call eos(eos_input_rh, eos_state)

                 ! store the thermodynamic state
                 sp(ii,jj,kk,tfromrh_comp) = eos_state%T
                 sp(ii,jj,kk,tfromrh_err_comp) = &
                      (eos_state%T - temp_zone)/temp_zone


                 !------------------------------------------------------------
                 ! call the EOS with T, p
                 ! input: eos_input_tp
                 !------------------------------------------------------------

                 ! change initial rho guess to make the root find do some work
                 eos_state%T     = temp_zone   
                 eos_state%rho   = THIRD*dens_zone
                 eos_state%p     = sp(ii,jj,kk,p_comp)
                 eos_state%xn(:) = xn_zone(:)

                 call eos(eos_input_tp, eos_state)

                 ! store the thermodynamic state
                 sp(ii,jj,kk,rfromtp_comp) = eos_state%rho
                 sp(ii,jj,kk,rfromtp_err_comp) = &
                      (eos_state%rho - dens_zone)/dens_zone


                 !------------------------------------------------------------
                 ! call the EOS with rho, p
                 ! input: eos_input_rp
                 !------------------------------------------------------------

                 ! change initial T guess to make the root find do some work
                 eos_state%T     = HALF*temp_zone   
                 eos_state%rho   = dens_zone
                 eos_state%p     = sp(ii,jj,kk,p_comp)
                 eos_state%xn(:) = xn_zone(:)

                 call eos(eos_input_rp, eos_state)

                 ! store the thermodynamic state
                 sp(ii,jj,kk,tfromrp_comp) = eos_state%T
                 sp(ii,jj,kk,tfromrp_err_comp) = &
                      (eos_state%T - temp_zone)/temp_zone


                 !------------------------------------------------------------
                 ! call the EOS with rho, e
                 ! input: eos_input_re
                 !------------------------------------------------------------

                 ! change initial T guess to make the root find do some work
                 eos_state%T     = HALF*temp_zone   
                 eos_state%rho   = dens_zone
                 eos_state%e     = sp(ii,jj,kk,e_comp)
                 eos_state%xn(:) = xn_zone(:)

                 call eos(eos_input_re, eos_state)

                 ! store the thermodynamic state
                 sp(ii,jj,kk,tfromre_comp) = eos_state%T
                 sp(ii,jj,kk,tfromre_err_comp) = &
                      (eos_state%T - temp_zone)/temp_zone


                 !------------------------------------------------------------
                 ! call the EOS with p, s
                 ! input: eos_input_ps
                 !------------------------------------------------------------

                 ! change initial T and rho guess to make the root find do 
                 ! some work
                 eos_state%T     = HALF*temp_zone   
                 eos_state%rho   = HALF*dens_zone
                 eos_state%p     = sp(ii,jj,kk,p_comp)
                 eos_state%s     = sp(ii,jj,kk,s_comp)
                 eos_state%xn(:) = xn_zone(:)

                 ! some EOSes don't have physically valid treatments
                 ! of entropy throughout the entire rho-T plane
                 if (eos_state%s > ZERO) then

                    call eos(eos_input_ps, eos_state)

                    ! store the thermodynamic state
                    sp(ii,jj,kk,tfromps_comp) = eos_state%T
                    sp(ii,jj,kk,tfromps_err_comp) = &
                         (eos_state%T - temp_zone)/temp_zone
                    
                 endif

              enddo
           enddo
        enddo
        !$OMP END PARALLEL DO

     enddo
  enddo


  ! write out a plotfile that contains the magnitude of the velocity
  ! field
  print *, mla%mba%rr(:,1)
  call fabio_ml_multifab_write_d(s,mla%mba%rr(:,1), &
                                 trim(run_prefix) // "eos_thermo", &
                                 names=varnames)



  ! clean-up
  do n = 1,nlevs
     call destroy(s(n))
  end do

  call destroy(mla)
  call destroy(mba)

  deallocate(s)

  call runtime_close()

end subroutine varden



