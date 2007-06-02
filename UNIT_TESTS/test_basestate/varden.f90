subroutine varden()

  use BoxLib
  use omp_module
  use f2kcli
  use list_box_module
  use ml_boxarray_module
  use layout_module
  use multifab_module
  use init_module
  use box_util_module
  use bl_IO_module
  use variables
  use geometry
  use network
  use eos_module

  implicit none

  integer :: i

  integer    :: narg, farg
  integer    :: max_step
  integer    :: dm,n_base
  integer    :: nscal, ntrac

  real(dp_t) :: dr_base
  real(dp_t) :: cflfac
  real(dp_t) :: stop_time
  real(dp_t) :: time,dt,half_time
  real(dp_t) :: prob_lo_x, prob_hi_x
  real(dp_t) :: anelastic_cutoff
  real(dp_t) :: smin,smax
  integer    :: spherical_in

  character(len=128) :: fname
  character(len=128) :: probin_env
  
  character(len=256) :: model_file

  integer :: un, ierr
  logical :: lexist
  logical :: need_inputs

  real(dp_t) :: y_0

  real(dp_t) :: dx(1)
  real(dp_t) :: prob_lo(1), prob_hi(1)

  real(dp_t), allocatable :: div_coeff_old(:)
  real(dp_t), allocatable :: div_coeff_new(:)
  real(dp_t), allocatable :: grav_cell(:)
  real(dp_t), allocatable :: gam1(:)
  real(dp_t), allocatable :: s0_old(:,:)
  real(dp_t), allocatable :: s0_new(:,:)
  real(dp_t), allocatable :: s0_1(:,:)
  real(dp_t), allocatable :: s0_2(:,:)
  real(dp_t), allocatable :: temp0(:)
  real(dp_t), allocatable :: p0_old(:)
  real(dp_t), allocatable :: p0_1(:)
  real(dp_t), allocatable :: p0_2(:)
  real(dp_t), allocatable :: p0_new(:)
  real(dp_t), allocatable :: w0(:)
  real(dp_t), allocatable :: Sbar_in(:)

  real(dp_t) :: coeff, Hbar

  namelist /probin/ model_file
  namelist /probin/ stop_time
  namelist /probin/ prob_lo_x
  namelist /probin/ prob_hi_x
  namelist /probin/ max_step
  namelist /probin/ cflfac
  namelist /probin/ anelastic_cutoff
  namelist /probin/ spherical_in
  namelist /probin/ n_base


  narg = command_argument_count()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize the runtime parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Defaults
  model_file = "model.hse"

  spherical_in = 0
  max_step  = 1

  prob_lo_x = ZERO
  prob_lo_x = ONE

  anelastic_cutoff = 3.e6

  ntrac = 1
  nscal = nspec + ntrac + 2

  n_base = 512

  need_inputs = .true.


  call get_environment_variable('PROBIN', probin_env, status = ierr)
  if ( need_inputs .AND. ierr == 0 ) then
     un = unit_new()
     open(unit=un, file = probin_env, status = 'old', action = 'read')
     read(unit=un, nml = probin)
     close(unit=un)
     need_inputs = .false.
  end if

  farg = 1
  if ( need_inputs .AND. narg >= 1 ) then
     call get_command_argument(farg, value = fname)
     inquire(file = fname, exist = lexist )
     if ( lexist ) then
        farg = farg + 1
        un = unit_new()
        open(unit=un, file = fname, status = 'old', action = 'read')
        read(unit=un, nml = probin)
        close(unit=un)
        need_inputs = .false.
     end if
  end if

  inquire(file = 'inputs_varden', exist = lexist)
  if ( need_inputs .AND. lexist ) then
     un = unit_new()
     open(unit=un, file = 'inputs_varden', status = 'old', action = 'read')
     read(unit=un, nml = probin)
     close(unit=un)
     need_inputs = .false.
  end if

  if ( .true. ) then
     do while ( farg <= narg )
        call get_command_argument(farg, value = fname)
        select case (fname)

        case ('--model_file')
           farg = farg + 1
           call get_command_argument(farg, value = model_file)

        case ('--prob_lo_x')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) prob_lo_x

        case ('--prob_hi_x')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) prob_hi_x

        case ('--cfl')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) cflfac

        case ('--stop_time')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) stop_time

        case ('--max_step')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) max_step

        case ('--spherical_in')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) spherical_in

        case ('--n_base')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) n_base

        case ('--')
           farg = farg + 1
           exit

        case default
           if ( .not. parallel_q() ) then
              write(*,*) 'UNKNOWN option = ', fname
              call bl_error("MAIN")
           end if
        end select

        farg = farg + 1
     end do
  end if

  dm = 1

  call init_spherical(spherical_in)

  center(1) = ZERO


  call init_variables(dm, nscal, nspec)
  call network_init()


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! define the grid spacing on all levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  prob_lo(1) = prob_lo_x
  prob_hi(1) = prob_hi_x

  dx(1) = (prob_hi_x-prob_lo_x) / real(n_base)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! allocate storage for the base state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dr_base = (prob_hi_x-prob_lo_x) / dble(n_base)

  allocate(div_coeff_old(n_base))
  allocate(div_coeff_new(n_base))

  allocate(grav_cell(n_base))

  allocate(   gam1(n_base  ))
  allocate( s0_old(n_base, nscal))
  allocate( s0_new(n_base, nscal))
  allocate( s0_1  (n_base, nscal))
  allocate( s0_2  (n_base, nscal))
  allocate(  temp0(n_base  ))
  allocate( p0_old(n_base  ))
  allocate( p0_new(n_base  ))
  allocate( p0_1  (n_base  ))
  allocate( p0_2  (n_base  ))
  allocate(     w0(0:n_base))
  allocate(Sbar_in(n_base))

  w0(:) = ZERO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in the base state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call init_geometry(center,n_base,dr_base)
    call init_base_state(model_file,n_base,s0_old,temp0,p0_old,gam1,dx,prob_lo,prob_hi)


    do i = 1, n_base
       print *, i, z(i), s0_old(i,rho_comp), p0_old(i), gam1(i)
    enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute the heating term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  y_0 = 4.e7

  do i = 1, n_base

     Hbar = 1.e16 * exp(-((z(i) - y_0)**2)/ 1.e14)
     
     ! (rho, T) --> p,h, etc
     input_flag = 1
     den_row(1)  = s0_old(i,rho_comp)
     temp_row(1) = temp0(i)
     xn_zone(:) = s0_old(i,spec_comp:spec_comp-1+nspec)/s0_old(i,rho_comp)

     call eos(input_flag, den_row, temp_row, NP, nspec, &
              xn_zone, aion, zion, &
              p_row, h_row, e_row, &
              cv_row, cp_row, xne_row, eta_row, pele_row, &
              dpdt_row, dpdr_row, dedt_row, dedr_row, &
              dpdX_row, dhdX_row, &
              gam1_row, cs_row, s_row, &
              dsdt_row, dsdr_row, &
              do_diag)

     s0_old(i,rhoh_comp) = den_row(1)*h_row(1)
     p0_old(i) = p_row(1)
     gam1(i) = gam1_row(1)

     coeff = dpdt_row(1)/ (den_row(1) * cp_row(1) * dpdr_row(1))

     Sbar_in(i) = coeff*Hbar
  enddo




  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute w_0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




  deallocate(div_coeff_old,div_coeff_new,grav_cell)
  deallocate(gam1,s0_old,s0_new,s0_1,s0_2,temp0,p0_old,p0_new,w0)
  

end subroutine varden
