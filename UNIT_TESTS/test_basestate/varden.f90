subroutine varden()

  use BoxLib
  use omp_module
  use f2kcli
  use list_box_module
  use ml_boxarray_module
  use layout_module
  use multifab_module
  use init_module
  use base_state_module
  use box_util_module
  use bl_IO_module
  use variables
  use geometry
  use network
  use eos_module
  use make_w0_module
  use advect_base_module
  use make_grav_module
  use make_div_coeff_module

  implicit none

  integer :: i

  integer    :: narg, farg
  integer    :: max_step
  integer    :: dm,n_base

  real(dp_t) :: dr_base(1)
  real(dp_t) :: cflfac
  real(dp_t) :: stop_time
  real(dp_t) :: time,dt,half_time,dtold
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

  real(dp_t) :: dx(1,1)
  real(dp_t) :: prob_lo(1), prob_hi(1)

  real(dp_t), allocatable :: div_coeff_old(:,:)
  real(dp_t), allocatable :: div_coeff(:,:)
  real(dp_t), allocatable :: grav_cell(:,:)
  real(dp_t), allocatable :: gam1(:,:)
  real(dp_t), allocatable :: s0_old(:,:,:)
  real(dp_t), allocatable :: s0(:,:,:)
  real(dp_t), allocatable :: p0_old(:,:)
  real(dp_t), allocatable :: p0(:,:)
  real(dp_t), allocatable :: w0(:,:)
  real(dp_t), allocatable :: w0_old(:,:)
  real(dp_t), allocatable :: eta(:,:,:)
  real(dp_t), allocatable :: f(:,:)
  real(dp_t), allocatable :: Sbar_in(:,:,:)

  real(dp_t) :: coeff, Hbar

  integer, parameter :: verbose = 0

  integer :: nlevs,n

  namelist /probin/ model_file
  namelist /probin/ stop_time
  namelist /probin/ prob_lo_x
  namelist /probin/ prob_hi_x
  namelist /probin/ max_step
  namelist /probin/ cflfac
  namelist /probin/ anelastic_cutoff
  namelist /probin/ spherical_in
  namelist /probin/ n_base

  nlevs = 1

  narg = command_argument_count()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize the runtime parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Defaults
  model_file = "model.hse"

  spherical_in = 0
  max_step  = 1

  prob_lo_x = ZERO
  prob_hi_x = 5.e8_dp_t

  anelastic_cutoff = 3.e6

  ntrac = 1
  nscal = nspec + ntrac + 2

  n_base = 512
  cflfac = 0.5_dp_t


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


  call init_variables(dm, nspec)
  call network_init()
  call helmeos_init()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! define the grid spacing on all levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  prob_lo(1) = prob_lo_x
  prob_hi(1) = prob_hi_x



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! allocate storage for the base state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dr_base(1) = (prob_hi_x-prob_lo_x) / dble(n_base)
  dx(1,1) = dr_base(1)

  allocate(div_coeff_old(nlevs,0:n_base-1))
  allocate(    div_coeff(nlevs,0:n_base-1))

  allocate(grav_cell(nlevs,0:n_base-1))

  allocate(   gam1(nlevs,0:n_base-1  ))
  allocate( s0_old(nlevs,0:n_base-1, nscal))
  allocate(     s0(nlevs,0:n_base-1, nscal))
  allocate( p0_old(nlevs,0:n_base-1  ))
  allocate(     p0(nlevs,0:n_base-1  ))
  allocate( w0_old(nlevs,0:n_base))
  allocate(     w0(nlevs,0:n_base))
  allocate(    eta(nlevs,0:n_base,   nscal))
  allocate(      f(nlevs,0:n_base))
  allocate(Sbar_in(nlevs,0:n_base-1,1))

  w0(:,:) = ZERO
  eta(:,:,:) = ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! read in the base state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call init_geometry(center,n_base,dr_base(1))

  do n = 1,nlevs
     call init_base_state(n,model_file,s0_old(n,:,:),p0_old(n,:),gam1(n,:),dx(n,:), &
                          prob_lo,prob_hi)
  enddo

  ! output
  open(unit=10,file="base.orig")
  do i = 0, n_base-1
     write(10,1000) z(1,i), s0(1,i,rho_comp), s0(1,i,temp_comp), p0(1,i)
  enddo
  close(unit=10)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! main timestepping loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
  time = ZERO
  dt = 0.0001_dp_t
  dtold = dt

  w0_old(:,:) = ZERO

  do while (time < stop_time)

     print *, 'time = ', time


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute the heating term
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     y_0 = 4.e7
       
     do i = 0, n_base-1

        Hbar = 1.e16 * exp(-((z(1,i) - y_0)**2)/ 1.e14)
     
        ! (rho, T) --> p,h, etc
        den_eos(1)  = s0(1,i,rho_comp)
        temp_eos(1) = s0(1,i,temp_comp)
        xn_eos(1,:) = s0(1,i,spec_comp:spec_comp-1+nspec)/s0(1,i,rho_comp)
        
        call eos(eos_input_rt, den_eos, temp_eos, NP, nspec, &
                 xn_eos, &
                 p_eos, h_eos, e_eos, &
                 cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                 dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                 dpdX_eos, dhdX_eos, &
                 gam1_eos, cs_eos, s_eos, &
                 dsdt_eos, dsdr_eos, &
                 do_diag)

        ! in the real Maestro code, we are updating the enthalpy by differencing
        ! the enthalpy equation with the heating term, rather than using the EOS.
        s0(1,i,rhoh_comp) = den_eos(1)*h_eos(1)

        p0(1,i) = p_eos(1)
        gam1(1,i) = gam1_eos(1)

        coeff = dpdt_eos(1)/ (den_eos(1) * cp_eos(1) * dpdr_eos(1))

        Sbar_in(1,i,1) = coeff*Hbar
     enddo


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute w_0
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     w0(:,:) = ZERO

     call make_w0(nlevs,w0,w0_old,f,Sbar_in(:,:,1),p0,s0(:,:,rho_comp), &
                  gam1,eta,dt,dtold,verbose)
  

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute the divergance coefficient (beta_0)
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     do n=1,nlevs
        call make_grav_cell(n,grav_cell(n,:),s0(n,:,rho_comp))
        call make_div_coeff(n,div_coeff(n,:),s0(n,:,rho_comp),p0(n,:), &
                            gam1(n,:),grav_cell(n,:),anelastic_cutoff)     
     enddo




     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute the new base state
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! here we set the old state to the current state -- advect_base will 
     ! update the state
     s0_old(:,:,:) = s0(:,:,:)
     p0_old(:,:) = p0(:,:)

     call advect_base(nlevs,w0,Sbar_in,p0_old,p0, &
                      s0_old,s0, &
                      gam1,div_coeff,eta, &
                      dr_base,dt,anelastic_cutoff)


     print *, 'new base pressure', p0(1,1)
     print *, 'new base density', s0(1,1,rho_comp)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! compute the new timestep
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     time = time + dt
     dtold = dt

     dt = min(1.1*dt,cflfac*dr_base(1)/maxval(abs(w0)))
     if (time+dt > stop_time) dt = stop_time - time

     ! store the old velocity
     w0_old(:,:) = w0(:,:)

  enddo

  ! output
  open(unit=10,file="base.new")
  do i = 0, n_base-1
     write(10,1000) z(1,i), s0(1,i,rho_comp), s0(1,i,temp_comp), p0(1,i)
  enddo
  close(unit=10)
1000 format(1x,5(g20.10))

  deallocate(div_coeff_old,div_coeff,grav_cell)
  deallocate(gam1,s0_old,s0,p0_old,p0,w0,f)
  

end subroutine varden
