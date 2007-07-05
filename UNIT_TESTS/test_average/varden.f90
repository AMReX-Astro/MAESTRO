subroutine varden()

  use BoxLib
  use f2kcli
  use list_box_module
  use ml_boxarray_module
  use layout_module
  use multifab_module
  use init_module
  use ml_restriction_module
  use bc_module
  use define_bc_module
  use bl_mem_stat_module
  use bl_timer_module
  use box_util_module
  use bl_IO_module
  use fabio_module
  use setbc_module
  use variables
  use geometry
  use network
  use average_module

  implicit none

  integer    :: narg, farg
  integer    :: dm,n_base
  real(dp_t) :: dr_base
  real(dp_t) :: lenx,leny,lenz,max_dist
  integer    :: bcx_lo,bcx_hi,bcy_lo,bcy_hi,bcz_lo,bcz_hi
  integer    :: k,ng_cell,ng_edge,ng_fill
  integer    :: i, j, d, n, nlevs, nscal, ntrac
  integer    :: comp,bc_comp
  logical    :: pmask_x,pmask_y,pmask_z
  logical    :: perturb_model

  real(dp_t) :: prob_lo_x,prob_lo_y,prob_lo_z
  real(dp_t) :: prob_hi_x,prob_hi_y,prob_hi_z

  real(dp_t) :: anelastic_cutoff
  integer    :: spherical_in

  integer     , allocatable :: domain_phys_bc(:,:)
  logical     :: pmask_xyz(MAX_SPACEDIM)
  logical     , allocatable :: pmask(:)
  real(dp_t)  , allocatable :: dx(:,:)
  real(dp_t)  , allocatable :: prob_hi(:)
  real(dp_t)  , allocatable :: prob_lo(:)
  type(ml_layout)           :: mla
  type(box)   , allocatable :: domain_boxes(:)

  type(multifab), allocatable ::       uold(:)
  type(multifab), allocatable ::       unew(:)
  type(multifab), allocatable ::       sold(:)
  type(multifab), allocatable ::       snew(:)
  type(multifab), allocatable ::     normal(:)

  real(kind=dp_t), pointer :: uop(:,:,:,:)
  real(kind=dp_t), pointer :: sop(:,:,:,:)
  real(kind=dp_t), pointer :: nrp(:,:,:,:)
  integer,allocatable      :: lo(:),hi(:)

  character(len=128) :: fname
  character(len=128) :: probin_env
  character(len=128) :: test_set
  character(len=256) :: model_file
  character(len=20), allocatable :: plot_names(:)
  integer :: un, ierr
  logical :: lexist
  logical :: need_inputs
  logical, allocatable :: nodal(:)
  logical, allocatable :: umac_nodal_flag(:)

  logical :: init_mode

  type(layout)    :: la
  type(box)       :: fine_domain
  type(ml_boxarray) :: mba

  real(dp_t), allocatable :: gam1(:)
  real(dp_t), allocatable :: s0_old(:,:)
  real(dp_t), allocatable :: s0_avg(:,:)
  real(dp_t), allocatable :: temp0(:)
  real(dp_t), allocatable :: p0_old(:)
  real(dp_t), allocatable :: w0(:)


  type(bc_tower) ::  the_bc_tower

  type(bc_level) ::  bc

  namelist /probin/ model_file
  namelist /probin/ prob_lo_x
  namelist /probin/ prob_lo_y
  namelist /probin/ prob_lo_z
  namelist /probin/ prob_hi_x
  namelist /probin/ prob_hi_y
  namelist /probin/ prob_hi_z
  namelist /probin/ test_set
  namelist /probin/ bcx_lo
  namelist /probin/ bcx_hi
  namelist /probin/ bcy_lo
  namelist /probin/ bcy_hi
  namelist /probin/ bcz_lo
  namelist /probin/ bcz_hi
  namelist /probin/ pmask_x
  namelist /probin/ pmask_y
  namelist /probin/ pmask_z
  namelist /probin/ pmask_xyz
  namelist /probin/ anelastic_cutoff
  namelist /probin/ spherical_in
  namelist /probin/ dr_base

  ng_edge = 2
  ng_cell = 3

  narg = command_argument_count()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize the runtime parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Defaults
  model_file = "model.hse"

  spherical_in = 0
  dr_base = -1.d0
  prob_lo_x = ZERO
  prob_lo_y = ZERO
  prob_lo_z = ZERO

  anelastic_cutoff = 3.e6

  ntrac = 1
  nscal = nspec + ntrac + 2

  need_inputs = .true.
  test_set = ''
  
  bcx_lo = SLIP_WALL
  bcy_lo = SLIP_WALL
  bcz_lo = SLIP_WALL
  bcx_hi = SLIP_WALL
  bcy_hi = SLIP_WALL
  bcz_hi = SLIP_WALL

  pmask_x = .false.
  pmask_y = .false.
  pmask_z = .false.
  
  perturb_model = .false.

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

  pmask_xyz = (/pmask_x, pmask_y, pmask_z/)

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

        case ('--prob_lo_y')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) prob_lo_y

        case ('--prob_lo_z')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) prob_lo_z

        case ('--prob_hi_x')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) prob_hi_x

        case ('--prob_hi_y')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) prob_hi_y

        case ('--prob_hi_z')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) prob_hi_z

        case ('--bcx_lo')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcx_lo
        case ('--bcy_lo')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcy_lo
        case ('--bcz_lo')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcz_lo
        case ('--bcx_hi')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcx_hi
        case ('--bcy_hi')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcy_hi
        case ('--bcz_hi')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) bcz_hi

        case ('--pmask_x')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) pmask_xyz(1)
        case ('--pmask_y')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) pmask_xyz(2)
        case ('--pmask_z')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) pmask_xyz(3)
        case ('--pmask_xyz')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) pmask_xyz(1)
           pmask_xyz = pmask_xyz(1)

        case ('--test_set')
           farg = farg + 1
           call get_command_argument(farg, value = test_set)

        case ('--anelastic_cutoff')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) anelastic_cutoff

        case ('--spherical_in')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) spherical_in

        case ('--dr_base')
           farg = farg + 1
           call get_command_argument(farg, value = fname)
           read(fname, *) dr_base

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

  call init_spherical(spherical_in)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the arrays and read the restart data if restart >= 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call read_a_hgproj_grid(mba, test_set)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the variable index pointers and the reaction network
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dm = get_dim(mba)
  allocate(lo(dm),hi(dm))
  call init_variables(dm, nscal, nspec)
  call network_init()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! allocate storage for the state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(prob_lo(dm), prob_hi(dm), pmask(dm))
  pmask = pmask_xyz(1:dm)

  if ( parallel_IOProcessor() ) &
    print *, 'pmask = ', pmask

  nlevs = mba%nlevel
  call ml_layout_build(mla,mba,pmask)
  allocate(uold(nlevs),sold(nlevs))

  allocate(nodal(dm), umac_nodal_flag(dm))
  nodal = .true.

  do n = 1,nlevs
     call multifab_build(      uold(n), mla%la(n),    dm, ng_cell)
     call multifab_build(      sold(n), mla%la(n), nscal, ng_cell)
     call setval( uold(n),0.0_dp_t, all=.true.)
     call setval( sold(n),0.0_dp_t, all=.true.)
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! define the grid spacing on all levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(dx(nlevs,dm))

  prob_lo(1) = prob_lo_x
  if (dm > 1) prob_lo(2) = prob_lo_y
  if (dm > 2) prob_lo(3) = prob_lo_z
  prob_hi(1) = prob_hi_x
  if (dm > 1) prob_hi(2) = prob_hi_y
  if (dm > 2) prob_hi(3) = prob_hi_z

  do i = 1, dm
    dx(1,i) = (prob_hi(i)-prob_lo(i)) / real(extent(mba%pd(1),i),kind=dp_t)
  end do
  do n = 2,nlevs
    dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! allocate storage for the base state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (spherical .eq. 1) then
    if (dr_base .gt. 0) then
      lenx = HALF * (prob_hi_x - prob_lo_x)
      leny = HALF * (prob_hi_y - prob_lo_y)
      lenz = HALF * (prob_hi_z - prob_lo_z)
      max_dist = sqrt(lenx**2 + leny**2 + lenz**2)
      n_base = int(max_dist / dr_base) + 1
      if ( parallel_IOProcessor() ) then
         print *,'DISTANCE FROM CENTER TO CORNER IS ',max_dist
         print *,'DR_BASE IS ',dr_base
         print *,'SETTING N_BASE TO ',n_base
      end if
    else
     if ( parallel_IOProcessor() ) &
       print *,'NEED TO DEFINE DR_BASE '
      stop
    endif
  else
    ! NOTE: WE ASSUME DR_BASE IS THE RESOLUTION OF THE FINEST LEVEL IN PLANE-PARALLEL!
    n_base = extent(mba%pd(nlevs),dm)
    dr_base = (prob_hi(dm)-prob_lo(dm)) / dble(n_base)
  end if


  allocate(   gam1(n_base  ))
  allocate( s0_old(n_base, nscal))
  allocate( s0_avg(n_base, nscal))
  allocate(  temp0(n_base  ))
  allocate( p0_old(n_base  ))
  allocate(     w0(0:n_base))

  s0_old(:,:) = ZERO
  s0_avg(:,:) = ZERO
  w0(:) = ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize all remaining arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(unew(nlevs),snew(nlevs))
  allocate(normal(nlevs))


  do n = nlevs,1,-1
     call multifab_build(   unew(n), mla%la(n),    dm, ng_cell)
     call multifab_build(   snew(n), mla%la(n), nscal, ng_cell)
     call multifab_build(normal(n), mla%la(n), dm, 1)

     call setval(  unew(n),ZERO, all=.true.)
     call setval(  snew(n),ZERO, all=.true.)
  end do

  la = mla%la(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocate the arrays for the boundary conditions at the physical boundaries.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(domain_phys_bc(dm,2))

  allocate(domain_boxes(nlevs))
  do n = 1,nlevs
     domain_boxes(n) = layout_get_pd(mla%la(n))
  end do

! Put the bc values from the inputs file into domain_phys_bc
  domain_phys_bc(1,1) = bcx_lo
  domain_phys_bc(1,2) = bcx_hi
  if (dm > 1) then
     domain_phys_bc(2,1) = bcy_lo
     domain_phys_bc(2,2) = bcy_hi
  end if
  if (dm > 2) then
     domain_phys_bc(3,1) = bcz_lo
     domain_phys_bc(3,2) = bcz_hi
  end if

  do i = 1, dm
     if ( pmask(i) ) domain_phys_bc(i,:) = BC_PER
  end do

! Build the arrays for each grid from the domain_bc arrays.
  call bc_tower_build( the_bc_tower,mla,domain_phys_bc,domain_boxes,nscal,nspec,ntrac)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now initialize the grid data, and do initial projection if restart < 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Initialize geometry (IMPT: dr is set in init_base_state)
  center(1:dm) = HALF * (prob_lo(1:dm) + prob_hi(1:dm))
  call init_geometry(center,n_base,dr_base)

  ! Initialize base state at finest level
  call init_base_state(model_file,n_base,s0_old,temp0,p0_old,gam1,dx(nlevs,:),prob_lo,prob_hi)

  ! Create the normal array once we have defined "center"
  if (spherical .eq. 1) then
    do n = 1,nlevs
       do i = 1, normal(n)%nboxes
         if ( multifab_remote(normal(n), i) ) cycle
         nrp => dataptr(normal(n), i)
          lo =  lwb(get_box(normal(n), i))
          hi =  upb(get_box(normal(n), i))
         call make_3d_normal(nrp(:,:,:,:),lo,hi,dx(n,:),1)
        end do
    end do
  end if

  do n = 1,nlevs
     call initveldata(uold(n),s0_old,p0_old,temp0,dx(n,:), &
                      prob_lo,prob_hi, &
                      the_bc_tower%bc_tower_array(n),nscal,ntrac)

     call initscalardata(sold(n),s0_old,p0_old,temp0,dx(n,:), &
                         perturb_model, &
                         prob_lo,prob_hi, &
                         the_bc_tower%bc_tower_array(n),nscal,ntrac)

  end do

  do n = nlevs,2,-1
     call ml_cc_restriction(uold(n-1),uold(n),mba%rr(n-1,:))
     call ml_cc_restriction(sold(n-1),sold(n),mba%rr(n-1,:))
  end do




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set the boundary conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do n = 1,nlevs
     bc = the_bc_tower%bc_tower_array(n)
     do i = 1, uold(n)%nboxes
        if ( multifab_remote(uold(n), i) ) cycle
        uop => dataptr(uold(n), i)
        sop => dataptr(sold(n), i)
        lo =  lwb(get_box(uold(n), i))
        hi =  upb(get_box(uold(n), i))
        select case (dm)
        case (2)
           do d = 1,dm
              call setbc_2d(uop(:,:,1,d), lo, ng_cell, &
                            bc%adv_bc_level_array(i,:,:,d), &
                            dx(n,:),d)
           end do
           do d = 1,nscal
              call setbc_2d(sop(:,:,1,d), lo, ng_cell, &
                            bc%adv_bc_level_array(i,:,:,dm+d), &
                            dx(n,:),dm+d)
           end do
        case (3)
           do d = 1,dm
              call setbc_3d(uop(:,:,:,d), lo, ng_cell, &
                            bc%adv_bc_level_array(i,:,:,d), &
                            dx(n,:),d)
           end do
           do d = 1, nscal
              call setbc_3d(sop(:,:,:,d), lo, ng_cell, &
                            bc%adv_bc_level_array(i,:,:,dm+d), &
                            dx(n,:),dm+d)
           end do
        end select
     end do

     call multifab_fill_boundary(uold(n))
     call multifab_fill_boundary(sold(n))

!    This is done to impose any Dirichlet bc's on unew or snew.
     call multifab_copy_c(unew(n),1,uold(n),1,dm   ,all=.true.)
     call multifab_copy_c(snew(n),1,sold(n),1,nscal,all=.true.)

  end do


  ! now that we are initialized, try averaging the state to 1-d
  ! and compare to the base state
  if ( parallel_IOProcessor() ) &
       print *, 'averaging...'
  call average(sold,s0_avg(:,:),dx,1,nscal)
  if ( parallel_IOProcessor() ) &
       print *, 'done'

  ! compute the error against the base state
  if ( parallel_IOProcessor() ) then
     open (unit=10, file="dens.error")
     do n = 1, n_base
        write (10,*) n, s0_old(n,rho_comp), s0_avg(n,rho_comp)
     enddo
     close (10)
  endif






  do n = 1,nlevs
     call destroy(uold(n))
     call destroy(unew(n))
     call destroy(sold(n))
     call destroy(snew(n))
     call destroy(normal(n))
  end do

  call destroy(mla)
  call destroy(mba)

  deallocate(uold,unew,sold,snew)
  deallocate(gam1,s0_old,s0_avg,temp0,p0_old,w0)

  deallocate(lo,hi)

  call bc_tower_destroy(the_bc_tower)

  call destroy_geometry()

end subroutine varden
