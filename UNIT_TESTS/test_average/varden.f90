subroutine varden()

  use BoxLib
  use f2kcli
  use list_box_module
  use ml_boxarray_module
  use layout_module
  use multifab_module
  use init_module
  use base_state_module
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
  use eos_module
  use fill_3d_module
  use probin_module
  use bl_constants_module

  implicit none

  integer    :: dm
  integer    :: nr_fine
  real(dp_t) :: lenx,leny,lenz,max_dist
  integer    :: k,ng_cell
  integer    :: i, j, d, n, nlevs
  integer    :: comp,bc_comp

  integer     , allocatable :: domain_phys_bc(:,:)

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

  real(dp_t), allocatable :: s0_old(:,:,:)
  real(dp_t), allocatable :: s0_avg(:,:,:)
  real(dp_t), allocatable :: p0_old(:,:)
  real(dp_t), allocatable :: w0(:,:)

  type(bc_tower) ::  the_bc_tower

  type(bc_level) ::  bc

  ng_cell = 3

  call probin_init()

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
  call init_variables(dm, nspec)
  call network_init()
  call eos_init(use_eos_coulomb=use_eos_coulomb)

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
        nr_fine = int(max_dist / dr_base) + 1
        if ( parallel_IOProcessor() ) then
           print *,'DISTANCE FROM CENTER TO CORNER IS ',max_dist
           print *,'DR_BASE IS ',dr_base
           print *,'SETTING NR_FINE TO ',nr_fine
        end if
     else
        if ( parallel_IOProcessor() ) &
             print *,'NEED TO DEFINE DR_BASE '
        stop
     endif
  else
     ! NOTE: WE ASSUME DR_BASE IS THE RESOLUTION OF THE FINEST LEVEL IN PLANE-PARALLEL!
     nr_fine = extent(mba%pd(nlevs),dm)
     dr_base = (prob_hi(dm)-prob_lo(dm)) / dble(nr_fine)
  end if

  allocate( s0_old(nlevs,0:nr_fine-1,nscal))
  allocate( s0_avg(nlevs,0:nr_fine-1,nscal))
  allocate( p0_old(nlevs,0:nr_fine-1))
  allocate(     w0(nlevs,0:nr_fine  ))

  s0_old(:,:,:) = ZERO
  s0_avg(:,:,:) = ZERO
  w0(:,:) = ZERO

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
  call bc_tower_build(the_bc_tower,mla,domain_phys_bc,domain_boxes,nspec)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now initialize the grid data, and do initial projection if restart < 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Initialize geometry (IMPT: dr is set in init_base_state)
  center(1:dm) = HALF * (prob_lo(1:dm) + prob_hi(1:dm))
  call init_geometry(center,dr_base,nlevs,mla)

  ! Initialize base state at finest level
  do n=1,nlevs
     call init_base_state(n,model_file,s0_old(n,:,:),p0_old(n,:),dx(n,:))
  enddo

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

  call initveldata(nlevs,uold,s0_old,p0_old,dx,the_bc_tower%bc_tower_array,mla)
  call initscalardata(nlevs,sold,s0_old,p0_old,dx,the_bc_tower%bc_tower_array,mla)

  do n = 1,nlevs
     ! This is done to impose any Dirichlet bc's on unew or snew.
     call multifab_copy_c(unew(n),1,uold(n),1,dm   ,ng=unew(n)%ng)
     call multifab_copy_c(snew(n),1,sold(n),1,nscal,ng=snew(n)%ng)
  end do

  ! now that we are initialized, try averaging the state to 1-d
  ! and compare to the base state
  if ( parallel_IOProcessor() ) &
       print *, 'averaging...'

  do i=1,nscal
     call average(mla,sold,s0_avg(:,:,i),dx,i)
  end do

  if ( parallel_IOProcessor() ) &
       print *, 'done'

  ! compute the error against the base state
  if ( parallel_IOProcessor() ) then
     open (unit=10, file="dens.error")
     do n=1,nlevs
        do i = 0, nr(n)-1
           write (10,*) r_cc_loc(n,i), s0_old(n,i,rho_comp), s0_avg(n,i,rho_comp)
        enddo
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
  deallocate(s0_old,s0_avg,p0_old,w0)

  deallocate(lo,hi)

  call bc_tower_destroy(the_bc_tower)

  call destroy_geometry()

end subroutine varden
