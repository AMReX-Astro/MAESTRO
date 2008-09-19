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
  use initialize_module

  implicit none

  real(dp_t) :: lenx,leny,lenz,max_dist
  integer    :: i,n

  real(dp_t)  , pointer :: dx(:,:)
  type(ml_layout)           :: mla

  type(multifab), allocatable ::       uold(:)
  type(multifab), allocatable ::       unew(:)
  type(multifab), allocatable ::       sold(:)
  type(multifab), allocatable ::       snew(:)
  type(multifab), allocatable ::     normal(:)

  real(kind=dp_t), pointer :: nrp(:,:,:,:)
  integer, allocatable      :: lo(:),hi(:)

  type(layout)    :: la
  type(ml_boxarray) :: mba

  real(dp_t), allocatable :: s0_old(:,:,:)
  real(dp_t), allocatable :: s0_avg(:,:,:)
  real(dp_t), allocatable :: p0_old(:,:)
  real(dp_t), allocatable :: w0(:,:)

  type(bc_tower) ::  the_bc_tower

  real(dp_t) :: dr_base

  call probin_init()
  call init_dm()
  call init_spherical()
  call init_center()

  call init_variables()

  call network_init()
  call eos_init(use_eos_coulomb=use_eos_coulomb)

  call read_a_hgproj_grid(mba, test_set)

  call ml_layout_build(mla,mba,pmask)

  ! check for proper nesting
  if (.not. ml_boxarray_properly_nested(mla%mba, 3, pmask)) then
     call bl_error('fixed_grids not properly nested')
  end if

  ! initialize nlevs
  nlevs = mla%nlevel

  ! initialize boundary conditions
  call initialize_bc(the_bc_tower,nlevs,pmask)
  do n = 1,nlevs
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  ! allocate states
  allocate(uold(nlevs),sold(nlevs))

  ! build states
  do n = 1,nlevs
     call multifab_build(      uold(n), mla%la(n),    dm, 3)
     call multifab_build(      sold(n), mla%la(n), nscal, 3)
     call setval( uold(n),0.0_dp_t, all=.true.)
     call setval( sold(n),0.0_dp_t, all=.true.)
  end do

  ! initialize_dx
  call initialize_dx(dx,mba,nlevs)

  ! now that we have dx we can initialize nr_fine and dr_fine
  if (spherical .eq. 1) then
     
     ! for spherical, we will now require that dr_fine = dx
     dr_fine = dx(1,nlevs) / dble(drdxfac)
     
     lenx = HALF * (prob_hi(1) - prob_lo(1))
     leny = HALF * (prob_hi(2) - prob_lo(2))
     lenz = HALF * (prob_hi(3) - prob_lo(3))
     
     max_dist = sqrt(lenx**2 + leny**2 + lenz**2)
     nr_fine = int(max_dist / dr_fine) + 1
     
  else
     
     nr_fine = extent(mla%mba%pd(nlevs),dm)
     dr_fine = (prob_hi(dm)-prob_lo(dm)) / dble(nr_fine)
     
  end if

  ! create numdisjointchunks, r_start_coord, r_end_coord
  call init_multilevel(sold)

  ! now that we have nr_fine and dr_fine we can create nr, dr, r_cc_loc, r_edge_loc
  call init_radial(nlevs,mba)

  if (spherical .eq. 1) then

     ! for spherical, we will now require that dr_base = dx 
     dr_base = dx(1,nlevs)

     lenx = HALF * (prob_hi(1) - prob_lo(1))
     leny = HALF * (prob_hi(2) - prob_lo(2))
     lenz = HALF * (prob_hi(3) - prob_lo(3))

     max_dist = sqrt(lenx**2 + leny**2 + lenz**2)
     nr_fine = int(max_dist / dr_base) + 1

     if ( parallel_IOProcessor() ) then
        print *,'DISTANCE FROM CENTER TO CORNER IS ',max_dist
        print *,'DR_BASE IS ',dr_base
        print *,'SETTING NR_FINE TO ',nr_fine
     end if
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
     call multifab_build(   unew(n), mla%la(n),    dm, 3)
     call multifab_build(   snew(n), mla%la(n), nscal, 3)
     call multifab_build( normal(n), mla%la(n),    dm, 1)

     call setval(  unew(n),ZERO, all=.true.)
     call setval(  snew(n),ZERO, all=.true.)
  end do

  la = mla%la(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now initialize the grid data, and do initial projection if restart < 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Initialize base state at finest level
  do n=1,nlevs
     call init_base_state(n,model_file,s0_old(n,:,:),p0_old(n,:),dx(n,:))
  enddo

  allocate(lo(dm))
  allocate(hi(dm))

  ! Create the normal array once we have defined "center"
  if (spherical .eq. 1) then
     do n = 1,nlevs
        do i = 1, normal(n)%nboxes
           if ( multifab_remote(normal(n), i) ) cycle
           nrp => dataptr(normal(n), i)
           lo =  lwb(get_box(normal(n), i))
           hi =  upb(get_box(normal(n), i))
           call make_normal_3d_sphr(nrp(:,:,:,:),lo,hi,dx(n,:),1)
        end do
     end do
  end if

  call initveldata(uold,s0_old,p0_old,dx,the_bc_tower%bc_tower_array,mla)
  call initscalardata(sold,s0_old,p0_old,dx,the_bc_tower%bc_tower_array,mla)

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

  call probin_close()

  deallocate(dr,r_cc_loc,r_edge_loc,r_start_coord,r_end_coord,nr,numdisjointchunks,dx)

end subroutine varden
