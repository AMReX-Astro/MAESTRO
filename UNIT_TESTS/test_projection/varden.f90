! initialize a divergence-free velocity field, add the gradient of a
! scalar, project, and check if we get the original field back.

subroutine varden()

  use variables
  use fabio_module
  use multifab_module, only: multifab, build, destroy
  use ml_layout_module, only: ml_layout
  use ml_boxarray_module, only: ml_boxarray
  use runtime_init_module, only: runtime_init, runtime_close
  use probin_module, only: test_set, pmask, max_levs, &
       prob_lo, prob_hi, drdxfac
  use geometry, only: spherical, init_spherical, &
       init_center, destroy_geometry, &
       dr_fine, nr_fine, &
       init_cutoff, init_radial, init_multilevel
  use box_util_module, only: read_a_hgproj_grid
  use initialize_module, only: initialize_bc, initialize_dx
  use bl_IO_module
  use bl_constants_module
  use define_bc_module, only: bc_tower, bc_tower_level_build, bc_tower_destroy

  implicit none

  integer :: i, n, nlevs, dm
  integer :: ng_s

  type(ml_boxarray) :: mba
  type(ml_layout)   :: mla
  type(bc_tower)    :: the_bc_tower

  type(multifab), allocatable :: uold(:)
  type(multifab), allocatable :: unew(:)

  real(kind=dp_t), pointer :: dx(:,:)

  real(kind=dp_t) :: lenx,leny,lenz,max_dist

  real(kind=dp_t) :: time

  real(dp_t), parameter :: SMALL = 1.d-13

  character (len=20), allocatable :: plot_names(:)


  !---------------------------------------------------------------------------
  ! initialize some things
  !---------------------------------------------------------------------------

  ! this reads in the inputs file and sets the values of the parameters
  ! in the probin_module
  call runtime_init()

  ! set spherical in the geometry module
  call init_spherical()     

  ! set the center() array in geometry with the coordinates
  ! of the center of the domain
  call init_center()

  ! initialize the integer keys used for variable names
  call init_variables()


  time = ZERO


  !---------------------------------------------------------------------------
  ! setup the grid
  !---------------------------------------------------------------------------

  ! we only use fixed grids
  if (test_set /= '') then
     call read_a_hgproj_grid(mba, test_set)
     call ml_layout_build(mla,mba,pmask)
  else
     call bl_error('test_projection only implemented with fixed grids!')
  endif


  ! check for proper nesting
  if (.not. ml_boxarray_properly_nested(mla%mba, 3, pmask)) then
     call bl_error('fixed_grids not properly nested')
  end if

  
  dm = mla%dim
  nlevs = mla%nlevel


  ! initialize dx
  call initialize_dx(dx,mba,nlevs)


  ! initialize boundary conditions
  call initialize_bc(the_bc_tower,nlevs,pmask)
  do n = 1,nlevs
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do


  !---------------------------------------------------------------------------
  ! sanity checks
  !---------------------------------------------------------------------------
  if (nlevs .ne. max_levs) then
     call bl_error('varden.f90: nlevs .ne. max_levs not supported yet')
  end if

  if (dm .ne. get_dim(mla%mba) .or. dm .ne. 2) then
     call bl_error('dm_in not properly set in inputs file')
  endif
  
  if (spherical .eq. 1 ) then
     call bl_error('spherical = 1 not supported')
  endif
  
  if (abs(dx(1,1) - dx(1,2)) > SMALL) then
     call bl_error('zones must be square')
  endif


  !---------------------------------------------------------------------------
  ! allocate arrays
  !---------------------------------------------------------------------------
  allocate(uold(nlevs), unew(nlevs))

  ng_s = 4

  do n = 1, nlevs
     call build(uold(n),  mla%la(n), dm, ng_s)
     call setval(uold(n), ZERO, all=.true.)
     
     call build(unew(n),  mla%la(n), dm, ng_s)
     call setval(unew(n), ZERO, all=.true.)
  enddo
  

  ! some stats
  if (parallel_IOProcessor()) then
     print *, 'number of processors = ', parallel_nprocs()
     print *, 'number of dimensions = ', dm
     do n = 1, nlevs
        print *, 'level: ', n
        print *, '   number of boxes = ', uold(n)%nboxes
        print *, '   maximum zones   = ', (extent(mla%mba%pd(n),i),i=1,dm)
     end do
     print *, ''
  end if


  !---------------------------------------------------------------------------
  ! some multilevel base state initialization
  !---------------------------------------------------------------------------

  ! setup some base state coordinate information: nr_fine and dr_fine in
  ! geometry
  if (spherical .eq. 1) then

     ! for spherical, we will now require that dr_fine = dx
     dr_fine = dx(nlevs,1) / dble(drdxfac)

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
  call init_multilevel(uold)


  ! now that we have nr_fine and dr_fine we can create nr, dr,
  ! r_cc_loc, r_edge_loc
  call init_radial(nlevs,mba)


  ! allocate the cutoff coordinate arrays
  call init_cutoff(nlevs)


  !---------------------------------------------------------------------------
  ! initialize velocity field
  !---------------------------------------------------------------------------



  ! fill ghostcells


  !---------------------------------------------------------------------------
  ! dump out the initial data
  !---------------------------------------------------------------------------

  ! plotnames
  allocate(plot_names(dm))
  plot_names(1) = "x-velocity"
  if (dm >= 2) plot_names(2) = "y-velocity"
  if (dm == 3) plot_names(3) = "z-velocity"

  call fabio_ml_write(uold, mla%mba%rr(:,1), "u_init", names=plot_names)


  !---------------------------------------------------------------------------
  ! clean-up
  !---------------------------------------------------------------------------

  do n = 1, nlevs
     call destroy(uold(n))
     call destroy(unew(n))
  enddo

  call destroy(mla)
  call destroy(mba)
  
  deallocate(uold, unew)
  deallocate(plot_names)
  
  call bc_tower_destroy(the_bc_tower)
  
  call runtime_close()
  
  call destroy_geometry()

end subroutine varden



  

