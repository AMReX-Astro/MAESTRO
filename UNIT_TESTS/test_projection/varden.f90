! initialize a divergence-free velocity field, add the gradient of a
! scalar, project, and check if we get the original field back.

subroutine varden()

  use variables
  use fabio_module
  use multifab_module, only: multifab, build, destroy
  use ml_layout_module, only: ml_layout
  use ml_boxarray_module, only: ml_boxarray
  use runtime_init_module, only: runtime_init, runtime_close
  use probin_module, only: test_set, pmask, max_levs
  use geometry, only: spherical, init_spherical, init_center, destroy_geometry
  use box_util_module, only: read_a_hgproj_grid
  use initialize_module, only: initialize_bc, initialize_dx
  use bl_IO_module
  use bl_constants_module
  use define_bc_module, only: bc_tower, bc_tower_level_build, bc_tower_destroy

  implicit none

  integer :: i, n, r, comp, nlevs, dm
  integer :: ng_s

  type(ml_boxarray) :: mba
  type(ml_layout)   :: mla
  type(bc_tower)    :: the_bc_tower

  type(multifab), allocatable :: uold(:)
  type(multifab), allocatable :: unew(:)

  real(kind=dp_t), pointer :: dx(:,:)

  real(kind=dp_t) :: time

  real(dp_t), parameter :: SMALL = 1.d-13

  character (len=20), allocatable :: plot_names(:)



  ! initialize some things
  call runtime_init()
  call init_spherical()
  call init_center()

  call init_variables()


  time = ZERO


 ! we only use fixed grids
 if (test_set /= '') then
    call read_a_hgproj_grid(mba, test_set)
    call ml_layout_build(mla,mba,pmask)
 else
    call bl_error('test_projection only implemented with fixed grids!')
 endif

 dm = mla%dim
 nlevs = mla%nlevel


 ! initialize dx
 call initialize_dx(dx,mba,nlevs)


 ! sanity checks
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



 ! initialize boundary conditions
 call initialize_bc(the_bc_tower,nlevs,pmask)
 do n = 1,nlevs
    call bc_tower_level_build(the_bc_tower,n,mla%la(n))
 end do


 ! allocate arrays
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


 ! plotnames
 allocate(plot_names(dm))
 plot_names(1) = "x-velocity"
 if (dm >= 2) plot_names(2) = "y-velocity"
 if (dm == 3) plot_names(3) = "z-velocity"


 ! initialize velocity


 ! fill ghostcells


 ! dump the initial data
 call fabio_ml_write(uold, mla%mba%rr(:,1), "u_init", &
                     names=plot_names)


 ! clean-up
 do n = 1, nlevs
    call destroy(uold(n))
    call destroy(unew(n))
 enddo

 call destroy(mla)
 call destroy(mba)

 deallocate(uold, unew)

 call bc_tower_destroy(the_bc_tower)

 call runtime_close()

 call destroy_geometry()

end subroutine varden



  

