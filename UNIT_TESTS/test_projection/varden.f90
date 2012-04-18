! initialize a divergence-free velocity field, add the gradient of a
! scalar, project, and check if we get the original field back.

subroutine varden()

  use variables
  use fabio_module
  use multifab_module, only: multifab, build, multifab_build_edge, destroy
  use ml_layout_module, only: ml_layout
  use ml_boxarray_module, only: ml_boxarray
  use runtime_init_module, only: runtime_init, runtime_close
  use probin_module, only: test_set, pmask, max_levs, &
       prob_lo, prob_hi, drdxfac, nodal, run_prefix, project_type
  use geometry, only: spherical, init_spherical, &
       init_center, destroy_geometry, &
       dr_fine, nr_fine, &
       init_cutoff, init_radial, init_multilevel
  use box_util_module, only: read_a_hgproj_grid
  use initialize_module, only: initialize_bc, initialize_dx
  use bl_IO_module
  use bl_constants_module
  use define_bc_module, only: bc_tower, bc_tower_level_build, bc_tower_destroy
  use test_projection_module
  use proj_parameters
  use hgproject_module
  use macproject_module

  implicit none

  integer :: i, n, nlevs, dm, comp
  integer :: ng_s

  type(ml_boxarray) :: mba
  type(ml_layout)   :: mla
  type(bc_tower)    :: the_bc_tower

  type(multifab), allocatable :: uold(:)
  type(multifab), allocatable :: umid(:)
  type(multifab), allocatable :: unew(:)
  type(multifab), allocatable :: gphi(:)

  type(multifab), allocatable :: umac_old(:,:)
  type(multifab), allocatable :: umac_mid(:,:)
  type(multifab), allocatable :: umac_new(:,:)
  type(multifab), allocatable :: gphi_mac(:,:)
  type(multifab), allocatable :: utemp(:)

  type(multifab), allocatable :: rhohalf(:)
  type(multifab), allocatable :: div_coeff(:)
  type(multifab), allocatable :: pi(:)
  type(multifab), allocatable :: macpi(:)
  type(multifab), allocatable :: gpi(:)

  real(kind=dp_t), pointer :: dx(:,:)

  real(kind=dp_t) :: lenx,leny,lenz,max_dist

  real(kind=dp_t) :: time, dt

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
  ng_s = 4

  if (project_type == 1) then
  
     ! HG projection.  Velocities are cell-centered

     allocate(uold(nlevs), umid(nlevs), unew(nlevs), gphi(nlevs))

     do n = 1, nlevs
        call build(uold(n),  mla%la(n), dm, ng_s)
        call setval(uold(n), ZERO, all=.true.)
        
        call build(umid(n),  mla%la(n), dm, ng_s)
        call setval(umid(n), ZERO, all=.true.)
        
        call build(unew(n),  mla%la(n), dm, ng_s)
        call setval(unew(n), ZERO, all=.true.)
        
        call build(gphi(n),  mla%la(n), dm, ng_s)
        call setval(gphi(n), ZERO, all=.true.)
     enddo
  
  else if (project_type == 2) then

     ! MAC projection.  Velocities are nodal in respective dimension

     allocate(umac_old(nlevs,dm), umac_mid(nlevs,dm), umac_new(nlevs,dm), &
              gphi_mac(nlevs,dm))

     do n = 1, nlevs
        do comp = 1, dm
           call multifab_build_edge(umac_old(n,comp),  mla%la(n), 1, ng_s, comp)
           call setval(umac_old(n,comp), ZERO, all=.true.)
        
           call multifab_build_edge(umac_mid(n,comp),  mla%la(n), 1, ng_s, comp)
           call setval(umac_mid(n,comp), ZERO, all=.true.)
        
           call multifab_build_edge(umac_new(n,comp),  mla%la(n), 1, ng_s, comp)
           call setval(umac_new(n,comp), ZERO, all=.true.)
        
           call multifab_build_edge(gphi_mac(n,comp),  mla%la(n), 1, ng_s, comp)
           call setval(gphi_mac(n,comp), ZERO, all=.true.)
        enddo
     enddo     

     ! some initialization stuff wants a cell-centered multifab, also,
     ! we will need a container to map the MAC velocity to a
     ! cell-centered velocity for output.
     allocate(utemp(nlevs))

     do n = 1, nlevs
        call build(utemp(n),  mla%la(n), dm, ng_s)
        call setval(utemp(n), ZERO, all=.true.)
     enddo

  else

     call bl_error("ERROR: invalid project_type")
     
  endif

  ! some stats
  if (parallel_IOProcessor()) then
     print *, 'number of processors = ', parallel_nprocs()
     print *, 'number of dimensions = ', dm
     do n = 1, nlevs
        print *, 'level: ', n
        print *, '   number of boxes = ', mla%la(n)%lap%nboxes
        print *, '   maximum zones   = ', (extent(mla%mba%pd(n),i),i=1,dm)
     end do
     print *, ''
  end if


  !------------------------------------------------------------------------
  ! some multilevel base state initialization
  !------------------------------------------------------------------------
  
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
  if (project_type == 1) then
     call init_multilevel(uold)
  else
     call init_multilevel(utemp)
  endif
  
  
  ! now that we have nr_fine and dr_fine we can create nr, dr,
  ! r_cc_loc, r_edge_loc
  call init_radial(nlevs,mba)
  
  
  ! allocate the cutoff coordinate arrays
  call init_cutoff(nlevs)
  
  
  !---------------------------------------------------------------------------
  ! initialize velocity field
  !---------------------------------------------------------------------------
  if (project_type == 1) then
     call init_velocity(uold, dx, mla, the_bc_tower%bc_tower_array)
  else
     call init_mac_velocity(umac_old, dx, mla, the_bc_tower%bc_tower_array)     
  endif

  ! output the initial velocity field
  allocate(plot_names(dm))
  plot_names(1) = "x-velocity"
  if (dm >= 2) plot_names(2) = "y-velocity"
  if (dm == 3) plot_names(3) = "z-velocity"

  
  if (project_type == 1) then

     call fabio_ml_write(uold, mla%mba%rr(:,1), trim(run_prefix) // "u_init", &
                         names=plot_names)

     ! copy the velocity field over to the intermediate state, umid
     do n = 1, nlevs
        call multifab_copy(umid(n), uold(n), nghost(uold(n)))
     enddo

  else
  
     ! cannot write out a MAC field -- convert to cell-centered
     call convert_MAC_to_cc(umac_old, utemp)

     call fabio_ml_write(utemp, mla%mba%rr(:,1), trim(run_prefix) // "u_init", &
                         names=plot_names)

     ! copy the velocity field over to the intermediate state, umid
     do n = 1, nlevs
        do comp = 1, dm
           call multifab_copy_c(umac_mid(n,comp), 1, umac_old(n,comp), 1, 1, &
                                nghost(umac_old(n,comp)))
        enddo
     enddo

  end if


  !---------------------------------------------------------------------------
  ! 'pollute' the velocity field by adding the gradient of a scalar
  !---------------------------------------------------------------------------

  if (project_type == 1) then
     call add_grad_scalar(umid, gphi, dx, mla, the_bc_tower%bc_tower_array)

     call fabio_ml_write(umid, mla%mba%rr(:,1), trim(run_prefix) // "u_plus_grad_phi", &
                         names=plot_names)

     call fabio_ml_write(gphi, mla%mba%rr(:,1), trim(run_prefix) // "grad_phi", &
                         names=plot_names)

     ! copy the velocity field over to the final state, unew
     do n = 1, nlevs
        call multifab_copy(unew(n), umid(n), nghost(uold(n)))
     enddo

  else
     call add_grad_scalar_mac(umac_mid, gphi_mac, dx, mla, the_bc_tower%bc_tower_array)

     ! cannot write out a MAC field -- convert to cell-centered
     call convert_MAC_to_cc(umac_mid, utemp)

     call fabio_ml_write(utemp, mla%mba%rr(:,1), trim(run_prefix) // "u_plus_grad_phi", &
                         names=plot_names)

     call convert_MAC_to_cc(gphi_mac, utemp)

     call fabio_ml_write(utemp, mla%mba%rr(:,1), trim(run_prefix) // "grad_phi", &
                         names=plot_names)

     ! copy the velocity field over to the final state, unew
     do n = 1, nlevs
        do comp = 1, dm
           call multifab_copy_c(umac_new(n,comp), 1, umac_mid(n,comp), 1, 1, &
                                nghost(umac_mid(n,comp)))
        enddo
     enddo
  endif
     

  !---------------------------------------------------------------------------
  ! project out the divergent portion of the velocity field
  !---------------------------------------------------------------------------
  if (project_type == 1) then

     ! hgprojection -- here pi is nodal and u is cell-centered

     allocate(rhohalf(nlevs), pi(nlevs), gpi(nlevs), div_coeff(nlevs))

     do n = 1, nlevs

        ! build the density used in the projection -- we are just doing
        ! constant density, so set it to 1
        call build(rhohalf(n),  mla%la(n), 1, ng_s)
        call setval(rhohalf(n), ONE, all=.true.)
        
        ! build pi
        call build(pi(n),  mla%la(n), 1, 1, nodal)
        call setval(pi(n), ZERO, all=.true.)
        
        ! build gpi
        call build(gpi(n),  mla%la(n), dm, 1)
        call setval(gpi(n), ZERO, all=.true.)
        
        ! build the coefficient in the divergence.  We are doing 
        ! divergence-free (incompressible), so set div_coeff = 1
        call build(div_coeff(n),  mla%la(n), 1, 1)
        call setval(div_coeff(n), ONE, all=.true.)
        
     enddo

     ! hgproject takes dt -- it has no meaning for our projection type
     dt = ONE

     call hgproject(initial_projection_comp, &
                    mla, &
                    unew, unew, rhohalf, &
                    pi, gpi, &
                    dx, dt, the_bc_tower, &
                    div_coeff)
     
     call fabio_ml_write(unew, mla%mba%rr(:,1), trim(run_prefix) // "u_new", &
                         names=plot_names)

  else

     ! mac projection -- here pi is cell-centered and u is MAC
     allocate(rhohalf(nlevs), macpi(nlevs))

     do n = 1, nlevs

        ! build the density used in the projection -- we are just doing
        ! constant density, so set it to 1
        call build(rhohalf(n),  mla%la(n), 1, ng_s)
        call setval(rhohalf(n), ONE, all=.true.)
        
        ! build macpi
        call build(macpi(n),  mla%la(n), 1, ng_s)
        call setval(macpi(n), ZERO, all=.true.)
        
     enddo

     call macproject(mla, &
                     umac_new, macpi, rhohalf, &
                     dx, the_bc_tower)


     ! convert to cell-centered for output
     call convert_MAC_to_cc(umac_new, utemp)
     
     call fabio_ml_write(utemp, mla%mba%rr(:,1), trim(run_prefix) // "u_new", &
                         names=plot_names)

  endif


  !---------------------------------------------------------------------------
  ! clean-up
  !---------------------------------------------------------------------------

  do n = 1, nlevs
     call destroy(uold(n))
     call destroy(umid(n))
     call destroy(unew(n))
     call destroy(gphi(n))
  enddo

  call destroy(mla)
  call destroy(mba)
  
  deallocate(uold, umid, unew)
  deallocate(plot_names)
  
  call bc_tower_destroy(the_bc_tower)
  
  call runtime_close()
  
  call destroy_geometry()

end subroutine varden



  

