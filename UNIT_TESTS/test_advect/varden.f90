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
  use variables, only: nscal, init_variables, rho_comp, spec_comp
  use geometry, only:  nlevs, nlevs_radial, spherical, dm, &
                       dr_fine, nr_fine, &
                       init_dm, init_spherical, init_center, init_multilevel, init_radial, &
                       init_cutoff, destroy_geometry
  use network, only: network_init, nspec
  use eos_module, only: eos_init
!  use fill_3d_module
  use probin_module, only: itest_dir, &
                           prob_lo, prob_hi, pmask, drdxfac, &
                           use_eos_coulomb, &
                           test_set, &
                           ppm_type, &
                           edge_nodal_flag, &
                           probin_init, probin_close
  use initialize_module, only: initialize_bc, initialize_dx
  use bl_constants_module
  use multifab_physbc_module
  use multifab_fill_ghost_module
  use test_advect_module, only: init_density_3d

  implicit none

  real(dp_t) :: lenx,leny,lenz,max_dist

  integer :: i,n,comp
  integer :: ng_s

  type(ml_layout) :: mla

  type(multifab), allocatable :: sold(:), snew(:)
  type(multifab), allocatable :: umac(:,:)

  real(kind=dp_t), pointer :: sp(:,:,:,:)

  real(dp_t), allocatable :: rho0(:,:), rhoh0(:,:), p0(:,:), w0(:,:)

  real(kind=dp_t), pointer :: dx(:,:)

  integer, allocatable :: lo(:),hi(:)

  type(ml_boxarray) :: mba

  type(bc_tower) ::  the_bc_tower

  ! general Maestro initializations
  call probin_init()
  call init_dm()
  call init_spherical()
  call init_center()

  call init_variables()

  call network_init()
  call eos_init(use_eos_coulomb=use_eos_coulomb)

  ! setup the grid
  call read_a_hgproj_grid(mba, test_set)

  call ml_layout_build(mla,mba,pmask)

  ! check for proper nesting
  if (.not. ml_boxarray_properly_nested(mla%mba, 3, pmask)) then
     call bl_error('fixed_grids not properly nested')
  end if

  ! initialize nlevs
  nlevs = mla%nlevel
  nlevs_radial = merge(1, nlevs, spherical .eq. 1)

  ! initialize boundary conditions
  call initialize_bc(the_bc_tower,nlevs,pmask)
  do n = 1,nlevs
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  ! allocate states
  allocate(sold(nlevs),snew(nlevs),umac(nlevs,dm))

  if (ppm_type .eq. 2) then
     ng_s = 4
  else
     ng_s = 3
  end if

  ! sanity checks
  if (itest_dir > dm) then
     call bl_error("ERROR: itest_dir > dm in test_advect")
  endif

  if (spherical == 1) then
     call bl_error("ERROR: test_advect not defined for spherical = 1")
  endif


  ! build states
  do n = 1,nlevs
     call multifab_build(sold(n), mla%la(n), nscal, ng_s)
     call multifab_build(snew(n), mla%la(n), nscal, ng_s)
     do comp=1,dm
        call multifab_build(umac(n,comp), mla%la(n),1,1,nodal=edge_nodal_flag(comp,:))
     end do
  end do

  ! initialize_dx
  call initialize_dx(dx,mba,nlevs)

  ! now that we have dx we can initialize nr_fine and dr_fine
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
  call init_multilevel(sold)

  ! now that we have nr_fine and dr_fine we can create nr, dr, r_cc_loc, r_edge_loc
  call init_radial(nlevs,mba)

  ! allocate the cutoff coordinate arrays
  call init_cutoff(nlevs)


  ! allocate the base state and set it all to 0
  allocate( rho0(nlevs,0:nr_fine-1))
  allocate(rhoh0(nlevs,0:nr_fine-1))
  allocate(   p0(nlevs,0:nr_fine-1))
  allocate(   w0(nlevs,0:nr_fine))

   rho0(:,:) = ZERO
  rhoh0(:,:) = ZERO
     p0(:,:) = ZERO
     w0(:,:) = ZERO


  ! other allocations
  allocate(lo(dm))
  allocate(hi(dm))
  

  ! initialize the velocity field -- it is unity in the direction of propagation
  do n = 1, nlevs

     select case (itest_dir)

     case (1)
        call setval(umac(n,1), ONE,  all=.true.)
        call setval(umac(n,2), ZERO, all=.true.)
        call setval(umac(n,3), ZERO, all=.true.)
     
     case (2)
        call setval(umac(n,1), ZERO, all=.true.)
        call setval(umac(n,2), ONE,  all=.true.)
        call setval(umac(n,3), ZERO, all=.true.)

     case (3)
        call setval(umac(n,1), ZERO, all=.true.)
        call setval(umac(n,2), ZERO, all=.true.)
        call setval(umac(n,3), ONE,  all=.true.)

     end select

  enddo


  ! initialize the density field
  do n=1,nlevs
     do i = 1, sold(n)%nboxes
        if ( multifab_remote(sold(n),i) ) cycle
        sp => dataptr(sold(n), i)
        lo = lwb(get_box(sold(n), i))
        hi = upb(get_box(sold(n), i))
        
        select case (dm)
        case (2)

        case (3)
           call init_density_3d(sp(:,:,:,rho_comp), sp(:,:,:,spec_comp:spec_comp-1+nspec), &
                                sold(n)%ng, lo, hi, dx(n,:))
        end select
     end do
  end do


  ! compute the initial timestep


  ! advance the density using the constant velocity field



  ! clean-up
  do n = 1,nlevs
     call destroy(sold(n))
     call destroy(snew(n))
     do comp=1,dm
        call destroy(umac(n,comp))
     end do
  end do

  call destroy(mla)
  call destroy(mba)

  deallocate(sold,snew,umac)

  call bc_tower_destroy(the_bc_tower)

  call probin_close()

  call destroy_geometry()

end subroutine varden



