! a driver that advects a Gaussian profile through a periodic domain
! in all possible coordinate directions.  This tests the advection
! scheme.  
!
! The norm of the error (initial profile compared to final profile)
! should be the same for all directions, otherwise there is some
! directional bias in the advection solver.  Note: this does not mean
! that the final density profile cannot have an imprint from the
! advection in a particular direction, just that advecting in any
! direction should show that same pattern of imprinting.
!
! Here we compute the norm level-by-level.  This driver assume that
! the grid is static, so the initial and final grids agree.  This
! allows us to compute a meaningful norm.

subroutine varden()

  use BoxLib
  use f2kcli
  use list_box_module
  use ml_boxarray_module
  use layout_module
  use multifab_module
  use base_state_module
  use ml_restriction_module
  use bc_module
  use define_bc_module
  use bl_mem_stat_module
  use bl_timer_module
  use box_util_module
  use bl_IO_module
  use fabio_module
  use variables, only: nscal, init_variables, rho_comp, spec_comp
  use fill_3d_module, only: make_normal
  use geometry, only:  nlevs_radial, spherical, &
                       dr_fine, nr_fine, &
                       init_spherical, init_center, init_multilevel, init_radial, &
                       init_cutoff, destroy_geometry
  use network, only: network_init, nspec
  use eos_module, only: eos_init
  use probin_module, only: dump_output, advect_test_tol, &
                           prob_lo, prob_hi, pmask, drdxfac, &
                           use_eos_coulomb, &
                           test_set, &
                           ppm_type, bds_type, &
                           cflfac, &
                           stop_time
  use runtime_init_module, only: runtime_init, runtime_close
  use initialize_module, only: initialize_bc, initialize_dx
  use bl_constants_module
  use multifab_physbc_module
  use multifab_fill_ghost_module
  use test_advect_module, only: init_density_2d, init_density_3d
  use density_advance_module, only: density_advance

  implicit none

  real(dp_t) :: lenx,leny,lenz,max_dist

  integer :: i,n,p,comp
  integer :: ng_s
  integer :: nlevs, dm

  integer :: idir, idim, itest_dir, index_t

  type(ml_layout) :: mla

  type(multifab), allocatable :: sold(:), snew(:)
  type(multifab), allocatable :: umac(:,:)
  type(multifab), allocatable :: normal(:)

  type(multifab), allocatable :: sedge(:,:), sflux(:,:)
  type(multifab), allocatable :: scal_force(:), etarhoflux(:)
  type(multifab), allocatable :: w0mac(:,:)

  type(multifab), allocatable :: dens_orig(:), dens_final(:)
  type(multifab), allocatable :: error(:)
  type(multifab), allocatable :: temporary(:)

  real(kind=dp_t), pointer :: sp(:,:,:,:)

  real(dp_t), allocatable :: rho0_old(:,:), rho0_new(:,:), rhoh0(:,:), p0(:,:), w0(:,:)
  real(dp_t), allocatable :: rho0_predicted_edge(:,:)

  real(kind=dp_t), pointer :: dx(:,:)

  integer, allocatable :: lo(:),hi(:)

  type(ml_boxarray) :: mba

  type(bc_tower) ::  the_bc_tower

  real(kind=dp_t) :: t, dt

  character (len=32) :: outname

  character (len=20) :: plot_names(3)

  real(kind=dp_t), allocatable :: abs_norm(:,:,:), rel_norm(:,:,:)
  real(kind=dp_t) :: max_dabs_error, max_abs_error
  real(kind=dp_t) :: max_drel_error, max_rel_error

  logical :: wrote_init_file = .false.

  ! general Maestro initializations
  call runtime_init()
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

!  if (nlevs .ne. max_levs) then
!     call bl_error('varden.f90: nlevs .ne. max_levs not supported yet')
!  end if

  nlevs_radial = merge(1, nlevs, spherical .eq. 1)

  ! initialize dm
  dm = mla%dim

  ! initialize boundary conditions
  call initialize_bc(the_bc_tower,nlevs,pmask)
  do n = 1,nlevs
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  ! allocate states
  allocate(sold(nlevs),snew(nlevs),umac(nlevs,dm))


  ! 4 ghostcells will work for all PPM types
  ng_s = 4


  if (spherical == 1) then
     call bl_error("ERROR: test_advect not defined for spherical = 1")
  endif


  ! build states
  do n = 1,nlevs
     call multifab_build(sold(n), mla%la(n), nscal, ng_s)
     call multifab_build(snew(n), mla%la(n), nscal, ng_s)
     do comp=1,dm
        call multifab_build_edge(umac(n,comp), mla%la(n),1,1,comp)
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


  ! allocate normal
  allocate (normal(nlevs))
  if (dm == 3) then
     do n = 1,nlevs
        call multifab_build(normal(n), mla%la(n),    dm, 1)
     enddo
  endif

  call make_normal(normal,dx)


  ! the initial and final densities
  allocate(dens_orig(nlevs))
  allocate(dens_final(nlevs))
  allocate(error(nlevs))
  allocate(temporary(nlevs))


  do n = 1,nlevs
     call multifab_build(dens_orig(n), mla%la(n), 1, ng_s)
     call multifab_build(dens_final(n), mla%la(n), 1, ng_s)
     call multifab_build(error(n), mla%la(n), 1, ng_s)
     call multifab_build(temporary(n), mla%la(n), 3, ng_s)
  end do

  plot_names(1) = "density"
  plot_names(2) = "rel. error"
  plot_names(3) = "abs. error"

  ! allocate the base state 
  allocate(           rho0_old(nlevs,0:nr_fine-1))
  allocate(           rho0_new(nlevs,0:nr_fine-1))
  allocate(              rhoh0(nlevs,0:nr_fine-1))
  allocate(                 p0(nlevs,0:nr_fine-1))
  allocate(                 w0(nlevs,0:nr_fine))
  allocate(rho0_predicted_edge(nlevs,0:nr_fine))


  ! other allocations
  allocate(lo(dm))
  allocate(hi(dm))
  
  allocate(sedge(nlevs,dm))
  allocate(sflux(nlevs,dm))

  allocate(scal_force(nlevs))
  allocate(etarhoflux(nlevs))

  allocate(w0mac(nlevs,dm))

  do n=1,nlevs
     do comp = 1,dm
        call multifab_build_edge(sedge(n,comp),mla%la(n),nscal,0,comp)
        call multifab_build_edge(sflux(n,comp),mla%la(n),nscal,0,comp)
        call multifab_build_edge(w0mac(n,comp),mla%la(n),1,1,comp)
        call setval(w0mac(n,comp),ZERO,all=.true.)
     end do

     call multifab_build(scal_force(n), mla%la(n), nscal, 1)
     call multifab_build_edge(etarhoflux(n), mla%la(n), 1, 0, dm)
     call setval(scal_force(n),ZERO,all=.true.)
     call setval(etarhoflux(n),ZERO,all=.true.)
  end do


  allocate(rel_norm(0:3,nlevs,2*dm))
  allocate(abs_norm(0:3,nlevs,2*dm))

  rel_norm = ZERO
  abs_norm = ZERO


  !---------------------------------------------------------------------------
  ! loop over all possible directions and all ppm_types
  !---------------------------------------------------------------------------
  do p = 0,3
     
     if (p <= 2) then
        ppm_type = p

        print *, ' '
        print *, '<<< ppm_type ', p, ' >>>'
        print *, ' '

     else
        bds_type = 1

        print *, ' '
        print *, '<<< bds_type ', 1, ' >>>'
        print *, ' '
     endif



     do idim = 1, dm
        do idir = -1, 1, 2

           itest_dir = idir * idim

           ! a unique index for indexing array by test
           index_t = idim
           if (idir < 0) then
              index_t = index_t + dm
           endif
        

           print *, ' '
           print *, '<<< advection test in direction ', itest_dir, ' >>>'
           print *, 'index_t = ', index_t

           ! initialize the velocity field -- it is unity in the
           ! direction of propagation a negative itest_dir indicates
           ! negative velocity
           do n = 1, nlevs

              select case (itest_dir)

              case (-1)
                 call setval(umac(n,1), -ONE,  all=.true.)
                 call setval(umac(n,2), ZERO, all=.true.)
                 if (dm == 3) call setval(umac(n,3), ZERO, all=.true.)

              case (1)
                 call setval(umac(n,1), ONE,  all=.true.)
                 call setval(umac(n,2), ZERO, all=.true.)
                 if (dm == 3) call setval(umac(n,3), ZERO, all=.true.)
     
              case (-2)
                 call setval(umac(n,1), ZERO, all=.true.)
                 call setval(umac(n,2), -ONE,  all=.true.)
                 if (dm == 3) call setval(umac(n,3), ZERO, all=.true.)

              case (2)
                 call setval(umac(n,1), ZERO, all=.true.)
                 call setval(umac(n,2), ONE,  all=.true.)
                 if (dm == 3) call setval(umac(n,3), ZERO, all=.true.)
                 
              case (-3)
                 call setval(umac(n,1), ZERO, all=.true.)
                 call setval(umac(n,2), ZERO, all=.true.)
                 call setval(umac(n,3), -ONE,  all=.true.)
                 
              case (3)
                 call setval(umac(n,1), ZERO, all=.true.)
                 call setval(umac(n,2), ZERO, all=.true.)
                 call setval(umac(n,3), ONE,  all=.true.)
                 
              end select

           enddo


           ! the base state will not carry any information in this test problem
           rho0_old(:,:) = ZERO
           rho0_new(:,:) = ZERO
           rhoh0(:,:) = ZERO
           p0(:,:) = ZERO
           w0(:,:) = ZERO
           rho0_predicted_edge(:,:) = ZERO


           ! initialize the density field and species
           do n=1,nlevs
              do i = 1, sold(n)%nboxes
                 if ( multifab_remote(sold(n),i) ) cycle
                 sp => dataptr(sold(n), i)
                 lo = lwb(get_box(sold(n), i))
                 hi = upb(get_box(sold(n), i))
                 
                 select case (dm)
                 case (2)
                    call init_density_2d(sp(:,:,1,rho_comp), &
                                         sp(:,:,1,spec_comp:spec_comp-1+nspec), &
                                         sold(n)%ng, lo, hi, dx(n,:))

                 case (3)
                    call init_density_3d(sp(:,:,:,rho_comp), &
                                         sp(:,:,:,spec_comp:spec_comp-1+nspec), &
                                         sold(n)%ng, lo, hi, dx(n,:))
                 end select
              end do
           end do


           ! ghost cell fill
           if (nlevs .eq. 1) then
           
              ! fill ghost cells for two adjacent grids at the same level
              ! this includes periodic domain boundary ghost cells
              call multifab_fill_boundary(sold(nlevs))
     
              ! fill non-periodic domain boundary ghost cells
              call multifab_physbc(sold(nlevs),rho_comp,dm+rho_comp,nscal, &
                                   the_bc_tower%bc_tower_array(nlevs))

           else

              ! the loop over nlevs must count backwards to make sure
              ! the finer grids are done first
              do n=nlevs,2,-1
        
                 ! set level n-1 data to be the average of the level n
                 ! data covering it
                 call ml_cc_restriction(sold(n-1),sold(n),mla%mba%rr(n-1,:))

                 ! fill level n ghost cells using interpolation from
                 ! level n-1 data note that multifab_fill_boundary and
                 ! multifab_physbc are called for both levels n-1 and
                 ! n
                 call multifab_fill_ghost_cells(sold(n),sold(n-1),sold(n)%ng, &
                                                mla%mba%rr(n-1,:), &
                                                the_bc_tower%bc_tower_array(n-1), &
                                                the_bc_tower%bc_tower_array(n), &
                                                rho_comp,dm+rho_comp,nscal, &
                                                fill_crse_input=.false.)
        
              enddo
     
           end if


           ! store the initial density field
           do n = 1,nlevs
              call multifab_copy_c(dens_orig(n),1,sold(n),rho_comp,1,ng_s)
           enddo


           ! write out the initial density field
           if (dump_output .and. .not. wrote_init_file) then

              if (dm == 2) then
                 outname = "dens_2d_orig"

              else if (dm == 3) then
                 outname = "dens_3d_orig"

              endif

              call fabio_ml_multifab_write_d(dens_orig,mla%mba%rr(:,1), &
                                             trim(outname),names=plot_names(1:1))

              print *, 'wrote file: ', trim(outname)
              print *, ' '
              
              wrote_init_file = .true.

           endif

           ! compute the initial timestep -- dt = dx / u
           dt = cflfac*dx(nlevs,1)/ONE


           ! advance the density using the constant velocity field
           t = ZERO
           do while (t < stop_time)

              print *, 't = ', t, 'dt = ', dt
     

              ! advance density according to rho_t + (rho U)_x = 0
              call density_advance(mla,1,sold,snew,sedge,sflux, &
                                   scal_force,umac,w0,w0mac,etarhoflux, &
                                   rho0_old,rho0_new,p0, &
                                   rho0_predicted_edge, &
                                   dx,dt,the_bc_tower%bc_tower_array)
              

              ! save the state for the next step
              do n = 1,nlevs
                 call multifab_copy_c(sold(n),1,snew(n),1,nscal,sold(n)%ng)
              enddo

              rho0_old = rho0_new

              ! update the time     
              t = t + dt

     
              ! adjust the timestep, if necessary
              if (t + dt > stop_time) then
                 dt = stop_time - t
              endif
              
           end do

           print *, 'finished evolution, t = ', t

           ! copy the final density field into a dummy multifab for
           ! output and analysis
           do n = 1,nlevs
              call multifab_copy_c(dens_final(n),1,snew(n),rho_comp,1,0)
           enddo


           ! compare the initial and final density

           ! compute dens_final - dens_orig
           do n = 1, nlevs
              call multifab_copy_c(error(n),1,dens_final(n),1,1,ng_s)
              call multifab_sub_sub(error(n), dens_orig(n))
           enddo
        
           do n = 1, nlevs
              abs_norm(p,n,index_t) = multifab_norm_l2(error(n))
           enddo
        
           ! now compute (dens_final - dens_orig)/dens_orig
           do n = 1, nlevs
              call multifab_div_div(error(n), dens_orig(n), 0)
           enddo

           do n = 1, nlevs
              rel_norm(p,n,index_t) = multifab_norm_l2(error(n))
           enddo


           if (dump_output) then

              do n = 1, nlevs
                 ! final density
                 call multifab_copy_c(temporary(n),1,dens_final(n),1,1,ng_s)

                 ! relative error 
                 call multifab_copy_c(temporary(n),2,error(n),1,1,ng_s)

                 ! absolute error (reconstructed)
                 call multifab_mult_mult(error(n), dens_orig(n), 0)
                 call multifab_copy_c(temporary(n),3,error(n),1,1,ng_s)
              enddo

              ! write out the initial density field the output name will
              ! include the dimensionality, ppm_type, and advection
              ! direction

              if (dm == 2) then
                 outname = "dens_2d_"
              else if (dm == 3) then
                 outname = "dens_3d_"
              endif
  
              if (bds_type == 1) then
                 outname = trim(outname) // "bds_"
              else
                 select case (p)
                 case (0)
                    outname = trim(outname) // "ppm0_"
                    
                 case (1)
                    outname = trim(outname) // "ppm1_"
                    
                 case (2)
                    outname = trim(outname) // "ppm2_"

                 case (3)
                    outname = trim(outname) // "bds_"
                 end select
              endif


              select case (itest_dir)
              case (-1)   
                 outname = trim(outname) // "xm_final"
                 
              case (1)
                 outname = trim(outname) // "xp_final"
                 
              case (-2)
                 outname = trim(outname) // "ym_final"
                 
              case (2)
                 outname = trim(outname) // "yp_final"
                 
              case (-3)
                 outname = trim(outname) // "zm_final"
                 
              case (3)
                 outname = trim(outname) // "zp_final"
              end select

              call fabio_ml_multifab_write_d(temporary,mla%mba%rr(:,1), &
                                             trim(outname),names=plot_names)

              print *, 'wrote file: ', trim(outname)
           
           endif


        enddo    ! idir loop
     enddo    ! idim loop
  enddo   ! ppm_type

  print *, ' '
  print *, ' '


  ! report
  do p = 0, 3   ! ppm_type

     if (p <= 2) then
        print *, 'ppm_type = ', p
     else
        print *, 'bds_type = ', 1
     endif

     do i = 1, 2*dm
        if (i <= dm) then
           itest_dir = i
        else
           itest_dir = -1*(i - dm)
        endif

        print *, 'advection direction: ', itest_dir
        do n = 1,nlevs
           print *, '  level ', n
           print *, '    | rho_final - rho_init |_2              = ', abs_norm(p,n,i)
           print *, '    | (rho_final - rho_init) / rho_init |_2 = ', rel_norm(p,n,i)
        enddo
        print *, ' '
     enddo
  enddo  ! ppm_type

  max_dabs_error = -1.d30
  max_abs_error  = -1.d30
  
  max_drel_error = -1.d30
  max_rel_error  = -1.d30
  
  do p = 0, 3
     do n = 1, nlevs
        max_dabs_error = max(abs(maxval(abs_norm(p,n,:)) - minval(abs_norm(p,n,:))), &
                             max_dabs_error)
        max_abs_error = max(maxval(abs_norm(p,n,:)), max_abs_error)

        max_drel_error = max(abs(maxval(rel_norm(p,n,:)) - minval(rel_norm(p,n,:))), &
                             max_drel_error)
        max_rel_error = max(maxval(rel_norm(p,n,:)), max_rel_error)
     enddo
  enddo  ! ppm_type

  if (max_dabs_error < advect_test_tol * max_abs_error .and. &
      max_drel_error < advect_test_tol * max_rel_error) then
     print *, "SUCCESS: advection errors for all directions agree with tolerance"
  else
     print *, "ERROR: advection errors across directions are too large"
     print *, "absolute error: ", max_dabs_error/max_abs_error
     print *, "relative error: ", max_drel_error/max_rel_error
  endif

  print *, " "


  ! clean-up
  do n = 1,nlevs
     call destroy(sold(n))
     call destroy(snew(n))
     do comp=1,dm
        call destroy(umac(n,comp))
     end do
     call destroy(dens_orig(n))
     call destroy(dens_final(n))
     call destroy(error(n))
     call destroy(temporary(n))
  end do

  call destroy(mla)
  call destroy(mba)

  deallocate(sold,snew,umac,dens_orig,dens_final,error,temporary)

  call bc_tower_destroy(the_bc_tower)

  call runtime_close()

  call destroy_geometry()

end subroutine varden



