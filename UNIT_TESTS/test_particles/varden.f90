! a driver that sets up a circular MAC velocity field and initializes 
! particles at various radii.  We then advect the particles around
! for one period.  We can check the accurary of the particle advection
! by comparing the initial radius to the final radius.


subroutine varden()

  use BoxLib
  use f2kcli
  use ml_boxarray_module
  use layout_module
  use multifab_module
  use ml_restriction_module
  use define_bc_module
  use bl_mem_stat_module
  use bl_timer_module
  use box_util_module
  use bl_IO_module
  use fabio_module
  use variables, only: nscal, init_variables, rho_comp, spec_comp
  use geometry, only:  spherical, &
                       init_spherical, init_center, &
                       init_cutoff
  use probin_module, only: prob_lo, prob_hi, pmask, &
                           test_set, cflfac, &
                           stop_time, vel_amp, &
                           probin_init, probin_close
  use initialize_module, only: initialize_bc, initialize_dx
  use bl_constants_module
  use particle_module
  use test_particles_module, only: init_umac_2d

  implicit none

  integer :: i,n,comp
  integer :: ng_s, ng_um
  integer :: nlevs, dm

  integer :: idir, idim, itest_dir, index_t

  type(ml_layout) :: mla

  type(multifab), allocatable :: umac(:,:)
  type(multifab), allocatable :: s(:)

  real(kind=dp_t), pointer :: ump(:,:,:,:), vmp(:,:,:,:), wmp(:,:,:,:)
  real(kind=dp_t), pointer :: sp(:,:,:,:)


  real(kind=dp_t), pointer :: dx(:,:)

  integer, allocatable :: lo(:),hi(:)

  type(ml_boxarray) :: mba

  type(bc_tower) ::  the_bc_tower

  real(kind=dp_t) :: t, dt, maxR

  type(particle_container) :: particles
  real(kind=dp_t), allocatable :: point(:)

  integer :: nparticles
  real(kind=dp_t) :: dr_part, dtheta_part, r, theta


  ! general Maestro initializations
  call probin_init()
  call init_spherical()
  call init_center()

  call init_variables()


  ! setup the grid
  call read_a_hgproj_grid(mba, test_set)

  call ml_layout_build(mla,mba,pmask)

  ! check for proper nesting
  if (.not. ml_boxarray_properly_nested(mla%mba, 3, pmask)) then
     call bl_error('fixed_grids not properly nested')
  end if

  ! initialize nlevs
  nlevs = mla%nlevel

  if (nlevs .ne. max_levs) then
     call bl_error('varden.f90: nlevs .ne. max_levs not supported yet')
  end if

  ! initialize dm
  dm = mla%dim

  ! initialize boundary conditions
  call initialize_bc(the_bc_tower,nlevs,pmask)
  do n = 1,nlevs
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  ! allocate states
  allocate(umac(nlevs,dm))
  allocate(s(nlevs))

  if (spherical == 1) then
     call bl_error("ERROR: test_particles not defined for spherical = 1")
  endif


  ! maximum length from center
  if (dm == 2) then
     maxR = max(HALF*(prob_hi(1)-prob_lo(1)), &
                HALF*(prob_hi(2)-prob_lo(2)))
  else
     maxR = max(HALF*(prob_hi(1)-prob_lo(1)), &
                HALF*(prob_hi(2)-prob_lo(2)), &
                HALF*(prob_hi(3)-prob_lo(3)))
  endif



  ! build MAC velocity
  do n = 1,nlevs
     call multifab_build(s(n), mla%la(n), 1, 1)

     do comp=1,dm
        call multifab_build_edge(umac(n,comp), mla%la(n),1,1,comp)
     end do
  end do


  ! initialize_dx
  call initialize_dx(dx,mba,nlevs)


  ! other allocations
  allocate(lo(dm))
  allocate(hi(dm))



  ! initialize the particles
  call build(particles)

  allocate(point(dm))

  nparticles = 20
  dr_part = maxR/nparticles
  dtheta_part = TWO*M_PI/nparticles

  do i = 1, nparticles
     theta = dble(i-ONE)*dtheta_part
     r = dble(i-HALF)*dr_part

     point(1) = r*cos(theta) + HALF*(prob_lo(1) + prob_hi(1))
     point(2) = r*sin(theta) + HALF*(prob_lo(2) + prob_hi(2))

     print *, point

     call add(particles,point,mla,dx,prob_lo)
  enddo

  call timestamp(particles, 'timestamp', s, &
                (/1/), (/"magvel"/), ZERO)


  ! initialize the MAC velocity field  ! negative velocity  
  ng_um = nghost(umac(1,1))
  ng_s  = nghost(s(1))

  do n = 1, nlevs
     do i = 1, umac(n,1)%nboxes
        if ( multifab_remote(umac(n,1),i) ) cycle

        ump  => dataptr(umac(n,1),i)
        vmp  => dataptr(umac(n,2),i)

        sp => dataptr(s(n), i)

        lo = lwb(get_box(umac(n,1), i))
        hi = upb(get_box(umac(n,1), i))
        
        select case (dm)
        case (2)
           call init_umac_2d(ump(:,:,1,1), vmp(:,:,1,1), ng_um, &
                              sp(:,:,1,1), ng_s, &
                             lo, hi, dx(n,:))

        case (3)
           wmp  => dataptr(umac(n,3),i)
           call bl_error("ERROR: 3-d not implemented")
        end select
     end do
  enddo


  ! write out a plotfile that contains the magnitude of the velocity
  ! field
  print *, mla%mba%rr(:,1)
  call fabio_ml_multifab_write_d(s,mla%mba%rr(:,1), &
                                 "magvel_field")


  ! compute the initial timestep -- dt = dx / |U|.  The velocity field
  ! has the magnitude |U| = A r, where A is the amplitude and r is the
  ! radius.
  dt = cflfac*dx(nlevs,1)/(vel_amp*maxR)


  ! advance the particls using the MAC velocity field
  t = ZERO
  do while (t < stop_time)

     print *, 't = ', t, 'dt = ', dt
     

     ! advance density according to rho_t + (rho U)_x = 0
     call particle_container_move_advect(particles,mla,umac,dx,dt, &
                                         prob_lo,prob_hi)
        
     
     ! update the time     
     t = t + dt


     ! write out a particle timestamp
     call timestamp(particles, 'timestamp', s, &
                    (/1/), (/"magvel"/), t)

        
     ! adjust the timestep, if necessary
     if (t + dt > stop_time) then
        dt = stop_time - t
     endif

  end do

  print *, 'finished evolution, t = ', t


  ! clean-up
  do n = 1,nlevs
     do comp=1,dm
        call destroy(umac(n,comp))
     end do
  end do

  call destroy(mla)
  call destroy(mba)

  deallocate(umac)

  call bc_tower_destroy(the_bc_tower)

  call probin_close()

end subroutine varden



