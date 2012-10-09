subroutine varden()

  use variables
  use network
  use geometry
  use init_module
  use base_state_module
  use base_io_module
  use estdt_module
  use firstdt_module
  use advance_timestep_module
  use box_util_module
  use fabio_module
  use bc_module
  use checkpoint_module
  use make_plotfile_module
  use restart_module
  use probin_module
  use bl_constants_module
  use average_module
  use make_grav_module
  use sponge_module
  use make_div_coeff_module
  use mk_vel_force_module
  use multifab_fill_ghost_module
  use fill_3d_module
  use multifab_physbc_module
  use eos_module
  use divu_iter_module
  use initial_proj_module
  use ml_restriction_module

  implicit none

  integer    :: init_step
  integer    :: istep_divu_iter,istep_init_iter,istep
  integer    :: ng_s,i,n,n_chk_comps
  integer    :: last_plt_written,last_chk_written
  real(dp_t) :: lenx,leny,lenz,max_dist,smin,smax
  real(dp_t) :: time,dt,dtold,dt_lev

  integer     , allocatable :: domain_phys_bc(:,:)
  integer     , allocatable :: lo(:), hi(:)
  logical     , allocatable :: nodal(:)
  real(dp_t)  , allocatable :: dx(:,:)
  type(box)   , allocatable :: domain_boxes(:)

  type(ml_layout) :: mla

  type(multifab), allocatable :: uold(:), unew(:), sold(:), snew(:)
  type(multifab), allocatable :: gpres(:), pres(:)
  type(multifab), allocatable :: vel_force(:)
  type(multifab), allocatable :: normal(:)
  type(multifab), allocatable :: sponge(:)
  type(multifab), allocatable :: rho_omegadot2(:)
  type(multifab), allocatable :: dSdt(:)         ! (S^n-S^{n-1})/t^{n-1}
  type(multifab), allocatable :: Source_old(:)   ! Source at timelevel n
  type(multifab), allocatable :: Source_new(:)   ! Source at timelevel n+1
  type(multifab), allocatable :: hgrhs(:)

  type(multifab), pointer     :: chkdata(:) => Null()
  type(multifab), pointer     :: chk_p(:) => Null()
  type(multifab), pointer     :: chk_dsdt(:) => Null()
  type(multifab), pointer     :: chk_src_old(:) => Null()
  type(multifab), pointer     :: chk_rho_omegadot2(:) => Null()

  real(kind=dp_t), pointer :: nop(:,:,:,:)

  character(len=8)               :: sd_name
  character(len=5)               :: plot_index
  character(len=256)             :: plot_file_name
  character(len=11)              :: base_state_name
  character(len=8)               :: base_w0_name
  character(len=9)               :: base_etarho_name
  character(len=20), allocatable :: plot_names(:)

  real(dp_t), parameter :: SMALL = 1.d-13

  logical :: init_mode

  type(ml_boxarray) :: mba

  real(dp_t)              :: r1, r2
  real(dp_t), allocatable :: div_coeff_old(:,:)
  real(dp_t), allocatable :: div_coeff_new(:,:)
  real(dp_t), allocatable :: grav_cell(:,:)
  real(dp_t), allocatable :: gamma1bar(:,:,:)
  real(dp_t), allocatable :: s0_old(:,:,:)
  real(dp_t), allocatable :: s0_new(:,:,:)
  real(dp_t), allocatable :: p0_old(:,:)
  real(dp_t), allocatable :: p0_new(:,:)
  real(dp_t), allocatable :: etarho(:,:)
  real(dp_t), allocatable :: psi(:,:) 
  real(dp_t), allocatable :: w0(:,:)
  real(dp_t), allocatable :: tempbar(:,:,:)

  type(bc_tower) ::  the_bc_tower

  type(box), allocatable :: boundingbox(:)

  type(boxarray), allocatable :: validboxarr(:)
  type(boxarray), allocatable :: diffboxarray(:)

  ng_s = 3

  last_plt_written = -1
  last_chk_written = -1

  call probin_init()

  if(parallel_IOProcessor() .and. do_initial_projection) then
     print*,'Warning: do_initial_projection = F'
     print*,'Setting init_divu_iter = init_iter = max_step = 0'
     print*,'So the output will be the loaded in data'

     init_divu_iter = 0
     init_iter = 0
     max_step = 0
  end if

  call init_spherical(spherical_in)
  call init_dm(dm_in)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize the arrays and read the restart data if restart >= 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (restart >= 0) then
     call fill_restart_data(restart, mba, chkdata, chk_p, chk_dsdt, &
                            chk_src_old, chk_rho_omegadot2, time, dt)
  else 
     call read_a_hgproj_grid(mba, test_set)
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize the variable index pointers and the reaction network
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(lo(dm),hi(dm))
  call init_variables(dm, nspec)
  call init_plot_variables(dm)
  call network_init()
  call eos_init(small_temp=small_temp,small_dens=small_dens)

  n_chk_comps = 2*dm + nscal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set up plot_names for writing plot files.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(plot_names(n_plot_comps))
  call get_plot_names(dm,plot_names)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! allocate storage for the state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if ( parallel_IOProcessor() ) &
       print *, 'pmask = ', pmask

  nlevs = mba%nlevel

  if (nlevs .ne. max_levs) then
     call bl_error('varden.f90: nlevs .ne. max_levs not supported yet')
  end if

  call ml_layout_build(mla,mba,pmask)

  allocate(uold(nlevs),sold(nlevs),gpres(nlevs),pres(nlevs),hgrhs(nlevs))
  allocate(dSdt(nlevs),Source_old(nlevs),Source_new(nlevs),rho_omegadot2(nlevs))
  allocate(unew(nlevs),snew(nlevs))
  allocate(vel_force(nlevs),normal(nlevs),sponge(nlevs))

  allocate(boundingbox(nlevs))

  allocate(nodal(dm))
  nodal = .true.

  do n = 1,nlevs
     call multifab_build(      uold(n),    mla%la(n),    dm, ng_s)
     call multifab_build(      sold(n),    mla%la(n), nscal, ng_s)
     call multifab_build(     gpres(n),    mla%la(n),    dm,    1)
     call multifab_build(      pres(n),    mla%la(n),     1,    1, nodal)
     call multifab_build(     hgrhs(n),    mla%la(n),     1,    0, nodal)
     call multifab_build(      dSdt(n),    mla%la(n),     1,    0)
     call multifab_build(Source_old(n),    mla%la(n),     1,    0)
     call multifab_build(Source_new(n),    mla%la(n),     1,    0)
     call multifab_build(rho_omegadot2(n), mla%la(n), nspec,    0)

     call setval(       uold(n), 0.0_dp_t, all=.true.)
     call setval(       sold(n), 0.0_dp_t, all=.true.)
     call setval(      gpres(n), 0.0_dp_t, all=.true.)
     call setval(       pres(n), 0.0_dp_t, all=.true.)
     call setval(      hgrhs(n), 0.0_dp_t, all=.true.)
     call setval( Source_old(n), 0.0_dp_t, all=.true.)
     call setval( Source_new(n), 0.0_dp_t, all=.true.)
     call setval(       dSdt(n), 0.0_dp_t, all=.true.)
     call setval(rho_omegadot2(n),0.0_dp_t,all=.true.)
  end do

  ! create a "bounding box" for each level
  do n=1,nlevs
     boundingbox(n) = get_box(sold(n),1)
     do i=2, sold(n)%nboxes
        boundingbox(n) = box_bbox(boundingbox(n),get_box(sold(n),i))
     end do
  end do

  ! compute diffboxarray
  ! each box in diffboxarray corresponds to an "empty space" between valid regions at 
  ! each level, excluding the coarsest level.
  ! I am going to use this to compute all of the intermediate r_start_coords and r_end_coords
  allocate(validboxarr(nlevs))
  allocate(diffboxarray(nlevs))
  do n=1,nlevs
     call boxarray_build_copy(validboxarr(n),get_boxarray(sold(n)))
     call boxarray_boxarray_diff(diffboxarray(n),boundingbox(n),validboxarr(n))
     call boxarray_simplify(diffboxarray(n))
  end do

  if (restart >= 0) then
     do n=1,nlevs
        call multifab_copy_c( uold(n),1,chkdata(n),1                ,dm)
        call multifab_copy_c( sold(n),1,chkdata(n),rho_comp+dm      ,nscal)
        call multifab_copy_c(gpres(n),1,chkdata(n),rho_comp+dm+nscal,dm)
        call destroy(chkdata(n))
     end do
        
     do n=1,nlevs
        call multifab_copy_c(pres(n),1,chk_p(n),1,1)       
        call destroy(chk_p(n))
     end do

     do n=1,nlevs
        call multifab_copy_c(dSdt(n),1,chk_dsdt(n),1,1)
        call destroy(chk_dsdt(n))
     end do

     do n=1,nlevs
        call multifab_copy_c(Source_old(n),1,chk_src_old(n),1,1)
        call destroy(chk_src_old(n))
     end do

     do n=1,nlevs
        call multifab_copy_c(rho_omegadot2(n),1,chk_rho_omegadot2(n),1,nspec)
        call destroy(chk_rho_omegadot2(n))
     end do

     deallocate(chkdata, chk_p, chk_dsdt, chk_src_old, chk_rho_omegadot2)
  end if

  if (parallel_IOProcessor()) then
     print *, 'number of processors = ', parallel_nprocs()
     print *, 'number of dimensions = ', dm

     do n = 1, nlevs
        print *, 'level: ', n
        print *, '   number of boxes = ', uold(n)%nboxes
        print *, '   maximum zones   = ', (extent(mba%pd(n),i),i=1,dm)
     end do

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! define the grid spacing on all levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(dx(nlevs,dm))

  do i = 1, dm
     dx(1,i) = (prob_hi(i)-prob_lo(i)) / real(extent(mba%pd(1),i),kind=dp_t)
  end do
  do n = 2,nlevs
     dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
  end do

  ! check to make sure our grid is square -- the solvers assume this
  if (dm == 2) then
     if (abs(dx(1,1) - dx(1,2)) > SMALL) then
        call bl_error('zones must be square')
     end if
  else if (dm == 3) then
     if (abs(dx(1,1) - dx(1,2)) > SMALL .OR. &
          abs(dx(1,1) - dx(1,3)) > SMALL) then
        call bl_error('zones must be square')
     end if
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! allocate storage for the base state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (spherical .eq. 1) then
     if (dr_base .gt. 0) then
        lenx = HALF * (prob_hi(1) - prob_lo(1))
        leny = HALF * (prob_hi(2) - prob_lo(2))
        lenz = HALF * (prob_hi(3) - prob_lo(3))
        max_dist = sqrt(lenx**2 + leny**2 + lenz**2)
        nr_fine = int(max_dist / dr_base) + 1
        if ( parallel_IOProcessor() ) then
           print *,'DISTANCE FROM CENTER TO CORNER IS ',max_dist
           print *,'DR_BASE IS ',dr_base
           print *,'SETTING nr_fine TO ',nr_fine
        end if
     else
        call bl_error('NEED TO DEFINE DR_BASE')
     end if
  else
     ! NOTE: WE ASSUME DR_BASE IS THE RESOLUTION OF THE FINEST LEVEL 
     ! IN PLANE-PARALLEL!
     nr_fine = extent(mba%pd(nlevs),dm)
     dr_base = (prob_hi(dm)-prob_lo(dm)) / dble(nr_fine)
  end if

  ! allocate the base state quantities -- these all use 0-based indexing
  allocate(div_coeff_old(nlevs,0:nr_fine-1))
  allocate(div_coeff_new(nlevs,0:nr_fine-1))
  allocate(gamma1bar    (nlevs,0:nr_fine-1,1))
  allocate(s0_old       (nlevs,0:nr_fine-1,nscal))
  allocate(s0_new       (nlevs,0:nr_fine-1,nscal))
  allocate(p0_old       (nlevs,0:nr_fine-1))
  allocate(p0_new       (nlevs,0:nr_fine-1))
  allocate(w0           (nlevs,0:nr_fine))
  allocate(etarho       (nlevs,0:nr_fine))
  allocate(psi          (nlevs,0:nr_fine))
  allocate(tempbar      (nlevs,0:nr_fine-1,1))

  div_coeff_old(:,:) = ZERO
  div_coeff_new(:,:) = ZERO

  gamma1bar(:,:,:) = ZERO

  s0_old(:,:,:) = ZERO
  s0_new(:,:,:) = ZERO

  etarho(:,:) = ZERO
  psi(:,:) = ZERO

  p0_old(:,:) = ZERO
  p0_new(:,:) = ZERO

  w0(:,:) = ZERO
  tempbar(:,:,:) = ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize all remaining arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do n = 1,nlevs
     call multifab_build(  unew(n), mla%la(n),    dm, ng_s)
     call multifab_build(  snew(n), mla%la(n), nscal, ng_s)
     call multifab_build(sponge(n), mla%la(n),     1, 0)

     call setval(  unew(n), ZERO, all=.true.)
     call setval(  snew(n), ZERO, all=.true.)
  end do

  if(spherical .eq. 1) then
     do n=1,nlevs
        call multifab_build(normal(n), mla%la(n),    dm, 1)
     end do
  end if

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
  if (pmask(1)) then
     if (bcx_lo .ne. -1 .or. bcx_hi .ne. -1) &
          call bl_error('MUST HAVE BCX = -1 if PMASK = T')
  end if
  if (dm > 1) then
     domain_phys_bc(2,1) = bcy_lo
     domain_phys_bc(2,2) = bcy_hi
     if (pmask(2)) then
        if (bcy_lo .ne. -1 .or. bcy_hi .ne. -1) &
             call bl_error('MUST HAVE BCY = -1 if PMASK = T') 
     end if
  end if
  if (dm > 2) then
     domain_phys_bc(3,1) = bcz_lo
     domain_phys_bc(3,2) = bcz_hi
     if (pmask(3)) then
        if (bcz_lo .ne. -1 .or. bcz_hi .ne. -1) &
             call bl_error('MUST HAVE BCZ = -1 if PMASK = T')
     end if
  end if

  do i = 1, dm
     if ( pmask(i) ) domain_phys_bc(i,:) = BC_PER
  end do

  ! Build the arrays for each grid from the domain_bc arrays.
  call bc_tower_build(the_bc_tower,mla,domain_phys_bc,domain_boxes)

  ! now that we have the_bc_tower, we can fill the ghost cells outside of the
  ! domain for the read-in checkpoint data
  if(restart >= 0) then
     
     if (nlevs .eq. 1) then
        
       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
        call multifab_fill_boundary(sold(nlevs))
        call multifab_fill_boundary(uold(nlevs))

        ! fill non-periodic domain boundary ghost cells
        call multifab_physbc(sold(nlevs),rho_comp,dm+rho_comp,nscal, &
                             the_bc_tower%bc_tower_array(nlevs))
        call multifab_physbc(uold(nlevs),       1,          1,   dm, &
                             the_bc_tower%bc_tower_array(nlevs))

     else

        ! the loop over nlevs must count backwards to make sure the finer grids are 
        ! done first
        do n = nlevs,2,-1

           ! set level n-1 data to be the average of the level n data covering it
           call ml_cc_restriction(uold(n-1),uold(n),mla%mba%rr(n-1,:))
           call ml_cc_restriction(sold(n-1),sold(n),mla%mba%rr(n-1,:))

           ! fill level n ghost cells using interpolation from level n-1 data
           ! note that multifab_fill_boundary and multifab_physbc are called for
           ! both levels n-1 and n
           call multifab_fill_ghost_cells(uold(n),uold(n-1), &
                                          ng_s,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1), &
                                          the_bc_tower%bc_tower_array(n  ), &
                                          1,1,dm,fill_crse_input=.false.)
           call multifab_fill_ghost_cells(sold(n),sold(n-1), &
                                          ng_s,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1), &
                                          the_bc_tower%bc_tower_array(n  ), &
                                          1,dm+rho_comp,nscal,fill_crse_input=.false.)
        end do

     end if

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Now initialize the grid data, and do initial projection if restart < 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(grav_cell(nlevs,nr_fine))

  ! Initialize geometry (IMPT: dr is set in init_base_state)
  center(1:dm) = HALF * (prob_lo(1:dm) + prob_hi(1:dm))
  if (parallel_IOProcessor()) & 
       print *,'DR_BASE ',dr_base
  call init_geometry(center,dr_base,nlevs,mla,boundingbox,diffboxarray)

  ! Initialize base state at each level independently
  if (restart < 0) then
     do n = 1,nlevs
        call init_base_state(n,model_file,s0_old(n,:,:),p0_old(n,:),dx(n,:))
     end do
  else
     write(unit=sd_name,fmt='("chk",i5.5)') restart
     write(unit=base_state_name,fmt='("model_",i5.5)') restart
     write(unit=base_w0_name,fmt='("w0_",i5.5)') restart
     write(unit=base_etarho_name,fmt='("eta_",i5.5)') restart

     call read_base_state(nlevs, base_state_name, base_w0_name, base_etarho_name, sd_name, &
                          s0_old, p0_old, gamma1bar(:,:,1), w0, etarho, div_coeff_old, psi)
  end if

  ! Create the normal array once we have defined "center"
  if (spherical .eq. 1) then
     do n = 1,nlevs
        do i = 1, nfabs(normal(n))
           nop => dataptr(normal(n), i)
           lo =  lwb(get_box(normal(n), i))
           hi =  upb(get_box(normal(n), i))
           call make_3d_normal(nop(:,:,:,:),lo,hi,dx(n,:),1)
        end do
     end do
  end if

  if (do_sponge) then
     call init_sponge(nlevs,s0_old(nlevs,:,:),prob_hi,dx(nlevs,:),prob_lo(dm))
  else
     do n=1,nlevs
        call setval(sponge(n), ONE, all=.true.)
     end do
  end if

  do n=1,nlevs
     call make_grav_cell(n,grav_cell(n,:),s0_old(n,:,rho_comp))
  end do

  if (restart < 0) then

     do n=1,nlevs
        call make_div_coeff(n,div_coeff_old(n,:),s0_old(n,:,rho_comp),p0_old(n,:), &
                            gamma1bar(n,:,1),grav_cell(n,:))
     end do

     time = ZERO
     dt = 1.d20

     call initveldata(nlevs,uold,s0_old,p0_old,dx,the_bc_tower%bc_tower_array,mla)
     call initscalardata(nlevs,sold,s0_old,p0_old,dx,the_bc_tower%bc_tower_array,mla)

     !----------------------------------------------------------------------
     ! Do an initial projection with omegadot = 0 and rho_Hext = 0
     !----------------------------------------------------------------------

     if(do_initial_projection) then
        call initial_proj(nlevs,uold,sold,pres,gpres,Source_old,hgrhs, &
                          div_coeff_old,s0_old,p0_old,gamma1bar(:,:,1),dx,the_bc_tower,mla)
     end if
    
     do n=1,nlevs
        call multifab_build(vel_force(n), mla%la(n), dm, 1)
     end do

     call mk_vel_force(nlevs,vel_force,gpres,sold,normal,s0_old,grav_cell,dx, &
                       the_bc_tower%bc_tower_array,mla)

     do n=1,nlevs
        call firstdt(n,uold(n),sold(n),vel_force(n),Source_old(n), &
                     p0_old(n,:),gamma1bar(n,:,1),dx(n,:),cflfac,dt_lev)

        dt_lev = min(dt_lev,0.05d0)

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
           print*,"Call to firstdt for level",n,"gives dt_lev =",dt_lev
        end if

        dt_lev = dt_lev*init_shrink
        if (parallel_IOProcessor() .and. verbose .ge. 1) then
           print*, "Multiplying dt_lev by init_shrink; dt_lev =",dt_lev
        end if

        dt = min(dt,dt_lev)
     end do

     do n=1,nlevs
        call destroy(vel_force(n))
     end do

     if (parallel_IOProcessor()) then
        print*,"Minimum firstdt over all levels =",dt
     end if

     if(fixed_dt .ne. -1.0d0) then
        dt = fixed_dt
        if (parallel_IOProcessor()) then
           print*, "Setting fixed dt =",dt
        end if
     end if

     !------------------------------------------------------------------------
     ! do the initial iterations to define a velocity field consistent with S
     !------------------------------------------------------------------------
     if ( parallel_IOProcessor() ) then
        print *, 'DOING',init_divu_iter ,'INITIAL DIVU ITERATIONS'
        print *, ' '
     end if

     do n=1,nlevs
        call multifab_build(vel_force(n), mla%la(n), dm, 1)
        call setval( vel_force(n),ZERO, all=.true.)
     end do

     do istep_divu_iter=1,init_divu_iter

        call divu_iter(nlevs,istep_divu_iter,uold,sold,pres,gpres,vel_force,normal, &
                       Source_old,hgrhs,dSdt,div_coeff_old,s0_old,p0_old,gamma1bar,w0, &
                       grav_cell,dx,dt,time,the_bc_tower,mla)

     end do

     do n=1,nlevs
        call destroy(vel_force(n))
     end do

  else if (restart >= 0) then

     if ( plot_int > 0 ) then

        ! tempbar is only used as an initial guess for eos calls
        call average(mla,sold,tempbar,dx,temp_comp,1,1)

        write(unit=plot_index,fmt='(i5.5)') restart
        plot_file_name = trim(plot_base_name) // plot_index
        call make_plotfile(plot_file_name,mla,uold,sold,gpres,rho_omegadot2,&
                           Source_old,sponge,mba,plot_names,time,dx,&
                           the_bc_tower,s0_old,p0_old,tempbar,plot_spec,plot_trac)

        call write_job_info(plot_file_name)
     end if

  end if

  do n = 1,nlevs

     ! This is done to impose any Dirichlet bc's on unew or snew.
     call multifab_copy_c(unew(n),1,uold(n),1,dm   ,ng=ng_s)
     call multifab_copy_c(snew(n),1,sold(n),1,nscal,ng=ng_s)

  end do

  if (restart < 0) then 
     istep = 0
  else 
     istep = restart
  end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Begin the initial iterations to define an initial pressure field.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (stop_time >= 0.d0) then
     if (time+dt > stop_time) dt = min(dt, stop_time - time)
  end if

  dtold = dt

  if (restart < 0) then
     ! initialize Source_new to the Source_old the first time through
     do n = 1,nlevs
        call multifab_copy_c(Source_new(n),1,Source_old(n),1,1)
     end do

     if (do_sponge) then
        call make_sponge(nlevs,sponge,dx,dt,mla)
     end if

     if (init_iter > 0) then
        if (parallel_IOProcessor() .and. verbose .ge. 1) then
           print*,'DOING',init_iter,'INITIAL PRESSURE ITERATIONS'
           print*,''
        end if

        !----------------------------------------------------------------------
        ! main initial iteration loop
        !----------------------------------------------------------------------

        do istep_init_iter = 1,init_iter

           ! Advance a single timestep at all levels.
           init_mode = .true.
           call advance_timestep(init_mode,mla,uold,sold,unew,snew, &
                                 gpres,pres,normal,s0_old, &
                                 s0_new,p0_old,p0_new,tempbar,gamma1bar,w0, &
                                 rho_omegadot2,div_coeff_old, &
                                 div_coeff_new,grav_cell,dx,time,dt,dtold,the_bc_tower, &
                                 dSdt,Source_old,Source_new,etarho,psi,sponge,hgrhs,istep)

        end do ! end do istep_init_iter = 1,init_iter

     end if ! end if (init_iter > 0)

     !------------------------------------------------------------------------
     ! write a checkpoint file
     !------------------------------------------------------------------------

     if ( chk_int > 0 ) then

        allocate(chkdata(nlevs))
        do n = 1,nlevs
           call multifab_build(chkdata(n), mla%la(n), n_chk_comps, 0)
           call multifab_copy_c(chkdata(n),1                ,uold(n), 1,dm)
           call multifab_copy_c(chkdata(n),rho_comp+dm      ,sold(n), 1,nscal)
           call multifab_copy_c(chkdata(n),rho_comp+dm+nscal,gpres(n),1,dm)
        end do
        write(unit=sd_name,fmt='("chk",i5.5)') istep

        call checkpoint_write(sd_name, chkdata, pres, dSdt, Source_old, &
                              rho_omegadot2, mba%rr, time, dt)

        last_chk_written = istep

        write(unit=base_state_name,fmt='("model_",i5.5)') istep
        write(unit=base_w0_name,fmt='("w0_",i5.5)') istep
        write(unit=base_etarho_name,fmt='("eta_",i5.5)') istep
        call write_base_state(nlevs,base_state_name,base_w0_name,base_etarho_name,sd_name, &
                              s0_old,p0_old,gamma1bar(:,:,1),w0,etarho, &
                              div_coeff_old,psi,prob_lo(dm))
        do n = 1,nlevs
           call destroy(chkdata(n))
        end do
        deallocate(chkdata)
     end if

     if ( plot_int > 0 ) then
        write(unit=plot_index,fmt='(i5.5)') istep
        plot_file_name = trim(plot_base_name) // plot_index
        call make_plotfile(plot_file_name,mla,uold,sold,gpres,rho_omegadot2,Source_new, &
                           sponge,mba,plot_names,time,dx,&
                           the_bc_tower, &
                           s0_old,p0_old,tempbar,plot_spec,plot_trac)

        call scalar_diags (istep,sold(1),s0_old(1,:,:),p0_old(1,:),dx(1,:))

        call write_job_info(plot_file_name)
        last_plt_written = istep
     end if

  end if ! end if (restart < 0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! main evolution loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (restart < 0) then
     init_step = 1
  else
     init_step = restart+1
  end if

  if ( parallel_IOProcessor()) then
     print*,""
     print*,"BEGIN MAIN EVOLUTION LOOP WITH dt =",dt
     print*,""
  end if

  if ( (max_step >= init_step) .and. &
       (time < stop_time .or. stop_time < 0.d0) ) then

     do istep = init_step,max_step

        if ( verbose .ge. 1 ) then
           if ( parallel_IOProcessor() ) then
              print *, 'MEMORY STATS AT START OF TIMESTEP ', istep
              print*, ' '
           end if
           call print(multifab_mem_stats(),    "    multifab")
           call print(fab_mem_stats(),         "         fab")
           call print(boxarray_mem_stats(),    "    boxarray")
           call print(layout_mem_stats(),      "      layout")
           call print(boxassoc_mem_stats(),    "    boxassoc")
           call print(fgassoc_mem_stats(),     "     fgassoc")
           call print(syncassoc_mem_stats(),   "   syncassoc")
           call print(copyassoc_mem_stats(),   "   copyassoc")
           call print(fluxassoc_mem_stats(),   "   fluxassoc")
           if ( parallel_IOProcessor() ) print*, ''
        end if

        do n = 1,nlevs

           if ( verbose .ge. 1 ) then
              smax = norm_inf(uold(n),1,1)
              if ( parallel_IOProcessor()) &
                   print *,'MAX OF UOLD ', smax,' AT LEVEL ',n

              smax = norm_inf(uold(n),2,1)
              if ( parallel_IOProcessor()) &
                   print *,'MAX OF VOLD ', smax,' AT LEVEL ',n

              if (dm > 2) then
                 smax = norm_inf(uold(n),3,1)
                 if ( parallel_IOProcessor()) &
                      print *,'MAX OF WOLD ', smax,' AT LEVEL ',n
              end if
           end if

        end do

        !---------------------------------------------------------------------
        ! get the new timestep
        !---------------------------------------------------------------------
        dtold = dt

        do n=1,nlevs
           call make_grav_cell(n,grav_cell(n,:),s0_old(n,:,rho_comp))
        end do

        if (istep > 1) then

           dt = 1.e20

           do n=1,nlevs
              call multifab_build(vel_force(n), mla%la(n), dm, 1)
              call setval( vel_force(n),ZERO, all=.true.)
           end do

           call mk_vel_force(nlevs,vel_force,gpres,sold,normal,s0_old,grav_cell,dx, &
                             the_bc_tower%bc_tower_array,mla)

           do n = 1,nlevs
              call estdt(n,uold(n),sold(n),vel_force(n),Source_old(n),dSdt(n),normal(n), &
                         w0(n,:),p0_old(n,:),gamma1bar(n,:,1),dx(n,:),cflfac,dt_lev)

              dt = min(dt,dt_lev)
           end do

           do n=1,nlevs
              call destroy(vel_force(n))
           end do

           if (parallel_IOProcessor() .and. verbose .ge. 1) then
              print*,''
              print*,"Call to estdt at beginning of step",istep
              print*,"gives dt =",dt
           end if

           if(dt .gt. max_dt_growth*dtold) then
              dt = max_dt_growth*dtold
              if (parallel_IOProcessor() .and. verbose .ge. 1) then
                 print*,'dt_growth factor limits the new dt =',dt
              end if
           end if

           if(fixed_dt .ne. -1.0d0) then
              dt = fixed_dt
              if (parallel_IOProcessor()) then
                 print*, "Setting fixed dt =",dt
              end if
           end if

           if (parallel_IOProcessor() .and. verbose .ge. 1) then
              print*,''
           end if

        end if

        if (stop_time >= 0.d0) then
           if (time+dt > stop_time) then
              dt = stop_time - time 
              if (parallel_IOProcessor()) &
                   print*, "Stop time limits dt =",dt
           end if
        end if

        !---------------------------------------------------------------------
        ! Advance a single timestep at all levels.
        !---------------------------------------------------------------------
        init_mode = .false.
        if (do_sponge) then
           call init_sponge(nlevs,s0_old(nlevs,:,:),prob_hi,dx(nlevs,:),prob_lo(dm))
           call make_sponge(nlevs,sponge,dx,dt,mla)
        end if
        r1 = parallel_wtime()
        call advance_timestep(init_mode,mla,uold,sold,unew,snew, &
                              gpres,pres,normal,s0_old, &
                              s0_new,p0_old,p0_new,tempbar,gamma1bar,w0, &
                              rho_omegadot2,div_coeff_old,div_coeff_new, &
                              grav_cell,dx,time,dt,dtold,the_bc_tower, &
                              dSdt,Source_old,Source_new,etarho,psi,&
                              sponge,hgrhs,istep)
        r2 = parallel_wtime() - r1
        call parallel_reduce(r1, r2, MPI_MAX, proc = parallel_IOProcessorNode())
        if (parallel_IOProcessor()) print*, 'Time to advance timestep: ', r1, ' seconds'

        call print_and_reset_fab_byte_spread()

        time = time + dt

        if ( verbose .ge. 1 ) then

           if ( parallel_IOProcessor() ) then
              print *, 'MEMORY STATS AT END OF TIMESTEP ', istep
              print*, ' '
           end if
           call print(multifab_mem_stats(),    "    multifab")
           call print(fab_mem_stats(),         "         fab")
           call print(boxarray_mem_stats(),    "    boxarray")
           call print(layout_mem_stats(),      "      layout")
           call print(boxassoc_mem_stats(),    "    boxassoc")
           call print(fgassoc_mem_stats(),     "     fgassoc")
           call print(syncassoc_mem_stats(),   "   syncassoc")
           call print(copyassoc_mem_stats(),   "   copyassoc")
           call print(fluxassoc_mem_stats(),   "   fluxassoc")
           if ( parallel_IOProcessor() ) print*, ''

           do n = 1,nlevs
              if (parallel_IOProcessor()) write(6,1100) n

              smin = multifab_min_c(unew(n),1) 
              smax = multifab_max_c(unew(n),1)
              if (parallel_IOProcessor()) write(6,1101) smin,smax

              smin = multifab_min_c(unew(n),2) 
              smax = multifab_max_c(unew(n),2)
              if (parallel_IOProcessor()) write(6,1102) smin,smax

              if (dm .eq. 3) then
                 smin = multifab_min_c(unew(n),3) 
                 smax = multifab_max_c(unew(n),3)
                 if (parallel_IOProcessor()) write(6,1103) smin,smax
              end if

              smin = multifab_min_c(gpres(n),1) 
              smax = multifab_max_c(gpres(n),1)
              if (parallel_IOProcessor()) write(6,1104) smin,smax

              smin = multifab_min_c(gpres(n),2) 
              smax = multifab_max_c(gpres(n),2)
              if (parallel_IOProcessor()) write(6,1105) smin,smax

              if (dm .eq. 3) then
                 smin = multifab_min_c(gpres(n),3) 
                 smax = multifab_max_c(gpres(n),3)
                 if (parallel_IOProcessor()) write(6,1106) smin,smax
              end if

              if (parallel_IOProcessor()) write(6,1107)
           end do
        end if

1100    format('At level ',i3)
1101    format('... min/max : x-velocity       ',e17.10,2x,e17.10)
1102    format('... min/max : y-velocity       ',e17.10,2x,e17.10)
1103    format('... min/max : z-velocity       ',e17.10,2x,e17.10)
1104    format('... min/max : gpresx              ',e17.10,2x,e17.10)
1105    format('... min/max : gpresy              ',e17.10,2x,e17.10)
1106    format('... min/max : gpresz              ',e17.10,2x,e17.10)
1107    format(' ')

        if (parallel_IOProcessor()) write(6,1000) istep,time,dt

        !---------------------------------------------------------------------
        ! save the old velocity, scalars and source terms for the next step
        !---------------------------------------------------------------------

        do n = 1,nlevs
           call multifab_copy_c(uold(n),      1,unew(n),      1,dm,   ng_s)
           call multifab_copy_c(sold(n),      1,snew(n),      1,nscal,ng_s)
           call multifab_copy_c(Source_old(n),1,Source_new(n),1,1)
        end do

        ! Set div_coeff_old equal to div_coeff_new from the last time step
        div_coeff_old = div_coeff_new

        ! Copy the base state
        s0_old = s0_new
        p0_old = p0_new

        !---------------------------------------------------------------------
        ! output
        !---------------------------------------------------------------------

        if (chk_int > 0) then
           if (mod(istep,chk_int) .eq. 0) then
              allocate(chkdata(nlevs))
              do n = 1,nlevs
                 call multifab_build(chkdata(n), mla%la(n), n_chk_comps, 0)
                 call multifab_copy_c(chkdata(n),1,unew(n),1,dm)
                 call multifab_copy_c(chkdata(n),rho_comp+dm,snew(n),1,nscal)
                 call multifab_copy_c(chkdata(n),rho_comp+dm+nscal,gpres(n),1,dm)
              end do
              write(unit=sd_name,fmt='("chk",i5.5)') istep
              call checkpoint_write(sd_name, chkdata, pres, dSdt,Source_old, &
                                    rho_omegadot2, mba%rr, time, dt)
              last_chk_written = istep

              write(unit=base_state_name,fmt='("model_",i5.5)') istep
              write(unit=base_w0_name,fmt='("w0_",i5.5)') istep
              write(unit=base_etarho_name,fmt='("eta_",i5.5)') istep
              call write_base_state(nlevs, base_state_name, base_w0_name, &
                                    base_etarho_name, sd_name, &
                                    s0_new, p0_new, gamma1bar(:,:,1), w0, etarho, &
                                    div_coeff_old,psi, &
                                    prob_lo(dm))
              do n = 1,nlevs
                 call destroy(chkdata(n))
              end do
              deallocate(chkdata)

           end if
        end if

        if (plot_int > 0) then
           if (mod(istep,plot_int) .eq. 0) then
              write(unit=plot_index,fmt='(i5.5)') istep
              plot_file_name = trim(plot_base_name) // plot_index
              call make_plotfile(plot_file_name,mla,unew,snew,gpres,&
                                 rho_omegadot2,Source_new,sponge, &
                                 mba,plot_names,time,dx,the_bc_tower, &
                                 s0_new,p0_new,tempbar,plot_spec,plot_trac)

              call write_job_info(plot_file_name)
              last_plt_written = istep

              call scalar_diags (istep,snew(1),s0_new(1,:,:),p0_new(1,:),dx(1,:))
           end if
        end if

        if (stop_time >= 0.d0) then
           if (time >= stop_time) goto 999
        end if

     end do

999  continue
     if (istep > max_step) istep = max_step

1000 format('STEP = ',i5,1x,' TIME = ',f16.10,1x,'DT = ',f14.9)

     if ( chk_int > 0 .and. last_chk_written .ne. istep ) then
        !       This writes a checkpoint file.
        allocate(chkdata(nlevs))
        do n = 1,nlevs
           call multifab_build(chkdata(n), mla%la(n), n_chk_comps, 0)
           call multifab_copy_c(chkdata(n),1,unew(n),1,dm)
           call multifab_copy_c(chkdata(n),rho_comp+dm,snew(n),1,nscal)
           call multifab_copy_c(chkdata(n),rho_comp+dm+nscal,gpres(n),1,dm)
        end do
        write(unit=sd_name,fmt='("chk",i5.5)') istep
        call checkpoint_write(sd_name, chkdata, pres, dSdt, Source_old, &
                              rho_omegadot2, mba%rr, time, dt)
        do n = 1,nlevs
           call destroy(chkdata(n))
        end do
        deallocate(chkdata)

        write(unit=base_state_name,fmt='("model_",i5.5)') istep
        write(unit=base_w0_name,fmt='("w0_",i5.5)') istep
        write(unit=base_etarho_name,fmt='("eta_",i5.5)') istep
        call write_base_state(nlevs, base_state_name, base_w0_name, &
                              base_etarho_name, sd_name, &
                              s0_new, p0_new, gamma1bar(:,:,1), w0, etarho, &
                              div_coeff_old, psi,prob_lo(dm))
     end if

     if ( plot_int > 0 .and. last_plt_written .ne. istep ) then
        write(unit=plot_index,fmt='(i5.5)') istep
        plot_file_name = trim(plot_base_name) // plot_index
        call make_plotfile(plot_file_name,mla,unew,snew,gpres, &
                           rho_omegadot2,Source_new,sponge, &
                           mba,plot_names,time,dx,the_bc_tower, &
                           s0_new,p0_new,tempbar,plot_spec,plot_trac)
        
        call write_job_info(plot_file_name)

        call scalar_diags (istep,snew(1),s0_new(1,:,:),p0_new(1,:),dx(1,:))
     end if
  end if

  do n=1,nlevs
     call destroy(uold(n))
     call destroy(unew(n))
     call destroy(sold(n))
     call destroy(snew(n))
     call destroy(pres(n))
     call destroy(gpres(n))
     call destroy(dSdt(n))
     call destroy(Source_old(n))
     call destroy(Source_new(n))
     call destroy(hgrhs(n))
     call destroy(rho_omegadot2(n))
     call destroy(sponge(n))
  end do

  do n=1,nlevs
     call destroy(validboxarr(n))
     call destroy(diffboxarray(n))
  end do

  if(spherical .eq. 1) then
     do n=1,nlevs
        call destroy(normal(n))
     end do
  end if

  call destroy(mla)
  call destroy(mba)

  deallocate(uold,unew,sold,snew)
  deallocate(div_coeff_old,div_coeff_new,grav_cell)
  deallocate(gamma1bar,s0_old,s0_new,p0_old,p0_new,w0,etarho)
  deallocate(dSdt,Source_old,Source_new)
  deallocate(rho_omegadot2)
  deallocate(vel_force,sponge,normal)
  deallocate(gpres,pres,hgrhs)
  deallocate(nodal)
  deallocate(dx)
  deallocate(domain_phys_bc,domain_boxes)

  deallocate(lo,hi,boundingbox)
  deallocate(plot_names)

  call bc_tower_destroy(the_bc_tower)

  call destroy_geometry()

end subroutine varden
