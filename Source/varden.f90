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
  use define_bc_module
  use fill_3d_module
  use eos_module
  use divu_iter_module
  use initial_proj_module
  use make_gamma_module
  use restrict_base_module
  use enforce_HSE_module
  use rhoh_vs_t_module
  use initialize_module
  use make_new_grids_module
  use regrid_module

  implicit none

  integer    :: init_step,istep
  integer    :: istep_divu_iter,istep_init_iter
  integer    :: i,n,r,comp
  integer    :: last_plt_written,last_chk_written
  real(dp_t) :: smin,smax
  real(dp_t) :: time,dt,dtold

  type(ml_layout) :: mla
  type(bc_tower)  :: the_bc_tower

  real(dp_t)  , pointer     :: dx(:,:)

  type(multifab), allocatable :: unew(:)
  type(multifab), allocatable :: snew(:)
  type(multifab), allocatable :: vel_force(:)
  type(multifab), allocatable :: normal(:)
  type(multifab), allocatable :: sponge(:)
  type(multifab), allocatable :: Source_new(:)
  type(multifab), allocatable :: hgrhs(:)
  type(multifab), allocatable :: gamma1(:)
  type(multifab), allocatable :: w0mac(:,:)

  ! these are pointers because they need to be allocated and built within another function
  type(multifab), pointer :: uold(:)
  type(multifab), pointer :: sold(:)
  type(multifab), pointer :: pres(:)
  type(multifab), pointer :: gpres(:)
  type(multifab), pointer :: dSdt(:)
  type(multifab), pointer :: Source_old(:)
  type(multifab), pointer :: rho_omegadot2(:)
  type(multifab), pointer :: rho_Hnuc2(:)

  type(multifab), pointer :: chkdata(:)

  character(len=5)               :: plot_index, check_index
  character(len=256)             :: plot_file_name, check_file_name
  character(len=11)              :: base_state_name
  character(len=8)               :: base_w0_name
  character(len=9)               :: base_etarho_name
  character(len=20), allocatable :: plot_names(:)

  real(dp_t), parameter :: SMALL = 1.d-13
  real(dp_t)            :: runtime1, runtime2

  logical :: init_mode

  real(dp_t), pointer :: div_coeff_old(:,:)
  real(dp_t), pointer :: div_coeff_new(:,:)
  real(dp_t), pointer :: gamma1bar(:,:)
  real(dp_t), pointer :: gamma1bar_hold(:,:)
  real(dp_t), pointer :: s0_init(:,:,:)
  real(dp_t), pointer :: rho0_old(:,:)
  real(dp_t), pointer :: rhoh0_old(:,:)
  real(dp_t), pointer :: rho0_new(:,:)
  real(dp_t), pointer :: rhoh0_new(:,:)
  real(dp_t), pointer :: p0_init(:,:)
  real(dp_t), pointer :: p0_old(:,:)
  real(dp_t), pointer :: p0_new(:,:)
  real(dp_t), pointer :: w0(:,:)
  real(dp_t), pointer :: etarho_ec(:,:)
  real(dp_t), pointer :: etarho_cc(:,:)
  real(dp_t), pointer :: psi(:,:)
  real(dp_t), pointer :: tempbar(:,:)
  real(dp_t), pointer :: grav_cell(:,:)

  real(dp_t), allocatable :: psi_temp(:,:)
  real(dp_t), allocatable :: etarho_cc_temp(:,:)
  real(dp_t), allocatable :: etarho_ec_temp(:,:)
  real(dp_t), allocatable :: w0_temp(:,:)

  last_plt_written = -1
  last_chk_written = -1

  call probin_init()
  call init_dm()
  call init_spherical()
  call init_center()
  call init_rotation()

  call init_variables()
  call init_plot_variables()

  call network_init()
  call eos_init(use_eos_coulomb=use_eos_coulomb,small_temp=small_temp,small_dens=small_dens)

  allocate(plot_names(n_plot_comps))
  call get_plot_names(plot_names)

  if (restart >= 0) then

     call initialize_from_restart(mla,restart,time,dt,pmask,dx,uold,sold,gpres,pres, &
                                  dSdt,Source_old,rho_omegadot2,rho_Hnuc2,the_bc_tower, &
                                  div_coeff_old,div_coeff_new,gamma1bar,gamma1bar_hold, &
                                  s0_init,rho0_old,rhoh0_old,rho0_new,rhoh0_new,p0_init, &
                                  p0_old,p0_new,w0,etarho_ec,etarho_cc,psi,tempbar,grav_cell)
     

  else if (test_set /= '') then

     call initialize_with_fixed_grids(mla,time,dt,pmask,dx,uold,sold,gpres,pres,dSdt, &
                                      Source_old,rho_omegadot2,rho_Hnuc2,the_bc_tower, &
                                      div_coeff_old,div_coeff_new,gamma1bar, &
                                      gamma1bar_hold,s0_init,rho0_old,rhoh0_old, &
                                      rho0_new,rhoh0_new,p0_init,p0_old,p0_new,w0, &
                                      etarho_ec,etarho_cc,psi,tempbar,grav_cell)

  else

     call initialize_with_adaptive_grids(mla,time,dt,pmask,dx,uold,sold,gpres,pres,dSdt, &
                                         Source_old,rho_omegadot2,rho_Hnuc2,the_bc_tower, &
                                         div_coeff_old,div_coeff_new,gamma1bar, &
                                         gamma1bar_hold,s0_init,rho0_old,rhoh0_old, &
                                         rho0_new,rhoh0_new,p0_init,p0_old,p0_new,w0, &
                                         etarho_ec,etarho_cc,psi,tempbar,grav_cell)

  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! error checking
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! check to make sure dimensionality is consistent in the inputs file
  if (dm .ne. get_dim(mla%mba)) then 
     call bl_error('dm_in not properly set in inputs file')
  end if

  ! check to make sure spherical is only used for 3d
  if (spherical .eq. 1 .and. dm .ne. 3) then
     call bl_error("spherical = 1 and dm != 3")
  end if
  
  ! check to make sure prob_lo for spherical is properly defined
  if(spherical .eq. 1) then
     do i=1,dm
        if(prob_lo(i) .ne. ZERO) then
           call bl_error('Error: prob_lo for spherical is not zero')
        end if
     end do
  end if

  ! check to make sure our grid is square -- the solvers assume this
  if (dm == 2) then
     if (abs(dx(1,1) - dx(1,2)) > SMALL) then
        call bl_error('zones must be square')
     end if
  else if (dm == 3) then
     if (abs(dx(1,1) - dx(1,2)) > SMALL .OR. abs(dx(1,1) - dx(1,3)) > SMALL) then
        call bl_error('zones must be square')
     end if
  end if

  if (parallel_IOProcessor()) then
     print *, 'number of processors = ', parallel_nprocs()
     print *, 'number of dimensions = ', dm

     do n = 1, nlevs
        print *, 'level: ', n
        print *, '   number of boxes = ', uold(n)%nboxes
        print *, '   maximum zones   = ', (extent(mla%mba%pd(n),i),i=1,dm)
     end do
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize all remaining arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(       psi_temp(max_levs,0:nr_fine-1))
  allocate( etarho_cc_temp(max_levs,0:nr_fine-1))
  allocate( etarho_ec_temp(max_levs,0:nr_fine))
  allocate(        w0_temp(max_levs,0:nr_fine))

  allocate(unew(nlevs),snew(nlevs),sponge(nlevs),hgrhs(nlevs),Source_new(nlevs))
  allocate(normal(nlevs))
  allocate(vel_force(nlevs))

  do n = 1,nlevs
     call multifab_build(      unew(n), mla%la(n),    dm, 3)
     call multifab_build(      snew(n), mla%la(n), nscal, 3)
     call multifab_build(    sponge(n), mla%la(n),     1, 0)
     call multifab_build(     hgrhs(n), mla%la(n),     1, 0, nodal)
     call multifab_build(Source_new(n), mla%la(n),     1, 1)
     if (dm .eq. 3) then
        call multifab_build(normal(n), mla%la(n),    dm, 1)
     end if

     call setval(      unew(n), ZERO, all=.true.)
     call setval(      snew(n), ZERO, all=.true.)
     call setval(    sponge(n), ONE,  all=.true.)
     call setval(     hgrhs(n), ZERO, all=.true.)
     call setval(Source_new(n), ZERO, all=.true.)

  end do

  ! Create normal now that we have defined center and dx
  call make_normal(normal,dx)

  call compute_cutoff_coords(rho0_old)

  if (do_sponge) then
     call init_sponge(rho0_old(1,:),dx(1,:),prob_lo(dm))
  end if

  call make_grav_cell(grav_cell,rho0_old)

  ! enforce HSE and then ensure state is thermodynamically consistens
  if (restart < 0) then

     ! this section is not needed for spherical since there is only
     ! one level of refinement of the base state, no matter what nlevs is
     if ( (nlevs .gt. 1 .or. perturb_model) .and. spherical .eq. 0) then

        ! enforce HSE
        call enforce_HSE(rho0_old,p0_old,grav_cell)

        ! compute full state T,h = T,h(rho,p0_old,X)
        call makeTHfromRhoP(sold,p0_old,the_bc_tower%bc_tower_array,mla,dx)

        ! force rhoh0 to be the average of rhoh
        call average(mla,sold,rhoh0_old,dx,rhoh_comp)

     end if

  end if

  if (restart < 0) then

     !----------------------------------------------------------------------
     ! Do an initial projection with omegadot = 0 and rho_Hext = 0
     !----------------------------------------------------------------------

     allocate(gamma1(nlevs))

     do n=1,nlevs
        call multifab_build(gamma1(n), mla%la(n), 1, 0)
     end do

     ! tempbar is only used as an initial guess for eos calls
     call average(mla,sold,tempbar,dx,temp_comp)
     tempbar = max(tempbar,small_temp)

     call make_gamma(mla,gamma1,sold,p0_old,tempbar,dx)
     call average(mla,gamma1,gamma1bar,dx,1)
     
     do n=1,nlevs
        call destroy(gamma1(n))
     end do

     call make_div_coeff(div_coeff_old,rho0_old,p0_old,gamma1bar,grav_cell)

     if(do_initial_projection > 0) then
        call initial_proj(uold,sold,pres,gpres,Source_old,hgrhs, &
                          div_coeff_old,p0_old,gamma1bar,dx,the_bc_tower,mla)
     end if

     !----------------------------------------------------------------------
     ! Compute the initial time step
     !----------------------------------------------------------------------
    
     do n=1,nlevs
        call multifab_build(vel_force(n), mla%la(n), dm, 1)
     end do

     call mk_vel_force(vel_force,uold,gpres,sold,normal,rho0_old,grav_cell,dx, &
                       the_bc_tower%bc_tower_array,mla)

     call firstdt(uold,sold,vel_force,Source_old,normal,p0_old,gamma1bar,dx,cflfac,dt)

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

        call divu_iter(istep_divu_iter,uold,sold,pres,gpres,vel_force,normal, &
                       Source_old,hgrhs,dSdt,div_coeff_old,rho0_old,p0_old, &
                       gamma1bar,tempbar,w0,grav_cell,dx,dt,time,the_bc_tower,mla)

     end do

     do n=1,nlevs
        call destroy(vel_force(n))
     end do

  else if (restart >= 0) then

     if ( plot_int > 0 ) then

        !------------------------------------------------------------------------
        ! write out a plotfile
        !------------------------------------------------------------------------

        ! tempbar is only used as an initial guess for eos calls
        ! Note: if we keep enthalpy_pred_type .eq. predict_Tprime_then_h, then
        ! we need to instead store tempbar (or t0 if we renname it) in the checkpoint file
        call average(mla,sold,tempbar,dx,temp_comp)
        tempbar = max(tempbar,small_temp)

        ! we generate this only for the initial plotfile
        if (do_sponge) then
           call make_sponge(sponge,dx,dt,mla)
        end if

        write(unit=plot_index,fmt='(i5.5)') restart
        plot_file_name = trim(plot_base_name) // plot_index
        call make_plotfile(plot_file_name,mla,uold,sold,gpres,rho_omegadot2,rho_Hnuc2, &
                           Source_old,sponge,mla%mba,plot_names,time,dx, &
                           the_bc_tower,w0,rho0_old,rhoh0_old,p0_old,tempbar,gamma1bar, &
                           div_coeff_old,normal)

        call write_job_info(plot_file_name)

     end if

  end if

  do n = 1,nlevs

     ! This is done to impose any Dirichlet bc's on unew or snew.
     call multifab_copy_c(unew(n),1,uold(n),1,dm   ,3)
     call multifab_copy_c(snew(n),1,sold(n),1,nscal,3)

  end do

  if (restart < 0) then 
     istep = 0
  else 
     istep = restart
  end if

  if (stop_time >= 0.d0) then
     if (time+dt > stop_time) dt = min(dt, stop_time - time)
  end if

  dtold = dt

  if (restart < 0) then

     !------------------------------------------------------------------------
     ! Begin the initial iterations to define an initial pressure field.
     !------------------------------------------------------------------------

     ! initialize Source_new to the Source_old the first time through
     do n = 1,nlevs
        call multifab_copy_c(Source_new(n),1,Source_old(n),1,1)
     end do

     if (do_sponge) then
        call make_sponge(sponge,dx,dt,mla)
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

           gamma1bar_hold = gamma1bar

           call advance_timestep(init_mode,mla,uold,sold,unew,snew,gpres,pres,normal, &
                                 rho0_old,rhoh0_old,rho0_new,rhoh0_new,p0_old,p0_new, &
                                 tempbar,gamma1bar,w0,rho_omegadot2,rho_Hnuc2, &
                                 div_coeff_old,div_coeff_new,grav_cell,dx,time,dt,dtold, &
                                 the_bc_tower,dSdt,Source_old,Source_new,etarho_ec, &
                                 etarho_cc,psi,sponge,hgrhs)

           gamma1bar = gamma1bar_hold

        end do ! end do istep_init_iter = 1,init_iter

     end if ! end if (init_iter > 0)

     if ( chk_int > 0 ) then

        !------------------------------------------------------------------------
        ! write a checkpoint file
        !------------------------------------------------------------------------

        allocate(chkdata(nlevs))
        do n = 1,nlevs
           call multifab_build(chkdata(n), mla%la(n), 2*dm+nscal, 0)
           call multifab_copy_c(chkdata(n),1                ,uold(n), 1,dm)
           call multifab_copy_c(chkdata(n),rho_comp+dm      ,sold(n), 1,nscal)
           call multifab_copy_c(chkdata(n),rho_comp+dm+nscal,gpres(n),1,dm)
        end do

        write(unit=check_index,fmt='(i5.5)') istep
        check_file_name = trim(check_base_name) // check_index

        call checkpoint_write(check_file_name, chkdata, &
                              pres, dSdt, Source_old, &
                              rho_omegadot2, rho_Hnuc2, mla%mba%rr, &
                              time, dt)

        last_chk_written = istep

        write(unit=base_state_name,fmt='("model_",i5.5)') istep
        write(unit=base_w0_name,fmt='("w0_",i5.5)') istep
        write(unit=base_etarho_name,fmt='("eta_",i5.5)') istep
        call write_base_state(base_state_name, base_w0_name, &
                              base_etarho_name, check_file_name, &
                              rho0_old, rhoh0_old, p0_old, gamma1bar, &
                              w0, etarho_ec, etarho_cc, &
                              div_coeff_old, psi, prob_lo(dm))

        do n = 1,nlevs
           call destroy(chkdata(n))
        end do
        deallocate(chkdata)
     end if

     if ( plot_int > 0 .or. plot_deltat > 0.0) then

        !------------------------------------------------------------------------
        ! write a plotfile
        !------------------------------------------------------------------------

        ! tempbar is only used as an initial guess for eos calls
        if (enthalpy_pred_type .ne. predict_Tprime_then_h) then
           call average(mla,sold,tempbar,dx,temp_comp)
           tempbar = max(tempbar,small_temp)
        end if

        write(unit=plot_index,fmt='(i5.5)') istep
        plot_file_name = trim(plot_base_name) // plot_index
        call make_plotfile(plot_file_name,mla,uold,sold,gpres,rho_omegadot2,rho_Hnuc2, &
                           Source_new,sponge,mla%mba,plot_names,time,dx,the_bc_tower,w0, &
                           rho0_old,rhoh0_old,p0_old,tempbar,gamma1bar,div_coeff_old,normal)

        call write_job_info(plot_file_name)
        last_plt_written = istep
     end if

  end if ! end if (restart < 0)

  !------------------------------------------------------------------------
  ! Main evolution loop
  !------------------------------------------------------------------------

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

  if ( (max_step >= init_step) .and. (time < stop_time .or. stop_time < 0.d0) ) then

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
           call print(aveassoc_mem_stats(),    "    aveassoc")
           if ( parallel_IOProcessor() ) print*, ''

           do n = 1,nlevs
              smax = norm_inf(uold(n),1,1)
              if ( parallel_IOProcessor()) then
                 print *,'MAX OF UOLD ', smax,' AT LEVEL ',n
              end if
                   
              smax = norm_inf(uold(n),2,1)
              if ( parallel_IOProcessor()) then
                 print *,'MAX OF VOLD ', smax,' AT LEVEL ',n
              end if
              
              if (dm > 2) then
                 smax = norm_inf(uold(n),3,1)
                 if ( parallel_IOProcessor()) then
                    print *,'MAX OF WOLD ', smax,' AT LEVEL ',n
                 end if
              end if
           end do
        end if

        !---------------------------------------------------------------------
        ! regrid
        !---------------------------------------------------------------------
        if (max_levs > 1 .and. regrid_int > 0 .and. (mod(istep-1,regrid_int) .eq. 0) ) then

           ! we do not regrid spherical base state arrays since there is only one
           ! level of refinement
           if (spherical .eq. 0) then
              
              ! first 'regrid' the 1d arrays 
              ! (other than rho0, rhoh0, and p0, which are handled below)
              ! copy the coarsest level only, interpolate to all the other levels
              ! and then copy the valid data from the old arrays onto the new
              ! this must be done before we call init_multilevel or else we lose track
              ! of where the old valid data was
              
              ! copy the coarsest level of the real arrays into the temp arrays
              psi_temp(1,:)        = psi(1,:)
              etarho_cc_temp(1,:)  = etarho_cc(1,:)
              etarho_ec_temp(1,:)  = etarho_ec(1,:)
              w0_temp(1,:)         = w0(1,:)

              ! piecewise constant interpolation to fill the cc temp arrays
              do n=2,max_levs
                 do r=0,nr(n)-1
                    psi_temp(n,r)        = psi_temp(n-1,r/2)
                    etarho_cc_temp(n,r)  = etarho_cc_temp(n-1,r/2)
                 end do
              end do

              ! piecewise linear interpolation to fill the edge-centered temp arrays
              do n=2,max_levs
                 do r=0,nr(n)
                    if (mod(r,2) .eq. 0) then
                       etarho_ec_temp(n,r) = etarho_ec_temp(n-1,r/2)
                       w0_temp(n,r)        = w0_temp(n-1,r/2)
                    else
                       etarho_ec_temp(n,r) = &
                            HALF*(etarho_ec_temp(n-1,r/2)+etarho_ec_temp(n-1,r/2+1))
                       w0_temp(n,r) = HALF*(w0_temp(n-1,r/2)+w0_temp(n-1,r/2+1))
                    end if
                 end do
              end do

              ! copy valid data into temp
              do n=2,nlevs_radial
                 do i=1,numdisjointchunks(n)
                    do r=r_start_coord(n,i),r_end_coord(n,i)
                       psi_temp(n,r)        = psi(n,r)
                       etarho_cc_temp(n,r)  = etarho_cc(n,r)
                    end do
                 end do
              end do
              do n=2,nlevs_radial
                 do i=1,numdisjointchunks(n)
                    do r=r_start_coord(n,i),r_end_coord(n,i)+1
                       etarho_ec_temp(n,r) = etarho_ec(n,r)
                       w0_temp(n,r)        = w0(n,r)
                    end do
                 end do
              end do

              ! copy temp array back into the real thing
              psi        = psi_temp
              etarho_cc  = etarho_cc_temp
              etarho_ec  = etarho_ec_temp
              w0         = w0_temp

           end if ! end regridding of base state

           ! create new grids and fill in data on those grids
           call regrid(mla,uold,sold,gpres,pres,dSdt,Source_old,rho_omegadot2,rho_Hnuc2, &
                       dx,the_bc_tower)

           call init_multilevel(sold)

           do n=1,nlevs
              call multifab_destroy(unew(n))
              call multifab_destroy(snew(n))
              call multifab_destroy(sponge(n))
              call multifab_destroy(hgrhs(n))
              call multifab_destroy(Source_new(n))
              if (dm .eq. 3) then
                 call multifab_destroy(normal(n))
              end if
           end do

           do n = 1,nlevs
              call multifab_build(      unew(n), mla%la(n),    dm, 3)
              call multifab_build(      snew(n), mla%la(n), nscal, 3)
              call multifab_build(    sponge(n), mla%la(n),     1, 0)
              call multifab_build(     hgrhs(n), mla%la(n),     1, 0, nodal)
              call multifab_build(Source_new(n), mla%la(n),     1, 1)
              if (dm .eq. 3) then
                 call multifab_build(normal(n), mla%la(n),    dm, 1)
              end if
              
              call setval(      unew(n), ZERO, all=.true.)
              call setval(      snew(n), ZERO, all=.true.)
              call setval(    sponge(n), ONE,  all=.true.)
              call setval(     hgrhs(n), ZERO, all=.true.)
              call setval(Source_new(n), ZERO, all=.true.)
           end do

           ! Create normal now that we have defined center and dx
           call make_normal(normal,dx)

           if (spherical .eq. 0) then
              call average(mla,sold,rho0_old,dx,rho_comp)
           end if
           
           call make_grav_cell(grav_cell,rho0_old)

           if (spherical .eq. 0) then
              ! enforce HSE
              call enforce_HSE(rho0_old,p0_old,grav_cell)
           end if

           ! compute full state T,h = T,h(rho,p0_old,X)
           call makeTHfromRhoP(sold,p0_old,the_bc_tower%bc_tower_array,mla,dx)

           if (spherical .eq. 0) then
              ! force rhoh0 to be the average of rhoh
              call average(mla,sold,rhoh0_old,dx,rhoh_comp)
           end if

           ! now that the grids and the full state have changed, recompute tempbar
           call average(mla,sold,tempbar,dx,temp_comp)
           tempbar = max(tempbar,small_temp)

           ! gamma1bar needs to be recomputed
           if (allocated(gamma1)) then
              deallocate(gamma1)
           end if
           allocate(gamma1(nlevs))
           
           do n=1,nlevs
              call multifab_build(gamma1(n), mla%la(n), 1, 0)
           end do
           
           call make_gamma(mla,gamma1,sold,p0_old,tempbar,dx)
           call average(mla,gamma1,gamma1bar,dx,1)
           
           do n=1,nlevs
              call destroy(gamma1(n))
           end do

           ! div_coeff_old needs to be recomputed
           call make_div_coeff(div_coeff_old,rho0_old,p0_old,gamma1bar,grav_cell)

        end if ! end regridding

        !---------------------------------------------------------------------
        ! get the new timestep
        !---------------------------------------------------------------------
        dtold = dt

        call make_grav_cell(grav_cell,rho0_old)

        if (istep > 1) then

           dt = 1.e20

           do n=1,nlevs
              call multifab_build(vel_force(n), mla%la(n), dm, 1)
              call setval( vel_force(n),ZERO, all=.true.)
           end do

           call mk_vel_force(vel_force,uold,gpres,sold,normal,rho0_old, &
                             grav_cell,dx,the_bc_tower%bc_tower_array,mla)

           allocate(w0mac(nlevs,dm))
           if (spherical .eq. 1) then
              do n=1,nlevs
                 do comp=1,dm
                    call multifab_build(w0mac(n,comp),mla%la(n),1,1, &
                                        nodal=edge_nodal_flag(comp,:))
                    call setval(w0mac(n,comp), ZERO, all=.true.)
                 end do
              end do
              call put_w0_on_edges(mla,w0,w0mac,dx,div_coeff_old,the_bc_tower)
           end if

           call estdt(uold,sold,vel_force,Source_old,dSdt,normal,w0,w0mac,p0_old, &
                      gamma1bar,dx,cflfac,dt)

           if (spherical .eq. 1) then
              do n=1,nlevs
                 do comp=1,dm
                    call destroy(w0mac(n,comp))
                 end do
              end do
           end if
           deallocate(w0mac)

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
              if (parallel_IOProcessor()) then
                 print*, "Stop time limits dt =",dt
              end if
           end if
        end if

        !---------------------------------------------------------------------
        ! Advance a single timestep at all levels.
        !---------------------------------------------------------------------
        init_mode = .false.
        if (do_sponge) then
           call init_sponge(rho0_old(1,:),dx(1,:),prob_lo(dm))
           call make_sponge(sponge,dx,dt,mla)
        end if
        runtime1 = parallel_wtime()

        call advance_timestep(init_mode,mla,uold,sold,unew,snew,gpres,pres,normal,rho0_old, &
                              rhoh0_old,rho0_new,rhoh0_new,p0_old,p0_new,tempbar,gamma1bar, &
                              w0,rho_omegadot2,rho_Hnuc2,div_coeff_old,div_coeff_new, &
                              grav_cell,dx,time,dt,dtold,the_bc_tower,dSdt,Source_old, &
                              Source_new,etarho_ec,etarho_cc,psi,sponge,hgrhs)

        runtime2 = parallel_wtime() - runtime1
        call parallel_reduce(runtime1, runtime2, MPI_MAX, proc = parallel_IOProcessorNode())
        if (parallel_IOProcessor()) print*,'Time to advance timestep: ',runtime1,' seconds'

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
           call print(aveassoc_mem_stats(),    "    aveassoc")
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
           call multifab_copy_c(uold(n),      1,unew(n),      1,dm,   3)
           call multifab_copy_c(sold(n),      1,snew(n),      1,nscal,3)
           call multifab_copy_c(Source_old(n),1,Source_new(n),1,1)
        end do

        ! Set div_coeff_old equal to div_coeff_new from the last time step
        div_coeff_old = div_coeff_new

        ! Copy the base state
        rho0_old = rho0_new
        rhoh0_old = rhoh0_new
        p0_old = p0_new

        !---------------------------------------------------------------------
        ! output
        !---------------------------------------------------------------------

        if (chk_int > 0) then
           if (mod(istep,chk_int) .eq. 0) then
              allocate(chkdata(nlevs))
              do n = 1,nlevs
                 call multifab_build(chkdata(n), mla%la(n), 2*dm+nscal, 0)
                 call multifab_copy_c(chkdata(n),1,unew(n),1,dm)
                 call multifab_copy_c(chkdata(n),rho_comp+dm,snew(n),1,nscal)
                 call multifab_copy_c(chkdata(n),rho_comp+dm+nscal,gpres(n),1,dm)
              end do

              write(unit=check_index,fmt='(i5.5)') istep
              check_file_name = trim(check_base_name) // check_index

              call checkpoint_write(check_file_name, chkdata, &
                                    pres, dSdt, Source_old, &
                                    rho_omegadot2, rho_Hnuc2, mla%mba%rr, &
                                    time, dt)

              last_chk_written = istep

              write(unit=base_state_name,fmt='("model_",i5.5)') istep
              write(unit=base_w0_name,fmt='("w0_",i5.5)') istep
              write(unit=base_etarho_name,fmt='("eta_",i5.5)') istep
              call write_base_state(base_state_name, base_w0_name, &
                                    base_etarho_name, check_file_name, &
                                    rho0_new, rhoh0_new, p0_new, gamma1bar(:,:), &
                                    w0, etarho_ec, etarho_cc, &
                                    div_coeff_old, psi, prob_lo(dm))

              do n = 1,nlevs
                 call destroy(chkdata(n))
              end do
              deallocate(chkdata)

           end if
        end if

        if (plot_int > 0 .or. plot_deltat > 0.0) then
           if ( (plot_int > 0 .and. mod(istep,plot_int) .eq. 0) .or. &
                (plot_deltat > 0.0 .and. &
                mod(time - dt,plot_deltat) > mod(time,plot_deltat))) then
              write(unit=plot_index,fmt='(i5.5)') istep
              plot_file_name = trim(plot_base_name) // plot_index
              call make_plotfile(plot_file_name,mla,unew,snew,gpres,rho_omegadot2,rho_Hnuc2,&
                                 Source_new,sponge,mla%mba,plot_names,time,dx,the_bc_tower, &
                                 w0,rho0_new,rhoh0_new,p0_new,tempbar,gamma1bar, &
                                 div_coeff_old,normal)

              call write_job_info(plot_file_name)
              last_plt_written = istep
           end if
        end if

        if (stop_time >= 0.d0) then
           if (time >= stop_time) goto 999
        end if

     end do

999  continue
     if (istep > max_step) istep = max_step

     !---------------------------------------------------------------------
     ! write the final checkpoint and plotfile
     !---------------------------------------------------------------------

1000 format('STEP = ',i5,1x,' TIME = ',f16.10,1x,'DT = ',f14.9)

     if ( chk_int > 0 .and. last_chk_written .ne. istep ) then
        !       This writes a checkpoint file.
        allocate(chkdata(nlevs))
        do n = 1,nlevs
           call multifab_build(chkdata(n), mla%la(n), 2*dm+nscal, 0)
           call multifab_copy_c(chkdata(n),1,unew(n),1,dm)
           call multifab_copy_c(chkdata(n),rho_comp+dm,snew(n),1,nscal)
           call multifab_copy_c(chkdata(n),rho_comp+dm+nscal,gpres(n),1,dm)
        end do

        write(unit=check_index,fmt='(i5.5)') istep
        check_file_name = trim(check_base_name) // check_index

        call checkpoint_write(check_file_name, chkdata, &
                              pres, dSdt, Source_old, &
                              rho_omegadot2, rho_Hnuc2, mla%mba%rr, &
                              time, dt)

        do n = 1,nlevs
           call destroy(chkdata(n))
        end do
        deallocate(chkdata)

        write(unit=base_state_name,fmt='("model_",i5.5)') istep
        write(unit=base_w0_name,fmt='("w0_",i5.5)') istep
        write(unit=base_etarho_name,fmt='("eta_",i5.5)') istep
        call write_base_state(base_state_name, base_w0_name, &
                              base_etarho_name, check_file_name, &
                              rho0_new, rhoh0_new, p0_new, gamma1bar, &
                              w0, etarho_ec, etarho_cc, &
                              div_coeff_old, psi, prob_lo(dm))
     end if

     if ( plot_int > 0 .and. last_plt_written .ne. istep ) then
        write(unit=plot_index,fmt='(i5.5)') istep
        plot_file_name = trim(plot_base_name) // plot_index
        call make_plotfile(plot_file_name,mla,unew,snew,gpres,rho_omegadot2,rho_Hnuc2, &
                           Source_new,sponge,mla%mba,plot_names,time,dx,the_bc_tower,w0, &
                           rho0_new,rhoh0_new,p0_new,tempbar,gamma1bar,div_coeff_old,normal)
        
        call write_job_info(plot_file_name)
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
     call destroy(rho_Hnuc2(n))
     call destroy(sponge(n))
  end do

  if(dm .eq. 3) then
     do n=1,nlevs
        call destroy(normal(n))
     end do
  end if

  call destroy(mla)

  call bc_tower_destroy(the_bc_tower)

  call destroy_geometry()

  call probin_close()

  deallocate(uold,sold,pres,gpres,dSdt,Source_old,rho_omegadot2,rho_Hnuc2,dx)
  deallocate(div_coeff_old,div_coeff_new,gamma1bar,gamma1bar_hold,s0_init,rho0_old)
  deallocate(rhoh0_old,rho0_new,rhoh0_new,p0_init,p0_old,p0_new,w0,etarho_ec,etarho_cc)
  deallocate(psi,tempbar,grav_cell)

  if (verbose .ge. 1) then
     if ( parallel_IOProcessor() ) then
        print *, 'MEMORY STATS AT END OF PROGRAM'
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
     call print(aveassoc_mem_stats(),    "    aveassoc")
  end if

end subroutine varden
