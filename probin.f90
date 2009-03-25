! This module stores the runtime parameters.  The probin_init() routine is
! used to initialize the runtime parameters

module probin_module

  use bl_types
  use bl_space
  use pred_parameters
  use multifab_module, only: multifab_set_alltoallv
  use cluster_module
  use layout_module

  implicit none

  character(len=256), save :: model_file
  real(dp_t), save         :: stop_time
  real(dp_t), save         :: prob_lo_x, prob_lo_y, prob_lo_z
  real(dp_t), save         :: prob_hi_x, prob_hi_y, prob_hi_z
  integer, save            :: max_step, plot_int
  real(dp_t), save         :: plot_deltat
  integer, save            :: chk_int, init_iter, init_divu_iter
  real(dp_t), save         :: cflfac, init_shrink
  character(len=128), save :: test_set
  integer, save            :: restart, do_initial_projection
  integer, save            :: bcx_lo, bcx_hi, bcy_lo, bcy_hi, bcz_lo, bcz_hi
  logical, save            :: pmask_x, pmask_y, pmask_z, pmask_xyz(MAX_SPACEDIM)
  integer, save            :: verbose, mg_verbose, cg_verbose
  integer, save            :: hg_bottom_solver, mg_bottom_solver
  logical, save            :: do_sponge, hg_dense_stencil
  real(dp_t), save         :: anelastic_cutoff, base_cutoff_density, dpdt_factor
  integer, save            :: spherical_in, dm_in
  logical, save            :: perturb_model, plot_spec, plot_trac, plot_base
  character(len=128), save :: plot_base_name, check_base_name
  logical, save            :: evolve_base_state
  logical, save            :: use_thermal_diffusion, temp_diffusion_formulation
  integer, save            :: thermal_diffusion_type
  logical, save            :: do_eos_h_above_cutoff
  integer, save            :: enthalpy_pred_type
  real(dp_t), save         :: max_dt_growth, fixed_dt
  real(dp_t), save         :: velpert_amplitude, velpert_radius, velpert_scale, velpert_steep
  logical, save            :: do_burning
  real(dp_t), save         :: grav_const
  real(dp_t), save         :: xrb_pert_factor, xrb_pert_size
  integer, save            :: xrb_pert_type
  logical, save            :: xrb_use_bottom_sponge
  logical, save            :: use_eos_coulomb
  real(dp_t), save         :: small_temp, small_dens
  real(dp_t), save         :: sponge_kappa
  logical, save            :: use_delta_gamma1_term, use_etarho
  integer, save            :: slope_order, interp_type_radial_bin_to_cart
  integer, save            :: w0mac_interp_type
  integer, save            :: nOutFiles
  logical, save            :: lUsingNFiles
  logical, save            :: use_tfromp, single_prec_plotfiles
  logical, save            :: use_soundspeed_firstdt, use_divu_firstdt
  logical, save            :: smallscale_beta, do_alltoallv, the_knapsack_verbosity
  integer, save            :: use_ppm
  integer, save            :: max_levs, max_grid_size, regrid_int, ref_ratio
  integer, save            :: n_cellx, n_celly, n_cellz
  integer, save            :: drdxfac, min_width, the_sfc_threshold
  real(dp_t), save         :: min_eff
  real(dp_t), save         :: burning_cutoff_density  ! note: presently not runtime parameter
  real(dp_t), save         :: rotational_frequency, co_latitude, radius
  character(len=256), save :: job_name

  ! These will be allocated and defined below
  logical,    allocatable, save :: edge_nodal_flag(:,:)
  logical,    allocatable, save :: nodal(:)
  logical,    allocatable, save :: pmask(:)
  real(dp_t), allocatable, save :: prob_lo(:)
  real(dp_t), allocatable, save :: prob_hi(:)

  namelist /probin/ model_file
  namelist /probin/ stop_time
  namelist /probin/ prob_lo_x
  namelist /probin/ prob_lo_y
  namelist /probin/ prob_lo_z
  namelist /probin/ prob_hi_x
  namelist /probin/ prob_hi_y
  namelist /probin/ prob_hi_z
  namelist /probin/ max_step
  namelist /probin/ plot_int
  namelist /probin/ plot_deltat
  namelist /probin/ chk_int
  namelist /probin/ init_iter
  namelist /probin/ init_divu_iter
  namelist /probin/ cflfac
  namelist /probin/ init_shrink
  namelist /probin/ test_set
  namelist /probin/ restart
  namelist /probin/ do_initial_projection
  namelist /probin/ bcx_lo
  namelist /probin/ bcx_hi
  namelist /probin/ bcy_lo
  namelist /probin/ bcy_hi
  namelist /probin/ bcz_lo
  namelist /probin/ bcz_hi
  namelist /probin/ pmask_x
  namelist /probin/ pmask_y
  namelist /probin/ pmask_z
  namelist /probin/ pmask_xyz
  namelist /probin/ verbose
  namelist /probin/ mg_verbose
  namelist /probin/ cg_verbose
  namelist /probin/ hg_bottom_solver
  namelist /probin/ mg_bottom_solver
  namelist /probin/ do_sponge
  namelist /probin/ hg_dense_stencil
  namelist /probin/ anelastic_cutoff
  namelist /probin/ base_cutoff_density
  namelist /probin/ dpdt_factor
  namelist /probin/ spherical_in
  namelist /probin/ dm_in
  namelist /probin/ perturb_model
  namelist /probin/ plot_spec
  namelist /probin/ plot_trac
  namelist /probin/ plot_base
  namelist /probin/ plot_base_name
  namelist /probin/ check_base_name
  namelist /probin/ evolve_base_state
  namelist /probin/ use_thermal_diffusion
  namelist /probin/ temp_diffusion_formulation
  namelist /probin/ thermal_diffusion_type
  namelist /probin/ do_eos_h_above_cutoff
  namelist /probin/ enthalpy_pred_type
  namelist /probin/ max_dt_growth
  namelist /probin/ fixed_dt
  namelist /probin/ velpert_amplitude
  namelist /probin/ velpert_radius
  namelist /probin/ velpert_scale
  namelist /probin/ velpert_steep
  namelist /probin/ do_burning
  namelist /probin/ grav_const
  namelist /probin/ xrb_pert_factor
  namelist /probin/ xrb_pert_size
  namelist /probin/ xrb_pert_type
  namelist /probin/ xrb_use_bottom_sponge
  namelist /probin/ use_eos_coulomb
  namelist /probin/ small_temp
  namelist /probin/ small_dens
  namelist /probin/ sponge_kappa
  namelist /probin/ use_delta_gamma1_term
  namelist /probin/ use_etarho
  namelist /probin/ slope_order
  namelist /probin/ interp_type_radial_bin_to_cart
  namelist /probin/ w0mac_interp_type
  namelist /probin/ nOutFiles
  namelist /probin/ lUsingNFiles
  namelist /probin/ use_tfromp
  namelist /probin/ single_prec_plotfiles
  namelist /probin/ use_soundspeed_firstdt
  namelist /probin/ use_divu_firstdt
  namelist /probin/ smallscale_beta
  namelist /probin/ do_alltoallv
  namelist /probin/ the_knapsack_verbosity
  namelist /probin/ use_ppm
  namelist /probin/ max_levs
  namelist /probin/ max_grid_size
  namelist /probin/ regrid_int
  namelist /probin/ ref_ratio
  namelist /probin/ n_cellx
  namelist /probin/ n_celly
  namelist /probin/ n_cellz
  namelist /probin/ drdxfac
  namelist /probin/ the_sfc_threshold
  namelist /probin/ min_eff
  namelist /probin/ min_width
  namelist /probin/ rotational_frequency
  namelist /probin/ co_latitude
  namelist /probin/ radius
  namelist /probin/ job_name

contains

  subroutine probin_init()

    use f2kcli
    use parallel
    use bc_module
    use bl_IO_module
    use bl_prof_module
    use bl_error_module
    use bl_constants_module
    use knapsack_module
    
    integer    :: narg, farg

    character(len=128) :: fname
    character(len=128) :: probin_env

    logical    :: lexist, need_inputs
    integer    :: i, natonce, myproc, nprocs, nsets, myset, iset, ibuff(1)
    integer    :: wakeuppid, waitforpid, tag, un, ierr
    real(dp_t) :: pistart, piend, pitotal, pistartall, piendall, pitotalall
    real(dp_t) :: pitotal_max, pitotalall_max

    type(bl_prof_timer), save :: bpt

    call build(bpt, "probin_init")

    narg = command_argument_count()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize the runtime parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Defaults
    model_file      = "model.hse"

    spherical_in = 0
    dm_in = 2

    max_step  = 1
    stop_time = -1.d0

    init_iter = 4
    init_divu_iter = 4
    plot_int  = 0
    plot_deltat = 0.0
    chk_int  = 0

    prob_lo_x = ZERO
    prob_lo_y = ZERO
    prob_lo_z = ZERO

    verbose = 0
    mg_verbose = 0
    cg_verbose = 0

    hg_bottom_solver = -1 ! valid values are >= 0
    mg_bottom_solver = -1 ! valid values are >= 0

    do_sponge = .false. 
    hg_dense_stencil = .true.
    anelastic_cutoff = 3.e6
    base_cutoff_density = 3.e6

    dpdt_factor = 0.0d0

    init_shrink = 1.0

    do_initial_projection  = 1

    need_inputs = .true.
    test_set = ''
    restart  = -1
  
    bcx_lo = SLIP_WALL
    bcy_lo = SLIP_WALL
    bcz_lo = SLIP_WALL
    bcx_hi = SLIP_WALL
    bcy_hi = SLIP_WALL
    bcz_hi = SLIP_WALL

    pmask_x = .false.
    pmask_y = .false.
    pmask_z = .false.
  
    perturb_model = .false.
    plot_spec     = .true.
    plot_trac     = .true.
    plot_base     = .false.
    plot_base_name = "plt"
    check_base_name = "chk"
    evolve_base_state = .true.
    use_thermal_diffusion = .false.
    temp_diffusion_formulation = .false.

    ! 1 = Crank-Nicolson, 2 = Backward Euler
    thermal_diffusion_type = 1

    do_eos_h_above_cutoff = .true.

    ! 1 = (rho h)', 2 = h, 3 = T into (rho h)prime, 4 = T into h
    enthalpy_pred_type = predict_rhohprime

    max_dt_growth = 1.1D0

    fixed_dt = -1.0d0

    slope_order = 4

    ! 1 = piecewise constant, 2 = piecewise linear
    interp_type_radial_bin_to_cart = 1

    ! 0.  Compute w0 at edges using a projection
    ! 1.  Interpolate w0 to cell centers, then average to edges
    ! 2.  Interpolate w0 to edges directly
    ! 3.  Interpolate w0 to nodes, then average to edges
    w0mac_interp_type = 1

    velpert_amplitude = 0.d0
    velpert_radius = 0.75d8
    velpert_scale = 0.8d8
    velpert_steep = 1.d0

    do_burning = .true.

    grav_const = -1.5d10

    xrb_pert_factor = 1.d-2
    xrb_pert_size = 5.0d1
    xrb_pert_type = 1
    xrb_use_bottom_sponge = .true.

    use_eos_coulomb = .false.

    small_temp = 5.d6
    small_dens = 1.d-5

    sponge_kappa = 10.d0

    use_delta_gamma1_term = .false.
    use_etarho = .true.

    nOutFiles    = 64
    lUsingNFiles = .true.

    use_tfromp = .false.

    single_prec_plotfiles = .false.

    use_soundspeed_firstdt = .false.
    use_divu_firstdt = .false.

    smallscale_beta = .false.

    do_alltoallv = .false.

    the_knapsack_verbosity = .false.

    the_sfc_threshold = 4  ! Same as the default in layout_module.

    ! 0 = no ppm
    ! 1 = 1985 ppm
    ! 2 = 2009 ppm
    use_ppm = 0

    max_levs = 1
    max_grid_size = 64
    regrid_int = -1
    ref_ratio = 2

    drdxfac = 1

    min_width = 16
    min_eff = 0.7d0

    rotational_frequency = ZERO
    co_latitude = ZERO
    radius = 1.0e6_dp_t    ! 10 km neutron star

    job_name = ""

    !
    ! Don't have more than 64 processes trying to read from disk at once.
    !
    natonce = min(64,parallel_nprocs())
    myproc  = parallel_myproc()
    nprocs  = parallel_nprocs()
    nsets   = ((nprocs + (natonce - 1)) / natonce)
    myset   = (myproc / natonce)

    pistartall = parallel_wtime()

    do iset = 0, nsets-1

       if (myset .eq. iset) then

          pistart = parallel_wtime()
          
          call get_environment_variable('PROBIN', probin_env, status = ierr)
          if ( need_inputs .AND. ierr == 0 ) then
             un = unit_new()
             open(unit=un, file = probin_env, status = 'old', action = 'read')
             read(unit=un, nml = probin)
             close(unit=un)
             need_inputs = .false.
          end if

          farg = 1
          if ( need_inputs .AND. narg >= 1 ) then
             call get_command_argument(farg, value = fname)
             inquire(file = fname, exist = lexist )
             if ( lexist ) then
                farg = farg + 1
                un = unit_new()
                open(unit=un, file = fname, status = 'old', action = 'read')
                read(unit=un, nml = probin)
                close(unit=un)
                need_inputs = .false.
             end if
          end if

          inquire(file = 'inputs_varden', exist = lexist)
          if ( need_inputs .AND. lexist ) then
             un = unit_new()
             open(unit=un, file = 'inputs_varden', status = 'old', action = 'read')
             read(unit=un, nml = probin)
             close(unit=un)
             need_inputs = .false.
          end if

          piend = parallel_wtime()

          ibuff(1)  = 0
          wakeuppid = myproc + natonce
          tag       = mod(myproc,natonce)
          
          if (wakeuppid < nprocs) call parallel_send(ibuff, wakeuppid, tag)

       end if

      if (myset .eq. (iset + 1)) then

         tag        = mod(myproc,natonce)
         waitforpid = myproc - natonce

         call parallel_recv(ibuff, waitforpid, tag)
      endif

    end do

    piendall   = parallel_wtime()
    pitotal    = piend - pistart
    pitotalall = piendall - pistartall

    call parallel_reduce(pitotal_max,    pitotal,    MPI_MAX, &
                         proc = parallel_IOProcessorNode())
    call parallel_reduce(pitotalall_max, pitotalall, MPI_MAX, &
                         proc = parallel_IOProcessorNode())

    if (parallel_IOProcessor()) then
      print*, "PROBINIT max time   = ", pitotal_max
      print*, "PROBINIT total time = ", pitotalall_max
    endif

    pmask_xyz = (/pmask_x, pmask_y, pmask_z/)
    
    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--model_file')
          farg = farg + 1
          call get_command_argument(farg, value = model_file)

       case ('--stop_time')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stop_time

       case ('--prob_lo_x')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_lo_x
       case ('--prob_lo_y')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_lo_y
       case ('--prob_lo_z')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_lo_z
       case ('--prob_hi_x')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_hi_x
       case ('--prob_hi_y')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_hi_y
       case ('--prob_hi_z')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) prob_hi_z

       case ('--max_step')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_step

       case ('--plot_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) plot_int

       case ('--plot_deltat')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) plot_deltat

       case ('--chk_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) chk_int

       case ('--init_iter')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) init_iter

       case ('--init_divu_iter')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) init_divu_iter

       case ('--cfl')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) cflfac

       case ('--init_shrink')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) init_shrink

       case ('--test_set')
          farg = farg + 1
          call get_command_argument(farg, value = test_set)

       case ('--restart')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) restart

       case ('--do_initial_projection')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) do_initial_projection

       case ('--bcx_lo')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) bcx_lo
       case ('--bcx_hi')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) bcx_hi
       case ('--bcy_lo')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) bcy_lo
       case ('--bcy_hi')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) bcy_hi
       case ('--bcz_lo')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) bcz_lo
       case ('--bcz_hi')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) bcz_hi

       case ('--pmask_x')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) pmask_xyz(1)
       case ('--pmask_y')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) pmask_xyz(2)
       case ('--pmask_z')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) pmask_xyz(3)
       case ('--pmask_xyz')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) pmask_xyz(1)
          pmask_xyz = pmask_xyz(1)

       case ('--verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) verbose

       case ('--mg_verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_verbose

       case ('--cg_verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) cg_verbose

       case ('--hg_bottom_solver')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) hg_bottom_solver

       case ('--mg_bottom_solver')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_bottom_solver

       case ('--do_sponge')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) do_sponge

       case ('--hg_dense_stencil')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) hg_dense_stencil

       case ('--anelastic_cutoff')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) anelastic_cutoff

       case ('--base_cutoff_density')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) base_cutoff_density

       case ('--dpdt_factor')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) dpdt_factor

       case ('--spherical_in')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) spherical_in

       case ('--dm_in')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) dm_in

       case ('--perturb_model')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) perturb_model

       case ('--plot_spec')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) plot_spec

       case ('--plot_trac')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) plot_trac

       case ('--plot_base')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) plot_base

       case ('--plot_base_name')
          farg = farg + 1
          call get_command_argument(farg, value = plot_base_name)

       case ('--check_base_name')
          farg = farg + 1
          call get_command_argument(farg, value = check_base_name)

       case ('--evolve_base_state')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) evolve_base_state

       case ('--use_thermal_diffusion')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_thermal_diffusion

       case ('--temp_diffusion_formulation')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) temp_diffusion_formulation

       case ('--thermal_diffusion_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) thermal_diffusion_type

       case ('--do_eos_h_above_cutoff')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) do_eos_h_above_cutoff

       case ('--enthalpy_pred_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) enthalpy_pred_type

       case ('--max_dt_growth')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_dt_growth

       case ('--fixed_dt')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) fixed_dt

       case ('--velpert_amplitude')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) velpert_amplitude

       case ('--velpert_radius')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) velpert_radius

       case ('--velpert_scale')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) velpert_scale

       case ('--velpert_steep')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) velpert_steep

       case ('--do_burning')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) do_burning

       case ('--grav_const')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) grav_const

       case ('--xrb_pert_factor')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) xrb_pert_factor

       case ('--xrb_pert_size')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) xrb_pert_size

       case ('--xrb_pert_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) xrb_pert_type

       case ('--xrb_use_bottom_sponge')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) xrb_use_bottom_sponge

       case ('--use_eos_coulomb')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_eos_coulomb

       case ('--small_temp')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) small_temp

       case ('--small_dens')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) small_dens

       case ('--sponge_kappa')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) sponge_kappa

       case ('--use_delta_gamma1_term')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_delta_gamma1_term

       case ('--use_etarho')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_etarho

       case ('--slope_order')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) slope_order

       case ('--interp_type_radial_bin_to_cart')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) interp_type_radial_bin_to_cart

       case ('--w0mac_interp_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) w0mac_interp_type

       case ('--nOutFiles')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) nOutFiles

       case ('--lUsingNFiles')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) lUsingNfiles

       case ('--use_tfromp')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_tfromp

       case ('--use_single_prec_plotfiles')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) single_prec_plotfiles

       case ('--use_soundspeed_firstdt')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_soundspeed_firstdt

       case ('--use_divu_firstdt')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_divu_firstdt

       case ('--smallscale_beta')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) smallscale_beta

       case ('--do_alltoallv')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) do_alltoallv

       case ('--the_knapsack_verbosity')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) the_knapsack_verbosity

       case ('--the_sfc_threshold')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) the_sfc_threshold

       case ('--use_ppm')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_ppm

       case ('--max_levs')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_levs

       case ('--max_grid_size')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) max_grid_size

       case ('--regrid_int')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) regrid_int

       case ('--ref_ratio')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) ref_ratio

       case ('--n_cellx')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) n_cellx
       case ('--n_celly')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) n_celly
       case ('--n_cellz')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) n_cellz

       case ('--drdxfac')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) drdxfac
       case ('--min_width')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) min_width

       case ('--rotational_frequency')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) rotational_frequency
       case ('--co_latitude')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) co_latitude
       case ('--radius')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) radius

       case ('--job_name')
          farg = farg + 1
          call get_command_argument(farg, value = job_name)

       case ('--')
          farg = farg + 1
          exit

       case default
          if ( .not. parallel_q() ) then
             write(*,*) 'UNKNOWN option = ', fname
             call bl_error("MAIN")
          end if
       end select

       farg = farg + 1
    end do

    ! for the moment, set the cutoff for burning to be base_cutoff_density
    burning_cutoff_density = base_cutoff_density

    ! initialize edge_nodal_flag
    allocate(edge_nodal_flag(dm_in,dm_in))
    edge_nodal_flag = .false.
    do i = 1,dm_in
       edge_nodal_flag(i,i) = .true.
    end do

    ! initialize nodal
    allocate(nodal(dm_in))
    nodal = .true.

    ! initialize pmask
    allocate(pmask(dm_in))
    pmask = .FALSE.
    pmask = pmask_xyz(1:dm_in)

    ! initialize prob_lo and prob_hi
    allocate(prob_lo(dm_in))
    prob_lo(1) = prob_lo_x
    if (dm_in > 1) prob_lo(2) = prob_lo_y
    if (dm_in > 2) prob_lo(3) = prob_lo_z
    allocate(prob_hi(dm_in))
    prob_hi(1) = prob_hi_x
    if (dm_in > 1) prob_hi(2) = prob_hi_y
    if (dm_in > 2) prob_hi(3) = prob_hi_z

    call cluster_set_min_eff(min_eff)
    call cluster_set_minwidth(min_width)

    if (do_alltoallv) call multifab_set_alltoallv(.true.)

    call knapsack_set_verbose(the_knapsack_verbosity)

    call layout_set_sfc_threshold(the_sfc_threshold)

    call destroy(bpt)
    
  end subroutine probin_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine probin_close()

    deallocate(edge_nodal_flag)
    deallocate(nodal)
    deallocate(pmask)
    deallocate(prob_lo)
    deallocate(prob_hi)

  end subroutine probin_close

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module probin_module
