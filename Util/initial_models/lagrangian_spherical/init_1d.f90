!!  The model is placed into HSE by the following differencing:
!!
!!   (1/dr) [ <P>_i - <P>_{i-1} ] = (1/2) [ <rho>_i + <rho>_{i-1} ] g
!!
!!  This will be iterated over in tandem with the EOS call,
!!  P(i-1) = P_eos(rho(i-1), T(i-1), X(i-1)
!!
!!  When smooth = T then also a smoothing of stiff convective boundaries is done by 
!!  minimizing 
!!
!!  f = p_EOS(rho, T, X) - p_HSE (rho)
!!  q = S_EOS(rho, T, X) - S_CORE
!!
!!  where S_CORE is the entropy at the convective boundary
!!
!!  after that the now changed density and composition profile get smoothed out by a moving average
!!  this also makes a renormalization of the composition necessary
!!
!!
!!
!!
!!
!!
!!  When change_strat = T then also the density stratification is changed to garantee a constant brunt-vaisala frequency outwards. 
!!  At default this will change the stratification after the maximum of N^2 is reached, but a arbitrary point can be given as well
!!
!!  f = p_EOS(rho, T, X) - p_HSE (rho)
!!  q = N_EOS(rho, T, X) - N_max
!!



!!
!!  We will integrate upwards only, because thats what MAESTRO assumes to be in HSE. 
!!

program init_1d
 
  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module, only: eos_input_rt, eos, eos_init
  use eos_type_module, only: eos_t
  use extern_probin_module, only: use_eos_coulomb
  use network, only : nspec, network_species_index, spec_names, network_init
  use fundamental_constants_module, only: Gconst
  use model_parser_module
  use conductivity_module

  implicit none

  integer :: i, n, j

  character(len=128) :: params_file

  real (kind=dp_t), DIMENSION(nspec) :: xn

  real (kind=dp_t), allocatable :: xzn_hse(:), xznl_hse(:), xznr_hse(:)
  real (kind=dp_t), allocatable :: model_hse(:,:),temp_model(:,:)

  real :: A

  integer, parameter :: var_names_size=80
  
  character(len=var_names_size), allocatable :: varnames_stored(:)
  integer :: lun1, lun2, lun3
  
  integer :: nx
  real (kind=dp_t) :: xmin, xmax, dCoord, xmin_smooth, xmax_smooth
  real (kind=dp_t), allocatable :: brunt(:), conductivity(:)

  real (kind=dp_t) :: sumx, sumrho, sumn
  real (kind=dp_t), DIMENSION(nspec) :: sumxn
  real (kind=dp_t) :: coreX, frhoT, qrhoT 
  integer :: comp
 
  real (kind=dp_t) :: dens_zone, temp_zone, pres_zone, entropy, s_zone
  real (kind=dp_t) :: dpd, dpdt, dsdt, dsdrho
  real (kind=dp_t) :: dndrho, dndt
  real (kind=dp_t) :: prev_mu, prev_p, prev_temp
  
  
  
  real (kind=dp_t) :: p_want, drho, dtemp, delx, s_want 
  real (kind=dp_t) :: brunt_want, brunt_zone, brunt_slope, prev_brunt
  real (kind=dp_t) :: change_strat_r = -1.0
  
  real (kind=dp_t) :: g_zone, g_const
  logical :: do_invsq_grav
    
  real (kind=dp_t), parameter :: TOL = 1.e-12

  integer, parameter :: MAX_ITER = 250

  integer :: iter

  logical :: converged_hse, fluff, smooth, converged_smooth, converged_strat
  logical :: neutral_strat, change_strat,put_in_hse, use_slope

  real (kind=dp_t) :: low_density_cutoff, temp_cutoff, smallx, max_T
  real (kind=dp_t) :: model_shift

  integer :: index_base, conv_base, min_base, max_base, cent_base, smoothness, csmooth, norm
  integer :: brunt_base
  real (kind=dp_t) :: wnorm

  character (len=256) :: outfile, outfile2, outfile3
  character (len=8) :: num
  character (len=32) :: dxstr
  character (len=32) :: num_to_unitstring
  character (len=64) :: model_file

  real (kind=dp_t) :: max_hse_error, dpdr, rhog

  character (len=128) :: model_prefix

  integer :: narg

  type (eos_t) :: eos_state
  

  namelist /params/ nx, model_file, xmin, xmax, g_const, &
                    temp_cutoff, do_invsq_grav, use_slope, &
                    low_density_cutoff, model_prefix, model_shift, smooth, &
                    xmin_smooth, xmax_smooth, smoothness, norm, &
                    change_strat_r, change_strat, put_in_hse, neutral_strat, index_base

  ! determine if we specified a runtime parameters file or use the default      
  narg = command_argument_count()

  if (narg == 0) then
     params_file = "_params"
  else
     call get_command_argument(1, value = params_file)
  endif


  ! define defaults for the parameters for this model
  nx = 640 

  model_file = "model.in"

  model_shift = 0.0_dp_t

  xmin = 0.0_dp_t
  xmax = 2.e3_dp_t
  xmin_smooth = 1.0d10
  xmax_smooth = 2.0d11
  
  do_invsq_grav = .false.
  g_const = -1.0

  model_prefix = "model"

  low_density_cutoff = 1.d-4
  temp_cutoff = 1.d6

  smallx = 1.d-10
  smooth = .false.
  smoothness = 4
  norm = 512
  change_strat = .false.
  put_in_hse = .true.
  neutral_strat = .false.
  use_slope = .false.
  ! this comes in via extern_probin_module -- override the default
  ! here if we want
  use_eos_coulomb = .true.
  index_base = -1


  ! initialize the EOS and network
  call eos_init()
  call network_init()

  ! check the namelist for any changed parameters
  open(unit=11, file=params_file, status="old", action="read")
  read(unit=11, nml=params)
  close(unit=11)


  ! start by reading in the Lagrangian initial model
  call read_model_file(model_file)

  ! apply the shift
  do i = 1, npts_model
     model_r(i) = model_r(i) - model_shift
  enddo

  
  ! check if all abundances sum to 1 and normalize if they don't 
  do i = 1, npts_model
     sumx = ZERO
     do comp = 1, nspec
	sumx = sumx + model_state(i,ispec_model-1+comp)
     enddo
     if (sumx .ne. ONE) then
	do comp = 1, nspec
	  model_state(i,ispec_model-1+comp) = model_state(i,ispec_model-1+comp) / sumx
	enddo 
     endif
  enddo
  
!-----------------------------------------------------------------------------
! Create a 1-d uniform grid that is identical to the mesh that we are
! mapping onto, and then we want to force it into HSE on that mesh.
!-----------------------------------------------------------------------------

  ! allocate storage
  allocate(xzn_hse(nx))
  allocate(xznl_hse(nx))
  allocate(xznr_hse(nx))
  allocate(model_hse(nx,nvars_model))
  allocate(temp_model(nx,nvars_model))
  allocate(brunt(nx))
  allocate(conductivity(nx))


  ! compute the coordinates of the new gridded function
  dCoord = (xmax - xmin) / dble(nx)

  do i = 1, nx
     xznl_hse(i) = xmin + (dble(i) - ONE)*dCoord
     xzn_hse(i)  = xmin + (dble(i) - HALF)*dCoord
     xznr_hse(i) = xmin + (dble(i))*dCoord
  enddo
  

!-----------------------------------------------------------------------------
! put the model onto our new uniform grid and make it thermodynamically consistent
!-----------------------------------------------------------------------------

  do i = 1, nx
     do n = 1, nvars_model
        if (n == itemp_model) then
           model_hse(i,n) = max(temp_cutoff, interpolate(xzn_hse(i), n, interpolate_top=.true.))
        else
           model_hse(i,n) = interpolate(xzn_hse(i), n)
        endif
     enddo
     
     ! make it all thermodynamically consistent
     eos_state%rho = model_hse(i,idens_model)
     eos_state%T = model_hse(i,itemp_model)
     eos_state%xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)

     call eos(eos_input_rt, eos_state)

     model_hse(i,ipres_model) = eos_state%p
  enddo
  

   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!make an output of the initial (interpolated) profiles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !compute the brunt-vaisala frequency at each zone
  brunt(1) = 0
  
  dens_zone = model_hse(1,idens_model)
  temp_zone = max(temp_cutoff, model_hse(1,itemp_model))
  xn(:) = model_hse(1,ispec_model:ispec_model-1+nspec)
  
  eos_state%T = temp_zone
  eos_state%rho = dens_zone
  eos_state%xn(:) = xn(:)
    
  call eos(eos_input_rt, eos_state)
  
  prev_p = eos_state%p
  prev_mu = eos_state%abar
  prev_temp = temp_zone
  
  do i = 2,nx 
     delx = xzn_hse(i) - xzn_hse(i-1)
     ! compute the gravitation acceleration at the lower edge
     call calc_grav_zone(i,xznr_hse,xznl_hse,xzn_hse,model_hse(:,idens_model),.True.,do_invsq_grav,g_zone)     


    dens_zone = model_hse(i,idens_model)
    temp_zone = max(temp_cutoff, model_hse(i,itemp_model))
    xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)
    
    eos_state%T = temp_zone
    eos_state%rho = dens_zone
    eos_state%xn(:) = xn(:)
    
    call eos(eos_input_rt, eos_state)
    call calc_brunt(g_zone,delx,dens_zone,temp_zone,xn,prev_p,prev_temp,prev_mu,brunt(i))
    write(*,*) eos_state%dpdr
    prev_p = eos_state%p
    prev_mu = eos_state%abar
    prev_temp = temp_zone

  enddo
  
  
  !!write everything into a file
  outfile = trim(model_prefix) // ".initial_interpolated_profile"
  open (newunit=lun1, file=outfile, status="unknown")


  write (lun1,2002), 'Initial interpolated profile: '
  write (lun1,2001), 'radius','entropy','pressure','density','X(H1)', 'N^2'
  do i = 1, nx
    dens_zone = model_hse(i,idens_model)
    temp_zone = max(temp_cutoff, model_hse(i,itemp_model))
    xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)
    
    eos_state%T = temp_zone
    eos_state%rho = dens_zone
    eos_state%xn(:) = xn(:)
    
    call eos(eos_input_rt, eos_state)
    
    write (lun1,2000), xzn_hse(i), eos_state%s, eos_state%p, eos_state%rho, eos_state%xn(1), brunt(i)
  enddo


2000 format (1x, 100(g26.16, 1x))
2001 format (1x,100(a26, 1x))
2002 format (a)

  close (unit=lun1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determine all the bases necessary to do the following operations:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!
  !index base
  !!!!
  index_base = 1

  !!!!
  ! conv_base
  !!!!
  if (neutral_strat) then
    !Find the position of the convective boundary. In garstec, convective zones
    !are full mixed -> first zone with varying composition (X) is just outside convective zone
    coreX = model_hse(1,ispec_model)
    conv_base = 1
    do i = 2,nx 
      if (model_hse(i,ispec_model) .ne. coreX) then
        conv_base = i-1
        exit
      endif 
    enddo
  
    print *, 'Convective boundary at zone ', conv_base
  else
    conv_base = -1
  endif
  
  if (index_base .lt. conv_base) then
    index_base = conv_base
  endif
  
  !!!!
  ! strat_base
  !!!!
  if (change_strat) then
    !Find the radius of the maximum of the brunt-vaisala frequency, or use a given radius. 
    !From this radius outward we will keep the brunt-vaisala frequency constant while also conserving the HSE
    
    dens_zone = model_hse(1,idens_model)
    temp_zone = max(temp_cutoff, model_hse(1,itemp_model))
    xn(:) = model_hse(1,ispec_model:ispec_model-1+nspec)
    
    eos_state%T = temp_zone
    eos_state%rho = dens_zone
    eos_state%xn(:) = xn(:)
      
    call eos(eos_input_rt, eos_state)
    
    prev_p = eos_state%p
    prev_mu = eos_state%abar
    prev_temp = temp_zone
    

    brunt_want = -1d30    

    if (change_strat_r .ge. 0) then
      brunt_base = 2 
      do i = 2,nx
        if (xzn_hse(i) .ge. change_strat_r) then
          brunt_base = i
          exit
        endif
      enddo           
    else 
      brunt_base = 2 
      do i = 2,nx
        delx = xzn_hse(i) - xzn_hse(i-1)
        ! compute the gravitation acceleration at the lower edge
        call calc_grav_zone(i,xznr_hse,xznl_hse,xzn_hse,model_hse(:,idens_model),.True.,do_invsq_grav,g_zone)     


        dens_zone = model_hse(i,idens_model)
        temp_zone = max(temp_cutoff, model_hse(i,itemp_model))
        xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)
        
        eos_state%T = temp_zone
        eos_state%rho = dens_zone
        eos_state%xn(:) = xn(:)
        
        call eos(eos_input_rt, eos_state)
        call calc_brunt(g_zone,delx,dens_zone,temp_zone,xn,prev_p,prev_temp,prev_mu,brunt_zone)
        
        if (brunt_zone .ge. brunt_want) then
          brunt_want = brunt_zone
          brunt_base = i
        endif

        prev_p = eos_state%p
        prev_mu = eos_state%abar
        prev_temp = temp_zone
      enddo
    endif   
  else
    brunt_base = nx
  endif
  
  
  !!!!
  ! smoothing
  !!!!
  if (smooth) then
    min_base = 0
    max_base = 0
    do i = 1, nx
      if (xzn_hse(i) .ge. xmin_smooth) then
          min_base = i
          exit
      endif
    enddo
    do i = 1, nx
      if (xzn_hse(i) .ge. xmax_smooth) then
          max_base = i
          exit
      endif
    enddo 
    
    
    if (min_base .eq. 0) then
      min_base = 1
    endif 
    
    if (max_base .eq. 0) then
      max_base = nx
    endif 

    cent_base = min_base + (max_base - min_base)/2
    
    smoothness = (max_base-min_base)/smoothness
    if ((smoothness .gt. min_base-1 ) .or. (smoothness .gt. nx-max_base)) then
      smoothness = min(min_base-1,nx-max_base)
    endif
    
    print *, 'min base for smoothing is ', min_base 
    print *, 'max base for smoothing is ', max_base 
    print *, 'centering around ', cent_base 
    print *, 'we will smooth over ', smoothness, ' grid points'
  else
    min_base = -1
    max_base = -1
    cent_base = -1
  endif
  
  
!-----------------------------------------------------------------------------
!Create a neutrally stratified core and smooth the density and composition profile afterwards
!-----------------------------------------------------------------------------

  if (neutral_strat) then  
    
    !keep the entropy at the core constant and enforce hse by integrating inwards from the convective boundary
    temp_zone = model_hse(conv_base,itemp_model)
    dens_zone = model_hse(conv_base,idens_model)
    xn(:) = model_hse(conv_base,ispec_model:ispec_model-1+nspec)
    
    eos_state%T     = temp_zone
    eos_state%rho   = dens_zone
    eos_state%xn(:) = xn(:)

    call eos(eos_input_rt, eos_state)
    
    s_want = eos_state%s
    
    
    
    
    !---------------------------------------------------------------------------
    ! integrate down -- using the temperature profile defined above
    ! we do this to set the i=1 consistent 
    ! if we set i=1 directlty to s_want the change can be too big 
    !---------------------------------------------------------------------------
    do i = conv_base, 1, -1

     delx = xzn_hse(i+1) - xzn_hse(i)
     ! compute the gravitation acceleration at the upper edge
     call calc_grav_zone(i,xznr_hse,xznl_hse,xzn_hse,model_hse(:,idens_model),.False.,do_invsq_grav,g_zone)     


      ! we already set the temperature and composition profiles
      temp_zone = max(temp_cutoff, model_hse(i,itemp_model))
      ! use our previous initial guess for density
      dens_zone = model_hse(i,idens_model)
      xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)

      !-----------------------------------------------------------------------
      ! iteration loop
      !-----------------------------------------------------------------------

      ! start off the Newton loop by saying that the zone has not converged
      converged_smooth = .FALSE.

      do iter = 1, MAX_ITER

          ! get the pressure we want from the HSE equation, just the
          ! zone below the current.  Note, we are using an average of
          ! the density of the two zones as an approximation of the
          ! interface value -- this means that we need to iterate for
          ! find the density and pressure that are consistent
          
          ! HSE differencing
          p_want = model_hse(i+1,ipres_model) - &
              delx*0.5*(dens_zone + model_hse(i+1,idens_model))*g_zone

          
          ! we need to zero:
          !   frhoT = p(rho) - p_want
          !   qrhoT = s_eos - s_want 
          
          ! (t, rho) -> (p)
          eos_state%T     = temp_zone
          eos_state%rho   = dens_zone
          eos_state%xn(:) = xn(:)

          call eos(eos_input_rt, eos_state)
          
          pres_zone = eos_state%p
          s_zone = eos_state%s
          
          dpd = eos_state%dpdr
          dpdt = eos_state%dpdt
          dsdt = eos_state%dsdt
          dsdrho = eos_state%dsdr
          
          frhoT = pres_zone - p_want
          qrhoT = s_zone - s_want
          
          dtemp = (qrhoT - frhoT * dsdrho / (dpd + 0.5*delx*g_zone) ) / (dsdrho * dpdt /(dpd + 0.5*delx*g_zone) - dsdt)
          drho = -(frhoT + dpdt * dtemp)/(dpd + 0.5*delx*g_zone)
          
          dens_zone = max(0.9_dp_t*dens_zone, &
              min(dens_zone + drho, 1.1_dp_t*dens_zone))
          temp_zone = max (0.9_dp_t*temp_zone, &
              min(temp_zone + dtemp, 1.1_dp_t*temp_zone))
                          
          if (abs(drho) < TOL*dens_zone .and. abs(dtemp) < TOL*temp_zone) then
            converged_smooth = .TRUE.
            exit
          endif
          

      enddo
          
      if (.NOT. converged_smooth) then
          
          print *, 'Error zone', i, ' did not converge while creating a neutrally stratified core'
          print *, dens_zone, temp_zone
          print *, p_want
          print *, drho, dtemp
          call bl_error('Error: HSE non-convergence')
            
      endif


      ! call the EOS one more time for this zone and then go on to the next
      ! (t, rho) -> (p)
      eos_state%T     = temp_zone
      eos_state%rho   = dens_zone
      eos_state%xn(:) = xn(:)

      call eos(eos_input_rt, eos_state)

      pres_zone = eos_state%p
      
      ! update the thermodynamics in this zone
      model_hse(i,idens_model) = dens_zone
      model_hse(i,itemp_model) = temp_zone
      model_hse(i,ipres_model) = pres_zone
      
    enddo 
    
    
    
!     !---------------------------------------------------------------------------
!     ! integrate i=1 to s_want and set p to the thermodynamically correct value
!     ! once this is set we can move forward 
!     !---------------------------------------------------------------------------
!     converged_smooth = .FALSE.
!     
!     do iter = 1, MAX_ITER
!           ! (t, rho) -> (s)
!           temp_zone = model_hse(conv_base,itemp_model)
!           dens_zone = model_hse(conv_base,idens_model)
!           xn(:) = model_hse(conv_base,ispec_model:ispec_model-1+nspec)
!     
!           call eos(eos_input_rt, eos_state)
!           
!           s_zone = eos_state%s
!           dsdrho = eos_state%dsdr
!           
!           qrhoT = s_zone - s_want
!           drho = -qrhoT/dsdrho
!           
!           dens_zone = max(0.9_dp_t*dens_zone, &
!               min(dens_zone + drho, 1.1_dp_t*dens_zone))
!                           
!           if (abs(drho) < TOL*dens_zone) then
!             converged_smooth = .TRUE.
!             exit
!           endif          
!     enddo    
! 
!     if (.NOT. converged_smooth) then
!         
!         print *, 'Error: did not converge while creating the inner most point of a neutrally stratified core'
!         print *, dens_zone, temp_zone
!         print *, s_want
!         print *, drho, dtemp
!         call bl_error('Error: HSE non-convergence')
!     
!     endif
! 
!     ! call the EOS one more time for this zone and then go on to the next
!     ! (t, rho) -> (p)
!     eos_state%T     = temp_zone
!     eos_state%rho   = dens_zone
!     eos_state%xn(:) = xn(:)
! 
!     call eos(eos_input_rt, eos_state)
! 
!     pres_zone = eos_state%p
!     
!     ! update the thermodynamics in this zone
!     model_hse(1,idens_model) = dens_zone
!     model_hse(1,itemp_model) = temp_zone
!     model_hse(1,ipres_model) = pres_zone

    !---------------------------------------------------------------------------
    ! integrate up -- using the temperature profile defined below
    !---------------------------------------------------------------------------
    do i = 2, conv_base, 1

     delx = xzn_hse(i) - xzn_hse(i-1)
     ! compute the gravitation acceleration at the lower edge
     call calc_grav_zone(i,xznr_hse,xznl_hse,xzn_hse,model_hse(:,idens_model),.True.,do_invsq_grav,g_zone)     


      ! we already set the temperature and composition profiles
      temp_zone = max(temp_cutoff, model_hse(i,itemp_model))
      ! use our previous initial guess for density
      dens_zone = model_hse(i,idens_model)


      !-----------------------------------------------------------------------
      ! iteration loop
      !-----------------------------------------------------------------------

      ! start off the Newton loop by saying that the zone has not converged
      converged_smooth = .FALSE.

      do iter = 1, MAX_ITER

          ! get the pressure we want from the HSE equation, just the
          ! zone below the current.  Note, we are using an average of
          ! the density of the two zones as an approximation of the
          ! interface value -- this means that we need to iterate for
          ! find the density and pressure that are consistent
          
          ! HSE differencing
          p_want = model_hse(i-1,ipres_model) + &
              delx*0.5*(dens_zone + model_hse(i-1,idens_model))*g_zone

          
          ! we need to zero:
          !   frhoT = p(rho) - p_want
          !   qrhoT = s_eos - s_want 
          
          ! (t, rho) -> (p)
          eos_state%T     = temp_zone
          eos_state%rho   = dens_zone
          eos_state%xn(:) = xn(:)

          call eos(eos_input_rt, eos_state)
          
          pres_zone = eos_state%p
          s_zone = eos_state%s
          
          dpd = eos_state%dpdr
          dpdt = eos_state%dpdt
          dsdt = eos_state%dsdt
          dsdrho = eos_state%dsdr
          
          frhoT = pres_zone - p_want
          qrhoT = s_zone - s_want
          
          dtemp = (qrhoT - frhoT * dsdrho / (dpd - 0.5*delx*g_zone) ) / (dsdrho * dpdt /(dpd - 0.5*delx*g_zone) - dsdt)
          drho = -(frhoT + dpdt * dtemp)/(dpd - 0.5*delx*g_zone)
          
          dens_zone = max(0.9_dp_t*dens_zone, &
              min(dens_zone + drho, 1.1_dp_t*dens_zone))
          temp_zone = max (0.9_dp_t*temp_zone, &
              min(temp_zone + dtemp, 1.1_dp_t*temp_zone))
                          
          if (abs(drho) < TOL*dens_zone .and. abs(dtemp) < TOL*temp_zone) then
            converged_smooth = .TRUE.
            exit
          endif
          

      enddo
          
      if (.NOT. converged_smooth) then
          
          print *, 'Error zone', i, ' did not converge while creating a neutrally stratified core'
          print *, dens_zone, temp_zone
          print *, p_want
          print *, drho, dtemp
          call bl_error('Error: HSE non-convergence')
            
      endif


      ! call the EOS one more time for this zone and then go on to the next
      ! (t, rho) -> (p)
      eos_state%T     = temp_zone
      eos_state%rho   = dens_zone
      eos_state%xn(:) = xn(:)

      call eos(eos_input_rt, eos_state)

      pres_zone = eos_state%p
      
      ! update the thermodynamics in this zone
      model_hse(i,idens_model) = dens_zone
      model_hse(i,itemp_model) = temp_zone
      model_hse(i,ipres_model) = pres_zone
      
    enddo
        
  endif  




!---------------------------------------------------------------------------
! integrate further up
!---------------------------------------------------------------------------
  do i = index_base+1, brunt_base

     delx = xzn_hse(i) - xzn_hse(i-1)
     ! compute the gravitation acceleration at the lower edge
     call calc_grav_zone(i,xznr_hse,xznl_hse,xzn_hse,model_hse(:,idens_model),.True.,do_invsq_grav,g_zone)     


     ! we've already set initial guesses for density, temperature, and
     ! composition
     dens_zone = model_hse(i,idens_model)
     temp_zone = max(temp_cutoff, model_hse(i,itemp_model))
     xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)

     
     !-----------------------------------------------------------------------
     ! iteration loop
     !-----------------------------------------------------------------------

     ! start off the Newton loop by saying that the zone has not converged
     converged_hse = .FALSE.

     do iter = 1, MAX_ITER

        ! what pressure does HSE say we want?
        p_want = model_hse(i-1,ipres_model) + &
             delx*0.5*(dens_zone + model_hse(i-1,idens_model))*g_zone
         
        ! (t, rho) -> (p)
        eos_state%T   = temp_zone
        eos_state%rho = dens_zone
        eos_state%xn(:) = xn(:)

        call eos(eos_input_rt, eos_state)
              
        entropy = eos_state%s
        pres_zone = eos_state%p

        dpd = eos_state%dpdr

        drho = (p_want - pres_zone)/(dpd - 0.5*delx*g_zone)
              
        dens_zone = max(0.9*dens_zone, &
             min(dens_zone + drho, 1.1*dens_zone))
       
        if (abs(drho) < TOL*dens_zone) then
           converged_hse = .TRUE.
           exit
        endif

        if (dens_zone < low_density_cutoff) then
           dens_zone = low_density_cutoff
           temp_zone = temp_cutoff
           converged_hse = .TRUE.
           exit
        endif

     enddo

        
     if (.NOT. converged_hse) then
        print *, 'Error zone', i, ' did not converge in init_1d'
        print *, 'integrate up'
        print *, dens_zone, temp_zone
        print *, p_want
        print *, drho, dtemp
        call bl_error('Error: HSE non-convergence')
     endif
        
  
     ! call the EOS one more time for this zone and then go on to the next
     ! (t, rho) -> (p)
     eos_state%T     = temp_zone
     eos_state%rho   = dens_zone
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state)

     pres_zone = eos_state%p
     
     ! update the thermodynamics in this zone
     model_hse(i,idens_model) = dens_zone
     model_hse(i,itemp_model) = temp_zone
     model_hse(i,ipres_model) = pres_zone

     ! to make this process converge faster, set the density in the
     ! next zone to the density in this zone
     ! model_hse(i+1,idens) = dens_zone

  enddo

!-----------------------------------------------------------------------------
! Change the stratification to a constanstant brunt-vaisala frequency 
!-----------------------------------------------------------------------------

  if (change_strat) then  

    
    !set up perv_p,prev_mu and prev_temp for N^2 calculation in brunt_base
    if (use_slope) then
      dens_zone = model_hse(brunt_base-2,idens_model)
      temp_zone = max(temp_cutoff, model_hse(brunt_base-2,itemp_model))
      xn(:) = model_hse(brunt_base-2,ispec_model:ispec_model-1+nspec)
      
      eos_state%T = temp_zone
      eos_state%rho = dens_zone
      eos_state%xn(:) = xn(:)

      call eos(eos_input_rt, eos_state)
      
      prev_p = eos_state%p
      prev_mu = eos_state%abar
      prev_temp = temp_zone
      
      !set up the rest and do the calculation. This is our brunt_want
      delx = xzn_hse(brunt_base-1) - xzn_hse(brunt_base-2)
      call calc_grav_zone(brunt_base-1,xznr_hse,xznl_hse,xzn_hse,model_hse(:,idens_model),.True.,do_invsq_grav,g_zone)     
    endif

    dens_zone = model_hse(brunt_base-1,idens_model)
    temp_zone = max(temp_cutoff, model_hse(brunt_base-1,itemp_model))
    xn(:) = model_hse(brunt_base-1,ispec_model:ispec_model-1+nspec)
    
    eos_state%T = temp_zone
    eos_state%rho = dens_zone
    eos_state%xn(:) = xn(:)

    call eos(eos_input_rt, eos_state)
    if (use_slope) then
      call calc_brunt(g_zone,delx,dens_zone,temp_zone,xn,prev_p,prev_temp,prev_mu,prev_brunt)
    endif
    
    prev_p = eos_state%p
    prev_mu = eos_state%abar
    prev_temp = temp_zone
    
    !set up the rest and do the calculation. This is our brunt_want
    delx = xzn_hse(brunt_base) - xzn_hse(brunt_base-1)
    call calc_grav_zone(brunt_base,xznr_hse,xznl_hse,xzn_hse,model_hse(:,idens_model),.True.,do_invsq_grav,g_zone)     

    dens_zone = model_hse(brunt_base,idens_model)
    temp_zone = max(temp_cutoff, model_hse(brunt_base,itemp_model))
    xn(:) = model_hse(brunt_base,ispec_model:ispec_model-1+nspec)
    
    eos_state%T = temp_zone
    eos_state%rho = dens_zone
    eos_state%xn(:) = xn(:)
      
    call eos(eos_input_rt, eos_state)
    !setting the wanted N^2
    call calc_brunt(g_zone,delx,dens_zone,temp_zone,xn,prev_p,prev_temp,prev_mu,brunt_want)
    
    prev_p = eos_state%p
    prev_mu = eos_state%abar
    prev_temp = temp_zone
    
    
    print *, 'We change the stratification outwards of zone ', brunt_base
    print *, 'this is at radius ', xzn_hse(brunt_base)
    if (use_slope) then
      print *, 'and we use a slope of ', (brunt_want-prev_brunt)/delx, ' (= dN^2 / dx)'    
    else
      print *, 'and the N^2 we want to reach is ', brunt_want   
    endif
    
    !---------------------------------------------------------------------------
    ! integrate up -- using the temperature profile defined above
    ! we need to do this to use the diff_brunt function, because with a downward 
    ! integration we would need to change the gravity of the cell when computing dndrho
    !---------------------------------------------------------------------------

    
    do i = brunt_base+1, nx

      if (use_slope) then
        brunt_slope = (brunt_want-prev_brunt)/delx
        prev_brunt = brunt_want
        delx = xzn_hse(i) - xzn_hse(i-1)
        brunt_want = brunt_want + brunt_slope * delx
      endif
      
      delx = xzn_hse(i) - xzn_hse(i-1)

      ! compute the gravitation acceleration at the lower edge
      call calc_grav_zone(i,xznr_hse,xznl_hse,xzn_hse,model_hse(:,idens_model),.True.,do_invsq_grav,g_zone)     

      ! we already set the temperature profiles
      temp_zone = max(temp_cutoff, model_hse(i,itemp_model))
      ! use our previous initial guess for density
      dens_zone = model_hse(i,idens_model)
      xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)


      !-----------------------------------------------------------------------
      ! iteration loop
      !-----------------------------------------------------------------------

      ! start off the Newton loop by saying that the zone has not converged
      converged_strat = .FALSE.

      do iter = 1, MAX_ITER

          ! get the pressure we want from the HSE equation, just the
          ! zone below the current.  Note, we are using an average of
          ! the density of the two zones as an approximation of the
          ! interface value -- this means that we need to iterate for
          ! find the density and pressure that are consistent
          
          ! HSE differencing
          p_want = model_hse(i-1,ipres_model) + &
              delx*HALF*(dens_zone + model_hse(i-1,idens_model))*g_zone
        
          
          ! we need to zero:
          !   frhoT = p(rho) - p_want
          !   qrhoT = N^2 - N^2_want 
          
          ! (t, rho) -> (p)
          eos_state%T     = temp_zone
          eos_state%rho   = dens_zone
          eos_state%xn(:) = xn(:)

          call eos(eos_input_rt, eos_state)
          call calc_brunt(g_zone,delx,dens_zone,temp_zone,xn,prev_p,prev_temp,prev_mu,brunt_zone)
          
          pres_zone = eos_state%p
          
          dpd = eos_state%dpdr
          dpdt = eos_state%dpdt

          frhoT = pres_zone - p_want
          qrhoT = brunt_zone - brunt_want
          
          call diff_brunt(g_zone,delx,dens_zone,temp_zone,xn,prev_p,prev_temp,prev_mu,.True.,dndt)
          call diff_brunt(g_zone,delx,dens_zone,temp_zone,xn,prev_p,prev_temp,prev_mu,.False.,dndrho)
          
          dtemp = (qrhoT - frhoT * dndrho / (dpd - 0.5*delx*g_zone) ) / (dndrho * dpdt /(dpd - 0.5*delx*g_zone) - dndt)
          drho = -(frhoT + dpdt * dtemp)/(dpd - 0.5*delx*g_zone)

          dens_zone = max(0.9_dp_t*dens_zone, &
              min(dens_zone + drho, 1.1_dp_t*dens_zone))
          temp_zone = max (0.9_dp_t*temp_zone, &
              min(temp_zone + dtemp, 1.1_dp_t*temp_zone))
              
          if (abs(drho) < TOL*dens_zone .and. abs(dtemp) < TOL*temp_zone) then
            converged_strat = .TRUE.
            exit
          endif
          

      enddo
          
      if (.NOT. converged_strat) then
          
          print *, 'Error zone', i, ' did not converge while changing the stratification'
          print *, dens_zone, temp_zone
          print *, p_want
          print *, drho, dtemp
          call bl_error('Error: HSE non-convergence')
            
      endif


      ! call the EOS one more time for this zone and then go on to the next
      ! (t, rho) -> (p)
      eos_state%T     = temp_zone
      eos_state%rho   = dens_zone
      eos_state%xn(:) = xn(:)

      call eos(eos_input_rt, eos_state)

      pres_zone = eos_state%p
      
      !set the needed previous zone values
      prev_p = pres_zone
      prev_temp = temp_zone
      prev_mu = eos_state%abar
      
      ! update the thermodynamics in this zone
      model_hse(i,idens_model) = dens_zone
      model_hse(i,itemp_model) = temp_zone
      model_hse(i,ipres_model) = pres_zone      
    enddo 
 endif

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! now smooth the density and composition profile by a moving average between xmin_smooth and xmax_smooth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (smooth) then
    temp_model = model_hse

    !one swipe upwards with weighted moving average ... this works fine, but does depend on the resolution 
    do i = min_base+1, max_base-1
      csmooth = smoothness
      
      sumn = ONE
      do n=2,csmooth+1
        wnorm = ONE+(dble(n)-ONE)*dble(norm)/dble(nx)
        sumn= sumn + TWO/wnorm**2
      enddo
      
      sumrho = model_hse(i,idens_model)
      do n=2,csmooth+1
        wnorm = ONE+(dble(n)-ONE)*dble(norm)/dble(nx)
	sumrho = sumrho + model_hse(i-n+1,idens_model)/wnorm**2 + model_hse(i+n-1,idens_model)/wnorm**2
      enddo
      
      do j = 1, nspec-1
	sumxn(j) = model_hse(i,ispec_model-1+j)
	do n=2,csmooth+1
          wnorm = ONE+(dble(n)-ONE)*dble(norm)/dble(nx)
	  sumxn(j) = sumxn(j) + model_hse(i-n+1,ispec_model-1+j)/wnorm**2 + model_hse(i+n-1,ispec_model-1+j)/wnorm**2
	enddo
      enddo
      
      temp_model(i,idens_model) = sumrho / sumn
      do j = 1, nspec-1
	temp_model(i,ispec_model-1+j) = sumxn(j) / sumn
      enddo    

    enddo


    
    model_hse = temp_model
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Normalize the composition again!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    
    do i=1,nx
      sumx = ZERO
      do j=1,nspec
        sumx = sumx + model_hse(i,ispec_model-1+j)
      enddo
      model_hse(i,ispec_model:ispec_model+nspec-1) = model_hse(i,ispec_model:ispec_model+nspec-1) / sumx
    enddo
  
  endif


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! write an intermediate output to show how the smoothing / stratification change worked
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  !compute the brunt-vaisala frequency at each zone
  brunt(1) = 0
  
  dens_zone = model_hse(1,idens_model)
  temp_zone = max(temp_cutoff, model_hse(1,itemp_model))
  xn(:) = model_hse(1,ispec_model:ispec_model-1+nspec)
  
  eos_state%T = temp_zone
  eos_state%rho = dens_zone
  eos_state%xn(:) = xn(:)
    
  call eos(eos_input_rt, eos_state)
  
  prev_p = eos_state%p
  prev_mu = eos_state%abar
  prev_temp = temp_zone

  do i = 2,nx 
     delx = xzn_hse(i) - xzn_hse(i-1)
     ! compute the gravitation acceleration at the lower edge
     call calc_grav_zone(i,xznr_hse,xznl_hse,xzn_hse,model_hse(:,idens_model),.True.,do_invsq_grav,g_zone)     


    dens_zone = model_hse(i,idens_model)
    temp_zone = max(temp_cutoff, model_hse(i,itemp_model))
    xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)

    eos_state%T = temp_zone
    eos_state%rho = dens_zone
    eos_state%xn(:) = xn(:)
    
    call eos(eos_input_rt, eos_state)
    call calc_brunt(g_zone,delx,dens_zone,temp_zone,xn,prev_p,prev_temp,prev_mu,brunt(i))
    
    prev_p = eos_state%p
    prev_mu = eos_state%abar
    prev_temp = temp_zone
  enddo

  !print the newly smoothed profiles into a new file 
  
  outfile = trim(model_prefix) // ".smoothed_profile"
  open (newunit=lun1, file=outfile, status="unknown")


  write (lun1,2002), 'Smoothed profile: '
  write (lun1,2001), 'radius','entropy','pressure','density','X(H1)','N^2'
  do i = 1, nx
      dens_zone = model_hse(i,idens_model)
      temp_zone = max(temp_cutoff, model_hse(i,itemp_model))
      xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)
      
      eos_state%T = temp_zone
      eos_state%rho = dens_zone
      eos_state%xn(:) = xn(:)
      
      call eos(eos_input_rt, eos_state)
      
      write (lun1,2000), xzn_hse(i), eos_state%s, eos_state%p, eos_state%rho, eos_state%xn(1),brunt(i)
  enddo


  close (unit=lun1)
  
  
  
  
  
!-----------------------------------------------------------------------------
! Do a final sweep, just ensuring HSE. This is required, if the composition profile was smoothed 
!-----------------------------------------------------------------------------
if (put_in_hse .or. smooth) then
  ! the HSE state will be done respecting the interpolated temperature 
  ! from the initial model.  When the temperature drops below T_lo,
  ! we floor it.
  
  ! make it all thermodynamically consistent in the centre
  eos_state%rho = model_hse(1,idens_model)
  eos_state%T = model_hse(1,itemp_model)
  eos_state%xn(:) = model_hse(1,ispec_model:ispec_model-1+nspec)

  call eos(eos_input_rt, eos_state)

  model_hse(1,ipres_model) = eos_state%p
  
  
  !---------------------------------------------------------------------------
  ! integrate up
  !---------------------------------------------------------------------------
  do i = 2, nx

     delx = xzn_hse(i) - xzn_hse(i-1)
     ! compute the gravitation acceleration at the lower edge
     call calc_grav_zone(i,xznr_hse,xznl_hse,xzn_hse,model_hse(:,idens_model),.True.,do_invsq_grav,g_zone)     


     ! we've already set initial guesses for density, temperature, and
     ! composition
     dens_zone = model_hse(i,idens_model)
     temp_zone = max(temp_cutoff, model_hse(i,itemp_model))
     xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)

     
     !-----------------------------------------------------------------------
     ! iteration loop
     !-----------------------------------------------------------------------

     ! start off the Newton loop by saying that the zone has not converged
     converged_hse = .FALSE.

     do iter = 1, MAX_ITER

        ! what pressure does HSE say we want?
        p_want = model_hse(i-1,ipres_model) + &
             delx*0.5*(dens_zone + model_hse(i-1,idens_model))*g_zone
         
        ! (t, rho) -> (p)
        eos_state%T   = temp_zone
        eos_state%rho = dens_zone
        eos_state%xn(:) = xn(:)

        call eos(eos_input_rt, eos_state)
              
        entropy = eos_state%s
        pres_zone = eos_state%p

        dpd = eos_state%dpdr

        drho = (p_want - pres_zone)/(dpd - 0.5*delx*g_zone)
              
        dens_zone = max(0.9*dens_zone, &
             min(dens_zone + drho, 1.1*dens_zone))
       
        if (abs(drho) < TOL*dens_zone) then
           converged_hse = .TRUE.
           exit
        endif

        if (dens_zone < low_density_cutoff) then
           dens_zone = low_density_cutoff
           temp_zone = temp_cutoff
           converged_hse = .TRUE.
           exit
        endif

     enddo

        
     if (.NOT. converged_hse) then
        print *, 'Error zone', i, ' did not converge in init_1d'
        print *, 'integrate up'
        print *, dens_zone, temp_zone
        print *, p_want
        print *, drho, dtemp
        call bl_error('Error: HSE non-convergence')
     endif
        
  
     ! call the EOS one more time for this zone and then go on to the next
     ! (t, rho) -> (p)
     eos_state%T     = temp_zone
     eos_state%rho   = dens_zone
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state)

     pres_zone = eos_state%p
     
     ! update the thermodynamics in this zone
     model_hse(i,idens_model) = dens_zone
     model_hse(i,itemp_model) = temp_zone
     model_hse(i,ipres_model) = pres_zone

     ! to make this process converge faster, set the density in the
     ! next zone to the density in this zone
     ! model_hse(i+1,idens) = dens_zone

  enddo

endif !put_in_hse  




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! final output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(num,'(i8)') nx

  dxstr = num_to_unitstring(dCoord)

  outfile = trim(model_prefix) // ".hse" // ".dx_" // trim(adjustl(dxstr))
  outfile2 = trim(outfile) // ".extras"
  outfile3 = trim(model_prefix) // '.binary.bin'

  open (newunit=lun1, file=outfile, status="unknown")
  open (newunit=lun2, file=outfile2, status="unknown")  
  
  open (newunit=lun3, file=outfile3, status="unknown", form="unformatted",access="stream")
  write (lun3) nx
  write (lun3) nvars_model

  write (lun3) xzn_hse(:)
  write (lun3) model_hse(:,:)
  
  write (lun3) var_names_size
  allocate(varnames_stored(nspec+3))
  varnames_stored(1) = "density"
  varnames_stored(2) = "temperature"
  varnames_stored(3) = "pressure"
  
  do n = 4, nspec+3
   varnames_stored(n) = spec_names(n-3)
  enddo
  write (lun3) varnames_stored(:)
  
  close (unit = lun3)  
  
  
  write (lun1,1001) "# npts = ", nx
  write (lun1,1001) "# num of variables = ", nvars_model
  write (lun1,1002) "# density"
  write (lun1,1002) "# temperature"
  write (lun1,1002) "# pressure"

  do n = 1, nspec
     write (lun1, 1003) "# ", spec_names(n)
  enddo

1000 format (1x, 100(g26.16, 1x))
1001 format (a, i5)
1002 format (a)
1003 format (a,a)

  do i = 1, nx
     write (lun1,1000) xzn_hse(i), model_hse(i,idens_model), &
                       model_hse(i,itemp_model), model_hse(i,ipres_model), &
                      (model_hse(i,ispec_model-1+n), n=1,nspec)
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! .extra output
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write (lun2,1001), "# npts = ", nx
  write (lun2,1001), "# num of variables = ", 3
  write (lun2,1002), "# entropy"
  write (lun2,1002), "# c_s"
  write (lun2,1002), "# N^2"
  
    !compute the brunt-vaisala frequency at each zone
  brunt(1) = 0
  
  dens_zone = model_hse(1,idens_model)
  temp_zone = max(temp_cutoff, model_hse(1,itemp_model))
  xn(:) = model_hse(1,ispec_model:ispec_model-1+nspec)
  
  eos_state%T = temp_zone
  eos_state%rho = dens_zone
  eos_state%xn(:) = xn(:)
    
  call eos(eos_input_rt, eos_state)
  
  prev_p = eos_state%p
  prev_mu = eos_state%abar
  prev_temp = temp_zone
  
  do i = 2,nx 
     delx = xzn_hse(i) - xzn_hse(i-1)
     ! compute the gravitation acceleration at the lower edge
     call calc_grav_zone(i,xznr_hse,xznl_hse,xzn_hse,model_hse(:,idens_model),.True.,do_invsq_grav,g_zone)     


    dens_zone = model_hse(i,idens_model)
    temp_zone = max(temp_cutoff, model_hse(i,itemp_model))
    xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)
    
    eos_state%T = temp_zone
    eos_state%rho = dens_zone
    eos_state%xn(:) = xn(:)
    
    call eos(eos_input_rt, eos_state)
    call conducteos(eos_input_rt, eos_state, .false., conductivity(i))  
    
    
    call calc_brunt(g_zone,delx,dens_zone,temp_zone,xn,prev_p,prev_temp,prev_mu,brunt(i))
    
    prev_p = eos_state%p
    prev_mu = eos_state%abar
    prev_temp = temp_zone
  enddo
  
  ! test: bulk EOS call -- Maestro will do this once we are mapped, so make
  ! sure that we are in HSE with updated thermodynamics
  do i = 1, nx
     eos_state%rho = model_hse(i,idens_model)
     eos_state%T = model_hse(i,itemp_model)
     eos_state%xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)

     call eos(eos_input_rt, eos_state)

     model_hse(i,ipres_model) = eos_state%p

     write (lun2,1000), xzn_hse(i), eos_state%s, eos_state%cs, brunt(i), conductivity(i)
  enddo
  
  
  ! compute the maximum HSE error
  max_hse_error = -1.d30

  do i = 2, nx-1
     delx = xzn_hse(i) - xzn_hse(i-1)
     ! compute the gravitation acceleration at the lower edge
     call calc_grav_zone(i,xznr_hse,xznl_hse,xzn_hse,model_hse(:,idens_model),.True.,do_invsq_grav,g_zone)     

     dpdr = (model_hse(i,ipres_model) - model_hse(i-1,ipres_model))/dCoord
     rhog = HALF*(model_hse(i,idens_model) + model_hse(i-1,idens_model))*g_zone
     if (dpdr /= ZERO .and. model_hse(i+1,idens_model) > low_density_cutoff) then
        max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(dpdr))
        write(lun2,'(g26.16)') abs(dpdr - rhog)/abs(dpdr)
     endif

  enddo

  print *, 'maximum HSE error = ', max_hse_error
  print *, ' '
  
  
  close (unit=lun1)
  close (unit=lun2)

 contains
 
 subroutine calc_brunt(g_zone,delx,dens_zone,temp_zone,xn,prev_p,prev_temp, prev_mu,brunt)
    use eos_module, only: eos_input_rt, eos, eos_init
    use eos_type_module, only: eos_t
    use bl_types
  
    real(kind=dp_t), intent(in) :: delx,g_zone,dens_zone,temp_zone,xn(:)
    real(kind=dp_t), intent(in) :: prev_p,prev_temp,prev_mu
    real(kind=dp_t), intent(out) :: brunt

    !local variables
    real(kind=dp_t) :: pres_zone,Hp
    real(kind=dp_t) :: dpd,dpdt,dpda,gam,chi_rho,chi_temp,chi_mu,adiab,dtdr
    real(kind=dp_t) :: grad_temp, grad_mu, grad_ad
    
    type (eos_t) :: eos_state
    
    
    !compute N^2

    eos_state%T = temp_zone
    eos_state%rho = dens_zone
    eos_state%xn(:) = xn(:)
    
    call eos(eos_input_rt, eos_state)    
    
    dpd = eos_state%dpdr
    dpdt = eos_state%dpdt
    dpda = eos_state%dpda
    gam = eos_state%gam1
    pres_zone = eos_state%p
    chi_rho = dens_zone * dpd / pres_zone
    chi_temp = temp_zone * dpdt / pres_zone
    chi_mu = eos_state%abar * dpda/ pres_zone
    Hp = -delx / (dlog(pres_zone)-dlog(prev_p))
    grad_temp = (dlog(temp_zone) - dlog(prev_temp)) / (dlog(pres_zone) - dlog(prev_p))
    grad_mu = (dlog(eos_state%abar) - dlog(prev_mu)) / (dlog(pres_zone) - dlog(prev_p))
    grad_ad = (gam - chi_rho)/(chi_temp*gam)
    
    dtdr = (temp_zone-prev_temp) / delx 
    adiab = (temp_zone-prev_temp) / delx - grad_temp
    brunt = - g_zone * chi_temp / (chi_rho * Hp) * (grad_ad - grad_temp - chi_mu / chi_temp * grad_mu) 
    
end subroutine calc_brunt


subroutine diff_brunt(g_zone,delx,dens,temp,xn,prev_p,prev_temp,prev_mu,temp_diff,brunt_diff)
    use bl_constants_module
    use bl_types  
    
    real(kind=dp_t), intent(in)  :: delx,g_zone,dens,temp,xn(:),prev_p,prev_temp,prev_mu
    real(kind=dp_t), intent(out) :: brunt_diff
    logical :: temp_diff
    
    !local variables
    real(kind=dp_t) :: dens_zone,temp_zone
    real(kind=dp_t) :: brunt_l,brunt_h,h
    
    !we compute the partial derivative of the brunt-vaisala frequency with temperature or density
    !we use a central differencing
    
    
    !first compute N(rho+h,T,X) or N(rho,T+h,X)
    if (temp_diff) then
      temp_zone = temp + temp/1000.0
      dens_zone = dens
    else
      temp_zone = temp
      dens_zone = dens + dens/1000.0
    endif 
    
    call calc_brunt(g_zone,delx,dens_zone,temp_zone,xn,prev_p,prev_temp,prev_mu,brunt_h)
    
    !second compute N(rho-h,T,X) or N(rho,T-h,X)
    if (temp_diff) then
      temp_zone = temp - temp/1000.0
      dens_zone = dens
    else
      temp_zone = temp
      dens_zone = dens - dens/1000.0
    endif 
    
    call calc_brunt(g_zone,delx,dens_zone,temp_zone,xn,prev_p,prev_temp,prev_mu,brunt_l)
   
   !third compute the differential quotient
    if (temp_diff) then
      h = temp/1000.0
    else
      h = dens/1000.0
    endif 
    
   brunt_diff = (brunt_h - brunt_l) / (TWO*h)
   
end subroutine diff_brunt

subroutine calc_grav_zone(i,r_r,r_l,r_m,dens,inner_edge,do_invsq_grav,g_zone)
  
  use fundamental_constants_module, only: Gconst
  use bl_constants_module
  use bl_types
  
  real(kind=dp_t), intent(in)  :: r_r(:),r_l(:),r_m(:),dens(:) 
  real(kind=dp_t), intent(out) :: g_zone
  integer, intent(in) :: i
  logical, intent(in) :: inner_edge, do_invsq_grav
  
  !local variables
  real(kind=dp_t) :: M_shell,M_enclosed
  integer :: j
  
  if (do_invsq_grav) then
    if (inner_edge) then
      ! compute the gravitation acceleration at the lower edge
      M_enclosed = four3rd*m_pi *  r_l(1)**3 * dens(1)  
        do j = 1, i-1
            M_shell = dens(j) * & 
                (four3rd*m_pi * (r_r(j)-r_l(j)) * (r_l(j)**2 + r_r(j)**2 + r_l(j)*r_r(j)) )
            M_enclosed = M_enclosed + M_shell
        enddo
        g_zone = -Gconst*M_enclosed/r_l(i)**2
    else 
      ! compute the gravitation acceleration at the upper edge
      M_enclosed = four3rd*m_pi *  r_l(1)**3 * dens(1)  
        do j = 1, i
            M_shell = dens(j) * & 
                (four3rd*m_pi * (r_r(j)-r_l(j)) * (r_l(j)**2 + r_r(j)**2 + r_l(j)*r_r(j)) )
            M_enclosed = M_enclosed + M_shell
        enddo
        g_zone = -Gconst*M_enclosed/r_r(i)**2
    endif
  
  else
    g_zone = g_const
  endif  
  
end subroutine calc_grav_zone
  
end program init_1d


function num_to_unitstring(value)

  use bl_types
  implicit none

  real (kind=dp_t) :: value
  character (len=32) :: num_to_unitstring
  character (len=16) :: temp

  if (value > 1.d5) then

     ! work in km
     write(temp,'(f6.3)') value/1.d5     
     num_to_unitstring = trim(temp) // "km"
  else

     ! work in cm
     if (value > 1.d3) then
        write(temp,'(f8.3)') value
        num_to_unitstring = trim(temp) // "cm"

     else
        write(temp,'(f6.3)') value
        num_to_unitstring = trim(temp) // "cm"
     endif

  endif

  return 
end function num_to_unitstring


