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
  real (kind=dp_t), allocatable :: brunt(:), s(:)

  real (kind=dp_t) :: sumx, sumrho, sumn
  real (kind=dp_t), DIMENSION(nspec) :: sumxn
  real (kind=dp_t) :: coreX, frhoT, qrhoT 
  integer :: comp
 
  real (kind=dp_t) :: dens_zone, temp_zone, pres_zone, entropy, s_zone
  real (kind=dp_t) :: dpd, dpdt, dsdt, dsdrho, dtdr, gam, Hp, dpda, adiab
  real (kind=dp_t) :: prev_mu, prev_p, prev_temp, prev_dtdr, prev_adiab
  real (kind=dp_t) :: grad_temp, grad_ad, grad_mu
  real (kind=dp_t) :: chi_rho, chi_temp, chi_mu    
  
  
  
  
  real (kind=dp_t) :: p_want, drho, dtemp, delx, s_want
  
  real (kind=dp_t) :: g_zone, g_const, M_enclosed, M_shell
  logical :: do_invsq_grav

  real (kind=dp_t), parameter :: TOL = 1.e-12

  integer, parameter :: MAX_ITER = 250

  integer :: iter

  logical :: converged_hse, fluff, smooth, converged_smooth

  real (kind=dp_t) :: low_density_cutoff, temp_cutoff, smallx, max_T
  real (kind=dp_t) :: model_shift

  integer :: index_base, conv_base, min_base, max_base, cent_base, smoothness, csmooth, norm
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
                    temp_cutoff, do_invsq_grav, &
                    low_density_cutoff, model_prefix, model_shift, smooth, &
                    xmin_smooth, xmax_smooth, smoothness, norm

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

  ! this comes in via extern_probin_module -- override the default
  ! here if we want
  use_eos_coulomb = .true.


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
  allocate(s(nx))


  ! compute the coordinates of the new gridded function
  dCoord = (xmax - xmin) / dble(nx)

  do i = 1, nx
     xznl_hse(i) = xmin + (dble(i) - ONE)*dCoord
     xzn_hse(i)  = xmin + (dble(i) - HALF)*dCoord
     xznr_hse(i) = xmin + (dble(i))*dCoord
  enddo
  

!-----------------------------------------------------------------------------
! put the model onto our new uniform grid
!-----------------------------------------------------------------------------

  fluff = .false.

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


  ! find the index to integrate from by looking for the peak temperature
  index_base = -1
  max_T = -1.0_dp_t

  do i = 1, nx
     if (model_hse(i,itemp_model) > max_T) then
        index_base = i
        max_T = model_hse(i,itemp_model)
     endif
  enddo
     
  if (index_base == -1) then
     call bl_error('ERROR: invalid base_height')
  endif
  
  print *, 'index_base = ', index_base

  ! make the base thermodynamics consistent for this base point -- that is
  ! what we will integrate from!
  eos_state%rho = model_hse(index_base,idens_model)
  eos_state%T = model_hse(index_base,itemp_model)
  eos_state%xn(:) = model_hse(index_base,ispec_model:ispec_model-1+nspec)

  call eos(eos_input_rt, eos_state)

  model_hse(index_base,ipres_model) = eos_state%p
  

   
  
  !make an output of the initial (interpolated) profiles
  
  
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
  prev_dtdr = huge(0.0)
  prev_adiab = huge(0.0)
  
  do i = 2,nx 
     delx = xzn_hse(i) - xzn_hse(i-1)
     ! compute the gravitation acceleration at the lower edge
     M_enclosed = four3rd*m_pi *  xznl_hse(1)**3 * model_hse(1,idens_model)     
     if (do_invsq_grav) then
        do j = 1, i-1
           M_shell = model_hse(j,idens_model) * & 
               (four3rd*m_pi * (xznr_hse(j)-xznl_hse(j)) * (xznl_hse(j)**2 + xznr_hse(j)**2 + xznl_hse(j)*xznr_hse(j)) )
           M_enclosed = M_enclosed + M_shell
        enddo
        g_zone = -Gconst*M_enclosed/xznl_hse(i)**2
     else
        g_zone = g_const
     endif

    dens_zone = model_hse(i,idens_model)
    temp_zone = max(temp_cutoff, model_hse(i,itemp_model))
    xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)
    
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
    s(i) = - prev_dtdr / dtdr * (adiab/prev_adiab)
    brunt(i) = - g_zone * chi_temp / (chi_rho * Hp) * (grad_ad - grad_temp - chi_mu / chi_temp * grad_mu) 
!     dtdr = (temp_zone - max(temp_cutoff, model_hse(i-1,itemp_model)))/delx
!     brunt(i) = - g_zone/temp_zone * (abs(dtdr) - g_zone/eos_state%cp)
    
    prev_adiab = adiab
    prev_dtdr = dtdr
    prev_p = eos_state%p
    prev_mu = eos_state%abar
    prev_temp = temp_zone

  enddo
  
  
  !!write everything into a file
  outfile = trim(model_prefix) // ".initial_interpolated_profile"
  open (newunit=lun1, file=outfile, status="unknown")


  write (lun1,2002), 'Initial interpolated profile: '
  write (lun1,2001), 'radius','entropy','pressure','density','X(H1)', 'N^2', 'Stiffness'
  do i = 1, nx
    dens_zone = model_hse(i,idens_model)
    temp_zone = max(temp_cutoff, model_hse(i,itemp_model))
    xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)
    
    eos_state%T = temp_zone
    eos_state%rho = dens_zone
    eos_state%xn(:) = xn(:)
    
    call eos(eos_input_rt, eos_state)
    
    write (lun1,2000), xzn_hse(i), eos_state%s, eos_state%p, eos_state%rho, eos_state%xn(1), brunt(i), s(i)
  enddo


2000 format (1x, 100(g26.16, 1x))
2001 format (1x,100(a26, 1x))
2002 format (a)

  close (unit=lun1)

  
!-----------------------------------------------------------------------------
!Create a neutrally stratified core and smooth the density and composition profile afterwards
!-----------------------------------------------------------------------------

  if (smooth) then  
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
    
    !keep the entropy at the core constant and enforce hse by integrating inwards from the convective boundary
    xn(:) = model_hse(conv_base,ispec_model:ispec_model-1+nspec)
    
    eos_state%T     = model_hse(conv_base+1,itemp_model)
    eos_state%rho   = model_hse(conv_base+1,idens_model)
    eos_state%xn(:) = xn(:)

    call eos(eos_input_rt, eos_state)
    
    s_want = eos_state%s
    !---------------------------------------------------------------------------
    ! integrate down -- using the temperature profile defined above
    !---------------------------------------------------------------------------
    do i = conv_base, 1, -1

      delx = xzn_hse(i+1) - xzn_hse(i)

      ! compute the gravitation acceleration at the upper edge
      M_enclosed = four3rd*m_pi *  xznl_hse(1)**3 * model_hse(1,idens_model)
      if (do_invsq_grav) then
	  do j = 1, i
	    M_shell = model_hse(j,idens_model) * & 
		(four3rd*m_pi * (xznr_hse(j)-xznl_hse(j)) * (xznl_hse(j)**2 + xznr_hse(j)**2 + xznl_hse(j)*xznr_hse(j)) )
	    M_enclosed = M_enclosed + M_shell
	  enddo
	  
	  g_zone = -Gconst*M_enclosed/xznr_hse(i)**2
      else
	  g_zone = g_const
      endif

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
	  p_want = model_hse(i+1,ipres_model) - &
	      delx*0.5*(dens_zone + model_hse(i+1,idens_model))*g_zone

	  
	  ! we need to zero:
	  !   frhoT = p_want - p(rho)
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
  
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
! now smooth the density and composition profile by a moving average between xmin_smooth and xmax_smooth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    temp_model = model_hse
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
    print *, ''
    
    
    
    
!     ! two swipes around the center of the smoothing region
!     !one swipe upwards to the center of the smoothing zone
!     do i = min_base+smoothness, cent_base
! !       if ((i-min_base)*2 - 1 .lt. smoothness) then
! !         csmooth = (i-min_base)*2 -1
! !       else
!         csmooth = smoothness
! !       endif
! 
!       sumrho = sum(model_hse(i-csmooth:i+csmooth,idens_model))
!       do j = 1, nspec-1
! 	sumxn(j) = sum(model_hse(i-csmooth:i+csmooth,ispec_model-1+j))
!       enddo
!       
!       temp_model(i,idens_model) = sumrho / DBLE(csmooth*2+1)
!       do j = 1, nspec-1
! 	temp_model(i,ispec_model-1+j) = sumxn(j) / DBLE(csmooth*2+1)
!       enddo    
! 
!     enddo
!     
!     model_hse = temp_model
!     
!     
!     !one swipe downwards to the center of the smoothing zone
!     do i = max_base-smoothness, cent_base, -1
! !       if ((max_base-i )*2 - 1 .lt. smoothness) then
! ! 	csmooth = (max_base-i )*2 -1
! !       else
!         csmooth = smoothness
! !       endif
!       
!       sumrho = sum(model_hse(i-csmooth:i+csmooth,idens_model))
!       do j = 1, nspec-1
! 	sumxn(j) = sum(model_hse(i-csmooth:i+csmooth,ispec_model-1+j))
!       enddo
!       
!       
!       
!       temp_model(i,idens_model) = sumrho / DBLE(csmooth*2+1)
!       do j = 1, nspec-1
! 	temp_model(i,ispec_model-1+j) = sumxn(j) / DBLE(csmooth*2+1)
!       enddo   
!     enddo  
!     
!     model_hse = temp_model
!      
!      ! one swipe upwards -> make effectively two swipes
!      do i = min_base+1, max_base-1
! !         if (((i-min_base)*2 - 1 .lt. smoothness) .or. ((max_base-i )*2 - 1 .lt. smoothness) ) then
! !           csmooth = min(i-min_base,max_base-i)*2 -1
! !         else
!          csmooth = smoothness
! !         endif
!  
!        sumrho = sum(model_hse(i-csmooth:i+csmooth,idens_model))
!        do j = 1, nspec-1
!  	sumxn(j) = sum(model_hse(i-csmooth:i+csmooth,ispec_model-1+j))
!        enddo
!        
!        temp_model(i,idens_model) = sumrho / DBLE(csmooth*2+1)
!        do j = 1, nspec-1
!  	temp_model(i,ispec_model-1+j) = sumxn(j) / DBLE(csmooth*2+1)
!        enddo    
!  
!      enddo
    
!     !one swipe upwards with weighted moving average ... this works fine, but does depend on the resolution 
!     do i = min_base+1, max_base-1
!       csmooth = smoothness
!       
!       sumn = ONE
!       do n=2,csmooth+1
!         sumn= sumn + TWO/dble(n)**2
!       enddo
!       
!       sumrho = model_hse(i,idens_model)
!       do n=2,csmooth+1
!         sumrho = sumrho + model_hse(i-n+1,idens_model)/dble(n)**2 + model_hse(i+n-1,idens_model)/dble(n)**2
!       enddo
!       
!       do j = 1, nspec-1
! 	sumxn(j) = model_hse(i,ispec_model-1+j)
! 	do n=2,csmooth+1
! 	  sumxn(j) = sumxn(j) + model_hse(i-n+1,ispec_model-1+j)/dble(n)**2 + model_hse(i+n-1,ispec_model-1+j)/dble(n)**2
! 	enddo
!       enddo
!       
!       temp_model(i,idens_model) = sumrho / sumn
!       do j = 1, nspec-1
! 	temp_model(i,ispec_model-1+j) = sumxn(j) / sumn
!       enddo    
! 
!     enddo



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

  endif
  
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
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
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
  prev_dtdr = huge(0.0)
  prev_adiab = huge(0.0)
  
  do i = 2,nx 
     delx = xzn_hse(i) - xzn_hse(i-1)
     ! compute the gravitation acceleration at the lower edge
     M_enclosed = four3rd*m_pi *  xznl_hse(1)**3 * model_hse(1,idens_model)     
     if (do_invsq_grav) then
        do j = 1, i-1
           M_shell = model_hse(j,idens_model) * & 
               (four3rd*m_pi * (xznr_hse(j)-xznl_hse(j)) * (xznl_hse(j)**2 + xznr_hse(j)**2 + xznl_hse(j)*xznr_hse(j)) )
           M_enclosed = M_enclosed + M_shell
        enddo
        g_zone = -Gconst*M_enclosed/xznl_hse(i)**2
     else
        g_zone = g_const
     endif

    dens_zone = model_hse(i,idens_model)
    temp_zone = max(temp_cutoff, model_hse(i,itemp_model))
    xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)
    
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
    s(i) = - prev_dtdr / dtdr * (adiab/prev_adiab)
    brunt(i) = - g_zone * chi_temp / (chi_rho * Hp) * (grad_ad - grad_temp - chi_mu / chi_temp * grad_mu) 
!     dtdr = (temp_zone - max(temp_cutoff, model_hse(i-1,itemp_model)))/delx
!     brunt(i) = - g_zone/temp_zone * (abs(dtdr) - g_zone/eos_state%cp)
    
    prev_adiab = adiab
    prev_dtdr = dtdr
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
! HSE + entropy solve
!-----------------------------------------------------------------------------

  ! the HSE state will be done respecting the interpolated temperature 
  ! from the initial model.  When the temperature drops below T_lo,
  ! we floor it.

  !---------------------------------------------------------------------------
  ! integrate up
  !---------------------------------------------------------------------------
  do i = index_base+1, nx

     delx = xzn_hse(i) - xzn_hse(i-1)

     ! compute the gravitation acceleration at the lower edge
     M_enclosed = four3rd*m_pi *  xznl_hse(1)**3 * model_hse(1,idens_model)     
     if (do_invsq_grav) then
        do j = 1, i-1
           M_shell = model_hse(j,idens_model) * & 
               (four3rd*m_pi * (xznr_hse(j)-xznl_hse(j)) * (xznl_hse(j)**2 + xznr_hse(j)**2 + xznl_hse(j)*xznr_hse(j)) )
           M_enclosed = M_enclosed + M_shell
        enddo
        g_zone = -Gconst*M_enclosed/xznl_hse(i)**2
     else
        g_zone = g_const
     endif

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


  !---------------------------------------------------------------------------
  ! integrate down -- using the temperature profile defined above
  !---------------------------------------------------------------------------
  do i = index_base-1, 1, -1

     delx = xzn_hse(i+1) - xzn_hse(i)

     ! compute the gravitation acceleration at the upper edge
     M_enclosed = four3rd*m_pi *  xznl_hse(1)**3 * model_hse(1,idens_model)
     if (do_invsq_grav) then
        do j = 1, i
           M_shell = model_hse(j,idens_model) * & 
               (four3rd*m_pi * (xznr_hse(j)-xznl_hse(j)) * (xznl_hse(j)**2 + xznr_hse(j)**2 + xznl_hse(j)*xznr_hse(j)) )
           M_enclosed = M_enclosed + M_shell
        enddo
        
        g_zone = -Gconst*M_enclosed/xznr_hse(i)**2
     else
        g_zone = g_const
     endif

     ! we already set the temperature and composition profiles
     temp_zone = max(temp_cutoff, model_hse(i,itemp_model))
     xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)

     ! use our previous initial guess for density
     dens_zone = model_hse(i,idens_model)


     !-----------------------------------------------------------------------
     ! iteration loop
     !-----------------------------------------------------------------------

     ! start off the Newton loop by saying that the zone has not converged
     converged_hse = .FALSE.

     do iter = 1, MAX_ITER

        ! get the pressure we want from the HSE equation, just the
        ! zone below the current.  Note, we are using an average of
        ! the density of the two zones as an approximation of the
        ! interface value -- this means that we need to iterate for
        ! find the density and pressure that are consistent
        
        ! HSE differencing
        p_want = model_hse(i+1,ipres_model) - &
             delx*0.5*(dens_zone + model_hse(i+1,idens_model))*g_zone

        
        ! we will take the temperature already defined in model_hse
        ! so we only need to zero:
        !   A = p_want - p(rho)
        
        ! (t, rho) -> (p)
        eos_state%T     = temp_zone
        eos_state%rho   = dens_zone
        eos_state%xn(:) = xn(:)

        call eos(eos_input_rt, eos_state)
        
        pres_zone = eos_state%p
        
        dpd = eos_state%dpdr
              
        A = p_want - pres_zone
              
        drho = A/(dpd + 0.5*delx*g_zone)
        
        dens_zone = max(0.9_dp_t*dens_zone, &
             min(dens_zone + drho, 1.1_dp_t*dens_zone))

                        
        if (abs(drho) < TOL*dens_zone) then
           converged_hse = .TRUE.
           exit
        endif
        

     enddo
        
     if (.NOT. converged_hse) then
        
        print *, 'Error zone', i, ' did not converge in init_1d'
        print *, 'integrate down'
        print *, dens_zone, temp_zone
        print *, p_want
        print *, drho
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

1000 format (1x, 100(g36.26, 1x))
1001 format (a, i5)
1002 format (a)
1003 format (a,a)

  do i = 1, nx
     write (lun1,1000) xzn_hse(i), model_hse(i,idens_model), &
                       model_hse(i,itemp_model), model_hse(i,ipres_model), &
                      (model_hse(i,ispec_model-1+n), n=1,nspec)
  enddo


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
  prev_dtdr = huge(0.0)
  prev_adiab = huge(0.0)
  
  do i = 2,nx 
     delx = xzn_hse(i) - xzn_hse(i-1)
     ! compute the gravitation acceleration at the lower edge
     M_enclosed = four3rd*m_pi *  xznl_hse(1)**3 * model_hse(1,idens_model)  
     if (do_invsq_grav) then
        do j = 1, i-1
           M_shell = model_hse(j,idens_model) * & 
               (four3rd*m_pi * (xznr_hse(j)-xznl_hse(j)) * (xznl_hse(j)**2 + xznr_hse(j)**2 + xznl_hse(j)*xznr_hse(j)) )
           M_enclosed = M_enclosed + M_shell
        enddo
        g_zone = -Gconst*M_enclosed/xznl_hse(i)**2
     else
        g_zone = g_const
     endif

    dens_zone = model_hse(i,idens_model)
    temp_zone = max(temp_cutoff, model_hse(i,itemp_model))
    xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)
    
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
    s(i) = - prev_dtdr / dtdr * (adiab/prev_adiab)
    brunt(i) = - g_zone * chi_temp / (chi_rho * Hp) * (grad_ad - grad_temp - chi_mu / chi_temp * grad_mu) 
!     dtdr = (temp_zone - max(temp_cutoff, model_hse(i-1,itemp_model)))/delx
!     brunt(i) = - g_zone/temp_zone * (abs(dtdr) - g_zone/eos_state%cp)
    
    prev_adiab = adiab
    prev_dtdr = dtdr
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

     write (lun2,1000), xzn_hse(i), eos_state%s, eos_state%cs, brunt(i)
  enddo
  
  
  ! compute the maximum HSE error
  max_hse_error = -1.d30

  do i = 2, nx-1

     ! compute the gravitation acceleration at the lower edge
     M_enclosed = four3rd*m_pi *  dCoord**3 * model_hse(1,idens_model)     
     if (do_invsq_grav) then
        do j = 2, i-1
           M_shell = model_hse(j,idens_model) * & 
               (four3rd*m_pi * (dCoord) * (xznl_hse(j)**2 + xznr_hse(j)**2 + xznl_hse(j)*xznr_hse(j)) )
           M_enclosed = M_enclosed + M_shell
        enddo
        g_zone = -Gconst*M_enclosed/xznl_hse(i)**2
     else
        g_zone = g_const
     endif
     

     dpdr = (model_hse(i,ipres_model) - model_hse(i-1,ipres_model))/dCoord
     rhog = HALF*(model_hse(i,idens_model) + model_hse(i-1,idens_model))*g_zone
  !   print *, abs(dpdr - rhog)/abs(dpdr)
     if (dpdr /= ZERO .and. model_hse(i+1,idens_model) > low_density_cutoff) then
        max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(dpdr))
        write(lun2,'(g26.16)') abs(dpdr - rhog)/abs(dpdr)
     endif

  enddo

  print *, 'maximum HSE error = ', max_hse_error
  print *, ' '
  
  
  close (unit=lun1)
  close (unit=lun2)
  
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


