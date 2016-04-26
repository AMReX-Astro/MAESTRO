!!  Create a 1-d hydrostatic, isoentropic atmosphere given the temperature,
!!  composition, and base density.  This version allows for an entropy 
!!  decrease below a height base_height, to make the lower region 
!!  convectively stable.
!!
!!  The model is placed into HSE by the following differencing:
!!
!!   (1/dr) [ <P>_i - <P>_{i-1} ] = (1/2) [ <rho>_i + <rho>_{i-1} ] g
!!
!!   We can take <P>_base and <rho>_base as given, and use the above
!!   to find all others.
!!
!!   This will be iterated over in tandem with the EOS call,
!!   P(i-1) = P_eos(rho(i-1), T(i-1), X(i-1)
!!

program init_1d
 
  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module, only: eos_input_rt, eos, eos_init
  use eos_type_module, only: eos_t
  use extern_probin_module, only: use_eos_coulomb
  use network, only : nspec, network_species_index, spec_names, network_init

  implicit none

  integer :: i, n

  character(len=128) :: params_file
      
  real (kind=dp_t) :: temp_base, dens_base
  real (kind=dp_t), DIMENSION(nspec) :: xn_base

  real (kind=dp_t), allocatable :: xzn_hse(:)
  real (kind=dp_t), allocatable :: model_hse(:,:)

  real :: A, B

  integer :: nx

  ! define convenient indices for the scalars
  integer, parameter :: nvar = 3 + nspec
  integer, parameter :: idens = 1, &
                        itemp = 2, &
                        ipres = 3, &
                        ispec = 4

  ! we'll get the composition from the network module
  integer :: ic12, io16, iash

  real (kind=dp_t) :: xmin, xmax, dCoord

  real (kind=dp_t) :: dens_zone, temp_zone, pres_zone, entropy
  real (kind=dp_t) :: dpd, dpt, dsd, dst

  real (kind=dp_t) :: p_want, drho, dtemp, delx
  real (kind=dp_t), allocatable :: entropy_store(:), entropy_want(:)
  
  real (kind=dp_t) :: g_zone

  real (kind=dp_t), parameter :: TOL = 1.e-10

  integer, parameter :: MAX_ITER = 250

  integer :: iter

  logical :: converged_hse, fluff

  real (kind=dp_t), dimension(nspec) :: xn

  real (kind=dp_t) :: low_density_cutoff, temp_fluff, smallx

  real (kind=dp_t) :: entropy_jump, entropy_end, base_height

  logical :: do_coulomb

  integer :: index_base

  real (kind=dp_t) :: slope

  logical :: isentropic 

  character (len=256) :: outfile
  character (len=8) :: num

  real (kind=dp_t) :: max_hse_error, dpdr, rhog

  integer :: narg

  type (eos_t) :: eos_state

  namelist /params/ nx, xmin, xmax, base_height, entropy_jump, entropy_end, &
       dens_base, temp_base, g_zone, low_density_cutoff, temp_fluff, do_coulomb

  ! determine if we specified a runtime parameters file or use the
  ! default
  narg = command_argument_count()

  if (narg == 0) then
     params_file = "_params"
  else
     call get_command_argument(1, value = params_file)
  endif


  ! define the defaults parameters for this model
  nx = 640

  xmin = 0_dp_t
  xmax = 3.6d8

  base_height = 1.d8

  entropy_jump = 3.d0
  entropy_end  = 6.d0

  dens_base = 2.6d9
  temp_base = 7.d8

  g_zone = -1.5d10

  low_density_cutoff = 1.d-4
  temp_fluff = 1.d7


  ! check the namelist for any changed parameters
  open(unit=11, file=trim(params_file), status="old", action="read")
  read(unit=11, nml=params)
  close(unit=11)




  ! EOS stuff
  smallx = 1.d-10

  ! this comes in from extern_probin_module -- override here if desired
  use_eos_coulomb = do_coulomb


  ! initialize the EOS and network
  call eos_init()
  call network_init()

  
  ! get the species indices
  ic12  = network_species_index("carbon-12")
  io16  = network_species_index("oxygen-16")

  ! ash can be either Mg24 or 'ash', depending on the network
  iash = network_species_index("magnesium-24")
  if (iash < 0) then
     iash = network_species_index("ash")
  endif

  if (ic12 < 0 .or. io16 < 0 .or. iash < 0) then
     call bl_error("ERROR: species not defined")
  endif

  xn_base(ic12) = 0.3_dp_t
  xn_base(io16) = 0.7_dp_t
  xn_base(iash) = 0.0_dp_t
  

!-----------------------------------------------------------------------------
! Create a 1-d uniform grid that is identical to the mesh that we are
! mapping onto, and then we want to force it into HSE on that mesh.
!-----------------------------------------------------------------------------

! allocate storage
  allocate(xzn_hse(nx))
  allocate(model_hse(nx,nvar))
  allocate(entropy_want(nx))
  allocate(entropy_store(nx))

! compute the coordinates of the new gridded function
  dCoord = (xmax - xmin) / dble(nx)

  do i = 1, nx
     xzn_hse(i) = xmin + (dble(i) - 0.5_dp_t)*dCoord
  enddo


  index_base = -1

! find the index of the base height 
  do i = 1, nx
     if (xzn_hse(i) >= base_height) then
        index_base = i
        exit
     endif
  enddo

  if (index_base == -1) then
     print *, 'ERROR: base_height not found on grid'
     call bl_error('ERROR: invalid base_height')
  endif



!-----------------------------------------------------------------------------
! put the model onto our new uniform grid
!-----------------------------------------------------------------------------

  fluff = .false.

  ! all the material above base_height will have constant entropy
  eos_state%T     = temp_base
  eos_state%rho   = dens_base
  eos_state%xn(:) = xn_base(:)

  print *, eos_state%T, eos_state%rho, eos_state%xn(:)

  call eos(eos_input_rt, eos_state)

  ! make the inital guess be completely uniform
  model_hse(:,idens) = eos_state%rho
  model_hse(:,itemp) = eos_state%T
  model_hse(:,ipres) = eos_state%p

  do i = 1, nspec
     model_hse(:,ispec-1+i) = eos_state%xn(i)
  enddo

  entropy = eos_state%s

  ! set the desired entropy profile
  do i = index_base, nx

     ! entropy is constant
     entropy_want(i) = entropy
     entropy_store(i) = entropy_want(i)
  enddo


  slope = (entropy/entropy_jump - entropy/entropy_end)/ &
       (xzn_hse(index_base) - xzn_hse(1))

  do i = index_base-1, 1, -1

     ! entropy gradient
     entropy_want(i) = slope*(xzn_hse(i) - xzn_hse(index_base-1)) + entropy/entropy_jump
     entropy_store(i) = entropy_want(i)
  enddo


!-----------------------------------------------------------------------------
! HSE + entropy solve
!-----------------------------------------------------------------------------

! the HSE state will be done putting creating an isentropic state until
! the temperature goes below temp_fluff -- then we will do isothermal.
! also, once the density goes below low_density_cutoff, we stop HSE

  isentropic = .true.

  !---------------------------------------------------------------------------
  ! integrate up
  !---------------------------------------------------------------------------
  do i = index_base+1, nx

     delx = xzn_hse(i) - xzn_hse(i-1)

     ! as the initial guess for the temperature and density, use the previous
     ! zone
     dens_zone = model_hse(i-1,idens)
     temp_zone = model_hse(i-1,itemp)
     xn(:) = model_hse(i,ispec:nvar)

     
     !-----------------------------------------------------------------------
     ! iteration loop
     !-----------------------------------------------------------------------

     ! start off the Newton loop by saying that the zone has not converged
     converged_hse = .FALSE.

     if (.not. fluff) then

        do iter = 1, MAX_ITER

           if (isentropic) then

              ! get the pressure we want from the HSE equation, just the
              ! zone below the current.  Note, we are using an average of
              ! the density of the two zones as an approximation of the
              ! interface value -- this means that we need to iterate for
              ! find the density and pressure that are consistent
              
              ! furthermore, we need to get the entropy that we need,
              ! which will come from adjusting the temperature in
              ! addition to the density.
        
              ! HSE differencing
              p_want = model_hse(i-1,ipres) + &
                   delx*0.5*(dens_zone + model_hse(i-1,idens))*g_zone

              ! now we have two functions to zero:
              !   A = p_want - p(rho,T)
              !   B = entropy_want - s(rho,T)
              ! We use a two dimensional Taylor expansion and find the deltas 
              ! for both density and temperature

              
              ! now we know the pressure and the entropy that we want, so we 
              ! need to find the temperature and density through a two 
              ! dimensional root find

              ! (t, rho) -> (p, s)
              eos_state%T     = temp_zone
              eos_state%rho   = dens_zone
              eos_state%xn(:) = xn(:)

              call eos(eos_input_rt, eos_state)
              
              entropy = eos_state%s
              pres_zone = eos_state%p

              dpt = eos_state%dpdt
              dpd = eos_state%dpdr
              dst = eos_state%dsdt
              dsd = eos_state%dsdr

              A = p_want - pres_zone
              B = entropy_want(i) - entropy
              
              dtemp = ((dsd/(dpd-0.5*delx*g_zone))*A - B)/ &
                   (dsd*dpt/(dpd -0.5*delx*g_zone) - dst)
              
              drho = (A - dpt*dtemp)/(dpd - 0.5*delx*g_zone)
              
              dens_zone = max(0.9_dp_t*dens_zone, &
                   min(dens_zone + drho, 1.1_dp_t*dens_zone))

              temp_zone = max(0.9_dp_t*temp_zone, &
                   min(temp_zone + dtemp, 1.1_dp_t*temp_zone))


              ! check if the density falls below our minimum cut-off -- 
              ! if so, floor it
              if (dens_zone < low_density_cutoff) then
                 
                 dens_zone = low_density_cutoff
                 temp_zone = temp_fluff
                 converged_hse = .TRUE.
                 fluff = .TRUE.
                 exit
              
              endif
                 
              ! if (A < TOL .and. B < ETOL) then
              if (abs(drho) < TOL*dens_zone .and. abs(dtemp) < TOL*temp_zone) then
                 converged_hse = .TRUE.
                 exit
              endif
        
           else

              ! do isothermal
              p_want = model_hse(i-1,ipres) + &
                   delx*0.5*(dens_zone + model_hse(i-1,idens))*g_zone
         
              temp_zone = temp_fluff

              ! (t, rho) -> (p)
              eos_state%T     = temp_zone
              eos_state%rho   = dens_zone
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
                 temp_zone = temp_fluff
                 converged_hse = .TRUE.
                 fluff = .TRUE.
                 exit
              
              endif

           endif

        enddo
        
        if (.NOT. converged_hse) then
           
           print *, 'Error zone', i, ' did not converge in init_1d'
           print *, 'integrate up'
           print *, dens_zone, temp_zone
           print *, p_want, entropy_want(i), entropy
           print *, drho, dtemp
           call bl_error('Error: HSE non-convergence')
           
        endif

        if (temp_zone < temp_fluff) then
           temp_zone = temp_fluff
           isentropic = .false.
        endif

     else
        dens_zone = low_density_cutoff
        temp_zone = temp_fluff
     endif


     ! call the EOS one more time for this zone and then go on to the next
     ! (t, rho) -> (p)
     eos_state%T     = temp_zone
     eos_state%rho   = dens_zone
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rt, eos_state)

     pres_zone = eos_state%p
     
     ! update the thermodynamics in this zone
     model_hse(i,idens) = dens_zone
     model_hse(i,itemp) = temp_zone
     model_hse(i,ipres) = pres_zone
     entropy_store(i) = entropy

     ! to make this process converge faster, set the density in the next zone to
     ! the density in this zone
     ! model_hse(i+1,idens) = dens_zone

  enddo


  !---------------------------------------------------------------------------
  ! integrate down
  !---------------------------------------------------------------------------
  do i = index_base-1, 1, -1

     delx = xzn_hse(i+1) - xzn_hse(i)

     ! as the initial guess for the temperature and density, use the previous
     ! zone
     dens_zone = model_hse(i+1,idens)
     temp_zone = model_hse(i+1,itemp)
     xn(:) = model_hse(i,ispec:nvar)


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
        
        ! furthermore, we need to get the entropy that we need,
        ! which will come from adjusting the temperature in
        ! addition to the density.
        
        ! HSE differencing
        p_want = model_hse(i+1,ipres) - &
             delx*0.5*(dens_zone + model_hse(i+1,idens))*g_zone

        ! now we have two functions to zero:
        !   A = p_want - p(rho,T)
        !   B = entropy_want - s(rho,T)
        ! We use a two dimensional Taylor expansion and find the deltas 
        ! for both density and temperature
        
        
        ! now we know the pressure and the entropy that we want, so we 
        ! need to find the temperature and density through a two 
        ! dimensional root find
        ! (t, rho) -> (p, s)
        eos_state%T     = temp_zone
        eos_state%rho   = dens_zone
        eos_state%xn(:) = xn(:)

        call eos(eos_input_rt, eos_state)
        
        entropy = eos_state%s
        pres_zone = eos_state%p
        
        dpt = eos_state%dpdt
        dpd = eos_state%dpdr
        dst = eos_state%dsdt
        dsd = eos_state%dsdr
              
        A = p_want - pres_zone
        B = entropy_want(i) - entropy
              
        dtemp = ((dsd/(dpd+0.5*delx*g_zone))*A - B)/ &
             (dsd*dpt/(dpd +0.5*delx*g_zone) - dst)
        
        drho = (A - dpt*dtemp)/(dpd + 0.5*delx*g_zone)
        
        dens_zone = max(0.9_dp_t*dens_zone, &
             min(dens_zone + drho, 1.1_dp_t*dens_zone))

        temp_zone = max(0.9_dp_t*temp_zone, &
             min(temp_zone + dtemp, 1.1_dp_t*temp_zone))
        
                        
        if (abs(drho) < TOL*dens_zone .and. abs(dtemp) < TOL*temp_zone) then
           converged_hse = .TRUE.
           exit
        endif
        

     enddo
        
     if (.NOT. converged_hse) then
        
        print *, 'Error zone', i, ' did not converge in init_1d'
        print *, 'integrate down'
        print *, dens_zone, temp_zone
        print *, p_want, entropy_want(i), entropy
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
     model_hse(i,idens) = dens_zone
     model_hse(i,itemp) = temp_zone
     model_hse(i,ipres) = pres_zone
     entropy_store(i) = entropy


  enddo

  write(num,'(i8)') nx
  outfile = "model.hse." // trim(adjustl(num))

  open (unit=50, file=outfile, status="unknown")

  write (50,1001) "# npts = ", nx
  write (50,1001) "# num of variables = ", nvar
  write (50,1002) "# density"
  write (50,1002) "# temperature"
  write (50,1002) "# pressure"

  do n = 1, nspec
     write (50, 1003) "# ", spec_names(n)
  enddo

1000 format (1x, 12(g26.16, 1x))
1001 format (a, i5)
1002 format (a)
1003 format (a,a)

  do i = 1, nx

     write (50,1000) xzn_hse(i), model_hse(i,idens), model_hse(i,itemp), model_hse(i,ipres), &
          (model_hse(i,ispec-1+n), n=1,nspec)

  enddo

  ! compute the maximum HSE error
  max_hse_error = -1.d30

  do i = 2, nx-1
     dpdr = (model_hse(i,ipres) - model_hse(i-1,ipres))/delx
     rhog = HALF*(model_hse(i,idens) + model_hse(i-1,idens))*g_zone

     if (dpdr /= ZERO .and. model_hse(i+1,idens) > low_density_cutoff) then
        max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(dpdr))
     endif
  enddo

  print *, 'maximum HSE error = ', max_hse_error
  print *, ' '

  close (unit=50)

end program init_1d






