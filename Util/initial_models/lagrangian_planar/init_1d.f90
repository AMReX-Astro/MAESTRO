!!  The model is placed into HSE by the following differencing:
!!
!!   (1/dr) [ <P>_i - <P>_{i-1} ] = (1/2) [ <rho>_i + <rho>_{i-1} ] g
!!
!!  This will be iterated over in tandem with the EOS call,
!!  P(i-1) = P_eos(rho(i-1), T(i-1), X(i-1)
!!

program init_1d
 
  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module, only: eos_input_rt, eos, eos_init
  use eos_type_module
  use extern_probin_module, only: use_eos_coulomb
  use network, only : nspec, network_species_index, spec_names, network_init
  use fundamental_constants_module, only: Gconst
  use model_parser_module

  implicit none

  integer :: i, n

  character(len=128) :: params_file

  real (kind=dp_t), DIMENSION(nspec) :: xn

  real (kind=dp_t), allocatable :: xzn_hse(:), xznl_hse(:), xznr_hse(:)
  real (kind=dp_t), allocatable :: model_hse(:,:)

  real :: A

  integer ::nx

  integer :: lun1, lun2

  real (kind=dp_t) :: xmin, xmax, dCoord

  real (kind=dp_t) :: dens_zone, temp_zone, pres_zone, entropy
  real (kind=dp_t) :: dpd

  real (kind=dp_t) :: p_want, drho, dtemp, delx
  
  real (kind=dp_t) :: g_zone, g_const, M_enclosed
  logical :: do_invsq_grav

  real (kind=dp_t), parameter :: TOL = 1.e-10

  integer, parameter :: MAX_ITER = 250

  integer :: iter

  logical :: converged_hse, fluff

  real (kind=dp_t) :: low_density_cutoff, temp_cutoff, smallx, max_T
  real (kind=dp_t) :: model_shift

  integer :: index_base

  character (len=256) :: outfile, outfile2
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
                    low_density_cutoff, model_prefix, model_shift

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

  do_invsq_grav = .false.
  g_const = -1.0

  model_prefix = "model"

  low_density_cutoff = 1.d-4
  temp_cutoff = 1.d6

  smallx = 1.d-10


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


!-----------------------------------------------------------------------------
! Create a 1-d uniform grid that is identical to the mesh that we are
! mapping onto, and then we want to force it into HSE on that mesh.
!-----------------------------------------------------------------------------

  ! allocate storage
  allocate(xzn_hse(nx))
  allocate(xznl_hse(nx))
  allocate(xznr_hse(nx))
  allocate(model_hse(nx,nvars_model))


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

     call eos(eos_input_rt, eos_state, .false.)

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

  call eos(eos_input_rt, eos_state, .false.)

  model_hse(index_base,ipres_model) = eos_state%p
  
  
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
     if (do_invsq_grav) then
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

        call eos(eos_input_rt, eos_state, .false.)
              
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

     call eos(eos_input_rt, eos_state, .false.)

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
     if (do_invsq_grav) then
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

        call eos(eos_input_rt, eos_state, .false.)
        
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

     call eos(eos_input_rt, eos_state, .false.)

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

  open (newunit=lun1, file=outfile, status="unknown")
  open (newunit=lun2, file=outfile2, status="unknown")  

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


  write (lun2,1001), "# npts = ", nx
  write (lun2,1001), "# num of variables = ", 2
  write (lun2,1002), "# entropy"
  write (lun2,1002), "# c_s"
  
  ! test: bulk EOS call -- Maestro will do this once we are mapped, so make
  ! sure that we are in HSE with updated thermodynamics
  do i = 1, nx
     eos_state%rho = model_hse(i,idens_model)
     eos_state%T = model_hse(i,itemp_model)
     eos_state%xn(:) = model_hse(i,ispec_model:ispec_model-1+nspec)

     call eos(eos_input_rt, eos_state, .false.)

     model_hse(i,ipres_model) = eos_state%p

     write (lun2,1000), xzn_hse(i), eos_state%s, eos_state%cs
  enddo
  
  ! compute the maximum HSE error
  max_hse_error = -1.d30

  do i = 2, nx-1

     ! compute the gravitation acceleration at the lower edge
     if (do_invsq_grav) then
        g_zone = -Gconst*M_enclosed/xznl_hse(i)**2
     else
        g_zone = g_const
     endif

     dpdr = (model_hse(i,ipres_model) - model_hse(i-1,ipres_model))/delx
     rhog = HALF*(model_hse(i,idens_model) + model_hse(i-1,idens_model))*g_zone

     if (dpdr /= ZERO .and. model_hse(i+1,idens_model) > low_density_cutoff) then
        max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(dpdr))
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


