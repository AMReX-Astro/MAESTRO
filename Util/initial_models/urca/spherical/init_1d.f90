!! generate an initial model for spherical geometry with a
!! uniform composition.  Here we take a base density and temperature
!! and use HSE and constant entropy to generate the model.

program init_1d

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module, only: eos, eos_input_rt, eos_init
  use eos_type_module, only: eos_t
  use extern_probin_module, only: use_eos_coulomb
  use network
  use fundamental_constants_module, only: Gconst
  use urca_composition_module

  implicit none

  integer :: i, j, n

  integer :: nx

  real (kind=dp_t) :: temp_base, dens_base
  real (kind=dp_t), DIMENSION(nspec) :: xn_base

  real (kind=dp_t), allocatable :: xzn_hse(:), xznl(:), xznr(:)
  real (kind=dp_t), allocatable :: model_hse(:,:), M_enclosed(:)

  real (kind=dp_t) :: A, B, dAdT, dAdrho, dBdT, dBdrho

  ! define convenient indices for the scalars
  integer, parameter :: nvar = 3 + nspec
  integer, parameter :: idens = 1, &
                        itemp = 2, &
                        ipres = 3, &
                        ispec = 4

  real (kind=dp_t), parameter :: M_sun = 1.9891e33

  ! we'll get the composition indices from the network module
  integer, save :: ihe4, ic12, io16, ine20, ine23, ina23, img23

  integer :: narg
  character(len=128) :: params_file

  real (kind=dp_t), save :: xmin, xmax, dCoord

  real (kind=dp_t) :: dens_zone, temp_zone, pres_zone, entropy
  real (kind=dp_t) :: dpd, dpt, dsd, dst

  real (kind=dp_t) :: p_want, drho, dtemp, delx
  real (kind=dp_t), allocatable :: entropy_store(:), entropy_want(:)

  real (kind=dp_t) :: g_zone

  real (kind=dp_t) :: TOL = 1.e-11

  integer :: MAX_ITER = 1000000

  integer :: iter

  logical :: converged_hse, fluff

  real (kind=dp_t), dimension(nspec) :: xn

  real (kind=dp_t), save :: low_density_cutoff, smallx, dens_conv_zone, M_conv_zone

  logical :: test_hse_convergence = .true.

  real (kind=dp_t), save :: temp_before_fluff, temp_fluff

  character (len=128) :: fluff_type = "continuous"

  logical :: isentropic

  character (len=256) :: outfile, prefix
  character (len=8) num

  real (kind=dp_t) :: max_hse_error, dpdr, rhog

  real (kind=dp_t) :: dtol_fac = 0.000001_dp_t

  integer :: i_conv, i_fluff

  type (eos_t) :: eos_state

  namelist /params/ nx, dens_base, temp_base, &
       low_density_cutoff, dens_conv_zone, M_conv_zone, temp_fluff, &
       test_hse_convergence, TOL, MAX_ITER, dtol_fac, &
       xmin, xmax, &
       fluff_type, &
       c12_in, c12_out, o16_in, o16_out, &
       ne23_in, ne23_out, na23_in, na23_out, &
       urca_23_dens, urca_shell_type, shell_atan_kappa, &
       prefix


  ! determine if we specified a runtime parameters file or use the default
  narg = command_argument_count()

  if (narg == 0) then
     params_file = "_params"
  else
     call get_command_argument(1, value = params_file)
  endif


  ! define the defaults parameters for this model
  nx = 1280

  smallx = 1.d-10

  xmin = 0_dp_t
  xmax = 5.d8

  dens_base = 2.6d9
  temp_base = 6.d8

  dens_conv_zone = -1.d0
  M_conv_zone = 2.0d0


  low_density_cutoff =1.d-4
  temp_fluff = 1.d7

  c12_in   = 0.0_dp_t
  c12_out  = 0.0_dp_t
  o16_in   = 0.0_dp_t
  o16_out  = 0.0_dp_t
  ne23_in  = 0.0_dp_t
  ne23_out = 0.0_dp_t
  na23_in  = 0.0_dp_t
  na23_out = 0.0_dp_t

  prefix = "spherical"


  ! check the namelist for any changed parameters
  open(unit=11, file=params_file, status="old", action="read")
  read(unit=11, nml=params)
  close(unit=11)

  ! Initialize the temperature interior to fluff
  ! (used to calculate fluff temperature for various values of fluff_type)
  temp_before_fluff = temp_fluff

  ! initialize the EOS and network

  ! use_eos_coulomb comes in from extern_probin_module -- override
  ! here if desired
  use_eos_coulomb = .true.
  call eos_init()
  call network_init()

  ! Initialize the composition module
  call init_urca_composition()

!-----------------------------------------------------------------------------
! Create a 1-d uniform grid that is identical to the mesh that we are
! mapping onto, and then we want to force it into HSE on that mesh.
!-----------------------------------------------------------------------------

! allocate storage
  allocate(xzn_hse(nx))
  allocate(xznl(nx))
  allocate(xznr(nx))
  allocate(model_hse(nx,nvar))
  allocate(M_enclosed(nx))
  allocate(entropy_want(nx))
  allocate(entropy_store(nx))

! compute the coordinates of the new gridded function
  dCoord = (xmax - xmin) / dble(nx)

  do i = 1, nx
     xznl(i) = xmin + (dble(i) - 1.0_dp_t)*dCoord
     xznr(i) = xmin + (dble(i))*dCoord
     xzn_hse(i) = 0.5_dp_t*(xznl(i) + xznr(i))
  enddo

  fluff = .false.

  ! call the EOS one more time for this zone and then go on to the next
  eos_state%T     = temp_base
  eos_state%rho   = dens_base
  call set_urca_composition(dens_base, xn_base)
  eos_state%xn(:) = xn_base(:)

  ! (t, rho) -> (p, s)
  call eos(eos_input_rt, eos_state)

  ! make the initial guess be completely uniform
  model_hse(:,idens) = eos_state%rho
  model_hse(:,itemp) = eos_state%T
  model_hse(:,ipres) = eos_state%p

  do i = 1, nspec
     model_hse(:,ispec-1+i) = eos_state%xn(i)
  enddo

  entropy_want(:) = eos_state%s


  ! keep track of the mass enclosed below the current zone
  M_enclosed(1) = FOUR3RD*M_PI*(xznr(1)**3 - xznl(1)**3)*model_hse(1,idens)


!-----------------------------------------------------------------------------
! HSE + entropy solve
!-----------------------------------------------------------------------------

  isentropic = .true.

  do i = 2, nx

     delx = xzn_hse(i) - xzn_hse(i-1)

     ! as the initial guess for the temperature and density, use the previous
     ! zone
     dens_zone = model_hse(i-1,idens)
     temp_zone = model_hse(i-1,itemp)
     call set_urca_composition(dens_zone, xn)

     g_zone = -Gconst*M_enclosed(i-1)/(xznl(i)*xznl(i))


     !-----------------------------------------------------------------------
     ! iteration loop
     !-----------------------------------------------------------------------

     ! start off the Newton loop by saying that the zone has not converged
     converged_hse = .FALSE.


     if (.not. fluff) then

        do iter = 1, MAX_ITER

           if (isentropic) then

              p_want = model_hse(i-1,ipres) + &
                   delx*0.5_dp_t*(dens_zone + model_hse(i-1,idens))*g_zone

              ! now we have two functions to zero:
              !   A = p_want - p(rho,T)
              !   B = entropy_want - s(rho,T)
              ! We use a two dimensional Taylor expansion and find the deltas
              ! for both density and temperature

              eos_state%T     = temp_zone
              eos_state%rho   = dens_zone
              call set_urca_composition(dens_zone, xn)
              eos_state%xn(:) = xn(:)

              ! (t, rho) -> (p, s)
              call eos(eos_input_rt, eos_state)

              entropy = eos_state%s
              pres_zone = eos_state%p

              dpt = eos_state%dpdt
              dpd = eos_state%dpdr
              dst = eos_state%dsdt
              dsd = eos_state%dsdr

              A = p_want - pres_zone
              B = entropy_want(i) - entropy

              dAdT = -dpt
              dAdrho = 0.5d0*delx*g_zone - dpd
              dBdT = -dst
              dBdrho = -dsd

              dtemp = (B - (dBdrho/dAdrho)*A)/ &
                   ((dBdrho/dAdrho)*dAdT - dBdT)

              drho = -(A + dAdT*dtemp)/dAdrho

              dens_zone = max((ONE-dtol_fac)*dens_zone, &
                              min(dens_zone + drho, (ONE+dtol_fac)*dens_zone))

              temp_zone = max((ONE-dtol_fac)*temp_zone, &
                              min(temp_zone + dtemp, (ONE+dtol_fac)*temp_zone))


              ! check if the density falls below our minimum cut-off --
              ! if so, floor it
              if (dens_zone < low_density_cutoff) then

                 i_fluff = i

                 dens_zone = low_density_cutoff

                 temp_before_fluff = model_hse(i-1,itemp)

                 call get_fluff_temperature(temp_zone, temp_fluff, temp_before_fluff, fluff_type)

                 converged_hse = .TRUE.
                 fluff = .TRUE.
                 exit

              endif

              if (dens_zone < dens_conv_zone .and. isentropic) then

                 i_conv = i
                 isentropic = .false.

              endif

              if (test_hse_convergence) then

                 ! Check to see if HSE is satisfied to relative tolerance TOL

                 ! Call EOS again using the new values of dens_zone and temp_zone
                 ! to get pres_zone ...
                 eos_state%T     = temp_zone
                 eos_state%rho   = dens_zone
                 call set_urca_composition(dens_zone, xn)
                 eos_state%xn(:) = xn(:)

                 ! (t, rho) -> (p, s)
                 call eos(eos_input_rt, eos_state)

                 pres_zone = eos_state%p

                 ! Compute dpdr
                 dpdr = (pres_zone - model_hse(i-1,ipres))/delx

                 ! Compute rho*g
                 rhog = HALF*(dens_zone + model_hse(i-1,idens))*g_zone

                 if (abs(dpdr - rhog) < TOL*abs(dpdr)) then
                    converged_hse = .TRUE.
                    exit
                 endif

              else

                 ! Check to see if drho and dtemp are under the relative tolerance TOL
                 if (abs(drho) < TOL*dens_zone .and. abs(dtemp) < TOL*temp_zone) then
                    converged_hse = .TRUE.
                    exit
                 endif

              endif

           else

              ! do isothermal
              p_want = model_hse(i-1,ipres) + &
                   delx*0.5*(dens_zone + model_hse(i-1,idens))*g_zone

              temp_zone = model_hse(i-1,itemp)


              eos_state%T     = temp_zone
              eos_state%rho   = dens_zone
              call set_urca_composition(dens_zone, xn)
              eos_state%xn(:) = xn(:)

              ! (t, rho) -> (p, s)
              call eos(eos_input_rt, eos_state)

              entropy = eos_state%s
              pres_zone = eos_state%p

              dpd = eos_state%dpdr

              drho = (p_want - pres_zone)/(dpd - 0.5*delx*g_zone)

              dens_zone = max(0.9*dens_zone, &
                   min(dens_zone + drho, 1.1*dens_zone))

              if (dens_zone < low_density_cutoff) then

                 i_fluff = i

                 dens_zone = low_density_cutoff

                 temp_before_fluff = model_hse(i-1,itemp)

                 call get_fluff_temperature(temp_zone, temp_fluff, temp_before_fluff, fluff_type)

                 converged_hse = .TRUE.
                 fluff = .TRUE.
                 exit

              endif

              if (test_hse_convergence) then

                 ! Check to see if HSE is satisfied to relative tolerance TOL

                 ! Call EOS again using the new values of dens_zone and temp_zone
                 ! to get pres_zone ...
                 eos_state%T     = temp_zone
                 eos_state%rho   = dens_zone
                 call set_urca_composition(dens_zone, xn)
                 eos_state%xn(:) = xn(:)

                 ! (t, rho) -> (p, s)
                 call eos(eos_input_rt, eos_state)

                 pres_zone = eos_state%p

                 ! Compute dpdr
                 dpdr = (pres_zone - model_hse(i-1,ipres))/delx

                 ! Compute rho*g
                 rhog = HALF*(dens_zone + model_hse(i-1,idens))*g_zone

                 if (abs(dpdr - rhog) < TOL*abs(dpdr)) then
                    converged_hse = .TRUE.
                    exit
                 endif

              else

                 ! Check to see if drho is under the relative tolerance TOL
                 if (abs(drho) < TOL*dens_zone) then
                    converged_hse = .TRUE.
                    exit
                 endif

              endif

           endif

        enddo

        if (.NOT. converged_hse) then

           print *, 'Error zone', i, ' did not converge in init_1d'
           print *, 'integrate up'
           print *, 'trial density (dens_zone): ', dens_zone
           print *, 'trial temperature (temp_zone): ', temp_zone
           print *, 'current pressure (pres_zone): ', pres_zone
           print *, 'desired pressure (p_want): ', p_want
           print *, 'desired pressure rel tolerance: ', abs(pres_zone-p_want)/pres_zone
           print *, 'delta density (drho): ', drho
           print *, 'delta temperature (dtemp): ', dtemp
           print *, 'pressure gradient (dpdr): ', dpdr
           print *, 'density * gravitational accel. (rhog): ', rhog
           print *, 'difference... abs(dpdr-rhog): ', abs(dpdr-rhog)
           print *, 'target... TOL*abs(dpdr): ', TOL*abs(dpdr)
           call bl_error('Error: HSE non-convergence')

        endif

        if (temp_zone < temp_fluff) then
           if (fluff_type .eq. "constant") then
              call get_fluff_temperature(temp_zone, temp_fluff, temp_before_fluff, fluff_type)
           endif
           isentropic = .false.
        endif

     else
        dens_zone = low_density_cutoff
        call get_fluff_temperature(temp_zone, temp_fluff, temp_before_fluff, fluff_type)
     endif


     ! call the EOS one more time for this zone and then go on to the next
     eos_state%T     = temp_zone
     eos_state%rho   = dens_zone
     call set_urca_composition(dens_zone, xn)
     eos_state%xn(:) = xn(:)

     ! (t, rho) -> (p, s)
     call eos(eos_input_rt, eos_state)

     pres_zone = eos_state%p

     dpd = eos_state%dpdr

     ! update the thermodynamics in this zone
     model_hse(i,idens) = dens_zone
     model_hse(i,itemp) = temp_zone
     model_hse(i,ipres) = pres_zone
     do j = 1, nspec
        model_hse(i,ispec-1+j) = eos_state%xn(j)
     enddo

     print *, i, dens_zone, temp_zone

     M_enclosed(i) = M_enclosed(i-1) + &
          FOUR3RD*M_PI*(xznr(i) - xznl(i))* &
            (xznr(i)**2 +xznl(i)*xznr(i) + xznl(i)**2)*model_hse(i,idens)

     if (M_enclosed(i) > M_conv_zone*M_sun .and. isentropic) then

        i_conv = i
        isentropic = .false.

     endif

  enddo


  print *, 'mass = ', M_enclosed(nx)/M_sun

  write(num,'(i8)') nx
  outfile = trim(prefix) // ".hse." // trim(adjustl(num))



  open (unit=50, file=outfile, status="unknown")

  write (50,1001) "# npts = ", nx
  write (50,1001) "# num of variables = ", nvar
  write (50,1002) "# density"
  write (50,1002) "# temperature"
  write (50,1002) "# pressure"

  do n = 1, nspec
     write (50,1003) "# ", spec_names(n)
  enddo

1000 format (1x, 12(g26.16, 1x))
1001 format(a, i5)
1002 format(a)
1003 format(a,a)

  do i = 1, nx

     write (50,1000) xzn_hse(i), model_hse(i,idens), model_hse(i,itemp), model_hse(i,ipres), &
          (model_hse(i,ispec-1+n), n=1,nspec)

  enddo


  ! compute the maximum HSE error
  max_hse_error = -1.d30

  do i = 2, nx-1
     g_zone = -Gconst*M_enclosed(i-1)/xznr(i-1)**2
     dpdr = (model_hse(i,ipres) - model_hse(i-1,ipres))/delx
     rhog = HALF*(model_hse(i,idens) + model_hse(i-1,idens))*g_zone

     if (dpdr /= ZERO .and. model_hse(i,idens) > low_density_cutoff) then
        max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(dpdr))
     endif

     if (dpdr /= 0) then
        print *, i, real(model_hse(i,idens)), real(model_hse(i,itemp)), real(dpdr), real(rhog), abs(dpdr - rhog)/abs(dpdr)
     endif
  enddo

  print *, 'maximum HSE error = ', max_hse_error
  print *, ' '


  print *, 'total mass = ', M_enclosed(i_fluff)/M_sun
  print *, 'convective zone mass = ', M_enclosed(i_conv)/M_sun



  close (unit=50)

end program init_1d


subroutine get_fluff_temperature(temp, temp_fluff, temp_previous, fluff_type)

  ! There are 2 kinds of fluff temperature handling:
  ! 1) fluff_type = "constant"   : fluff is at the temperature temp_fluff in the inputs
  ! 2) fluff_type = "continuous" : fluff is at the same temperature as the material
  !                                immediately interior to the fluff.

  use bl_types, only: dp_t
  use bl_error_module

  implicit none

  real (kind=dp_t), intent(out)   :: temp
  real (kind=dp_t), intent(in)    :: temp_fluff, temp_previous
  character (len=128), intent(in) :: fluff_type

  if (fluff_type .eq. "constant") then
     temp = temp_fluff
  else if (fluff_type .eq. "continuous") then
     temp = temp_previous
  else
     call bl_error("ERROR: invalid fluff_type")
  end if

end subroutine get_fluff_temperature
