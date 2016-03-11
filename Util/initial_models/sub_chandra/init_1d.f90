!! generate an initial model for an arbitrary-mass, isothermal C WD 
!! with an isentropic He envelope on the surface.

program init_1d

  use bl_types
  use bl_constants_module
  use bl_error_module
  use extern_probin_module, only: use_eos_coulomb
  use eos_module, only: eos_input_rt, eos, eos_init
  use eos_type_module
  use network
  use fundamental_constants_module, only: Gconst
  use f2kcli

  implicit none

  integer :: i, n

  character(len=128) :: params_file
      
  real (kind=dp_t) :: temp_core, temp_base, delta
  real (kind=dp_t), DIMENSION(nspec) :: xn_core, xn_he

  logical :: mixed_co_wd

  real (kind=dp_t), allocatable :: xzn_hse(:), xznl(:), xznr(:)
  real (kind=dp_t), allocatable :: model_hse(:,:), M_enclosed(:)
  real (kind=dp_t), allocatable :: cs_hse(:), s_hse(:)

  real (kind=dp_t) :: rho_c, rho_c_old, mass_wd, mass_wd_old, drho_c
  real (kind=dp_t) :: rho_he, rho_he_old, mass_he, mass_he_old, drho_he

  real (kind=dp_t) :: slope_T, slope_xn(nspec)

  real (kind=dp_t) :: A, B, dAdT, dAdrho, dBdT, dBdrho
  logical :: isentropic

  real (kind=dp_t) :: test

  integer :: nx

  ! define convenient indices for the scalars
  integer, parameter :: nvar = 3 + nspec
  integer, parameter :: idens = 1, &
                        itemp = 2, &
                        ipres = 3, &
                        ispec = 4

  ! we'll get the composition indices from the network module
  integer, save :: ihe4, ic12, io16

  real (kind=dp_t), save :: xmin, xmax, dCoord

  real (kind=dp_t) :: dens_zone, temp_zone, pres_zone, entropy
  real (kind=dp_t) :: dpd, dpt, dsd, dst

  real (kind=dp_t) :: p_want, drho, dtemp, delx
  real (kind=dp_t) :: entropy_base

  real (kind=dp_t) :: g_zone

  ! TOL_HSE is the tolerance used when iterating over a zone to force
  ! it into HSE by adjusting the current density (and possibly
  ! temperature).  TOL_HSE should be very small (~ 1.e-10).
  real (kind=dp_t), parameter :: TOL_HSE = 1.d-10

  ! TOL_WD_MASS is tolerance used for getting the total WD mass equal
  ! to M_tot (defined below).  It can be reasonably small, since there
  ! will always be a central density value that can give the desired
  ! WD mass on the grid we use
  real (kind=dp_t), parameter :: TOL_WD_MASS = 1.d-4

  ! TOL_HE_MASS is the tolerance used for getting the mass of the He
  ! envelope.  This should not be too small, since the values of the
  ! He envelope mass we can achieve will depend on our grid spacing.
  real (kind=dp_t), parameter :: TOL_HE_MASS = 2.d-2


  integer, parameter :: MAX_ITER = 250

  integer :: iter, iter_mass

  integer :: icutoff, ihe_layer, ihe_entropy

  logical :: converged_hse, fluff, mass_converged

  real (kind=dp_t), dimension(nspec) :: xn

  real (kind=dp_t) :: low_density_cutoff, temp_fluff, smallx, smallt

  real (kind=dp_t) :: M_tot, M_He
  real (kind=dp_t) :: solar_mass = 1.98892d33

  character (len=256) :: outfile
  character (len=8) num, mass_wd_str, mass_he_str

  real (kind=dp_t) :: max_hse_error, dpdr, rhog

  integer :: narg

  type (eos_t) :: eos_state

  namelist /params/ nx, M_tot, M_He, delta, xmin, xmax, &
       temp_core, temp_base, mixed_co_wd, low_density_cutoff, temp_fluff, smallt


  ! determine if we specified a runtime parameters file or use the default
  narg = command_argument_count()

  if (narg == 0) then
     params_file = "_params"
  else
     call get_command_argument(1, value = params_file)
  endif
  


  ! define the defaults parameters for this model
  nx = 2560

  M_tot = 0.6
  M_He  = 0.2

  delta = 1.d-6

  xmin = 0_dp_t
  xmax = 1.6e9_dp_t

  temp_core = 1.d7
  temp_base = 4.d8

  mixed_co_wd = .true.

  low_density_cutoff =1.d-4
  temp_fluff = 1.d5
  smallt = 1.d5

  ! check the namelist for any changed parameters  
  open(unit=11, file=trim(params_file), status="old", action="read")
  read(unit=11, nml=params)
  close(unit=11)


  ! convert the envelope and WD mass into solar masses
  M_tot = M_tot * solar_mass
  M_He = M_He * solar_mass


  ! this comes in via extern_probin_module, override the default if desired.
  use_eos_coulomb = .true.

  smallx = 1.d-10


  ! initialize the EOS and network
  call eos_init()
  call network_init()


  ! get the species indices
  ihe4 = network_species_index("helium-4")
  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")

  if (ihe4 < 0 .or. ic12 < 0 .or. io16 < 0) then
     call bl_error("ERROR: species not defined")
  endif

  if (mixed_co_wd) then
     xn_core(:)    = smallx
     xn_core(ic12) = 0.5_dp_t - 0.5*(nspec - 1)*smallx
     xn_core(io16) = 0.5_dp_t - 0.5*(nspec - 1)*smallx
  else
     xn_core(:)    = smallx
     xn_core(ic12) = 1.0_dp_t - (nspec - 1)*smallx
  endif


  xn_he(:)    = smallx
  xn_he(ihe4) = 1.0_dp_t - (nspec - 1)*smallx
  

  !----------------------------------------------------------------------------
  ! Create a 1-d uniform grid that is identical to the mesh that we are
  ! mapping onto, and then we want to force it into HSE on that mesh.
  !----------------------------------------------------------------------------

  ! allocate storage
  allocate(xzn_hse(nx))
  allocate(xznl(nx))
  allocate(xznr(nx))
  allocate(model_hse(nx,nvar))
  allocate(M_enclosed(nx))
  allocate(cs_hse(nx))
  allocate(s_hse(nx))

  ! compute the coordinates of the new gridded function
  dCoord = (xmax - xmin) / dble(nx)

  do i = 1, nx
     xznl(i) = xmin + (dble(i) - 1.0_dp_t)*dCoord
     xznr(i) = xmin + (dble(i))*dCoord
     xzn_hse(i) = 0.5_dp_t*(xznl(i) + xznr(i))
  enddo


  ! We don't know what WD central density will give the desired total
  ! mass, so we need to iterate over central density

  ! we will do a secant iteration.  rho_c_old is the 'old' guess for
  ! the central density and rho_c is the current guess.  After 2
  ! loops, we can start estimating the density required to yield our
  ! desired mass
  rho_c_old = -1.0_dp_t
  rho_c     = 1.e9_dp_t     ! 1.e9 is a reasonable starting WD central density

  ! rho_he_old is the old guess for the density to transition to He,
  ! where we will be isentropic, and rho_he is the currrent guess.  
  rho_he_old = -1.0_dp_t
  rho_he     = 0.5*rho_c

  mass_converged = .false.

  
  do iter_mass = 1, MAX_ITER

     print *, 'mass iter = ', iter_mass, rho_c, temp_core

     fluff = .false.

     ! we start at the center of the WD and integrate outward.  Initialize
     ! the central conditions.
     eos_state%T     = temp_core
     eos_state%rho   = rho_c
     eos_state%xn(:) = xn_core(:)

     ! (t, rho) -> (p, s)    
     call eos(eos_input_rt, eos_state)

     ! make the initial guess be completely uniform
     model_hse(:,idens) = eos_state%rho
     model_hse(:,itemp) = eos_state%T
     model_hse(:,ipres) = eos_state%p

     do i = 1, nspec
        model_hse(:,ispec-1+i) = eos_state%xn(i)
     enddo


     ! keep track of the mass enclosed below the current zone
     M_enclosed(1) = FOUR3RD*M_PI*(xznr(1)**3 - xznl(1)**3)*model_hse(1,idens)

     ihe_layer = -1
     ihe_entropy = -1

     !-------------------------------------------------------------------------
     ! HSE + entropy solve
     !-------------------------------------------------------------------------
     do i = 2, nx

        delx = xzn_hse(i) - xzn_hse(i-1)

        ! as the initial guess for the density, use the previous zone
        dens_zone = model_hse(i-1,idens)

        if (dens_zone > rho_he) then
           temp_zone = temp_core
           xn(:) = xn_core(:)

           isentropic = .false.

        else

           if (ihe_layer == -1) then
              ihe_layer = i
           endif

           ! determine whether we are starting the ramp up.  We will
           ! use a tanh profile, centered at (xzn_hse(ihe_layer) +
           ! FOUR*delta).  The "+ FOUR*delta" enables us to capture
           ! the leading edge of the profile.  Since rho_he is
           ! computed by considering the integral of He on the grid,
           ! shifting the profile by FOUR*delta doesn't affect the
           ! overall mass.

           test = HALF*(ONE + tanh((xzn_hse(i) - xzn_hse(ihe_layer) - FOUR*delta)/delta))

           if (test < 0.999d0) then
 
              ! small tanh ramp up regime
              xn(:) = xn_core(:) + HALF*(xn_he(:) - xn_core(:))* &
                   (ONE + tanh((xzn_hse(i) - xzn_hse(ihe_layer) - FOUR*delta)/delta))

              temp_zone = temp_core + HALF*(temp_base - temp_core)* &
                   (ONE + tanh((xzn_hse(i) - xzn_hse(ihe_layer) - FOUR*delta)/delta))
              
              isentropic = .false.

           else

              ! fully isentropic
              if (ihe_entropy == -1) then
                 ihe_entropy = i
                 temp_zone = temp_base
                 isentropic = .false.
              else
                 temp_zone = model_hse(i-1,itemp)
                 isentropic = .true.
              endif

              xn(:) = xn_he(:)

           endif

        endif

        g_zone = -Gconst*M_enclosed(i-1)/(xznl(i)*xznl(i))


        !----------------------------------------------------------------------
        ! iteration loop
        !----------------------------------------------------------------------

        ! start off the Newton loop by saying that the zone has not converged
        converged_hse = .FALSE.

        if (.not. fluff) then

           do iter = 1, MAX_ITER


              if (isentropic) then

                 p_want = model_hse(i-1,ipres) + &
                      delx*0.5_dp_t*(dens_zone + model_hse(i-1,idens))*g_zone


                 ! now we have two functions to zero:
                 !   A = p_want - p(rho,T)
                 !   B = entropy_base - s(rho,T)
                 ! We use a two dimensional Taylor expansion and find the deltas
                 ! for both density and temperature   
                 
                 eos_state%T     = temp_zone
                 eos_state%rho   = dens_zone
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
                 B = entropy_base - entropy

                 dAdT = -dpt
                 dAdrho = 0.5d0*delx*g_zone - dpd
                 dBdT = -dst
                 dBdrho = -dsd

                 dtemp = (B - (dBdrho/dAdrho)*A)/ &
                      ((dBdrho/dAdrho)*dAdT - dBdT)

                 drho = -(A + dAdT*dtemp)/dAdrho

                 dens_zone = max(0.9_dp_t*dens_zone, &
                                 min(dens_zone + drho, 1.1_dp_t*dens_zone))

                 temp_zone = max(0.9_dp_t*temp_zone, &
                                 min(temp_zone + dtemp, 1.1_dp_t*temp_zone))

                 ! check if the density falls below our minimum
                 ! cut-off -- if so, floor it
                 if (dens_zone < low_density_cutoff) then

                    dens_zone = low_density_cutoff
                    temp_zone = temp_fluff
                    converged_hse = .TRUE.
                    fluff = .TRUE.
                    exit

                 endif

                 if ( abs(drho) < TOL_HSE*dens_zone .and. &
                      abs(dtemp) < TOL_HSE*temp_zone) then
                    converged_hse = .TRUE.
                    exit
                 endif

              else
                 ! the core is isothermal, so we just need to constrain
                 ! the density and pressure to agree with the EOS and HSE

                 ! We difference HSE about the interface between the current
                 ! zone and the one just inside.
                 p_want = model_hse(i-1,ipres) + &
                      delx*0.5*(dens_zone + model_hse(i-1,idens))*g_zone

                 eos_state%T     = temp_zone
                 eos_state%rho   = dens_zone
                 eos_state%xn(:) = xn(:)
                 
                 ! (t, rho) -> (p, s)
                 call eos(eos_input_rt, eos_state)
        
                 entropy = eos_state%s
                 pres_zone = eos_state%p
                 
                 dpd = eos_state%dpdr
              
                 drho = (p_want - pres_zone)/(dpd - 0.5*delx*g_zone)
              
                 dens_zone = max(0.9*dens_zone, &
                      min(dens_zone + drho, 1.1*dens_zone))
              
                 if (abs(drho) < TOL_HSE*dens_zone) then
                    converged_hse = .TRUE.
                    exit
                 endif
              
                 if (dens_zone < low_density_cutoff) then
                 
                    icutoff = i
                    dens_zone = low_density_cutoff
                    temp_zone = temp_fluff
                    converged_hse = .TRUE.
                    fluff = .TRUE.
                    exit
                 
                 endif
              endif

              if (temp_zone < temp_fluff .and. isentropic) then
                 temp_zone = temp_fluff
                 isentropic = .false.
              endif


           enddo

           if (.NOT. converged_hse) then
           
              print *, 'Error zone', i, ' did not converge in init_1d'
              print *, dens_zone, temp_zone
              print *, p_want
              print *, drho
              call bl_error('Error: HSE non-convergence')
           
           endif

        else
           dens_zone = low_density_cutoff
           temp_zone = temp_fluff
        endif
           

        ! call the EOS one more time for this zone and then go on
        ! to the next
        eos_state%T     = temp_zone
        eos_state%rho   = dens_zone
        eos_state%xn(:) = xn(:)

        ! (t, rho) -> (p, s)    
        call eos(eos_input_rt, eos_state)

        pres_zone = eos_state%p

        ! determine the entropy that we want to constrain to, if
        ! this is the first zone of the He layer
        if (i == ihe_entropy) then
           entropy_base = entropy
        endif


        ! update the thermodynamics in this zone
        model_hse(i,idens) = dens_zone
        model_hse(i,itemp) = temp_zone
        model_hse(i,ipres) = pres_zone

        model_hse(i,ispec:ispec-1+nspec) = xn(:)

        M_enclosed(i) = M_enclosed(i-1) + &
             FOUR3RD*M_PI*(xznr(i) - xznl(i))* &
             (xznr(i)**2 +xznl(i)*xznr(i) + xznl(i)**2)*model_hse(i,idens)

        cs_hse(i) = eos_state%cs
        s_hse(i) = eos_state%s

     enddo  ! end loop over zones


     ! compute the total mass of the He layer and C/O WD
     mass_he = FOUR3RD*M_PI*(xznr(1)**3 - xznl(1)**3)* &
          model_hse(1,idens)*model_hse(1,ispec-1+ihe4)

     mass_wd = FOUR3RD*M_PI*(xznr(1)**3 - xznl(1)**3)*model_hse(1,idens)* &
          (model_hse(1,ispec-1+ic12) + model_hse(1,ispec-1+io16))

     do i = 2, icutoff
        mass_he = mass_he + &
             FOUR3RD*M_PI*(xznr(i) - xznl(i))* &
             (xznr(i)**2 +xznl(i)*xznr(i) + xznl(i)**2)*model_hse(i,idens)* &
             model_hse(i,ispec-1+ihe4)

        mass_wd = mass_wd + &
             FOUR3RD*M_PI*(xznr(i) - xznl(i))* &
             (xznr(i)**2 +xznl(i)*xznr(i) + xznl(i)**2)*model_hse(i,idens)* &
             (model_hse(i,ispec-1+ic12) + model_hse(i,ispec-1+io16))
     enddo


     if (rho_c_old < 0.0_dp_t) then
        ! not enough iterations yet -- store the old central density and
        ! mass and pick a new value
        rho_c_old = rho_c
        mass_wd_old = mass_wd

        rho_he_old = rho_he
        mass_he_old = mass_he

        rho_c = 0.5*rho_c_old
        rho_he = 0.5*rho_he_old

     else
        ! have we converged
        if ( abs(mass_wd - M_tot)/M_tot < TOL_WD_MASS .and. &
             abs(mass_he - M_He)/M_He < TOL_HE_MASS) then
           mass_converged = .true.
           exit
        endif

        ! do a secant iteration:
        ! M_tot = M(rho_c) + dM/drho |_rho_c x drho + ...        
        drho_c = (M_tot - mass_wd)/ &
             ( (mass_wd  - mass_wd_old)/(rho_c - rho_c_old) )

        rho_c_old = rho_c
        mass_wd_old = mass_wd

        rho_c = min(1.1_dp_t*rho_c_old, &
                    max((rho_c + drho_c), 0.9_dp_t*rho_c_old))


        drho_he = (M_He - mass_he)/ &
             ( (mass_he  - mass_he_old)/(rho_he - rho_he_old) )

        rho_he_old = rho_he
        mass_he_old = mass_he

        rho_he = min(1.1_dp_t*rho_he_old, &
                    max((rho_he + drho_he), 0.9_dp_t*rho_he_old))

        print *, 'current mass = ', mass_wd/solar_mass, mass_he/solar_mass

     endif     

  enddo  ! end mass constraint loop

  if (.not. mass_converged) then
     print *, 'ERROR: WD mass did not converge'
     call bl_error("ERROR: mass did not converge")
  endif

  print *, 'final masses: '
  print *, ' mass WD: ', mass_wd/solar_mass
  print *, ' mass He: ', mass_He/solar_mass
  print *, ihe_layer


  ! store the model
  write(num,'(i8)') nx
  write(mass_wd_str,'(f4.2)') mass_wd/solar_mass
  write(mass_he_str,'(f5.3)') mass_He/solar_mass
  !if (mass_He/solar_mass > 0.01) then
  !   write(mass_he_str,'(f4.2)') mass_He/solar_mass
  !else
  !   write(mass_he_str,'(f6.4)') mass_He/solar_mass
  !endif

  if (mixed_co_wd) then
     outfile = "sub_chandra.M_WD-" // trim(adjustl(mass_wd_str)) // &
          ".M_He-" // trim(adjustl(mass_he_str)) // &
          ".hse.CO." // trim(adjustl(num))
  else
     outfile = "sub_chandra.M_WD-" // trim(adjustl(mass_wd_str)) // &
          ".M_He-" // trim(adjustl(mass_he_str)) // &
          ".hse.C." // trim(adjustl(num))
  endif

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

  close (50)


  ! extra info
  if (mixed_co_wd) then
     outfile = "sub_chandra.M_WD-" // trim(adjustl(mass_wd_str)) // &
          ".M_He-" // trim(adjustl(mass_he_str)) // &
          ".extras.CO." // trim(adjustl(num))
  else
     outfile = "sub_chandra.M_WD-" // trim(adjustl(mass_wd_str)) // &
          ".M_He-" // trim(adjustl(mass_he_str)) // &
          ".extras.C." // trim(adjustl(num))
  endif

  open (unit=51, file=outfile, status="unknown")

  write (51,1001) "# npts = ", nx
  write (51,1002) "# cs"
  write (51,1002) "# entropy"

  do i = 1, nx
     write (51,1000) xzn_hse(i), cs_hse(i), s_hse(i) 
  enddo

  close (51)

  ! compute the maximum HSE error
  max_hse_error = -1.d30

  do i = 2, nx-1
     g_zone = -Gconst*M_enclosed(i-1)/xznr(i-1)**2
     dpdr = (model_hse(i,ipres) - model_hse(i-1,ipres))/delx
     rhog = HALF*(model_hse(i,idens) + model_hse(i-1,idens))*g_zone

     print *, xzn_hse(i), g_zone*xzn_hse(i)

     if (dpdr /= ZERO .and. model_hse(i+1,idens) > low_density_cutoff) then
        max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(dpdr))
     endif

  enddo

  print *, 'maximum HSE error = ', max_hse_error
  print *, ' '

end program init_1d

