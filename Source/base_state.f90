module base_state_module

  use bl_types
  use network, only: nspec
  
  implicit none

  real(dp_t), save :: base_cutoff_density_loc, prob_hi_r
  real(dp_t), save :: rho_above_cutoff, rhoh_above_cutoff
  real(dp_t), save :: spec_above_cutoff(nspec), p_above_cutoff
  real(dp_t), save :: temp_above_cutoff
  real(dp_t), save :: trac_above_cutoff(10) ! if ntrac=0 we'd have a problem if we used ntrac

  private

  public :: init_base_state, get_model_npts

contains

  subroutine init_base_state(n,model_file,s0_init,p0_init,dx)

    use bl_prof_module
    use parallel
    use bl_error_module
    use bl_constants_module
    use eos_module
    use network, only: spec_names
    use probin_module, only: base_cutoff_density, anelastic_cutoff, prob_lo, prob_hi, &
                             small_temp, small_dens, grav_const
    use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp, trac_comp, ntrac
    use geometry, only: dr, spherical, nr, dm
    use inlet_bc_module, only: set_inlet_bcs
    use fundamental_constants_module, only: Gconst

    
    integer           , intent(in   ) :: n
    character(len=256), intent(in   ) :: model_file
    real(kind=dp_t)   , intent(inout) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(inout) :: p0_init(0:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    ! local
    integer         :: i,j,r,comp
    real(kind=dp_t) :: rloc,dr_in,rmax,starting_rad,mod_dr
    real(kind=dp_t) :: d_ambient,t_ambient,p_ambient,xn_ambient(nspec)
    real(kind=dp_t) :: sumX

    ! these indices define how the initial model is stored in the 
    ! base_state array
    integer, parameter :: nvars_model = 3 + nspec
    integer, parameter :: idens_model = 1
    integer, parameter :: itemp_model = 2
    integer, parameter :: ipres_model = 3
    integer, parameter :: ispec_model = 4

    integer, parameter :: MAX_VARNAME_LENGTH=80
    integer :: npts_model, nvars_model_file

    real(kind=dp_t) :: min_dens, max_dens, min_temp, max_temp

    real(kind=dp_t), allocatable :: base_state(:,:), base_r(:)
    real(kind=dp_t), allocatable :: vars_stored(:)
    character(len=MAX_VARNAME_LENGTH), allocatable :: varnames_stored(:)
    logical :: found

    integer :: ipos
    character (len=256) :: header_line

    type(bl_prof_timer), save :: bpt
    
    real(kind=dp_t), parameter :: TINY = 1.0d-10

    real(kind=dp_t) :: mencl, g, r_l, r_r, dpdr, rhog
    real(kind=dp_t) :: max_hse_error

    call build(bpt, "init_base_state")

    do_diag = .false.

    ! open the model file and read in the header
    ! the model file is assumed to be of the follow form:
    ! # npts = 896
    ! # num of variables = 6
    ! # density
    ! # temperature
    ! # pressure
    ! # carbon-12
    ! # oxygen-16
    ! # magnesium-24
    ! 195312.5000  5437711139.  8805500.952   .4695704813E+28  0.3  0.7  0
    ! 585937.5000  5410152416.  8816689.836  0.4663923963E+28  0.3  0.7  0

    ! we read in the number of variables and their order and use this to map 
    ! them into the base_state array.  We ignore anything other than density, 
    ! temperature, pressure and composition.  

    ! Presently, we take density, temperature, and composition as the 
    ! independent variables and use them to define the thermodynamic state.

    ! composition is assumed to be in terms of mass fractions

    open(99,file=model_file)

    ! the first line has the number of points in the model
    read (99, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) npts_model

    if ( parallel_IOProcessor() ) &
         print *, npts_model, '    points found in the initial model file'

    ! now read in the number of variables
    read (99, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) nvars_model_file

    if ( parallel_IOProcessor() ) &
         print *, nvars_model_file, ' variables found in the initial model file'

    allocate (vars_stored(nvars_model_file))
    allocate (varnames_stored(nvars_model_file))

    ! now read in the names of the variables
    do i = 1, nvars_model_file
       read (99, '(a256)') header_line
       ipos = index(header_line, '#') + 1
       varnames_stored(i) = trim(adjustl(header_line(ipos:)))
    enddo

    ! allocate storage for the model data
    allocate (base_state(npts_model, nvars_model))
    allocate (base_r(npts_model))

    do i = 1, npts_model
       read(99,*) base_r(i), (vars_stored(j), j = 1, nvars_model_file)

       base_state(i,:) = ZERO

       do j = 1,nvars_model_file

          found = .false.

          if (trim(varnames_stored(j)) == "density") then
             base_state(i,idens_model) = vars_stored(j)
             found = .true.

          else if (trim(varnames_stored(j)) == "temperature") then
             base_state(i,itemp_model) = vars_stored(j)
             found = .true.

          else if (trim(varnames_stored(j)) == "pressure") then
             base_state(i,ipres_model) = vars_stored(j)
             found = .true.

          else
             do comp = 1, nspec
                if (trim(varnames_stored(j)) == spec_names(comp)) then
                   base_state(i,ispec_model-1+comp) = vars_stored(j)
                   found = .true.
                   exit
                endif
             enddo
          endif

          if (.NOT. found) then
             if ( parallel_IOProcessor() ) then
                print *, 'ERROR: variable not found: ', varnames_stored(j)
             end if
          endif

       enddo

    end do

    close(99)

    max_dens = maxval(base_state(:,idens_model))
    min_dens = minval(base_state(:,idens_model))

    max_temp = maxval(base_state(:,itemp_model))
    min_temp = minval(base_state(:,itemp_model))

    if ( parallel_IOProcessor() ) then
       print *, ' '
       print *, 'model read in:'
       print *, '    maximum density of model =                   ', max_dens
       print *, '    minimum density of model =                   ', min_dens
       print *, '    maximum temperature of model =               ', max_temp
       print *, '    minimum temperature of model =               ', min_temp
       print *, '    anelastic cutoff =                           ', anelastic_cutoff
       print *, '    low density cutoff (for mapping the model) = ', base_cutoff_density
       print *, ' '
    end if

    if (min_dens < base_cutoff_density .OR. min_dens < anelastic_cutoff) then
       if ( parallel_IOProcessor() ) then
          print *, ' '
          print *, 'WARNING: minimum model density is lower than one of the cutoff densities'
          print *, '         make sure that the cutoff densities are lower than any density'
          print *, '         of dynamical interest'
          print *, ' '
       end if
    endif

    if (min_temp < small_temp) then
       if ( parallel_IOProcessor() ) then
          print *, ' '
          print *, 'WARNING: minimum model temperature is lower than the EOS cutoff'
          print *, '         temperature, small_temp'
          print *, ' '
       endif
    endif

    if (min_dens < small_dens) then
       if ( parallel_IOProcessor() ) then
          print *, ' '
          print *, 'WARNING: minimum model density is lower than the EOS cutoff'
          print *, '         density, small_dens'
          print *, ' '
       endif
    endif


    if (anelastic_cutoff < base_cutoff_density) then
       print *, 'ERROR: anelastic cutoff should be at a higher density than the base state'
       print *, '       cutoff density.'
       call bl_error("anelastic cutoff < base_cutoff_density")
    endif

    if ( parallel_IOProcessor() ) then
        print *, ' '
        print *, ' '
    end if

    dr_in = (base_r(npts_model) - base_r(1)) / dble(npts_model-1)
    rmax = base_r(npts_model)

    if ( parallel_IOProcessor() ) then
       print *,'DR , RMAX OF MODEL     ',dr_in, rmax

       if (dr(n) .lt. dr_in) then
          mod_dr = mod(dr_in,dr(n))
       else
          mod_dr = mod(dr(n),dr_in)
       endif

       if (mod_dr .gt. TINY) then
          print *, ''
          print *, "WARNING: resolution of base state array is not an integer"
          print *, "         multiple of the initial model's resolution.     "
          print *, "         make sure this is a desired property as this    "
          print *, "         could lead to aliasing when performing the      "
          print *, "         interpolation.                                  "
          print *, " "
          print *, "modulus = ", mod_dr
       endif
    end if

    if (spherical .eq. 0) then
       starting_rad = prob_lo(dm)
    else
       starting_rad = ZERO
    endif

    if (n .eq. 1) then
       base_cutoff_density_loc = prob_hi(dm)
       prob_hi_r = base_cutoff_density_loc
    end if
    
    do r=0,nr(n)-1

       rloc = starting_rad + (dble(r) + HALF)*dr(n)

       ! here we account for r > rmax of the model.hse array, assuming
       ! that the state stays constant beyond rmax
       rloc = min(rloc, rmax)

       ! also, if we've falled below the cutoff density, just keep the
       ! model constant
       if (rloc .gt. base_cutoff_density_loc) then

          s0_init(r,rho_comp) = rho_above_cutoff
          s0_init(r,rhoh_comp) = rhoh_above_cutoff
          s0_init(r,spec_comp:spec_comp+nspec-1) = spec_above_cutoff(1:nspec)
          p0_init(r) = p_above_cutoff
          s0_init(r,temp_comp) = temp_above_cutoff
          if(ntrac .gt. 0) then
             s0_init(r,trac_comp:trac_comp+ntrac-1) = trac_above_cutoff(1:ntrac)
          end if

       else

          d_ambient = interpolate(rloc, npts_model, base_r, base_state(:,idens_model))
          t_ambient = interpolate(rloc, npts_model, base_r, base_state(:,itemp_model))
          p_ambient = interpolate(rloc, npts_model, base_r, base_state(:,ipres_model))

          sumX = ZERO
          do comp = 1, nspec
             xn_ambient(comp) = max(ZERO,min(ONE, &
                  interpolate(rloc, npts_model, base_r, base_state(:,ispec_model-1+comp))))
             sumX = sumX + xn_ambient(comp)
          enddo
          xn_ambient = xn_ambient/sumX

          ! use the EOS to make the state consistent
          temp_eos(1) = t_ambient
          den_eos(1)  = d_ambient
          p_eos(1)    = p_ambient
          xn_eos(1,:) = xn_ambient(:)

          ! (rho,T) --> p,h
          call eos(eos_input_rt, den_eos, temp_eos, &
               npts, nspec, &
               xn_eos, &
               p_eos, h_eos, e_eos, &
               cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
               dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
               dpdX_eos, dhdX_eos, &
               gam1_eos, cs_eos, s_eos, &
               dsdt_eos, dsdr_eos, &
               do_diag)

          s0_init(r, rho_comp ) = d_ambient
          s0_init(r,rhoh_comp ) = d_ambient * h_eos(1)
          s0_init(r,spec_comp:spec_comp+nspec-1) = d_ambient * xn_ambient(1:nspec)
          p0_init(r) = p_eos(1)

          s0_init(r,temp_comp) = t_ambient

          if (ntrac .gt. 0) then
             s0_init(r,trac_comp:trac_comp+ntrac-1) = ZERO
          end if

          ! keep track of the height where we drop below the cutoff density
          if (s0_init(r,rho_comp) .le. base_cutoff_density .and. &
               base_cutoff_density_loc .eq. prob_hi_r .and. n .eq. 1 ) then

             if ( parallel_IOProcessor() ) then
                print *,'SETTING R_CUTOFF TO ',r
             end if

             base_cutoff_density_loc = rloc

             rho_above_cutoff = s0_init(r,rho_comp)
             rhoh_above_cutoff = s0_init(r,rhoh_comp)
             spec_above_cutoff(1:nspec) = s0_init(r,spec_comp:spec_comp+nspec-1)
             p_above_cutoff = p0_init(r)
             temp_above_cutoff = s0_init(r,temp_comp)
             if (ntrac .gt. 0) then
                trac_above_cutoff(1:ntrac) = s0_init(r,trac_comp:trac_comp+ntrac-1)
             end if

          end if

       end if

    end do

    ! check whether we are in HSE
    
    if (spherical .eq. 1) then
       mencl = four3rd*m_pi*dr(n)**3*s0_init(0,rho_comp)
    endif

    max_hse_error = -1.d30

    do r=1,nr(n)-1
       
       rloc = starting_rad + (dble(r) + HALF)*dr(n)
       rloc = min(rloc, rmax)

       if (rloc .lt. base_cutoff_density_loc) then

          r_r = starting_rad + dble(r+1)*dr(n)
          r_l = starting_rad + dble(r)*dr(n)

          if (spherical .eq. 1) then
             g = -Gconst*mencl/r_l**2
             mencl = mencl &
                  + four3rd*m_pi*dr(n)*(r_l**2+r_l*r_r+r_r**2)*s0_init(r,rho_comp)
          else
             g = grav_const
          endif

          dpdr = (p0_init(r) - p0_init(r-1))/dr(n)
          rhog = HALF*(s0_init(r,rho_comp) + s0_init(r-1,rho_comp))*g

          !write(*,1000) r, dpdr, rhog, abs(dpdr - rhog)/abs(dpdr), s0_init(r,rho_comp)
1000      format(1x,6(g20.10))

          max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(dpdr))

       end if

    enddo

    if ( parallel_IOProcessor() ) then
       print *, " " 
       print *, "Level: ", n
       print *, "Maximum HSE Error = ", max_hse_error
       print *, "   (after putting initial model into base state arrays, and"
       print *, "    for density < base_cutoff_density)"
       print *, " "
    endif


    ! initialize any inlet BC parameters
    call set_inlet_bcs()

    deallocate(vars_stored,varnames_stored)
    deallocate(base_state,base_r)

    call destroy(bpt)

  end subroutine init_base_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function interpolate(r, npts, model_r, model_var)

    use bl_constants_module

    ! given the array of model coordinates (model_r), and variable (model_var),
    ! find the value of model_var at point r using linear interpolation.
    ! Eventually, we can do something fancier here.

    real(kind=dp_t) :: interpolate
    real(kind=dp_t), intent(in) :: r
    integer :: npts
    real(kind=dp_t), dimension(npts) :: model_r, model_var

    real(kind=dp_t) :: slope
    real(kind=dp_t) :: minvar, maxvar

    integer :: i, id

    ! find the location in the coordinate array where we want to interpolate
    do i = 1, npts
       if (model_r(i) >= r) exit
    enddo
    if(i .gt. 1 .and. i .lt. npts+1) then
       if(abs(r-model_r(i-1)) .lt. abs(r-model_r(i))) then
          i = i-1
       end if
    end if
    if (i == npts+1) then
       i = npts
    end if

    id = i

    if (id == 1) then

       slope = (model_var(id+1) - model_var(id))/(model_r(id+1) - model_r(id))
       interpolate = slope*(r - model_r(id)) + model_var(id)

       ! safety check to make sure interpolate lies within the bounding points
       minvar = min(model_var(id+1),model_var(id))
       maxvar = max(model_var(id+1),model_var(id))
       interpolate = max(interpolate,minvar)
       interpolate = min(interpolate,maxvar)

    else if (id == npts) then

       slope = (model_var(id) - model_var(id-1))/(model_r(id) - model_r(id-1))
       interpolate = slope*(r - model_r(id)) + model_var(id)

       ! safety check to make sure interpolate lies within the bounding points
       minvar = min(model_var(id),model_var(id-1))
       maxvar = max(model_var(id),model_var(id-1))
       interpolate = max(interpolate,minvar)
       interpolate = min(interpolate,maxvar)

    else

       if (r .ge. model_r(id)) then

          ! we should not wind up in here

          slope = (model_var(id+1) - model_var(id))/(model_r(id+1) - model_r(id))
          interpolate = slope*(r - model_r(id)) + model_var(id)
          
          ! safety check to make sure interpolate lies within the bounding points
          minvar = min(model_var(id+1),model_var(id))
          maxvar = max(model_var(id+1),model_var(id))
          interpolate = max(interpolate,minvar)
          interpolate = min(interpolate,maxvar)
          
       else

          slope = (model_var(id) - model_var(id-1))/(model_r(id) - model_r(id-1))
          interpolate = slope*(r - model_r(id)) + model_var(id)
          
          ! safety check to make sure interpolate lies within the bounding points
          minvar = min(model_var(id),model_var(id-1))
          maxvar = max(model_var(id),model_var(id-1))
          interpolate = max(interpolate,minvar)
          interpolate = min(interpolate,maxvar)
          
       end if

    end if

    return

  end function interpolate


  function get_model_npts(model_file)

    ! look in the model file and return the number of points
    real(kind=dp_t) :: get_model_npts

    character(len=256), intent(in   ) :: model_file

    character (len=256) :: header_line
    integer :: ipos

    open(99,file=model_file)

    ! the first line has the number of points in the model
    read (99, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) get_model_npts

    close(99)

    return

  end function get_model_npts

end module base_state_module
