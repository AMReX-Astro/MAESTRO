module base_state_module

  use bl_types

  implicit none

  private

  public :: init_base_state

contains

  subroutine init_base_state(n,model_file,s0_init,p0_init,dx)

    use bc_module
    use setbc_module
    use bl_prof_module
    use multifab_module
    use define_bc_module
    use bl_constants_module
    use eos_module
    use network, only: nspec, spec_names
    use probin_module, only: base_cutoff_density, anelastic_cutoff, prob_lo
    use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp, trac_comp, ntrac
    use geometry, only: dr, nr, spherical

    integer           , intent(in   ) :: n
    character(len=256), intent(in   ) :: model_file
    real(kind=dp_t)   , intent(inout) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(inout) :: p0_init(0:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    real(kind=dp_t)   :: max_speed

    ! local
    integer         :: i,j,j_cutoff,comp
    real(kind=dp_t) :: r,dr_in,rmax,starting_rad
    real(kind=dp_t) :: d_ambient,t_ambient,p_ambient,xn_ambient(nspec)
    real(kind=dp_t) :: sum

    ! these indices define how the initial model is stored in the 
    ! base_state array
    integer, parameter :: nvars_model = 3 + nspec
    integer, parameter :: idens_model = 1
    integer, parameter :: itemp_model = 2
    integer, parameter :: ipres_model = 3
    integer, parameter :: ispec_model = 4

    integer, parameter :: MAX_VARNAME_LENGTH=80
    integer :: npts_model, nvars_model_file

    real(kind=dp_t) :: min_dens, max_dens

    real(kind=dp_t), allocatable :: base_state(:,:), base_r(:)
    real(kind=dp_t), allocatable :: vars_stored(:)
    character(len=MAX_VARNAME_LENGTH), allocatable :: varnames_stored(:)
    logical :: found

    integer :: ipos
    character (len=256) :: header_line

    type(bl_prof_timer), save :: bpt

    call build(bpt, "init_base_state")

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

       do j = 1, nvars_model_file

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
             print *, 'ERROR: variable not found: ', varnames_stored(j)
          endif

       enddo

    end do

    close(99)

    max_dens = maxval(base_state(:,idens_model))
    min_dens = minval(base_state(:,idens_model))

    if ( parallel_IOProcessor() ) then
       print *, ' '
       print *, 'model read in'
       print *, ' *** maximum density of model =                   ', max_dens
       print *, ' *** minimum density of model =                   ', min_dens
       print *, ' *** anelastic cutoff =                           ', anelastic_cutoff
       print *, ' *** low density cutoff (for mapping the model) = ', base_cutoff_density
       print *, ' '
    end if

    if (min_dens < base_cutoff_density .OR. min_dens < anelastic_cutoff) then
       if ( parallel_IOProcessor() ) then
          print *, 'WARNING: minimum model density is lower than one of the cutoff densities'
          print *, '         make sure that the cutoff densities are lower than any density'
          print *, '         of dynamical interest'
       end if
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
       print *,'DR , RMAX OF BASE ARRAY',dr(n), dble(nr(n)) * dr(n)
    end if

    starting_rad = prob_lo(size(dx))

    j_cutoff = nr(n)
    max_speed = 0.d0
    do j = 0,nr(n)-1

       if (j .ge. j_cutoff) then

          s0_init(j, rho_comp ) = s0_init(j_cutoff, rho_comp )
          s0_init(j,rhoh_comp ) = s0_init(j_cutoff,rhoh_comp )
          s0_init(j,spec_comp:spec_comp+nspec-1) = s0_init(j_cutoff,spec_comp:spec_comp+nspec-1)
          p0_init(j)            = p0_init(j_cutoff)
          s0_init(j,temp_comp)  = s0_init(j_cutoff,temp_comp)

       else

          ! compute the coordinate height at this level
          ! NOTE: we are assuming that the basestate is in the y-direction
          ! and that ymin = 0.0
          r = starting_rad + (dble(j) + HALF)*dr(n)

          ! here we account for r > rmax of the model.hse array, assuming
          ! that the state stays constant beyond rmax
          r = min(r, rmax)

          d_ambient = interpolate(r, npts_model, base_r, base_state(:,idens_model))
          t_ambient = interpolate(r, npts_model, base_r, base_state(:,itemp_model))
          p_ambient = interpolate(r, npts_model, base_r, base_state(:,ipres_model))

          sum = ZERO
          do comp = 1, nspec
             xn_ambient(comp) = max(ZERO,min(ONE, &
                  interpolate(r, npts_model, base_r, base_state(:,ispec_model-1+comp))))
             sum = sum + xn_ambient(comp)
          enddo
          xn_ambient(:) = xn_ambient(:)/sum

          ! use the EOS to make the state consistent
          temp_eos  = t_ambient
          den_eos   = d_ambient
          p_eos     = p_ambient
          xn_eos(:) = xn_ambient(:)

          ! (rho,T) --> p,h
          call eos(eos_input_rt, den_eos, temp_eos, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, &
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   .false.)
          max_speed = max(max_speed,cs_eos)
!         print *,'CS ',j,cs_eos(1)

          s0_init(j, rho_comp ) = d_ambient
          s0_init(j,rhoh_comp ) = d_ambient * h_eos
          s0_init(j,spec_comp:spec_comp+nspec-1) = d_ambient * xn_ambient(1:nspec)
          p0_init(j)    = p_eos

          s0_init(j,temp_comp) = t_ambient

          s0_init(j,trac_comp:trac_comp+ntrac-1) = ZERO

          ! keep track of the height where we drop below the cutoff density
          if (s0_init(j,rho_comp) .lt. base_cutoff_density .and. j_cutoff .eq. nr(n)) then
             if ( parallel_IOProcessor() ) print *,'SETTING J_CUTOFF TO ',j
             j_cutoff = j
          end if

       end if

    end do
    print *,'MAX CS ',max_speed

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

    real(kind=dp_t) :: slope, xi, dr_model
    real(kind=dp_t) :: minvar, maxvar

    integer :: i, id
    logical :: use_quadratic_interpolation

    use_quadratic_interpolation = .false.

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

       if (use_quadratic_interpolation) then

          ! do a quadratic interpolation
          dr_model = model_r(id+1) - model_r(id)
          xi = r - model_r(id)
          interpolate = &
               (model_var(id+1) - 2*model_var(id) + model_var(id-1))*xi**2/(2*dr_model**2) + &
               (model_var(id+1) - model_var(id-1))*xi/(2*dr_model) + &
               (-model_var(id+1) + 26*model_var(id) - model_var(id-1))/24.0_dp_t
          
          minvar = min(model_var(id), min(model_var(id-1), model_var(id+1)))
          maxvar = max(model_var(id), max(model_var(id-1), model_var(id+1)))
          
          ! if we overshoot or undershoot, use linear interpolation
          if (interpolate > maxvar .OR. interpolate < minvar) then
             
             slope = (model_var(id+1) - model_var(id-1))/(model_r(id+1) - model_r(id-1))
             interpolate = slope*(r - model_r(id)) + model_var(id)      
             interpolate = max(interpolate,minvar)
             interpolate = min(interpolate,maxvar)
             
          endif
          
       else
          
          if (r .ge. model_r(id)) then

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

    endif

    return

  end function interpolate

end module base_state_module
