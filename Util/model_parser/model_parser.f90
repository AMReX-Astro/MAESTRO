module model_parser_module

  ! read in an initial model and return arrays with the model data.
  ! take care to match the species available in the model file to
  ! those defined by the network

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
  !
  ! we read in the number of variables and their order and use this to
  ! map them into the model_state array.  We ignore anything other than
  ! density, temperature, pressure and composition.
  !
  ! composition is assumed to be in terms of mass fractions     

  use parallel, only: parallel_IOProcessor
  use network
  use bl_types
  use simple_log_module

  implicit none

  ! integer keys for indexing the model_state array
  integer, parameter :: nvars_model = 3 + nspec
  integer, parameter :: idens_model = 1
  integer, parameter :: itemp_model = 2
  integer, parameter :: ipres_model = 3
  integer, parameter :: ispec_model = 4

  ! number of points in the model file
  integer, save :: npts_model

  ! arrays for storing the model data
  real (kind=dp_t), allocatable, save :: model_state(:,:)
  real (kind=dp_t), allocatable, save :: model_r(:)

  ! model_initialized will be .true. once the model is read in and the
  ! model data arrays are initialized and filled
  logical, save :: model_initialized = .false.

  integer, parameter :: MAX_VARNAME_LENGTH=80

  public :: read_model_file, model_parser_finalize, interpolate

contains

  subroutine read_model_file(model_file, grav_const_in, base_cutoff_density_in)

    use bl_constants_module
    use bl_error_module
    use probin_module, only: use_binary_input  

    character(len=*), intent(in   ) :: model_file

    real(kind=dp_t), optional :: grav_const_in, base_cutoff_density_in

    ! local variables
    integer :: nvars_model_file
    integer :: ierr
    integer :: varname_length_binary
    integer :: i, j, comp

    real(kind=dp_t) :: grav_const, base_cutoff_density

    double precision :: rhog, dpdr, max_hse_error

    real(kind=dp_t), allocatable :: vars_stored(:)
    real(kind=dp_t), allocatable :: vars_stored_binary(:,:)
    character(len=MAX_VARNAME_LENGTH), allocatable :: varnames_stored(:)
    logical :: found_model, found_dens, found_temp, found_pres
    logical :: found_spec(nspec)
    integer :: ipos
    character (len=256) :: header_line

    if (present(grav_const_in)) then
       grav_const = grav_const_in
    else
       grav_const = 0.0_dp_t
    endif

    if (present(base_cutoff_density_in)) then
       base_cutoff_density = base_cutoff_density_in
    else
       base_cutoff_density = 0.0_dp_t
    endif
       
    
    ! open the model file
    if (use_binary_input) then
      open(99, file=trim(model_file),status='old',iostat=ierr,form="unformatted",access="stream")
    else
      open(99,file=trim(model_file),status='old',iostat=ierr)
    endif
    if (ierr .ne. 0) then
       print *,'Couldnt open model_file: ',model_file
       call bl_error('Aborting now -- please supply model_file')
    end if
    
    if (use_binary_input) then
        read(99) npts_model
	read(99) nvars_model_file
    else
      ! the first line has the number of points in the model
      read (99, '(a256)') header_line
      ipos = index(header_line, '=') + 1
      read (header_line(ipos:),*) npts_model

      ! now read in the number of variables
      read (99, '(a256)') header_line
      ipos = index(header_line, '=') + 1
      read (header_line(ipos:),*) nvars_model_file
    endif
    
    if (use_binary_input) then
      allocate (vars_stored_binary(npts_model,nvars_model_file))
      allocate (varnames_stored(nvars_model_file))
      !allocate storage for model_r 
      allocate (model_r(npts_model))
      read(99) model_r
      read(99) vars_stored_binary
      read(99) varname_length_binary
      
      if (varname_length_binary .ne. MAX_VARNAME_LENGTH) then 
	call bl_error ('string size in binary files does not match expected size')
      endif
      read(99) varnames_stored
    else
      allocate (vars_stored(nvars_model_file))
      allocate (varnames_stored(nvars_model_file))
      !allocate storage for model_r
      allocate (model_r(npts_model))
      ! now read in the names of the variables
      do i = 1, nvars_model_file
	read (99, '(a256)') header_line
	ipos = index(header_line, '#') + 1
	varnames_stored(i) = trim(adjustl(header_line(ipos:)))
      enddo
    endif
    
    
    ! allocate storage for the remaining model data
    allocate (model_state(npts_model, nvars_model))


887 format(78('-'))
889 format(a60)

    if ( parallel_IOProcessor()) then
       write (*,889) ' '
       write (*,887)
       write (*,*)   'reading initial model'
       write (*,*)   npts_model, 'points found in the initial model file'
       write (*,*)   nvars_model_file, ' variables found in the initial model file'
    endif
    
    if (use_binary_input) then
      ! ensure that the data is put to the right position in model_state and 
      ! check wether all variables that MASETRO cares about are there
      do j = 1,nvars_model_file
        found_dens = .false.
	found_temp = .false.
	found_pres = .false.
	found_spec(:) = .false.
	
	if (varnames_stored(j) == "density") then
	      model_state(:,idens_model) = vars_stored_binary(:,j)
	      found_model = .true.
	      found_dens  = .true.

	    else if (varnames_stored(j) == "temperature") then
	      model_state(:,itemp_model) = vars_stored_binary(:,j)
	      found_model = .true.
	      found_temp  = .true.

	    else if (varnames_stored(j) == "pressure") then
	      model_state(:,ipres_model) = vars_stored_binary(:,j)
	      found_model = .true.
	      found_pres  = .true.

	    else
	      do comp = 1, nspec
		  if (varnames_stored(j) == spec_names(comp) .or. &
		      varnames_stored(j) == short_spec_names(comp)) then
		    model_state(:,ispec_model-1+comp) = vars_stored_binary(:,j)
		    found_model = .true.
		    found_spec(comp) = .true.
		    exit
		  endif
	      enddo
	    endif

	    ! is the current variable from the model file one that we
	    ! care about?
	    if (.NOT. found_model .and. i == 1) then
	      if ( parallel_IOProcessor() ) then
		  call log('WARNING: variable not found: ' // &
		      trim(varnames_stored(j)))
	      end if
	    endif
      enddo
      
      ! were all the variables we care about provided?
	if (i == 1) then
	    if (.not. found_dens) then
	      if ( parallel_IOProcessor() ) then
		  call log('WARNING: density not provided in inputs file')
	      end if
	    endif

	    if (.not. found_temp) then
	      if ( parallel_IOProcessor() ) then
		  call log('WARNING: temperature not provided in inputs file')
	      end if
	    endif

	    if (.not. found_pres) then
	      if ( parallel_IOProcessor() ) then
		  call log('WARNING: pressure not provided in inputs file')
	      end if
	    endif

	    do comp = 1, nspec
	      if (.not. found_spec(comp)) then
		  if ( parallel_IOProcessor() ) then
		    call log('WARNING: ' // trim(spec_names(comp)) // &
			  ' not provided in inputs file')
		  end if
	      endif
	    enddo
	endif
      
    else ! Use a textfile as input
      ! start reading in the data
      do i = 1, npts_model
	read(99,*) model_r(i), (vars_stored(j), j = 1, nvars_model_file)
	
	model_state(i,:) = ZERO

	! make sure that each of the variables that MAESTRO cares about
	! are found
	found_dens = .false.
	found_temp = .false.
	found_pres = .false.
	found_spec(:) = .false.

	do j = 1,nvars_model_file

	    ! keep track of whether the current variable from the model
	    ! file is one that MAESTRO cares about
	    found_model = .false.


	    if (varnames_stored(j) == "density") then
	      model_state(i,idens_model) = vars_stored(j)
	      found_model = .true.
	      found_dens  = .true.

	    else if (varnames_stored(j) == "temperature") then
	      model_state(i,itemp_model) = vars_stored(j)
	      found_model = .true.
	      found_temp  = .true.

	    else if (varnames_stored(j) == "pressure") then
	      model_state(i,ipres_model) = vars_stored(j)
	      found_model = .true.
	      found_pres  = .true.

	    else
	      do comp = 1, nspec
		  if (varnames_stored(j) == spec_names(comp) .or. &
		      varnames_stored(j) == short_spec_names(comp)) then
		    model_state(i,ispec_model-1+comp) = vars_stored(j)
		    found_model = .true.
		    found_spec(comp) = .true.
		    exit
		  endif
	      enddo
	    endif

	    ! is the current variable from the model file one that we
	    ! care about?
	    if (.NOT. found_model .and. i == 1) then
	      if ( parallel_IOProcessor() ) then
		  call log('WARNING: variable not found: ' // &
		      trim(varnames_stored(j)))
	      end if
	    endif

	enddo   ! end loop over nvars_model_file

	! were all the variables we care about provided?
	if (i == 1) then
	    if (.not. found_dens) then
	      if ( parallel_IOProcessor() ) then
		  call log('WARNING: density not provided in inputs file')
	      end if
	    endif

	    if (.not. found_temp) then
	      if ( parallel_IOProcessor() ) then
		  call log('WARNING: temperature not provided in inputs file')
	      end if
	    endif

	    if (.not. found_pres) then
	      if ( parallel_IOProcessor() ) then
		  call log('WARNING: pressure not provided in inputs file')
	      end if
	    endif

	    do comp = 1, nspec
	      if (.not. found_spec(comp)) then
		  if ( parallel_IOProcessor() ) then
		    call log('WARNING: ' // trim(spec_names(comp)) // &
			  ' not provided in inputs file')
		  end if
	      endif
	    enddo
	endif

      end do   ! end loop over npts_model
    endif
    
    model_initialized = .true.

    !max_hse_error = -1.e33
    !do i = 2, npts_model
    !   dpdr = (model_state(i,ipres_model) - model_state(i-1,ipres_model))/ &
    !        (model_r(i) - model_r(i-1))
    !   rhog = HALF*(model_state(i,idens_model) + model_state(i-1,idens_model))*grav_const
    !
    !   if (model_state(i,idens_model) > base_cutoff_density) then
    !      max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(rhog))
    !   endif
    !enddo

    
    if ( parallel_IOProcessor() ) then
       call log("initial model read...")
       call log(" ")
    endif
    
    close(99)
    
    if (use_binary_input) then
      deallocate(vars_stored_binary)
    else
      deallocate(vars_stored)
    endif
    deallocate(varnames_stored)

  end subroutine read_model_file


  function get_model_npts(model_file)

    integer :: get_model_npts
  
    ! look in the model file and return the number of points
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

  function interpolate(r, ivar, interpolate_top)

    ! use the module's array of model coordinates (model_r), and
    ! variables (model_state), to find the value of model_var at point
    ! r using linear interpolation.
    !
    ! Eventually, this should all be made a class.

    real(kind=dp_t) :: interpolate
    real(kind=dp_t), intent(in) :: r
    integer, intent(in) :: ivar
    logical, intent(in), optional :: interpolate_top

    real(kind=dp_t) :: slope
    real(kind=dp_t) :: minvar, maxvar

    logical :: interp_top

    integer :: i, id

    if (present(interpolate_top)) then
       interp_top = interpolate_top
    else
       interp_top = .false.
    endif

    ! find the location in the coordinate array where we want to interpolate
    do i = 1, npts_model
       if (model_r(i) >= r) exit
    enddo
    if (i > 1 .and. i < npts_model+1) then
       if(abs(r-model_r(i-1)) < abs(r-model_r(i))) then
          i = i-1
       end if
    end if
    if (i == npts_model+1) then
       i = npts_model
    end if

    id = i

    if (id == 1) then

       slope = (model_state(id+1,ivar) - model_state(id,ivar))/(model_r(id+1) - model_r(id))
       interpolate = slope*(r - model_r(id)) + model_state(id,ivar)

       ! safety check to make sure interpolate lies within the bounding points
       minvar = min(model_state(id+1,ivar), model_state(id,ivar))
       maxvar = max(model_state(id+1,ivar), model_state(id,ivar))
       interpolate = max(interpolate,minvar)
       interpolate = min(interpolate,maxvar)

    else if (id == npts_model) then

       slope = (model_state(id,ivar) - model_state(id-1,ivar))/(model_r(id) - model_r(id-1))
       interpolate = slope*(r - model_r(id)) + model_state(id,ivar)

       ! safety check to make sure interpolate lies within the bounding points
       if (.not. interp_top) then
          minvar = min(model_state(id,ivar),model_state(id-1,ivar))
          maxvar = max(model_state(id,ivar),model_state(id-1,ivar))
          interpolate = max(interpolate,minvar)
          interpolate = min(interpolate,maxvar)
       endif

    else

       if (r >= model_r(id)) then

          ! we should not wind up in here

          slope = (model_state(id+1,ivar) - model_state(id,ivar))/(model_r(id+1) - model_r(id))
          interpolate = slope*(r - model_r(id)) + model_state(id,ivar)
          
          ! safety check to make sure interpolate lies within the bounding points
          minvar = min(model_state(id+1,ivar),model_state(id,ivar))
          maxvar = max(model_state(id+1,ivar),model_state(id,ivar))
          interpolate = max(interpolate,minvar)
          interpolate = min(interpolate,maxvar)
          
       else

          slope = (model_state(id,ivar) - model_state(id-1,ivar))/(model_r(id) - model_r(id-1))
          interpolate = slope*(r - model_r(id)) + model_state(id,ivar)
          
          ! safety check to make sure interpolate lies within the bounding points
          minvar = min(model_state(id,ivar),model_state(id-1,ivar))
          maxvar = max(model_state(id,ivar),model_state(id-1,ivar))
          interpolate = max(interpolate,minvar)
          interpolate = min(interpolate,maxvar)
          
       end if

    end if

    return

  end function interpolate

  subroutine model_parser_finalize()
    if (model_initialized) then
       deallocate (model_state)
       deallocate (model_r)
    endif
  end subroutine model_parser_finalize

end module model_parser_module

