! print out the radius range of cells above a user-specified temperature
! along with the radius of the CO/He interface as determined by C mass fraction
!
! Author:        Adam Jacobs (adam.jacobs@stonybrook.edu)
! Creation Date: Feb. 5th, 2013

! Usage: $ fhotspots -t <T_hot>  <pltfile> [<pltfile2> ...]
! Flags:
!     -t|--temp <T_hot> : specifies the hot temperature cutoff as
!                         the real number <T_hot> in Kelvin 

! TODO:
!   1) Fill out usage section above
!   2) Turn "sections" into functions
!   3) Check that t_hot is reasonable

program fhotspots
  !============================================================
  ! fhotspots - data section                                  |
  !============================================================
  ! modules
  use f2kcli
  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module

  implicit none

  ! argument variables
  character(len=256) :: pltfile, outputfile, t_hot_char
  integer :: idir, nvars_in
  logical :: idir_specified, do_favre
  character(len=80), allocatable :: vars_in(:)
  logical :: quiet
  real(kind=dp_t) :: t_hot
  
  ! auxilliary variables for arguments
  logical :: do_specific_vars
  integer, allocatable :: icomps_in(:)
  logical, allocatable :: skip(:)
  character :: dir

  ! fk2cli variables
  integer :: narg, farg
  character(len=256) :: fname

  ! local constants
  real(kind=dp_t), parameter :: XC_TOL = 0.5

  ! local variables
  type(plotfile) :: pf
  integer :: unit, uout
  integer :: dim, irho_comp, irho_comp_pass, itemp_comp, ixc_comp
  integer :: i, j
  real(kind=dp_t) :: dx(MAX_SPACEDIM), r_lo, r_hi, r_interface
  real(kind=dp_t) :: cur_temp, cur_r, cur_xc
  integer, dimension(MAX_SPACEDIM) :: flo, fhi
  logical, allocatable :: imask(:,:,:)
  real(kind=dp_t), allocatable :: avg(:,:), rms(:,:), cnt(:)
  real(kind=dp_t), pointer :: p(:,:,:,:)

  ! variables to be passed into averaging routines
  ! nvars_pass contains the total number of valid variables to be analyzed
  ! icomps_pass contains the indices (within the pf) of the valid variables
  !     being analyzed
  integer :: nvars_pass
  integer, allocatable :: icomps_pass(:)

  integer :: r1, max_points, ii, jj, kk, index, index_lo, index_hi
  real(kind=dp_t) :: weight

  ! make formatting flexible
  character(len=128) :: column_format, data_format, columns

  unit = unit_new()

  ! defaults
  pltfile = ''
  outputfile = 'stdout'
  idir_specified = .false.
  idir = 1
  dir = "x"
  do_favre = .false.
  do_specific_vars = .false.
  irho_comp = -10
  itemp_comp = -10
  ixc_comp = -10
  quiet = .false.

  !============================================================
  ! fhotspots - execution section                             |
  !============================================================
  !------------------------------------------------------------
  ! parse command line arguments                              |
  !------------------------------------------------------------
  narg = command_argument_count()

  farg = 1
  do while (farg <= narg)
     call get_command_argument(farg, value = fname)

     select case(fname)

     case('-t', '--temp')
        farg = farg + 1
        call get_command_argument(farg, value = t_hot_char)
        read(t_hot_char, *) t_hot

     case('-p', '--pltfile')
        farg = farg + 1
        call get_command_argument(farg, value = pltfile)

     case default
        exit

     end select
     farg = farg + 1
  enddo

  !------------------------------------------------------------
  ! sanity checks on specified options                        |
  !------------------------------------------------------------
  !TODO: check that t_hot is reasonable
  if (pltfile == '') then
     call print_usage()
     stop
  endif

  ! build the plotfile datatype
  call build(pf, pltfile, unit)

  dim = pf%dim

  ! make sure we don't specify z-direction for a 2d pltfile
  if (idir > dim) call bl_error("error: idir too large for pltfile")

  ! if the idir was not specified with the -d flag, then use the direction
  ! perpendicular to gravity
  if (.not. idir_specified) then
     if (dim == 2) idir = 2
     if (dim == 3) idir = 3
  endif

  if (idir == 2) dir = "y"
  if (idir == 3) dir = "z"

  ! make sure there weren't too many variables specified
  if (do_specific_vars) then
     if (nvars_in > pf%nvars) &
          call bl_error("error: -v option set with too many vars for pltfile")
  endif

  ! grab the density index in the plotfile
  ! if it is not in the plotfile and we are doing favre averaging then we abort
  do i = 1, pf%nvars
     if (var_name(pf, i) == "density") then
        irho_comp = i
        exit
     endif
  enddo

  if (do_favre .and. irho_comp < 0) &
       call bl_error("-f option set but density is not in pltfile")

  ! grab the temperature index in the plotfile
  ! if it is not in the plotfile then we abort
  do i = 1, pf%nvars
     if (var_name(pf, i) == "tfromh") then
        itemp_comp = i
        exit
     endif
  enddo

  if (itemp_comp < 0) call bl_error("no temperature component in plotfile!")

  ! grab the Carbon mass fraction index in the plotfile
  ! if it is not in the plotfile then we abort
  do i = 1, pf%nvars
     if (var_name(pf, i) == "X(C12)") then
        ixc_comp = i
        exit
     endif
  enddo

  if (ixc_comp < 0) call bl_error("no X(C12) component in plotfile!")





  !------------------------------------------------------------
  ! Allocate storage for the averaging                        |
  !------------------------------------------------------------
  ! grab the finest level's dx
  dx = plotfile_get_dx(pf, pf%flevel)

  ! get the index bounds for the finest level
  flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
  fhi = upb(plotfile_get_pd_box(pf, pf%flevel))

  ! imask will be set to false if we've already output the data.
  ! Note, imask is defined in terms of the finest level.  As we loop 
  ! over levels, we will compare to the finest level index space to
  ! determine if we've already output here
  if (dim == 2) then
     allocate(imask(flo(1):fhi(1), &
                    flo(2):fhi(2), &
                    1              ))
  else if (dim == 3) then
     allocate(imask(flo(1):fhi(1), &
                    flo(2):fhi(2), &
                    flo(3):fhi(3) ))
  endif

  ! this is the maximum number of points along the specified direction if
  ! the entire domain were at the finest level resolution
  max_points = fhi(idir) - flo(idir) + 1

  ! allocate storage for the data
  allocate(avg(0:max_points-1,nvars_pass), &
           rms(0:max_points-1,nvars_pass), &
           cnt(0:max_points-1))

  avg = ZERO
  rms = ZERO

  !------------------------------------------------------------
  ! calculate hotspot radius ranges                           |
  !------------------------------------------------------------
  ! dump some info for the user
  if (.not. quiet) then
     print *, 'calculating...'
     print *, 'pltfile = "', trim(pltfile), '"'
  endif

  ! r1 is the refinement factor between the current level and the FINEST level
  r1 = 1

  imask = .true.
  cnt = 0
  r_hi = 0.0
  r_lo = huge(r_lo)
  r_interface = huge(r_interface)


  ! loop over the data starting at the finest grid, and if we havn't already
  ! checked the data in that grid location (according to imask), check it.
  do i = pf%flevel, 1, -1
     ! grab current level's dx
     dx = plotfile_get_dx(pf, i)

     do j = 1, nboxes(pf, i)

        ! bind to a comp vector
        call fab_bind(pf, i, j)
        !call fab_bind_comp_vec(pf, i, j, icomps_pass)

        ! get the data
        p => dataptr(pf, i, j)

        ! loop over all the zones in the current box
        ! Here, we convert the cell-centered indices at the current level
        ! into the corresponding RANGE on the finest level, and test if 
        ! we've stored data in any of those locations.  If we havn't then
        ! we store this level's data and mark that range as filled.
        do kk = lbound(p,dim=3), ubound(p,dim=3)

           do jj = lbound(p,dim=2), ubound(p,dim=2)

              do ii = lbound(p,dim=1), ubound(p,dim=1)

                 ! calculate radii of interest, if neeed
                 if (dim == 2) then
                    if (any(imask(ii*r1:(ii+1)*r1-1, &
                                  jj*r1:(jj+1)*r1-1, &
                                  1))) then
                       ! get temp and radius
                       cur_temp = p(ii,jj,kk,itemp_comp)
                       cur_r = sqrt( (ii*dx(1))**2 + (jj*dx(2))**2 )

                       ! if hot, check radius 
                       if(cur_temp > t_hot) then
                         if(cur_r < r_lo) r_lo = cur_r
                         if(cur_r > r_hi) r_hi = cur_r
                       endif

                       ! get current Carbon mass fraction
                       cur_xc =  p(ii,jj,kk,ixc_comp)

                       ! record minimum radius of CO/He interface
                       if(cur_xc < XC_TOL) then
                         if(cur_r < r_interface) r_interface = cur_r
                       endif

                       imask(ii*r1:(ii+1)*r1-1, &
                             jj*r1:(jj+1)*r1-1, &
                             1) = .false.
                    endif

                 else if (dim == 3) then
                    if (any(imask(ii*r1:(ii+1)*r1-1, &
                                  jj*r1:(jj+1)*r1-1, &
                                  kk*r1:(kk+1)*r1-1))) then
                    
                       ! get temp and radius
                       cur_temp = p(ii,jj,kk,itemp_comp)
                       cur_r = sqrt( (ii*dx(1))**2 + (jj*dx(2))**2 + (kk*dx(3))**2 )

                       ! if hot, check radius 
                       if(cur_temp > t_hot) then
                         if(cur_r < r_lo) r_lo = cur_r
                         if(cur_r > r_hi) r_hi = cur_r
                       endif

                       ! get current Carbon mass fraction
                       cur_xc =  p(ii,jj,kk,ixc_comp)

                       ! record minimum radius of CO/He interface
                       if(cur_xc < XC_TOL) then
                         if(cur_r < r_interface) r_interface = cur_r
                       endif

                       imask(ii*r1:(ii+1)*r1-1, &
                             jj*r1:(jj+1)*r1-1, &
                             kk*r1:(kk+1)*r1-1) = .false.
                    endif
                 endif

              enddo
           enddo
        enddo

        call fab_unbind(pf, i, j)

     enddo

     ! adjust r1 for the next level
     if (i /= 1) r1 = r1*pf%refrat(i-1,1)
     
  enddo

  ! error checking; this should never happen with non-corrupted data
  !if (any(cnt == 0)) call bl_error("pltfile contains zones with empty data!")

  ! normalize
  !do i = 0, max_points-1
  !      avg(i,:) = avg(i,:) / cnt(i)
  !enddo

  !------------------------------------------------------------
  ! calculate CO/He interface radius                          |
  !------------------------------------------------------------
  !imask = .true.
  !r1 = 1
  !cnt = 0

  !do i = pf%flevel, 1, -1

  !   do j = 1, nboxes(pf, i)

  !      call fab_bind_comp_vec(pf, i, j, icomps_pass)

  !      p => dataptr(pf, i, j)

  !      do kk = lbound(p,dim=3), ubound(p,dim=3)

  !         if (idir == 3) then
  !            index_lo = kk*r1
  !            index_hi = (kk+1)*r1 - 1
  !         endif

  !         do jj = lbound(p,dim=2), ubound(p,dim=2)

  !            if (idir == 2) then
  !               index_lo = jj*r1
  !               index_hi = (jj+1)*r1 - 1
  !            endif

  !            do ii = lbound(p,dim=1), ubound(p,dim=1)

  !               if (idir == 1) then
  !                  index_lo = ii*r1
  !                  index_hi = (ii+1)*r1 - 1
  !               endif

  !               if (dim == 2) then

  !                  if(any(imask(ii*r1:(ii+1)*r1-1, &
  !                               jj*r1:(jj+1)*r1-1, &
  !                               1))) then

  !                     weight = r1**2
  !                  
  !                     do index = index_lo, index_hi
  !                        rms(index,:) = rms(index,:) + &
  !                             (p(ii,jj,kk,:) - avg(index,:))**2 * weight

  !                        cnt(index) = cnt(index) + weight
  !                     enddo

  !                     imask(ii*r1:(ii+1)*r1-1, &
  !                           jj*r1:(jj+1)*r1-1, &
  !                           1) = .false.

  !                  endif

  !               else if (dim == 3) then

  !                  if(any(imask(ii*r1:(ii+1)*r1-1, &
  !                               jj*r1:(jj+1)*r1-1, &
  !                               kk*r1:(kk+1)*r1-1))) then
  !                     weight = r1**3

  !                     do index = index_lo, index_hi
  !                        rms(index,:) = rms(index,:) + &
  !                             (p(ii,jj,kk,:) - avg(index,:))**2 * weight
  !                  
  !                        cnt(index) = cnt(index) + weight
  !                     enddo

  !                     imask(ii*r1:(ii+1)*r1-1, &
  !                           jj*r1:(jj+1)*r1-1, &
  !                           kk*r1:(kk+1)*r1-1) = .false.

  !                  endif

  !               endif

  !            enddo
  !         enddo
  !      enddo

  !      call fab_unbind(pf, i,j)

  !   enddo

  !   if (i/=1) r1 = r1*pf%refrat(i-1,1)

  !enddo

  ! error checking
  !if (any(cnt == 0)) call bl_error("pltfile contains zones with empty data!")

  ! normalize
  !do i = 0, max_points-1
  !   rms(i,:) = sqrt(rms(i,:)/cnt(i))
  !enddo

  !------------------------------------------------------------
  ! output                                                    |
  !------------------------------------------------------------
  
100  format("# time:", 1x, g24.12)

  print *, 'r_lo: ', r_lo
  print *, 'r_hi: ', r_hi
  print *, 'r_interface: ', r_interface
  ! flexible formats based on the number of variables
  ! number of columns = 2*nvars_pass + 1
  !write(columns,'(i10)') 2*nvars_pass+1
  !column_format = '("#",5x,a,5x,' // trim(columns) &
  !     // '(a,"_avg",5x,a,"_rms",5x))'
  !data_format = '(1x,' // trim(columns) // '(g24.12e3,1x))'


  !if (trim(outputfile) == "stdout") then
  !   uout = 6
  !else
  !   uout = unit_new()
  !   open(unit=uout, file=outputfile, status='replace')
  !endif

  !write(uout, 100) pf%tm

  !write(uout,trim(column_format)) dir, &
  !     (trim(pf%names(icomps_pass(j))), &
  !     trim(pf%names(icomps_pass(j))), j = 1, nvars_pass)

  !do i = flo(idir), fhi(idir)
  !   write(uout,trim(data_format)) pf%plo(idir) + (i+HALF)*dx(idir), &
  !        (avg(i,j), rms(i,j), j=1, nvars_pass)
  !enddo

  !if (uout .ne. 6) close(unit=uout)

  call destroy(pf)
  
contains
  !============================================================
  ! fhotspots - subroutine, functions section                 |
  !============================================================
  subroutine print_usage()
    implicit none

    print *, ''
    print *, 'Description: '
    print *, '  This program takes a specified pltfile and calculates the '
    print *, '  average and rms quantities as a function of "height" along a '
    print *, '  specified direction for a specified set of variables. '
    print *, ''
    print *, 'Usage: '
    print *, '  faverage <args> '
    print *, ''
    print *, 'Arguments: '
    print *, '  [-p|--pltfile]    <filename> : '
    print *, '      specifies the input pltfile; (required) '
    print *, ''
    print *, '  [-o|--outputfile] <filename> : '
    print *, '      specifies the output file; default is "output.dat"'
    print *, ''
    print *, '  [-d|--dir]        <integer>  : '
    print *, '      specifies "height" direction; default is 2 (y) for 2-d'
    print *, '      and 3 (z) for 3-d'
    print *, ''
    print *, '  [-v|--vars] <integer> <list> : '
    print *, '      specifies a specific <list> of variables to analyze. '
    print *, '      <integer> is the number of variables in the whitespace '
    print *, '      separated <list> of variable names.  any variables that '
    print *, '      are specified but are not in the pltfile will be skipped; '
    print *, '      default behaviour is to analyze all variables. '
    print *, ''
    print *, '  [-f|--favre]                 : '
    print *, '      toggles favre (density-weighted) averaging; default is to '
    print *, '      do normal averaging '
    print *, ''
    print *, '  [-q|--quiet]                 : '
    print *, '      perform averaging operations quietly'
    print *, ''
  end subroutine print_usage
end program fhotspots
