! Analysis routines for the 3-d WD convect MAESTRO problem 
! These are based on fsedov3d_sph.f90
!
! Here we compute angle-averaged quantities, <q>, and RMS quantities,
! q' = sqrt { < (q - <q>)**2 > }, where the averaging is done at constant
! radius.
!
! For density, rho_0 = <rho>, and in the plotfiles we store rhopert =
! rho - rho0, so we can compute the RMS density directly from this.
!
! Similarly, for temperature, we store tpert = T - <T>, so we will use
! tpert directly from the plotfiles for computing the RMS value.
!
! we also compute some global quantities, including T_peak and its
! location, enuc_peak and its location, and the total kinetic energy
!
! for the total kinetic energy inside the star, we only consider those
! zones where rho > rho_cutoff
!
! --globals_only will just print T_max, it's location, and kinetic energy
!
! This is essentially fwdconvect modified to angle-average all variables
! in the supplied plotfile.
program central_angle_average

  use bl_space
  use bl_error_module
  use bl_constants_module
  use bl_IO_module
  use plotfile_module

  implicit none

  integer, parameter :: MAX_VAR_NAME  = 128
  ! Variable averaging type to use:
  ! 1 = standard volume-weighted averaging
  ! 2 = mass-weighted averaging
  ! 3 = root-mean-square (perturbations)
  integer, parameter :: volume_weighted  = 1
  integer, parameter :: mass_weighted    = 2
  integer, parameter :: rms_perturbation = 3
  
  type(plotfile) pf
  integer :: unit
  character(len=50) :: rfmt ! Slice file row format
  character(len=50) :: hfmt ! Slice file header format
  integer :: i, j, ii, jj, kk, ivar
  real(kind=dp_t) :: xx, yy, zz
  integer :: rr, r1
  integer :: uno

  integer :: nbins
  real(kind=dp_t), allocatable :: r(:)
  real(kind=dp_t) :: maxdist, x_maxdist, y_maxdist, z_maxdist
  real(kind=dp_t) :: xctr, yctr, zctr

  real(kind=dp_t) :: dx(MAX_SPACEDIM)
  real(kind=dp_t) :: dx_fine

  real(kind=dp_t) :: r_zone
  integer :: irbin

  real(kind=dp_t), pointer :: p(:,:,:,:)

  integer, allocatable :: ncount(:)
  ! Profile bin data array. Indices are (variable #, bin #)
  real(kind=dp_t), allocatable :: profile_bins(:,:)
  character(len=MAX_VAR_NAME), allocatable :: var_names(:)
  character(len=MAX_VAR_NAME), allocatable :: column_names(:)
  real(kind=dp_t), allocatable :: row_data(:)
  integer, allocatable :: var_avg_type(:)
  integer, allocatable :: var_indices(:)
  real(kind=dp_t) :: kinetic_energy

  integer :: numvars ! Number of variables in the plotfile
  ! Indexes in the plotfile
  integer :: dens_comp, temp_comp
  integer :: magvel_comp, enuc_comp
  integer :: xvel_comp, yvel_comp, zvel_comp
  ! Index in the binned profile
  integer :: idens 

  logical, allocatable :: imask(:,:,:)
  integer :: lo(MAX_SPACEDIM), hi(MAX_SPACEDIM)
  integer :: flo(MAX_SPACEDIM), fhi(MAX_SPACEDIM)

  character(len=256) :: slicefile
  character(len=256) :: pltfile
  real(kind=dp_t) :: rho_cutoff

  integer :: narg, farg
  character(len=256) :: fname

  real(kind=dp_t) :: time

  real(kind=dp_t) :: T_peak
  real(kind=dp_t) :: xloc_Tpeak, yloc_Tpeak, zloc_Tpeak, R_Tpeak
  real(kind=dp_t) :: vx_Tpeak, vy_Tpeak, vz_Tpeak, vr_Tpeak

  real(kind=dp_t) :: enuc_peak
  real(kind=dp_t) :: xloc_enucpeak, yloc_enucpeak, zloc_enucpeak, R_enucpeak

  logical :: globals_only

  unit = unit_new()
  uno =  unit_new()


  ! set the defaults
  slicefile = ''
  pltfile  = ''
  rho_cutoff = 3.d6
  globals_only = .false.

  narg = command_argument_count()

  farg = 1
  do while ( farg <= narg )
     call get_command_argument(farg, value = fname)

     select case (fname)

     case ('-p', '--pltfile')
        farg = farg + 1
        call get_command_argument(farg, value = pltfile)

     case ('-s', '--slicefile')
        farg = farg + 1
        call get_command_argument(farg, value = slicefile)

     case ('--rho_cutoff')
        farg = farg + 1
        call get_command_argument(farg, value = fname)
        read(fname,*) rho_cutoff

     case ('--globals_only')
        globals_only = .true.

     case default
        exit

     end select
     farg = farg + 1
  end do

  if ( len_trim(pltfile) == 0 .OR. (len_trim(slicefile) == 0 .and. .not. globals_only)) then
     print *, "usage: central_angle_average args"
     print *, "args [-p|--pltfile]   plotfile   : plot file directory (required)"
     print *, "     [-s|--slicefile] slice file : slice file (required if not globals_only)"
     print *, "     --rho_cutoff cutoff value   : low density cutoff (optional)"
     print *, "     --globals_only              : only output global quantities (optional) "
     stop
  end if

  if (.not. globals_only) then
     print *, 'pltfile    = "', trim(pltfile), '"'
     print *, 'slicefile  = "', trim(slicefile), '"'
     print *, 'rho_cutoff = "', rho_cutoff
  endif

  call build(pf, pltfile, unit)

  time = pf%tm

!  do i = 1, pf%flevel
!     call fab_bind_level(pf, i)
!  end do

  ! Get the number of variables in the plotfile
  numvars = plotfile_nvars(pf)

  ! figure out the variable indices
  
  ! density
  dens_comp = plotfile_var_index(pf, "density")

  ! magvel
  magvel_comp = plotfile_var_index(pf, "magvel") 

  ! enuc
  enuc_comp = plotfile_var_index(pf, "enucdot")

  ! temperature (here we used T from p0)
  temp_comp = plotfile_var_index(pf, "tfromp") 

  ! xvel
  xvel_comp = plotfile_var_index(pf, "x_vel")

  ! yvel
  yvel_comp = plotfile_var_index(pf, "y_vel") 

  ! zvel
  zvel_comp = plotfile_var_index(pf, "z_vel") 

  if ( dens_comp < 0 .or. temp_comp < 0 .or. &
       magvel_comp < 0 .or. enuc_comp < 0 .or. &
       xvel_comp < 0 .or. yvel_comp < 0 .or. zvel_comp < 0 ) then
     call bl_error("ERROR: varaible(s) not defined")
  endif
  
  ! get dx for the coarse level.  
  dx = plotfile_get_dx(pf, 1)


  ! get the index bounds for the finest level.  Note, lo and hi are
  ! ZERO based indicies
  flo = lwb(plotfile_get_pd_box(pf, pf%flevel))
  fhi = upb(plotfile_get_pd_box(pf, pf%flevel))

  if (.not. globals_only) then
     print *, 'Size of domain (zones): ', fhi(1)-flo(1)+1, fhi(2)-flo(2)+1, fhi(3)-flo(3)+1
  endif

  ! the default for the center of the star will be the geometric center
  ! of the domain
  xctr = HALF*(pf%phi(1) + pf%plo(1))
  yctr = HALF*(pf%phi(2) + pf%plo(2))
  zctr = HALF*(pf%phi(3) + pf%plo(3))

  if (.not. globals_only) then
     print *, 'Center of the star: ', xctr, yctr, zctr
  endif

  ! compute the size of the radially-binned array -- we'll do it to
  ! the furtherest corner of the domain
  x_maxdist = max(abs(pf%phi(1) - xctr), abs(pf%plo(1) - xctr))
  y_maxdist = max(abs(pf%phi(2) - yctr), abs(pf%plo(2) - yctr))
  z_maxdist = max(abs(pf%phi(3) - zctr), abs(pf%plo(3) - zctr))
  
  maxdist = sqrt(x_maxdist**2 + y_maxdist**2 + z_maxdist**2)

  dx_fine = minval(plotfile_get_dx(pf, pf%flevel))
  nbins = int(maxdist/dx_fine)

  ! Allocate binning data
  allocate(r(0:nbins-1))
  
  ! Construct radius bin centers
  do i = 0, nbins-1
     r(i) = (dble(i) + HALF)*dx_fine
  enddo

  ! imask will be set to false if we've already output the data.
  ! Note, imask is defined in terms of the finest level.  As we loop
  ! over levels, we will compare to the finest level index space to
  ! determine if we've already output here
  allocate(imask(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3)))

  ! ncount will keep track of how many fine zones were written into
  ! each bin
  allocate(   ncount(0:nbins-1))
  
  !----------------------------------------------------------------------------
  ! get plotfile variable names and determine averaging type for each
  !----------------------------------------------------------------------------  
  allocate(var_names(0:numvars-1))
  allocate(var_indices(0:numvars-1))
  allocate(var_avg_type(0:numvars-1))
  do i = 0, numvars-1
     var_names(i) = pf%names(i+1)
     if (var_names(i) == "density") then
        idens = i
     endif
     var_indices(i) = plotfile_var_index(pf, var_names(i))
     if (var_indices(i) == -1) then
        call bl_error("ERROR: variable " // var_names(i) // " not indexed.")
     endif
     ! Is the variable a perturbation? (use RMS)
     jj = index(var_names(i), "pert")
     kk = index(var_names(i), "X(")
     if (jj /= 0) then
        var_avg_type(i) = rms_perturbation
     elseif (kk /= 0) then
        var_avg_type(i) = mass_weighted
     else
        var_avg_type(i) = volume_weighted
     endif
  enddo

  !----------------------------------------------------------------------------
  ! compute the angle averaged quantities
  !----------------------------------------------------------------------------

  ! allocate storage for the data
  allocate(profile_bins(0:numvars-1, 0:nbins-1))
  
  profile_bins(:,:) = ZERO
  
  ncount(:) = 0

  kinetic_energy = ZERO

  imask(:,:,:) = .true.

  T_peak = 0.0
  xloc_Tpeak = 0.0
  yloc_Tpeak = 0.0
  zloc_Tpeak = 0.0

  vx_Tpeak = 0.0
  vy_Tpeak = 0.0
  vz_Tpeak = 0.0
  vr_Tpeak = 0.0

  enuc_peak = 0.0
  xloc_enucpeak = 0.0
  yloc_enucpeak = 0.0
  zloc_enucpeak = 0.0

  ! loop over the data, starting at the finest grid, and if we haven't
  ! already stored data in that grid location (according to imask),
  ! store it.  


  ! r1 is the factor between the current level grid spacing and the
  ! FINEST level
  r1  = 1

  do i = pf%flevel, 1, -1

     ! rr is the factor between the COARSEST level grid spacing and
     ! the current level
     rr = product(pf%refrat(1:i-1,1))

     if (.not. globals_only) then
        print *, 'processing level ', i, ' rr = ', rr
     endif

     do j = 1, nboxes(pf, i)
        
        ! read in the data 1 patch at a time -- read in all the variables
        call fab_bind(pf, i, j)

        lo(:) = 1
        hi(:) = 1
        lo = lwb(get_box(pf, i, j))
        hi = upb(get_box(pf, i, j))

        ! get a pointer to the current patch
        p => dataptr(pf, i, j)

        
        ! loop over all of the zones in the patch.  Here, we convert
        ! the cell-centered indices at the current level into the
        ! corresponding RANGE on the finest level, and test if we've
        ! stored data in any of those locations.  If we haven't then
        ! we store this level's data and mark that range as filled.
        do kk = lo(3), hi(3)
           zz = (kk + HALF)*dx(3)/rr

           do jj = lo(2), hi(2)
              yy = (jj + HALF)*dx(2)/rr

              do ii = lo(1), hi(1)
                 xx = (ii + HALF)*dx(1)/rr

                 if ( any(imask(ii*r1:(ii+1)*r1-1, &
                                jj*r1:(jj+1)*r1-1, &
                                kk*r1:(kk+1)*r1-1) ) ) then

                    r_zone = sqrt((xx-xctr)**2 + (yy-yctr)**2 + (zz-zctr)**2)

                    irbin = r_zone/dx_fine

                    ! weight the zone's data by its size

                    ! note, for p(:,:,:,n), n refers to index of the
                    ! variable as found via plotfile_var_index

                    ! Bin the variables in the plotfile according to their averaging type
                    do ivar = 0, numvars-1
                       if (var_avg_type(ivar) == volume_weighted) then
                          profile_bins(ivar, irbin) = profile_bins(ivar, irbin) + &
                               p(ii,jj,kk,var_indices(ivar)) * r1**3
                       else if (var_avg_type(ivar) == mass_weighted) then
                          profile_bins(ivar, irbin) = profile_bins(ivar, irbin) + &
                               p(ii,jj,kk,dens_comp)*p(ii,jj,kk,var_indices(ivar)) * r1**3
                       else if (var_avg_type(ivar) == rms_perturbation) then
                          profile_bins(ivar, irbin) = profile_bins(ivar, irbin) + &
                               p(ii,jj,kk,var_indices(ivar)) * &
                               p(ii,jj,kk,var_indices(ivar)) * r1**3
                       endif
                    enddo

                    ! kinetic energy is integral of rho U U dV
                    if (p(ii,jj,kk,1) >= rho_cutoff) then
                       kinetic_energy = kinetic_energy + &
                            p(ii,jj,kk,dens_comp) * p(ii,jj,kk,magvel_comp)**2 * &
                            (dx(1)/rr)*(dx(2)/rr)*(dx(3)/rr)
                    endif

                    ncount(irbin) = ncount(irbin) + r1**3


                    ! store the location and value of the peak temperature
                    if (p(ii,jj,kk,temp_comp) > T_peak) then
                       T_peak = p(ii,jj,kk,temp_comp)
                       xloc_Tpeak = xx
                       yloc_Tpeak = yy
                       zloc_Tpeak = zz

                       R_Tpeak = sqrt( (xctr - xloc_Tpeak)**2 + &
                                       (yctr - yloc_Tpeak)**2 + &
                                       (zctr - zloc_Tpeak)**2 )

                       vx_Tpeak = p(ii,jj,kk,xvel_comp)
                       vy_Tpeak = p(ii,jj,kk,yvel_comp)
                       vz_Tpeak = p(ii,jj,kk,zvel_comp)

                       vr_Tpeak = vx_Tpeak*(xx - xctr)/R_Tpeak + &
                                  vy_Tpeak*(yy - yctr)/R_Tpeak + &
                                  vz_Tpeak*(zz - zctr)/R_Tpeak
                       

                    endif


                    ! store the location and value of the peak enucdot
                    if (p(ii,jj,kk,enuc_comp) > enuc_peak) then
                       enuc_peak = p(ii,jj,kk,enuc_comp)

                       xloc_enucpeak = xx
                       yloc_enucpeak = yy
                       zloc_enucpeak = zz

                       R_enucpeak = sqrt( (xctr - xloc_enucpeak)**2 + &
                                          (yctr - yloc_enucpeak)**2 + &
                                          (zctr - zloc_enucpeak)**2 )
                    endif

                    imask(ii*r1:(ii+1)*r1-1, &
                          jj*r1:(jj+1)*r1-1, &
                          kk*r1:(kk+1)*r1-1) = .false.
                 end if

              end do
           enddo
        enddo

        call fab_unbind(pf, i, j)
                
     end do

     ! adjust r1 for the next lowest level
     if ( i /= 1 ) r1 = r1*pf%refrat(i-1,1)
  end do


  ! normalize
  do i = 0, nbins-1
     if (ncount(i) /= 0) then

        ! Normalize the variables in the plotfile according to their averaging type
        do ivar = 0, numvars-1
           if (var_avg_type(ivar) == volume_weighted) then
              profile_bins(ivar, i) = profile_bins(ivar, i)/ncount(i)
           else if (var_avg_type(ivar) == mass_weighted) then
              profile_bins(ivar, i) = (profile_bins(ivar, i)/ncount(i)) / &
                   profile_bins(idens, i)
           else if (var_avg_type(ivar) == rms_perturbation) then
              profile_bins(ivar, i) = sqrt(profile_bins(ivar, i)/ncount(i))
           endif
        enddo
        
        ! we already normalized the kinetic energy by multiplying 
        ! each zone by dV

     endif
  enddo



990  format("# time = ",g24.12)
991  format("# ---------------------------------------------------------------------------")
992  format("# total kinetic energy = ", g24.12)
994  format("# peak temperature = ", g24.12)
995  format("# peak temp loc (x,y,z) = ", 3(g24.12,1x))
996  format("# peak temp radius = ", g24.12)
9961 format("# velocity @ peak T loc (vx, vy, vz) = ", 3(g24.12,1x))
9962 format("# radial velocity @ peak T loc = ", g24.12,1x)
997  format("# peak enucdot = ", g24.12)
998  format("# peak enucdot loc (x,y,z) = ", 3(g24.12,1x))
999  format("# peak enucdot radius = ", g24.12)
1000 format("#",100(a24,1x))
1001 format(1x, 100(g24.12,1x))

  ! slicefile
  if (.not. globals_only) then
     open(unit=uno, file=slicefile, status = 'replace')

     ! write the header
     write(uno,990) time
     write(uno,991)

     write(uno,992) kinetic_energy
     write(uno,991)

     write(uno,994) T_peak
     write(uno,995) xloc_Tpeak, yloc_Tpeak, zloc_Tpeak     
     write(uno,996) R_Tpeak
     write(uno,9961) vx_Tpeak, vy_Tpeak, vz_Tpeak
     write(uno,9962) vr_Tpeak

     write(uno,997) enuc_peak
     write(uno,998) xloc_enucpeak, yloc_enucpeak, zloc_enucpeak
     write(uno,999) R_enucpeak
     write(uno,991)

     ! Set format for writing column entries
     write(hfmt,'(A,I5,A)') '(', numvars+1, '(A24))'
     write(rfmt,'(A,I5,A)') '(', numvars+1, '(g24.12))'

     ! Create column header names
     allocate(column_names(0:numvars))
     column_names(0) = "r"
     column_names(1:numvars) = var_names(0:numvars-1)

     ! Allocate row data
     allocate(row_data(0:numvars))
     row_data = ZERO
     
     write(uno, fmt=hfmt) (column_names(j), j = 0, numvars)

     ! write the data in columns
     do i = 0, nbins-1
        ! Fill row
        ! Use this to protect against a number being xx.e-100
        !   which will print without the "e"
        row_data(0) = r(i)
        do j = 0, numvars-1
           if (abs(profile_bins(j, i)) .lt. 1.d-99) profile_bins(j, i) = 0.d0
           row_data(j+1) = profile_bins(j, i)
        enddo
        
        ! Write row data
        write(uno, fmt=rfmt) (row_data(j), j = 0, numvars)
     end do

     close(unit=uno)
  endif

  if (.not. globals_only) then
     print *, 'Total kinetic energy = ', kinetic_energy
     print *, 'Peak temperature = ', T_peak
     print *, 'Peak temperature location (x,y,z) = ', xloc_Tpeak, yloc_Tpeak, zloc_Tpeak
     print *, 'Peak temperature radius = ', R_Tpeak
     print *, 'Peak enucdot = ', enuc_peak
     print *, 'Peak enucdot location (x,y,z) = ', xloc_enucpeak, yloc_enucpeak, zloc_enucpeak
     print *, 'Peak enucdot radius = ', R_enucpeak
  else
     write (*,1000) "time", "T_peak", "x(T_peak)", "y(T_peak)", "z(T_peak)", "R(T_peak)", &
          "enuc_peak", "x(enuc_peak)", "y(enuc_peak)", "z(enuc_peak)", "R(enuc_peak)", "KE"
     write (*,1001) time, T_peak, xloc_Tpeak, yloc_Tpeak, zloc_Tpeak, R_Tpeak, &
          enuc_peak, xloc_enucpeak, yloc_enucpeak, zloc_enucpeak, R_enucpeak, kinetic_energy
  endif

  do i = 1, pf%flevel
     call fab_unbind_level(pf, i)
  end do

  call destroy(pf)

end program central_angle_average
