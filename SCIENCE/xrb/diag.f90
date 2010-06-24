! xrb-specific diagnostic routine
!
! output files:
!
!    xrb_enuc_diag.out:
!        peak nuclear energy generation rate (erg / g / s)
!        x/y/z location of peak 
!        total mass of C12
!        total mass of O16 (if using the triple_alpha_plus_cago network
!
!    xrb_temp_diag.out:
!        peak temperature in the helium layer
!            helium layer defined by X(He4) >= diag_define_he_layer where 
!            diag_define_he_layer is a xrb _parameter
!        x/y/z location of peak
!
!    xrb_vel_diag.out:
!        peak velocity
!        x/y/z location of peak velocity
!        peak Mach number
!        x/y/z location of peak Mach number
!
! optional output files:
!
!  if do_deltap_diag .eq. T:
!    xrb_deltap_diag.out:
!        peak deltap
!        location of peak deltap
!        average deltap
!        
!


module diag_module

  use bl_types
  use bl_IO_module
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use network, only: network_species_index
  use probin_module, only: do_deltap_diag

  implicit none

  private

  public :: diag, flush_diag

  logical, save :: network_contains_oxygen = .false.

contains

  subroutine diag(time,dt,dx,s,rho_Hnuc,rho_Hext,thermal,rho_omegadot, &
                  rho0,rhoh0,p0,tempbar, &
                  gamma1bar,div_coeff, &
                  u,w0,normal, &
                  mla,the_bc_tower)

    use bl_prof_module
    use geometry, only: dm, nlevs
    use bl_constants_module

    real(kind=dp_t), intent(in   ) :: dt,dx(:,:),time
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: rho_Hnuc(:)
    type(multifab) , intent(in   ) :: rho_Hext(:)
    type(multifab) , intent(in   ) :: thermal(:)
    type(multifab) , intent(in   ) :: rho_omegadot(:)    
    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) ::      rho0(:,0:)
    real(kind=dp_t), intent(in   ) ::     rhoh0(:,0:)
    real(kind=dp_t), intent(in   ) ::        p0(:,0:)
    real(kind=dp_t), intent(in   ) ::   tempbar(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: div_coeff(:,0:)
    real(kind=dp_t), intent(in   ) ::        w0(:,0:)
    type(ml_layout), intent(in   ) :: mla
    type(bc_tower) , intent(in   ) :: the_bc_tower


    ! Local
    real(kind=dp_t), pointer::  sp(:,:,:,:)
    real(kind=dp_t), pointer::  rhnp(:,:,:,:)
    real(kind=dp_t), pointer::  rhep(:,:,:,:)
    real(kind=dp_t), pointer::  up(:,:,:,:)
    logical,         pointer :: mp(:,:,:,:)


    real(kind=dp_t) :: T_max, T_max_level, T_max_local
    real(kind=dp_t) :: coord_T_max(dm), coord_T_max_level(dm), &
                       coord_T_max_local(dm)

    real(kind=dp_t) :: enuc_max, enuc_max_level, enuc_max_local
    real(kind=dp_t) :: coord_enuc_max(dm), coord_enuc_max_level(dm), &
                       coord_enuc_max_local(dm)

    real(kind=dp_t) :: total_c12_mass, total_c12_mass_level, &
                       total_c12_mass_local
    real(kind=dp_t) :: total_o16_mass, total_o16_mass_level, &
                       total_o16_mass_local

    real(kind=dp_t) :: vel_max, vel_max_level, vel_max_local
    real(kind=dp_t) :: coord_vel_max(dm), coord_vel_max_level(dm), &
                       coord_vel_max_local(dm)

    real(kind=dp_t) :: Machno_max, Machno_max_level, Machno_max_local
    real(kind=dp_t) :: coord_Machno_max(dm), coord_Machno_max_level(dm), &
                       coord_Machno_max_local(dm)

    real(kind=dp_t) :: deltap_max, deltap_max_level, deltap_max_local
    real(kind=dp_t) :: coord_deltap_max(dm), coord_deltap_max_level(dm), &
                       coord_deltap_max_local(dm)

    real(kind=dp_t) :: deltap_avg, deltap_avg_level, deltap_avg_local
    real(kind=dp_t) :: nzones, nzones_level, nzones_local

    ! buffers
    real(kind=dp_t) :: T_max_data_local(1), T_max_coords_local(dm)
    real(kind=dp_t), allocatable :: T_max_data(:), T_max_coords(:)

    real(kind=dp_t) :: enuc_max_data_local(1), enuc_max_coords_local(dm)
    real(kind=dp_t), allocatable :: enuc_max_data(:), enuc_max_coords(:)

    real(kind=dp_t), allocatable :: sum_data_level(:), sum_data_local(:)
    integer :: nsums

    real(kind=dp_t) :: vel_max_data_local(1), vel_max_coords_local(dm)
    real(kind=dp_t), allocatable :: vel_max_data(:), vel_max_coords(:)
    
    real(kind=dp_t) :: Machno_max_data_local(1), Machno_max_coords_local(dm)
    real(kind=dp_t), allocatable :: Machno_max_data(:), Machno_max_coords(:)

    real(kind=dp_t) :: deltap_max_data_local(1), deltap_max_coords_local(dm)
    real(kind=dp_t), allocatable :: deltap_max_data(:), deltap_max_coords(:)

    integer :: lo(dm),hi(dm),ng_s,ng_u,ng_rhn,ng_rhe
    integer :: i,n, index_max, isum
    integer :: un, un2, un3, un4
    logical :: lexist

    logical, save :: firstCall = .true.

    type(bl_prof_timer), save :: bpt


    call build(bpt, "diagnostics")

    ! find out if we are using the triple_alpha_plus_cago network
    if (firstCall .and. network_species_index("oxygen-16") .ne. -1) &
         network_contains_oxygen = .true.

    ng_s = s(1)%ng
    ng_u = u(1)%ng
    ng_rhn = rho_Hnuc(1)%ng
    ng_rhe = rho_Hext(1)%ng

    ! initialize
    ! note that T_max corresponds to the maximum temperature in the 
    ! helium layer defined by X(He4) >= diag_define_he_layer
    T_max          = ZERO
    coord_T_max(:) = ZERO

    enuc_max          = ZERO
    coord_enuc_max(:) = ZERO

    total_c12_mass = ZERO

    nsums = 1

    if (network_contains_oxygen) then
       total_o16_mass = ZERO
       nsums = nsums + 1
    endif

    vel_max          = ZERO
    coord_vel_max(:) = ZERO

    Machno_max          = ZERO
    coord_Machno_max(:) = ZERO

    if (do_deltap_diag) then
       deltap_max = ZERO
       coord_deltap_max(:) = ZERO
       
       deltap_avg = ZERO

       nzones = 0

       ! one for the avg and one for the number of zones
       nsums = nsums + 2
    endif

    ! loop over the levels and calculate global quantities
    do n = 1, nlevs

       ! initialize local and level quantities
       T_max_local = ZERO
       T_max_level = ZERO

       coord_T_max_level(:) = ZERO
       coord_T_max_local(:) = ZERO

       enuc_max_local = ZERO
       enuc_max_level = ZERO

       coord_enuc_max_level(:) = ZERO
       coord_enuc_max_local(:) = ZERO
       
       total_c12_mass_local = ZERO
       total_c12_mass_level = ZERO

       if (network_contains_oxygen) then
          total_o16_mass_local = ZERO
          total_o16_mass_level = ZERO
       endif
       
       vel_max_local = ZERO
       vel_max_level = ZERO
       
       coord_vel_max_level(:) = ZERO
       coord_vel_max_local(:) = ZERO
       
       Machno_max_level = ZERO
       Machno_max_local = ZERO

       coord_Machno_max_level(:) = ZERO
       coord_Machno_max_local(:) = ZERO

       if (do_deltap_diag) then
          deltap_max_level = ZERO
          deltap_max_local = ZERO
          
          coord_deltap_max_level(:) = ZERO
          coord_deltap_max_local(:) = ZERO

          deltap_avg_level = ZERO
          deltap_avg_local = ZERO

          nzones_level = 0
          nzones_local = 0
       endif

       allocate(sum_data_local(nsums), sum_data_level(nsums))

       ! loop over the boxes at the current level
       do i = 1, s(n)%nboxes

          if ( multifab_remote(s(n), i) ) cycle

          sp => dataptr(s(n) , i)
          rhnp => dataptr(rho_Hnuc(n), i)
          rhep => dataptr(rho_Hext(n), i)
          up => dataptr(u(n) , i)
          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))

          ! build the local data
          select case (dm)
          case (2)
             ! only do those boxes that aren't masked
             if (n .eq. nlevs) then

                call diag_2d(n,time,dt,dx(n,:), &
                             sp(:,:,1,:),ng_s, &
                             rhnp(:,:,1,1),ng_rhn, &
                             rhep(:,:,1,1),ng_rhe, &
                             rho0(n,:),rhoh0(n,:), &
                             p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                             up(:,:,1,:),ng_u, &
                             w0(n,:), &
                             lo,hi, &
                             T_max_local, &
                             coord_T_max_local, &
                             enuc_max_local, &
                             coord_enuc_max_local, &
                             total_c12_mass_local, &
                             vel_max_local, &
                             coord_vel_max_local, &
                             Machno_max_local, &
                             coord_Machno_max_local, &
                             total_o16_mass_local, &
                             deltap_max_local, &
                             coord_deltap_max_local, &
                             deltap_avg_local, &
                             nzones_local)

             else
                mp => dataptr(mla%mask(n),i)

                call diag_2d(n,time,dt,dx(n,:), &
                             sp(:,:,1,:),ng_s, &
                             rhnp(:,:,1,1),ng_rhn, &
                             rhep(:,:,1,1),ng_rhe, &
                             rho0(n,:),rhoh0(n,:), &
                             p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                             up(:,:,1,:),ng_u, &
                             w0(n,:), &
                             lo,hi, &
                             T_max_local, &
                             coord_T_max_local, &
                             enuc_max_local, &
                             coord_enuc_max_local, &
                             total_c12_mass_local, &
                             vel_max_local, &
                             coord_vel_max_local, &
                             Machno_max_local, &
                             coord_Machno_max_local, &
                             total_o16_mass_local, &
                             deltap_max_local, &
                             coord_deltap_max_local, &
                             deltap_avg_local, &
                             nzones_local, &
                             mask=mp(:,:,1,1))

             endif
          case (3)
             ! only do those boxes that aren't masked
             if (n .eq. nlevs) then

                call diag_3d(n,time,dt,dx(n,:), &
                             sp(:,:,:,:),ng_s, &
                             rhnp(:,:,:,1),ng_rhn, &
                             rhep(:,:,:,1),ng_rhe, &
                             rho0(n,:),rhoh0(n,:), &
                             p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                             up(:,:,:,:),ng_u, &
                             w0(n,:), &
                             lo,hi, &
                             T_max_local, &
                             coord_T_max_local, &
                             enuc_max_local, &
                             coord_enuc_max_local, &
                             total_c12_mass_local, &
                             vel_max_local, &
                             coord_vel_max_local, &
                             Machno_max_local, &
                             coord_Machno_max_local, &
                             total_o16_mass_local, &
                             deltap_max_local, &
                             coord_deltap_max_local, &
                             deltap_avg_local, &
                             nzones_local)

             else
                mp => dataptr(mla%mask(n),i)

                call diag_3d(n,time,dt,dx(n,:), &
                             sp(:,:,:,:),ng_s, &
                             rhnp(:,:,:,1),ng_rhn, &
                             rhep(:,:,:,1),ng_rhe, &
                             rho0(n,:),rhoh0(n,:), &
                             p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                             up(:,:,:,:),ng_u, &
                             w0(n,:), &
                             lo,hi, &
                             T_max_local, &
                             coord_T_max_local, &
                             enuc_max_local, &
                             coord_enuc_max_local, &
                             total_c12_mass_local, &
                             vel_max_local, &
                             coord_vel_max_local, &
                             Machno_max_local, &
                             coord_Machno_max_local, &
                             total_o16_mass_local, &
                             deltap_max_local, &
                             coord_deltap_max_local, &
                             deltap_avg_local, &
                             nzones_local, &
                             mask=mp(:,:,:,1))

             endif

          end select

       end do

       ! build the level data
       
       ! pack the summed quantities to reduce communications
       isum = 1
       sum_data_local(isum) = total_c12_mass_local
       if (network_contains_oxygen) then
          isum = isum + 1
          sum_data_local(isum) = total_o16_mass_local
       endif
       if (do_deltap_diag) then
          isum = isum + 1
          sum_data_local(isum) = deltap_avg_local
          isum = isum + 1
          sum_data_local(isum) = nzones_local
       endif

       call parallel_reduce(sum_data_level, sum_data_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       isum = 1
       total_c12_mass_level = sum_data_level(isum)
       if (network_contains_oxygen) then
          isum = isum + 1
          total_o16_mass_level = sum_data_level(isum)
       endif
       if (do_deltap_diag) then
          isum = isum + 1
          deltap_avg_level = sum_data_level(isum)
          isum = isum + 1
          nzones_level = sum_data_level(isum)
       endif

       ! get the maximum temperature and its location
       ! gather all the T_max data into an array; find the index of the maximum
       allocate(T_max_data(parallel_nprocs()))
       T_max_data_local(1) = T_max_local

       call parallel_gather(T_max_data_local, T_max_data, 1, &
                            root = parallel_IOProcessorNode())
       
       index_max = maxloc(T_max_data, dim=1)

       ! gather all the T_max_coords into an array and use index_max to 
       ! get the correct location
       allocate(T_max_coords(dm*parallel_nprocs()))
       T_max_coords_local(:) = coord_T_max_local(:)

       call parallel_gather(T_max_coords_local , T_max_coords, dm, &
                            root = parallel_IOProcessorNode())

       T_max_level = T_max_data(index_max)
       
       coord_T_max_level(1) = T_max_coords(dm*(index_max-1) + 1)
       coord_T_max_level(2) = T_max_coords(dm*(index_max-1) + 2)
       if (dm>2) coord_T_max_level(3) = T_max_coords(dm*(index_max-1) + 3)

       deallocate(T_max_data, T_max_coords)

       ! get the maximum enuc and its location
       allocate(enuc_max_data(parallel_nprocs()))
       enuc_max_data_local(1) = enuc_max_local

       call parallel_gather(enuc_max_data_local, enuc_max_data, 1, &
                            root = parallel_IOProcessorNode())

       index_max = maxloc(enuc_max_data, dim=1)

       allocate(enuc_max_coords(dm*parallel_nprocs()))
       enuc_max_coords_local(:) = coord_enuc_max_local(:)

       call parallel_gather(enuc_max_coords_local, enuc_max_coords, dm, &
                            root = parallel_IOProcessorNode())

       enuc_max_level = enuc_max_data(index_max)

       coord_enuc_max_level(1) = enuc_max_coords(dm*(index_max-1) + 1)
       coord_enuc_max_level(2) = enuc_max_coords(dm*(index_max-1) + 2)
       if(dm>2) coord_enuc_max_level(3) = enuc_max_coords(dm*(index_max-1) + 3)

       deallocate(enuc_max_data, enuc_max_coords)

       ! get the maximum of vel and its location
       allocate(vel_max_data(parallel_nprocs()))
       vel_max_data_local(1) = vel_max_local

       call parallel_gather(vel_max_data_local, vel_max_data, 1, &
                            root = parallel_IOProcessorNode())
       
       index_max = maxloc(vel_max_data, dim=1)

       allocate(vel_max_coords(dm*parallel_nprocs()))
       vel_max_coords_local(:) = coord_vel_max_local(:)

       call parallel_gather(vel_max_coords_local, vel_max_coords, dm, &
                            root = parallel_IOProcessorNode())

       vel_max_level = vel_max_data(index_max)

       coord_vel_max_level(1) = vel_max_coords(dm*(index_max-1) + 1)
       coord_vel_max_level(2) = vel_max_coords(dm*(index_max-1) + 2)
       if (dm>2) coord_vel_max_level(3) = vel_max_coords(dm*(index_max-1) + 3)

       deallocate(vel_max_data, vel_max_coords)

       ! get the maximum of Mach number and its location
       allocate(Machno_max_data(parallel_nprocs()))
       Machno_max_data_local(1) = Machno_max_local
       
       call parallel_gather(Machno_max_data_local, Machno_max_data, 1, &
                            root = parallel_IOProcessorNode())

       index_max = maxloc(Machno_max_data, dim=1)

       allocate(Machno_max_coords(dm*parallel_nprocs()))
       Machno_max_coords_local(:) = coord_Machno_max_local(:)

       call parallel_gather(Machno_max_coords_local, Machno_max_coords, dm, &
                            root = parallel_IOProcessorNode())

       Machno_max_level = Machno_max_data(index_max)

       coord_Machno_max_level(1) = Machno_max_coords(dm*(index_max-1) + 1)
       coord_Machno_max_level(2) = Machno_max_coords(dm*(index_max-1) + 2)
       if (dm>2) coord_Machno_max_level(3) = Machno_max_coords(dm*(index_max-1) + 3)

       deallocate(Machno_max_data, Machno_max_coords)

       ! get the maximum of deltap and its location
       if (do_deltap_diag) then
          allocate(deltap_max_data(parallel_nprocs()))
          deltap_max_data_local(1) = deltap_max_local

          call parallel_gather(deltap_max_data_local, deltap_max_data, 1, &
                               root = parallel_IOProcessorNode())

          index_max = maxloc(deltap_max_data, dim=1)

          allocate(deltap_max_coords(dm*parallel_nprocs()))
          deltap_max_coords_local(:) = coord_deltap_max_local(:)

          call parallel_gather(deltap_max_coords_local, deltap_max_coords, &
                               dm, &
                               root = parallel_IOProcessorNode())

          deltap_max_level = deltap_max_data(index_max)

          coord_deltap_max_level(1) = deltap_max_coords(dm*(index_max-1) + 1)
          coord_deltap_max_level(2) = deltap_max_coords(dm*(index_max-1) + 2)
          if (dm>2) coord_deltap_max_level(3) = deltap_max_coords(dm*(index_max-1) + 3)

          deallocate(deltap_max_data, deltap_max_coords)

       endif

       ! reduce the current level's data with the global data
       if (parallel_IOProcessor()) then

          if (T_max_level > T_max) then
             T_max = T_max_level

             coord_T_max(:) = coord_T_max_level(:)

          endif

          if (enuc_max_level > enuc_max) then
             enuc_max = enuc_max_level

             coord_enuc_max(:) = coord_enuc_max_level(:)

          endif

          total_c12_mass = total_c12_mass + total_c12_mass_level
          if (network_contains_oxygen) &
               total_o16_mass = total_o16_mass + total_o16_mass_level

          if (vel_max_level > vel_max) then
             vel_max = vel_max_level

             coord_vel_max(:) = coord_vel_max_level(:)

          endif

          if (Machno_max_level > Machno_max) then
             Machno_max = Machno_max_level

             coord_Machno_max(:) = coord_Machno_max_level(:)

          endif

          if (do_deltap_diag) then

             nzones     = nzones + nzones_level
             deltap_avg = deltap_avg + deltap_avg_level
             
             if (deltap_max_level > deltap_max) then

                deltap_max = deltap_max_level

                coord_deltap_max(:) = coord_deltap_max_level(:)
                
             endif

          endif

       endif

    end do

    ! normalize 
    ! we weight things in the loop over zones by the coarse level resolution
    total_c12_mass = total_c12_mass * dx(1,1) * dx(1,2)

    if (network_contains_oxygen) &
         total_o16_mass = total_o16_mass * dx(1,1) * dx(1,2)

    if (dm>2) then
       total_c12_mass = total_c12_mass * dx(1,3)

       if (network_contains_oxygen) &
            total_o16_mass = total_o16_mass * dx(1,3)

    endif

    if (do_deltap_diag) deltap_avg = deltap_avg / nzones

1000 format(1x,10(g20.10,1x))
1001 format("#",10(a18,1x))

    if (parallel_IOProcessor()) then

       ! open the diagnostic files for output, taking care not to overwrite
       ! an existing file
       un = unit_new()
       inquire(file="xrb_temp_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un, file="xrb_temp_diag.out", &
               status="old", position="append")
       else
          open(unit=un, file="xrb_temp_diag.out", status="new")
       endif

       un2 = unit_new()
       inquire(file="xrb_enuc_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un2, file="xrb_enuc_diag.out", &
               status="old", position="append")
       else
          open(unit=un2, file="xrb_enuc_diag.out", status="new")
       endif

       un3 = unit_new()
       inquire(file="xrb_vel_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un3, file="xrb_vel_diag.out", &
               status="old", position="append")
       else
          open(unit=un3, file="xrb_vel_diag.out", status="new")
       endif

       if (do_deltap_diag) then
          un4 = unit_new()

          inquire(file="xrb_deltap_diag.out", exist=lexist)
          if (lexist) then
             open(unit=un4, file="xrb_deltap_diag.out", &
                  status="old", position="append")
          else
             open(unit=un4, file="xrb_deltap_diag.out", status="new")
          endif

       endif

       ! print the headers
       if(firstCall) then
          
          if (dm > 2) then
             write(un , 1001) "time", "max{T}", "x_loc", "y_loc", "z_loc"

             if (network_contains_oxygen) then
                write(un2, 1001) "time", "max{enuc}", &
                                 "x_loc", "y_loc", "z_loc", &
                                 "mass_c12", "mass_o16"
             else
                write(un2, 1001) "time", "max{enuc}", &
                                 "x_loc", "y_loc", "z_loc", &
                                 "mass_c12"
             endif
                
             write(un3, 1001) "time", "max{vel}", "x_loc", "y_loc", "z_loc", &
                  "max{Machno}", "x_loc", "y_loc", "z_loc"

             if (do_deltap_diag) &
                write(un4, 1001) "time", "max{deltap}", &
                     "x_loc", "y_loc", "z_loc", "avg{deltap}"

          else
             write(un , 1001) "time", "max{T}", "x_loc", "y_loc"
             
             if (network_contains_oxygen) then
                write(un2, 1001) "time", "max{enuc}", "x_loc", "y_loc", &
                                 "mass_c12", "mass_o16"
             else
                write(un2, 1001) "time", "max{enuc}", "x_loc", "y_loc", &
                                 "mass_c12"
             endif

             write(un3, 1001) "time", "max{vel}", "x_loc", "y_loc", &
                  "max{Machno}", "x_loc", "y_loc"

             if (do_deltap_diag) &
                  write(un4, 1001) "time", "max{deltap}", "x_loc", "y_loc", &
                                   "avg{deltap}"
          endif

          firstCall = .false.
       endif

       ! print the data
       write(un,  1000) time, T_max, (coord_T_max(i), i=1,dm)

       if (network_contains_oxygen) then
          write(un2, 1000) time, enuc_max, (coord_enuc_max(i), i=1,dm), &
                           total_c12_mass, total_o16_mass
       else
          write(un2, 1000) time, enuc_max, (coord_enuc_max(i), i=1,dm), &
                           total_c12_mass
       endif

       write(un3, 1000) time, vel_max, (coord_vel_max(i), i=1,dm), &
            Machno_max, (coord_Machno_max(i), i=1,dm)

       if (do_deltap_diag) &
            write(un4, 1000) time, deltap_max, (coord_deltap_max(i), i=1,dm), &
                             deltap_avg

       close(un )
       close(un2)
       close(un3)
       if (do_deltap_diag) close(un4)
    endif

    call destroy(bpt)

  end subroutine diag


  subroutine flush_diag()
    ! flush_diag is called immediately before checkpointing.  If an
    ! implementation of these diagnostic routines wants to buffer the
    ! data a write out a lot of timestep's worth of information all 
    ! at once, flush_diag() is the routine that should do the writing.

  end subroutine flush_diag


  subroutine diag_2d(n,time,dt,dx, &
                     s,ng_s, &
                     rho_Hnuc,ng_rhn, &
                     rho_Hext,ng_rhe, &
                     rho0,rhoh0,p0,tempbar,gamma1bar, &
                     u,ng_u, &
                     w0, &
                     lo,hi, &
                     T_max, &
                     coord_T_max, &
                     enuc_max, &
                     coord_enuc_max, &
                     c12_mass, &
                     vel_max, &
                     coord_vel_max, &
                     Machno_max, &
                     coord_Machno_max, &
                     o16_mass, &
                     deltap_max, &
                     coord_deltap_max, &
                     deltap_avg, &
                     nzones, &
                     mask)

    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp
    use network, only: nspec
    use probin_module, only: diag_define_he_layer, prob_lo, base_cutoff_density
    use bl_constants_module, only: HALF, ONE, FOUR
    use eos_module

    integer, intent(in) :: n, lo(:), hi(:), ng_s, ng_u, ng_rhn, ng_rhe
    real (kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho_Hnuc(lo(1)-ng_rhn:,lo(2)-ng_rhn:)
    real (kind=dp_t), intent(in   ) :: rho_Hext(lo(1)-ng_rhe:,lo(2)-ng_rhe:)
    real (kind=dp_t), intent(in   ) :: rho0(0:), rhoh0(0:), &
                                         p0(0:),tempbar(0:),gamma1bar(0:)
    real (kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: w0(0:)
    real (kind=dp_t), intent(in   ) :: time, dt, dx(:)
    real (kind=dp_t), intent(inout) :: T_max, enuc_max
    real (kind=dp_t), intent(inout) :: coord_T_max(2), coord_enuc_max(2)
    real (kind=dp_t), intent(inout) :: c12_mass
    real (kind=dp_t), intent(inout) :: vel_max, Machno_max
    real (kind=dp_t), intent(inout) :: coord_vel_max(2), coord_Machno_max(2)
    real (kind=dp_t), intent(inout) :: o16_mass
    real (kind=dp_t), intent(inout) :: deltap_max, coord_deltap_max(2)
    real (kind=dp_t), intent(inout) :: deltap_avg, nzones
    logical,          intent(in   ), optional :: mask(lo(1):,lo(2):)

    !     Local variables
    integer                     :: i, j
    real (kind=dp_t) :: x, y 
    real (kind=dp_t) :: enuc_local
    logical :: cell_valid
    real (kind=dp_t) :: weight
    real (kind=dp_t) :: w0_cent, vel, Mach, deltap

    ! weight is the volume of a cell at the current level divided by the
    ! volume of a cell at the COARSEST level
    weight = ONE / FOUR**(n-1)

    do j = lo(2), hi(2)
       y = prob_lo(2) + (dble(j) + HALF) * dx(2)

       ! recall that w0 is edge centered
       w0_cent = HALF * (w0(j) + w0(j+1))

       do i = lo(1), hi(1)
          x = prob_lo(1) + (dble(i) + HALF) * dx(1)

          cell_valid = .true.
          if (present(mask)) then
             if ( (.not. mask(i,j)) ) cell_valid = .false.
          endif

          if (cell_valid .and. s(i,j,rho_comp) > base_cutoff_density) then

             ! temperature diagnostic
             ! check to see if we are in the helium layer
             ! if we are, then get T_max and its loc
             if ( s(i,j,spec_comp) .ge. &
		  diag_define_he_layer * s(i,j,rho_comp) ) then

                if (s(i,j,temp_comp) > T_max) then
                
                   T_max = s(i,j,temp_comp)
                
                   coord_T_max(1) = x
                   coord_T_max(2) = y

                endif
             endif

             ! enuc diagnostic
             enuc_local = rho_Hnuc(i,j)/s(i,j,rho_comp)
             if (enuc_local > enuc_max) then

                enuc_max = enuc_local
             
                coord_enuc_max(1) = x
                coord_enuc_max(2) = y

             endif

             ! c12 mass diagnostic
             c12_mass = c12_mass + weight*s(i,j,spec_comp+1)

             ! o16 mass diagnostic
             if (network_contains_oxygen) &
                  o16_mass = o16_mass + weight*s(i,j,spec_comp+2)

             ! vel diagnostic
             vel = sqrt(u(i,j,1)**2 + (u(i,j,2) + w0_cent)**2)
             if (vel > vel_max) then
                
                vel_max = vel

                coord_vel_max(1) = x
                coord_vel_max(2) = y
             endif
             
             ! Mach number diagnostic
             ! call the EOS to get the sound speed
             den_eos(1) = s(i,j,rho_comp)
             temp_eos(1) = s(i,j,temp_comp)
             h_eos(1) = s(i,j,rhoh_comp) / s(i,j,rho_comp)

             xn_eos(1,:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)

             call eos(eos_input_rh, den_eos, temp_eos, &
                      npts, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      .false.)

             Mach = vel/cs_eos(1)
             if (Mach > Machno_max) then

                Machno_max = Mach

                coord_Machno_max(1) = x
                coord_Machno_max(2) = y
             endif

             ! deltap diagnostic
             if (do_deltap_diag) then
                ! here we mirror what we do in make_tfromH
                ! we already have the pressure from the above eos call
                deltap = abs(p_eos(1) - p0(j)) / p0(j)

                if (deltap > deltap_max) then
                   
                   deltap_max = deltap

                   coord_deltap_max(1) = x
                   coord_deltap_max(2) = y

                endif

                deltap_avg = deltap_avg + weight*deltap
                nzones = nzones + weight
                
             endif

          endif

       enddo
    enddo

  end subroutine diag_2d

  subroutine diag_3d(n,time,dt,dx, &
                     s,ng_s, &
                     rho_Hnuc,ng_rhn, &
                     rho_Hext,ng_rhe, &
                     rho0,rhoh0,p0,tempbar,gamma1bar, &
                     u,ng_u,w0, &
                     lo,hi, &
                     T_max, &
                     coord_T_max, &
                     enuc_max, &
                     coord_enuc_max, &
                     c12_mass, &
                     vel_max, &
                     coord_vel_max, &
                     Machno_max, &
                     coord_Machno_max, &
                     o16_mass, &
                     deltap_max, &
                     coord_deltap_max, &
                     deltap_avg, &
                     nzones, &
                     mask)

    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp
    use network, only: nspec
    use bl_constants_module, only: HALF, ONE, EIGHT
    use probin_module, only: diag_define_he_layer, prob_lo, base_cutoff_density
    use eos_module

    integer, intent(in) :: n,lo(:), hi(:), ng_s, ng_u, ng_rhn, ng_rhe
    real (kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho_Hnuc(lo(1)-ng_rhn:,lo(2)-ng_rhn:,lo(3)-ng_rhn:)
    real (kind=dp_t), intent(in   ) :: rho_Hext(lo(1)-ng_rhe:,lo(2)-ng_rhe:,lo(3)-ng_rhe:)
    real (kind=dp_t), intent(in   ) :: rho0(0:), rhoh0(0:), &
                                         p0(0:),tempbar(0:),gamma1bar(0:)
    real (kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: w0(0:)
    real (kind=dp_t), intent(in   ) :: time, dt, dx(:)
    real (kind=dp_t), intent(inout) :: T_max, enuc_max
    real (kind=dp_t), intent(inout) :: coord_T_max(3), coord_enuc_max(3)
    real (kind=dp_t), intent(inout) :: c12_mass
    real (kind=dp_t), intent(inout) :: vel_max, Machno_max
    real (kind=dp_t), intent(inout) :: coord_vel_max(3), coord_Machno_max(3)
    real (kind=dp_t), intent(inout) :: o16_mass
    real (kind=dp_t), intent(inout) :: deltap_max, coord_deltap_max(3)
    real (kind=dp_t), intent(inout) :: deltap_avg, nzones
    logical,          intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    !     Local variables
    integer                     :: i, j, k
    real (kind=dp_t) :: x, y , z
    real (kind=dp_t) :: enuc_local
    logical :: cell_valid
    real (kind=dp_t) :: weight
    real (kind=dp_t) :: w0_cent, vel, Mach, deltap


    ! weight is the volume of a cell at the current level divided by the
    ! volume of a cell at the COARSEST level
    weight = ONE / EIGHT**(n-1)

    do k = lo(3), hi(3)
       z = prob_lo(3) + (dble(k) + HALF) * dx(3)       
       
       ! recall w0 is edge centered
       w0_cent = HALF * (w0(k) + w0(k+1))

       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j) + HALF) * dx(2)
          do i = lo(1), hi(1)
             x = prob_lo(1) + (dble(i) + HALF) * dx(1)

             cell_valid = .true.
             if (present(mask)) then
                if ( (.not. mask(i,j,k)) ) cell_valid = .false.
             endif

             if (cell_valid .and. s(i,j,k,rho_comp) > base_cutoff_density) then

                ! temperature diagnostic
                ! check to see if we are in the helium layer
                ! if we are, then get T_max
                if ( s(i,j,k,spec_comp) .ge. &
                     diag_define_he_layer * s(i,j,k,rho_comp) ) then

                   if (s(i,j,k,temp_comp) > T_max) then
                   
                      T_max = s(i,j,k,temp_comp)

                      coord_T_max(1) = x
                      coord_T_max(2) = y
                      coord_T_max(3) = z
                   
                   endif
                endif

                ! enuc diagnostic
                enuc_local = rho_Hnuc(i,j,k)/s(i,j,k,rho_comp)
                if (enuc_local > enuc_max) then
                
                   enuc_max = enuc_local

                   coord_enuc_max(1) = x
                   coord_enuc_max(2) = y
                   coord_enuc_max(3) = z

                endif

                ! c12 mass diagnostic
                c12_mass = c12_mass + weight*s(i,j,k,spec_comp+1)

                ! o16 mass diagnostic
                if (network_contains_oxygen) &
                     o16_mass = o16_mass + weight*s(i,j,k,spec_comp+2)

                ! vel diagnostic
                vel = sqrt(u(i,j,k,1)**2 + u(i,j,k,2)**2 + &
                     (u(i,j,k,3) + w0_cent)**2)
                if (vel > vel_max) then
                   
                   vel_max = vel
                   
                   coord_vel_max(1) = x
                   coord_vel_max(2) = y
                   coord_vel_max(3) = z

                endif

                ! Mach number diagnostic
                ! call the EOS to get the sound speed
                temp_eos(1) = s(i,j,k,temp_comp)
                den_eos(1)  = s(i,j,k,rho_comp)
                xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

                call eos(eos_input_rt, den_eos, temp_eos, &
                         npts, &
                         xn_eos, &
                         p_eos, h_eos, e_eos, &
                         cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                         dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                         dpdX_eos, dhdX_eos, &
                         gam1_eos, cs_eos, s_eos, &
                         dsdt_eos, dsdr_eos, &
                         .false.)

                Mach = vel/cs_eos(1)
                if (Mach > Machno_max) then

                   Machno_max = Mach

                   coord_Machno_max(1) = x
                   coord_Machno_max(2) = y
                   coord_Machno_max(3) = z

                endif

                ! deltap diagnostic
                if (do_deltap_diag) then
                   ! here we mirror what we do in make_tfromH
                   ! we already have the pressure from the above eos call
                   deltap = abs(p_eos(1) - p0(k)) / p0(k)

                   if (deltap > deltap_max) then
                   
                      deltap_max = deltap

                      coord_deltap_max(1) = x
                      coord_deltap_max(2) = y
                      coord_deltap_max(3) = z

                   endif

                   deltap_avg = deltap_avg + weight*deltap
                   nzones = nzones + weight
                   
                endif

             endif

          enddo
       enddo
    enddo

  end subroutine diag_3d

end module diag_module
