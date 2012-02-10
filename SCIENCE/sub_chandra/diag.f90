! sub-Chandra-specific diagnostic routine
! 
! currently, there are 3 output files (we mirror the files structure 
! from the wdconvect problem)
!
!   subchandra_enuc_diag.out:
!          peak nuc energy generation rate (erg / g / s)
!          x/y/z location of peak enuc
!          velocity components at location of peak enuc
!          radius of peak enuc
!
!   subchandra_temp_diag.out:          
!          peak temperature
!          x/y/z location of peak temperature
!          velocity components at location of peak temperature
!          radius of peak temperature
!
!   subchandra_vel_diag.out:
!          peak total velocity
!          peak Mach number (both in the star and full domain)
!          timestep
!
!
! We hold many timesteps-worth of diagnostic information in a buffer
! and output to the files only when flush_diag() is called.  This
! gives better performance on large machines with slow filesystems.
!

module diag_module

  use bl_types, only: dp_t

  implicit none

  private


  ! buffers
  real(kind=dp_t), allocatable, save :: time_data(:)
  real(kind=dp_t), allocatable, save :: file1_data(:,:), file2_data(:,:), &
                                        file3_data(:,:), file4_data(:,:)

  integer, save :: nstored = 0

  public :: diag, flush_diag

contains

  subroutine diag(time,dt,dx,s,rho_Hnuc,rho_Hext,thermal,rho_omegadot, &
                  rho0,rhoh0,p0,tempbar, &
                  gamma1bar,div_coeff, &
                  u,w0,normal, &
                  mla,the_bc_tower)

    use parallel
    use bl_constants_module, only: ZERO, FOUR, FOUR3RD, M_PI
    use bl_prof_module, only: bl_prof_timer, build
    use bl_error_module, only: bl_error

    use fundamental_constants_module, only: Gconst

    use fab_module, only: lwb, upb
    use multifab_module, only: multifab, multifab_build, multifab_build_edge, destroy, &
                               multifab_remote, dataptr, setval, get_box
    use ml_layout_module, only: ml_layout
    use define_bc_module, only: bc_tower

    use geometry, only: spherical, nr_fine, &
                        r_cc_loc, r_edge_loc, dr, center
    use variables, only: foextrap_comp
    use fill_3d_module, only: put_1d_array_on_cart, make_w0mac
    use probin_module, only: prob_lo_x, prob_lo_y, prob_lo_z, &
                             prob_hi_x, prob_hi_y, prob_hi_z, &
                             base_cutoff_density, &
                             diag_buf_size, octant

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
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: rhnp(:,:,:,:)
    real(kind=dp_t), pointer :: rhep(:,:,:,:)
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: nop(:,:,:,:)
    real(kind=dp_t), pointer :: w0rp(:,:,:,:)
    real(kind=dp_t), pointer :: w0xp(:,:,:,:)
    real(kind=dp_t), pointer :: w0yp(:,:,:,:)
    real(kind=dp_t), pointer :: w0zp(:,:,:,:)
    logical,         pointer :: mp(:,:,:,:)

    type(multifab) :: w0r_cart(mla%nlevel)
    type(multifab) ::    w0mac(mla%nlevel,mla%dim)

    real(kind=dp_t) :: T_max,     T_max_level,     T_max_local
    real(kind=dp_t) :: enuc_max,  enuc_max_level,  enuc_max_local

    real(kind=dp_t) :: U_max,           U_max_level,           U_max_local
    real(kind=dp_t) :: Mach_max,        Mach_max_level,        Mach_max_local
    real(kind=dp_t) :: Mach_max_domain, Mach_max_domain_level, Mach_max_domain_local

    real(kind=dp_t) :: coord_Tmax_local(mla%dim), coord_Tmax_level(mla%dim), coord_Tmax(mla%dim)
    real(kind=dp_t) :: Rloc_Tmax
    real(kind=dp_t) :: vel_Tmax_local(mla%dim), vel_Tmax_level(mla%dim), vel_Tmax(mla%dim)
    real(kind=dp_t) :: vr_Tmax
    
    real(kind=dp_t) :: coord_enucmax_local(mla%dim), coord_enucmax_level(mla%dim), coord_enucmax(mla%dim)
    real(kind=dp_t) :: Rloc_enucmax
    real(kind=dp_t) :: vel_enucmax_local(mla%dim), vel_enucmax_level(mla%dim), vel_enucmax(mla%dim)
    real(kind=dp_t) :: vr_enucmax


    ! buffers
    real(kind=dp_t) :: T_max_data_local(1), T_max_coords_local(2*mla%dim)
    real(kind=dp_t), allocatable :: T_max_data(:), T_max_coords(:)

    real(kind=dp_t) :: enuc_max_data_local(1), enuc_max_coords_local(2*mla%dim)
    real(kind=dp_t), allocatable :: enuc_max_data(:), enuc_max_coords(:)

    real(kind=dp_t) :: max_data_level(3), max_data_local(3)

    integer :: index

    integer :: index_max

    integer :: lo(mla%dim),hi(mla%dim)
    integer :: ng_s,ng_u,ng_n,ng_w,ng_wm,ng_rhn,ng_rhe
    integer :: i,n, comp, dm, nlevs

    type(bl_prof_timer), save :: bpt

    ! the maximum number of quantities to store in a size file -- for the 
    ! buffering
    integer, parameter :: MAX_FIELDS_PER_FILE = 32

    logical, save :: firstCall = .true.

    call build(bpt, "diagnostics")

    dm = mla%dim
    nlevs = mla%nlevel

    if (firstCall) then
       
       ! allocate the storage space for the buffers -- diag_buf_size
       ! is a runtime parameter that specifies how many steps we
       ! should go between outputting.  We need to make sure that we
       ! call flush_diag() before (or when) we reach diag_buf_size
       ! timesteps stored.
       allocate(time_data(diag_buf_size))
       allocate(file1_data(diag_buf_size, MAX_FIELDS_PER_FILE))
       allocate(file2_data(diag_buf_size, MAX_FIELDS_PER_FILE))
       allocate(file3_data(diag_buf_size, MAX_FIELDS_PER_FILE))
       allocate(file4_data(diag_buf_size, MAX_FIELDS_PER_FILE))

       nstored = 0
       time_data(:) = ZERO
       file1_data(:,:) = ZERO
       file2_data(:,:) = ZERO
       file3_data(:,:) = ZERO
       file4_data(:,:) = ZERO

       firstCall = .false.
    endif

       
    if (spherical .eq. 1) then

       do n=1,nlevs

          do comp=1,dm
             ! w0mac will contain an edge-centered w0 on a Cartesian grid,   
             ! for use in computing divergences.                            
             call multifab_build_edge(w0mac(n,comp), mla%la(n),1,1,comp)
             call setval(w0mac(n,comp), ZERO, all=.true.)
          enddo

          ! w0r_cart is w0 but onto a Cartesian grid in cell-centered as
          ! a scalar.  Since w0 is the radial expansion velocity, w0r_cart
          ! is the radial w0 in a zone
          call multifab_build(w0r_cart(n), mla%la(n),1,0)
          call setval(w0r_cart(n), ZERO, all=.true.)
       end do

       ! put w0 on Cartesian edges as a vector  
       call make_w0mac(mla,w0,w0mac,dx,the_bc_tower%bc_tower_array)


       ! put w0 in Cartesian cell-centers as a scalar (the radial 
       ! expansion velocity)
       call put_1d_array_on_cart(w0,w0r_cart,1,.true.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla)

    else

       call bl_error("ERROR: wdconvect/diag.f90: geometry not spherical")

    endif

    ng_s = s(1)%ng
    ng_u = u(1)%ng
    ng_n = normal(1)%ng
    ng_w = w0r_cart(1)%ng
    ng_wm = w0mac(1,1)%ng
    ng_rhn = rho_Hnuc(1)%ng
    ng_rhe = rho_Hext(1)%ng


    !=========================================================================
    ! initialize
    !=========================================================================
    T_max    = ZERO
    enuc_max = ZERO

    U_max           = ZERO
    Mach_max        = ZERO
    Mach_max_domain = ZERO

    coord_Tmax(:) = ZERO
    vel_Tmax(:) = ZERO

    coord_enucmax(:) = ZERO
    vel_enucmax(:) = ZERO


    !=========================================================================
    ! loop over the levels and compute the global quantities
    !=========================================================================
    do n = 1, nlevs

       ! initialize the local (processor's version) and level quantities to 0
       T_max_level = ZERO
       T_max_local = ZERO
       
       enuc_max_level = ZERO
       enuc_max_local = ZERO

       U_max_level = ZERO
       U_max_local = ZERO

       Mach_max_level = ZERO
       Mach_max_local = ZERO

       Mach_max_domain_level = ZERO
       Mach_max_domain_local = ZERO

       coord_Tmax_local(:) = ZERO
       coord_Tmax_level(:) = ZERO

       vel_Tmax_local(:) = ZERO
       vel_Tmax_level(:) = ZERO

       coord_enucmax_local(:) = ZERO
       coord_enucmax_level(:) = ZERO

       vel_enucmax_local(:) = ZERO
       vel_enucmax_level(:) = ZERO


       !----------------------------------------------------------------------
       ! loop over boxes in a given level
       !----------------------------------------------------------------------
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          sp => dataptr(s(n) , i)
          rhnp => dataptr(rho_Hnuc(n), i)
          rhep => dataptr(rho_Hext(n), i)
          up => dataptr(u(n) , i)
          nop => dataptr(normal(n) , i)
          w0rp => dataptr(w0r_cart(n), i)
          w0xp => dataptr(w0mac(n,1), i)
          w0yp => dataptr(w0mac(n,2), i)
          w0zp => dataptr(w0mac(n,3), i)

          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))

          select case (dm)
          case (2)
             call bl_error("ERROR: 2-d diag not implmented")
          case (3)
             if (n .eq. nlevs) then
                call diag_3d(n,time,dt,dx(n,:), &
                             sp(:,:,:,:),ng_s, &
                             rhnp(:,:,:,1), ng_rhn, &
                             rhep(:,:,:,1), ng_rhe, &
                             up(:,:,:,:),ng_u, &
                             w0rp(:,:,:,1), ng_w, &
                             w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1),ng_wm, & 
                             nop(:,:,:,:),ng_n, &
                             lo,hi, &
                             T_max_local, coord_Tmax_local, vel_Tmax_local, &
                             enuc_max_local, coord_enucmax_local, vel_enucmax_local, &
                             U_max_local, Mach_max_local, Mach_max_domain_local)

             else
                mp => dataptr(mla%mask(n), i)
                call diag_3d(n,time,dt,dx(n,:), &
                             sp(:,:,:,:),ng_s, &
                             rhnp(:,:,:,1), ng_rhn, &
                             rhep(:,:,:,1), ng_rhe, &
                             up(:,:,:,:),ng_u, &
                             w0rp(:,:,:,1), ng_w, &
                             w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1),ng_wm, & 
                             nop(:,:,:,:),ng_n, &
                             lo,hi, &
                             T_max_local, coord_Tmax_local, vel_Tmax_local, &
                             enuc_max_local, coord_enucmax_local, vel_enucmax_local, &
                             U_max_local, Mach_max_local, Mach_max_domain_local, &
                             mp(:,:,:,1))
             end if
          end select
       end do


       !----------------------------------------------------------------------
       ! do the appropriate parallel reduction for the current level
       !----------------------------------------------------------------------

       ! NOTE: only the I/O Processor will have the correct reduced value


       ! pack the quantities that we are taking the max of into a vector
       ! to reduce communication
       max_data_local(1) = U_max_local
       max_data_local(2) = Mach_max_local
       max_data_local(3) = Mach_max_domain_local

       call parallel_reduce(max_data_level, max_data_local, MPI_MAX, &
                            proc = parallel_IOProcessorNode())

       U_max_level           = max_data_level(1)
       Mach_max_level        = max_data_level(2)
       Mach_max_domain_level = max_data_level(3)


       ! for T_max, we want to know where the hot spot is, so we do a
       ! gather on the temperature and find the index corresponding to
       ! the maxiumum.  We then pack the coordinates and velocities
       ! into a local array and gather that to the I/O processor and
       ! pick the values corresponding to the maximum.
       allocate(T_max_data(parallel_nprocs()))
       T_max_data_local(1) = T_max_local

       call parallel_gather(T_max_data_local, T_max_data, 1, &
                            root = parallel_IOProcessorNode())


       index_max = maxloc(T_max_data, dim=1)
       
       ! T_max_coords will contain both the coordinate information and
       ! the velocity information, so there are 2*dm values on each
       ! proc
       allocate(T_max_coords(2*dm*parallel_nprocs()))
       T_max_coords_local(1) = coord_Tmax_local(1)
       T_max_coords_local(2) = coord_Tmax_local(2)
       T_max_coords_local(3) = coord_Tmax_local(3)
       T_max_coords_local(4) = vel_Tmax_local(1)
       T_max_coords_local(5) = vel_Tmax_local(2)
       T_max_coords_local(6) = vel_Tmax_local(3)
       
       call parallel_gather(T_max_coords_local, T_max_coords, 2*dm, &
                            root = parallel_IOProcessorNode())

       
       T_max_level = T_max_data(index_max)

       coord_Tmax_level(1) = T_max_coords(2*dm*(index_max-1)+1)
       coord_Tmax_level(2) = T_max_coords(2*dm*(index_max-1)+2)
       coord_Tmax_level(3) = T_max_coords(2*dm*(index_max-1)+3)
       vel_Tmax_level(1)   = T_max_coords(2*dm*(index_max-1)+4)
       vel_Tmax_level(2)   = T_max_coords(2*dm*(index_max-1)+5)
       vel_Tmax_level(3)   = T_max_coords(2*dm*(index_max-1)+6)


       deallocate(T_max_data)
       deallocate(T_max_coords)


       ! for enuc_max, we also want to know where the hot spot is, so
       ! we do the same gather procedure as with the temperature
       ! (above).
       allocate(enuc_max_data(parallel_nprocs()))
       enuc_max_data_local(1) = enuc_max_local

       call parallel_gather(enuc_max_data_local, enuc_max_data, 1, &
                            root = parallel_IOProcessorNode())


       index_max = maxloc(enuc_max_data, dim=1)
       
       ! enuc_max_coords will contain both the coordinate information
       ! and the velocity information, so there are 2*dm values on
       ! each proc
       allocate(enuc_max_coords(2*dm*parallel_nprocs()))
       enuc_max_coords_local(1) = coord_enucmax_local(1)
       enuc_max_coords_local(2) = coord_enucmax_local(2)
       enuc_max_coords_local(3) = coord_enucmax_local(3)
       enuc_max_coords_local(4) = vel_enucmax_local(1)
       enuc_max_coords_local(5) = vel_enucmax_local(2)
       enuc_max_coords_local(6) = vel_enucmax_local(3)
       
       call parallel_gather(enuc_max_coords_local, enuc_max_coords, 2*dm, &
                            root = parallel_IOProcessorNode())

       
       enuc_max_level = enuc_max_data(index_max)

       coord_enucmax_level(1) = enuc_max_coords(2*dm*(index_max-1)+1)
       coord_enucmax_level(2) = enuc_max_coords(2*dm*(index_max-1)+2)
       coord_enucmax_level(3) = enuc_max_coords(2*dm*(index_max-1)+3)
       vel_enucmax_level(1)   = enuc_max_coords(2*dm*(index_max-1)+4)
       vel_enucmax_level(2)   = enuc_max_coords(2*dm*(index_max-1)+5)
       vel_enucmax_level(3)   = enuc_max_coords(2*dm*(index_max-1)+6)


       deallocate(enuc_max_data)
       deallocate(enuc_max_coords)


       !----------------------------------------------------------------------
       ! reduce the current level's data with the global data
       !----------------------------------------------------------------------
       if (parallel_IOProcessor()) then
          U_max           = max(U_max,    U_max_level)
          Mach_max        = max(Mach_max, Mach_max_level)
          Mach_max_domain = max(Mach_max_domain, Mach_max_domain_level)
          
          ! if T_max_level is the new max, then copy the location as well
          if (T_max_level > T_max) then
             T_max = T_max_level

             coord_Tmax(:) = coord_Tmax_level(:)

             ! compute the radius of the bubble from the center
             Rloc_Tmax = sqrt( (coord_Tmax(1) - center(1))**2 + &
                               (coord_Tmax(2) - center(2))**2 + &
                               (coord_Tmax(3) - center(3))**2 )

             vel_Tmax(:) = vel_Tmax_level(:)

             ! use the coordinates of the hot spot and the velocity components
             ! to compute the radial velocity at the hotspot
             vr_Tmax = &
                  ((coord_Tmax(1) - center(1))/Rloc_Tmax)*vel_Tmax(1) + &
                  ((coord_Tmax(2) - center(2))/Rloc_Tmax)*vel_Tmax(2) + &
                  ((coord_Tmax(3) - center(3))/Rloc_Tmax)*vel_Tmax(3)

          endif

          ! if enuc_max_level is the new max, then copy the location as well
          if (enuc_max_level > enuc_max) then
             enuc_max = enuc_max_level

             coord_enucmax(:) = coord_enucmax_level(:)

             ! compute the radius of the bubble from the center
             Rloc_enucmax = sqrt( (coord_enucmax(1) - center(1))**2 + &
                                  (coord_enucmax(2) - center(2))**2 + &
                                  (coord_enucmax(3) - center(3))**2 )

             vel_enucmax(:) = vel_enucmax_level(:)

             ! use the coordinates of the hot spot and the velocity components
             ! to compute the radial velocity at the hotspot
             vr_enucmax = &
                  ((coord_enucmax(1) - center(1))/Rloc_enucmax)*vel_enucmax(1) + &
                  ((coord_enucmax(2) - center(2))/Rloc_enucmax)*vel_enucmax(2) + &
                  ((coord_enucmax(3) - center(3))/Rloc_enucmax)*vel_enucmax(3)


          endif

       endif

    end do


    !=========================================================================
    ! normalize
    !=========================================================================

    ! any normalization of integral quantities would go here


    !=========================================================================
    ! store the current step's data in the buffers
    !=========================================================================

    ! get the index into the buffer arrays for the current step's information.
    index = get_next_buffer_index()

    ! time
    time_data(index) = time

    ! for the file information, we need to coordinate with the header
    ! information printing out in flush_diag() to make sure that we are storing
    ! the right information in the right order.

    if (parallel_IOProcessor()) then
       ! file1 -- we are not storing a file1, to maintain consistency with
       ! the wdconvect diag.f90
       

       ! file2 -- subchandra_temp_diag.out
       file2_data(index, 1) = T_max
       file2_data(index, 2) = coord_Tmax(1)
       file2_data(index, 3) = coord_Tmax(2)
       file2_data(index, 4) = coord_Tmax(3)
       file2_data(index, 5) = vel_Tmax(1)
       file2_data(index, 6) = vel_Tmax(2)
       file2_data(index, 7) = vel_Tmax(3)
       file2_data(index, 8) = Rloc_Tmax
       file2_data(index, 9) = vr_Tmax


       ! file3 -- subchandra_enuc_diag.out
       file3_data(index, 1) = enuc_max
       file3_data(index, 2) = coord_enucmax(1)
       file3_data(index, 3) = coord_enucmax(2)
       file3_data(index, 4) = coord_enucmax(3)
       file3_data(index, 5) = vel_enucmax(1)
       file3_data(index, 6) = vel_enucmax(2)
       file3_data(index, 7) = vel_enucmax(3)
       file3_data(index, 8) = Rloc_enucmax
       file3_data(index, 9) = vr_enucmax


       ! file4 -- subchandra_vel_diag.out
       file4_data(index, 1) = U_max
       file4_data(index, 2) = Mach_max
       file4_data(index, 3) = Mach_max_domain
       file4_data(index, 4) = dt
    endif


    !=========================================================================
    ! output, if needed
    !=========================================================================

    ! if we've filled the buffers, flush them
    if (index == diag_buf_size) then
       call flush_diag()
    endif


    !=========================================================================
    ! clean-up
    !=========================================================================
    if (spherical .eq. 1) then
       do n=1,nlevs
          call destroy(w0r_cart(n))
          do comp=1,dm
             call destroy(w0mac(n,comp))
          enddo
       end do
    end if

    call destroy(bpt)

  end subroutine diag



  !===========================================================================
  ! flush_diag -- the output routine.  When this routine is called, it 
  ! outputs all the stored information in the buffers and resets them.
  !===========================================================================
  subroutine flush_diag()

    use parallel
    use bl_constants_module, only: ZERO
    use bl_IO_module, only: unit_new
    use bl_system_module, only: BL_CWD_SIZE, get_cwd 
    use probin_module, only: job_name

    integer :: un1,un2,un3,un4
    logical :: lexist

    character (len=16) :: date_str, time_str
    integer, dimension(8) :: values
    character (len=BL_CWD_SIZE) :: cwd

    integer :: n

    logical, save :: firstCall = .true.


    ! if the buffers are empty, move on
    if (nstored == 0) return


    ! IMPORTANT: make sure that there are enough entries in the format
    ! statement to write out all of the data in each file.
999 format("# job name: ",a)
1000 format(1x,16(g20.10,1x))
1001 format("#",16(a20,1x))
800 format("# ",a,i4.4,'-',i2.2,'-',i2.2)
801 format("# ",a,i2.2,':',i2.2,':',i2.2)
802 format("# ",a,a)

    if (parallel_IOProcessor()) then

       ! open the diagnostic files for output, taking care not to overwrite
       ! an existing file
       !un1 = unit_new()
       !inquire(file="subchandra_radvel_diag.out", exist=lexist)
       !if (lexist) then
       !   open(unit=un1, file="subchandra_radvel_diag.out", &
       !        status="old", position="append")
       !else
       !   open(unit=un1, file="subchandra_radvel_diag.out", status="new")
       !endif

       un2 = unit_new()
       inquire(file="subchandra_temp_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un2, file="subchandra_temp_diag.out", &
               status="old", position="append")
       else
          open(unit=un2, file="subchandra_temp_diag.out", status="new")
       endif

       un3 = unit_new()
       inquire(file="subchandra_enuc_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un3, file="subchandra_enuc_diag.out", &
               status="old", position="append")
       else
          open(unit=un3, file="subchandra_enuc_diag.out", status="new")
       endif

       un4 = unit_new()
       inquire(file="subchandra_vel_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un4, file="subchandra_vel_diag.out", &
               status="old", position="append")
       else
          open(unit=un4, file="subchandra_vel_diag.out", status="new")
       endif


       ! write out the headers
       if (firstCall) then

          ! get the data and time
          call date_and_time(date_str, time_str, VALUES=values)

          ! get the output directory
          call get_cwd(cwd)

          ! temp
          write (un2, *) " "
          write (un2, 800) "output date: ", values(1), values(2), values(3)
          write (un2, 801) "output time: ", values(5), values(6), values(7)
          write (un2, 802) "output dir:  ", trim(cwd)
          write (un2, 999) trim(job_name)
          write (un2,1001) "time", "max{T}", "x(max{T})", "y(max{T})", "z(max{T})", &
                           "vx(max{T})", "vy(max{T})", "vz(max{T})", &
                           "R(max{T})", "vr(max{T})"

          ! enuc
          write (un3, *) " "
          write (un3, 800) "output date: ", values(1), values(2), values(3)
          write (un3, 801) "output time: ", values(5), values(6), values(7)
          write (un3, 802) "output dir:  ", trim(cwd)
          write (un3, 999) trim(job_name)
          write (un3,1001) "time", "max{enuc}", &
                           "x(max{enuc})", "y(max{enuc})", "z(max{enuc})", &
                           "vx(max{enuc})", "vy(max{enuc})", "vz(max{enuc})", &
                           "R(max{enuc})", "vr(max{enuc})"

          ! vel
          write (un4, *) " "
          write (un4, 800) "output date: ", values(1), values(2), values(3)
          write (un4, 801) "output time: ", values(5), values(6), values(7)
          write (un4, 802) "output dir:  ", trim(cwd)
          write (un4, 999) trim(job_name)
          write (un4,1001) "time", "max{|U + w0|}", &
                           "max{Mach #}", "max{Mach #} (domain)", &
                           "dt"

          firstCall = .false.
       endif

       do n = 1, nstored

          ! no data to un1

          write (un2,1000) time_data(n), &
               file2_data(n,1), &
               file2_data(n,2), file2_data(n,3), file2_data(n,4), &
               file2_data(n,5), file2_data(n,6), file2_data(n,7), &
               file2_data(n,8), file2_data(n,9)

          write (un3,1000) time_data(n), &
               file3_data(n,1), &
               file3_data(n,2), file3_data(n,3), file3_data(n,4), &
               file3_data(n,5), file3_data(n,6), file3_data(n,7), &
               file3_data(n,8), file3_data(n,9)

          write (un4,1000) time_data(n), &
               file4_data(n,1), file4_data(n,2), file4_data(n,3), file4_data(n,4)
          
       enddo

       !close(un1)
       close(un2)
       close(un3)
       close(un4)

    endif


    ! reset the buffers
    nstored = 0
    time_data(:) = ZERO
    file1_data(:,:) = ZERO
    file2_data(:,:) = ZERO
    file3_data(:,:) = ZERO
    file4_data(:,:) = ZERO

  end subroutine flush_diag



  !===========================================================================
  ! get_next_buffer_index returns the next index into the main buffers for
  ! storing one timestep's data.  It increments the nstored index to account
  ! for the new data
  !===========================================================================
  function get_next_buffer_index() result (next_index)
    integer :: next_index

    ! we will use 1-based indexing
    nstored = nstored + 1
    next_index = nstored

  end function get_next_buffer_index



  !===========================================================================
  ! the actual diagnostic routine
  !===========================================================================
  subroutine diag_3d(n,time,dt,dx, &
                     s,ng_s, &
                     rho_Hnuc,ng_rhn, &
                     rho_Hext,ng_rhe, &
                     u,ng_u, &
                     w0r,ng_w, &
                     w0macx,w0macy,w0macz,ng_wm, &
                     normal,ng_n, &
                     lo,hi, &
                     T_max,coord_Tmax,vel_Tmax, &
                     enuc_max,coord_enucmax,vel_enucmax, &
                     U_max,Mach_max,Mach_max_domain, &
                     mask)

    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp
    use bl_constants_module, only: HALF
    use bl_error_module, only: bl_error
    use network, only: nspec
    use geometry, only: spherical, center
    use probin_module, only: base_cutoff_density, prob_lo, sponge_start_factor, &
         sponge_center_density
    use eos_module

    integer,          intent(in   ) :: n,lo(:),hi(:),ng_s,ng_u,ng_n,ng_w,ng_wm,ng_rhn,ng_rhe
    real (kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:  ,lo(2)-ng_s:  ,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho_Hnuc(lo(1)-ng_rhn:,lo(2)-ng_rhn:,lo(3)-ng_rhn:)
    real (kind=dp_t), intent(in   ) :: rho_Hext(lo(1)-ng_rhe:,lo(2)-ng_rhe:,lo(3)-ng_rhe:)
    real (kind=dp_t), intent(in   ) ::        u(lo(1)-ng_u:  ,lo(2)-ng_u:  ,lo(3)-ng_u:,:)
    real (kind=dp_t), intent(in   ) ::      w0r(lo(1)-ng_w:  ,lo(2)-ng_w:  ,lo(3)-ng_w:)
    real (kind=dp_t), intent(in   ) ::   w0macx(lo(1)-ng_wm: ,lo(2)-ng_wm: ,lo(3)-ng_wm:)
    real (kind=dp_t), intent(in   ) ::   w0macy(lo(1)-ng_wm: ,lo(2)-ng_wm: ,lo(3)-ng_wm:)
    real (kind=dp_t), intent(in   ) ::   w0macz(lo(1)-ng_wm: ,lo(2)-ng_wm: ,lo(3)-ng_wm:)
    real (kind=dp_t), intent(in   ) ::   normal(lo(1)-ng_n:  ,lo(2)-ng_n:  ,lo(3)-ng_n:,:)
    real (kind=dp_t), intent(in   ) :: time, dt, dx(:)
    real (kind=dp_t), intent(inout) :: T_max, coord_Tmax(:)
    real (kind=dp_t), intent(inout) :: vel_Tmax(:)
    real (kind=dp_t), intent(inout) :: enuc_max, coord_enucmax(:)
    real (kind=dp_t), intent(inout) :: vel_enucmax(:)
    real (kind=dp_t), intent(inout) :: U_max, Mach_max, Mach_max_domain
    logical,          intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    !     Local variables
    integer            :: i, j, k
    real (kind=dp_t)   :: velr, vel, weight
    logical            :: cell_valid
    real (kind=dp_t)   :: x, y, z

    ! weight is the factor by which the volume of a cell at the
    ! current level relates to the volume of a cell at the coarsest
    ! level of refinement.
    weight = 1.d0 / 8.d0**(n-1)

    if (.not. spherical == 1) then
       call bl_error("ERROR: geometry not spherical in diag")
    endif

!$omp parallel do private(i,j,k,x,y,z,cell_valid,velr,vel) &
!$omp reduction(max:vr_max,U_max,Mach_max,Mach_max_domain) &
    do k = lo(3), hi(3)
       z = prob_lo(3) + (dble(k)+HALF) * dx(3)

       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j)+HALF) * dx(2)

          do i = lo(1), hi(1) 
             x = prob_lo(1) + (dble(i)+HALF) * dx(1)
                
             cell_valid = .true.
             if (present(mask)) then
                if ( (.not. mask(i,j,k)) ) cell_valid = .false.
             end if

             ! we only consider cells inside of where the sponging begins
             if (cell_valid) then

                ! vel is the magnitude of the velocity, including w0
                vel = sqrt( (u(i,j,k,1)+HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)))**2 + &
                            (u(i,j,k,2)+HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)))**2 + &
                            (u(i,j,k,3)+HALF*(w0macz(i,j,k)+w0macz(i,j,k+1)))**2)

                if (s(i,j,k,rho_comp) >= sponge_start_factor*sponge_center_density) then

!$omp critical
                   ! max T, location, and velocity at that location (including w0)
                   if (s(i,j,k,temp_comp) > T_max) then
                      T_max = s(i,j,k,temp_comp)
                      coord_Tmax(1) = x
                      coord_Tmax(2) = y
                      coord_Tmax(3) = z
                      vel_Tmax(1) = u(i,j,k,1)+HALF*(w0macx(i,j,k)+w0macx(i+1,j,k))
                      vel_Tmax(2) = u(i,j,k,2)+HALF*(w0macy(i,j,k)+w0macy(i,j+1,k))
                      vel_Tmax(3) = u(i,j,k,3)+HALF*(w0macz(i,j,k)+w0macz(i,j,k+1))
                   endif

                   ! max enuc
                   if (rho_Hnuc(i,j,k)/s(i,j,k,rho_comp) > enuc_max) then
                      enuc_max = rho_Hnuc(i,j,k)/s(i,j,k,rho_comp)
                      coord_enucmax(1) = x
                      coord_enucmax(2) = y
                      coord_enucmax(3) = z
                      vel_enucmax(1) = u(i,j,k,1)+HALF*(w0macx(i,j,k)+w0macx(i+1,j,k))
                      vel_enucmax(2) = u(i,j,k,2)+HALF*(w0macy(i,j,k)+w0macy(i,j+1,k))
                      vel_enucmax(3) = u(i,j,k,3)+HALF*(w0macz(i,j,k)+w0macz(i,j,k+1))
                   endif
!$omp end critical

                endif   ! end density check

                ! call the EOS to get the sound speed and internal energy
                temp_eos = s(i,j,k,temp_comp)
                den_eos  = s(i,j,k,rho_comp)
                xn_eos(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos

                call eos(eos_input_rt, den_eos, temp_eos, &
                         xn_eos, &
                         p_eos, h_eos, e_eos, &
                         cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                         dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                         dpdX_eos, dhdX_eos, &
                         gam1_eos, cs_eos, s_eos, &
                         dsdt_eos, dsdr_eos, &
                         .false.)

                ! max vel and Mach number
                U_max = max(U_max,vel)

                if (s(i,j,k,rho_comp) >= sponge_start_factor*sponge_center_density) then
                   Mach_max = max(Mach_max,vel/cs_eos)
                endif

                Mach_max_domain = max(Mach_max_domain,vel/cs_eos)

             endif  ! end cell_valid and density check

          enddo
       enddo
    enddo

  end subroutine diag_3d

end module diag_module
