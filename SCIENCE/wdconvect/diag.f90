! wdconvect-specific diagnostic routine
! 
! currently, there are 4 output files:
!
!   wdconvect_enuc_diag.out:
!          peak nuc energy generation rate (erg / g / s)
!          x/y/z location of peak enuc
!          velocity components at location of peak enuc
!          total nuclear energy release (erg / s)
!          radius of peak enuc
!
!   wdconvect_radvel_diag.out:
!          radial velocity components (std. average & Favre average)
!          peak radial velocity
!          total mass
!
!   wdconvect_temp_diag.out:          
!          peak temperature
!          x/y/z location of peak temperature
!          velocity components at location of peak temperature
!          radius of peak temperature
!
!   wdconvect_vel_diag.out:
!          peak total velocity
!          peak Mach number
!          total kinetic energy
!          gravitational potential energy
!          total internal energy
!

module diag_module

  implicit none

  private

  public :: diag

contains

  subroutine diag(time,dt,dx,s,rho_Hnuc,rho_Hext,thermal,rho_omegadot, &
                  rho0,rhoh0,p0,tempbar, &
                  gamma1bar,div_coeff, &
                  u,w0,normal, &
                  mla,the_bc_tower)

    use parallel                        

    use bl_types, only: dp_t
    use bl_IO_module, only: unit_new
    use bl_constants_module, only: ZERO, FOUR, FOUR3RD, M_PI
    use bl_prof_module, only: bl_prof_timer, build
    use bl_error_module, only: bl_error
    use bl_system_module, only: BL_CWD_SIZE, get_cwd 

    use fundamental_constants_module, only: Gconst

    use fab_module, only: lwb, upb
    use multifab_module, only: multifab, multifab_build, destroy, &
                               multifab_remote, dataptr, setval, get_box
    use ml_layout_module, only: ml_layout
    use define_bc_module, only: bc_tower

    use geometry, only: dm, nlevs, spherical, nr_fine, &
                        r_cc_loc, r_edge_loc, dr, center
    use variables, only: foextrap_comp
    use fill_3d_module, only: put_1d_array_on_cart, make_w0mac
    use probin_module, only: prob_lo_x, prob_lo_y, prob_lo_z, &
                             prob_hi_x, prob_hi_y, prob_hi_z, &
                             job_name, &
                             edge_nodal_flag, &
                             base_cutoff_density

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
    real(kind=dp_t), pointer :: np(:,:,:,:)
    real(kind=dp_t), pointer :: w0rp(:,:,:,:)
    real(kind=dp_t), pointer :: w0xp(:,:,:,:)
    real(kind=dp_t), pointer :: w0yp(:,:,:,:)
    real(kind=dp_t), pointer :: w0zp(:,:,:,:)
    logical,         pointer :: mp(:,:,:,:)

    type(multifab) :: w0r_cart(mla%nlevel)
    type(multifab) ::    w0mac(mla%nlevel,dm)

    real(kind=dp_t) :: vr(dm),    vr_level(dm),    vr_local(dm)
    real(kind=dp_t) :: vr_max,    vr_max_level,    vr_max_local
    real(kind=dp_t) :: rhovr(dm), rhovr_level(dm), rhovr_local(dm)

    real(kind=dp_t) :: mass,      mass_level,      mass_local
    real(kind=dp_t) :: nzones,    nzones_level,    nzones_local

    real(kind=dp_t) :: T_max,     T_max_level,     T_max_local
    real(kind=dp_t) :: enuc_max,  enuc_max_level,  enuc_max_local

    real(kind=dp_t) :: kin_ener,  kin_ener_level,  kin_ener_local
    real(kind=dp_t) :: int_ener,  int_ener_level,  int_ener_local
    real(kind=dp_t) :: nuc_ener,  nuc_ener_level,  nuc_ener_local

    real(kind=dp_t) :: U_max,     U_max_level,     U_max_local
    real(kind=dp_t) :: Mach_max,  Mach_max_level,  Mach_max_local

    real(kind=dp_t) :: coord_Tmax_local(dm), coord_Tmax_level(dm), coord_Tmax(dm)
    real(kind=dp_t) :: Rloc_Tmax
    real(kind=dp_t) :: vel_Tmax_local(dm), vel_Tmax_level(dm), vel_Tmax(dm)
    real(kind=dp_t) :: vr_Tmax
    
    real(kind=dp_t) :: coord_enucmax_local(dm), coord_enucmax_level(dm), coord_enucmax(dm)
    real(kind=dp_t) :: Rloc_enucmax
    real(kind=dp_t) :: vel_enucmax_local(dm), vel_enucmax_level(dm), vel_enucmax(dm)
    real(kind=dp_t) :: vr_enucmax

    real(kind=dp_t) :: vel_center_local(dm), vel_center_level(dm), vel_center(dm)
    real(kind=dp_t) :: T_center_local, T_center_level, T_center

    integer         :: ncenter_local, ncenter_level, ncenter

    ! buffers
    real(kind=dp_t) :: T_max_data_local(1), T_max_coords_local(2*dm)
    real(kind=dp_t), allocatable :: T_max_data(:), T_max_coords(:)

    real(kind=dp_t) :: enuc_max_data_local(1), enuc_max_coords_local(2*dm)
    real(kind=dp_t), allocatable :: enuc_max_data(:), enuc_max_coords(:)

    real(kind=dp_t) :: max_data_level(3), max_data_local(3)
    real(kind=dp_t) :: sum_data_level(2*dm+9), sum_data_local(2*dm+9)


    integer :: index_max

    real(kind=dp_t) :: vr_favre(dm)

    real(kind=dp_t) :: grav_ener, term1, term2
    real(kind=dp_t), allocatable :: m(:)


    integer :: lo(dm),hi(dm)
    integer :: ng_s,ng_u,ng_n,ng_w,ng_wm,ng_rhn,ng_rhe
    integer :: i,n, comp, r
    integer :: un1,un2,un3,un4
    logical :: lexist

    character (len=16) :: date_str, time_str
    integer, dimension(8) :: values
    character (len=BL_CWD_SIZE) :: cwd

    logical, save :: firstCall = .true.

    type(bl_prof_timer), save :: bpt

    call build(bpt, "diagnostics")

    if (spherical .eq. 1) then

       do n=1,nlevs

          do comp=1,dm
             ! w0mac will contain an edge-centered w0 on a Cartesian grid,   
             ! for use in computing divergences.                            
             call multifab_build(w0mac(n,comp), mla%la(n),1,1, &
                                 nodal=edge_nodal_flag(comp,:))
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
       call put_1d_array_on_cart(w0,w0r_cart,foextrap_comp,.true.,.false.,dx, &
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
    vr(:)    = ZERO
    rhovr(:) = ZERO

    mass     = ZERO
    nzones   = ZERO

    vr_max   = ZERO

    T_max    = ZERO
    enuc_max = ZERO

    kin_ener = ZERO
    int_ener = ZERO
    nuc_ener = ZERO

    U_max    = ZERO
    Mach_max  = ZERO

    coord_Tmax(:) = ZERO
    vel_Tmax(:) = ZERO

    coord_enucmax(:) = ZERO
    vel_enucmax(:) = ZERO

    vel_center(:) = ZERO
    T_center = ZERO
    
    ncenter = 0


    !=========================================================================
    ! loop over the levels and compute the global quantities
    !=========================================================================
    do n = 1, nlevs

       ! initialize the local (processor's version) and level quantities to 0
       vr_level(:) = ZERO
       vr_local(:) = ZERO
       
       rhovr_level(:) = ZERO
       rhovr_local(:) = ZERO
       
       mass_level = ZERO
       mass_local = ZERO
       
       nzones_level = ZERO
       nzones_local = ZERO
       
       vr_max_level = ZERO
       vr_max_local = ZERO
       
       T_max_level = ZERO
       T_max_local = ZERO
       
       enuc_max_level = ZERO
       enuc_max_local = ZERO

       kin_ener_level = ZERO
       kin_ener_local = ZERO

       int_ener_level = ZERO
       int_ener_local = ZERO

       nuc_ener_level = ZERO
       nuc_ener_local = ZERO

       U_max_level = ZERO
       U_max_local = ZERO

       Mach_max_level = ZERO
       Mach_max_local = ZERO

       coord_Tmax_local(:) = ZERO
       coord_Tmax_level(:) = ZERO

       vel_Tmax_local(:) = ZERO
       vel_Tmax_level(:) = ZERO

       coord_enucmax_local(:) = ZERO
       coord_enucmax_level(:) = ZERO

       vel_enucmax_local(:) = ZERO
       vel_enucmax_level(:) = ZERO

       vel_center_local(:) = ZERO
       vel_center_level(:) = ZERO
       
       T_center_local = ZERO
       T_center_level = ZERO

       ncenter_local = 0
       ncenter_level = 0


       !----------------------------------------------------------------------
       ! loop over boxes in a given level
       !----------------------------------------------------------------------
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          sp => dataptr(s(n) , i)
          rhnp => dataptr(rho_Hnuc(n), i)
          rhep => dataptr(rho_Hext(n), i)
          up => dataptr(u(n) , i)
          np => dataptr(normal(n) , i)
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
                             np(:,:,:,:),ng_n, &
                             lo,hi, &
                             nzones_local, &
                             vr_local(1),vr_local(2),vr_local(3),vr_max_local, &
                             rhovr_local(1), rhovr_local(2), rhovr_local(3), mass_local, &
                             T_max_local, coord_Tmax_local, vel_Tmax_local, &
                             enuc_max_local, coord_enucmax_local, vel_enucmax_local, &
                             kin_ener_local, int_ener_local, nuc_ener_local, &
                             U_max_local, Mach_max_local, &
                             ncenter_local,T_center_local, &
                             vel_center_local(1),vel_center_local(2),vel_center_local(3))
             else
                mp => dataptr(mla%mask(n), i)
                call diag_3d(n,time,dt,dx(n,:), &
                             sp(:,:,:,:),ng_s, &
                             rhnp(:,:,:,1), ng_rhn, &
                             rhep(:,:,:,1), ng_rhe, &
                             up(:,:,:,:),ng_u, &
                             w0rp(:,:,:,1), ng_w, &
                             w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1),ng_wm, & 
                             np(:,:,:,:),ng_n, &
                             lo,hi, &
                             nzones_local, &
                             vr_local(1),vr_local(2),vr_local(3),vr_max_local, &
                             rhovr_local(1), rhovr_local(2), rhovr_local(3), mass_local, &
                             T_max_local, coord_Tmax_local, vel_Tmax_local, &
                             enuc_max_local, coord_enucmax_local, vel_enucmax_local, &
                             kin_ener_local, int_ener_local, nuc_ener_local, &
                             U_max_local, Mach_max_local, &
                             ncenter_local,T_center_local, &
                             vel_center_local(1),vel_center_local(2),vel_center_local(3), &
                             mp(:,:,:,1))
             end if
          end select
       end do

       !----------------------------------------------------------------------
       ! do the appropriate parallel reduction for the current level
       !----------------------------------------------------------------------

       ! NOTE: only the I/O Processor will have the correct reduced value

       ! pack the quantities that we are summing into a vector to reduce
       ! communication

       ! start with the real quantities...
       sum_data_local(1:dm)      = vr_local(:)
       sum_data_local(dm+1:2*dm) = rhovr_local(:)
       sum_data_local(2*dm+1)    = mass_local
       sum_data_local(2*dm+2)    = nzones_local
       sum_data_local(2*dm+3)    = kin_ener_local
       sum_data_local(2*dm+4)    = int_ener_local
       sum_data_local(2*dm+5)    = nuc_ener_local
       sum_data_local(2*dm+6)    = vel_center_local(1)
       sum_data_local(2*dm+7)    = vel_center_local(2)
       sum_data_local(2*dm+8)    = vel_center_local(3)
       sum_data_local(2*dm+9)    = T_center_local

       call parallel_reduce(sum_data_level, sum_data_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       vr_level(:)         = sum_data_level(1:dm)
       rhovr_level(:)      = sum_data_level(dm+1:2*dm)
       mass_level          = sum_data_level(2*dm+1)
       nzones_level        = sum_data_level(2*dm+2)
       kin_ener_level      = sum_data_level(2*dm+3)
       int_ener_level      = sum_data_level(2*dm+4)
       nuc_ener_level      = sum_data_level(2*dm+5)
       vel_center_level(1) = sum_data_level(2*dm+6)   
       vel_center_level(2) = sum_data_level(2*dm+7)    
       vel_center_level(3) = sum_data_level(2*dm+8)   
       T_center_level      = sum_data_level(2*dm+9)   

       ! ...and the integer quantities
       call parallel_reduce(ncenter_level, ncenter_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())


       ! pack the quantities that we are taking the max of into a vector
       ! to reduce communication
       max_data_local(1) = vr_max_local
       max_data_local(2) = U_max_local
       max_data_local(3) = Mach_max_local

       call parallel_reduce(max_data_level, max_data_local, MPI_MAX, &
                            proc = parallel_IOProcessorNode())

       vr_max_level   = max_data_level(1)
       U_max_level    = max_data_level(2)
       Mach_max_level = max_data_level(3)


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
          vr       = vr     + vr_level
          rhovr    = rhovr  + rhovr_level

          mass     = mass   + mass_level
          nzones   = nzones + nzones_level

          vr_max   = max(vr_max,   vr_max_level)

          kin_ener = kin_ener + kin_ener_level
          int_ener = int_ener + int_ener_level
          nuc_ener = nuc_ener + nuc_ener_level

          U_max    = max(U_max,    U_max_level)
          Mach_max = max(Mach_max, Mach_max_level)
          
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

          T_center = T_center + T_center_level

          vel_center = vel_center + vel_center_level


          ncenter = ncenter + ncenter_level

       endif

    end do

    !-------------------------------------------------------------------------
    ! compute the gravitational potential energy too.
    !-------------------------------------------------------------------------
    allocate(m(0:nr_fine-1))
    grav_ener = ZERO

    ! m(r) will contain mass enclosed by the center
    m(0) = FOUR3RD*M_PI*rho0(1,0)*r_cc_loc(1,0)**3

    ! dU = - G M dM / r;  dM = 4 pi r**2 rho dr  -->  dU = - 4 pi G r rho dr
    grav_ener = -FOUR*M_PI*Gconst*m(0)*r_cc_loc(1,0)*rho0(1,0)*dr(1)

    do r=1,nr_fine-1

       ! the mass is defined at the cell-centers, so to compute the
       ! mass at the current center, we need to add the contribution
       ! of the upper half of the zone below us and the lower half of
       ! the current zone.
       
       ! don't add any contributions from outside the star -- i.e.
       ! rho < base_cutoff_density
       if (rho0(1,r-1) > base_cutoff_density) then
          term1 = FOUR3RD*M_PI*rho0(1,r-1) * &
               (r_edge_loc(1,r) - r_cc_loc(1,r-1)) * &
               (r_edge_loc(1,r)**2 + &
                r_edge_loc(1,r)*r_cc_loc(1,r-1) + &
                r_cc_loc(1,r-1)**2)
       else
          term1 = ZERO
       endif

       if (rho0(1,r) > base_cutoff_density) then
          term2 = FOUR3RD*M_PI*rho0(1,r  )*&
               (r_cc_loc(1,r) - r_edge_loc(1,r  )) * &
               (r_cc_loc(1,r)**2 + &
                r_cc_loc(1,r)*r_edge_loc(1,r  ) + &
                r_edge_loc(1,r  )**2)          
       else
          term2 = ZERO
       endif

       m(r) = m(r-1) + term1 + term2
          
       ! dU = - G M dM / r;  
       ! dM = 4 pi r**2 rho dr  -->  dU = - 4 pi G r rho dr
       grav_ener = grav_ener - &
            FOUR*M_PI*Gconst*m(r)*r_cc_loc(1,r)*rho0(1,r)*dr(1)

    enddo


    

    !=========================================================================
    ! normalize
    !=========================================================================
    if (parallel_IOProcessor()) then

       vr(:) = vr(:)/nzones
       vr_favre(:) = rhovr(:)/mass    ! note: common dV normalization cancels

       ! the volume we normalize with is that of a single coarse-level
       ! zone.  This is because the weight used in the loop over cells
       ! was with reference to the coarse level

       mass = mass*dx(1,1)*dx(1,2)*dx(1,3)
       kin_ener = kin_ener*dx(1,1)*dx(1,2)*dx(1,3)
       int_ener = int_ener*dx(1,1)*dx(1,2)*dx(1,3)
       nuc_ener = nuc_ener*dx(1,1)*dx(1,2)*dx(1,3)

       
       ! ncenter should be 8 -- there are only 8 zones that have a
       ! vertex at the center of the star
       if (ncenter /= 8) then
          call bl_error("ERROR: ncenter /= 8 in diag")
       else
          T_center = T_center/ncenter
          vel_center(:) = vel_center(:)/ncenter
       endif
    endif


    !=========================================================================
    ! output
    !=========================================================================
 999 format("# job name: ",a)
1000 format(1x,10(g20.10,1x))
1001 format("#",10(a20,1x))
 800 format("# ",a,i4.4,'-',i2.2,'-',i2.2)
 801 format("# ",a,i2.2,':',i2.2,':',i2.2)
 802 format("# ",a,a)

    if (parallel_IOProcessor()) then

       ! open the diagnostic files for output, taking care not to overwrite
       ! an existing file
       un1 = unit_new()
       inquire(file="wdconvect_radvel_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un1, file="wdconvect_radvel_diag.out", &
               status="old", position="append")
       else
          open(unit=un1, file="wdconvect_radvel_diag.out", status="new")
       endif

       un2 = unit_new()
       inquire(file="wdconvect_temp_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un2, file="wdconvect_temp_diag.out", &
               status="old", position="append")
       else
          open(unit=un2, file="wdconvect_temp_diag.out", status="new")
       endif

       un3 = unit_new()
       inquire(file="wdconvect_enuc_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un3, file="wdconvect_enuc_diag.out", &
               status="old", position="append")
       else
          open(unit=un3, file="wdconvect_enuc_diag.out", status="new")
       endif

       un4 = unit_new()
       inquire(file="wdconvect_vel_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un4, file="wdconvect_vel_diag.out", &
               status="old", position="append")
       else
          open(unit=un4, file="wdconvect_vel_diag.out", status="new")
       endif


       ! write out the headers
       if (firstCall) then

          ! get the data and time
          call date_and_time(date_str, time_str, VALUES=values)

          ! get the output directory
          call get_cwd(cwd)

          ! radvel
          write (un1, *) " "
          write (un1, 800) "output date: ", values(1), values(2), values(3)
          write (un1, 801) "output time: ", values(5), values(6), values(7)
          write (un1, 802) "output dir:  ", trim(cwd)
          write (un1, 999) trim(job_name)
          write (un1, 1001) "time", "<vr_x>", "<vr_y>", "<vr_z>", "<vr>", &
                            "max{|vr|}", &
                            "int{rhovr_x}/mass", "int{rhovr_y}/mass", "int{rhovr_z}/mass", &
                            "mass"

          ! temp
          write (un2, *) " "
          write (un2, 800) "output date: ", values(1), values(2), values(3)
          write (un2, 801) "output time: ", values(5), values(6), values(7)
          write (un2, 802) "output dir:  ", trim(cwd)
          write (un2, 999) trim(job_name)
          write (un2,1001) "time", "max{T}", "x(max{T})", "y(max{T})", "z(max{T})", &
                           "vx(max{T})", "vy(max{T})", "vz(max{T})", &
                           "R(max{T})", "vr(max{T})", "T_center"

          ! enuc
          write (un3, *) " "
          write (un3, 800) "output date: ", values(1), values(2), values(3)
          write (un3, 801) "output time: ", values(5), values(6), values(7)
          write (un3, 802) "output dir:  ", trim(cwd)
          write (un3, 999) trim(job_name)
          write (un3,1001) "time", "max{enuc}", &
                           "x(max{enuc})", "y(max{enuc})", "z(max{enuc})", &
                           "vx(max{enuc})", "vy(max{enuc})", "vz(max{enuc})", &
                           "R(max{enuc})", "vr(max{enuc})", "tot nuc ener (erg/s)"

          ! vel
          write (un4, *) " "
          write (un4, 800) "output date: ", values(1), values(2), values(3)
          write (un4, 801) "output time: ", values(5), values(6), values(7)
          write (un4, 802) "output dir:  ", trim(cwd)
          write (un4, 999) trim(job_name)
          write (un4,1001) "time", "max{|U + w0|}", "max{Mach #}", &
                           "tot kin energy", "grav pot energy", "tot int energy", &
                           "velx_center", "vely_center", "velz_center", "dt"

          firstCall = .false.
       endif

       ! write out the data
       write (un1,1000) time, vr(1), vr(2), vr(3), &
            sqrt(vr(1)**2 + vr(2)**2 + vr(3)**2), vr_max, &
            vr_favre(1), vr_favre(2), vr_favre(3), mass
       
       write (un2,1000) time, T_max, &
            coord_Tmax(1), coord_Tmax(2), coord_Tmax(3), &
            vel_Tmax(1), vel_Tmax(2), vel_Tmax(3), &
            Rloc_Tmax, vr_Tmax, &
            T_center

       write (un3,1000) time, enuc_max, &
            coord_enucmax(1), coord_enucmax(2), coord_enucmax(3), &
            vel_enucmax(1), vel_enucmax(2), vel_enucmax(3), &
            Rloc_enucmax, vr_enucmax, &
            nuc_ener

       write (un4,1000) time, U_max, Mach_max, kin_ener, grav_ener, int_ener, &
            vel_center(1), vel_center(2), vel_center(3), dt

       close(un1)
       close(un2)
       close(un3)
       close(un4)
    endif

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
                     nzones, &
                     vr_x,vr_y,vr_z,vr_max, &
                     rhovr_x,rhovr_y,rhovr_z,mass, &
                     T_max,coord_Tmax,vel_Tmax, &
                     enuc_max,coord_enucmax,vel_enucmax, &
                     kin_ener,int_ener,nuc_ener, &
                     U_max,Mach_max, &
                     ncenter,T_center, &
                     velx_center,vely_center,velz_center, &
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
    real (kind=dp_t), intent(inout) :: vr_x, vr_y, vr_z, vr_max
    real (kind=dp_t), intent(inout) :: rhovr_x, rhovr_y, rhovr_z, mass, nzones
    real (kind=dp_t), intent(inout) :: T_max, coord_Tmax(:)
    real (kind=dp_t), intent(inout) :: vel_Tmax(:)
    real (kind=dp_t), intent(inout) :: enuc_max, coord_enucmax(:)
    real (kind=dp_t), intent(inout) :: vel_enucmax(:)
    real (kind=dp_t), intent(inout) :: kin_ener, int_ener, nuc_ener
    real (kind=dp_t), intent(inout) :: U_max, Mach_max
    integer         , intent(inout) :: ncenter
    real (kind=dp_t), intent(inout) :: T_center, velx_center, vely_center, velz_center
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

             if (cell_valid .and. &
                  s(i,j,k,rho_comp) >= sponge_start_factor*sponge_center_density) then
                   

                ! is it one of the 8 zones surrounding the center?
                if (abs(x - center(1)) < dx(1)  .and. &
                    abs(y - center(2)) < dx(2)  .and. &
                    abs(z - center(3)) < dx(3)) then

                   ncenter = ncenter + 1

                   T_center = T_center + s(i,j,k,temp_comp)

                   velx_center = velx_center + &
                        u(i,j,k,1) + HALF*(w0macx(i,j,k)+w0macx(i+1,j,k))

                   vely_center = vely_center + &
                        u(i,j,k,2) + HALF*(w0macy(i,j,k)+w0macy(i,j+1,k))

                   velz_center = velz_center + &
                        u(i,j,k,3) + HALF*(w0macz(i,j,k)+w0macz(i,j,k+1))
                endif


                ! velr is the projection of the velocity (including w0) onto 
                ! the radial unit vector 
                velr = u(i,j,k,1)*normal(i,j,k,1) + &
                       u(i,j,k,2)*normal(i,j,k,2) + &
                       u(i,j,k,3)*normal(i,j,k,3) + w0r(i,j,k)

                ! vel is the magnitude of the velocity, including w0
                vel = sqrt( (u(i,j,k,1)+HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)))**2 + &
                            (u(i,j,k,2)+HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)))**2 + &
                            (u(i,j,k,3)+HALF*(w0macz(i,j,k)+w0macz(i,j,k+1)))**2)
                

                ! radial velocity diagnostics
                vr_max = max(vr_max,abs(velr))
                
                vr_x = vr_x + weight*velr*normal(i,j,k,1)
                vr_y = vr_y + weight*velr*normal(i,j,k,2)
                vr_z = vr_z + weight*velr*normal(i,j,k,3)
                
                rhovr_x = rhovr_x + weight*s(i,j,k,rho_comp)*velr*normal(i,j,k,1)
                rhovr_y = rhovr_y + weight*s(i,j,k,rho_comp)*velr*normal(i,j,k,2)
                rhovr_z = rhovr_z + weight*s(i,j,k,rho_comp)*velr*normal(i,j,k,3)
                

                ! normalization quantities
                mass = mass + weight*s(i,j,k,rho_comp)
                nzones = nzones + weight


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


                ! call the EOS to get the sound speed and internal energy
                temp_eos(1) = s(i,j,k,temp_comp)
                den_eos(1)  = s(i,j,k,rho_comp)
                xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

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


                ! kinetic, internal, and nuclear energies
                kin_ener = kin_ener + weight*s(i,j,k,rho_comp)*vel**2
                int_ener = int_ener + weight*s(i,j,k,rho_comp)*e_eos(1)
                nuc_ener = nuc_ener + weight*rho_Hnuc(i,j,k)


                ! max vel and Mach number
                U_max = max(U_max,vel)
                Mach_max = max(Mach_max,vel/cs_eos(1))

             endif  ! end cell_valid and density check

          enddo
       enddo
    enddo

  end subroutine diag_3d

end module diag_module
