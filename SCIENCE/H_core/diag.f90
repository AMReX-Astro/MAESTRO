! hcore-specific diagnostic routine
! 
! currently, there are 4 output files:
!
!
!   hcore_cz_diag.out:          
!          radius of convective boundary 
!
!   hcore_sphrvel_diag.out:
!          radial velocity components (std. average & Favre average)
!          peak radial velocity
!          x/y/z location of peak rad vel
!          circumferential velocity components (std. average & Favre average)
!          peak circumferential velocity
!          x/y/z location 
!          average (v_rad/v_tot)
!          total mass
!
!   hcore_vel_diag.out:
!          peak total velocity
!          x/y/z location of peak vel
!          peak Mach number
!
!   hcore_energy_diag.out:
!          total kinetic energy
!          total internal energy
!          gravitational potential energy
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

  integer, parameter :: n_file1 = 18 
  integer, parameter :: n_file2 = 4
  integer, parameter :: n_file3 = 1
  integer, parameter :: n_file4 = 7

  integer, save :: nstored = 0

  public :: diag, flush_diag

contains

  subroutine diag(time,dt,dx,s,rho_Hnuc,rho_Hext,thermal,rho_omegadot, &
                  rho0,rhoh0,p0,tempbar, &
                  gamma1bar,div_coeff, &
                  u,w0,normal, &
                  mla,the_bc_tower)

    use parallel
    use bl_constants_module, only: ZERO, HALF, TWO, FOUR, FOUR3RD, M_PI
    use bl_prof_module, only: bl_prof_timer, build
    use bl_error_module, only: bl_error

    use fundamental_constants_module, only: Gconst

    use fab_module, only: lwb, upb
    use multifab_module, only: multifab, multifab_build, destroy, &
                               multifab_remote, dataptr, setval, get_box, &
                               multifab_div_div_c, multifab_copy_c
    use ml_layout_module, only: ml_layout
    use define_bc_module, only: bc_tower

    use geometry, only: dm, nlevs, spherical, nr_fine, &
                        r_cc_loc, r_edge_loc, dr, center
    use variables, only: foextrap_comp, rho_comp, spec_comp
    use fill_3d_module, only: put_1d_array_on_cart, make_w0mac
    use probin_module, only: prob_lo_x, prob_lo_y, prob_lo_z, &
                             prob_hi_x, prob_hi_y, prob_hi_z, &
                             edge_nodal_flag, &
                             base_cutoff_density, &
                             diag_buf_size, octant, &
                             sponge_start_factor, sponge_center_density
    use average_module, only: average

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
    real(kind=dp_t) :: vrvt,    vrvt_level,    vrvt_local

    real(kind=dp_t) :: vc(dm),    vc_level(dm),    vc_local(dm)
    real(kind=dp_t) :: vc_max,    vc_max_level,    vc_max_local
    real(kind=dp_t) :: rhovc(dm), rhovc_level(dm), rhovc_local(dm)

    real(kind=dp_t) :: mass,      mass_level,      mass_local
    real(kind=dp_t) :: nzones,    nzones_level,    nzones_local

    real(kind=dp_t) :: kin_ener,  kin_ener_level,  kin_ener_local
    real(kind=dp_t) :: int_ener,  int_ener_level,  int_ener_local
    real(kind=dp_t) :: nuc_ener,  nuc_ener_level,  nuc_ener_local

    real(kind=dp_t) :: U_max,     U_max_level,     U_max_local
    real(kind=dp_t) :: coord_Umax_local(dm), coord_Umax_level(dm), coord_Umax(dm)
    real(kind=dp_t) :: Mach_max,  Mach_max_level,  Mach_max_local

    real(kind=dp_t) :: r_cz

    ! buffers
    real(kind=dp_t) :: max_data_level(3), max_data_local(3)
    real(kind=dp_t) :: sum_data_level(2*dm+9), sum_data_local(2*dm+9)

    real(kind=dp_t) :: U_max_data_local(1), U_max_coords_local(2*dm)
    real(kind=dp_t), allocatable :: U_max_data(:), U_max_coords(:)

    type(multifab)  :: XH(mla%nlevel)
    real(kind=dp_t) :: XH_avg(1,0:nr_fine-1)
    real(kind=dp_t) :: grad_XH(0:nr_fine-1)

    integer :: i_cz

    integer :: index

    integer :: index_max

    real(kind=dp_t) :: vr_favre(dm)
    real(kind=dp_t) :: vc_favre(dm)

    real(kind=dp_t) :: grav_ener, term1, term2
    real(kind=dp_t), allocatable :: m(:)


    integer :: lo(dm),hi(dm)
    integer :: ng_s,ng_u,ng_n,ng_w,ng_wm,ng_rhn,ng_rhe
    integer :: i,n, comp, r

    type(bl_prof_timer), save :: bpt

    ! the maximum number of quantities to store in a size file -- for the 
    ! buffering
    integer, parameter :: MAX_FIELDS_PER_FILE = 32

    logical, save :: firstCall = .true.

    call build(bpt, "diagnostics")


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

       call bl_error("ERROR: hcore/diag.f90: geometry not spherical")

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
    vrvt    = ZERO

    vc(:)    = ZERO
    rhovc(:) = ZERO

    mass     = ZERO
    nzones   = ZERO

    vr_max   = ZERO

    vc_max   = ZERO

    kin_ener = ZERO
    int_ener = ZERO
    nuc_ener = ZERO

    U_max    = ZERO
    coord_Umax(:) = ZERO

    Mach_max  = ZERO

    r_cz = ZERO



    !=========================================================================
    ! loop over the levels and compute the global quantities
    !=========================================================================
    do n = 1, nlevs

       ! initialize the local (processor's version) and level quantities to 0
       vr_level(:) = ZERO
       vr_local(:) = ZERO

       vrvt_level = ZERO
       vrvt_local = ZERO
       
       rhovr_level(:) = ZERO
       rhovr_local(:) = ZERO

       vc_level(:) = ZERO
       vc_local(:) = ZERO
       
       rhovc_level(:) = ZERO
       rhovc_local(:) = ZERO
       
       mass_level = ZERO
       mass_local = ZERO
       
       nzones_level = ZERO
       nzones_local = ZERO
       
       vr_max_level = ZERO
       vr_max_local = ZERO

       vc_max_level = ZERO
       vc_max_local = ZERO
       
       kin_ener_level = ZERO
       kin_ener_local = ZERO

       int_ener_level = ZERO
       int_ener_local = ZERO

       U_max_level = ZERO
       U_max_local = ZERO

       coord_Umax_local(:) = ZERO
       coord_Umax_level(:) = ZERO

       Mach_max_level = ZERO
       Mach_max_local = ZERO

       r_cz = ZERO
       r_cz = ZERO


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
                             rhovr_local(1), rhovr_local(2), rhovr_local(3),  &
                             vrvt_local, &
                             vc_local(1),vc_local(2),vc_local(3),vc_max_local, &
                             rhovc_local(1), rhovc_local(2), rhovc_local(3), &
                             mass_local, &
                             kin_ener_local, int_ener_local, &
                             U_max_local, coord_Umax_local, & 
                             Mach_max_local)
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
                             rhovr_local(1), rhovr_local(2), rhovr_local(3), &
                             vrvt_local, &
                             vc_local(1),vc_local(2),vc_local(3),vc_max_local, &
                             rhovc_local(1), rhovc_local(2), rhovc_local(3), &
                             mass_local, &
                             kin_ener_local, int_ener_local, &
                             U_max_local,  coord_Umax_local, & 
                             Mach_max_local, mp(:,:,:,1))
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
       sum_data_local(2*dm+1:3*dm) = vc_local(:)
       sum_data_local(3*dm+1:4*dm) = rhovc_local(:)
       sum_data_local(4*dm+1)    = mass_local
       sum_data_local(4*dm+2)    = nzones_local
       sum_data_local(4*dm+3)    = kin_ener_local
       sum_data_local(4*dm+4)    = int_ener_local
       sum_data_local(4*dm+5)    = vrvt_local


       call parallel_reduce(sum_data_level, sum_data_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       vr_level(:)         = sum_data_level(1:dm)
       rhovr_level(:)      = sum_data_level(dm+1:2*dm)
       vc_level(:)         = sum_data_level(2*dm+1:3*dm)
       rhovc_level(:)      = sum_data_level(3*dm+1:4*dm)
       mass_level          = sum_data_level(4*dm+1)
       nzones_level        = sum_data_level(4*dm+2)
       kin_ener_level      = sum_data_level(4*dm+3)
       int_ener_level      = sum_data_level(4*dm+4)
       vrvt_level          = sum_data_level(4*dm+5)


       ! pack the quantities that we are taking the max of into a vector
       ! to reduce communication
       max_data_local(1) = vr_max_local
       max_data_local(2) = vc_max_local
       max_data_local(3) = Mach_max_local

       call parallel_reduce(max_data_level, max_data_local, MPI_MAX, &
                            proc = parallel_IOProcessorNode())

       vr_max_level   = max_data_level(1)
       vc_max_level   = max_data_level(2)
       Mach_max_level = max_data_level(3)


       ! for U_max, we want to know where the spot is, so we do a
       ! gather on U and find the index corresponding to
       ! the maxiumum.  We then pack the coordinates 
       ! into a local array and gather that to the I/O processor and
       ! pick the values corresponding to the maximum.
       allocate(U_max_data(parallel_nprocs()))
       U_max_data_local(1) = U_max_local

       call parallel_gather(U_max_data_local, U_max_data, 1, &
                            root = parallel_IOProcessorNode())


       index_max = maxloc(U_max_data, dim=1)
       
       ! U_max_coords will contain both the coordinate information and
       ! the velocity information, so there are 2*dm values on each
       ! proc
       allocate(U_max_coords(dm*parallel_nprocs()))
       U_max_coords_local(1) = coord_Umax_local(1)
       U_max_coords_local(2) = coord_Umax_local(2)
       U_max_coords_local(3) = coord_Umax_local(3)
       
       call parallel_gather(U_max_coords_local, U_max_coords, dm, &
                            root = parallel_IOProcessorNode())

       
       U_max_level = U_max_data(index_max)

       coord_Umax_level(1) = U_max_coords(dm*(index_max-1)+1)
       coord_Umax_level(2) = U_max_coords(dm*(index_max-1)+2)
       coord_Umax_level(3) = U_max_coords(dm*(index_max-1)+3)


       deallocate(U_max_data)
       deallocate(U_max_coords)



       !----------------------------------------------------------------------
       ! reduce the current level's data with the global data
       !----------------------------------------------------------------------
       if (parallel_IOProcessor()) then
          vr       = vr     + vr_level
          rhovr    = rhovr  + rhovr_level
          vrvt     = vrvt     + vrvt_level

          vc       = vc     + vc_level
          rhovc    = rhovc  + rhovc_level

          mass     = mass   + mass_level
          nzones   = nzones + nzones_level

          vr_max   = max(vr_max,   vr_max_level)

          vc_max   = max(vc_max,   vc_max_level)

          kin_ener = kin_ener + kin_ener_level
          int_ener = int_ener + int_ener_level

          Mach_max = max(Mach_max, Mach_max_level)
          
          ! if U_max_level is the new max, then copy the location as well
          if (U_max_level > U_max) then
             U_max = U_max_level

             coord_Umax(:) = coord_Umax_level(:)

          endif

       endif

    end do


    !-------------------------------------------------------------------------
    ! compute the location of convection zone boundary
    !-------------------------------------------------------------------------
    ! get X_H by itself
    do n = 1, nlevs
       call multifab_build(XH(n), mla%la(n), 1, s(n)%ng)
       call multifab_copy_c(XH(n),1,s(n),spec_comp, 1,s(n)%ng)
       call multifab_div_div_c(XH(n), 1, s(n), rho_comp, 1, s(n)%ng)
    enddo

    call average(mla,XH,XH_avg,dx,spec_comp)

    ! compute the radial gradient of X_H
    do n = 1,nr_fine-2
       grad_XH(n) = abs((XH_avg(1,n+1) - XH_avg(1,n-1))/(TWO*dr(1))) 
       write(11,*)n, XH_avg(1,n), grad_XH(n)
    end do

    r_cz = ZERO
    i_cz = -1
    do n = 1,nr_fine-2
       if ( grad_XH(n) > r_cz ) then
          r_cz = grad_XH(n)
          i_cz = n
       end if
    end do
    write(*,*) i, r_cz
    r_cz = dr(1)*(dble(i_cz) + HALF)
    write(*,*)r_cz
    stop
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
       
       ! don't add any contributions from the sponged region
       if (rho0(1,r-1) > sponge_start_factor*sponge_center_density) then
          term1 = FOUR3RD*M_PI*rho0(1,r-1) * &
               (r_edge_loc(1,r) - r_cc_loc(1,r-1)) * &
               (r_edge_loc(1,r)**2 + &
                r_edge_loc(1,r)*r_cc_loc(1,r-1) + &
                r_cc_loc(1,r-1)**2)
       else
          term1 = ZERO
       endif

       if (rho0(1,r) > sponge_start_factor*sponge_center_density) then
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
       vrvt = vrvt/nzones

       vc(:) = vc(:)/nzones
       vc_favre(:) = rhovc(:)/mass    ! note: common dV normalization cancels

       ! the volume we normalize with is that of a single coarse-level
       ! zone.  This is because the weight used in the loop over cells
       ! was with reference to the coarse level

       mass = mass*dx(1,1)*dx(1,2)*dx(1,3)
       kin_ener = kin_ener*dx(1,1)*dx(1,2)*dx(1,3)
       int_ener = int_ener*dx(1,1)*dx(1,2)*dx(1,3)

    endif


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
       ! file1 -- hcore_sphrvel_diag.out
       file1_data(index, 1) = vr(1)
       file1_data(index, 2) = vr(2)
       file1_data(index, 3) = vr(3)
       file1_data(index, 4) = sqrt(vr(1)**2 + vr(2)**2 + vr(3)**2)
       file1_data(index, 5) = vr_max
       file1_data(index, 6) = vrvt
       file1_data(index, 7) = vr_favre(1)
       file1_data(index, 8) = vr_favre(2)
       file1_data(index, 9) = vr_favre(3)
       file1_data(index, 10) = vc(1)
       file1_data(index, 11) = vc(2)
       file1_data(index, 12) = vc(3)
       file1_data(index, 13) = sqrt(vc(1)**2 + vc(2)**2 + vc(3)**2)
       file1_data(index, 14) = vc_max
       file1_data(index, 15) = vc_favre(1)
       file1_data(index, 16) = vc_favre(2)
       file1_data(index, 17) = vc_favre(3)
       file1_data(index, 18) = mass


       ! file2 -- hcore_ener_diag.out
       file2_data(index, 1) = kin_ener
       file2_data(index, 2) = grav_ener
       file2_data(index, 3) = int_ener
       file2_data(index, 4) = dt


       ! file3 -- hcore_cz_diag.out
       file3_data(index, 1) = r_cz


       ! file4 -- hcore_vel_diag.out
       file4_data(index, 1) = U_max
       file4_data(index, 2) = coord_Umax(1)
       file4_data(index, 3) = coord_Umax(2)
       file4_data(index, 4) = coord_Umax(3)
       file4_data(index, 5) = sqrt( (coord_Umax(1) - center(1))**2 + &
                                    (coord_Umax(2) - center(2))**2 + &
                                    (coord_Umax(3) - center(3))**2 )
       file4_data(index, 6) = Mach_max
       file4_data(index, 7) = dt
    end if

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

    integer :: n,i

    logical, save :: firstCall = .true.


    ! if the buffers are empty, move on
    if (nstored == 0) return


    ! IMPORTANT: make sure that there are enough entries in the format
    ! statement to write out all of the data in each file.
999 format("# job name: ",a)
1000 format(1x,20(g20.10,1x))
1001 format("#",20(a20,1x))
800 format("# ",a,i4.4,'-',i2.2,'-',i2.2)
801 format("# ",a,i2.2,':',i2.2,':',i2.2)
802 format("# ",a,a)

    if (parallel_IOProcessor()) then

       ! open the diagnostic files for output, taking care not to overwrite
       ! an existing file
       un1 = unit_new()
       inquire(file="hcore_sphrvel_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un1, file="hcore_sphrvel_diag.out", &
               status="old", position="append")
       else
          open(unit=un1, file="hcore_sphrvel_diag.out", status="new")
       endif

       un2 = unit_new()
       inquire(file="hcore_ener_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un2, file="hcore_ener_diag.out", &
               status="old", position="append")
       else
          open(unit=un2, file="hcore_ener_diag.out", status="new")
       endif

       un3 = unit_new()
       inquire(file="hcore_cz_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un3, file="hcore_cz_diag.out", &
               status="old", position="append")
       else
          open(unit=un3, file="hcore_cz_diag.out", status="new")
       endif

       un4 = unit_new()
       inquire(file="hcore_vel_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un4, file="hcore_vel_diag.out", &
               status="old", position="append")
       else
          open(unit=un4, file="hcore_vel_diag.out", status="new")
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
                            "max{|vr|}", "<|vr|/|vtot|>", &
                            "int{rhovr_x}/mass", "int{rhovr_y}/mass", "int{rhovr_z}/mass", &
                            "<vc_x>", "<vc_y>", "<vc_z>", "<vc>", &
                            "max{|vc|}", &
                            "int{rhovc_x}/mass", "int{rhovc_y}/mass", "int{rhovc_z}/mass", &
                            "mass"

          ! energy
          write (un2, *) " "
          write (un2, 800) "output date: ", values(1), values(2), values(3)
          write (un2, 801) "output time: ", values(5), values(6), values(7)
          write (un2, 802) "output dir:  ", trim(cwd)
          write (un2, 999) trim(job_name)
          write (un2,1001) "time", "tot kin energy", "grav pot energy", &
               "tot int energy", "dt"

          ! convective boundary
          write (un3, *) " "
          write (un3, 800) "output date: ", values(1), values(2), values(3)
          write (un3, 801) "output time: ", values(5), values(6), values(7)
          write (un3, 802) "output dir:  ", trim(cwd)
          write (un3, 999) trim(job_name)
          write (un3,1001) "time", "radius of convective boundary"
               
          ! vel
          write (un4, *) " "
          write (un4, 800) "output date: ", values(1), values(2), values(3)
          write (un4, 801) "output time: ", values(5), values(6), values(7)
          write (un4, 802) "output dir:  ", trim(cwd)
          write (un4, 999) trim(job_name)
          write (un4,1001) "time", "max{|U + w0|}",  &
               "x(max{V})", "y(max{V})", "z(max{V})", "R{max{V})", &
               "max{Mach #}", "dt"
  
          firstCall = .false.
       endif

       do n = 1, nstored

          ! write out the data
          write (un1,1000) time_data(n), &
               (file1_data(n,i), i=1,n_file1)
       
          write (un2,1000) time_data(n), &
               (file2_data(n,i), i=1,n_file2)

          write (un3,1000) time_data(n), &
               (file3_data(n,i), i=1,n_file3)

          write (un4,1000) time_data(n), &
               (file4_data(n,i), i=1,n_file4)

       enddo

       close(un1)
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
                     nzones, &
                     vr_x,vr_y,vr_z,vr_max, &
                     rhovr_x,rhovr_y,rhovr_z, vrvt, &
                     vc_x,vc_y,vc_z,vc_max, &
                     rhovc_x,rhovc_y,rhovc_z, &
                     mass, &
                     kin_ener,int_ener, &
                     U_max, coord_Umax, Mach_max, &
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
    real (kind=dp_t), intent(inout) :: rhovr_x, rhovr_y, rhovr_z, vrvt
    real (kind=dp_t), intent(inout) :: vc_x, vc_y, vc_z, vc_max
    real (kind=dp_t), intent(inout) :: rhovc_x, rhovc_y, rhovc_z, mass, nzones
    real (kind=dp_t), intent(inout) :: kin_ener, int_ener
    real (kind=dp_t), intent(inout) :: U_max, coord_Umax(:)
    real (kind=dp_t), intent(inout) :: Mach_max
    logical,          intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    !     Local variables
    integer            :: i, j, k
    real (kind=dp_t)   :: velr, velc, vel, weight
    logical            :: cell_valid
    real (kind=dp_t)   :: x, y, z

    ! weight is the factor by which the volume of a cell at the
    ! current level relates to the volume of a cell at the coarsest
    ! level of refinement.
    weight = 1.d0 / 8.d0**(n-1)

    if (.not. spherical == 1) then
       call bl_error("ERROR: geometry not spherical in diag")
    endif

!$omp parallel do private(i,j,k,x,y,z,cell_valid,velr,vel,vrvt) &
!$omp reduction(max:vr_max,U_max,Mach_max) &
!$omp reduction(vr_x,vr_y,vr_z,rhovr_x,rhovr_y,rhovr_z, &
!$omp           vc_x,vc_y,vc_z,rhovc_x,rhovc_y,rhovc_z,mass,nzones,kin_ener,int_ener)
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
             if (cell_valid .and. &
                  s(i,j,k,rho_comp) > sponge_start_factor*sponge_center_density) then
                   
                ! velr is the projection of the velocity (including w0 below) 
                ! onto the radial unit vector 
                velr = u(i,j,k,1)*normal(i,j,k,1) + &
                       u(i,j,k,2)*normal(i,j,k,2) + &
                       u(i,j,k,3)*normal(i,j,k,3) 

                ! vel is the magnitude of the velocity, including w0
                vel = sqrt( (u(i,j,k,1)+HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)))**2 + &
                            (u(i,j,k,2)+HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)))**2 + &
                            (u(i,j,k,3)+HALF*(w0macz(i,j,k)+w0macz(i,j,k+1)))**2)
                velc = sqrt( (u(i,j,k,1) - velr*normal(i,j,k,1))**2 + &
                             (u(i,j,k,2) - velr*normal(i,j,k,2))**2 + &
                             (u(i,j,k,3) - velr*normal(i,j,k,3))**2)


                ! diagnostics
                vc_max = max(vc_max,abs(velc))

                vc_x = vc_x + weight*(u(i,j,k,1)-velr*normal(i,j,k,1))
                vc_y = vc_y + weight*(u(i,j,k,2)-velr*normal(i,j,k,2))
                vc_z = vc_z + weight*(u(i,j,k,3)-velr*normal(i,j,k,3))
                
                rhovc_x = rhovc_x + weight*s(i,j,k,rho_comp)*(u(i,j,k,1)-velr*normal(i,j,k,1))
                rhovc_y = rhovc_y + weight*s(i,j,k,rho_comp)*(u(i,k,j,2)-velr*normal(i,j,k,2))
                rhovc_z = rhovc_z + weight*s(i,j,k,rho_comp)*(u(i,j,k,3)-velr*normal(i,j,k,3))


                ! add in w0 before computing radial velocity diagnostics
                velr = velr + w0r(i,j,k)                

                vr_max = max(vr_max,abs(velr))
                
                vrvt = vrvt + weight*velr/vel

                vr_x = vr_x + weight*velr*normal(i,j,k,1)
                vr_y = vr_y + weight*velr*normal(i,j,k,2)
                vr_z = vr_z + weight*velr*normal(i,j,k,3)
                
                rhovr_x = rhovr_x + weight*s(i,j,k,rho_comp)*velr*normal(i,j,k,1)
                rhovr_y = rhovr_y + weight*s(i,j,k,rho_comp)*velr*normal(i,j,k,2)
                rhovr_z = rhovr_z + weight*s(i,j,k,rho_comp)*velr*normal(i,j,k,3)
                
                ! normalization quantities
                mass = mass + weight*s(i,j,k,rho_comp)
                nzones = nzones + weight

!$omp critical
                ! Mach number
                Mach_max = max(Mach_max,vel/cs_eos(1))

                ! max U and  location
                if (vel > U_max) then
                   U_max = vel
                   coord_Umax(1) = x
                   coord_Umax(2) = y
                   coord_Umax(3) = z
                endif

!$omp end critical

                ! call the EOS to get the sound speed and internal energy
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


                ! kinetic, internal, and nuclear energies
                kin_ener = kin_ener + weight*s(i,j,k,rho_comp)*vel**2
                int_ener = int_ener + weight*s(i,j,k,rho_comp)*e_eos(1)

             endif  ! end cell_valid and density check

          enddo
       enddo
    enddo

  end subroutine diag_3d

end module diag_module
