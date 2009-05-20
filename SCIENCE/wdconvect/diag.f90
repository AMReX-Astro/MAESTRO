! wdconvect-specific diagnostic routine
! 
! currently, there are 4 output files:
!
!   wdconvect_enuc_diag.out:
!          peak nuc energy / g / s
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
!
!   wdconvect_vel_diag.out:
!          peak total velocity
!          peak Mach number
!          total kinetic energy
!          gravitational potential energy
!          total internal energy
!

module diag_module

  use bl_types
  use bl_IO_module
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: diag

contains

  subroutine diag(time,dt,dx,s,rho_Hnuc,rho_Hext, &
                  rho0,rhoh0,p0,tempbar,gamma1bar,div_coeff, &
                  u,w0,normal, &
                  mla,the_bc_tower)

    use bl_prof_module
    use geometry, only: dm, nlevs, spherical, nr_fine, r_cc_loc, r_edge_loc, dr
    use fundamental_constants_module, only: Gconst
    use bl_constants_module
    use variables, only: foextrap_comp
    use fill_3d_module
    use probin_module, only: prob_lo_x, prob_lo_y, prob_lo_z, &
                             prob_hi_x, prob_hi_y, prob_hi_z, &
                             job_name, &
                             edge_nodal_flag, &
                             base_cutoff_density

    real(kind=dp_t), intent(in   ) :: dt,dx(:,:),time
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: rho_Hnuc(:)
    type(multifab) , intent(in   ) :: rho_Hext(:)    
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

    real(kind=dp_t) :: U_max,     U_max_level,     U_max_local
    real(kind=dp_t) :: Mach_max,  Mach_max_level,  Mach_max_local

    real(kind=dp_t) :: xloc_Tmax_local, yloc_Tmax_local, zloc_Tmax_local
    real(kind=dp_t) :: xloc_Tmax_level, yloc_Tmax_level, zloc_Tmax_level
    real(kind=dp_t) :: xloc_Tmax,       yloc_Tmax,       zloc_Tmax

    real(kind=dp_t) :: vx_Tmax_local, vy_Tmax_local, vz_Tmax_local
    real(kind=dp_t) :: vx_Tmax_level, vy_Tmax_level, vz_Tmax_level
    real(kind=dp_t) :: vx_Tmax,       vy_Tmax,       vz_Tmax

    real(kind=dp_t) :: T_max_data_local(1), T_max_coords_local(dm)
    real(kind=dp_t), allocatable :: T_max_data(:), T_max_coords(:)

    integer :: index_max

    real(kind=dp_t) :: vr_favre(dm)

    real(kind=dp_t) :: grav_ener, term1, term2
    real(kind=dp_t), allocatable :: m(:)


    integer :: lo(dm),hi(dm)
    integer :: ng_s,ng_u,ng_n,ng_w,ng_wm,ng_rhn,ng_rhe
    integer :: i,n, comp, r
    integer :: un,un2,un3,un4
    logical :: lexist

    logical, save :: firstCall = .true.

    type(bl_prof_timer), save :: bpt

    call build(bpt, "diagnostics")

    if (spherical .eq. 1) then

       do n=1,nlevs

          do comp=1,dm
             ! w0mac will contain an edge-centered w0 on a Cartesian grid,   
             ! for use in computing divergences.                            
             call multifab_build(w0mac(n,comp), mla%la(n),1,1,nodal=edge_nodal_flag(comp,:))
             call setval(w0mac(n,comp), ZERO, all=.true.)
          enddo

          ! w0r_cart is w0 but onto a Cartesian grid in cell-centered as
          ! a scalar.  Since w0 is the radial expansion velocity, w0r_cart
          ! is the radial w0 in a zone
          call multifab_build(w0r_cart(n), mla%la(n),1,0)
          call setval(w0r_cart(n), ZERO, all=.true.)
       end do

       ! put w0 on Cartesian edges as a vector  
       call put_w0_on_edges(mla,w0,w0mac,dx,div_coeff,the_bc_tower)


       ! put w0 in Cartesian cell-centers as a scalar (the radial 
       ! expansion velocity)
       call put_1d_array_on_cart(w0,w0r_cart,foextrap_comp,.true.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla,normal=normal)

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

    U_max    = ZERO
    Mach_max  = ZERO

    xloc_Tmax = ZERO
    yloc_Tmax = ZERO
    zloc_Tmax = ZERO

    vx_Tmax = ZERO
    vy_Tmax = ZERO
    vz_Tmax = ZERO


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

       U_max_level = ZERO
       U_max_local = ZERO

       Mach_max_level = ZERO
       Mach_max_local = ZERO

       xloc_Tmax_local = ZERO
       yloc_Tmax_local = ZERO
       zloc_Tmax_local = ZERO

       xloc_Tmax_level = ZERO
       yloc_Tmax_level = ZERO
       zloc_Tmax_level = ZERO

       vx_Tmax_local = ZERO
       vy_Tmax_local = ZERO
       vz_Tmax_local = ZERO

       vx_Tmax_level = ZERO
       vy_Tmax_level = ZERO
       vz_Tmax_level = ZERO
       

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
                             T_max_local, xloc_Tmax_local, yloc_Tmax_local, zloc_Tmax_local, &
                             vx_Tmax_local, vy_Tmax_local, vz_Tmax_local, &
                             enuc_max_local, kin_ener_local, int_ener_local, &
                             U_max_local, Mach_max_local)
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
                             T_max_local, xloc_Tmax_local, yloc_Tmax_local, zloc_Tmax_local, &
                             vx_Tmax_local, vy_Tmax_local, vz_Tmax_local, &
                             enuc_max_local, kin_ener_local, int_ener_local, &
                             U_max_local, Mach_max_local, &
                             mp(:,:,:,1))
             end if
          end select
       end do

       !----------------------------------------------------------------------
       ! do the appropriate parallel reduction for the current level
       !----------------------------------------------------------------------

       ! NOTE: only the I/O Processor will have the correct reduced value

       call parallel_reduce(vr_level, vr_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(rhovr_level, rhovr_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(mass_level, mass_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(nzones_level, nzones_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(vr_max_level, vr_max_local, MPI_MAX, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(enuc_max_level, enuc_max_local, MPI_MAX, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(kin_ener_level, kin_ener_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(int_ener_level, int_ener_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(U_max_level, U_max_local, MPI_MAX, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(Mach_max_level, Mach_max_local, MPI_MAX, &
                            proc = parallel_IOProcessorNode())

       ! for T_max, we want to know where the hot spot is, so we do a gather on
       ! the temperature and find the index corresponding to the maxiumum.  We
       ! then pack the coordinates and velocities into a local array and gather 
       ! that to the I/O processor and pick the values corresponding to the 
       ! maximum.
       allocate(T_max_data(parallel_nprocs()))
       T_max_data_local(1) = T_max_local

       call parallel_gather(T_max_data_local, T_max_data, 1, &
                            root = parallel_IOProcessorNode())


       index_max = maxloc(T_max_data, dim=1)
       
       ! T_max_coords will contain both the coordinate information and 
       ! the velocity information, so there are 2*dm values on each proc
       allocate(T_max_coords(2*dm*parallel_nprocs()))
       T_max_coords_local(1) = xloc_Tmax_local
       T_max_coords_local(2) = yloc_Tmax_local
       T_max_coords_local(3) = zloc_Tmax_local
       T_max_coords_local(4) = vx_Tmax_local
       T_max_coords_local(5) = vy_Tmax_local
       T_max_coords_local(6) = vz_Tmax_local
       
       call parallel_gather(T_max_coords_local, T_max_coords, 2*dm, &
                            root = parallel_IOProcessorNode())

       
       T_max_level = T_max_data(index_max)

       xloc_Tmax_level = T_max_coords(2*dm*(index_max-1)+1)
       yloc_Tmax_level = T_max_coords(2*dm*(index_max-1)+2)
       zloc_Tmax_level = T_max_coords(2*dm*(index_max-1)+3)
       vx_Tmax_level   = T_max_coords(2*dm*(index_max-1)+4)
       vy_Tmax_level   = T_max_coords(2*dm*(index_max-1)+5)
       vz_Tmax_level   = T_max_coords(2*dm*(index_max-1)+6)


       deallocate(T_max_data)
       deallocate(T_max_coords)

       ! reduce the current level's data with the global data
       if (parallel_IOProcessor()) then
          vr       = vr     + vr_level
          rhovr    = rhovr  + rhovr_level
          mass     = mass   + mass_level
          nzones   = nzones + nzones_level
          vr_max   = max(vr_max,   vr_max_level)
          enuc_max = max(enuc_max, enuc_max_level)     
          kin_ener = kin_ener + kin_ener_level
          int_ener = int_ener + int_ener_level
          U_max    = max(U_max,    U_max_level)
          Mach_max = max(Mach_max, Mach_max_level)
          
          ! if T_max_level is the new max, then copy the location as well
          if (T_max_level > T_max) then
             T_max = T_max_level

             xloc_Tmax = xloc_Tmax_level
             yloc_Tmax = yloc_Tmax_level
             zloc_Tmax = zloc_Tmax_level

             vx_Tmax = vx_Tmax_level
             vy_Tmax = vy_Tmax_level
             vz_Tmax = vz_Tmax_level
          endif

       endif

    end do

    !-------------------------------------------------------------------------
    ! compute the gravitational potential energy too.
    !-------------------------------------------------------------------------
    allocate(m(0:nr_fine-1))
    grav_ener = 0.0

    ! m(r) will contain mass enclosed by the center
    m(0) = FOUR3RD*M_PI*rho0(1,0)*r_cc_loc(1,0)**3

    ! dU = - G M dM / r;  dM = 4 pi r**2 rho dr  -->  dU = - 4 pi G r rho dr
    grav_ener = -FOUR*M_PI*Gconst*m(0)*r_cc_loc(1,0)*rho0(1,0)*dr(1)

    do r=1,nr_fine-1

       ! the mass is defined at the cell-centers, so to compute the
       ! mass at the current center, we need to add the contribution of
       ! the upper half of the zone below us and the lower half of the
       ! current zone.
       
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
          
       ! dU = - G M dM / r;  dM = 4 pi r**2 rho dr  -->  dU = - 4 pi G r rho dr
       grav_ener = grav_ener - FOUR*M_PI*Gconst*m(r)*r_cc_loc(1,r)*rho0(1,r)*dr(1)

    enddo


    

    !=========================================================================
    ! normalize
    !=========================================================================
    vr(:) = vr(:)/nzones
    vr_favre(:) = rhovr(:)/mass    ! note, the common dV normalization cancels

    ! the volume we normalize with is that of a single coarse-level zone.
    ! This is because the weight used in the loop over cells was with reference
    ! to the coarse level

    mass = mass*dx(1,1)*dx(1,2)*dx(1,3)
    kin_ener = kin_ener*dx(1,1)*dx(1,2)*dx(1,3)
    int_ener = int_ener*dx(1,1)*dx(1,2)*dx(1,3)


    !=========================================================================
    ! output
    !=========================================================================
 999 format("# job name: ",a)
1000 format(1x,10(g20.10,1x))
1001 format("#",10(a20,1x))

    if (parallel_IOProcessor()) then

       ! open the diagnostic files for output, taking care not to overwrite
       ! an existing file
       un = unit_new()
       inquire(file="wdconvect_radvel_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un, file="wdconvect_radvel_diag.out", &
               status="old", position="append")
       else
          open(unit=un, file="wdconvect_radvel_diag.out", status="new")
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
          
          write (un, *) " "
          write (un, 999) trim(job_name)
          write (un, 1001) "time", "<vr_x>", "<vr_y>", "<vr_z>", "<vr>", &
                           "max{|vr|}", &
                           "int{rhovr_x}/mass", "int{rhovr_y}/mass", "int{rhovr_z}/mass", &
                           "mass"

          write (un2, *) " "
          write (un2, 999) trim(job_name)
          write (un2,1001) "time", "max{T}", 'x(max{T})', 'y(max{T})', 'z(max{T})', &
               'vx(max{T})', 'vy(max{T})', 'vz(max{T})'

          write (un3, *) " "
          write (un3, 999) trim(job_name)
          write (un3,1001) "time", "max{enuc}"

          write (un4, *) " "
          write (un4, 999) trim(job_name)
          write (un4,1001) "time", "max{|U + w0|}", "max{Mach #}", "tot. kin. energy", "grav. pot. energy", "tot. int. energy"

          firstCall = .false.
       endif

       ! write out the data
       write (un,1000) time, vr(1), vr(2), vr(3), &
            sqrt(vr(1)**2 + vr(2)**2 + vr(3)**2), vr_max, &
            vr_favre(1), vr_favre(2), vr_favre(3), mass
       
       write (un2,1000) time, T_max, xloc_Tmax, yloc_Tmax, zloc_Tmax, vx_Tmax, vy_Tmax, vz_Tmax

       write (un3,1000) time, enuc_max

       write (un4,1000) time, U_max, Mach_max, kin_ener, grav_ener, int_ener

       close(un)
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
                     T_max,xloc_Tmax,yloc_Tmax,zloc_Tmax, &
                     vx_Tmax, vy_Tmax, vz_Tmax, &
                     enuc_max,kin_ener,int_ener, &
                     U_max,Mach_max, &
                     mask)

    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp
    use bl_constants_module
    use network, only: nspec
    use geometry, only: spherical
    use probin_module, only: base_cutoff_density, prob_lo
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
    real (kind=dp_t), intent(inout) :: T_max, xloc_Tmax, yloc_Tmax, zloc_Tmax
    real (kind=dp_t), intent(inout) :: vx_Tmax, vy_Tmax, vz_Tmax
    real (kind=dp_t), intent(inout) :: enuc_max, kin_ener, int_ener
    real (kind=dp_t), intent(inout) :: U_max, Mach_max
    logical,          intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    !     Local variables
    integer            :: i, j, k
    real (kind=dp_t)   :: velr, vel, weight
    logical            :: cell_valid
    real (kind=dp_t)   :: x, y, z

    ! weight is the factor by which the volume of a cell at the current level 
    ! relates to the volume of a cell at the coarsest level of refinement.
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

             if (cell_valid .and. s(i,j,k,rho_comp) > base_cutoff_density) then
                   

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
                   xloc_Tmax = x
                   yloc_Tmax = y
                   zloc_Tmax = z
                   vx_Tmax = u(i,j,k,1)+HALF*(w0macx(i,j,k)+w0macx(i+1,j,k))
                   vy_Tmax = u(i,j,k,2)+HALF*(w0macy(i,j,k)+w0macy(i,j+1,k))
                   vz_Tmax = u(i,j,k,3)+HALF*(w0macz(i,j,k)+w0macz(i,j,k+1))
                endif


                ! max enuc
                enuc_max = max(enuc_max,rho_Hnuc(i,j,k)/s(i,j,k,rho_comp))


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


                ! kinetic and internal energies
                kin_ener = kin_ener + weight*s(i,j,k,rho_comp)*vel**2
                int_ener = int_ener + weight*s(i,j,k,rho_comp)*e_eos(1)               


                ! max vel and Mach number
                U_max = max(U_max,vel)
                Mach_max = max(Mach_max,vel/cs_eos(1))

             endif

          enddo
       enddo
    enddo

  end subroutine diag_3d

end module diag_module
