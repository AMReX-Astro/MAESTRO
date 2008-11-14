! This is a general interface for doing runtime diagnostics on the state.
! It is called at the end of advance

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
                  rho0,rhoh0,p0,tempbar,gamma1bar, &
                  u,w0,normal, &
                  mla,the_bc_tower)

    use bl_prof_module
    use geometry, only: dm, nlevs, spherical
    use bl_constants_module
    use variables, only: foextrap_comp
    use fill_3d_module
    use probin_module, only: prob_lo_x, prob_lo_y, prob_lo_z, &
                             prob_hi_x, prob_hi_y, prob_hi_z

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
    real(kind=dp_t), intent(in   ) ::        w0(:,0:)
    type(ml_layout), intent(in   ) :: mla
    type(bc_tower) , intent(in   ) :: the_bc_tower


    ! Local
    real(kind=dp_t), pointer::  sp(:,:,:,:)
    real(kind=dp_t), pointer::  rhnp(:,:,:,:)
    real(kind=dp_t), pointer::  rhep(:,:,:,:)
    real(kind=dp_t), pointer::  up(:,:,:,:)
    real(kind=dp_t), pointer::  np(:,:,:,:)
    real(kind=dp_t), pointer::  w0rp(:,:,:,:)

    type(multifab) ::    w0r_cart(mla%nlevel)

    real(kind=dp_t) :: vr(dm), vr_local(dm), vr_max, vr_max_local
    real(kind=dp_t) :: rhovr(dm), rhovr_local(dm), mass, mass_local
    integer         :: nzones, nzones_local
    real(kind=dp_t) :: T_max, T_max_local
    real(kind=dp_t) :: enuc_max, enuc_max_local

    integer :: lo(dm),hi(dm),ng_s,ng_u,ng_n,ng_w,ng_rhn,ng_rhe
    integer :: i,n
    integer :: un, un2, un3
    logical :: lexist

    logical, save :: firstCall = .true.

    type(bl_prof_timer), save :: bpt

    call build(bpt, "diagnostics")

    if (spherical .eq. 1) then

       ! even though we are looping over levels, this is not multilevel
       ! aware.  We are assuming here that there is only 1 level.
       do n=1,nlevs

          ! w0r_cart is w0 but onto a Cartesian grid in cell-centered as
          ! a scalar.  Since w0 is the radial expansion velocity, w0r_cart
          ! is the radial w0 in a zone
          call multifab_build(w0r_cart(n), mla%la(n),1,0)
          call setval(w0r_cart(n), ZERO, all=.true.)
       end do

       ! put w0 in Cartesian cell-centers as a scalar (the radial 
       ! expansion velocity)
       call put_1d_array_on_cart(w0,w0r_cart,foextrap_comp,.true.,.false.,dx, &
            the_bc_tower%bc_tower_array,mla,normal=normal)

    else
       call bl_error("ERROR: geometry not spherical")
    endif


    ng_s = s(1)%ng
    ng_u = u(1)%ng
    ng_n = normal(1)%ng
    ng_w = w0r_cart(1)%ng
    ng_rhn = rho_Hnuc(1)%ng
    ng_rhe = rho_Hext(1)%ng

    vr(:) = ZERO
    vr_local(:) = ZERO

    rhovr(:) = ZERO
    rhovr_local(:) = ZERO

    mass = ZERO
    mass_local = ZERO

    nzones = 0
    nzones_local = 0

    vr_max = ZERO
    vr_max_local = ZERO

    T_max = ZERO
    T_max_local = ZERO

    enuc_max = ZERO
    enuc_max_local = ZERO


    do n = 1, nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          sp => dataptr(s(n) , i)
          rhnp => dataptr(rho_Hnuc(n), i)
          rhep => dataptr(rho_Hext(n), i)
          up => dataptr(u(n) , i)
          np => dataptr(normal(n) , i)
          w0rp => dataptr(w0r_cart(n), i)

          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))

          select case (dm)
          case (2)
             call bl_error("ERROR: 2-d diag not implmented")
          case (3)
             call diag_3d(time,dt,dx(n,:), &
                          sp(:,:,:,:),ng_s, &
                          rhnp(:,:,:,1), ng_rhn, &
                          rhep(:,:,:,1), ng_rhe, &
                          up(:,:,:,:),ng_u, &
                          w0rp(:,:,:,1), ng_w, &
                          np(:,:,:,:),ng_n, &
                          lo,hi, &
                          nzones_local, &
                          vr_local(1),vr_local(2),vr_local(3),vr_max_local, &
                          rhovr_local(1), rhovr_local(2), rhovr_local(3), mass_local, &
                          T_max_local, enuc_max_local)
          end select
       end do
    end do

    ! we now have vr_local on each processor -- do a reduce on
    ! all dm components
    call parallel_reduce(vr,vr_local, MPI_SUM, &
                         proc = parallel_IOProcessorNode())

    call parallel_reduce(rhovr,rhovr_local, MPI_SUM, &
                         proc = parallel_IOProcessorNode())

    call parallel_reduce(mass,mass_local, MPI_SUM, &
                         proc = parallel_IOProcessorNode())

    call parallel_reduce(nzones,nzones_local, MPI_SUM, &
                         proc = parallel_IOProcessorNode())

    call parallel_reduce(vr_max,vr_max_local, MPI_MAX, &
                         proc = parallel_IOProcessorNode())

    call parallel_reduce(T_max,T_max_local, MPI_MAX, &
                         proc = parallel_IOProcessorNode())

    call parallel_reduce(enuc_max,enuc_max_local, MPI_MAX, &
                         proc = parallel_IOProcessorNode())

    
    ! normalize
    vr = vr/nzones
    
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

       ! write out the headers
       if (firstCall) then
          write (un, *) " "
          write (un, 1001) "time", "<vr_x>", "<vr_y>", "<vr_z>", "<vr>", &
                           "max{|vr|}", &
                           "int{rhovr_x}/mass", "int{rhovr_y}/mass", "int{rhovr_z}/mass", &
                           "mass"

          write (un2, *) " "
          write (un2,1001) "time", "max{T}"

          write (un3, *) " "
          write (un3,1001) "time", "max{enuc}"
          firstCall = .false.
       endif

       ! write out the data
       write (un,1000) time, vr(1), vr(2), vr(3), &
            sqrt(vr(1)**2 + vr(2)**2 + vr(3)**2), vr_max, &
            rhovr(1)/mass, rhovr(2)/mass, rhovr(3)/mass, mass*dx(1,1)*dx(1,2)*dx(1,3)
       
       write (un2,1000) time, T_max

       write (un3,1000) time, enuc_max

       close(un)
       close(un2)
       close(un3)
    endif

    if (spherical .eq. 1) then
       do n=1,nlevs
          call destroy(w0r_cart(n))
       end do
    end if




    call destroy(bpt)

  end subroutine diag

  subroutine diag_3d(time,dt,dx, &
                     s,ng_s, &
                     rho_Hnuc,ng_rhn, &
                     rho_Hext,ng_rhe, &
                     u,ng_u, &
                     w0r,ng_w, &
                     normal,ng_n, &
                     lo,hi, &
                     nzones, &
                     vr_x,vr_y,vr_z,vr_max, &
                     rhovr_x,rhovr_y,rhovr_z,mass, &
                     T_max,enuc_max)

    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp
    use network, only: nspec
    use geometry, only: spherical
    use probin_module, only: base_cutoff_density

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_u, ng_n, ng_w, ng_rhn, ng_rhe
    real (kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho_Hnuc(lo(1)-ng_rhn:,lo(2)-ng_rhn:,lo(3)-ng_rhn:)
    real (kind=dp_t), intent(in   ) :: rho_Hext(lo(1)-ng_rhe:,lo(2)-ng_rhe:,lo(3)-ng_rhe:)
    real (kind=dp_t), intent(in   ) ::        u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
    real (kind=dp_t), intent(in   ) ::      w0r(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:)
    real (kind=dp_t), intent(in   ) ::   normal(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)
    real (kind=dp_t), intent(in   ) :: time, dt, dx(:)
    real (kind=dp_t), intent(inout) :: vr_x, vr_y, vr_z, vr_max
    real (kind=dp_t), intent(inout) :: T_max, enuc_max
    real (kind=dp_t), intent(inout) :: rhovr_x, rhovr_y, rhovr_z, mass
    integer         , intent(inout) :: nzones

    !     Local variables
    integer            :: i, j, k
    real (kind=dp_t) :: velr

    if (.not. spherical == 1) then
       call bl_error("ERROR: geometry not spherical in diag")
    endif


    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)           

             if (s(i,j,k,rho_comp) > base_cutoff_density) then

                velr = u(i,j,k,1)*normal(i,j,k,1) + &
                       u(i,j,k,2)*normal(i,j,k,2) + &
                       u(i,j,k,3)*normal(i,j,k,3) + w0r(i,j,k)

                vr_max = max(vr_max,abs(velr))

                vr_x = vr_x + velr*normal(i,j,k,1)
                vr_y = vr_y + velr*normal(i,j,k,2)
                vr_z = vr_z + velr*normal(i,j,k,3)

                rhovr_x = rhovr_x + s(i,j,k,rho_comp)*velr*normal(i,j,k,1)
                rhovr_y = rhovr_y + s(i,j,k,rho_comp)*velr*normal(i,j,k,2)
                rhovr_z = rhovr_z + s(i,j,k,rho_comp)*velr*normal(i,j,k,3)

                mass = mass + s(i,j,k,rho_comp)
                nzones = nzones + 1

                T_max = max(T_max,s(i,j,k,temp_comp))
                enuc_max = max(enuc_max,rho_Hnuc(i,j,k)/s(i,j,k,rho_comp))

             endif

          enddo
       enddo
    enddo

  end subroutine diag_3d

end module diag_module
