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

  subroutine diag(time,dt,dx,s,rho0,rhoh0,p0,tempbar,gamma1bar,u,w0,normal, &
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
    real(kind=dp_t), pointer::  up(:,:,:,:)
    real(kind=dp_t), pointer::  np(:,:,:,:)
    real(kind=dp_t), pointer::  w0rp(:,:,:,:)

    type(multifab) ::    w0r_cart(mla%nlevel)

    real(kind=dp_t) :: vr_x, vr_y, vr_z, vr_max
    real(kind=dp_t) :: vr_x_local, vr_y_local, vr_z_local, vr_max_local
    integer         :: nzones, nzones_local

    integer :: lo(dm),hi(dm),ng_s,ng_u,ng_n,ng_w
    integer :: i,n
    integer :: un
    logical :: lexist

    type(bl_prof_timer), save :: bpt

    call build(bpt, "diagnostics")

    if (spherical .eq. 1) then
       do n=1,nlevs

          ! w0r_cart is w0 but onto a Cartesian grid in cell-centered as
          ! a scalar.  Since w0 is the radial expansion velocity, w0r_cart
          ! is the radial w0 in a zone
          call multifab_build(w0r_cart(n), mla%la(n),1,0)
          call setval(w0r_cart(n), ZERO, all=.true.)
       end do

       ! put w0 in Cartesian cell-centers as a scalar (the radial expansion velocity)
       call put_1d_array_on_cart(w0,w0r_cart,foextrap_comp,.true.,.false.,dx, &
            the_bc_tower%bc_tower_array,mla,normal=normal)

    endif


    ng_s = s(1)%ng
    ng_u = u(1)%ng
    ng_n = normal(1)%ng
    ng_w = w0r_cart(1)%ng

    vr_x = ZERO
    vr_y = ZERO
    vr_z = ZERO
    vr_max = ZERO

    vr_x_local = ZERO
    vr_y_local = ZERO
    vr_z_local = ZERO
    vr_max_local = ZERO

    nzones = 0
    nzones_local = 0


    do n = 1, nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          sp => dataptr(s(n) , i)
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
                          rho0(n,:),rhoh0(n,:), &
                          p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                          up(:,:,:,:),ng_u, &
                          w0(n,:), &
                          w0rp(:,:,:,1), ng_w, &
                          np(:,:,:,:),ng_n, &
                          lo,hi,vr_x_local,vr_y_local,vr_z_local,vr_max_local,nzones_local)
          end select
       end do
    end do

    ! we now have vr_n on each processor -- do a reduce
    call parallel_reduce(vr_x,vr_x_local, MPI_SUM)
    call parallel_reduce(vr_y,vr_y_local, MPI_SUM)
    call parallel_reduce(vr_z,vr_z_local, MPI_SUM)
    call parallel_reduce(nzones,nzones_local, MPI_SUM)

    call parallel_reduce(vr_max,vr_max_local,MPI_MAX)

    vr_x = vr_x/nzones
    vr_y = vr_y/nzones
    vr_z = vr_z/nzones

1000 format(1x,6(g20.10,1x))
1001 format(1x,6(a20,1x))

    if (parallel_IOProcessor()) then
       un = unit_new()
       inquire(file="wdconvect_radvel_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un, file="wdconvect_radvel_diag.out", &
               status="old", position="append")
       else
          open(unit=un, file="wdconvect_radvel_diag.out", status="new")
          write (un,1001) "time", "<vr_x>", "<vr_y>", "<vr_z>", "<vr>", "max{|vr|}"
       endif

       write (un,1000) time, vr_x, vr_y, vr_z, &
            sqrt(vr_x**2 + vr_y**2 + vr_z**2), vr_max
       close(un)
    endif

    call destroy(bpt)

  end subroutine diag

  subroutine diag_3d(time,dt,dx,s,ng_s,rho0,rhoh0,p0,tempbar,gamma1bar, &
                     u,ng_u,w0,w0r,ng_w,normal,ng_n,lo,hi,vr_x,vr_y,vr_z,vr_max,nzones)

    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp
    use network, only: nspec
    use geometry, only: spherical
    use probin_module, only: base_cutoff_density

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_u, ng_n, ng_w
    real (kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:), rhoh0(0:), &
                                         p0(0:),tempbar(0:),gamma1bar(0:)
    real (kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: w0(0:)
    real (kind=dp_t), intent(in   ) ::    w0r(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:)
    real (kind=dp_t), intent(in   ) :: normal(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)
    real (kind=dp_t), intent(in   ) :: time, dt, dx(:)
    real (kind=dp_t), intent(inout) :: vr_x, vr_y, vr_z, vr_max
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

                nzones = nzones + 1

             endif

          enddo
       enddo
    enddo

  end subroutine diag_3d

end module diag_module
