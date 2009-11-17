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

  subroutine diag(time,dt,dx,s,rho_Hnuc,rho_Hext,thermal,rho_omegadot, &
                  rho0,rhoh0,p0,tempbar, &
                  gamma1bar,div_coeff, &
                  u,w0,normal, &
                  mla,the_bc_tower)

    use bl_prof_module
    use geometry, only: dm, nlevs, spherical
    use bl_constants_module

    real(kind=dp_t), intent(in   ) :: dt,dx(:,:),time
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: rho_Hnuc(:)
    type(multifab) , intent(in   ) :: rho_Hext(:)
    type(multifab) , intent(in   ) :: thermal(:)
    type(multifab),  intent(in   ) :: rho_omegadot(:)
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
    real(kind=dp_t), pointer::  np(:,:,:,:)

    integer :: lo(dm),hi(dm),ng_s,ng_u,ng_n,ng_rhn,ng_rhe
    integer :: i,n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "diagnostics")

    ng_s = s(1)%ng
    ng_u = u(1)%ng
    ng_n = normal(1)%ng
    ng_rhn = rho_Hnuc(1)%ng
    ng_rhe = rho_Hext(1)%ng

    do n = 1, nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          sp => dataptr(s(n) , i)
          rhnp => dataptr(rho_Hnuc(n), i)
          rhep => dataptr(rho_Hext(n), i)
          up => dataptr(u(n) , i)
          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))
          print *,'DM ',dm
          select case (dm)
          case (1)
             call diag_1d(time,dt,dx(n,:), &
                          sp(:,1,1,:),ng_s, &
                          rhnp(:,1,1,1),ng_rhn, &
                          rhep(:,1,1,1),ng_rhe, &
                          rho0(n,:),rhoh0(n,:), &
                          p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                          up(:,1,1,:),ng_u, &
                          w0(n,:), &
                          lo,hi)
          case (2)
             call diag_2d(time,dt,dx(n,:), &
                          sp(:,:,1,:),ng_s, &
                          rhnp(:,:,1,1),ng_rhn, &
                          rhep(:,:,1,1),ng_rhe, &
                          rho0(n,:),rhoh0(n,:), &
                          p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                          up(:,:,1,:),ng_u, &
                          w0(n,:), &
                          lo,hi)
          case (3)
             np => dataptr(normal(n) , i)
             if (spherical .eq. 0) then
                call diag_3d(time,dt,dx(n,:), &
                             sp(:,:,:,:),ng_s, &
                             rhnp(:,:,:,1),ng_rhn, &
                             rhep(:,:,:,1),ng_rhe, &
                             rho0(n,:),rhoh0(n,:), &
                             p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                             up(:,:,:,:),ng_u, &
                             w0(n,:), &
                             np(:,:,:,:),ng_n, &
                             lo,hi)
             else
                call diag_3d_sphr(time,dt,dx(n,:), &
                                 sp(:,:,:,:),ng_s, &
                                 rhnp(:,:,:,1),ng_rhn, &
                                 rhep(:,:,:,1),ng_rhe, &
                                 rho0(1,:),rhoh0(1,:), &
                                 p0(1,:),tempbar(1,:),gamma1bar(1,:), &
                                 up(:,:,:,:),ng_u, &
                                 w0(1,:), &
                                 np(:,:,:,:),ng_n, &
                                 lo,hi)
             end if
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine diag

  subroutine diag_1d(time,dt,dx, &
                     s,ng_s, &
                     rho_Hnuc,ng_rhn, &
                     rho_Hext,ng_rhe, &
                     rho0,rhoh0,p0,tempbar,gamma1bar, &
                     u,ng_u, &
                     w0, &
                     lo,hi)

    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp
    use network, only: nspec

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_u, ng_rhn, ng_rhe
    real (kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho_Hnuc(lo(1)-ng_rhn:)
    real (kind=dp_t), intent(in   ) :: rho_Hext(lo(1)-ng_rhe:)
    real (kind=dp_t), intent(in   ) :: rho0(0:), rhoh0(0:), &
                                         p0(0:),tempbar(0:),gamma1bar(0:)
    real (kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: w0(0:)
    real (kind=dp_t), intent(in   ) :: time, dt, dx(:)

    !     Local variables
    integer            :: i

    do i = lo(1), hi(1)

       ! do diagnostics here
       !
       ! access state variables as:
       !   s(i,rho_comp)
       !
       ! access velocity components as:
       !   u(i,1)
       !

    enddo

  end subroutine diag_1d

  subroutine diag_2d(time,dt,dx, &
                     s,ng_s, &
                     rho_Hnuc,ng_rhn, &
                     rho_Hext,ng_rhe, &
                     rho0,rhoh0,p0,tempbar,gamma1bar, &
                     u,ng_u, &
                     w0, &
                     lo,hi)

    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp
    use network, only: nspec

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_u, ng_rhn, ng_rhe
    real (kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho_Hnuc(lo(1)-ng_rhn:,lo(2)-ng_rhn:)
    real (kind=dp_t), intent(in   ) :: rho_Hext(lo(1)-ng_rhe:,lo(2)-ng_rhe:)
    real (kind=dp_t), intent(in   ) :: rho0(0:), rhoh0(0:), &
                                         p0(0:),tempbar(0:),gamma1bar(0:)
    real (kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: w0(0:)
    real (kind=dp_t), intent(in   ) :: time, dt, dx(:)

    !     Local variables
    integer            :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          ! do diagnostics here
          !
          ! access state variables as:
          !   s(i,j,rho_comp)
          !
          ! access velocity components as:
          !   u(i,j,1), u(i,j,2), u(i,j,3)
          !

       enddo
    enddo

  end subroutine diag_2d

  subroutine diag_3d(time,dt,dx, &
                     s,ng_s, &
                     rho_Hnuc,ng_rhn, &
                     rho_Hext,ng_rhe, &
                     rho0,rhoh0,p0,tempbar,gamma1bar, &
                     u,ng_u,w0,normal,ng_n,lo,hi)

    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp
    use network, only: nspec

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_u, ng_n, ng_rhn, ng_rhe
    real (kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho_Hnuc(lo(1)-ng_rhn:,lo(2)-ng_rhn:,lo(3)-ng_rhn:)
    real (kind=dp_t), intent(in   ) :: rho_Hext(lo(1)-ng_rhe:,lo(2)-ng_rhe:,lo(3)-ng_rhe:)
    real (kind=dp_t), intent(in   ) :: rho0(0:), rhoh0(0:), &
                                         p0(0:),tempbar(0:),gamma1bar(0:)
    real (kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: w0(0:)
    real (kind=dp_t), intent(in   ) :: normal(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)
    real (kind=dp_t), intent(in   ) :: time, dt, dx(:)

    !     Local variables
    integer            :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! do diagnostics here
             !
             ! access state variables as:
             !   s(i,j,k,rho_comp)
             !
             ! access velocity components as:
             !   u(i,j,k,1), u(i,j,k,2), u(i,j,k,3)
             !
             ! access normal vector as:
             !   normal(i,j,k,1), normal(i,j,k,2), normal(i,j,k,3)

          enddo
       enddo
    enddo

  end subroutine diag_3d

  subroutine diag_3d_sphr(time,dt,dx, &
                          s,ng_s, &
                          rho_Hnuc,ng_rhn, &
                          rho_Hext,ng_rhe, &
                          rho0,rhoh0,p0,tempbar,gamma1bar, &
                          u,ng_u,w0,normal,ng_n,lo,hi)

    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp
    use network, only: nspec

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_u, ng_n, ng_rhn, ng_rhe
    real (kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho_Hnuc(lo(1)-ng_rhn:,lo(2)-ng_rhn:,lo(3)-ng_rhn:)
    real (kind=dp_t), intent(in   ) :: rho_Hext(lo(1)-ng_rhe:,lo(2)-ng_rhe:,lo(3)-ng_rhe:)
    real (kind=dp_t), intent(in   ) :: rho0(0:), rhoh0(0:), &
                                         p0(0:),tempbar(0:),gamma1bar(0:)
    real (kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: w0(0:)
    real (kind=dp_t), intent(in   ) :: normal(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)
    real (kind=dp_t), intent(in   ) :: time, dt, dx(:)

    !     Local variables
    integer            :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! do diagnostics here
             !
             ! access state variables as:
             !   s(i,j,k,rho_comp)
             !
             ! access velocity components as:
             !   u(i,j,k,1), u(i,j,k,2), u(i,j,k,3)
             !
             ! access normal vector as:
             !   normal(i,j,k,1), normal(i,j,k,2), normal(i,j,k,3)

          enddo
       enddo
    enddo

  end subroutine diag_3d_sphr

end module diag_module
