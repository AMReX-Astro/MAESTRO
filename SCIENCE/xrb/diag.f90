! This is an interface for doing runtime diagnostics on the state.
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
    use geometry, only: dm, nlevs
    use bl_constants_module

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

    real(kind=dp_t) :: T_max, T_max_local
    real(kind=dp_t) :: enuc_max, enuc_max_local

    integer :: lo(dm),hi(dm),ng_s,ng_u,ng_n,ng_rhn,ng_rhe
    integer :: i,n
    integer :: un, un2
    logical :: lexist

    logical, save :: firstCall = .true.

    type(bl_prof_timer), save :: bpt

    call build(bpt, "diagnostics")

    ng_s = s(1)%ng
    ng_u = u(1)%ng
    ng_n = normal(1)%ng
    ng_rhn = rho_Hnuc(1)%ng
    ng_rhe = rho_Hext(1)%ng

    T_max          = ZERO
    T_max_local    = ZERO
    enuc_max       = ZERO
    enuc_max_local = ZERO

    do n = 1, nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          sp => dataptr(s(n) , i)
          rhnp => dataptr(rho_Hnuc(n), i)
          rhep => dataptr(rho_Hext(n), i)
          up => dataptr(u(n) , i)
          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))
          select case (dm)
          case (2)
             call diag_2d(time,dt,dx(n,:), &
                          sp(:,:,1,:),ng_s, &
                          rhnp(:,:,1,1),ng_rhn, &
                          rhep(:,:,1,1),ng_rhe, &
                          rho0(n,:),rhoh0(n,:), &
                          p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                          up(:,:,1,:),ng_u, &
                          w0(n,:), &
                          lo,hi, &
                          T_max_local, enuc_max_local)
          case (3)
             np => dataptr(normal(n) , i)
             call diag_3d(time,dt,dx(n,:), &
                          sp(:,:,:,:),ng_s, &
                          rhnp(:,:,:,1),ng_rhn, &
                          rhep(:,:,:,1),ng_rhe, &
                          rho0(n,:),rhoh0(n,:), &
                          p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                          up(:,:,:,:),ng_u, &
                          w0(n,:), &
                          np(:,:,:,:),ng_n, &
                          lo,hi, &
                          T_max, enuc_max_local)
          end select
       end do
    end do

    call parallel_reduce(T_max, T_max_local, MPI_MAX, &
                         proc = parallel_IOProcessorNode())

    call parallel_reduce(enuc_max, enuc_max_local, MPI_MAX, &
                         proc = parallel_IOProcessorNode())

1000 format(1x,10(g20.10,1x))
1001 format("#",10(a20,1x))

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

       ! print the headers
       if(firstCall) then
          write(un ,    *) " "
          write(un , 1001) "time", "max{T}"
          write(un2, 1001) "time", "max{enuc}"

          firstCall = .false.
       endif

       ! print the data
       write(un , 1000) time, T_max
       write(un2, 1000) time, enuc_max

       close(un )
       close(un2)
    endif

    call destroy(bpt)

  end subroutine diag

  subroutine diag_2d(time,dt,dx, &
                     s,ng_s, &
                     rho_Hnuc,ng_rhn, &
                     rho_Hext,ng_rhe, &
                     rho0,rhoh0,p0,tempbar,gamma1bar, &
                     u,ng_u, &
                     w0, &
                     lo,hi, &
                     T_max, enuc_max)

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
    real (kind=dp_t), intent(inout) :: T_max, enuc_max

    !     Local variables
    integer            :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          T_max = max(T_max,s(i,j,temp_comp))

          enuc_max = max(enuc_max,rho_Hnuc(i,j)/s(i,j,rho_comp))

       enddo
    enddo

  end subroutine diag_2d

  subroutine diag_3d(time,dt,dx, &
                     s,ng_s, &
                     rho_Hnuc,ng_rhn, &
                     rho_Hext,ng_rhe, &
                     rho0,rhoh0,p0,tempbar,gamma1bar, &
                     u,ng_u,w0,normal,ng_n,lo,hi,T_max,enuc_max)

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
    real (kind=dp_t), intent(inout) :: T_max, enuc_max

    !     Local variables
    integer            :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             T_max = max(T_max,s(i,j,k,temp_comp))

             enuc_max = max(enuc_max,rho_Hnuc(i,j,k)/s(i,j,k,rho_comp))

          enddo
       enddo
    enddo

  end subroutine diag_3d

end module diag_module
