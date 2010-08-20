module addw0_module

  use bl_types
  use multifab_module

  implicit none

  private

  public :: addw0, addw0_3d_sphr

contains

  subroutine addw0(umac,w0,w0mac,mult)

    use bl_prof_module
    use create_umac_grown_module
    use geometry, only: spherical, dm, nlevs
    
    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    real(kind=dp_t), intent(in   ) :: mult

    ! Local variables
    integer :: i,lo(dm),hi(dm),n,ng_um,ng_w0
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: w0xp(:,:,:,:)
    real(kind=dp_t), pointer :: w0yp(:,:,:,:)
    real(kind=dp_t), pointer :: w0zp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "addw0")

    ng_um = nghost(umac(1,1)) ! here we are assuming all components have the 
                              ! same # of ghostcells

    do n = 1, nlevs
       do i = 1, nboxes(umac(n,dm))
          if ( multifab_remote(umac(n,dm), i) ) cycle
          wmp  => dataptr(umac(n,dm), i)
          lo =  lwb(get_box(umac(n,dm), i))
          hi =  upb(get_box(umac(n,dm), i))
          select case(dm)
          case(1)
             call addw0_1d(wmp(:,1,1,1),ng_um,w0(n,:),lo,hi,mult)
          case(2)
             call addw0_2d(wmp(:,:,1,1),ng_um,w0(n,:),lo,hi,mult)
          case(3)
             if (spherical .eq. 0) then
                call addw0_3d(wmp(:,:,:,1),ng_um,w0(n,:),lo,hi,mult)
             else
                ump  => dataptr(umac(n,1), i)
                vmp  => dataptr(umac(n,2), i)
                wmp  => dataptr(umac(n,3), i)
                w0xp  => dataptr(w0mac(n,1), i)
                w0yp  => dataptr(w0mac(n,2), i)
                w0zp  => dataptr(w0mac(n,3), i)
                ng_w0 = nghost(w0mac(n,1))
                call addw0_3d_sphr(ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),ng_um, &
                                   w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1),ng_w0, &
                                   lo,hi,mult)
             end if
          end select
       end do
    end do

    if (nlevs .gt. 1) then
       do n=2,nlevs
          call create_umac_grown(n,umac(n,:),umac(n-1,:))
       end do
    else
       do i=1,dm
          call multifab_fill_boundary(umac(1,i))
       enddo
    end if

    call destroy(bpt)

  end subroutine addw0

  subroutine addw0_1d(umac,ng_um,w0,lo,hi,mult)

    integer        , intent(in   ) :: lo(:),hi(:),ng_um
    real(kind=dp_t), intent(inout) :: umac(lo(1)-ng_um:)
    real(kind=dp_t), intent(in   ) ::   w0(0:)
    real(kind=dp_t), intent(in   ) :: mult

    integer :: i

    do i = lo(1),hi(1)+1
       umac(i) = umac(i) + mult * w0(i)
    end do

  end subroutine addw0_1d

  subroutine addw0_2d(vmac,ng_um,w0,lo,hi,mult)

    integer        , intent(in   ) :: lo(:),hi(:),ng_um
    real(kind=dp_t), intent(inout) :: vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::   w0(0:)
    real(kind=dp_t), intent(in   ) :: mult

    integer :: i,j

    do j = lo(2),hi(2)+1
       do i = lo(1)-1,hi(1)+1
          vmac(i,j) = vmac(i,j) + mult * w0(j)
       end do
    end do

  end subroutine addw0_2d

  subroutine addw0_3d(wmac,ng_um,w0,lo,hi,mult)

    integer        , intent(in   ) :: lo(:),hi(:),ng_um
    real(kind=dp_t), intent(inout) :: wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::   w0(0:)
    real(kind=dp_t), intent(in   ) :: mult

    integer :: i,j,k

    do k = lo(3),hi(3)+1
       do j = lo(2)-1,hi(2)+1
          do i = lo(1)-1,hi(1)+1
             wmac(i,j,k) = wmac(i,j,k) + mult * w0(k)
          end do
       end do
    end do

  end subroutine addw0_3d

  subroutine addw0_3d_sphr(umac,vmac,wmac,ng_um,w0macx,w0macy,w0macz,ng_w0,lo,hi,mult)

    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_um,ng_w0
    real(kind=dp_t), intent(inout) ::    umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(inout) ::    vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(inout) ::    wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) ::  w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) ::  w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: mult

    integer :: i,j,k

    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             umac(i,j,k) = umac(i,j,k) + mult * w0macx(i,j,k)
          end do
       end do
    end do
    !$OMP END DO NOWAIT
    !$OMP DO
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
             vmac(i,j,k) = vmac(i,j,k) + mult * w0macy(i,j,k)
          end do
       end do
    end do
    !$OMP END DO NOWAIT
    !$OMP DO
    do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             wmac(i,j,k) = wmac(i,j,k) + mult * w0macz(i,j,k)
          end do
       end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine addw0_3d_sphr

end module addw0_module
