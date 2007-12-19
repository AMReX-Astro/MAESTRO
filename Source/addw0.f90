module addw0_module

  use bl_types
  use multifab_module

  implicit none

  private

  public :: addw0, addw0_3d_sphr

contains

  subroutine addw0(nlevs,umac,w0,w0_cart,dx,mult)

    use geometry

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0_cart(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),mult

    ! Local variables
    integer :: i,lo(umac(1,1)%dim),hi(umac(1,1)%dim),dm,n
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: w0p(:,:,:,:)

    dm = umac(1,1)%dim

    do n = 1, nlevs
       do i = 1, umac(n,dm)%nboxes
          if ( multifab_remote(umac(n,dm), i) ) cycle
          wmp  => dataptr(umac(n,dm), i)
          lo =  lwb(get_box(umac(n,dm), i))
          hi =  upb(get_box(umac(n,dm), i))
          select case(dm)
          case(2)
             call addw0_2d(wmp(:,:,1,1),w0(n,:),lo,hi,dx(n,:),mult)
          case(3)
             if (spherical .eq. 0) then
                call addw0_3d(wmp(:,:,:,1),w0(n,:),lo,hi,dx(n,:),mult)
             else
                ump  => dataptr(umac(n,1), i)
                vmp  => dataptr(umac(n,2), i)
                wmp  => dataptr(umac(n,3), i)
                w0p  => dataptr(w0_cart(n), i)
                lo =  lwb(get_box(w0_cart(n), i))
                hi =  upb(get_box(w0_cart(n), i))
                call addw0_3d_sphr(ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),w0p(:,:,:,:), &
                                   lo,hi,mult)
             end if
          end select
       end do
    end do

  end subroutine addw0

  subroutine addw0_2d(vmac,w0,lo,hi,dx,mult)

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: vmac(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(in   ) ::   w0(0:)
    real(kind=dp_t), intent(in   ) :: mult
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer :: i,j

    do j = lo(2)  ,hi(2)
       do i = lo(1)-1,hi(1)+1
          vmac(i,j) = vmac(i,j) + mult * w0(j)
       end do
    end do

  end subroutine addw0_2d

  subroutine addw0_3d(wmac,w0,lo,hi,dx,mult)

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) :: wmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::   w0(0:)
    real(kind=dp_t), intent(in   ) :: mult
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer :: i,j,k

    do k = lo(3),hi(3)
       do j = lo(2)-1,hi(2)+1
          do i = lo(1)-1,hi(1)+1
             wmac(i,j,k) = wmac(i,j,k) + mult * w0(k)
          end do
       end do
    end do

  end subroutine addw0_3d

  subroutine addw0_3d_sphr(umac,vmac,wmac,w0_cart,lo,hi,mult)

    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:)
    real(kind=dp_t), intent(inout) ::   umac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(inout) ::   vmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(inout) ::   wmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: w0_cart(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)
    real(kind=dp_t), intent(in   ) :: mult

    integer :: i,j,k

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             umac(i,j,k) = umac(i,j,k) + mult * &
                  HALF * ( w0_cart(i  ,j,k,1) +w0_cart(i-1,j,k,1) )
          end do
       end do
    end do
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
             vmac(i,j,k) = vmac(i,j,k) + mult * &
                  HALF * ( w0_cart(i,j  ,k,2) + w0_cart(i,j-1,k,2) )
          end do
       end do
    end do
    do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             wmac(i,j,k) = wmac(i,j,k) + mult * &
                  HALF * ( w0_cart(i,j,k  ,3) + w0_cart(i,j,k-1,3) )
          end do
       end do
    end do

  end subroutine addw0_3d_sphr

end module addw0_module
