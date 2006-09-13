module addw0_module

  use bl_types
  use bl_constants_module

  implicit none

contains

      subroutine addw0_2d(vmac,w0,lo,hi,mult)

      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: vmac(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(in   ) ::   w0(          lo(2)   :)
      real(kind=dp_t), intent(in   ) :: mult

      integer :: i,j

      do j = lo(2),hi(2)+1
      do i = lo(1),hi(1)
         vmac(i,j) = vmac(i,j) + mult * w0(j)
      end do
      end do

      end subroutine addw0_2d

      subroutine addw0_3d(wmac,w0,lo,hi,mult)

      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: wmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
      real(kind=dp_t), intent(in   ) ::   w0(                  lo(3)  :)
      real(kind=dp_t), intent(in   ) :: mult

      integer :: i,j,k

      do k = lo(3),hi(3)+1
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
         wmac(i,j,k) = wmac(i,j,k) + mult * w0(k)
      end do
      end do
      end do

      end subroutine addw0_3d


end module addw0_module
