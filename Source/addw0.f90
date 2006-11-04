module addw0_module

  use bl_types
  use bl_constants_module
  use multifab_module

  implicit none

contains

      subroutine addw0(umac,w0,mult)

      type(multifab) , intent(inout) :: umac(:)
      real(kind=dp_t), intent(in   ) ::   w0(:)
      real(kind=dp_t), intent(in   ) :: mult

      integer :: i,lo(umac(1)%dim),hi(umac(1)%dim),dm
      real(kind=dp_t), pointer :: wmp(:,:,:,:)

      dm = umac(1)%dim

      do i = 1, umac(dm)%nboxes
         if ( multifab_remote(umac(dm), i) ) cycle
         wmp  => dataptr(umac(dm), i)
         lo =  lwb(get_box(umac(dm), i))
         hi =  upb(get_box(umac(dm), i))
         select case(dm)
         case(2)
           call addw0_2d(wmp(:,:,1,1),w0,lo,hi,mult)
         case(3)
           call addw0_3d(wmp(:,:,:,1),w0,lo,hi,mult)
         end select
      end do

      end subroutine addw0

      subroutine addw0_2d(vmac,w0,lo,hi,mult)

      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: vmac(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(in   ) ::   w0(          lo(2)   :)
      real(kind=dp_t), intent(in   ) :: mult

      integer :: i,j

      do j = lo(2),hi(2)
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

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
         wmac(i,j,k) = wmac(i,j,k) + mult * w0(k)
      end do
      end do
      end do

      end subroutine addw0_3d


end module addw0_module
