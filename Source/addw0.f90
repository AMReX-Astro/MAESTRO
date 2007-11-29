module addw0_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use geometry
  use fill_3d_module

  implicit none

  private
  public :: addw0, addw0_3d_sphr

contains

      subroutine addw0(umac,w0,w0_cart,dx,mult)

      type(multifab) , intent(inout) :: umac(:)
      real(kind=dp_t), intent(in   ) :: w0(0:)
      type(multifab) , intent(in   ) :: w0_cart
      real(kind=dp_t), intent(in   ) :: dx(:),mult

      ! Local variables
      integer :: i,lo(umac(1)%dim),hi(umac(1)%dim),dm
      real(kind=dp_t), pointer :: ump(:,:,:,:)
      real(kind=dp_t), pointer :: vmp(:,:,:,:)
      real(kind=dp_t), pointer :: wmp(:,:,:,:)
      real(kind=dp_t), pointer :: w0p(:,:,:,:)

      dm = umac(1)%dim

      do i = 1, umac(dm)%nboxes
         if ( multifab_remote(umac(dm), i) ) cycle
         wmp  => dataptr(umac(dm), i)
         lo =  lwb(get_box(umac(dm), i))
         hi =  upb(get_box(umac(dm), i))
         select case(dm)
         case(2)
           call addw0_2d(wmp(:,:,1,1),w0,lo,hi,dx,mult)
         case(3)
           if (spherical .eq. 0) then
             call addw0_3d(wmp(:,:,:,1),w0,lo,hi,dx,mult)
           else
             ump  => dataptr(umac(1), i)
             vmp  => dataptr(umac(2), i)
             wmp  => dataptr(umac(3), i)
             w0p  => dataptr(w0_cart, i)
               lo =  lwb(get_box(w0_cart, i))
               hi =  upb(get_box(w0_cart, i))
             call addw0_3d_sphr(ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),w0p(:,:,:,:),lo,hi,mult)
           end if
         end select
      end do

      end subroutine addw0

      subroutine addw0_2d(vmac,w0,lo,hi,dx,mult)

      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: vmac(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(in   ) ::   w0(0:)
      real(kind=dp_t), intent(in   ) :: mult
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j
      integer :: rr

      rr = int( dx(2) / dr + 1.d-12)

      do j = lo(2)  ,hi(2)
      do i = lo(1)-1,hi(1)+1
         vmac(i,j) = vmac(i,j) + mult * w0(rr*j)
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
      integer :: rr

      rr = int( dx(3) / dr + 1.d-12)

      do k = lo(3),hi(3)
      do j = lo(2)-1,hi(2)+1
      do i = lo(1)-1,hi(1)+1
         wmac(i,j,k) = wmac(i,j,k) + mult * w0(rr*k)
      end do
      end do
      end do

      end subroutine addw0_3d

      subroutine addw0_3d_sphr(umac,vmac,wmac,w0_cart,lo,hi,mult)

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
