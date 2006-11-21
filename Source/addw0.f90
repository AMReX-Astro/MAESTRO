module addw0_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use geometry
  use fill_3d_module

  implicit none

contains

      subroutine addw0(umac,normal,w0,dx,mult)

      type(multifab) , intent(inout) :: umac(:)
      type(multifab) , intent(in   ) :: normal
      real(kind=dp_t), intent(in   ) :: w0(:)
      real(kind=dp_t), intent(in   ) :: dx(:),mult

      ! Local variables
      integer :: i,lo(umac(1)%dim),hi(umac(1)%dim),dm
      real(kind=dp_t), pointer :: ump(:,:,:,:)
      real(kind=dp_t), pointer :: vmp(:,:,:,:)
      real(kind=dp_t), pointer :: wmp(:,:,:,:)
      real(kind=dp_t), pointer ::  np(:,:,:,:)

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
           if (spherical .eq. 0) then
             call addw0_3d(wmp(:,:,:,1),w0,lo,hi,mult)
           else
             ump  => dataptr(umac(1), i)
             vmp  => dataptr(umac(2), i)
              np  => dataptr(normal, i)
             call addw0_3d_sphr(ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),np(:,:,:,:),w0,dx,mult)
           end if
         end select
      end do

      end subroutine addw0

      subroutine addw0_2d(vmac,w0,lo,hi,mult)

      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(inout) :: vmac(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(in   ) ::   w0(0:)
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
      real(kind=dp_t), intent(in   ) ::   w0(0:)
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

      subroutine addw0_3d_sphr(umac,vmac,wmac,normal,w0,dx,mult)

      real(kind=dp_t), intent(inout) ::   umac(0:,0:,0:)
      real(kind=dp_t), intent(inout) ::   vmac(0:,0:,0:)
      real(kind=dp_t), intent(inout) ::   wmac(0:,0:,0:)
      real(kind=dp_t), intent(in   ) :: normal(0:,0:,0:,:)
      real(kind=dp_t), intent(in   ) ::   w0(:)
      real(kind=dp_t), intent(in   ) :: dx(:),mult

      integer :: i,j,k
      integer :: nx,ny,nz,nr
      real(kind=dp_t), allocatable :: w0_rad(:)
      real(kind=dp_t), allocatable :: w0_cart(:,:,:)

      nx = size(vmac,dim=1) - 2
      ny = size(wmac,dim=2) - 2
      nz = size(umac,dim=3) - 2
      nr = size(w0,dim=1) - 1

      allocate(w0_rad (nr))
      allocate(w0_cart(0:nx+1,0:ny+1,0:nz+1))

      ! Put w0 on centers of radial grid first
      do k = 1,nr
        w0_rad(k) = HALF * (w0(k) + w0(k+1))
      end do

      ! Then put w0 on centers of Cartesian grid
      call fill_3d_data(w0_cart,w0_rad,dx,1)

      do k = 1,nz
      do j = 1,ny
         w0_cart(   0,j,k) = w0_cart( 1,j,k)
         w0_cart(nx+1,j,k) = w0_cart(nx,j,k)
      end do
      end do
      do k = 1,nz
      do i = 1,nx
         w0_cart(i,   0,k) = w0_cart(i, 1,k)
         w0_cart(i,ny+1,k) = w0_cart(i,ny,k)
      end do
      end do
      do j = 1,ny
      do i = 1,nx
         w0_cart(i,j,   0) = w0_cart(i,j, 1)
         w0_cart(i,j,nz+1) = w0_cart(i,j,nz)
      end do
      end do

      do k = 1,nz
      do j = 1,ny
      do i = 1,nx+1
         umac(i,j,k) = umac(i,j,k) + mult * &
                       HALF * ( w0_cart(i  ,j,k) * normal(i  ,j,k,1) &
                               +w0_cart(i-1,j,k) * normal(i-1,j,k,1) )
      end do
      end do
      end do
      do k = 1,nz
      do j = 1,ny+1
      do i = 1,nx
         vmac(i,j,k) = vmac(i,j,k) + mult * &
                       HALF * ( w0_cart(i,j  ,k) * normal(i,j  ,k,2) &
                               +w0_cart(i,j-1,k) * normal(i,j-1,k,2) )
      end do
      end do
      end do
      do k = 1,nz+1
      do j = 1,ny
      do i = 1,nx
         wmac(i,j,k) = wmac(i,j,k) + mult * &
                       HALF * ( w0_cart(i,j  ,k) * normal(i,j,k  ,3) &
                               +w0_cart(i,j-1,k) * normal(i,j,k-1,3) )
      end do
      end do
      end do

      deallocate(w0_rad,w0_cart)

      end subroutine addw0_3d_sphr

end module addw0_module
