module average_module

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_module

  implicit none

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine average (phi,phibar)

      type(multifab) , intent(inout) :: phi
      real(kind=dp_t), intent(  out) :: phibar(:,:)

      real(kind=dp_t), pointer:: pp(:,:,:,:)
      integer :: lo(phi%dim),hi(phi%dim),ng,dm
      integer :: i

      dm = phi%dim
      ng = phi%ng

      do i = 1, phi%nboxes
         if ( multifab_remote(phi, i) ) cycle
         pp => dataptr(phi, i)
         lo =  lwb(get_box(phi, i))
         hi =  upb(get_box(phi, i))
         select case (dm)
            case (2)
              call average_2d(pp(:,:,1,:),phibar,lo,hi,ng)
            case (3)
              call average_3d(pp(:,:,:,:),phibar,lo,hi,ng)
         end select
      end do

   end subroutine average

   subroutine average_2d (phi,phibar,lo,hi,ng)

      integer         , intent(in   ) :: lo(:), hi(:), ng
      real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,:)
      real (kind=dp_t), intent(  out) :: phibar(lo(2):,:)

!     Local variables
      integer          :: i, j, n
      real (kind=dp_t) :: denom

      denom = dble(hi(1)-lo(1)+1)

      do n = 1,size(phibar,dim=2)
      do j = lo(2),hi(2)
        phibar(j,n) = zero
        do i = lo(1),hi(1)
          phibar(j,n) = phibar(j,n) + phi(i,j,n)
        end do
        phibar(j,n) = phibar(j,n) / denom
      end do
      end do
 
   end subroutine average_2d

   subroutine average_3d (phi,phibar,lo,hi,ng)

      integer         , intent(in   ) :: lo(:), hi(:), ng
      real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind=dp_t), intent(  out) :: phibar(lo(3):,:)

!     Local variables
      integer          :: i, j, k, n
      real (kind=dp_t) :: denom

      denom = dble((hi(1)-lo(1)+1)*(hi(2)-lo(2)+1))

      do n = 1,size(phibar,dim=2)
      do k = lo(3),hi(3)
        phibar(k,n) = zero
        do j = lo(2),hi(3)
        do i = lo(1),hi(1)
          phibar(k,n) = phibar(k,n) + phi(i,j,k,n)
        end do
        end do
        phibar(k,n) = phibar(k,n) / denom
      end do
      end do
 
   end subroutine average_3d

end module average_module
