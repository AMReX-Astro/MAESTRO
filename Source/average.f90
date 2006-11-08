module average_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use geometry

  implicit none

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine average (phi,phibar,dx)

      type(multifab) , intent(inout) :: phi
      real(kind=dp_t), intent(  out) :: phibar(:,:)
      real(kind=dp_t), intent(in   ) :: dx(:)

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
              if (spherical .eq. 1) then
                call average_3d_sphr(pp(:,:,:,:),phibar,ng,dx)
              else
                call average_3d(pp(:,:,:,:),phibar,lo,hi,ng)
              end if
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

   subroutine average_3d_sphr (phi,phibar,ng,dx)

      integer         , intent(in   ) :: ng
      real (kind=dp_t), intent(in   ) :: phi(1-ng:,1-ng:,1-ng:,:)
      real (kind=dp_t), intent(  out) :: phibar(:,:)
      real (kind=dp_t), intent(in   ) :: dx(:)

!     Local variables
      integer                       :: i, j, k, n, index
      integer                       :: nx, ny, nz, nr, nc
      real (kind=dp_t)              :: x,y,z,radius,vol
      real (kind=dp_t), allocatable :: sum(:)
      real (kind=dp_t), parameter   :: fourthirdspi = 4.18879020478639098400_dp_t

      nx = size(phi,dim=1)-2*ng
      ny = size(phi,dim=2)-2*ng
      nz = size(phi,dim=3)-2*ng

      nr = size(phibar,dim=1)
      nc = size(phibar,dim=2)

      allocate(sum(nr))

      phibar = ZERO
      sum    = ZERO

      do k = 1,nz
        z = (dble(k)-HALF)*dx(3) - center(3)
        do j = 1,nz
          y = (dble(j)-HALF)*dx(2) - center(2)
          do i = 1,nz
            x = (dble(i)-HALF)*dx(1) - center(1)

            radius = sqrt(x**2 + y**2 + z**2)
            index = radius / dr + 1

            if (index .lt. 1 .or. index .gt. nr) then
              print *,'RADIUS ',radius
              print *,'BOGUS INDEX ',index
              print *,'NOT IN RANGE 0 TO ',nr
              print *,'I J K ',i,j,k
              print *,'X Y Z ',x,y,z
              stop
            end if

            vol = fourthirdspi * (zl(index+1)**3  - zl(index)**3)

            do n = 1,nc
              phibar(index,n) = phibar(index,n) + vol * phi(i,j,k,n)
            end do

            sum(index) = sum(index) + vol

          end do
        end do
      end do

      do n = 1,nc
      do i = 1,nr-1
         phibar(i,n) = phibar(i,n) / sum(i)
      end do
      end do

      deallocate(sum)
 
   end subroutine average_3d_sphr

end module average_module
