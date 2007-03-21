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
      real(kind=dp_t), intent(  out) :: phibar(0:,:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      real(kind=dp_t), pointer:: pp(:,:,:,:)
      integer :: lo(phi%dim),hi(phi%dim),ng,dm,nr
      integer :: i,k
      integer        , allocatable :: npts_grid(:), npts_proc(:), npts_tot(:)
      real(kind=dp_t), allocatable :: vol_grid(:), vol_proc(:), vol_tot(:)

      dm = phi%dim
      ng = phi%ng
      nr = size(phibar,dim=1)

      if (spherical .eq. 1) then
        allocate(vol_grid(0:nr-1),vol_proc(0:nr-1),vol_tot(0:nr-1))
        vol_proc(:) = ZERO
        vol_tot(:)  = ZERO
      else
        allocate(npts_grid(0:nr-1),npts_proc(0:nr-1),npts_tot(0:nr-1))
        npts_proc(:) = 0
        npts_tot(:)  = 0
      end if

      phibar(:,:) = ZERO

      do i = 1, phi%nboxes
         if ( multifab_remote(phi, i) ) cycle
         pp => dataptr(phi, i)
         lo =  lwb(get_box(phi, i))
         hi =  upb(get_box(phi, i))
         select case (dm)
            case (2)
              npts_grid(:) = 0
              call average_2d(pp(:,:,1,:),phibar(lo(2):,:),lo,hi,ng,npts_grid(lo(2):))
              npts_proc(lo(2):hi(2)) = npts_proc(lo(2):hi(2)) + npts_grid(lo(2):hi(2))
            case (3)
              if (spherical .eq. 1) then
                vol_grid(:) = ZERO
                call average_3d_sphr(pp(:,:,:,:),phibar,lo,hi,ng,dx,vol_grid)
                vol_proc = vol_proc + vol_grid
              else
                npts_grid(:) = 0
                call average_3d(pp(:,:,:,:),phibar(lo(3):,:),lo,hi,ng,npts_grid(lo(3):))
                npts_proc(lo(3):hi(3)) = npts_proc(lo(3):hi(3)) + npts_grid(lo(3):hi(3))
              end if
         end select
      end do

      if (dm .eq. 2 .or. (dm.eq.3.and.spherical.eq.0)) then
        call parallel_reduce(npts_tot,npts_proc,MPI_SUM)
        do k = 0,nr-1
          phibar(k,:) = phibar(k,:) / dble(npts_tot(k))
        end do
        deallocate(npts_grid,npts_proc,npts_tot)
      else
        do k = 0,nr-1
          call parallel_reduce(vol_tot(k),vol_proc(k),MPI_SUM)
          if (vol_tot(k) .gt. ZERO) then
            phibar(k,:) = phibar(k,:) / vol_tot(k)
          else
            phibar(k,:) = ZERO
          end if
        end do
        deallocate(vol_grid,vol_proc,vol_tot)
      end if

   end subroutine average

   subroutine average_2d (phi,phibar,lo,hi,ng,npts)

      integer         , intent(in   ) :: lo(:), hi(:), ng
      real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,:)
      real (kind=dp_t), intent(  out) :: phibar(lo(2):,:)
      integer         , intent(  out) :: npts(lo(2):)

!     Local variables
      integer          :: i, j, n


      do n = 1,size(phibar,dim=2)
      do j = lo(2),hi(2)
        npts(j) = hi(1)-lo(1)+1
        do i = lo(1),hi(1)
          phibar(j,n) = phibar(j,n) + phi(i,j,n)
        end do
      end do
      end do
 
   end subroutine average_2d

   subroutine average_3d (phi,phibar,lo,hi,ng,npts)

      integer         , intent(in   ) :: lo(:), hi(:), ng
      real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind=dp_t), intent(  out) :: phibar(lo(3):,:)
      integer         , intent(  out) :: npts(lo(3):)

!     Local variables
      integer          :: i, j, k, n

      do n = 1,size(phibar,dim=2)
      do k = lo(3),hi(3)
        npts(k) = (hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          phibar(k,n) = phibar(k,n) + phi(i,j,k,n)
        end do
        end do
      end do
      end do
 
   end subroutine average_3d

   subroutine average_3d_sphr (phi,phibar,lo,hi,ng,dx,sum)

      integer         , intent(in   ) :: lo(:),hi(:),ng
      real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind=dp_t), intent(  out) :: phibar(0:,:)
      real (kind=dp_t), intent(in   ) :: dx(:)
      real (kind=dp_t), intent(  out) :: sum(:)

!     Local variables
      integer                       :: i, j, k, n, index
      integer                       :: nr, nc
      real (kind=dp_t)              :: x,y,z,radius,vol
      real (kind=dp_t), parameter   :: fourthirdspi = 4.18879020478639098400_dp_t

      nr = size(phibar,dim=1)
      nc = size(phibar,dim=2)

      do k = lo(3),hi(3)
        z = (dble(k)-HALF)*dx(3) - center(3)
        do j = lo(2),hi(2)
          y = (dble(j)-HALF)*dx(2) - center(2)
          do i = lo(3),hi(3)
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
 
   end subroutine average_3d_sphr

end module average_module
