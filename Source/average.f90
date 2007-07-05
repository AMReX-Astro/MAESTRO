module average_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use geometry

  implicit none

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine average (phi,phibar,dx,comp,ncomp)

      integer        , intent(in   ) :: comp,ncomp
      type(multifab) , intent(inout) :: phi(:)
      real(kind=dp_t), intent(inout) :: phibar(0:,:)
      real(kind=dp_t), intent(in   ) :: dx(:,:)

      real(kind=dp_t), pointer:: pp(:,:,:,:)
      integer :: lo(phi(1)%dim),hi(phi(1)%dim),ng,dm,nr
      integer :: i,k,n,nlevs
      real(kind=dp_t), allocatable :: vol_grid(:), vol_proc(:), vol_tot(:), phibar_proc(:,:)

      dm = phi(1)%dim
      ng = phi(1)%ng
      nr = size(phibar,dim=1)
      nlevs = size(dx,dim=1)

      if (spherical .eq. 1) then
        allocate(vol_grid(0:nr-1))
      end if
      allocate(vol_proc(0:nr-1),vol_tot(0:nr-1))
      allocate(phibar_proc(0:nr-1,ncomp))

      vol_proc(:) = ZERO
      vol_tot(:)  = ZERO

      phibar(:,:) = ZERO
      phibar_proc(:,:) = ZERO

      do n = 1, nlevs
       do i = 1, phi(n)%nboxes
         if ( multifab_remote(phi(n), i) ) cycle
         pp => dataptr(phi(n), i)
         lo =  lwb(get_box(phi(n), i))
         hi =  upb(get_box(phi(n), i))
         select case (dm)
            case (2)
              call average_2d(pp(:,:,1,:),phibar_proc,lo,hi,ng,comp,ncomp,dx(n,:))
              vol_proc(lo(2):hi(2)) = vol_proc(lo(2):hi(2)) + (hi(1)-lo(1)+1)*dx(n,1)*dx(n,2)
            case (3)
              if (spherical .eq. 1) then
                vol_grid(:) = ZERO
                call average_3d_sphr(pp(:,:,:,:),phibar_proc,lo,hi,ng,dx(n,:),vol_grid,comp,ncomp)
                vol_proc = vol_proc + vol_grid
              else
                call average_3d(pp(:,:,:,:),phibar_proc,lo,hi,ng,comp,ncomp,dx(n,:))
                vol_proc(lo(3):hi(3)) = vol_proc(lo(3):hi(3)) + (hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)*dx(n,1)*dx(n,2)*dx(n,3)
              end if
         end select
       end do
      end do

      if (dm .eq. 2 .or. (dm.eq.3.and.spherical.eq.0)) then
        do k = 0,nr-1
           call parallel_reduce(vol_tot(k),  vol_proc(k),      MPI_SUM)
           call parallel_reduce(phibar(k,:), phibar_proc(k,:), MPI_SUM)

           phibar(k,:) = phibar(k,:) / vol_tot(k)
        end do

      else
        do k = 0,nr-1
          call parallel_reduce(vol_tot(k),  vol_proc(k),      MPI_SUM)
          call parallel_reduce(phibar(k,:), phibar_proc(k,:), MPI_SUM)

          if (vol_tot(k) .gt. ZERO) then
            phibar(k,:) = phibar(k,:) / vol_tot(k)
          else
            phibar(k,:) = ZERO
          end if

        end do
        deallocate(vol_grid)

      end if
      deallocate(vol_proc,vol_tot,phibar_proc)

   end subroutine average

   subroutine average_2d (phi,phibar,lo,hi,ng,comp,ncomp,dx)

      integer         , intent(in   ) :: lo(:), hi(:), ng, comp, ncomp
      real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,:)
      real (kind=dp_t), intent(inout) :: phibar(0:,:)
      real (kind=dp_t), intent(in   ) :: dx(:)

!     Local variables
      integer          :: i, j, n
      real (kind=dp_t) :: vol

      vol = dx(1)*dx(2)

      do n = comp,comp+ncomp-1
      do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          phibar(j,n) = phibar(j,n) + phi(i,j,n)
        end do
        phibar(j,n) = phibar(j,n) * vol
      end do
      end do
 
   end subroutine average_2d

   subroutine average_3d (phi,phibar,lo,hi,ng,comp,ncomp,dx)

      integer         , intent(in   ) :: lo(:), hi(:), ng, comp, ncomp
      real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind=dp_t), intent(inout) :: phibar(0:,:)
      real (kind=dp_t), intent(in   ) :: dx(:)

!     Local variables
      integer          :: i, j, k, n
      real (kind=dp_t) :: vol

      vol = dx(1)*dx(2)*dx(3)

      do n = comp,comp+ncomp-1
      do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          phibar(k,n) = phibar(k,n) + phi(i,j,k,n)
        end do
        end do
        phibar(k,n) = phibar(k,n) * vol
      end do
      end do

 
   end subroutine average_3d

   subroutine average_3d_sphr (phi,phibar,lo,hi,ng,dx,sum,comp,ncomp)

      integer         , intent(in   ) :: lo(:), hi(:), ng, comp, ncomp
      real (kind=dp_t), intent(in   ) :: phi(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind=dp_t), intent(inout) :: phibar(0:,:)
      real (kind=dp_t), intent(in   ) :: dx(:)
      real (kind=dp_t), intent(inout) :: sum(0:)

!     Local variables
      integer                       :: i, j, k, n, index
      integer                       :: nr, nc
      real (kind=dp_t)              :: x,y,z,radius,vol

      nr = size(phibar,dim=1)
      nc = size(phibar,dim=2)

      do k = lo(3),hi(3)
        z = (dble(k)+HALF)*dx(3) - center(3)
        do j = lo(2),hi(2)
          y = (dble(j)+HALF)*dx(2) - center(2)
          do i = lo(1),hi(1)
            x = (dble(i)+HALF)*dx(1) - center(1)

            radius = sqrt(x**2 + y**2 + z**2)
            index = radius / dr 

            if (index .lt. 0 .or. index .gt. nr-1) then
              print *,'RADIUS ',radius
              print *,'BOGUS INDEX IN AVERAGE ',index
              print *,'NOT IN RANGE 0 TO ',nr-1
              print *,'I J K ',i,j,k
              print *,'X Y Z ',x,y,z
              stop
            end if

            vol = FOUR3RD*M_PI * (zl(index+1)  - zl(index)) * &
                                 (zl(index+1)**2 + zl(index+1)*zl(index) + zl(index)**2)

            do n = comp,comp+ncomp-1
              phibar(index,n) = phibar(index,n) + vol * phi(i,j,k,n)
            end do

            sum(index) = sum(index) + vol

          end do
        end do
      end do
 
   end subroutine average_3d_sphr

end module average_module
