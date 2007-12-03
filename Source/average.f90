module average_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use geometry

  implicit none

  ! if we set phibar to 0, then the evolution of the base state is effectively 
  ! turned off.  Here we set a parameter that allows us to do that.  We will 
  ! use the routines enable_base_evolution() and disable_base_evolution() to 
  ! control this.  By default, we have base state evolution enabled.
  logical, private, save :: evolve_base = .true.

contains

  subroutine enable_base_evolution()
    
    implicit none

    evolve_base = .true.

  end subroutine enable_base_evolution
  

  subroutine disable_base_evolution()

    implicit none

    evolve_base = .false.

  end subroutine disable_base_evolution



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
      real(kind=dp_t), allocatable :: source_buffer(:), target_buffer(:)

      dm = phi(1)%dim
      ng = phi(1)%ng
      nr = size(phibar,dim=1)
      nlevs = size(dx,dim=1)

      phibar(:,:) = ZERO

      if (evolve_base) then

         if (spherical .eq. 1) then
            allocate(vol_grid(0:nr-1))
         end if

         allocate(vol_proc(0:nr-1),vol_tot(0:nr-1))
         allocate(phibar_proc(0:nr-1,ncomp))
         allocate(source_buffer(ncomp), target_buffer(ncomp))

         vol_proc(:) = ZERO
         vol_tot(:)  = ZERO

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
                  vol_proc(lo(2):hi(2)) = vol_proc(lo(2):hi(2)) + &
                       (hi(1)-lo(1)+1)*dx(n,1)*dx(n,2)

               case (3)
                  if (spherical .eq. 1) then
                     vol_grid(:) = ZERO
                     call average_3d_sphr(pp(:,:,:,:),phibar_proc,lo,hi,ng,dx(n,:),vol_grid,comp,ncomp)
                     vol_proc = vol_proc + vol_grid
                  else
                     call average_3d(pp(:,:,:,:),phibar_proc,lo,hi,ng,comp,ncomp,dx(n,:))
                     vol_proc(lo(3):hi(3)) = vol_proc(lo(3):hi(3)) + &
                          (hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)*dx(n,1)*dx(n,2)*dx(n,3)
                  end if

               end select
            end do
         end do

         if (dm .eq. 2 .or. (dm.eq.3.and.spherical.eq.0)) then
            do k = 0,nr-1
               call parallel_reduce(vol_tot(k),  vol_proc(k),      MPI_SUM)

               ! put all the components for the current k into a buffer array
               ! and reduce
               source_buffer(:) = phibar_proc(k,:)
               call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
               phibar(k,:) = target_buffer(:)

               phibar(k,:) = phibar(k,:) / vol_tot(k)
            end do

         else
            do k = 0,nr-1
               call parallel_reduce(vol_tot(k),  vol_proc(k),      MPI_SUM)

               ! put all the components for the current k into a buffer array
               ! and reduce
               source_buffer(:) = phibar_proc(k,:)
               call parallel_reduce(target_buffer, source_buffer, MPI_SUM)
               phibar(k,:) = target_buffer(:)

               if (vol_tot(k) .gt. ZERO) then
                  phibar(k,:) = phibar(k,:) / vol_tot(k)
               else
                  phibar(k,:) = ZERO
               end if

            end do
            deallocate(vol_grid)

         end if
         deallocate(vol_proc,vol_tot,phibar_proc)
         deallocate(source_buffer, target_buffer)

      endif

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
      integer                       :: ii, jj, kk
      integer                       :: nr, nc
      real (kind=dp_t)              :: x,y,z,radius,vol
      real (kind=dp_t)              :: xx, yy, zz
      real (kind=dp_t)              :: xmin, ymin, zmin
      integer :: nsub

      ! compute nsub such that we are always guaranteed to fill each of
      ! the base state radial bins
      nsub = int(dx(1)/dr(1)) + 1
      nr = size(phibar,dim=1)

      do k = lo(3),hi(3)
        zmin = dble(k)*dx(3) - center(3)

        do j = lo(2),hi(2)
          ymin = dble(j)*dx(2) - center(2)

          do i = lo(1),hi(1)
            xmin = dble(i)*dx(1) - center(1)


            do kk = 0, nsub-1
               zz = zmin + (dble(kk) + HALF)*dx(3)/nsub

               do jj = 0, nsub-1
                  yy = ymin + (dble(jj) + HALF)*dx(2)/nsub

                  do ii = 0, nsub-1
                     xx = xmin + (dble(ii) + HALF)*dx(1)/nsub

                     radius = sqrt(xx**2 + yy**2 + zz**2)
                     index = radius / dr(1) 

                     if (index .lt. 0 .or. index .gt. nr-1) then
                        print *,'RADIUS ',radius
                        print *,'BOGUS INDEX IN AVERAGE ',index
                        print *,'NOT IN RANGE 0 TO ',nr-1
                        print *,'I J K ',i,j,k
                        print *,'X Y Z ',x,y,z
                        stop
                     end if

                     vol = FOUR3RD*M_PI * dr(1) * &
                          (zl(index+1)**2 + zl(index+1)*zl(index) + zl(index)**2)

                     do n = comp,comp+ncomp-1
                        phibar(index,n) = phibar(index,n) + vol * phi(i,j,k,n)
                     end do

                     sum(index) = sum(index) + vol

                  enddo
               enddo
            enddo

          end do
        end do
      end do

   end subroutine average_3d_sphr

end module average_module
