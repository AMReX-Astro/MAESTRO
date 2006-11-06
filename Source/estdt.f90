module estdt_module

  use bl_types
  use bl_constants_module
  use multifab_module

  implicit none

contains

   subroutine estdt (istep, u, force, w0, dx, cflfac, dtold, dt)

      integer        , intent(in ) :: istep
      type(multifab) , intent(in ) :: u,force
      real(kind=dp_t), intent(in ) :: w0(:)
      real(kind=dp_t), intent(in ) :: dx(:)
      real(kind=dp_t), intent(in ) :: cflfac, dtold
      real(kind=dp_t), intent(out) :: dt

      real(kind=dp_t), pointer:: uop(:,:,:,:)
      real(kind=dp_t), pointer:: fp(:,:,:,:)
      integer :: lo(u%dim),hi(u%dim),ng,dm
      real(kind=dp_t) :: dt_hold
      real(kind=dp_t) :: dtchange
      integer         :: i

      ng = u%ng
      dm = u%dim

      dtchange = 1.1d0
      dt_hold  = 1.d20

      do i = 1, u%nboxes
         if ( multifab_remote(u, i) ) cycle
         uop => dataptr(u, i)
          fp => dataptr(force, i)
         lo =  lwb(get_box(u, i))
         hi =  upb(get_box(u, i))
         select case (dm)
            case (2)
              call estdt_2d(uop(:,:,1,:), fp(:,:,1,:),&
                            w0, lo, hi, ng, dx, dt)
            case (3)
              call estdt_3d(uop(:,:,:,:), fp(:,:,:,:),&
                            w0, lo, hi, ng, dx, dt)
         end select
         dt_hold = min(dt_hold,dt)
      end do

      dt = dt_hold
      dt = dt * cflfac

      if (dtold .gt. 0.0D0 ) dt = min(dt,dtchange*dtold)

      if (istep.le.10)  dt = min(0.005_dp_t,dt)

      print *,'Computing dt at istep ',istep,' to be ',dt

   end subroutine estdt

   subroutine estdt_2d (u, force, w0, lo, hi, ng, dx, dt)

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in ) ::     u(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in ) :: force(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent( in) ::   w0(lo(2):)
      real (kind = dp_t), intent(in ) :: dx(:)
      real (kind = dp_t), intent(out) :: dt

!     Local variables
      real (kind = dp_t)  :: spdx, spdy
      real (kind = dp_t)  :: pforcex, pforcey
      real (kind = dp_t)  :: eps
      integer             :: i,j

      eps = 1.0e-8

      spdx  = 0.0D0 
      spdy  = 0.0D0 
      pforcex = 0.0D0 
      pforcey = 0.0D0 

      do j = lo(2), hi(2)
        spdy = max(spdy ,abs(w0(j)))
      enddo

      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          spdx    = max(spdx ,abs(u(i,j,1)))
          spdy    = max(spdy ,abs(u(i,j,2)))
          pforcex = max(pforcex,abs(force(i,j,1)))
          pforcey = max(pforcey,abs(force(i,j,2)))
        enddo
      enddo

      spdx = spdx / dx(1)
      spdy = spdy / dx(2)

      if (spdx < eps .and. spdy < eps) then

        dt = min(dx(1),dx(2))

      else

        dt = 1.0D0  / max(spdx,spdy)

      endif

      if (pforcex > eps) then
        dt = min(dt,sqrt(2.0D0 *dx(1)/pforcex))
      endif

      if (pforcey > eps) then
        dt = min(dt,sqrt(2.0D0 *dx(2)/pforcey))
      endif

   end subroutine estdt_2d

   subroutine estdt_3d (u, force, w0, lo, hi, ng, dx, dt)

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in ) ::     u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(in ) :: force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent( in) ::   w0(lo(3):)
      real (kind = dp_t), intent(in ) :: dx(:)
      real (kind = dp_t), intent(out) :: dt

!     Local variables
      real (kind = dp_t)  :: spdx, spdy, spdz
      real (kind = dp_t)  :: pforcex, pforcey, pforcez
      real (kind = dp_t)  :: eps
      integer             :: i,j,k

      eps = 1.0e-8

      spdx    = ZERO
      spdy    = ZERO 
      spdz    = ZERO 
      pforcex = ZERO 
      pforcey = ZERO 
      pforcez = ZERO 

      do k = lo(3), hi(3)
        spdz = max(spdz ,abs(w0(k)))
      enddo

      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          spdx    = max(spdx ,abs(u(i,j,k,1)))
          spdy    = max(spdy ,abs(u(i,j,k,2)))
          spdz    = max(spdz ,abs(u(i,j,k,3)))
          pforcex = max(pforcex,abs(force(i,j,k,1)))
          pforcey = max(pforcey,abs(force(i,j,k,2)))
          pforcez = max(pforcez,abs(force(i,j,k,3)))
        enddo
      enddo
      enddo

      spdx = spdx / dx(1)
      spdy = spdy / dx(2)
      spdz = spdz / dx(3)

      if (spdx < eps .and. spdy < eps .and. spdz < eps) then

        dt = min(dx(1),dx(2),dx(3))

      else

        dt = 1.0D0  / max(spdx,spdy,spdz)

      endif

      if (pforcex > eps) then
        dt = min(dt,sqrt(2.0D0*dx(1)/pforcex))
      endif

      if (pforcey > eps) then
        dt = min(dt,sqrt(2.0D0*dx(2)/pforcey))
      endif

      if (pforcez > eps) then
        dt = min(dt,sqrt(2.0D0*dx(3)/pforcez))
      endif

   end subroutine estdt_3d

end module estdt_module
