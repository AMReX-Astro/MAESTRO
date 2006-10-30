module estdt_module

  use bl_types
  use multifab_module

  implicit none

contains

   subroutine estdt (istep, u, s, gp, force, rho0, dx, cflfac, dtold, dt)

      integer        , intent(in ) :: istep
      type(multifab) , intent(in ) :: u,s,gp,force
      real(kind=dp_t), intent(in ) :: rho0(:)
      real(kind=dp_t), intent(in ) :: dx(:)
      real(kind=dp_t), intent(in ) :: cflfac, dtold
      real(kind=dp_t), intent(out) :: dt

      real(kind=dp_t), pointer:: uop(:,:,:,:), sop(:,:,:,:)
      real(kind=dp_t), pointer:: gpp(:,:,:,:),  fp(:,:,:,:)
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
         sop => dataptr(s, i)
         gpp => dataptr(gp, i)
          fp => dataptr(force, i)
         lo =  lwb(get_box(u, i))
         hi =  upb(get_box(u, i))
         select case (dm)
            case (2)
              call estdt_2d(uop(:,:,1,:), sop(:,:,1,1), gpp(:,:,1,:), fp(:,:,1,:),&
                            rho0, lo, hi, ng, dx, dt)
            case (3)
              call estdt_3d(uop(:,:,:,:), sop(:,:,:,1), gpp(:,:,:,:), fp(:,:,:,:),&
                            rho0, lo, hi, ng, dx, dt)
         end select
         dt_hold = min(dt_hold,dt)
      end do

      dt = dt_hold
      dt = dt * cflfac

      if (dtold .gt. 0.0D0 ) dt = min(dt,dtchange*dtold)

      print *,'Computing dt at istep ',istep,' to be ',dt

   end subroutine estdt

   subroutine estdt_2d (u, s, gp, force, rho0, lo, hi, ng, dx, dt)

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in ) ::     u(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in ) ::     s(lo(1)-ng:,lo(2)-ng:)  
      real (kind = dp_t), intent(in ) ::    gp(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent(in ) :: force(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent( in) :: rho0(lo(2):)
      real (kind = dp_t), intent(in ) :: dx(:)
      real (kind = dp_t), intent(out) :: dt

!     Local variables
      real (kind = dp_t)  spdx, spdy
      real (kind = dp_t)  pforcex, pforcey
      real (kind = dp_t)  eps
      real (kind = dp_t)  vert_force
      integer :: i, j

      eps = 1.0e-8

      spdx  = 0.0D0 
      spdy  = 0.0D0 
      pforcex = 0.0D0 
      pforcey = 0.0D0 

      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          spdx    = max(spdx ,abs(u(i,j,1))/dx(1))
          spdy    = max(spdy ,abs(u(i,j,2))/dx(2))
          pforcex = max(pforcex,abs(gp(i,j,1)/s(i,j)-force(i,j,1)))
          vert_force = (s(i,j)-rho0(j))/s(i,j) * force(i,j,2)
          pforcey = max(pforcey,abs(gp(i,j,2)/s(i,j)-vert_force))
        enddo
      enddo

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

   subroutine estdt_3d (u, s, gp, force, rho0, lo, hi, ng, dx, dt)

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in ) ::     u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(in ) ::     s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)  
      real (kind = dp_t), intent(in ) ::    gp(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent(in ) :: force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent( in) :: rho0(lo(3):)
      real (kind = dp_t), intent(in ) :: dx(:)
      real (kind = dp_t), intent(out) :: dt

!     Local variables
      real (kind = dp_t)  spdx,spdy, spdz
      real (kind = dp_t)  pforcex,pforcey,pforcez
      real (kind = dp_t)  eps
      real (kind = dp_t)  vert_force
      integer :: i, j, k

      eps = 1.0e-8

      spdx  = 0.0D0 
      spdy  = 0.0D0 
      spdz  = 0.0D0 
      pforcex = 0.0D0 
      pforcey = 0.0D0 
      pforcez = 0.0D0 

      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          spdx    = max(spdx ,abs(u(i,j,k,1))/dx(1))
          spdy    = max(spdy ,abs(u(i,j,k,2))/dx(2))
          spdz    = max(spdz ,abs(u(i,j,k,3))/dx(3))
          pforcex = max(pforcex,abs(gp(i,j,k,1)/s(i,j,k)-force(i,j,k,1)))
          pforcey = max(pforcey,abs(gp(i,j,k,2)/s(i,j,k)-force(i,j,k,2)))
          vert_force = (s(i,j,k)-rho0(k))/s(i,j,k) * force(i,j,k,3)
          pforcez = max(pforcez,abs(gp(i,j,k,3)/s(i,j,k)-vert_force))
        enddo
      enddo
      enddo

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
