module make_w0_module

  ! adjust the base state quantities in response to the heating.
  ! This is step 3 of ABRZ2.

  use bl_types
  use bl_constants_module
  use multifab_module
  use variables
  use geometry
  use mkflux_module
  use make_grav_module

  implicit none

contains

   subroutine make_w0(vel,vel_old,f,Sbar_in,p0,rho0,gam1,dz,dt,verbose)

      real(kind=dp_t), intent(  out) :: vel(:)
      real(kind=dp_t), intent(in   ) :: vel_old(:)
      real(kind=dp_t), intent(inout) ::   f(:)
      real(kind=dp_t), intent(in   ) :: p0(:),rho0(:),gam1(:)
      real(kind=dp_t), intent(in   ) :: Sbar_in(:)
      real(kind=dp_t), intent(in   ) :: dz,dt
      integer        , intent(in   ) :: verbose

      integer         :: j,nz
      real(kind=dp_t) :: max_vel

      nz = size(vel,dim=1)

      print *, '<<< integrating to get w0>>>'

      f = ZERO

      if (spherical .eq. 0) then

        call make_w0_planar(vel,vel_old,Sbar_in,f,dz,dt)

      else

        call make_w0_spherical(vel,Sbar_in,p0,rho0,gam1)

      endif

      max_vel = zero
      do j = 1,nz
         max_vel = max(max_vel, abs(vel(j)))
      end do
      if (parallel_IOProcessor() .and. verbose .ge. 1) &
        write(6,*) '... max CFL of w0: ',max_vel * dt / dr

   end subroutine make_w0

   subroutine make_w0_planar (vel,vel_old,Sbar_in,f,dz,dt)

      implicit none
      real(kind=dp_t), intent(  out) :: vel(:)
      real(kind=dp_t), intent(in   ) :: vel_old(:)
      real(kind=dp_t), intent(in   ) :: Sbar_in(:)
      real(kind=dp_t), intent(inout) ::   f(:)
      real(kind=dp_t), intent(in   ) :: dz,dt

!     Local variables
      integer         :: j,nz
      real(kind=dp_t), allocatable :: vel_old_cen(:)
      real(kind=dp_t), allocatable :: vel_new_cen(:)
      real(kind=dp_t), allocatable ::   force(:)
      real(kind=dp_t), allocatable ::    edge(:)

      nz = size(vel,dim=1)

      allocate(edge(nz))
      allocate(vel_old_cen(nz-1),vel_new_cen(nz-1),force(nz-1))

      ! Initialize new velocity to zero.
      vel(1) = ZERO
      do j = 2,nz
         vel(j) = vel(j-1) + Sbar_in(j-1) * dz
      end do

      ! Compute the 1/rho0 grad pi0 term.

      do j = 1,nz-1
         vel_old_cen(j) = HALF * (vel_old(j) + vel_old(j+1))
         vel_new_cen(j) = HALF * (vel    (j) + vel    (j+1))
      end do

      force = ZERO
      call mkflux_1d(vel_old_cen,edge,vel_old,force,1,dz,dt)

      do j = 1,nz-1
         f(j) = (vel_new_cen(j)-vel_old_cen(j)) / dt + vel_old_cen(j) * (edge(j+1)-edge(j)) / dz
      end do

      deallocate(edge)
      deallocate(vel_old_cen,vel_new_cen,force)

   end subroutine make_w0_planar

   subroutine make_w0_spherical (vel,Sbar_in,p0,rho0,gam1)

      implicit none
      real(kind=dp_t), intent(  out) :: vel(:)
      real(kind=dp_t), intent(in   ) :: p0(:),rho0(:),gam1(:)
      real(kind=dp_t), intent(in   ) :: Sbar_in(:)

!     Local variables
      integer         :: j, k, n, nz
      real(kind=dp_t) :: mencl,rhohalf,integral,velmax
      real(kind=dp_t), allocatable :: c(:),d(:),e(:),u(:),rhs(:)
      real(kind=dp_t), allocatable :: m(:),grav_edge(:)

      nz = size(vel,dim=1)-1

      ! Cell-centered
      allocate(m(nz))

      ! Edge-centered
      allocate(c(nz+1),d(nz+1),e(nz+1),rhs(nz+1),u(nz+1))
      allocate(grav_edge(nz+1))
   
     call make_grav_edge(grav_edge,rho0)

     do j = 2,nz+1
       c(j) = gam1(j-1) * p0(j-1) * zl(j-1)**2 / z(j-1)**2
       c(j) = c(j) / dr**2
     end do

     do j = 2,nz
!       rhohalf = half * (rho0(j) + rho0(j-1))
        if (j == 2) then
           rhohalf = half * (rho0(1) + rho0(2))
        else if (j == nz) then
           rhohalf = half * (rho0(nz-1) + rho0(nz))
        else
           rhohalf = 7.d0/12.d0 * (rho0(j) + rho0(j-1)) - &
                     1.d0/12.d0 * (rho0(j+1) + rho0(j-2))
        endif

       d(j) = -( gam1(j-1) * p0(j-1) / z(j-1)**2 &
                +gam1(j  ) * p0(j  ) / z(j  )**2 ) * (zl(j)**2/dr**2) &
              + four * rhohalf * grav_edge(j) / zl(j)
     end do

     do j = 2,nz
       rhs(j) = ( gam1(j  )*p0(j  )*Sbar_in(j) - gam1(j-1)*p0(j-1)*Sbar_in(j-1) ) 
       rhs(j) = rhs(j) / dr
     end do

     do j = 1,nz
       e(j) = gam1(j) * p0(j) * zl(j+1)**2 / z(j)**2
       e(j) = e(j) / dr**2
     end do

     ! Lower boundary
       d(1) = one
       e(1) = zero
     rhs(1) = zero

     ! Upper boundary
       c(nz+1) = zero
       d(nz+1) = one
     rhs(nz+1) = zero

     ! Call the tridiagonal solver
     call tridiag(c, d, e, rhs, u, nz+1)

     velmax = zero
     do j = 1,nz+1
       vel(j) = u(j)
       velmax = max(velmax,abs(vel(j)))
     end do

     deallocate(m,grav_edge)

   end subroutine make_w0_spherical

   subroutine tridiag(a,b,c,r,u,n)

   integer, intent(in) ::  n

   real(kind=dp_t), intent(in   ) :: a(:), b(:), c(:), r(:)
   real(kind=dp_t), intent(  out) :: u(:)

   integer, parameter :: nmax = 4098
 
   real(kind=dp_t) :: bet, gam(nmax)
   integer         :: j
 
   if (n .gt. nmax ) then
     print *,'tridiag: size exceeded'
     stop
   end if
   if (b(1) .eq. 0) then
    print *,'tridiag: CANT HAVE B(1) = ZERO'
     stop
   end if
 
   bet = b(1)
   u(1) = r(1)/bet
 
   do j = 2,n
     gam(j) = c(j-1)/bet
     bet = b(j) - a(j)*gam(j)
     if (bet .eq. 0) then
       print *,'tridiag: TRIDIAG FAILED'
       stop
     end if
     u(j) = (r(j)-a(j)*u(j-1))/bet
   end do
 
   do j = n-1,1,-1
     u(j) = u(j) - gam(j+1)*u(j+1)
   end do
 
   return
 
   end subroutine tridiag

end module make_w0_module
