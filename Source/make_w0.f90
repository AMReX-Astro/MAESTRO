module make_w0_module

  ! adjust the base state quantities in response to the heating.
  ! This is step 3 of ABRZ2.

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_module
  use heating_module
  use variables
  use make_grav_module

  implicit none

contains

   subroutine make_w0(vel,Sbar_in,p0,rho0,temp0,gam1,dr,dt,spherical)

      integer        , intent(in   ) :: spherical
      real(kind=dp_t), intent(  out) :: vel(:)
      real(kind=dp_t), intent(in   ) :: p0(:),rho0(:),temp0(:),gam1(:)
      real(kind=dp_t), intent(in   ) :: Sbar_in(:)
      real(kind=dp_t), intent(in   ) :: dr,dt

      integer         :: j,nz
      real(kind=dp_t) :: max_vel

      nz = size(vel,dim=1)

      print *, '<<< integrating to get w0>>>'

      if (spherical .eq. 0) then

        call make_w0_planar(vel,Sbar_in,dr)

      else

        call make_w0_spherical(vel,Sbar_in,p0,rho0,temp0,gam1,dr)

      endif

      max_vel = zero
      do j = 1,nz
         max_vel = max(max_vel, abs(vel(j)))
      end do
      print *,'MAX CFL FRAC OF DISPL ',max_vel * dt / dr

   end subroutine make_w0

   subroutine make_w0_planar (vel,Sbar_in,dr)

      implicit none
      real(kind=dp_t), intent(  out) :: vel(:)
      real(kind=dp_t), intent(in   ) :: Sbar_in(:)
      real(kind=dp_t), intent(in   ) :: dr

!     Local variables
      integer         :: j,nz

      nz = size(vel,dim=1)

      ! Initialize velocity to zero.
      vel(1) = ZERO
      do j = 2,nz
         vel(j) = vel(j-1) + Sbar_in(j-1) * dr
      end do

   end subroutine make_w0_planar

   subroutine make_w0_spherical (vel,Sbar_in,p0,rho0,temp0,gam1,dr)

      implicit none
      real(kind=dp_t), intent(  out) :: vel(:)
      real(kind=dp_t), intent(in   ) :: p0(:),rho0(:),temp0(:),gam1(:)
      real(kind=dp_t), intent(in   ) :: Sbar_in(:)
      real(kind=dp_t), intent(in   ) :: dr

!     Local variables
      integer         :: j, k, n, nz
      real(kind=dp_t) :: mencl,rhohalf,integral,velmax
      real(kind=dp_t), allocatable :: z(:),zl(:),c(:),d(:),e(:),u(:),rhs(:)
      real(kind=dp_t), allocatable :: m(:),grav_edge(:),grav_cell(:),beta(:)

      nz = size(vel,dim=1)-1

      ! Cell-centered
      allocate(z(nz))
      allocate(m(nz),c(nz),d(nz),e(nz),rhs(nz))
      allocate(grav_cell(nz))

      ! Edge-centered
      allocate(zl(nz+1))
      allocate(grav_edge(nz+1),beta(nz+1))
   
      ! z(j)  is the location of the cell center of cell (j)
      ! zl(j) is the location of the lower edge of cell (j) 
      z(1) = 0.5_dp_t*dr
      do j = 2,nz
         z(j) = z(j-1) + dr
      enddo
      zl(1) = zero
      do j = 2,nz+1
         zl(j) = zl(j-1) + dr
      enddo

     call make_grav_edge(grav_edge,rho0,dr,1)
     call make_grav_cell(grav_cell,rho0,dr,1)

     ! Define beta_0 at cell edges using the gravity above
     beta(1) = 1.5d0 * rho0(1) - 0.5d0 * rho0(2)
     do j = 2,nz+1
        integral  = rho0(j-1) * grav_cell(j) * dr / (gam1(j-1) * p0(j-1))
        beta(j) = beta(j-1) * exp(-integral)
     end do 

     do j = 1,nz+1
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
!    c(nz) = -beta(nz  ) * zl(nz  )**2
!    d(nz) =  beta(nz+1) * zl(nz+1)**2
     c(nz) = zero
     d(nz) = one
     rhs(nz) = zero

     ! Call the tridiagonal solver
     call tridiag(c, d, e, rhs, u, nz+1)

     velmax = zero
     do j = 1,nz+1
       vel(j) = u(j)
       velmax = max(velmax,abs(vel(j)))
     end do

     deallocate(z,zl,m,grav_cell,grav_edge,beta)

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
