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
  use cell_to_edge_module

  implicit none

contains

  subroutine make_w0(vel,vel_old,f,Sbar_in,p0,rho0,gam1,dt,dtold,verbose)

    real(kind=dp_t), intent(  out) :: vel(0:)
    real(kind=dp_t), intent(in   ) :: vel_old(0:)
    real(kind=dp_t), intent(inout) ::   f(0:)
    real(kind=dp_t), intent(in   ) :: p0(0:),rho0(0:),gam1(0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold
    integer        , intent(in   ) :: verbose

    integer         :: j,nz
    real(kind=dp_t) :: max_vel

    ! nz is the dimension of a cell-centered quantity
    nz = size(vel,dim=1)-1

    f = ZERO

    if (spherical .eq. 0) then

       call make_w0_planar(vel,vel_old,Sbar_in,f,dt,dtold)

    else

       call make_w0_spherical(vel,Sbar_in,p0,rho0,gam1)

    endif

    max_vel = zero
    do j = 0,nz
       max_vel = max(max_vel, abs(vel(j)))
    end do

    if (parallel_IOProcessor() .and. verbose .ge. 1) &
         write(6,*) '... max CFL of w0: ',max_vel * dt / dr

  end subroutine make_w0

  subroutine make_w0_planar (vel,vel_old,Sbar_in,f,dt,dtold)

    implicit none
    real(kind=dp_t), intent(  out) :: vel(0:)
    real(kind=dp_t), intent(in   ) :: vel_old(0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(0:)
    real(kind=dp_t), intent(inout) ::   f(0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold

!     Local variables
    integer         :: j,nz
    real(kind=dp_t), allocatable :: vel_old_cen(:)
    real(kind=dp_t), allocatable :: vel_new_cen(:)
    real(kind=dp_t), allocatable ::   force(:)
    real(kind=dp_t), allocatable ::    edge(:)

    ! nz is the dimension of a cell-centered quantity
    nz = size(vel,dim=1)-1

    ! edge-centered
    allocate(edge(0:nz))

    ! cell-centered
    allocate(vel_old_cen(0:nz-1),vel_new_cen(0:nz-1),force(0:nz-1))

    ! Initialize new velocity to zero.
    vel(0) = ZERO
    do j = 1,nz
       vel(j) = vel(j-1) + Sbar_in(j-1) * dr
    end do

    ! Compute the 1/rho0 grad pi0 term.

    do j = 0,nz-1
       vel_old_cen(j) = HALF * (vel_old(j) + vel_old(j+1))
       vel_new_cen(j) = HALF * (vel    (j) + vel    (j+1))
    end do

    force = ZERO
    call mkflux_1d(vel_old_cen,edge,vel_old,force,1,dr,dt)

    do j = 0,nz-1
       f(j) = (vel_new_cen(j)-vel_old_cen(j)) / (HALF*(dt+dtold)) + &
            HALF*(vel_old_cen(j)+vel_new_cen(j)) * (edge(j+1)-edge(j)) / dr
    end do

    deallocate(edge)
    deallocate(vel_old_cen,vel_new_cen,force)

  end subroutine make_w0_planar

  subroutine make_w0_spherical (vel,Sbar_in,p0,rho0,gam1)

    implicit none
    real(kind=dp_t), intent(  out) :: vel(0:)
    real(kind=dp_t), intent(in   ) :: p0(0:),rho0(0:),gam1(0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(0:)

!     Local variables
    integer         :: j, nz
    real(kind=dp_t), allocatable :: c(:),d(:),e(:),u(:),rhs(:)
    real(kind=dp_t), allocatable :: m(:),grav_edge(:),rho0_edge(:)
    
    ! nz is the dimension of an cell-centered quantity
    nz = size(vel,dim=1)-1

    ! Cell-centered
    allocate(m(0:nz-1))

    ! Edge-centered
    allocate(c(0:nz),d(0:nz),e(0:nz),rhs(0:nz),u(0:nz))
    allocate(grav_edge(0:nz),rho0_edge(0:nz))

    c(:)   = ZERO
    d(:)   = ZERO
    e(:)   = ZERO
    rhs(:) = ZERO
    u(:)   = ZERO
   
    call make_grav_edge(grav_edge,rho0)

    do j = 1,nz
       c(j) = gam1(j-1) * p0(j-1) * zl(j-1)**2 / z(j-1)**2
       c(j) = c(j) / dr**2
    end do

    call cell_to_edge(rho0,rho0_edge)

    do j = 1,nz-1

       d(j) = -( gam1(j-1) * p0(j-1) / z(j-1)**2 &
                +gam1(j  ) * p0(j  ) / z(j  )**2 ) * (zl(j)**2/dr**2) &
                - four * rho0_edge(j) * grav_edge(j) / zl(j)
    end do

    do j = 1,nz-1
       rhs(j) = ( gam1(j  )*p0(j  )*Sbar_in(j) - gam1(j-1)*p0(j-1)*Sbar_in(j-1) ) 
       rhs(j) = rhs(j) / dr
    end do

    do j = 0,nz-1
       e(j) = gam1(j) * p0(j) * zl(j+1)**2 / z(j)**2
       e(j) = e(j) / dr**2
    end do

    ! Lower boundary
       d(0) = one
       e(0) = zero
     rhs(0) = zero

    ! Upper boundary
       c(nz) = zero
       d(nz) = one
     rhs(nz) = zero

    ! Call the tridiagonal solver
    call tridiag(c, d, e, rhs, u, nz+1)

    do j = 0,nz
       vel(j) = u(j)
    end do

    deallocate(c,d,e,rhs,u)
    deallocate(m,grav_edge,rho0_edge)

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
