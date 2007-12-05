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
  use probin_module, only: grav_const

  implicit none

  private
  public :: make_w0

contains

  subroutine make_w0(nlevs,vel,vel_old,f,Sbar_in,p0,rho0,gam1,eta,dt,dtold,verbose)

    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(  out) :: vel(:,0:)
    real(kind=dp_t), intent(in   ) :: vel_old(:,0:)
    real(kind=dp_t), intent(in   ) :: eta(:,0:,:)
    real(kind=dp_t), intent(inout) :: f(:,0:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:),rho0(:,0:),gam1(:,0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(:,0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold
    integer        , intent(in   ) :: verbose

    integer         :: j,n
    real(kind=dp_t) :: max_vel

    f(:,:) = ZERO

    do n=1,nlevs
       if (spherical .eq. 0) then
          call make_w0_planar(n,vel(n,:),vel_old(n,:),Sbar_in(n,:),p0(n,:),gam1(n,:), &
                              eta(n,:,:),f(n,:),dt,dtold)
       else
          call make_w0_spherical(n,vel(n,:),Sbar_in(n,:),p0(n,:),rho0(n,:),gam1(n,:))
       endif

       max_vel = zero
       do j = 0,nr(n)
          max_vel = max(max_vel, abs(vel(n,j)))
       end do

       if (parallel_IOProcessor() .and. verbose .ge. 1) &
            write(6,*) '... max CFL of w0: ',max_vel * dt / dr(n)
    enddo

  end subroutine make_w0

  subroutine make_w0_planar(n,vel,vel_old,Sbar_in,p0,gam1,eta,f,dt,dtold)

    implicit none
    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(  out) :: vel(0:)
    real(kind=dp_t), intent(in   ) :: vel_old(0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(0:)
    real(kind=dp_t), intent(in   ) :: p0(0:),gam1(0:),eta(0:,:)
    real(kind=dp_t), intent(inout) ::   f(0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold

    ! Local variables
    integer                      :: j
    real(kind=dp_t), allocatable :: vel_old_cen(:)
    real(kind=dp_t), allocatable :: vel_new_cen(:)
    real(kind=dp_t), allocatable ::   force(:)
    real(kind=dp_t), allocatable ::    edge(:)
    real(kind=dp_t)              :: eta_avg

    ! edge-centered
    allocate(edge(0:nr(n)))

    ! cell-centered
    allocate(vel_old_cen(0:nr(n)-1))
    allocate(vel_new_cen(0:nr(n)-1))
    allocate(      force(0:nr(n)-1))

    ! Initialize new velocity to zero.
    vel(0) = ZERO
    do j = 1,nr(n)
       eta_avg = HALF * (eta(j,rho_comp)+eta(j-1,rho_comp))
       vel(j) = vel(j-1) + Sbar_in(j-1) * dr(1) - &
                         ( eta_avg * abs(grav_const) / (gam1(j-1)*p0(j-1)) ) * dr(1)
    end do

    ! Compute the 1/rho0 grad pi0 term.

    do j = 0,nr(n)-1
       vel_old_cen(j) = HALF * (vel_old(j) + vel_old(j+1))
       vel_new_cen(j) = HALF * (vel    (j) + vel    (j+1))
    end do

    force = ZERO
    call mkflux_1d(vel_old_cen,edge,vel_old,force,1,dr(1),dt)

    do j = 0,nr(n)-1
       f(j) = (vel_new_cen(j)-vel_old_cen(j)) / (HALF*(dt+dtold)) + &
            HALF*(vel_old_cen(j)+vel_new_cen(j)) * (edge(j+1)-edge(j)) / dr(1)
    end do

    deallocate(edge)
    deallocate(vel_old_cen,vel_new_cen,force)

  end subroutine make_w0_planar

  subroutine make_w0_spherical(n,vel,Sbar_in,p0,rho0,gam1)
    
    implicit none
    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(  out) :: vel(0:)
    real(kind=dp_t), intent(in   ) :: p0(0:),rho0(0:),gam1(0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(0:)

    ! Local variables
    integer                      :: j
    real(kind=dp_t), allocatable :: c(:),d(:),e(:),u(:),rhs(:)
    real(kind=dp_t), allocatable :: m(:),grav_edge(:),rho0_edge(:)
    
    ! Cell-centered
    allocate(m(0:nr(n)-1))

    ! Edge-centered
    allocate(c(0:nr(n)),d(0:nr(n)),e(0:nr(n)),rhs(0:nr(n)),u(0:nr(n)))
    allocate(grav_edge(0:nr(n)),rho0_edge(0:nr(n)))

    c(:)   = ZERO
    d(:)   = ZERO
    e(:)   = ZERO
    rhs(:) = ZERO
    u(:)   = ZERO
   
    call make_grav_edge(grav_edge,rho0)

    do j = 1,nr(n)
       c(j) = gam1(j-1) * p0(j-1) * zl(j-1)**2 / z(j-1)**2
       c(j) = c(j) / dr(1)**2
    end do

    call cell_to_edge(n,rho0,rho0_edge)

    do j = 1,nr(n)-1

       d(j) = -( gam1(j-1) * p0(j-1) / z(j-1)**2 &
                +gam1(j  ) * p0(j  ) / z(j  )**2 ) * (zl(j)**2/dr(1)**2) &
                - four * rho0_edge(j) * grav_edge(j) / zl(j)
    end do

    do j = 1,nr(n)-1
       rhs(j) = ( gam1(j  )*p0(j  )*Sbar_in(j) - gam1(j-1)*p0(j-1)*Sbar_in(j-1) ) 
       rhs(j) = rhs(j) / dr(1)
    end do

    do j = 0,nr(n)-1
       e(j) = gam1(j) * p0(j) * zl(j+1)**2 / z(j)**2
       e(j) = e(j) / dr(1)**2
    end do

    ! Lower boundary
       d(0) = one
       e(0) = zero
     rhs(0) = zero

    ! Upper boundary
       c(nr(n)) = zero
       d(nr(n)) = one
     rhs(nr(n)) = zero

    ! Call the tridiagonal solver
    call tridiag(c, d, e, rhs, u, nr(n)+1)

    do j = 0,nr(n)
       vel(j) = u(j)
    end do

    deallocate(c,d,e,rhs,u)
    deallocate(m,grav_edge,rho0_edge)

  end subroutine make_w0_spherical

  subroutine tridiag(a,b,c,r,u,n)

    real(kind=dp_t), intent(in   ) :: a(:), b(:), c(:), r(:)
    real(kind=dp_t), intent(  out) :: u(:)
    integer, intent(in)            :: n
    
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
