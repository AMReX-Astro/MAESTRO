! compute w0 -- the base state velocity.  This is based on the average 
! heating in a layer (Sbar) and the mixing (the eta quantities).  The
! computation of w0 for plane-parallel atmospheres was first described
! in paper II, with modifications due to mixing in paper III.  For 
! spherical geometry, it was first described in paper III.

module make_w0_module

  use bl_types

  implicit none

  private

  public :: make_w0

contains

  subroutine make_w0(nlevs,vel,vel_old,f,Sbar_in,rho0,p0_old,p0_new, &
                     gamma1bar_old,gamma1bar_new,psi,dt,dtold)

    use parallel
    use bl_prof_module
    use geometry, only: spherical, nr, dr
    use bl_constants_module
    use probin_module, only: verbose

    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(  out) :: vel(:,0:)
    real(kind=dp_t), intent(in   ) :: vel_old(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    real(kind=dp_t), intent(inout) :: f(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:), p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar_old(:,0:), gamma1bar_new(:,0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(:,0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold

    integer         :: r,n
    real(kind=dp_t) :: max_vel

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_w0")

    f = ZERO

    do n=1,nlevs
       if (spherical .eq. 0) then
          call make_w0_planar(n,vel(n,0:),vel_old(n,0:),Sbar_in(n,0:), &
                              p0_old(n,0:),p0_new(n,0:), &
                              gamma1bar_old(n,0:),gamma1bar_new(n,0:), &
                              psi(n,0:),f(n,0:),dt,dtold)
       else
          call make_w0_spherical(n,vel(n,:),vel_old(n,0:),Sbar_in(n,:), &
                                 rho0(n,:),p0_old(n,:), &
                                 gamma1bar_old(n,:), &
                                 f(n,0:),dt,dtold)
       endif

       max_vel = zero
       do r = 0,nr(n)
          max_vel = max(max_vel, abs(vel(n,r)))
       end do

       if (parallel_IOProcessor() .and. verbose .ge. 1) &
            write(6,*) '... max CFL of w0: ',max_vel * dt / dr(n)
    enddo

    call destroy(bpt)

  end subroutine make_w0

  subroutine make_w0_planar(n,vel,vel_old,Sbar_in,p0_old,p0_new, &
                            gamma1bar_old,gamma1bar_new,psi,f,dt,dtold)

    use geometry, only: nr, dr
    use variables, only: rho_comp
    use bl_constants_module
    use probin_module, only: grav_const

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(  out) :: vel(0:)
    real(kind=dp_t), intent(in   ) :: vel_old(0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(0:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:), p0_new(0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar_old(0:), gamma1bar_new(0:)
    real(kind=dp_t), intent(in   ) :: psi(0:)
    real(kind=dp_t), intent(inout) ::   f(0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold

    ! Local variables
    integer                      :: r
    real(kind=dp_t), allocatable :: vel_old_cen(:)
    real(kind=dp_t), allocatable :: vel_new_cen(:)
    real(kind=dp_t), allocatable ::   force(:)
    real(kind=dp_t)              :: vel_avg, div_avg, dt_avg, gamma1bar_p0_avg

    ! Cell-centered
    allocate(vel_old_cen(0:nr(n)-1))
    allocate(vel_new_cen(0:nr(n)-1))
    allocate(      force(0:nr(n)-1))

    ! Initialize new velocity to zero.
    vel(0) = ZERO
    
    do r = 1,nr(n)
       gamma1bar_p0_avg = (gamma1bar_old(r-1)+gamma1bar_new(r-1))*(p0_old(r-1)+p0_new(r-1)) &
            / 4.0d0

       vel(r) = vel(r-1) + Sbar_in(r-1) * dr(n) &
          - ( psi(r-1) / gamma1bar_p0_avg ) * dr(n)
    end do

    ! Compute the forcing term in the base state velocity equation, - 1/rho0 grad pi0 
    dt_avg = HALF * (dt + dtold)
    do r = 0,nr(n)-1
       vel_old_cen(r) = HALF * (vel_old(r) + vel_old(r+1))
       vel_new_cen(r) = HALF * (vel    (r) + vel    (r+1))
       vel_avg = HALF * (dt *  vel_old_cen(r)           + dtold *  vel_new_cen(r)  ) / dt_avg
       div_avg = HALF * (dt * (vel_old(r+1)-vel_old(r)) + dtold * (vel(r+1)-vel(r))) / dt_avg
       f(r) = (vel_new_cen(r)-vel_old_cen(r)) / dt_avg + &
               vel_avg * div_avg / dr(n)
    end do

    deallocate(vel_old_cen,vel_new_cen,force)

  end subroutine make_w0_planar

  subroutine make_w0_spherical(n,vel,vel_old,Sbar_in,rho0,p0,gamma1bar,f,dt,dtold)

    use geometry, only: base_cc_loc, nr, base_loedge_loc, dr
    use make_grav_module
    use cell_to_edge_module
    use bl_constants_module
    
    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(  out) :: vel(0:)
    real(kind=dp_t), intent(in   ) :: vel_old(0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(0:)
    real(kind=dp_t), intent(in   ) :: rho0(0:),p0(0:),gamma1bar(0:)
    real(kind=dp_t), intent(inout) ::   f(0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold

    ! Local variables
    integer                      :: r
    real(kind=dp_t), allocatable :: vel_old_cen(:)
    real(kind=dp_t), allocatable :: vel_new_cen(:)
    real(kind=dp_t), allocatable :: c(:),d(:),e(:),u(:),rhs(:)
    real(kind=dp_t), allocatable :: m(:),grav_edge(:),rho0_edge(:),rho0star(:)
    real(kind=dp_t)              :: vel_avg, div_avg, dt_avg

    real(kind=dp_t), parameter :: eps = 1.d-8

    ! Cell-centered
    allocate(m(0:nr(n)-1))
    allocate(vel_old_cen(0:nr(n)-1))
    allocate(vel_new_cen(0:nr(n)-1))
    allocate(rho0star(0:nr(n)-1))

    ! Edge-centered
    allocate(c(0:nr(n)),d(0:nr(n)),e(0:nr(n)),rhs(0:nr(n)),u(0:nr(n)))
    allocate(grav_edge(0:nr(n)),rho0_edge(0:nr(n)))

    c   = ZERO
    d   = ZERO
    e   = ZERO
    rhs = ZERO
    u   = ZERO
   
    ! rho0star is defined as (dp_0/dr)/g -- i.e. it should be rho0 except
    ! outside the star where we hold density constant (and therefore are not
    ! in HSE).  In that region, rho0star = 0
    rho0star(0) = rho0(0)
    rho0star(nr(n)-1) = ZERO
    
    do r = 1, nr(n)-2
       if ( (p0(r-1)-p0(r+1))/p0(r) < eps) then
          rho0star(r) = ZERO
       else
          rho0star(r) = rho0(r)
       endif
    enddo

    call make_grav_edge(n,grav_edge,rho0)

    ! Note that we are solving for (r^2 w0), not just w0. 

    do r = 1,nr(n)
!      c(r) = gamma1bar(r-1) * p0(r-1) * base_loedge_loc(n,r-1)**2 / base_cc_loc(n,r-1)**2
       c(r) = gamma1bar(r-1) * p0(r-1)                             / base_cc_loc(n,r-1)**2
       c(r) = c(r) / dr(n)**2
    end do

    call cell_to_edge(n,rho0star,rho0_edge)

    do r = 1,nr(n)-1

!      d(r) = -( gamma1bar(r-1) * p0(r-1) / base_cc_loc(n,r-1)**2 &
!               +gamma1bar(r  ) * p0(r  ) / base_cc_loc(n,r  )**2 ) &
!               * (base_loedge_loc(n,r)**2/dr(n)**2) &
!               - four * rho0_edge(r) * grav_edge(r) / base_loedge_loc(n,r)

       d(r) = -( gamma1bar(r-1) * p0(r-1) / base_cc_loc(n,r-1)**2 &
                +gamma1bar(r  ) * p0(r  ) / base_cc_loc(n,r  )**2 ) / dr(n)**2 &
                - four * rho0_edge(r) * grav_edge(r) / (base_loedge_loc(n,r))**3
    end do

    do r = 0,nr(n)-1
!      e(r) = gamma1bar(r) * p0(r) * base_loedge_loc(n,r+1)**2 / base_cc_loc(n,r)**2
       e(r) = gamma1bar(r) * p0(r)                             / base_cc_loc(n,r)**2
       e(r) = e(r) / dr(n)**2
    end do

    do r = 1,nr(n)-1
       rhs(r) = ( gamma1bar(r  )*p0(r  )*Sbar_in(r) - gamma1bar(r-1)*p0(r-1)*Sbar_in(r-1) ) 
       rhs(r) = rhs(r) / dr(n)
    end do

    ! Lower boundary
       d(0) = one
       e(0) = zero
     rhs(0) = zero

    ! Upper boundary
       c(nr(n)) = -one
       d(nr(n)) =  one
     rhs(nr(n)) = zero

    ! Call the tridiagonal solver
    call tridiag(c, d, e, rhs, u, nr(n)+1)

    vel(0) = ZERO
    do r = 1,nr(n)
       vel(r) = u(r) / base_loedge_loc(n,r)**2
    end do

    ! Compute the forcing term in the base state velocity equation, - 1/rho0 grad pi0 
    dt_avg = HALF * (dt + dtold)
    do r = 0,nr(n)-1
       vel_old_cen(r) = HALF * (vel_old(r) + vel_old(r+1))
       vel_new_cen(r) = HALF * (vel    (r) + vel    (r+1))
       vel_avg = HALF * (dt *  vel_old_cen(r)           + dtold *  vel_new_cen(r)  ) / dt_avg
       div_avg = HALF * (dt * (vel_old(r+1)-vel_old(r)) + dtold * (vel(r+1)-vel(r))) / dt_avg
       f(r) = (vel_new_cen(r)-vel_old_cen(r)) / dt_avg + &
               vel_avg * div_avg / dr(n)
    end do

    deallocate(c,d,e,rhs,u)
    deallocate(m,grav_edge,rho0_edge,rho0star)
    deallocate(vel_old_cen,vel_new_cen)

  end subroutine make_w0_spherical

  subroutine tridiag(a,b,c,r,u,n)

    use bl_error_module

    real(kind=dp_t), intent(in   ) :: a(:), b(:), c(:), r(:)
    real(kind=dp_t), intent(  out) :: u(:)
    integer, intent(in)            :: n

    real(kind=dp_t)              :: bet
    real(kind=dp_t), allocatable :: gam(:)
    integer                      :: j

    allocate(gam(n))
    
    if ( b(1) .eq. 0 ) call bl_error('tridiag: CANT HAVE B(1) = ZERO')
    
    bet = b(1)
    u(1) = r(1)/bet
    
    do j = 2,n
       gam(j) = c(j-1)/bet
       bet = b(j) - a(j)*gam(j)
       if ( bet .eq. 0 ) call bl_error('tridiag: TRIDIAG FAILED')
       u(j) = (r(j)-a(j)*u(j-1))/bet
    end do
    
    do j = n-1,1,-1
       u(j) = u(j) - gam(j+1)*u(j+1)
    end do
    
  end subroutine tridiag
  
end module make_w0_module
