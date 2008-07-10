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

  subroutine make_w0(nlevs,vel,vel_old,f,Sbar_in,rho0_old,rho0_new,p0_old,p0_new, &
                     gamma1bar_old,gamma1bar_new,delta_p0_ptherm_bar,psi,etarho,etarho_cc,dt,dtold)

    use parallel
    use bl_prof_module
    use geometry, only: spherical, dr, r_start_coord, r_end_coord
    use bl_constants_module
    use probin_module, only: verbose

    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(  out) :: vel(:,0:)
    real(kind=dp_t), intent(in   ) :: vel_old(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    real(kind=dp_t), intent(in   ) :: etarho(:,0:)
    real(kind=dp_t), intent(in   ) :: etarho_cc(:,0:)
    real(kind=dp_t), intent(inout) :: f(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:), rho0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:), p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar_old(:,0:), gamma1bar_new(:,0:)
    real(kind=dp_t), intent(in   ) :: delta_p0_ptherm_bar(:,0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(:,0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold

    integer         :: r,n
    real(kind=dp_t) :: max_vel

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_w0")

    f = ZERO

    if (spherical .eq. 0) then

       call make_w0_planar(nlevs,vel,vel_old,Sbar_in,p0_old,p0_new,gamma1bar_old, &
                           gamma1bar_new,delta_p0_ptherm_bar,psi,f,dt,dtold)

    else

       do n=1,nlevs
          call make_w0_spherical(n,vel(n,:),vel_old(n,0:),Sbar_in(n,0:), &
                                 rho0_old(n,:),rho0_new(n,:),p0_old(n,0:),p0_new(n,0:), &
                                 gamma1bar_old(n,0:),gamma1bar_new(n,0:), &
                                 delta_p0_ptherm_bar(n,0:), &
                                 etarho(n,0:),etarho_cc(n,0:), &
                                 f(n,0:),dt,dtold)
       end do

    end if

    do n=1,nlevs
       max_vel = zero
       do r=r_start_coord(n),r_end_coord(n)+1
          max_vel = max(max_vel, abs(vel(n,r)))
       end do
       
       if (parallel_IOProcessor() .and. verbose .ge. 1) &
            write(6,*) '... max CFL of w0: ',max_vel * dt / dr(n)
    end do

    call destroy(bpt)

  end subroutine make_w0

  subroutine make_w0_planar(nlevs,vel,vel_old,Sbar_in,p0_old,p0_new, &
                            gamma1bar_old,gamma1bar_new,delta_p0_ptherm_bar, &
                            psi,f,dt,dtold)

    use geometry, only: nr_fine, r_start_coord, r_end_coord, dr, base_cutoff_density_coord
    use variables, only: rho_comp
    use bl_constants_module
    use probin_module, only: grav_const, dpdt_factor, base_cutoff_density

    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(  out) :: vel(:,0:)
    real(kind=dp_t), intent(in   ) :: vel_old(:,0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:), p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar_old(:,0:), gamma1bar_new(:,0:)
    real(kind=dp_t), intent(in   ) :: delta_p0_ptherm_bar(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    real(kind=dp_t), intent(inout) ::   f(:,0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold

    ! Local variables
    integer                      :: r, n, i, refrat
    real(kind=dp_t), allocatable :: vel_old_cen(:,:)
    real(kind=dp_t), allocatable :: vel_new_cen(:,:)
    real(kind=dp_t)              :: vel_avg, div_avg, dt_avg, gamma1bar_p0_avg
    real(kind=dp_t)              :: volume_discrepancy, offset

    ! Cell-centered
    allocate(vel_old_cen(nlevs,0:nr_fine-1))
    allocate(vel_new_cen(nlevs,0:nr_fine-1))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Multilevel Outline
    !
    ! Compute vel at level 1 only
    ! do n=2,nlevs
    !   Compute vel on edges at level n
    !   Obtain the starting value of vel from the coarser grid
    !   do i=n-1,1,-1
    !     Restrict vel from level n to level i
    !     Compare the difference between vel at top of level n to the corresponding point
    !       on level i
    !     Offset the vel on level i above this point
    !   end do
    ! end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Compute vel on edges at level n
    do n=1,nlevs

       if (n .eq. 1) then
          ! Initialize new velocity at bottom of coarse base array to zero.
          vel(1,0) = ZERO
       else
          ! Obtain the starting value of vel from the coarser grid
          vel(n,r_start_coord(n)) = vel(n-1,r_start_coord(n)/2)
       end if

       do r=r_start_coord(n)+1,r_end_coord(n)+1
          gamma1bar_p0_avg = (gamma1bar_old(n,r-1)+gamma1bar_new(n,r-1)) * &
               (p0_old(n,r-1)+p0_new(n,r-1))/4.0d0
       
          if (r-1 .lt. base_cutoff_density_coord(n)) then
             volume_discrepancy = dpdt_factor * delta_p0_ptherm_bar(n,r-1)/dt
          else
             volume_discrepancy = 0.0d0
          end if
          
          vel(n,r) = vel(n,r-1) + Sbar_in(n,r-1) * dr(n) &
               - ( (psi(n,r-1)+volume_discrepancy) / gamma1bar_p0_avg ) * dr(n)
       end do

       do i=n-1,1,-1

          refrat = 2**(n-i)

          ! Compare the difference between vel at top of level n to the corresponding point
          !   on level i
          offset = vel(n,r_end_coord(n)+1) - vel(i,(r_end_coord(n)+1)/refrat)

          ! Restrict vel from level n to level i
          do r=r_start_coord(n),r_end_coord(n)+1,refrat
             vel(i,r/refrat) = vel(n,r)
          end do

          ! Offset the vel on level i above this point
          do r=(r_end_coord(n)+1)/refrat+1,r_end_coord(i)+1
             vel(i,r) = vel(i,r) + offset
          end do

       end do

    end do

    do n=1,nlevs
       
       ! Compute the forcing term in the base state velocity equation, - 1/rho0 grad pi0 
       dt_avg = HALF * (dt + dtold)
       do r=r_start_coord(n),r_end_coord(n)
          vel_old_cen(n,r) = HALF * (vel_old(n,r) + vel_old(n,r+1))
          vel_new_cen(n,r) = HALF * (vel(n,r) + vel(n,r+1))
          vel_avg = HALF * (dt * vel_old_cen(n,r) + dtold *  vel_new_cen(n,r)) / dt_avg
          div_avg = HALF * (dt * (vel_old(n,r+1)-vel_old(n,r)) + &
               dtold * (vel(n,r+1)-vel(n,r))) / dt_avg
          f(n,r) = (vel_new_cen(n,r)-vel_old_cen(n,r)) / dt_avg + vel_avg * div_avg / dr(n)
       end do

    end do

    deallocate(vel_old_cen,vel_new_cen)

  end subroutine make_w0_planar

  subroutine make_w0_spherical(n,vel,vel_old,Sbar_in,rho0,rho0_new,p0,p0_new, &
                               gamma1bar,gamma1bar_new,delta_p0_ptherm_bar, &
                               etarho,etarho_cc,f,dt,dtold)

    use geometry, only: r_cc_loc, nr_fine, r_edge_loc, dr, r_end_coord
    use make_grav_module
    use cell_to_edge_module
    use bl_constants_module
    use fundamental_constants_module, only: Gconst
    use probin_module, only: dpdt_factor, base_cutoff_density

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(  out) :: vel(0:)
    real(kind=dp_t), intent(in   ) :: vel_old(0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(0:)
    real(kind=dp_t), intent(in   ) :: rho0(0:),rho0_new(0:)
    real(kind=dp_t), intent(in   ) :: p0(0:),p0_new(0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(0:),gamma1bar_new(0:),delta_p0_ptherm_bar(0:)
    real(kind=dp_t), intent(in   ) :: etarho(0:),etarho_cc(0:)
    real(kind=dp_t), intent(inout) ::   f(0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold

    ! Local variables
    integer                      :: r
    real(kind=dp_t), allocatable :: vel_old_cen(:)
    real(kind=dp_t), allocatable :: vel_new_cen(:)
    real(kind=dp_t), allocatable :: c(:),d(:),e(:),u(:),rhs(:)
    real(kind=dp_t), allocatable :: m(:)
    real(kind=dp_t)              :: vel_avg, div_avg, dt_avg
    real(kind=dp_t), allocatable :: vel_bar(:)
    real(kind=dp_t), allocatable :: grav_edge(:)

    real(kind=dp_t), parameter :: eps = 1.d-8

    real(kind=dp_t) :: dpdr

    real(kind=dp_t) :: volume_discrepancy

    real(kind=dp_t), allocatable :: gamma1bar_nph(:), rho0_nph(:), p0_nph(:)

    ! Cell-centered
    allocate(m          (0:nr_fine-1))
    allocate(vel_old_cen(0:nr_fine-1))
    allocate(vel_new_cen(0:nr_fine-1))

    allocate(gamma1bar_nph(0:nr_fine-1))
    allocate(     rho0_nph(0:nr_fine-1))
    allocate(       p0_nph(0:nr_fine-1))


    ! Edge-centered
    allocate(c(0:nr_fine),d(0:nr_fine),e(0:nr_fine),rhs(0:nr_fine),u(0:nr_fine))
    allocate(vel_bar(0:nr_fine))
    allocate(grav_edge(0:nr_fine))


    ! create time-centered base-state quantities
    do r = 0, r_end_coord(n)
       p0_nph(r)        = HALF*(p0(r)        + p0_new(r))
       rho0_nph(r)      = HALF*(rho0(r)      + rho0_new(r))
       gamma1bar_nph(r) = HALF*(gamma1bar(r) + gamma1bar_new(r))       
    enddo


    ! NOTE:  we first solve for the w0 resulting only from Sbar -- then we will
    ! solve for the update to w0.  We integrate d/dr (r^2 w0) = (r^2 Sbar)

    vel_bar = ZERO
    do r=1,r_end_coord(n)+1

       if (rho0(r-1) .gt. base_cutoff_density) then
          volume_discrepancy = dpdt_factor * delta_p0_ptherm_bar(r-1)/dt
       else
          volume_discrepancy = ZERO
       endif

       vel_bar(r) = vel_bar(r-1) + dr(n) * Sbar_in(r-1) * r_cc_loc(n,r-1)**2 - &
            dr(n)* volume_discrepancy * r_cc_loc(n,r-1)**2 / &
            (gamma1bar_nph(r-1)*p0_nph(r-1))

    end do

    do r = 1,r_end_coord(n)+1
       vel_bar(r) = vel_bar(r) / r_edge_loc(n,r)**2
    end do


    ! make the edge-centered gravity
    call make_grav_edge(n,grav_edge,rho0_nph)

    ! NOTE:  now we solve for the remainder of (r^2 * w0)

    c   = ZERO
    d   = ZERO
    e   = ZERO
    rhs = ZERO
    u   = ZERO
   
    ! Note that we are solving for (r^2 w0), not just w0. 

    do r=1,r_end_coord(n)+1
       c(r) = gamma1bar_nph(r-1) * p0_nph(r-1) / r_cc_loc(n,r-1)**2
       c(r) = c(r) / dr(n)**2
    end do

    do r=1,r_end_coord(n)

       d(r) = -( gamma1bar_nph(r-1) * p0_nph(r-1) / r_cc_loc(n,r-1)**2 &
                +gamma1bar_nph(r  ) * p0_nph(r  ) / r_cc_loc(n,r  )**2 ) / dr(n)**2 

       dpdr = (p0_nph(r)-p0_nph(r-1))/dr(n)
       d(r) = d(r) - four * dpdr / (r_edge_loc(n,r))**3
    end do

    do r = 0,r_end_coord(n)
       e(r) = gamma1bar_nph(r) * p0_nph(r) / r_cc_loc(n,r)**2
       e(r) = e(r) / dr(n)**2
    end do

    do r = 1,r_end_coord(n)
       dpdr = (p0_nph(r)-p0_nph(r-1))/dr(n)
       rhs(r) = four * dpdr * vel_bar(r) / r_edge_loc(n,r) - &
            grav_edge(r) * (r_cc_loc(n,r  )**2 * etarho_cc(r  ) - &
                            r_cc_loc(n,r-1)**2 * etarho_cc(r-1)) / &
                           (dr(n) * r_edge_loc(n,r)**2) - &
            four * M_PI * Gconst * HALF * (rho0_nph(r) + rho0_nph(r-1)) * etarho(r)
    end do

    ! Lower boundary
       d(0) = one
       e(0) = zero
     rhs(0) = zero

    ! Upper boundary
!      c(r_end_coord(n)+1) = -one
       c(r_end_coord(n)+1) = zero
       d(r_end_coord(n)+1) =  one
     rhs(r_end_coord(n)+1) = zero

    ! Call the tridiagonal solver
    call tridiag(c, d, e, rhs, u, r_end_coord(n)+2)

    vel(0) = ZERO
    do r=1,r_end_coord(n)+1
       vel(r) = u(r) / r_edge_loc(n,r)**2
    end do

    do r=0,r_end_coord(n)+1
       vel(r) = vel(r) + vel_bar(r)
    end do

    ! Compute the forcing term in the base state velocity equation, - 1/rho0 grad pi0 
    dt_avg = HALF * (dt + dtold)
    do r = 0,r_end_coord(n)
       vel_old_cen(r) = HALF * (vel_old(r) + vel_old(r+1))
       vel_new_cen(r) = HALF * (vel    (r) + vel    (r+1))
       vel_avg = HALF * (dt *  vel_old_cen(r)           + dtold *  vel_new_cen(r)  ) / dt_avg
       div_avg = HALF * (dt * (vel_old(r+1)-vel_old(r)) + dtold * (vel(r+1)-vel(r))) / dt_avg
       f(r) = (vel_new_cen(r)-vel_old_cen(r)) / dt_avg + &
               vel_avg * div_avg / dr(n)
    end do

    deallocate(c,d,e,rhs,u)
    deallocate(m)
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
