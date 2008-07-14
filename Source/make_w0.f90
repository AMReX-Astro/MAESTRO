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

  subroutine make_w0(nlevs,w0,w0_old,w0_force,Sbar_in,rho0_old,rho0_new,p0_old,p0_new, &
                     gamma1bar_old,gamma1bar_new,delta_p0_ptherm_bar,psi,etarho,etarho_cc, &
                     dt,dtold)

    use parallel
    use bl_prof_module
    use geometry, only: spherical, dr, r_start_coord, r_end_coord
    use bl_constants_module
    use probin_module, only: verbose
    use restrict_base_module, only: fill_ghost_base

    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(  out) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: w0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    real(kind=dp_t), intent(in   ) :: etarho(:,0:)
    real(kind=dp_t), intent(in   ) :: etarho_cc(:,0:)
    real(kind=dp_t), intent(inout) :: w0_force(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:), rho0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:), p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar_old(:,0:), gamma1bar_new(:,0:)
    real(kind=dp_t), intent(in   ) :: delta_p0_ptherm_bar(:,0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(:,0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold

    integer         :: r,n
    real(kind=dp_t) :: max_w0

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_w0")

    w0_force = ZERO

    if (spherical .eq. 0) then

       call make_w0_planar(nlevs,w0,w0_old,Sbar_in,p0_old,p0_new,gamma1bar_old, &
                           gamma1bar_new,delta_p0_ptherm_bar,psi,w0_force,dt,dtold)

    else

       do n=1,nlevs
          ! NOTE: need to fix this to put loop over nlevs within so we can call
          ! fill_ghost_base on w0 before computing w0_force
          call make_w0_spherical(n,w0(n,:),w0_old(n,0:),Sbar_in(n,0:), &
                                 rho0_old(n,:),rho0_new(n,:),p0_old(n,0:),p0_new(n,0:), &
                                 gamma1bar_old(n,0:),gamma1bar_new(n,0:), &
                                 delta_p0_ptherm_bar(n,0:), &
                                 etarho(n,0:),etarho_cc(n,0:), &
                                 w0_force(n,0:),dt,dtold)
       end do

    end if

    call fill_ghost_base(nlevs,w0_force,.true.)

    do n=1,nlevs
       max_w0 = zero
       do r=r_start_coord(n),r_end_coord(n)+1
          max_w0 = max(max_w0, abs(w0(n,r)))
       end do
       if (parallel_IOProcessor() .and. verbose .ge. 1) &
            write(6,*) '... max CFL of w0: ',max_w0 * dt / dr(n)
    end do

    call destroy(bpt)

  end subroutine make_w0

  subroutine make_w0_planar(nlevs,w0,w0_old,Sbar_in,p0_old,p0_new, &
                            gamma1bar_old,gamma1bar_new,delta_p0_ptherm_bar, &
                            psi,w0_force,dt,dtold)

    use geometry, only: nr_fine, r_start_coord, r_end_coord, dr, base_cutoff_density_coord
    use variables, only: rho_comp
    use bl_constants_module
    use probin_module, only: grav_const, dpdt_factor, base_cutoff_density
    use restrict_base_module, only: fill_ghost_base

    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(  out) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: w0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:), p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar_old(:,0:), gamma1bar_new(:,0:)
    real(kind=dp_t), intent(in   ) :: delta_p0_ptherm_bar(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    real(kind=dp_t), intent(inout) :: w0_force(:,0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold

    ! Local variables
    integer                      :: r, n, i, refrat
    real(kind=dp_t), allocatable :: w0_old_cen(:,:)
    real(kind=dp_t), allocatable :: w0_new_cen(:,:)
    real(kind=dp_t)              :: w0_avg, div_avg, dt_avg, gamma1bar_p0_avg
    real(kind=dp_t)              :: volume_discrepancy, offset

    ! Cell-centered
    allocate(w0_old_cen(nlevs,0:nr_fine-1))
    allocate(w0_new_cen(nlevs,0:nr_fine-1))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Multilevel Outline
    !
    ! Compute w0 at level 1 only
    ! do n=2,nlevs
    !   Compute w0 on edges at level n
    !   Obtain the starting value of w0 from the coarser grid
    !   do i=n-1,1,-1
    !     Restrict w0 from level n to level i
    !     Compare the difference between w0 at top of level n to the corresponding point
    !       on level i
    !     Offset the w0 on level i above this point
    !   end do
    ! end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Compute w0 on edges at level n
    do n=1,nlevs

       if (n .eq. 1) then
          ! Initialize new w0 at bottom of coarse base array to zero.
          w0(1,0) = ZERO
       else
          ! Obtain the starting value of w0 from the coarser grid
          w0(n,r_start_coord(n)) = w0(n-1,r_start_coord(n)/2)
       end if

       do r=r_start_coord(n)+1,r_end_coord(n)+1
          gamma1bar_p0_avg = (gamma1bar_old(n,r-1)+gamma1bar_new(n,r-1)) * &
               (p0_old(n,r-1)+p0_new(n,r-1))/4.0d0
       
          if (r-1 .lt. base_cutoff_density_coord(n)) then
             volume_discrepancy = dpdt_factor * delta_p0_ptherm_bar(n,r-1)/dt
          else
             volume_discrepancy = 0.0d0
          end if
          
          w0(n,r) = w0(n,r-1) + Sbar_in(n,r-1) * dr(n) &
               - ( (psi(n,r-1)+volume_discrepancy) / gamma1bar_p0_avg ) * dr(n)
       end do

       do i=n-1,1,-1

          refrat = 2**(n-i)

          ! Compare the difference between w0 at top of level n to the corresponding point
          !   on level i
          offset = w0(n,r_end_coord(n)+1) - w0(i,(r_end_coord(n)+1)/refrat)

          ! Restrict w0 from level n to level i
          do r=r_start_coord(n),r_end_coord(n)+1,refrat
             w0(i,r/refrat) = w0(n,r)
          end do

          ! Offset the w0 on level i above this point
          do r=(r_end_coord(n)+1)/refrat+1,r_end_coord(i)+1
             w0(i,r) = w0(i,r) + offset
          end do

       end do

    end do

    call fill_ghost_base(nlevs,w0,.false.)

    do n=1,nlevs
       
       ! Compute the forcing term in the base state velocity equation, - 1/rho0 grad pi0 
       dt_avg = HALF * (dt + dtold)
       do r=r_start_coord(n),r_end_coord(n)
          w0_old_cen(n,r) = HALF * (w0_old(n,r) + w0_old(n,r+1))
          w0_new_cen(n,r) = HALF * (w0(n,r) + w0(n,r+1))
          w0_avg = HALF * (dt * w0_old_cen(n,r) + dtold *  w0_new_cen(n,r)) / dt_avg
          div_avg = HALF * (dt * (w0_old(n,r+1)-w0_old(n,r)) + &
               dtold * (w0(n,r+1)-w0(n,r))) / dt_avg
          w0_force(n,r) = (w0_new_cen(n,r)-w0_old_cen(n,r))/dt_avg + w0_avg*div_avg/dr(n)
       end do

    end do

    deallocate(w0_old_cen,w0_new_cen)

  end subroutine make_w0_planar

  subroutine make_w0_spherical(n,w0,w0_old,Sbar_in,rho0,rho0_new,p0,p0_new, &
                               gamma1bar,gamma1bar_new,delta_p0_ptherm_bar, &
                               etarho,etarho_cc,w0_force,dt,dtold)

    use geometry, only: r_cc_loc, nr_fine, r_edge_loc, dr, r_end_coord
    use make_grav_module
    use cell_to_edge_module
    use bl_constants_module
    use fundamental_constants_module, only: Gconst
    use probin_module, only: dpdt_factor, base_cutoff_density

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(  out) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: w0_old(0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(0:)
    real(kind=dp_t), intent(in   ) :: rho0(0:),rho0_new(0:)
    real(kind=dp_t), intent(in   ) :: p0(0:),p0_new(0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(0:),gamma1bar_new(0:),delta_p0_ptherm_bar(0:)
    real(kind=dp_t), intent(in   ) :: etarho(0:),etarho_cc(0:)
    real(kind=dp_t), intent(inout) :: w0_force(0:)
    real(kind=dp_t), intent(in   ) :: dt,dtold

    ! Local variables
    integer                      :: r
    real(kind=dp_t), allocatable :: w0_old_cen(:)
    real(kind=dp_t), allocatable :: w0_new_cen(:)
    real(kind=dp_t), allocatable :: c(:),d(:),e(:),u(:),rhs(:)
    real(kind=dp_t), allocatable :: m(:)
    real(kind=dp_t)              :: w0_avg, div_avg, dt_avg
    real(kind=dp_t), allocatable :: w0_bar(:)
    real(kind=dp_t), allocatable :: grav_edge(:)

    real(kind=dp_t), parameter :: eps = 1.d-8

    real(kind=dp_t) :: dpdr

    real(kind=dp_t) :: volume_discrepancy

    real(kind=dp_t), allocatable :: gamma1bar_nph(:), rho0_nph(:), p0_nph(:)

    ! Cell-centered
    allocate(m         (0:nr_fine-1))
    allocate(w0_old_cen(0:nr_fine-1))
    allocate(w0_new_cen(0:nr_fine-1))

    allocate(gamma1bar_nph(0:nr_fine-1))
    allocate(     rho0_nph(0:nr_fine-1))
    allocate(       p0_nph(0:nr_fine-1))


    ! Edge-centered
    allocate(c(0:nr_fine),d(0:nr_fine),e(0:nr_fine),rhs(0:nr_fine),u(0:nr_fine))
    allocate(w0_bar(0:nr_fine))
    allocate(grav_edge(0:nr_fine))


    ! create time-centered base-state quantities
    do r = 0, r_end_coord(n)
       p0_nph(r)        = HALF*(p0(r)        + p0_new(r))
       rho0_nph(r)      = HALF*(rho0(r)      + rho0_new(r))
       gamma1bar_nph(r) = HALF*(gamma1bar(r) + gamma1bar_new(r))       
    enddo


    ! NOTE:  we first solve for the w0 resulting only from Sbar -- then we will
    ! solve for the update to w0.  We integrate d/dr (r^2 w0) = (r^2 Sbar)

    w0_bar = ZERO
    do r=1,r_end_coord(n)+1

       if (rho0(r-1) .gt. base_cutoff_density) then
          volume_discrepancy = dpdt_factor * delta_p0_ptherm_bar(r-1)/dt
       else
          volume_discrepancy = ZERO
       endif

       w0_bar(r) = w0_bar(r-1) + dr(n) * Sbar_in(r-1) * r_cc_loc(n,r-1)**2 - &
            dr(n)* volume_discrepancy * r_cc_loc(n,r-1)**2 / &
            (gamma1bar_nph(r-1)*p0_nph(r-1))

    end do

    do r = 1,r_end_coord(n)+1
       w0_bar(r) = w0_bar(r) / r_edge_loc(n,r)**2
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
       rhs(r) = four * dpdr * w0_bar(r) / r_edge_loc(n,r) - &
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

    w0(0) = ZERO
    do r=1,r_end_coord(n)+1
       w0(r) = u(r) / r_edge_loc(n,r)**2
    end do

    do r=0,r_end_coord(n)+1
       w0(r) = w0(r) + w0_bar(r)
    end do

    ! Compute the forcing term in the base state velocity equation, - 1/rho0 grad pi0 
    dt_avg = HALF * (dt + dtold)
    do r = 0,r_end_coord(n)
       w0_old_cen(r) = HALF * (w0_old(r) + w0_old(r+1))
       w0_new_cen(r) = HALF * (w0    (r) + w0    (r+1))
       w0_avg = HALF * (dt *  w0_old_cen(r)           + dtold *  w0_new_cen(r)  ) / dt_avg
       div_avg = HALF * (dt * (w0_old(r+1)-w0_old(r)) + dtold * (w0(r+1)-w0(r))) / dt_avg
       w0_force(r) = (w0_new_cen(r)-w0_old_cen(r)) / dt_avg + w0_avg * div_avg / dr(n)
    end do

    deallocate(c,d,e,rhs,u)
    deallocate(m)
    deallocate(w0_old_cen,w0_new_cen)

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
