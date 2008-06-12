! compute beta_0, the coefficient in our constraint equation,
! div{beta_0 U} = beta_0 S

module make_div_coeff_module

  use bl_types

  implicit none

  private

  public :: make_div_coeff

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_div_coeff(nlevs,div_coeff,rho0,p0,gamma1bar,grav_center)

    use bl_constants_module
    use geometry, only: nr_fine, dr, anelastic_cutoff_coord, r_start_coord, r_end_coord

    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(  out) :: div_coeff(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0(:,0:), p0(:,0:), gamma1bar(:,0:), grav_center(:,0:)

    ! local
    integer :: r, n, i
    real(kind=dp_t) :: integral
    real(kind=dp_t) :: beta0_edge(nlevs,0:nr_fine)
    real(kind=dp_t) :: lambda, mu, nu
    real(kind=dp_t) :: denom, coeff1, coeff2
    real(kind=dp_t) :: del,dpls,dmin,slim,sflag

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Compute beta0 on the edges and average to the center      
    !
    ! Multilevel Outline:
    !
    ! First, compute beta0 on edges and centers at level 1 only
    ! do n=2,nlevs
    !   Compute beta0 on edges and centers at level n
    !   Obtain the starting value of beta0_edge_lo from the coarser grid
    !   Modify the slope calculation at the level edges to look at coarser data
    !   do i=n,2,-1
    !     Restrict beta0 at edges from level i to level i-1
    !     Recompute beta0 at centers at level i-1 for cells that are covered by level i data
    !     Compare the difference between beta0 at the top of level i to the corresponding
    !      point on level i-1
    !     Offset the centered beta on level i-1 above this point so the total integral 
    !      is consistent
    !   end do
    ! end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Compute beta0 on edges and centers at level 1 only

    beta0_edge(1,0) = rho0(1,0)

    do r = r_start_coord(1),r_end_coord(1)

       ! compute the slopes
       if (r == r_start_coord(1) .or. r == r_end_coord(1)) then

          lambda = ZERO
          mu = ZERO
          nu = ZERO

       else

          del    = HALF* (rho0(1,r+1) - rho0(1,r-1))/dr(1)
          dpls   = TWO * (rho0(1,r+1) - rho0(1,r  ))/dr(1)
          dmin   = TWO * (rho0(1,r  ) - rho0(1,r-1))/dr(1)
          slim   = min(abs(dpls), abs(dmin))
          slim   = merge(slim, zero, dpls*dmin.gt.ZERO)
          sflag  = sign(ONE,del)
          lambda = sflag*min(slim,abs(del))

          del   = HALF* (gamma1bar(1,r+1) - gamma1bar(1,r-1))/dr(1)
          dpls  = TWO * (gamma1bar(1,r+1) - gamma1bar(1,r  ))/dr(1)
          dmin  = TWO * (gamma1bar(1,r  ) - gamma1bar(1,r-1))/dr(1)
          slim  = min(abs(dpls), abs(dmin))
          slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
          sflag = sign(ONE,del)
          mu    = sflag*min(slim,abs(del))

          del   = HALF* (p0(1,r+1) - p0(1,r-1))/dr(1)
          dpls  = TWO * (p0(1,r+1) - p0(1,r  ))/dr(1)
          dmin  = TWO * (p0(1,r  ) - p0(1,r-1))/dr(1)
          slim  = min(abs(dpls), abs(dmin))
          slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
          sflag = sign(ONE,del)
          nu    = sflag*min(slim,abs(del))

       endif

       if (r == r_start_coord(1) .or. r == r_end_coord(1)) then

          integral = abs(grav_center(1,r))*rho0(1,r)*dr(1)/(p0(1,r)*gamma1bar(1,r))

       else if (nu .eq. ZERO .or. mu .eq. ZERO .or. &
               (nu*gamma1bar(1,r) - mu*p0(1,r)) .eq. ZERO .or. &
               ((gamma1bar(1,r) + HALF*mu*dr(1))/ &
               (gamma1bar(1,r) - HALF*mu*dr(1))) .le. ZERO .or. &
               ((p0(1,r) + HALF*nu*dr(1))/ &
               (p0(1,r) - HALF*nu*dr(1))) .le. ZERO) then

          integral = abs(grav_center(1,r))*rho0(1,r)*dr(1)/(p0(1,r)*gamma1bar(1,r))

       else 

          denom = nu*gamma1bar(1,r) - mu*p0(1,r)
          coeff1 = lambda*gamma1bar(1,r)/mu - rho0(1,r)
          coeff2 = lambda*p0(1,r)/nu - rho0(1,r)

          integral = (abs(grav_center(1,r))/denom)* &
               (coeff1*log( (gamma1bar(1,r) + HALF*mu*dr(1))/ &
               (gamma1bar(1,r) - HALF*mu*dr(1))) - &
               coeff2*log( (p0(1,r) + HALF*nu*dr(1))/ &
               (p0(1,r) - HALF*nu*dr(1))) )

       endif

       beta0_edge(1,r+1) = beta0_edge(1,r) * exp(-integral)
       div_coeff(1,r) = HALF*(beta0_edge(1,r) + beta0_edge(1,r+1))

    end do

    do r = anelastic_cutoff_coord(1),r_end_coord(1)
       div_coeff(1,r) = div_coeff(1,r-1) * (rho0(1,r)/rho0(1,r-1))
    end do

    do n=2,nlevs

       ! Compute beta0 on edges and centers at level n
       ! Obtain the starting value of beta0_edge_lo from the coarser grid
       ! Modify the slope calculation at the level edges to look at coarser data



       do i=n,2,-1

          ! Restrict beta0 at edges from level i to level i-1
          ! Recompute beta0 at centers at level i-1 for cells that are covered by level i data
          ! Compare the difference between beta0 at the top of level i to the corresponding
          !  point on level i-1
          ! Offset the centered beta on level i-1 above this point so the total integral 
          !  is consistent



       end do

    end do

  end subroutine make_div_coeff

end module make_div_coeff_module
