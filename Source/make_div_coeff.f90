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
    integer :: r, n, i, refrat
    real(kind=dp_t) :: integral
    real(kind=dp_t) :: beta0_edge(nlevs,0:nr_fine)
    real(kind=dp_t) :: lambda, mu, nu
    real(kind=dp_t) :: denom, coeff1, coeff2
    real(kind=dp_t) :: del,dpls,dmin,slim,sflag
    real(kind=dp_t) :: offset,ghostval

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
    !   do i=n-1,1,-1
    !     Restrict beta0 at edges from level n to level i
    !     Recompute beta0 at centers at level i for cells that are covered by level n data
    !     Compare the difference between beta0 at the top of level n to the corresponding
    !      point on level i
    !     Offset the centered beta on level i above this point so the total integral 
    !      is consistent
    !     Redo the anelastic cutoff part
    !   end do
    ! end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Compute beta0 on edges and centers at level 1 only

    beta0_edge(1,0) = rho0(1,0)

    do r = r_start_coord(1),r_end_coord(1)

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
       beta0_edge(n,r_start_coord(n)) = beta0_edge(n-1,r_start_coord(n)/2)

       do r=r_start_coord(n),r_end_coord(n)

          ! Modify the slope calculation at the level edges to look at coarser data
          if (r .eq. r_start_coord(n)) then

             ghostval = (1.d0/3.d0)*rho0(n,r) + (2.d0/3.d0)*rho0(n-1,r/2-1)
             del    = HALF* (rho0(n,r+1) - ghostval   )/dr(n)
             dpls   = TWO * (rho0(n,r+1) - rho0(n,r  ))/dr(n)
             dmin   = TWO * (rho0(n,r  ) - ghostval   )/dr(n)
             slim   = min(abs(dpls), abs(dmin))
             slim   = merge(slim, zero, dpls*dmin.gt.ZERO)
             sflag  = sign(ONE,del)
             lambda = sflag*min(slim,abs(del))

             ghostval = (1.d0/3.d0)*gamma1bar(n,r) + (2.d0/3.d0)*gamma1bar(n-1,r/2-1)
             del   = HALF* (gamma1bar(n,r+1) - ghostval        )/dr(n)
             dpls  = TWO * (gamma1bar(n,r+1) - gamma1bar(n,r  ))/dr(n)
             dmin  = TWO * (gamma1bar(n,r  ) - ghostval        )/dr(n)
             slim  = min(abs(dpls), abs(dmin))
             slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
             sflag = sign(ONE,del)
             mu    = sflag*min(slim,abs(del))

             ghostval = (1.d0/3.d0)*p0(n,r) + (2.d0/3.d0)*p0(n-1,r/2-1)
             del   = HALF* (p0(n,r+1) - ghostval )/dr(n)
             dpls  = TWO * (p0(n,r+1) - p0(n,r  ))/dr(n)
             dmin  = TWO * (p0(n,r  ) - ghostval )/dr(n)
             slim  = min(abs(dpls), abs(dmin))
             slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
             sflag = sign(ONE,del)
             nu    = sflag*min(slim,abs(del))

          else if (r .eq. r_end_coord(n)) then

             ghostval = (1.d0/3.d0)*rho0(n,r) + (2.d0/3.d0)*rho0(n-1,(r+1)/2)
             del    = HALF* (ghostval    - rho0(n,r-1))/dr(n)
             dpls   = TWO * (ghostval    - rho0(n,r  ))/dr(n)
             dmin   = TWO * (rho0(n,r  ) - rho0(n,r-1))/dr(n)
             slim   = min(abs(dpls), abs(dmin))
             slim   = merge(slim, zero, dpls*dmin.gt.ZERO)
             sflag  = sign(ONE,del)
             lambda = sflag*min(slim,abs(del))

             ghostval = (1.d0/3.d0)*gamma1bar(n,r) + (2.d0/3.d0)*gamma1bar(n-1,(r+1)/2)
             del   = HALF* (ghostval         - gamma1bar(n,r-1))/dr(n)
             dpls  = TWO * (ghostval         - gamma1bar(n,r  ))/dr(n)
             dmin  = TWO * (gamma1bar(n,r  ) - gamma1bar(n,r-1))/dr(n)
             slim  = min(abs(dpls), abs(dmin))
             slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
             sflag = sign(ONE,del)
             mu    = sflag*min(slim,abs(del))

             ghostval = (1.d0/3.d0)*p0(n,r) + (2.d0/3.d0)*p0(n-1,(r+1)/2)
             del   = HALF* (ghostval  - p0(n,r-1))/dr(n)
             dpls  = TWO * (ghostval  - p0(n,r  ))/dr(n)
             dmin  = TWO * (p0(n,r  ) - p0(n,r-1))/dr(n)
             slim  = min(abs(dpls), abs(dmin))
             slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
             sflag = sign(ONE,del)
             nu    = sflag*min(slim,abs(del))

          else

             del    = HALF* (rho0(n,r+1) - rho0(n,r-1))/dr(n)
             dpls   = TWO * (rho0(n,r+1) - rho0(n,r  ))/dr(n)
             dmin   = TWO * (rho0(n,r  ) - rho0(n,r-1))/dr(n)
             slim   = min(abs(dpls), abs(dmin))
             slim   = merge(slim, zero, dpls*dmin.gt.ZERO)
             sflag  = sign(ONE,del)
             lambda = sflag*min(slim,abs(del))
             
             del   = HALF* (gamma1bar(n,r+1) - gamma1bar(n,r-1))/dr(n)
             dpls  = TWO * (gamma1bar(n,r+1) - gamma1bar(n,r  ))/dr(n)
             dmin  = TWO * (gamma1bar(n,r  ) - gamma1bar(n,r-1))/dr(n)
             slim  = min(abs(dpls), abs(dmin))
             slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
             sflag = sign(ONE,del)
             mu    = sflag*min(slim,abs(del))
             
             del   = HALF* (p0(n,r+1) - p0(n,r-1))/dr(n)
             dpls  = TWO * (p0(n,r+1) - p0(n,r  ))/dr(n)
             dmin  = TWO * (p0(n,r  ) - p0(n,r-1))/dr(n)
             slim  = min(abs(dpls), abs(dmin))
             slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
             sflag = sign(ONE,del)
             nu    = sflag*min(slim,abs(del))

          end if

          if (nu .eq. ZERO .or. mu .eq. ZERO .or. &
             (nu*gamma1bar(n,r) - mu*p0(n,r)) .eq. ZERO .or. &
             ((gamma1bar(n,r) + HALF*mu*dr(n))/ &
             (gamma1bar(n,r) - HALF*mu*dr(n))) .le. ZERO .or. &
             ((p0(n,r) + HALF*nu*dr(n))/ &
             (p0(n,r) - HALF*nu*dr(n))) .le. ZERO) then

             integral = abs(grav_center(n,r))*rho0(n,r)*dr(n)/(p0(n,r)*gamma1bar(n,r))

          else 

             denom = nu*gamma1bar(n,r) - mu*p0(n,r)
             coeff1 = lambda*gamma1bar(n,r)/mu - rho0(n,r)
             coeff2 = lambda*p0(n,r)/nu - rho0(n,r)
             
             integral = (abs(grav_center(n,r))/denom)* &
                  (coeff1*log( (gamma1bar(n,r) + HALF*mu*dr(n))/ &
                  (gamma1bar(n,r) - HALF*mu*dr(n))) - &
                  coeff2*log( (p0(n,r) + HALF*nu*dr(n))/ &
                  (p0(n,r) - HALF*nu*dr(n))) )
             
          endif

          beta0_edge(n,r+1) = beta0_edge(n,r) * exp(-integral)
          div_coeff(n,r) = HALF*(beta0_edge(n,r) + beta0_edge(n,r+1))

       end do

       do r = anelastic_cutoff_coord(n),r_end_coord(n)
          div_coeff(n,r) = div_coeff(n,r-1) * (rho0(n,r)/rho0(n,r-1))
       end do

       do i=n-1,1,-1

          refrat = 2**(n-i)

          ! Restrict beta0 at edges from level n to level i
          do r=r_start_coord(n),r_end_coord(n)+1,refrat
             beta0_edge(i,r/refrat) = beta0_edge(n,r)
          end do

          ! Recompute beta0 at centers at level i for cells that are covered by level n data
          do r=r_start_coord(n),r_end_coord(n)-refrat+1,refrat
             div_coeff(i,r/refrat) = HALF*(beta0_edge(i,r/refrat) + beta0_edge(i,r/refrat+1))
          end do

          ! Compare the difference between beta0 at the top of level n to the corresponding
          !  point on level i
          offset = beta0_edge(n,r_end_coord(n)+1) - beta0_edge(i,(r_end_coord(n)+1)/refrat)

          ! Offset the centered beta on level i above this point so the total integral 
          !  is consistent
          do r=(r_end_coord(n)+1)/refrat,r_end_coord(i)
             div_coeff(i,r) = div_coeff(i,r) + offset
          end do

          ! Redo the anelastic cutoff part
          do r = anelastic_cutoff_coord(i),r_end_coord(i)
             div_coeff(i,r) = div_coeff(i,r-1) * (rho0(i,r)/rho0(i,r-1))
          end do

       end do

    end do

  end subroutine make_div_coeff

end module make_div_coeff_module
