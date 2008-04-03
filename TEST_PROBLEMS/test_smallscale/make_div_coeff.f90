module make_div_coeff_module

  use bl_types

  implicit none

  private

  public :: make_div_coeff

contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_div_coeff(n,div_coeff,rho0,p0,gamma1bar,grav_center)

    use bl_constants_module
    use geometry, only: dr
    use probin_module, only: anelastic_cutoff

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(  out) :: div_coeff(0:)
    real(kind=dp_t), intent(in   ) :: rho0(0:), p0(0:), gamma1bar(0:)
    real(kind=dp_t), intent(in   ) :: grav_center(0:)

    integer :: r,nr,r_anel
    real(kind=dp_t) :: integral

    real(kind=dp_t) :: beta0_edge_lo, beta0_edge_hi

    real(kind=dp_t) :: lambda, mu, nu
    real(kind=dp_t) :: denom, coeff1, coeff2
    real(kind=dp_t) :: del,dpls,dmin,slim,sflag

    nr = size(div_coeff,dim=1)
    r_anel = nr-1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! compute beta0 on the edges and average to the center      
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    beta0_edge_lo = rho0(0)
    do r = 0,nr-1

       ! compute the slopes
       if (r == 0 .or. r == nr-1) then
          lambda = ZERO
          mu = ZERO
          nu = ZERO
       else

          del    = HALF* (rho0(r+1) - rho0(r-1))/dr(1)
          dpls   = TWO * (rho0(r+1) - rho0(r  ))/dr(1)
          dmin   = TWO * (rho0(r  ) - rho0(r-1))/dr(1)
          slim   = min(abs(dpls), abs(dmin))
          slim   = merge(slim, zero, dpls*dmin.gt.ZERO)
          sflag  = sign(ONE,del)
          lambda = sflag*min(slim,abs(del))

          del   = HALF* (gamma1bar(r+1) - gamma1bar(r-1))/dr(1)
          dpls  = TWO * (gamma1bar(r+1) - gamma1bar(r  ))/dr(1)
          dmin  = TWO * (gamma1bar(r  ) - gamma1bar(r-1))/dr(1)
          slim  = min(abs(dpls), abs(dmin))
          slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
          sflag = sign(ONE,del)
          mu    = sflag*min(slim,abs(del))

          del   = HALF* (  p0(r+1) -   p0(r-1))/dr(1)
          dpls  = TWO * (  p0(r+1) -   p0(r  ))/dr(1)
          dmin  = TWO * (  p0(r  ) -   p0(r-1))/dr(1)
          slim  = min(abs(dpls), abs(dmin))
          slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
          sflag = sign(ONE,del)
          nu    = sflag*min(slim,abs(del))

       endif

       if (r == 0 .or. r == nr-1) then

          integral = abs(grav_center(r))*rho0(r)*dr(1)/(p0(r)*gamma1bar(r))

       else if (nu .eq. ZERO .or. mu .eq. ZERO .or. &
            ((gamma1bar(r) + HALF*mu*dr(1))/(gamma1bar(r) - HALF*mu*dr(1))) .le. ZERO .or. &
            ((p0(r) + HALF*nu*dr(1))/(p0(r) - HALF*nu*dr(1))) .le. ZERO) then

          integral = abs(grav_center(r))*rho0(r)*dr(1)/(p0(r)*gamma1bar(r))

       else 
          denom = nu*gamma1bar(r) - mu*p0(r)

          coeff1 = lambda*gamma1bar(r)/mu - rho0(r)
          coeff2 = lambda*p0(r)/nu - rho0(r)

          integral = (abs(grav_center(r))/denom)* &
               (coeff1*log( (gamma1bar(r) + HALF*mu*dr(1))/ &
               (gamma1bar(r) - HALF*mu*dr(1))) - &
               coeff2*log( (p0(r) + HALF*nu*dr(1))/ &
               (p0(r) - HALF*nu*dr(1))) )

       endif

       beta0_edge_hi = beta0_edge_lo * exp(-integral)

       div_coeff(r) = HALF*(beta0_edge_lo + beta0_edge_hi)


       if (rho0(r) .lt. anelastic_cutoff .and. r_anel .eq. nr-1) then
          r_anel = r
          exit
       end if

       beta0_edge_lo = beta0_edge_hi

    end do

    !      do r = r_anel,nr-1
    !        div_coeff(r) = div_coeff(r-1) * (rho0(r)/rho0(r-1))
    !      end do 

    ! HACK HACK HACK
    do r = 0,nr-1
       div_coeff(r) = 1.d0
    end do

  end subroutine make_div_coeff

end module make_div_coeff_module

