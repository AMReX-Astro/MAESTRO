! compute beta_0, the coefficient in our constraint equation,
! div{beta_0 U} = beta_0 S

module make_div_coeff_module

  use bl_types

  implicit none

  private

  public :: make_div_coeff

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine make_div_coeff(n,div_coeff,rho0,p0,gamma1bar,grav_center)

     use bl_constants_module
     use geometry, only: dr, anelastic_cutoff_coord, r_start_coord, r_end_coord

      integer        , intent(in   ) :: n
      real(kind=dp_t), intent(  out) :: div_coeff(0:)
      real(kind=dp_t), intent(in   ) :: rho0(0:), p0(0:), gamma1bar(0:)
      real(kind=dp_t), intent(in   ) :: grav_center(0:)

      integer :: r
      real(kind=dp_t) :: integral

      real(kind=dp_t) :: beta0_edge_lo, beta0_edge_hi

      real(kind=dp_t) :: lambda, mu, nu
      real(kind=dp_t) :: denom, coeff1, coeff2
      real(kind=dp_t) :: del,dpls,dmin,slim,sflag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute beta0 on the edges and average to the center      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      beta0_edge_lo = rho0(0)
      do r = r_start_coord(n),r_end_coord(n)

         ! compute the slopes
         if (r == r_start_coord(n) .or. r == r_end_coord(n)) then
            lambda = ZERO
            mu = ZERO
            nu = ZERO
         else

            del    = HALF* (rho0(r+1) - rho0(r-1))/dr(n)
            dpls   = TWO * (rho0(r+1) - rho0(r  ))/dr(n)
            dmin   = TWO * (rho0(r  ) - rho0(r-1))/dr(n)
            slim   = min(abs(dpls), abs(dmin))
            slim   = merge(slim, zero, dpls*dmin.gt.ZERO)
            sflag  = sign(ONE,del)
            lambda = sflag*min(slim,abs(del))

            del   = HALF* (gamma1bar(r+1) - gamma1bar(r-1))/dr(n)
            dpls  = TWO * (gamma1bar(r+1) - gamma1bar(r  ))/dr(n)
            dmin  = TWO * (gamma1bar(r  ) - gamma1bar(r-1))/dr(n)
            slim  = min(abs(dpls), abs(dmin))
            slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
            sflag = sign(ONE,del)
            mu    = sflag*min(slim,abs(del))

            del   = HALF* (  p0(r+1) -   p0(r-1))/dr(n)
            dpls  = TWO * (  p0(r+1) -   p0(r  ))/dr(n)
            dmin  = TWO * (  p0(r  ) -   p0(r-1))/dr(n)
            slim  = min(abs(dpls), abs(dmin))
            slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
            sflag = sign(ONE,del)
            nu    = sflag*min(slim,abs(del))

         endif

         if (r == r_start_coord(n) .or. r == r_end_coord(n)) then

            integral = abs(grav_center(r))*rho0(r)*dr(n)/(p0(r)*gamma1bar(r))

         else if (nu .eq. ZERO .or. mu .eq. ZERO .or. &
                  (nu*gamma1bar(r) - mu*p0(r)) .eq. ZERO .or. &
                  ((gamma1bar(r) + HALF*mu*dr(n))/ &
                   (gamma1bar(r) - HALF*mu*dr(n))) .le. ZERO .or. &
                  ((p0(r) + HALF*nu*dr(n))/ &
                   (p0(r) - HALF*nu*dr(n))) .le. ZERO) then

            integral = abs(grav_center(r))*rho0(r)*dr(n)/(p0(r)*gamma1bar(r))

         else 
            denom = nu*gamma1bar(r) - mu*p0(r)

            coeff1 = lambda*gamma1bar(r)/mu - rho0(r)
            coeff2 = lambda*p0(r)/nu - rho0(r)
 
            integral = (abs(grav_center(r))/denom)* &
                 (coeff1*log( (gamma1bar(r) + HALF*mu*dr(n))/ &
                              (gamma1bar(r) - HALF*mu*dr(n))) - &
                  coeff2*log( (p0(r) + HALF*nu*dr(n))/ &
                              (p0(r) - HALF*nu*dr(n))) )

         endif

         beta0_edge_hi = beta0_edge_lo * exp(-integral)

         div_coeff(r) = HALF*(beta0_edge_lo + beta0_edge_hi)

         beta0_edge_lo = beta0_edge_hi

      end do
      
      do r = anelastic_cutoff_coord(n),r_end_coord(n)
        div_coeff(r) = div_coeff(r-1) * (rho0(r)/rho0(r-1))
      end do 
   end subroutine make_div_coeff

end module make_div_coeff_module
