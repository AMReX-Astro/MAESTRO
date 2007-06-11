module make_div_coeff_module

  use bl_types
  use bl_constants_module
  use geometry

  implicit none

contains


!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine make_div_coeff (div_coeff,rho0,p0,gam1,grav_center,anelastic_cutoff)

      real(kind=dp_t), intent(  out) :: div_coeff(:)
      real(kind=dp_t), intent(in   ) :: rho0(:), p0(:), gam1(:)
      real(kind=dp_t), intent(in   ) :: grav_center(:), anelastic_cutoff

      integer :: j,ny,j_anel
      real(kind=dp_t) :: integral

      real(kind=dp_t) :: beta0_edge_lo, beta0_edge_hi

      real(kind=dp_t) :: lambda, mu, nu
      real(kind=dp_t) :: denom, coeff1, coeff2
      real(kind=dp_t) :: del,dpls,dmin,slim,sflag

      ny = size(div_coeff,dim=1)
      j_anel = ny

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     compute beta0 on the edges and average to the center      
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      beta0_edge_lo = rho0(1)
      do j = 1,ny

         ! compute the slopes
         if (j == 1 .or. j == ny) then
            lambda = ZERO
            mu = ZERO
            nu = ZERO
         else

            del    = HALF* (rho0(j+1) - rho0(j-1))/dr
            dpls   = TWO * (rho0(j+1) - rho0(j  ))/dr
            dmin   = TWO * (rho0(j  ) - rho0(j-1))/dr
            slim   = min(abs(dpls), abs(dmin))
            slim   = merge(slim, zero, dpls*dmin.gt.ZERO)
            sflag  = sign(ONE,del)
            lambda = sflag*min(slim,abs(del))

            del   = HALF* (gam1(j+1) - gam1(j-1))/dr
            dpls  = TWO * (gam1(j+1) - gam1(j  ))/dr
            dmin  = TWO * (gam1(j  ) - gam1(j-1))/dr
            slim  = min(abs(dpls), abs(dmin))
            slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
            sflag = sign(ONE,del)
            mu    = sflag*min(slim,abs(del))

            del   = HALF* (  p0(j+1) -   p0(j-1))/dr
            dpls  = TWO * (  p0(j+1) -   p0(j  ))/dr
            dmin  = TWO * (  p0(j  ) -   p0(j-1))/dr
            slim  = min(abs(dpls), abs(dmin))
            slim  = merge(slim, zero, dpls*dmin.gt.ZERO)
            sflag = sign(ONE,del)
            nu    = sflag*min(slim,abs(del))

         endif

         if (j == 1 .or. j == ny) then

            integral = abs(grav_center(j))*rho0(j)*dr/(p0(j)*gam1(j))

         else if (nu .eq. ZERO .or. mu .eq. ZERO .or. &
                  ((gam1(j) + HALF*mu*dr)/(gam1(j) - HALF*mu*dr)) .le. ZERO .or. &
                  ((p0(j) + HALF*nu*dr)/(p0(j) - HALF*nu*dr)) .le. ZERO) then

            integral = abs(grav_center(j))*rho0(j)*dr/(p0(j)*gam1(j))

         else 
            denom = nu*gam1(j) - mu*p0(j)

            coeff1 = lambda*gam1(j)/mu - rho0(j)
            coeff2 = lambda*p0(j)/nu - rho0(j)
 
            integral = (abs(grav_center(j))/denom)* &
                 (coeff1*log( (gam1(j) + HALF*mu*dr)/ &
                              (gam1(j) - HALF*mu*dr)) - &
                  coeff2*log( (p0(j) + HALF*nu*dr)/ &
                              (p0(j) - HALF*nu*dr)) )

         endif

         beta0_edge_hi = beta0_edge_lo * exp(-integral)

         div_coeff(j) = HALF*(beta0_edge_lo + beta0_edge_hi)


         if (rho0(j) .lt. anelastic_cutoff .and. j_anel .eq. ny) then
            j_anel = j
            exit
         end if
         
         beta0_edge_lo = beta0_edge_hi

      end do
      
      do j = j_anel,ny
        div_coeff(j) = div_coeff(j-1) * (rho0(j)/rho0(j-1))
      end do 
   end subroutine make_div_coeff

   subroutine put_1d_beta_on_edges (div_coeff_cell,div_coeff_edge)

      real(kind=dp_t), intent(in   ) :: div_coeff_cell(:)
      real(kind=dp_t), intent(  out) :: div_coeff_edge(:)

      integer :: j,ny
      ny = size(div_coeff_cell,dim=1)

      div_coeff_edge(   1) = div_coeff_cell(1)
      div_coeff_edge(   2) = HALF*(div_coeff_cell( 1) + div_coeff_cell(2))
      div_coeff_edge(ny  ) = HALF*(div_coeff_cell(ny) + div_coeff_cell(ny-1))
      div_coeff_edge(ny+1) = div_coeff_cell(ny)
      do j = 3,ny-1
        div_coeff_edge(j) = 7.d0/12.d0 * (div_coeff_cell(j  ) + div_coeff_cell(j-1)) &
                           -1.d0/12.d0 * (div_coeff_cell(j+1) + div_coeff_cell(j-2))
      end do

   end subroutine put_1d_beta_on_edges

end module make_div_coeff_module

