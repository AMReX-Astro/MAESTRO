module make_div_coeff_module

  use bl_types
  use bl_constants_module

  implicit none

contains


!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine make_div_coeff (div_coeff,rho0,p0,gam1,grav_edge,dy,anelastic_cutoff)

      real(kind=dp_t), intent(  out) :: div_coeff(:)
      real(kind=dp_t), intent(in   ) :: rho0(:), p0(:), gam1(:)
      real(kind=dp_t), intent(in   ) :: grav_edge(:), dy, anelastic_cutoff

      integer :: j,ny,j_anel
      real(kind=dp_t) :: integral,rho0_lo,rho0_hi
      real(kind=dp_t) :: rho0_edge,gam1_edge,p0_edge

      ny = size(div_coeff,dim=1)
      j_anel = ny

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     CREATING CELL FIRST THEN AVERAGING ONTO HALF - DIFFERENT INTEGRAL
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      div_coeff(1) = rho0(1)
      do j = 2,ny
         rho0_edge = HALF * (rho0(j) + rho0(j-1))
         gam1_edge = HALF * (gam1(j) + gam1(j-1))
           p0_edge = HALF * (  p0(j) +   p0(j-1))
         integral  = rho0_edge * abs(grav_edge(j)) * dy / (gam1_edge * p0_edge)

         div_coeff(j) = div_coeff(j-1) * exp(-integral)
         if (rho0(j) .lt. anelastic_cutoff .and. j_anel .eq. ny) j_anel = j
      end do

      do j = j_anel,ny
        div_coeff(j) = div_coeff(j-1) * (rho0(j)/rho0(j-1))
      end do
!     print *,'SETTING BETA TO RHO0 STARTING AT  ',j_anel

   end subroutine make_div_coeff

   subroutine put_beta_on_edges (div_coeff_cell,div_coeff_edge)

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
      div_coeff_edge(ny+1) = div_coeff_edge(ny)

   end subroutine put_beta_on_edges

end module make_div_coeff_module

