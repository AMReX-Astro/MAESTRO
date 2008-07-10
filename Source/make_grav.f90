module make_grav_module
  
  use bl_types

  implicit none

  private

  public :: make_grav_cell, make_grav_edge

contains

  subroutine make_grav_cell(n,grav_cell,rho0)

    use bl_constants_module
    use geometry, only: spherical, nr_fine, r_cc_loc, r_edge_loc, r_end_coord
    use probin_module, only: grav_const, base_cutoff_density
    use fundamental_constants_module, only: Gconst

    ! compute the base state gravitational acceleration at the cell
    ! centers.  The base state uses 0-based indexing, so grav_cell 
    ! does too.
    
    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(  out) :: grav_cell(0:)
    real(kind=dp_t), intent(in   ) :: rho0(0:)

    ! Local variables
    integer                      :: r
    real(kind=dp_t), allocatable :: m(:)
    real(kind=dp_t) :: term1, term2

    if (spherical .eq. 0) then

       grav_cell = grav_const
       
    else

       allocate(m(0:nr_fine-1))

       m(0) = FOUR3RD*M_PI*rho0(0)*r_cc_loc(n,0)**3
       grav_cell(0) = -Gconst * m(0) / r_cc_loc(n,0)**2

       do r = 1, r_end_coord(n)
          ! the mass is defined at the cell-centers, so to compute the
          ! mass at the current center, we need to add the contribution of
          ! the upper half of the zone below us and the lower half of the
          ! current zone.

          ! don't add any contributions from outside the star -- i.e.
          ! rho < base_cutoff_density
          if (rho0(r-1) > base_cutoff_density) then
             term1 = FOUR3RD*M_PI*rho0(r-1) * &
               (r_edge_loc(n,r) - r_cc_loc(n,r-1)) * &
               (r_edge_loc(n,r)**2 + &
                r_edge_loc(n,r)*r_cc_loc(n,r-1) + &
                r_cc_loc(n,r-1)**2)
          else
             term1 = ZERO
          endif

          if (rho0(r) > base_cutoff_density) then
             term2 = FOUR3RD*M_PI*rho0(r  )*&
               (r_cc_loc(n,r) - r_edge_loc(n,r  )) * &
               (r_cc_loc(n,r)**2 + &
                r_cc_loc(n,r)*r_edge_loc(n,r  ) + &
                r_edge_loc(n,r  )**2)          
          else
             term2 = ZERO
          endif

          m(r) = m(r-1) + term1 + term2

          grav_cell(r) = -Gconst * m(r) / r_cc_loc(n,r)**2
       enddo

       deallocate(m)

    end if

  end subroutine make_grav_cell

  subroutine make_grav_edge(n,grav_edge,rho0)

  use bl_constants_module
  use geometry, only: spherical, r_edge_loc, r_end_coord
  use probin_module, only: grav_const, base_cutoff_density
  use fundamental_constants_module, only: Gconst

    ! compute the base state gravity at the cell edges (grav_edge(1)
    ! is the gravitational acceleration at the left edge of zone 1).
    ! The base state uses 0-based indexing, so grav_edge does too.

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(  out) :: grav_edge(0:)
    real(kind=dp_t), intent(in   ) :: rho0(0:)

    ! Local variables
    integer                      :: r
    real(kind=dp_t)              :: mencl
    
    if (spherical .eq. 0) then
       
       grav_edge = grav_const
       
    else
       
       grav_edge(0) = zero 
       mencl = ZERO

       do r = 1, r_end_coord(n)

          ! only add to the enclosed mass if the density is 
          ! > base_cutoff_density
          if (rho0(r-1) > base_cutoff_density) then
             mencl = mencl + FOUR3RD*M_PI * &
                  (r_edge_loc(n,r) - r_edge_loc(n,r-1)) * &
                  (r_edge_loc(n,r)**2 + &
                   r_edge_loc(n,r)*r_edge_loc(n,r-1) + &
                   r_edge_loc(n,r-1)**2) * rho0(r-1)
          endif
          
          grav_edge(r) = -Gconst * mencl / r_edge_loc(n,r)**2
       end do
       
    end if
    
  end subroutine make_grav_edge
  
end module make_grav_module
