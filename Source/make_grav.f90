module make_grav_module
  
  use bl_types

  implicit none

  real(kind=dp_t), parameter :: Gconst = 6.6725985E-8_dp_t

  private

  public :: make_grav_cell, make_grav_edge

contains

  subroutine make_grav_cell(n,grav_cell,rho0)

    use bl_constants_module
    use geometry, only: spherical, nr, base_cc_loc, base_loedge_loc
    use probin_module, only: grav_const


    ! compute the base state gravitational acceleration at the cell
    ! centers.  The base state uses 0-based indexing, so grav_cell 
    ! does too.
    
    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(  out) :: grav_cell(0:)
    real(kind=dp_t), intent(in   ) :: rho0(0:)

    ! Local variables
    integer                      :: r
    real(kind=dp_t), allocatable :: m(:)

    if (spherical .eq. 0) then

       grav_cell(:) = grav_const
       
    else

       allocate(m(0:nr(n)-1))

       m(0) = FOUR3RD*M_PI*rho0(0)*base_cc_loc(n,0)**3
       grav_cell(0) = -Gconst * m(0) / base_cc_loc(n,0)**2

       do r = 1, nr(n)-1
          ! the mass is defined at the cell-centers, so to compute the
          ! mass at the current center, we need to add the contribution of
          ! the upper half of the zone below us and the lower half of the
          ! current zone.
          m(r) = m(r-1) + FOUR3RD*M_PI*rho0(r-1)*(base_loedge_loc(n,r) - &
               base_cc_loc(n,r-1))*(base_loedge_loc(n,r)**2 + &
               base_loedge_loc(n,r)* base_cc_loc(n,r-1) +  base_cc_loc(n,r-1)**2) &
               + FOUR3RD*M_PI*rho0(r  )*&
               ( base_cc_loc(n,r) - base_loedge_loc(n,r  ))*( base_cc_loc(n,r)**2 + &
               base_cc_loc(n,r)*base_loedge_loc(n,r  ) + base_loedge_loc(n,r  )**2)
          grav_cell(r) = -Gconst * m(r) / base_cc_loc(n,r)**2
       enddo

       deallocate(m)

    end if

  end subroutine make_grav_cell

  subroutine make_grav_edge(n,grav_edge,rho0)

  use bl_constants_module
  use geometry, only: spherical, nr, base_loedge_loc
  use probin_module, only: grav_const


    ! compute the base state gravity at the cell edges (grav_edge(1)
    ! is the gravitational acceleration at the left edge of zone 1).
    ! The base state uses 0-based indexing, so grav_edge does too.

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(  out) :: grav_edge(0:)
    real(kind=dp_t), intent(in   ) :: rho0(0:)

    ! Local variables
    integer                      :: r,r2
    real(kind=dp_t)              :: mencl
    
    if (spherical .eq. 0) then
       
       grav_edge(:) = grav_const
       
    else
       
       grav_edge(0) = zero 
       do r = 1,nr(n)-1
          
          mencl = zero 
          do r2 = 1,r
             mencl = mencl + FOUR3RD*M_PI * &
                  (base_loedge_loc(n,r2) - base_loedge_loc(n,r2-1)) &
                  * (base_loedge_loc(n,r2)**2 &
                  + base_loedge_loc(n,r2)*base_loedge_loc(n,r2-1) &
                  + base_loedge_loc(n,r2-1)**2) * rho0(r2-1)
          end do
          
          grav_edge(r) = -Gconst * mencl / base_loedge_loc(n,r)**2
       end do
       
    end if
    
  end subroutine make_grav_edge
  
end module make_grav_module
