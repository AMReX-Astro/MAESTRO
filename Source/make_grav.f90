module make_grav_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use variables
  use geometry
  use probin_module, only: grav_const

  implicit none

  real(kind=dp_t), parameter :: Gconst = 6.6725985E-8_dp_t

contains

  subroutine make_grav_cell(grav_cell,rho0)

    ! compute the base state gravitational acceleration at the cell
    ! centers.  The base state uses 0-based indexing, so grav_cell 
    ! does too.

    real(kind=dp_t), intent(  out) :: grav_cell(0:)
    real(kind=dp_t), intent(in   ) :: rho0(0:)

    ! Local variables
    integer                      :: k,nz
    real(kind=dp_t), allocatable :: m(:)

    nz = size(grav_cell,dim=1)

    if (spherical .eq. 0) then

       grav_cell(:) = grav_const
       
    else

       allocate(m(0:nz-1))

       m(0) = FOUR3RD*M_PI*rho0(0)*z(0)**3
       grav_cell(0) = -Gconst * m(0) / z(0)**2

       do k = 1, nz-1
          ! the mass is defined at the cell-centers, so to compute the
          ! mass at the current center, we need to add the contribution of
          ! the upper half of the zone below us and the lower half of the
          ! current zone.
          m(k) = m(k-1) + FOUR3RD*M_PI*rho0(k-1)*(zl(k) -  z(k-1))*(zl(k)**2 + zl(k)* z(k-1) +  z(k-1)**2) &
                        + FOUR3RD*M_PI*rho0(k  )*( z(k) - zl(k  ))*( z(k)**2 +  z(k)*zl(k  ) + zl(k  )**2)
          grav_cell(k) = -Gconst * m(k) / z(k)**2
       enddo

       deallocate(m)

    end if

  end subroutine make_grav_cell

  subroutine make_grav_edge(grav_edge,rho0)

    ! compute the base state gravity at the cell edges (grav_edge(1)
    ! is the gravitational acceleration at the left edge of zone 1).
    ! The base state uses 0-based indexing, so grav_edge does too.

    real(kind=dp_t), intent(  out) :: grav_edge(0:)
    real(kind=dp_t), intent(in   ) :: rho0(0:)

      ! Local variables
      integer                      :: j,k,nz
      real(kind=dp_t)              :: mencl

      nz = size(grav_edge,dim=1)

      if (spherical .eq. 0) then

        grav_edge(:) = grav_const

      else

        grav_edge(0) = zero 
        do k = 1,nz-1

          mencl = zero 
          do j = 1, k
            mencl = mencl + FOUR3RD*M_PI * (zl(j) - zl(j-1)) * (zl(j)**2 + zl(j)*zl(j-1) + zl(j-1)**2) * rho0(j-1)
          end do

          grav_edge(k) = -Gconst * mencl / zl(k)**2
        end do

      end if

   end subroutine make_grav_edge

end module make_grav_module
