module make_grav_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use variables
  use geometry

  implicit none

contains

   subroutine make_grav_cell(grav_cell,rho0)

      real(kind=dp_t), intent(  out) :: grav_cell(:)
      real(kind=dp_t), intent(in   ) :: rho0(:)

      ! Local variables
      integer                      :: k,nz
      real(kind=dp_t), allocatable :: m(:)

      real(kind=dp_t), parameter :: Gconst = 6.6725985E-8_dp_t
      real(kind=dp_t), parameter ::     pi = 3.141592653589793238_dp_t
      real(kind=dp_t), parameter :: fourthirds = 4.0_dp_t/3.0_dp_t

      nz = size(grav_cell,dim=1)

      if (spherical .eq. 0) then

        grav_cell(:) = -1.5d10

      else

        allocate(m(nz))

        m(1) = fourthirds*pi*rho0(1)*z(1)**3
        grav_cell(1) = -Gconst * m(1) / z(1)**2
        do k = 2, nz
           ! the mass is defined at the cell-centers, so to compute the
           ! mass at the current center, we need to add the contribution of
           ! the upper half of the zone below us and the lower half of the
           ! current zone.
           m(k) = m(k-1) + fourthirds*pi*rho0(k-1)*(zl(k)**3 -  z(k-1)**3) &
                         + fourthirds*pi*rho0(k  )*( z(k)**3 - zl(k  )**3)
           grav_cell(k) = -Gconst * m(k) / z(k)**2
        enddo

        deallocate(m)

      end if

   end subroutine make_grav_cell

   subroutine make_grav_edge(grav_edge,rho0)

      real(kind=dp_t), intent(  out) :: grav_edge(:)
      real(kind=dp_t), intent(in   ) :: rho0(:)

      ! Local variables
      integer                      :: j,k,nz
      real(kind=dp_t)              :: mencl

      real(kind=dp_t), parameter :: Gconst = 6.6725985E-8_dp_t
      real(kind=dp_t), parameter ::     pi = 3.141592653589793238_dp_t
      real(kind=dp_t), parameter :: fourthirds = 4.0_dp_t/3.0_dp_t

      nz = size(grav_edge,dim=1)

      if (spherical .eq. 0) then

        grav_edge(:) = -1.5d10

      else

        grav_edge(1) = zero 
        do k = 2,nz
          mencl = zero 
          do j = 2, k
            mencl = mencl + fourthirds * pi * (zl(j)**3 - zl(j-1)**3) * rho0(j-1)
          end do
          grav_edge(k) = -Gconst * mencl / zl(k)**2
        end do

      end if

   end subroutine make_grav_edge

end module make_grav_module
