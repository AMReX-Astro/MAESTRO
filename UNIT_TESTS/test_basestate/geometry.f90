! a module for storing the geometric information so we don't have to pass it
!
! This module provides the coordinate value for the left edge of a base-state
! zone (base_loedge_loc) and the zone center (base_cc_loc).  As always, it is assumed that 
! the base state arrays begin with index 0, not 1.

module geometry

  use bl_types
  use ml_layout_module

  implicit none

  integer   , save :: spherical
  real(dp_t), save :: center(3)
  real(dp_t), allocatable, save :: dr(:)
  integer   , allocatable, save :: nr(:)
  real(dp_t), allocatable, save :: base_cc_loc(:,:), base_loedge_loc(:,:)

  private

  public :: spherical, center, dr, base_cc_loc, base_loedge_loc, nr
  public :: init_spherical, init_geometry, destroy_geometry

contains

  subroutine init_spherical(spherical_in)

    integer   , intent(in) :: spherical_in

    spherical = spherical_in

  end subroutine init_spherical

  subroutine init_geometry(center_in, nr_in, dr_in)

    use bl_constants_module

    real(dp_t)     , intent(in)    :: center_in(:)
    integer        , intent(in)    :: nr_in
    real(dp_t)     , intent(in)    :: dr_in

    ! Local variables
    integer :: i,n,nlevs

    nlevs = 1

    allocate(    base_cc_loc(nlevs,0:nr_in-1))
    allocate(base_loedge_loc(nlevs,0:nr_in))

    allocate(dr(nlevs))
    allocate(nr(nlevs))

    nr(1) = nr_in
    dr(1) = dr_in

    center(:) = center_in(:)

    do n=1,nlevs
       do i = 0,nr(n)-1
          base_cc_loc(n,i) = (dble(i)+HALF)*dr(n)
       end do
       do i = 0,nr(n)
          base_loedge_loc(n,i) = (dble(i))*dr(n)
       end do
    enddo

  end subroutine init_geometry

  subroutine destroy_geometry()

    deallocate(base_cc_loc,base_loedge_loc)
    deallocate(dr)

  end subroutine destroy_geometry

end module geometry
