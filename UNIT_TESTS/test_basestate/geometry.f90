! a module for storing the geometric information so we don't have to pass it
!
! This module provides the coordinate value for the left edge of a base-state
! zone (r_edge_loc) and the zone center (r_cc_loc).  As always, it is assumed that 
! the base state arrays begin with index 0, not 1.
!
! Note: this version is specific to test_basestate -- it must deal with a 1-d
! grid.

module geometry

  use bl_types
  use ml_layout_module

  implicit none

  integer   , save :: spherical
  real(dp_t), save :: center(3)
  integer   , save :: nr_fine
  real(dp_t), allocatable, save :: dr(:), r_cc_loc(:,:), r_edge_loc(:,:)
  integer   , allocatable, save :: r_start_coord(:), r_end_coord(:), nr(:)
  integer   , allocatable, save :: anelastic_cutoff_coord(:)
  integer   , allocatable, save :: base_cutoff_density_coord(:)
  integer   , allocatable, save :: burning_cutoff_density_coord(:)


  private

  public :: spherical, center, nr_fine
  public :: dr, r_cc_loc, r_edge_loc
  public :: r_start_coord, r_end_coord, nr
  public :: anelastic_cutoff_coord
  public :: base_cutoff_density_coord
  public :: burning_cutoff_density_coord
  
  public :: init_spherical, init_geometry, destroy_geometry

contains

  subroutine init_spherical(spherical_in)

    integer   , intent(in) :: spherical_in

    spherical = spherical_in

  end subroutine init_spherical

  subroutine init_geometry(center_in, nr_in, dr_in)

    use bl_constants_module
    use probin_module, only: prob_lo_x

    real(dp_t)     , intent(in)    :: center_in(:)
    integer        , intent(in)    :: nr_in
    real(dp_t)     , intent(in)    :: dr_in

    ! Local variables
    integer :: i,n,nlevs
    real(dp_t) :: base_lo

    nlevs = 1

    base_lo = prob_lo_x

    allocate(dr(nlevs))
    allocate(r_cc_loc(nlevs,0:nr_in-1))
    allocate(r_edge_loc(nlevs,0:nr_in))
    allocate(r_start_coord(nlevs))
    allocate(r_end_coord(nlevs))
    allocate(nr(nlevs))
    allocate(anelastic_cutoff_coord(nlevs))
    allocate(base_cutoff_density_coord(nlevs))
    allocate(burning_cutoff_density_coord(nlevs))

    nr(1) = nr_in
    dr(1) = dr_in

    center(:) = center_in(:)

    do n=1,nlevs
       do i = 0,nr(n)-1
          r_cc_loc(n,i) = base_lo + (dble(i)+HALF)*dr(n)
       end do
       do i = 0,nr(n)
          r_edge_loc(n,i) = base_lo + (dble(i))*dr(n)
       end do
    enddo

    r_start_coord(nlevs) = 0
    r_end_coord(nlevs) = nr_in-1
    
  end subroutine init_geometry

  subroutine destroy_geometry()

    deallocate(r_cc_loc,r_edge_loc)
    deallocate(dr,anelastic_cutoff_coord)

  end subroutine destroy_geometry

end module geometry
