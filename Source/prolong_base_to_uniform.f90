! Take a multilevel base state quantity (for plane-parallel only),
! and convert it to a uniformly gridded, single-level array at
! the finest base state resolution.

! this works for cell-centered data only

module prolong_base_to_uniform_module

  use bl_types

  implicit none

  private

  public :: prolong_base_to_uniform

contains

  subroutine prolong_base_to_uniform(base_ml, base_fine)

    use bl_prof_module
    use bl_error_module
    use geometry, only: r_start_coord, r_end_coord, numdisjointchunks, nlevs_radial, nr_fine

    real(kind=dp_t), intent(in   ) :: base_ml(:,0:)
    real(kind=dp_t), intent(inout) :: base_fine(0:)

    ! local
    integer :: n, j, r
    logical, allocatable :: imask_fine(:)
    integer :: r1
    
    type(bl_prof_timer), save :: bpt

    call build(bpt, "prolong_base")

    ! the mask array will keep track of whether we've filled in data
    ! in a corresponding radial bin.  .false. indicates that we've
    ! already output there.
    allocate(imask_fine(0:nr_fine-1))
    imask_fine(:) = .true.

    ! r1 is the factor between the current level grid spacing and the
    ! FINEST level
    r1 = 1

    do n = nlevs_radial, 1, -1
       do j = 1,numdisjointchunks(n)
          do r = r_start_coord(n,j), r_end_coord(n,j)

             if (any(imask_fine(r*r1:(r+1)*r1-1) ) ) then
                 base_fine(r*r1:(r+1)*r1-1) = base_ml(n,r)
                imask_fine(r*r1:(r+1)*r1-1) = .false.
             endif
             
          enddo
       enddo
       
       ! update r1 for the next coarsest level -- assume a jump by
       ! factor of 2
       r1 = r1*2

    enddo

    ! check to make sure that no mask values are still true    
    if (any(imask_fine)) then
       call bl_error("ERROR: unfilled cells in prolong_base_to_uniform")
    endif

    call destroy(bpt)

  end subroutine prolong_base_to_uniform

end module prolong_base_to_uniform_module
    
