module correct_base_module

  use bl_types

  implicit none

  private

  public :: correct_base

contains

  subroutine correct_base(s0_new,sprimebar)

    use bl_prof_module
    use geometry, only: spherical
    use restrict_base_module
    use bl_error_module

    real(kind=dp_t), intent(inout) :: s0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: sprimebar(:,0:)
    
    ! local
    type(bl_prof_timer), save :: bpt

    call build(bpt, "correct_base")
    
    if (spherical .eq. 1) then
       call correct_base_state_spherical(s0_new(1,:),sprimebar(1,:))
    else
       call bl_error("correct_base not defined for plane-parallel")
    end if

    call restrict_base(s0_new,.true.)
    call fill_ghost_base(s0_new,.true.)

    call destroy(bpt)
       
  end subroutine correct_base

  subroutine correct_base_state_spherical(s0_new,sprimebar)

    use geometry, only: anelastic_cutoff_coord, r_cc_loc, r_edge_loc, dr

    real(kind=dp_t), intent(inout) :: s0_new(0:)
    real(kind=dp_t), intent(in   ) :: sprimebar(0:)
    
    ! Local variables
    integer :: r

    do r=0,anelastic_cutoff_coord(1)-1
       s0_new(r) = s0_new(r) + sprimebar(r)
    end do
    
  end subroutine correct_base_state_spherical
  
end module correct_base_module
