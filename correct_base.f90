module correct_base_module

  use bl_types

  implicit none

  private

  public :: correct_base

contains

  subroutine correct_base(rho0_new,sprimebar)

    use bl_prof_module
    use geometry, only: spherical, nlevs_radial
    use restrict_base_module
    use bl_error_module

    real(kind=dp_t), intent(inout) :: rho0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: sprimebar(:,0:)
    
    ! local
    integer :: n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "correct_base")
    
    do n=1,nlevs_radial
       if (spherical .eq. 1) then
          call correct_base_state_spherical(n,rho0_new(n,0:),sprimebar(n,0:))
       else
          call bl_error("correct_base not defined for plane-parallel")
       end if
    enddo

    call restrict_base(rho0_new,.true.)
    call fill_ghost_base(rho0_new,.true.)

    call destroy(bpt)
       
  end subroutine correct_base

  subroutine correct_base_state_spherical(n,rho0_new,sprimebar)

    use geometry, only: anelastic_cutoff_coord, r_cc_loc, r_edge_loc, dr

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(inout) :: rho0_new(0:)
    real(kind=dp_t), intent(in   ) :: sprimebar(0:)
    
    ! Local variables
    integer :: r
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE RHO0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do r=0,anelastic_cutoff_coord(n)-1
       rho0_new(r) = rho0_new(r) + sprimebar(r)
    end do
    
  end subroutine correct_base_state_spherical
  
end module correct_base_module
