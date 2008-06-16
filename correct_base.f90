module correct_base_module

  use bl_types

  implicit none

  private

  public :: correct_base

contains

  subroutine correct_base(nlevs,rho0_new,etarho,dz,dt)

    use bl_prof_module
    use geometry, only: spherical
    use restrict_base_module

    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(inout) :: rho0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: etarho(:,0:)
    real(kind=dp_t), intent(in   ) :: dz(:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! local
    integer :: n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "correct_base")
    
    do n=1,nlevs
       if (spherical .eq. 1) then
          ! spherical is not updated.
       else
          call correct_base_state_planar(n,rho0_new(n,0:),etarho(n,0:),dz(n),dt)
       end if
    enddo

    call restrict_base(nlevs,rho0_new,.true.)

    call destroy(bpt)
       
  end subroutine correct_base

  subroutine correct_base_state_planar(n,rho0_new,etarho,dz,dt)

    use geometry, only: anelastic_cutoff_coord, r_start_coord, r_end_coord

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(inout) :: rho0_new(0:)
    real(kind=dp_t), intent(in   ) :: etarho(0:)
    real(kind=dp_t), intent(in   ) :: dz,dt
    
    ! Local variables
    integer :: r
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE RHO0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do r=r_start_coord(n),anelastic_cutoff_coord(n)-1
       rho0_new(r) = rho0_new(r) - dt/dz*(etarho(r+1) - etarho(r))
    end do
    
  end subroutine correct_base_state_planar
  
end module correct_base_module
