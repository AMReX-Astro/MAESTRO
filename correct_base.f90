module correct_base_module

  use bl_types

  implicit none

  private

  public :: correct_base

contains

  subroutine correct_base(rho0_new,div_etarho,dt)

    use bl_prof_module
    use geometry, only: spherical, nlevs
    use restrict_base_module

    real(kind=dp_t), intent(inout) :: rho0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: div_etarho(:,0:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! local
    integer :: n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "correct_base")
    
    do n=1,nlevs
       if (spherical .eq. 1) then
          call correct_base_state_spherical(n,rho0_new(n,0:),div_etarho(n,0:),dt)
       else
          call correct_base_state_planar(n,rho0_new(n,0:),div_etarho(n,0:),dt)
       end if
    enddo

    call restrict_base(rho0_new,.true.)
    call fill_ghost_base(rho0_new,.true.)

    call destroy(bpt)
       
  end subroutine correct_base

  subroutine correct_base_state_planar(n,rho0_new,div_etarho,dt)

    use geometry, only: anelastic_cutoff_coord, r_start_coord, r_end_coord, numdisjointchunks

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(inout) :: rho0_new(0:)
    real(kind=dp_t), intent(in   ) :: div_etarho(0:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer :: r,i
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE RHO0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i=1,numdisjointchunks(n)
       do r=r_start_coord(n,i),r_end_coord(n,i)
          if (r .lt. anelastic_cutoff_coord(n)) then
             rho0_new(r) = rho0_new(r) - dt*div_etarho(r)
          end if
       end do
    end do
    
  end subroutine correct_base_state_planar

  subroutine correct_base_state_spherical(n,rho0_new,div_etarho,dt)

    use geometry, only: anelastic_cutoff_coord, r_cc_loc, r_edge_loc, dr

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(inout) :: rho0_new(0:)
    real(kind=dp_t), intent(in   ) :: div_etarho(0:)
    real(kind=dp_t), intent(in   ) :: dt
    
    ! Local variables
    integer :: r
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE RHO0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do r=0,anelastic_cutoff_coord(n)-1
       rho0_new(r) = rho0_new(r) - dt*div_etarho(r)
    end do
    
  end subroutine correct_base_state_spherical
  
end module correct_base_module
