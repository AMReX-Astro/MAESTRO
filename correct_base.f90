module correct_base_module

  use bl_types

  implicit none

  private

  public :: correct_base

contains

  subroutine correct_base(nlevs,s0_old,s0_new,etarho,dz,dt)

    ! NOTE: the convection at the moment here is that s0_new
    ! is updated.  s0_old is only used to find the anelastic
    ! cutoff.  

    use bl_prof_module
    use geometry, only: spherical

    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(in   ) :: s0_old(:,0:,:)
    real(kind=dp_t), intent(inout) :: s0_new(:,0:,:)
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
          call correct_base_state_planar(n,s0_old(n,0:,:),s0_new(n,0:,:),etarho(n,0:), &
                                         dz(n),dt)
       end if
    enddo

    call destroy(bpt)
       
  end subroutine correct_base

  subroutine correct_base_state_planar(n,s0_old,s0_new,etarho,dz,dt)

    use bl_constants_module
    use eos_module
    use variables, only: spec_comp, rho_comp, rhoh_comp
    use geometry, only: nr
    use probin_module, only: grav_const, anelastic_cutoff, enthalpy_pred_type

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(in   ) :: s0_old(0:,:)
    real(kind=dp_t), intent(inout) :: s0_new(0:,:)
    real(kind=dp_t), intent(in   ) :: etarho(0:)
    real(kind=dp_t), intent(in   ) :: dz,dt
    
    ! Local variables
    integer :: r,comp
    integer :: r_anel
   
    ! This is used to zero the etarho contribution above the anelastic_cutoff
    r_anel = nr(1)-1
    do r = 0,nr(1)-1
       if (s0_old(r,rho_comp) .lt. anelastic_cutoff .and. r_anel .eq. nr(1)-1) then
          r_anel = r
          exit
       end if
    end do
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! UPDATE RHO0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do r = 0, r_anel-1
      s0_new(r,rho_comp) = s0_new(r,rho_comp) - dt/dz*(etarho(r+1) - etarho(r))
    end do
    
  end subroutine correct_base_state_planar
  
end module correct_base_module
