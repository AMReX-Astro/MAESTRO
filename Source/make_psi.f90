module make_psi_module

  use bl_types
  implicit none
  private
  public :: make_psi

contains

  subroutine make_psi(eta,psi,s0_old)

    use bl_prof_module
    use geometry, only: spherical

    real(kind=dp_t), intent(in   ) :: eta(:,0:,:)
    real(kind=dp_t), intent(inout) :: psi(:,0:)
    real(kind=dp_t), intent(in   ) :: s0_old(:,0:,:)
    
    ! local
    integer :: n,nlevs

    type(bl_prof_timer), save :: bpt
    call build(bpt, "make_psi")

    nlevs = size(psi,dim=1)
    
    if (spherical .eq. 0) then
       do n = 1,nlevs
          call make_psi_planar(n,eta(n,0:,:),psi(n,0:),s0_old(n,0:,:))
       end do
    end if

    call destroy(bpt)
       
  end subroutine make_psi

  subroutine make_psi_planar(n,eta,psi,s0_old)

    use bl_constants_module
    use variables, only: rho_comp
    use geometry, only: nr
    use probin_module, only: grav_const, anelastic_cutoff

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(in   ) :: eta(0:,:)
    real(kind=dp_t), intent(inout) :: psi(0:)
    real(kind=dp_t), intent(in   ) :: s0_old(0:,:)
    
    ! Local variables
    integer         :: r,r_anel
    real(kind=dp_t) :: eta_avg
   
    ! This is used to zero the eta contribution above the anelastic_cutoff
    r_anel = nr(n)-1
    do r = 0,nr(n)-1
       if (s0_old(r,rho_comp) .lt. anelastic_cutoff .and. r_anel .eq. nr(n)-1) then
          r_anel = r
          exit
       end if
    end do

    psi(:) = ZERO

    do r = 0, r_anel-1
      eta_avg = HALF * (eta(r,rho_comp)+eta(r+1,rho_comp))
      psi(r) = eta_avg * abs(grav_const)
    end do
    
  end subroutine make_psi_planar
  
end module make_psi_module
