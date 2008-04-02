! Create the psi term, where psi = D_0 p_0/Dt 

module make_psi_module

  use bl_types
  implicit none
  private
  public :: make_psi

contains

  subroutine make_psi(nlevs,etarho,psi,s0,w0,gamma1bar,p0_old,p0_new,Sbar_in)

    use bl_prof_module
    use geometry, only: spherical

    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(in   ) :: etarho(:,0:)
    real(kind=dp_t), intent(inout) :: psi(:,0:)
    real(kind=dp_t), intent(in   ) :: s0(:,0:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:), p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(:,0:)

    
    ! local
    integer :: n

    type(bl_prof_timer), save :: bpt
    call build(bpt, "make_psi")
    
    if (spherical .eq. 0) then 
       do n = 1,nlevs
          call make_psi_planar(n,etarho(n,0:),psi(n,0:),s0(n,0:,:))
       end do
    else
       do n = 1,nlevs
          call make_psi_spherical(n,psi(n,0:),s0(n,0:,:),w0(n,0:),gamma1bar(n,0:), &
                                  p0_old(n,0:),p0_new(n,0:),Sbar_in(n,0:))
       end do
    endif

    
    call destroy(bpt)
       
  end subroutine make_psi

  subroutine make_psi_planar(n,etarho,psi,s0)

    use bl_constants_module
    use variables, only: rho_comp
    use geometry, only: nr
    use probin_module, only: grav_const, anelastic_cutoff

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(in   ) :: etarho(0:)
    real(kind=dp_t), intent(inout) :: psi(0:)
    real(kind=dp_t), intent(in   ) :: s0(0:,:)
    
    ! Local variables
    integer         :: r,r_anel
    real(kind=dp_t) :: etarho_avg
   
    ! This is used to zero the etarho contribution above the anelastic_cutoff
    r_anel = nr(n)-1
    do r = 0,nr(n)-1
       if (s0(r,rho_comp) .lt. anelastic_cutoff .and. r_anel .eq. nr(n)-1) then
          r_anel = r
          exit
       end if
    end do

    psi(:) = ZERO

    do r = 0, r_anel-1
      etarho_avg = HALF * (etarho(r)+etarho(r+1))
      psi(r) = etarho_avg * abs(grav_const)
    end do
    
  end subroutine make_psi_planar

  subroutine make_psi_spherical(n,psi,s0,w0,gamma1bar,p0_old,p0_new,Sbar_in)

    use bl_constants_module
    use variables, only: rho_comp
    use geometry, only: nr, dr, base_cc_loc, base_loedge_loc
    use probin_module, only: grav_const, anelastic_cutoff

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(inout) :: psi(0:)
    real(kind=dp_t), intent(in   ) :: s0(0:,:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:), p0_new(0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(0:)
    
    ! local variables
    integer :: r
    real(kind=dp_t) :: div_w0_sph

    do r = 0, nr(n)-1

       div_w0_sph = one/(base_cc_loc(n,r)**2)* &
            (base_loedge_loc(n,r+1)**2 * w0(r+1) - &
             base_loedge_loc(n,r  )**2 * w0(r  )) / dr(n)

       psi(r) = gamma1bar(r) * HALF*(p0_old(r) + p0_new(r)) * (Sbar_in(r) - div_w0_sph)

    enddo


  end subroutine make_psi_spherical

  
end module make_psi_module
