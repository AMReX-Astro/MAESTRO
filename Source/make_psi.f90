! Create the psi term, where psi = D_0 p_0/Dt 

module make_psi_module

  use bl_types
  implicit none
  private
  public :: make_psi

contains

  subroutine make_psi(etarho_cc,psi,w0,gamma1bar,p0_old,p0_new,Sbar_in)

    use bl_prof_module
    use geometry, only: spherical
    use restrict_base_module
    use probin_module, only: nlevs

    real(kind=dp_t), intent(in   ) :: etarho_cc(:,0:)
    real(kind=dp_t), intent(inout) :: psi(:,0:)
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
          call make_psi_planar(n,etarho_cc(n,0:),psi(n,0:))
       end do
    else
       do n = 1,nlevs
          call make_psi_spherical(n,psi(n,0:),w0(n,0:),gamma1bar(n,0:), &
                                  p0_old(n,0:),p0_new(n,0:),Sbar_in(n,0:))
       end do
    endif

    call fill_ghost_base(psi,.true.)
    call restrict_base(psi,.true.)

    call destroy(bpt)
       
  end subroutine make_psi

  subroutine make_psi_planar(n,etarho_cc,psi)

    use bl_constants_module
    use geometry, only: anelastic_cutoff_coord, r_start_coord, r_end_coord, numdisjointchunks
    use probin_module, only: grav_const

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(in   ) :: etarho_cc(0:)
    real(kind=dp_t), intent(inout) :: psi(0:)
    
    ! Local variables
    integer         :: r,i
   
    psi = ZERO

    do i=1,numdisjointchunks(n)
       do r = r_start_coord(n,i), r_end_coord(n,i)
          if (r .lt. anelastic_cutoff_coord(n)) then
             psi(r) = etarho_cc(r) * abs(grav_const)
          end if
       end do
    end do
    
  end subroutine make_psi_planar

  subroutine make_psi_spherical(n,psi,w0,gamma1bar,p0_old,p0_new,Sbar_in)

    use bl_constants_module
    use geometry, only: dr, r_cc_loc, r_edge_loc, nr_fine

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(inout) :: psi(0:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:), p0_new(0:)
    real(kind=dp_t), intent(in   ) :: Sbar_in(0:)
    
    ! local variables
    integer :: r
    real(kind=dp_t) :: div_w0_sph

    do r=0,nr_fine-1

       div_w0_sph = one/(r_cc_loc(n,r)**2)* &
            (r_edge_loc(n,r+1)**2 * w0(r+1) - &
             r_edge_loc(n,r  )**2 * w0(r  )) / dr(n)

       psi(r) = gamma1bar(r) * HALF*(p0_old(r) + p0_new(r)) * (Sbar_in(r) - div_w0_sph)

    enddo

  end subroutine make_psi_spherical
  
end module make_psi_module
