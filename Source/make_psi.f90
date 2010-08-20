! Create the psi term, where psi = D_0 p_0/Dt 

module make_psi_module

  use bl_types
  implicit none
  private
  public :: make_psi_planar, make_psi_spherical

contains

  subroutine make_psi_planar(etarho_cc,psi)

    use bl_constants_module
    use geometry, only: base_cutoff_density_coord, r_start_coord, r_end_coord, &
         numdisjointchunks, nlevs_radial
    use probin_module, only: grav_const
    use restrict_base_module

    real(kind=dp_t), intent(in   ) :: etarho_cc(:,0:)
    real(kind=dp_t), intent(inout) ::       psi(:,0:)
    
    ! Local variables
    integer         :: r,i,n
   
    psi = ZERO

    do n=1,nlevs_radial
       do i=1,numdisjointchunks(n)
          !$OMP PARALLEL DO PRIVATE(r)
          do r = r_start_coord(n,i), r_end_coord(n,i)
             if (r .lt. base_cutoff_density_coord(n)) then
                psi(n,r) = etarho_cc(n,r) * abs(grav_const)
             end if
          end do
          !$OMP END PARALLEL DO
       end do
    end do

    call restrict_base(psi,.true.)
    call fill_ghost_base(psi,.true.)
    
  end subroutine make_psi_planar

  subroutine make_psi_spherical(psi,w0,gamma1bar,p0_avg,Sbar_in)

    use bl_constants_module
    use geometry, only: dr, r_cc_loc, r_edge_loc, nr_fine, base_cutoff_density_coord

    real(kind=dp_t), intent(inout) ::       psi(:,0:)
    real(kind=dp_t), intent(in   ) ::        w0(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) ::    p0_avg(:,0:)
    real(kind=dp_t), intent(in   ) ::   Sbar_in(:,0:)
    
    ! local variables
    integer :: r
    real(kind=dp_t) :: div_w0_sph

    psi = ZERO

    !$OMP PARALLEL DO PRIVATE(r,div_w0_sph)
    do r=0,base_cutoff_density_coord(1)-1

       div_w0_sph = one/(r_cc_loc(1,r)**2)* &
            (r_edge_loc(1,r+1)**2 * w0(1,r+1) - &
             r_edge_loc(1,r  )**2 * w0(1,r  )) / dr(1)

       psi(1,r) = gamma1bar(1,r) * p0_avg(1,r) * (Sbar_in(1,r) - div_w0_sph)

    enddo
    !$OMP END PARALLEL DO

  end subroutine make_psi_spherical
  
end module make_psi_module
