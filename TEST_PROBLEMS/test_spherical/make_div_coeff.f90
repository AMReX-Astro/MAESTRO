! compute beta_0, the coefficient in our constraint equation,
! div{beta_0 U} = beta_0 S

module make_div_coeff_module

  use bl_types

  implicit none

  private

  public :: make_div_coeff

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_div_coeff(div_coeff,rho0,p0,gamma1bar,grav_center)

    use bl_constants_module
    use geometry, only: nr_fine, dr, anelastic_cutoff_coord, r_start_coord, r_end_coord, &
         nr, numdisjointchunks, nlevs_radial
    use restrict_base_module
    use probin_module, only: smallscale_beta

    real(kind=dp_t), intent(  out) :: div_coeff(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0(:,0:), p0(:,0:), gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: grav_center(:,0:)

    ! local
    integer :: r, n, i, refrat, j
    real(kind=dp_t) :: integral
    real(kind=dp_t) :: beta0_edge(nlevs_radial,0:nr_fine)
    real(kind=dp_t) :: lambda, mu, nu
    real(kind=dp_t) :: denom, coeff1, coeff2
    real(kind=dp_t) :: del,dpls,dmin,slim,sflag
    real(kind=dp_t) :: offset

    if (smallscale_beta) then

       div_coeff = 1.d0

    else

       ! for this test problem, we want beta_0 = rho_0
       
       do n=1,nlevs_radial

          do j=1,numdisjointchunks(n)

             do r=r_start_coord(n,j),r_end_coord(n,j)

                div_coeff(i,r) = rho0(i,r)

             enddo

          end do ! end loop over disjoint chunks

       end do ! end loop over levels

       call fill_ghost_base(div_coeff,.true.)
       call restrict_base(div_coeff,.true.)

    end if

  end subroutine make_div_coeff

end module make_div_coeff_module
