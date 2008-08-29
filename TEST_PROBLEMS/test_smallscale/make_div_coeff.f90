module make_div_coeff_module

  use bl_types

  implicit none

  private

  public :: make_div_coeff

contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_div_coeff(div_coeff,rho0,p0,gamma1bar,grav_center)

    use geometry, only: r_start_coord, r_end_coord
    use probin_module, only: nlevs

    real(kind=dp_t), intent(  out) :: div_coeff(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0(:,0:), p0(:,0:), gamma1bar(:,0:), grav_center(:,0:)

    div_coeff = 1.d0

  end subroutine make_div_coeff

end module make_div_coeff_module

