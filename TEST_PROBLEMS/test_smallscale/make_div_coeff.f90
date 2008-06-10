module make_div_coeff_module

  use bl_types

  implicit none

  private

  public :: make_div_coeff

contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine make_div_coeff(n,div_coeff,rho0,p0,gamma1bar,grav_center)

    use geometry, only: r_start_coord, r_end_coord

    integer        , intent(in   ) :: n
    real(kind=dp_t), intent(  out) :: div_coeff(0:)
    real(kind=dp_t), intent(in   ) :: rho0(0:), p0(0:), gamma1bar(0:)
    real(kind=dp_t), intent(in   ) :: grav_center(0:)

    integer :: r

    do r=r_start_coord(n),r_end_coord(n)
       div_coeff(r) = 1.d0
    end do

  end subroutine make_div_coeff

end module make_div_coeff_module

