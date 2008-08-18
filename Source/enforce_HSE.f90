module enforce_HSE_module

  use bl_types

  implicit none

  private

  public :: enforce_HSE

contains

  subroutine enforce_HSE(rho0,p0,grav)

    use geometry, only: dr, r_start_coord, r_end_coord, numdisjointchunks, spherical

    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    real(kind=dp_t), intent(inout) ::   p0(:,0:)
    real(kind=dp_t), intent(in   ) :: grav(:,0:)

    if (spherical .eq. 0) then





    else

    end if

  end subroutine enforce_HSE
  
end module enforce_HSE_module
