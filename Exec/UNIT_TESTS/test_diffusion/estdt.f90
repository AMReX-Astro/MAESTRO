! Compute the timestep
module estdt_module

contains

  subroutine estdt(dx,diffusion_coefficient,dt)

    use bl_constants_module, only: HALF
    use bl_types

    implicit none

    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: diffusion_coefficient
    real(kind=dp_t), intent(  out) :: dt

    ! calculate the timestep
    dt = minval(dx*dx / diffusion_coefficient)

  end subroutine estdt
  
end module estdt_module
