module heating_module

  use bl_types

  implicit none

  private
  public :: get_H_2d, get_H_3d

contains

   subroutine get_H_2d (H,lo,hi,dx,time)

     use bl_constants_module

      integer, intent(in) :: lo(:), hi(:)
      real(kind=dp_t), intent(inout) :: H(lo(1):,lo(2):)
      real(kind=dp_t), intent(in   ) :: dx(:),time

      H = 0.0_dp_t

   end subroutine get_H_2d

   subroutine get_H_3d (H,lo,hi,dx,time)

     use bl_constants_module

      integer, intent(in) :: lo(:), hi(:)
      real(kind=dp_t), intent(inout) :: H(lo(1):,lo(2):,lo(3):)
      real(kind=dp_t), intent(in   ) :: dx(:),time

      H = 0.0_dp_t

   end subroutine get_H_3d
end module heating_module
