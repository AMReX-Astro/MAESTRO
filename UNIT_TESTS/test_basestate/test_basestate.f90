module test_basestate_module

  use bl_types
  implicit none

  public get_heating

contains

  subroutine get_heating(Hbar)
  
    use geometry, ONLY : nr, spherical, r_cc_loc

    implicit none

    real(dp_t), intent(inout) :: Hbar(0:)

    real(dp_t) :: y_0
    integer :: r
    
    do r = 0, nr(1)-1
       if (spherical .eq. 0) then
          ! plane-parallel -- do the heating term in paper II (section 4)
          y_0 = 4.d7
          Hbar(r) = 1.d17 * exp(-((r_cc_loc(1,r) - y_0)**2)/ 1.d14)
       else
          ! spherical -- lower amplitude heating term
          y_0 = 4.d7
          Hbar(r) = 1.d16 * exp(-((r_cc_loc(1,r) - y_0)**2)/ 1.d14)
       endif
       
    enddo

    return
  end subroutine get_heating

end module test_basestate_module

