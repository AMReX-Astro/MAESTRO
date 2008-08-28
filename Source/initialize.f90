module initialize_module


  implicit none

  private

  public :: initialize_from_restart, initialize_with_fixed_grids, &
       initialize_with_adaptive_grids, initialize_bc

  contains

    subroutine initialize_from_restart()

    end subroutine initialize_from_restart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine initialize_with_fixed_grids()
      
    end subroutine initialize_with_fixed_grids
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine initialize_with_adaptive_grids()

    end subroutine initialize_with_adaptive_grids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine initialize_bc()

    end subroutine initialize_bc
  
end module initialize_module
