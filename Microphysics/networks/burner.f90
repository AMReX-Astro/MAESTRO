module burner_module

  implicit none

  logical :: burner_initialized = .false.

contains

  subroutine burner_init()

    use actual_burner_module, only: actual_burner_init

    implicit none

    call actual_burner_init()

    burner_initialized = .true.

  end subroutine burner_init


  subroutine burner(state_in, state_out, dt, rhowdot, rhoH)

    !$acc routine seq
    
    use bl_types, only: dp_t
    use network, only: nspec
    use burn_type_module, only: burn_t
    use actual_burner_module, only: actual_burner

    implicit none
    
    type (burn_t),   intent(in   ) :: state_in
    type (burn_t),   intent(inout) :: state_out
    real(kind=dp_t), intent(in   ) :: dt
    real(kind=dp_t), intent(  out) :: rhowdot(nspec), rhoH
    real(kind=dp_t) :: time
    integer :: n

    time = 0.0d0

    call actual_burner(state_in, state_out, dt, time)

  end subroutine burner

end module burner_module
