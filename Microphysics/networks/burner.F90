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


  subroutine burner(state_in, state_out, dt)

    !$acc routine seq
    
    use bl_types, only: dp_t
#ifndef SDC
    use burn_type_module, only: burn_t
#else
    use sdc_type_module, only: sdc_t
#endif
    use actual_burner_module, only: actual_burner

    implicit none

#ifndef SDC
    type (burn_t),   intent(in   ) :: state_in
    type (burn_t),   intent(inout) :: state_out
#else
    type (sdc_t),   intent(in   )  :: state_in
    type (sdc_t),   intent(inout)  :: state_out
#endif
    real(kind=dp_t), intent(in   ) :: dt
    real(kind=dp_t) :: time

    time = 0.0d0

    call actual_burner(state_in, state_out, dt, time)

  end subroutine burner

end module burner_module
