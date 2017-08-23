module actual_conductivity_module

  use bl_types
  implicit none

contains

  subroutine actual_conductivity_init()
    implicit none
  end subroutine actual_conductivity_init

  subroutine actual_conductivity(eos_state, conductivity)

    use eos_type_module
    use extern_probin_module, only: conductivity_constant

    type (eos_t)    , intent(inout) :: eos_state
    real (kind=dp_t), intent(inout) :: conductivity

    ! fill the conductivity
    conductivity = conductivity_constant

  end subroutine actual_conductivity

end module actual_conductivity_module

