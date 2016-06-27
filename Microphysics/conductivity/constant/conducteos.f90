module conductivity_module

  use bl_types
  implicit none

contains

  subroutine conductivity_init()
  end subroutine conductivity_init

  subroutine conducteos(input, eos_state, do_diag, conductivity)

    use eos_module
    use eos_type_module
    use network, only: nspec
    use extern_probin_module, only: conductivity_constant

    integer         , intent(in   ) :: input
    type (eos_t)    , intent(inout) :: eos_state
    logical         , intent(in   ) :: do_diag
    real (kind=dp_t), intent(inout) :: conductivity

    call eos(input, eos_state)

    ! fill the conductivity
    conductivity = conductivity_constant

  end subroutine conducteos

end module conductivity_module

