module conductivity_module

  use bl_types
  implicit none

  real (kind=dp_t), save, private :: conductivity_constant

contains

  subroutine conductivity_init(cond_const)

    real (kind=dp_t), optional :: cond_const

    if (present(cond_const)) then
       conductivity_constant = cond_const
    else
       conductivity_constant = 1.0_dp_t
    endif

  end subroutine conductivity_init

  subroutine conducteos(input, eos_state, do_diag, conductivity)

    use eos_module
    use eos_type_module
    use network, only: nspec

    integer         , intent(in   ) :: input
    type (eos_t)    , intent(inout) :: eos_state
    logical         , intent(in   ) :: do_diag
    real (kind=dp_t), intent(inout) :: conductivity

    call eos(input, eos_state, do_diag)

    ! fill the conductivity
    conductivity = conductivity_constant

  end subroutine conducteos

end module conductivity_module

