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



  subroutine burner(rho, T_in, x_in, dt, x_out, rhowdot, rhoH)

     !$acc routine seq

     use network, only: nspec
     use burn_type_module, only: burn_t
     use actual_burner_module, only: actual_burner

     implicit none

     real(kind=dp_t), intent(in   ) :: rho, T_in, x_in(nspec), dt
     real(kind=dp_t), intent(inout) :: x_out(nspec), rhowdot(nspec), rhoH

     real(kind=dp_t) :: time
     integer :: n
     type (burn_t) :: state_in, state_out

     time = 0.0d0

     state_in % rho = rho
     state_in % T   = T_in
     do n = 1, nspec
        state_in % xn(n) = x_in(n)
     enddo

     state_in % e = 0.0d0

     ! Initialize the outgoing state to be equal to the incoming state.

     state_out = state_in

     call actual_burner(state_in, state_out, dt, time)

     do n = 1, nspec
        x_out(n) = state_out % xn(n)
        rhowdot(n) = rho * (x_out(n) - x_in(n)) / dt
     enddo

     rhoH = rho * (state_out % e - state_in % e) / dt

   end subroutine burner

end module burner_module
