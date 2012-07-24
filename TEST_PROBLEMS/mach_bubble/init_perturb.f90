module init_perturb_module

  use bl_constants_module
  use variables
  use network, only: nspec
  use eos_module, only: eos_input_tp, eos
  use eos_type_module
  use bl_error_module

  implicit none

  private
  public :: perturb_2d, perturb_3d, perturb_3d_sphr

contains

  subroutine perturb_2d(x, y, p0_init, s0_init, dens_pert, rhoh_pert, rhoX_pert, &
                        temp_pert, trac_pert)

    use geometry, only: center
    use probin_module, only: pert_factor, y_pert_center, pert_width

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp,t0
    real(kind=dp_t) :: x1, y1, r1, x2, y2, r2, x3, y3, r3

    type (eos_t) :: eos_state

    t0 = s0_init(temp_comp)

    x1 = center(1)
    y1 = y_pert_center
    r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / pert_width
    
    ! temperature perturbation -- pert_factor sets the amplitude
    ! of the perturbation (higher values means larger perturbation)
    temp = t0 * (1.d0 + (pert_factor * (1.d0 + tanh(2.d0-r1))))
          
    ! Use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    eos_state%T     = temp
    eos_state%p     = p0_init
    eos_state%rho   = s0_init(rho_comp)
    eos_state%xn(:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_tp, eos_state, .false.)

    dens_pert = eos_state%rho
    rhoh_pert = eos_state%rho * eos_state%h
    rhoX_pert = dens_pert*eos_state%xn(:)

    temp_pert = temp

    trac_pert = ZERO

  end subroutine perturb_2d

  subroutine perturb_3d(x, y, z, p0_init, s0_init, dens_pert, rhoh_pert, &
                        rhoX_pert, temp_pert, trac_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y, z
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    call bl_error("ERROR: perturb_3d not implemented")

  end subroutine perturb_3d

  subroutine perturb_3d_sphr(x, y, z, p0_init, s0_init, dens_pert, rhoh_pert, &
                             rhoX_pert, temp_pert, trac_pert)

    real(kind=dp_t), intent(in ) :: x, y, z
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    call bl_error("ERROR: perturb_3d_sphr not implemented")

  end subroutine perturb_3d_sphr

end module init_perturb_module
