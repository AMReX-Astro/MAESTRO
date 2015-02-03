module init_perturb_module

  use bl_constants_module
  use variables
  use network, only: nspec
  use eos_module, only: eos_input_rp, eos
  use eos_type_module
  use bl_error_module

  implicit none

  private
  public :: perturb_2d, perturb_3d, perturb_3d_sphr

contains

  subroutine perturb_2d(x, y, p0_init, s0_init, &
                        dens_pert, rhoh_pert, rhoX_pert, &
                        temp_pert, trac_pert)

    use geometry, only: center
    use probin_module, only: pert_factor, y_pert_center, pert_width, &
         prob_lo_x, prob_hi_x, single

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: xn(nspec)
    real(kind=dp_t) :: dens
    real(kind=dp_t) :: x1, y1, r1, x2, y2, r2

    type (eos_t) :: eos_state

    if (.not. single) then
       x1 = prob_lo_x + (prob_hi_x-prob_lo_x)/3.d0
       x2 = prob_lo_x + 2.d0*(prob_hi_x-prob_lo_x)/3.d0

       y1 = y_pert_center
       y2 = y_pert_center

       r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / pert_width
       r2 = sqrt( (x-x2)**2 +(y-y2)**2 ) / pert_width
    
       if (r1 < 2.0d0) then
          dens = s0_init(rho_comp) * (1.d0 - (pert_factor * (1.d0 + tanh(2.d0-r1))))
          xn(:) = 0.d0
          xn(2) = 1.d0

       else if (r2 < 2.0d0) then
          dens = s0_init(rho_comp) * (1.d0 - (pert_factor * (1.d0 + tanh(2.d0-r2))))
          xn(:) = 0.d0
          xn(3) = 1.d0

       else
          dens = s0_init(rho_comp)
          xn(:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)
       endif

    else

       x1 = prob_lo_x + 0.5d0*(prob_hi_x-prob_lo_x)
       y1 = y_pert_center
       r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / pert_width
    
       if (r1 < 2.0d0) then
          dens = s0_init(rho_comp) * (1.d0 - (pert_factor * (1.d0 + tanh(2.d0-r1))))
          xn(:) = 0.d0
          xn(2) = 1.d0

       else
          dens = s0_init(rho_comp)
          xn(:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)
       endif

    endif

    ! Use the EOS to make this density perturbation occur at constant 
    ! pressure
    eos_state%T     = 10000.d0   ! guess
    eos_state%p     = p0_init
    eos_state%rho   = dens
    eos_state%xn(:) = xn(:)

    call eos(eos_input_rp, eos_state, .false.)

    dens_pert = dens
    rhoh_pert = dens_pert * eos_state%h
    rhoX_pert = dens_pert * xn(:)

    temp_pert = eos_state % T

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
