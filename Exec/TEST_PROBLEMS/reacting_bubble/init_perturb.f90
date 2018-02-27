module init_perturb_module

  use variables
  use network, only: nspec
  use eos_module, only: eos, eos_input_tp
  use eos_type_module
  use bl_constants_module

  implicit none

  private
  public :: perturb_2d, perturb_3d, perturb_3d_sphr

contains

  subroutine perturb_2d(x, y, p0_init, s0_init, dens_pert, rhoh_pert, rhoX_pert, &
                        temp_pert, trac_pert)

    use probin_module, only: pert_temp_factor, pert_rad_factor, do_small_domain

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp,t0
    real(kind=dp_t) :: x1, y1, r1, x2, y2, r2, x3, y3, r3, x4, y4, r4

    type (eos_t) :: eos_state

    t0 = s0_init(temp_comp)

    x1 = 5.0d7
    y1 = 6.5d7
    r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / (2.5d6*pert_rad_factor)
    
    x2 = 1.2d8
    y2 = 8.5d7
    r2 = sqrt( (x-x2)**2 +(y-y2)**2 ) / (2.5d6*pert_rad_factor)
    
    x3 = 2.0d8
    y3 = 7.5d7
    r3 = sqrt( (x-x3)**2 +(y-y3)**2 ) / (2.5d6*pert_rad_factor)

    ! this is a tiny bubble for inputs_2d_smalldomain
    x4 = 0.5d0
    y4 = 85218750.25d0
    r4 = sqrt( (x-x4)**2 +(y-y4)**2 ) / (2.5d-2*pert_rad_factor)

    if (do_small_domain) then
       temp = t0 * (1.d0 + pert_temp_factor* &
            (0.150d0 * (1.d0 + tanh(2.d0-r1)) + &
            0.300d0 * (1.d0 + tanh(2.d0-r2)) + &
            0.225d0 * (1.d0 + tanh(2.d0-r3)) + &
            0.300d0 * (1.d0 + tanh(2.d0-r4))))
    else
       temp = t0 * (1.d0 + pert_temp_factor* &
            (0.150d0 * (1.d0 + tanh(2.d0-r1)) + &
            0.300d0 * (1.d0 + tanh(2.d0-r2)) + &
            0.225d0 * (1.d0 + tanh(2.d0-r3))))
    end if
          
    ! Use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    eos_state%T     = temp
    eos_state%p     = p0_init
    eos_state%rho   = s0_init(rho_comp)
    eos_state%xn(:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_tp, eos_state)

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

    real(kind=dp_t) :: temp, t0
    real(kind=dp_t) :: x0, y0, z0, r0

    type (eos_t) :: eos_state

    t0 = s0_init(temp_comp)

    x0 = 1.35d7
    y0 = 1.35d7
    z0 = 8.5d7

    r0 = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 ) / 2.5d6

    temp = t0 * (ONE + TWO * (0.15d0 * (1.d0 + tanh((2.d0-r0)))))

    ! Use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    eos_state%T     = temp
    eos_state%p     = p0_init
    eos_state%rho   = s0_init(rho_comp)
    eos_state%xn(:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_tp, eos_state)

    dens_pert = eos_state%rho
    rhoh_pert = eos_state%rho * eos_state%h
    rhoX_pert = dens_pert*eos_state%xn(:)

    temp_pert = temp

    trac_pert = ZERO

  end subroutine perturb_3d

  subroutine perturb_3d_sphr(x, y, z, p0_init, s0_init, dens_pert, rhoh_pert, &
                             rhoX_pert, temp_pert, trac_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles
    real(kind=dp_t), intent(in ) :: x, y, z
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp, t0
    real(kind=dp_t) :: x0, y0, z0, r0

    type (eos_t) :: eos_state

    t0 = s0_init(temp_comp)

    ! center of star is at 2.5d8
    x0 = 2.5d8
    y0 = 2.5d8
    z0 = 3.0d8

    ! Tanh bubbles
    r0 = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 ) / 2.5d6
    
    ! This case works
    temp = t0 * (ONE + TWO*(.150_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0_dp_t-r0)))))

    ! Use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    eos_state%T     = temp
    eos_state%p     = p0_init
    eos_state%rho   = s0_init(rho_comp)
    eos_state%xn(:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_tp, eos_state)

    dens_pert = eos_state%rho
    rhoh_pert = eos_state%rho * eos_state%h
    rhoX_pert = dens_pert*eos_state%xn(:)

    temp_pert = temp

    trac_pert = ZERO

  end subroutine perturb_3d_sphr

end module init_perturb_module
