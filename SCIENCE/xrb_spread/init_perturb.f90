module init_perturb_module

  use variables
  use network, only: nspec
  use eos_module, only: eos_input_tp, eos
  use eos_type_module
  use bl_constants_module

  implicit none

  private
  public :: perturb_2d, perturb_3d, perturb_3d_sphr

contains

  subroutine perturb_2d(x, y, p0_init, s0_init, &
                        dens_pert, rhoh_pert, rhoX_pert, &
                        temp_pert, trac_pert)

    use probin_module, only: prob_lo, prob_hi, &
                             xrb_pert_temp, xrb_pert_size, xrb_pert_height

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)
    
    real(kind=dp_t) :: temp,t0
    real(kind=dp_t) :: xmid, dist

    type (eos_t) :: eos_state

    xmid = 0.5_dp_t * (prob_lo(1) + prob_hi(1))

    dist = sqrt((x - xmid)**2 + (y - xrb_pert_height)**2)

    

    if (dist < xrb_pert_size) then
       temp = xrb_pert_temp

    else
       temp = s0_init(temp_comp)

    endif

          
    ! Use the EOS to make this temperature perturbation occur at
    ! constant pressure
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

    real(kind=dp_t) :: temp, t0
    real(kind=dp_t) :: x0, y0, z0, r0

    

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

    call eos(eos_input_tp, eos_state, .false.)

    dens_pert = eos_state%rho
    rhoh_pert = eos_state%rho * eos_state%h
    rhoX_pert = dens_pert*eos_state%xn(:)

    temp_pert = temp

    trac_pert = ZERO

  end subroutine perturb_3d_sphr

end module init_perturb_module

