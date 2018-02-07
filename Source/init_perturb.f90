! apply an optional perturbation to the scalar state.  This routine 
! is called on a zone-by-zone basis from init_scalar_data.  It is 
! assumed that the perturbation is done at constant pressure.

module init_perturb_module

  use variables
  use network, only: nspec
  use eos_module, only: eos_input_tp, eos
  use eos_type_module
  use bl_constants_module

  implicit none

  private
  public :: perturb_2d, perturb_3d, perturb_3d_sphr, perturb_2d_polar

contains

  subroutine perturb_2d(x, y, p0_init, s0_init, dens_pert, rhoh_pert, rhoX_pert, &
                        temp_pert, trac_pert)

    real(kind=dp_t), intent(in ) :: x, y
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp

    type (eos_t) :: eos_state
    

    temp = s0_init(temp_comp)

    ! apply some perturbation to density here
    ! temp = ...

    ! use the EOS to make this temperature perturbation occur at
    ! constant pressure
    eos_state%T     = temp
    eos_state%p     = p0_init
    eos_state%rho   = s0_init(rho_comp)
    eos_state%xn(:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_tp, eos_state)

    dens_pert = eos_state%rho
    rhoh_pert = eos_state%rho*eos_state%h
    rhoX_pert = dens_pert*eos_state%xn(:)

    temp_pert = temp

    trac_pert = ZERO

  end subroutine perturb_2d

  subroutine perturb_2d_polar(x, y, p0_init, s0_init, dens_pert, rhoh_pert, &
                             rhoX_pert, temp_pert, trac_pert)

    real(kind=dp_t), intent(in ) :: x, y
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp

    type (eos_t) :: eos_state

    temp = s0_init(temp_comp)

    ! apply some perturbation to density here
    ! temp = ...

    ! use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    eos_state%T     = temp
    eos_state%p     = p0_init
    eos_state%rho   = s0_init(rho_comp)
    eos_state%xn(:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_tp, eos_state)

    dens_pert = eos_state%rho
    rhoh_pert = eos_state%rho*eos_state%h
    rhoX_pert = dens_pert*eos_state%xn(:)

    temp_pert = temp

    trac_pert = ZERO

  end subroutine perturb_2d_polar  
  
  subroutine perturb_3d(x, y, z, p0_init, s0_init, dens_pert, rhoh_pert, &
                        rhoX_pert, temp_pert, trac_pert)

    real(kind=dp_t), intent(in ) :: x, y, z
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp

    type (eos_t) :: eos_state

    temp = s0_init(temp_comp)

    ! apply some perturbation to density here
    ! temp = ...

    ! use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    eos_state%T     = temp
    eos_state%p     = p0_init
    eos_state%rho   = s0_init(rho_comp)
    eos_state%xn(:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_tp, eos_state)

    dens_pert = eos_state%rho
    rhoh_pert = eos_state%rho*eos_state%h
    rhoX_pert = dens_pert*eos_state%xn(:)

    temp_pert = temp

    trac_pert = ZERO

  end subroutine perturb_3d

  subroutine perturb_3d_sphr(x, y, z, p0_init, s0_init, dens_pert, rhoh_pert, &
                             rhoX_pert, temp_pert, trac_pert)

    real(kind=dp_t), intent(in ) :: x, y, z
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp

    type (eos_t) :: eos_state

    temp = s0_init(temp_comp)

    ! apply some perturbation to density here
    ! temp = ...

    ! use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    eos_state%T     = temp
    eos_state%p     = p0_init
    eos_state%rho   = s0_init(rho_comp)
    eos_state%xn(:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_tp, eos_state)

    dens_pert = eos_state%rho
    rhoh_pert = eos_state%rho*eos_state%h
    rhoX_pert = dens_pert*eos_state%xn(:)

    temp_pert = temp

    trac_pert = ZERO

  end subroutine perturb_3d_sphr

end module init_perturb_module

