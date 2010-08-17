module init_perturb_module

  use variables
  use network, only: nspec
  use eos_module
  use bl_constants_module

  implicit none

  private
  public :: perturb_2d, perturb_3d

contains

  subroutine perturb_2d(x, y, p0_init, s0_init, dens_pert, rhoh_pert, rhoX_pert, &
                        temp_pert, trac_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp,t0
    real(kind=dp_t) :: x0, y0, x1, y1, x2, y2
    real(kind=dp_t) :: r0, r1, r2

    t0 = s0_init(temp_comp)

    x0 = 5.0d7
    y0 = 6.5d7
    
    x1 = 1.2d8
    y1 = 8.5d7
    
    x2 = 2.0d8
    y2 = 7.5d7

    ! Tanh bubbles
    r0 = sqrt( (x-x0)**2 +(y-y0)**2 ) / 2.5e6
    r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / 2.5e6
    r2 = sqrt( (x-x2)**2 +(y-y2)**2 ) / 2.5e6
    
    ! This case works
    temp = t0 * (ONE + TWO * ( &
         .15_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r0))) + &
         .3_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r1))) + &
         .225_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r2)))  ) )

    ! This case breaks
!   temp = t0 * (ONE + FOUR * ( &
!        .15_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r0))) + &
!        .3_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r1))) + &
!        .225_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r2)))  ) )
          
    ! Use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    temp_eos(1) = temp
    p_eos(1) = p0_init
    den_eos(1) = s0_init(rho_comp)
    xn_eos(1,:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_tp, den_eos, temp_eos, &
             npts, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             .false.)

    dens_pert = den_eos(1)
    rhoh_pert = den_eos(1)*h_eos(1)
    rhoX_pert(:) = dens_pert*xn_eos(1,:)

    temp_pert = temp
    
!   if ( (r0 .lt. 2.0) .or. (r1 .lt. 2.0) .or. (r2 .lt. 2.0) ) then
!     trac_pert(:) = ONE
!   else
      trac_pert(:) = ZERO
!   end if

  end subroutine perturb_2d

  subroutine perturb_3d(x, y, z, p0_init, s0_init, dens_pert, rhoh_pert, rhoX_pert, &
                        temp_pert, trac_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y, z
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp, t0
    real(kind=dp_t) :: x0, y0, z0, x1, y1, z1, x2, y2, z2
    real(kind=dp_t) :: r0, r1, r2

    t0 = s0_init(temp_comp)

    x0 = 5.0d7
    y0 = 5.0d7
    z0 = 6.5d7
    
    x1 = 1.2d8
    y1 = 1.2d8
    z1 = 8.5d7
    
    x2 = 2.0d8
    y2 = 2.0d8
    z2 = 7.5d7

!   temp = t0 * (ONE + TWO * ( &
!        .0625_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r0))) + &
!        .1875_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r1))) + &
!        .1250_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r2)))  ) )

    ! Tanh bubbles from perturb_2d
    r0 = sqrt( (y-y0)**2 +(z-z0)**2 ) / 2.5e6
    r1 = sqrt( (y-y1)**2 +(z-z1)**2 ) / 2.5e6
    r2 = sqrt( (y-y2)**2 +(z-z2)**2 ) / 2.5e6
    
    ! This case works
    temp = t0 * (ONE + TWO * ( &
         .150_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r0))) + &
         .300_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r1))) + &
         .225_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r2)))  ) )

    ! Use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    temp_eos(1) = temp
    p_eos(1) = p0_init
    den_eos(1) = s0_init(rho_comp)
    xn_eos(1,:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_tp, den_eos, temp_eos, &
             npts, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             .false.)

    dens_pert = den_eos(1)
    rhoh_pert = den_eos(1)*h_eos(1)
    rhoX_pert(:) = dens_pert*xn_eos(1,:)

    temp_pert = temp
    
!   if (r1 .lt. 2.0) then
!     trac_pert(:) = ONE
!   else
      trac_pert(:) = ZERO
!   end if

  end subroutine perturb_3d

end module init_perturb_module
