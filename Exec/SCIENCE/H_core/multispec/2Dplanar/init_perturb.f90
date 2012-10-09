! apply an optional perturbation to the scalar state.  This routine 
! is called on a zone-by-zone basis from init_scalar_data.  It is 
! assumed that the perturbation is done at constant pressure.

module init_perturb_module

  use variables
  use network, only: nspec
  use eos_module
  use bl_constants_module

  implicit none

  private
  public :: perturb_2d, perturb_3d, perturb_3d_sphr

contains

  subroutine perturb_2d(x, y, p0_init, s0_init, dens_pert, rhoh_pert, rhoX_pert, &
                        temp_pert, trac_pert)

    use mt19937_module
    use probin_module, only: octant, &
         scalpert_amplitude, scalpert_radius, scalpert_steep, scalpert_scale
    use geometry, only: center

    real(kind=dp_t), intent(in ) :: x, y
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: rho, temp
    real(kind=dp_t) :: x0, y0, r0

    ! random numbers between -1 and 1
    real(kind=dp_t) :: alpha(3,3,3), beta(3,3,3), gamma(3,3,3)

    ! random numbers between 0 and 2*pi
    real(kind=dp_t) :: phix(3,3,3), phiy(3,3,3), phiz(3,3,3)

    ! L2 norm of k
    real(kind=dp_t) :: normk(3,3,3)

    ! random number
    real(kind=dp_t) :: rand
    
    ! Local variables
    integer :: i, j, k

    ! cos and sin of (2*pi*kx/L + phix), etc
    real(kind=dp_t) :: cx(3,3,3), cy(3,3,3), cz(3,3,3)
    real(kind=dp_t) :: sx(3,3,3), sy(3,3,3), sz(3,3,3)

    real(kind=dp_t) :: theta,phi


    ! load in random numbers alpha, beta, gamma, phix, phiy, and phiz
!     call init_genrand(20908)
!     do i=1,3
!        do j=1,3
!           do k=1,3
!              rand = genrand_real1()
!              rand = 2.0d0*rand - 1.0d0
!              alpha(i,j,k) = rand
!              rand = genrand_real1()
!              rand = 2.0d0*rand - 1.0d0
!              beta(i,j,k) = rand
!              rand = genrand_real1()
!              rand = 2.0d0*rand - 1.0d0
!              gamma(i,j,k) = rand
!              rand = genrand_real1()
!              rand = 2.0d0*M_PI*rand
!              phix(i,j,k) = rand
!              rand = genrand_real1()
!              rand = 2.0d0*M_PI*rand
!              phiy(i,j,k) = rand
!              rand = genrand_real1()
!              rand = 2.0d0*M_PI*rand
!              phiz(i,j,k) = rand
!           enddo
!        enddo
!     enddo

!     ! compute the norm of k
!     do i=1,3
!        do j=1,3
!           do k=1,3
!              normk(i,j,k) = sqrt(dble(i)**2+dble(j)**2+dble(k)**2)
!           enddo
!        enddo
!     enddo


! ! random temperature fluctuations
!     temp = ZERO

!     ! loop over the 27 combinations of fourier components
!     do i=1,3
!        do j=1,3
!           do k=1,3
!              ! compute cosines and sines
!              cx(i,j,k) = cos(2.0d0*M_PI*dble(i)*x/scalpert_scale + phix(i,j,k))
!              cy(i,j,k) = cos(2.0d0*M_PI*dble(j)*y/scalpert_scale + phiy(i,j,k))
!              cz(i,j,k) = cos(phiz(i,j,k))
!              sx(i,j,k) = sin(2.0d0*M_PI*dble(i)*x/scalpert_scale + phix(i,j,k))
!              sy(i,j,k) = sin(2.0d0*M_PI*dble(j)*y/scalpert_scale + phiy(i,j,k))
!              sz(i,j,k) = sin(phiz(i,j,k))
!           enddo
!        enddo
!     enddo

!     ! loop over the 27 combinations of fourier components
!     do i=1,3
!        do j=1,3
!           do k=1,3
!              ! compute contribution from perturbation scalar from each mode
!              temp = temp + &
!                     (-gamma(i,j,k)*dble(j)*cx(i,j,k)*cz(i,j,k)*sy(i,j,k) &
!                      +beta(i,j,k)*dble(k)*cx(i,j,k)*cy(i,j,k)*sz(i,j,k)) &
!                     / normk(i,j,k)
!           enddo
!        enddo
!     enddo
    
!     ! apply the cutoff function to the perturbational scalar
!     ! with 2D hack y is like radius
!     temp = scalpert_amplitude * temp &
!            *(0.5d0+0.5d0*tanh((scalpert_radius - y)/scalpert_steep))
    
!     ! add perturbational velocity to background velocity
!     temp = temp + s0_init(temp_comp)


! tanh density perturbation
    rho = s0_init(rho_comp)

!    x0 = center(1) + 5.d10
!    y0 = 7.35d9

    x0 = center(1) !+ 2.d10
    y0 = 3.5d10

    ! Tanh bubbles
    r0 = sqrt( (x-x0)**2 + (y-y0)**2 ) / scalpert_scale !1.d9!8.e8
!    r0 = (y-y0) / 1.e10
    
!    temp = sin(TWO*M_PI*y/scalpert_scale) 
!    temp = s0_init(temp_comp) + temp * scalpert_amplitude *  &
!           (0.5d0+0.5d0*tanh((scalpert_radius - y)/scalpert_steep))
!    temp = temp + 3.d3*(1.d0 + tanh(2.0_dp_t-r0))

    ! This case works for a hot bubble - used for 2D v full test
!    rho = rho - 3.d-2*(1.d0 + tanh(2.0_dp_t-r0))
!    temp = s0_init(temp_comp) + 3.d5*(1.d0 + tanh(2.0_dp_t-r0))
     temp = s0_init(temp_comp) + scalpert_amplitude *  &
           (1.d0 + tanh((2.d0 - r0)/scalpert_steep))
   

    ! use the EOS to make this temperature perturbation occur at
    ! constant pressure
    temp_eos(1) = temp
    p_eos(1) = p0_init
    den_eos(1) = rho
    xn_eos(1,:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

!    write(*,*) 'temp ', temp
!    write (*,*)'rho ', rho

    pt_index_eos(:) = (/i, j, -1/)

    call eos(eos_input_tp, den_eos, temp_eos, &
!    call eos(eos_input_rp, den_eos, temp_eos, &
             npts, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             .false., pt_index_eos)

    dens_pert = den_eos(1)
    rhoh_pert = den_eos(1)*h_eos(1)
    rhoX_pert = dens_pert*xn_eos(1,:)

    temp_pert = temp_eos(1)

    trac_pert = ZERO

  end subroutine perturb_2d

  subroutine perturb_3d(x, y, z, p0_init, s0_init, dens_pert, rhoh_pert, &
                        rhoX_pert, temp_pert, trac_pert)

    real(kind=dp_t), intent(in ) :: x, y, z
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp


    temp = s0_init(temp_comp)

    ! apply some perturbation to density here
    ! temp = ...

    ! use the EOS to make this temperature perturbation occur at constant 
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
    rhoX_pert = dens_pert*xn_eos(1,:)

    temp_pert = temp

    trac_pert = ZERO

  end subroutine perturb_3d

  subroutine perturb_3d_sphr(x, y, z, p0_init, s0_init, dens_pert, rhoh_pert, &
                             rhoX_pert, temp_pert, trac_pert)

    use geometry, only: center
    
    real(kind=dp_t), intent(in ) :: x, y, z
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp, r0, x0, y0, z0


    temp = s0_init(temp_comp)

    x0 = center(1) 
    y0 = center(2) + 1.04d10
    z0 = center(3) 

    ! Tanh bubbles
    r0 = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 ) / 1.15d10
    
    ! This case works
    ! temp = t0 * (ONE + TWO*(.150_dp_t * 0.5_dp_t * & 
!                             (1.0_dp_t + tanh((2.0_dp_t-r0)))))
    temp = temp - 3.d-6 * tanh(2.0_dp_t-r0)

    ! use the EOS to make this temperature perturbation occur at constant 
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
    rhoX_pert = dens_pert*xn_eos(1,:)

    temp_pert = temp

    trac_pert = ZERO

  end subroutine perturb_3d_sphr

end module init_perturb_module

