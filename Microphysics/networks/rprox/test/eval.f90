! a simple code to check the analytic Jacobian via numerical
! differencing

program eval

  use bl_types
  use bl_constants_module
  use bl_error_module
  use network
  use eos_module
  use eos_type_module
  use burner_module
  use rpar_indices

  implicit none

  real(kind=dp_t) :: dens, temp
  integer, parameter :: neq = nspec+1
  real(kind=dp_t), dimension(nspec) :: Xin
  real(kind=dp_t), dimension(neq) :: y, ydot
  real(kind=dp_t), dimension(neq) :: yp, ym
  real(kind=dp_t), dimension(neq) :: ydotp, ydotm
  real(kind=dp_t), dimension(neq,neq) :: pd
  real(kind=dp_t) :: enucdot

  real(kind=dp_t), allocatable :: rpar(:)
  integer :: ipar

  type(eos_t) :: eos_state

  integer :: i, j, n

  call network_init()
  call eos_init()

  allocate(rpar(n_rpar_comps))

  dens = 6.558e5_dp_t
  temp = 7.108e8_dp_t

  Xin(:) = ZERO
  Xin(io14) = 1.027e-3_dp_t
  Xin(io15) = 3.558e-3_dp_t
  Xin(if17) = 2.788e-2_dp_t
  Xin(is30) = 1.5735e-2_dp_t
  Xin(ihe4) = 0.2624e0_dp_t
  Xin(ih1) = 0.6894e0_dp_t

  rpar(irp_dens) = dens
  rpar(irp_T9_eos) = temp/1.e9

  eos_state%rho = dens
  eos_state%T   = temp
  eos_state%xn  = Xin

  call eos(eos_input_rt, eos_state)

  rpar(irp_cp) = eos_state%cp
  rpar(irp_dhdX:irp_dhdX+nspec-1) = eos_state%dhdX

  print *, 'evaluating the RHS...'

  ! load the state
  y(1:nspec) = Xin / aion
  y(neq) = temp/1.e9

  ! load the common block

  call f_rhs(neq, ZERO, y, ydot, rpar, ipar)

  do j = 1, nspec
     print *, Xin(j), ydot(j)
  enddo

end program eval
