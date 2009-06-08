program testburn

  use bl_types
  use bl_constants_module
  use network
  use eos_module
  use burner_module

  implicit none

  real(kind=dp_t) :: dens, temp
  real(kind=dp_t), dimension(nspec) :: Xin
  real(kind=dp_t), dimension(nspec+1) :: y, ydot
  real(kind=dp_t) :: enucdot

  real(kind=dp_t) :: rpar
  integer :: ipar

  real(kind=dp_t) :: c_p, dhdx(nspec), T_eos, dT_crit

  integer :: ic12, io16, img24
  integer :: n

  common /zone_state/ dens, c_p, dhdx, T_eos, dT_crit

  call network_init()
  call eos_init()

  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")
  img24 = network_species_index("magnesium-24")

  dens = 2.e9_dp_t
  temp = 7.e8_dp_t

  Xin(ic12) = 0.5_dp_t
  Xin(io16) = 0.5_dp_t
  Xin(img24) = 0.0_dp_t


  den_eos(1) = dens
  temp_eos(1) = temp
  xn_eos(1,:) = Xin(:)
  
  call eos(eos_input_rt, den_eos, temp_eos, &
           npts, nspec, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           do_diag)


  print *, 'evaluating the RHS...'

  ! load the state
  y(1) = Xin(ic12)
  y(2) = Xin(io16)
  y(3) = Xin(img24)
  y(4) = temp

  ! load the common block
  T_eos = temp
  dT_crit = 0.01d0
  dhdx(:) = dhdX_eos(1,:)
  c_p = cp_eos(1)
  

  call f_rhs(nspec+1, ZERO, y, ydot, rpar, ipar)

  
  print *, 'done!'

  print *, 'Xin:  ', Xin
  print *, 'rhs:  ', ydot

  ! compute the energy release/s (erg/g/s)
  enucdot = ZERO
  do n = 1, nspec
     enucdot = enucdot - ebin(n)*ydot(n)
  enddo
  print *, 'enucdot = ', enucdot

end program testburn
