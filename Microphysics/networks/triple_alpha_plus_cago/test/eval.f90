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

  integer :: n

  integer :: ihe4, ic12, io16, ife56
  common /species_index/ ihe4, ic12, io16, ife56

  integer :: ir3a, irg3a, ircago, irogac
  common /reaction_index/ ir3a, irg3a, ircago, irogac

  common /zone_state/ dens, c_p, dhdx, T_eos, dT_crit

  call network_init()
  call eos_init()

  dens = 2.e9_dp_t
  temp = 7.e8_dp_t

  ihe4  = network_species_index("helium-4")
  ic12  = network_species_index("carbon-12")
  io16  = network_species_index("oxygen-16")
  ife56 = network_species_index("iron-56")
  
  ir3a   = network_reaction_index("3agc")
  irg3a  = network_reaction_index("cg3a")
  ircago = network_reaction_index("cago")
  irogac = network_reaction_index("ogac")

  Xin(ihe4)  = 0.5_dp_t
  Xin(ic12)  = 0.25_dp_t
  Xin(io16)  = 0.25_dp_t
  Xin(ife56) = 0.0_dp_t


  den_eos = dens
  temp_eos = temp
  xn_eos(:) = Xin(:)

  call eos(eos_input_rt, den_eos, temp_eos, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           .false.)


  print *, 'evaluating the RHS...'

  ! load the state
  y(1) = Xin(ihe4)
  y(2) = Xin(ic12)
  y(3) = Xin(io16)
  y(4) = Xin(ife56)
  y(5) = temp

  ! load the common block
  T_eos = temp
  dT_crit = 0.01d0
  dhdx(:) = dhdX_eos(:)
  c_p = cp_eos
  

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
