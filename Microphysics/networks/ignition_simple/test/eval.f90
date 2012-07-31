program testburn

  use bl_types
  use bl_constants_module
  use bl_error_module
  use network
  use eos_module, only: eos_input_rt, eos, eos_init
  use eos_type_module
  use burner_module
  use burner_aux_module

  implicit none

  real(kind=dp_t) :: dens, temp
  real(kind=dp_t), dimension(nspec) :: Xin
  real(kind=dp_t), dimension(nspec_advance+1) :: y, ydot
  real(kind=dp_t) :: enucdot

  real(kind=dp_t) :: rpar
  integer :: ipar

  integer :: ic12, io16, img24
  integer :: n

  type (eos_t) :: eos_state

  call network_init()
  call eos_init()

  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")
  img24 = network_species_index("magnesium-24")

  if (ic12 < 0 .or. io16 < 0 .or. img24 < 0) then
     call bl_error("ERROR: species index not defined")
  endif

  dens = 2.6e9_dp_t
  temp = 6.e8_dp_t

  Xin(ic12) = 0.5_dp_t
  Xin(io16) = 0.5_dp_t
  Xin(img24) = 0.0_dp_t


  eos_state%rho = dens
  eos_state%T = temp
  eos_state%xn(:) = Xin(:)
  
  call eos(eos_input_rt, eos_state, .false.)

  print *, 'evaluating the RHS...'

  ! load the state
  y(1) = Xin(ic12)
  y(nspec_advance+1) = temp

  ! set the burner_aux variables
  dens_pass = dens
  c_p_pass = eos_state%cp
  dhdx_pass(:) = eos_state%dhdX(:)
  X_O16_pass = Xin(io16)
  

  call f_rhs(nspec_advance+1, ZERO, y, ydot, rpar, ipar)

  
  print *, 'done!'

  print *, 'Xin:  ', Xin
  print *, 'rhs:  ', ydot

  ! compute the energy release/s (erg/g/s)
  enucdot = (ebin(img24) - ebin(ic12))*ydot(ic12)

  print *, 'enucdot = ', enucdot

end program testburn
