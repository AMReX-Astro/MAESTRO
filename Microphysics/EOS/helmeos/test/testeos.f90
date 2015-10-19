program testeos

  use bl_types
  use network
  use eos_module
  use eos_type_module

  implicit none

  real(kind=dp_t) :: dens, temp, pres, entr
  real(kind=dp_t), dimension(nspec) :: Xin
  type(eos_t) :: eos_state
  
  integer :: ic12, io16, img24

  logical :: do_diag

  call network_init()
  call eos_init()

  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")
  img24 = network_species_index("magnesium-24")

  dens = 2.6e9_dp_t
  temp = 1.e9_dp_t

  pres = 1.7e27_dp_t
  entr = 6.5e7_dp_t

  Xin(ic12) = 0.5_dp_t
  Xin(io16) = 0.5_dp_t
  Xin(img24) = 0.0_dp_t


  eos_state%rho = dens
  eos_state%T = temp
  eos_state%xn(:) = Xin(:)

  do_diag = .false.

  call eos(eos_input_rt, eos_state)

  print *, 'eos_input_rt:'
  print *, 'dens: ', dens, ' temp: ', temp
  print *, 'X: ', Xin
  print *, 'pres: ', eos_state%p,  ' ener: ', eos_state%e
  print *, 'h:    ', eos_state%h,  ' entr: ', eos_state%s
  print *, 'c_v:  ', eos_state%cv, ' c_p : ', eos_state%cp
  print *, 'dpdT: ', eos_state%dpdt, ' dpdr: ', eos_state%dpdr
  print *, 'dedT: ', eos_state%dedt, ' dedr: ', eos_state%dedr
  print *, 'dpdX: ', eos_state%dpdX(:)
  print *, 'dhdX: ', eos_state%dhdX(:)

  eos_state%p = pres
  eos_state%s = entr

  call eos(eos_input_ps, eos_state, do_diag)

  print *, 'eos_input_ps:'
  print *, 'dens: ', eos_state%rho, ' temp: ', eos_state%T
  print *, 'X: ', Xin
  print *, 'pres: ', pres,  ' entr: ', entr
  print *, 'p_eos: ', eos_state%p, ' s_eos: ', eos_state%s
  print *, 'h:    ', eos_state%h,  ' ener: ', eos_state%e
  print *, 'c_v:  ', eos_state%cv, ' c_p : ', eos_state%cp
  print *, 'dpdT: ', eos_state%dpdt, ' dpdr: ', eos_state%dpdr
  print *, 'dedT: ', eos_state%dedt, ' dedr: ', eos_state%dedr
  print *, 'dpdX: ', eos_state%dpdX(:)
  print *, 'dhdX: ', eos_state%dhdX(:)
  

end program testeos
