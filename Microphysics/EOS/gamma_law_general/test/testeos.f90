program testburn

  use bl_types
  use network
  use eos_module

  implicit none

  real(kind=dp_t) :: dens, temp, pres, entr
  real(kind=dp_t), dimension(nspec) :: Xin
  
  integer :: ic12, io16, img24

  logical :: do_diag

  call network_init()
  call eos_init(gamma_in=5.0d0/3.0d0)

  do_diag = .false.

  ic12 = network_species_index("carbon-12")
  io16 = network_species_index("oxygen-16")
  img24 = network_species_index("magnesium-24")

  dens = 1.e4_dp_t
  temp = 1.e8_dp_t

  Xin(ic12) = 0.5_dp_t
  Xin(io16) = 0.5_dp_t
  Xin(img24) = 0.0_dp_t


  den_eos = dens
  temp_eos = temp
  xn_eos(:) = Xin(:)

  
  ! test eos_input_rt
  call eos(eos_input_rt, den_eos, temp_eos, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           do_diag)

  print *, 'input: eos_input_rt'
  print *, 'dens: ', dens, ' temp: ', temp
  print *, 'X: ', Xin
  print *, 'pres: ', p_eos,  ' ener: ', e_eos
  print *, 'h:    ', h_eos
  print *, 'c_v:  ', cv_eos, ' c_p : ', cp_eos
  print *, 'dpdT: ', dpdt_eos, ' dpdr: ', dpdr_eos
  print *, 'dedT: ', dedt_eos, ' dedr: ', dedr_eos
  print *, 'dpdX: ', dpdX_eos(:)
  print *, 'dhdX: ', dhdX_eos(:)

  print *, ' '
  print *, 'setting pres = ', p_eos
  pres = p_eos
  print *, 'setting entr = ', s_eos
  entr = s_eos

  
  ! test eos_input_rh
  temp_eos = 0.d0
  call eos(eos_input_rh, den_eos, temp_eos, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           do_diag)

  print *, ' '
  print *, 'input: eos_input_rh'
  print *, 'dens: ', dens, ' temp: ', temp
  print *, 'X: ', Xin
  print *, 'pres: ', p_eos,  ' ener: ', e_eos
  print *, 'entropy: ', s_eos
  print *, 'h:    ', h_eos
  print *, 'c_v:  ', cv_eos, ' c_p : ', cp_eos
  print *, 'dpdT: ', dpdt_eos, ' dpdr: ', dpdr_eos
  print *, 'dedT: ', dedt_eos, ' dedr: ', dedr_eos
  print *, 'dpdX: ', dpdX_eos(:)
  print *, 'dhdX: ', dhdX_eos(:)

  ! test eos_input_tp
  den_eos = 0.d0
  call eos(eos_input_tp, den_eos, temp_eos, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           do_diag)

  print *, ' '
  print *, 'input: eos_input_tp'
  print *, 'dens: ', dens, ' temp: ', temp
  print *, 'X: ', Xin
  print *, 'pres: ', p_eos,  ' ener: ', e_eos
  print *, 'h:    ', h_eos
  print *, 'c_v:  ', cv_eos, ' c_p : ', cp_eos
  print *, 'dpdT: ', dpdt_eos, ' dpdr: ', dpdr_eos
  print *, 'dedT: ', dedt_eos, ' dedr: ', dedr_eos
  print *, 'dpdX: ', dpdX_eos(:)
  print *, 'dhdX: ', dhdX_eos(:)

  ! test eos_input_rp
  temp_eos = 0.d0
  call eos(eos_input_rp, den_eos, temp_eos, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           do_diag)

  print *, ' '
  print *, 'input: eos_input_rp'
  print *, 'dens: ', dens, ' temp: ', temp
  print *, 'X: ', Xin
  print *, 'pres: ', p_eos,  ' ener: ', e_eos
  print *, 'h:    ', h_eos
  print *, 'c_v:  ', cv_eos, ' c_p : ', cp_eos
  print *, 'dpdT: ', dpdt_eos, ' dpdr: ', dpdr_eos
  print *, 'dedT: ', dedt_eos, ' dedr: ', dedr_eos
  print *, 'dpdX: ', dpdX_eos(:)
  print *, 'dhdX: ', dhdX_eos(:)

  ! test eos_input_re
  temp_eos = 0.d0
  call eos(eos_input_re, den_eos, temp_eos, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           do_diag)

  print *, ' '
  print *, 'input: eos_input_re'
  print *, 'dens: ', dens, ' temp: ', temp
  print *, 'X: ', Xin
  print *, 'pres: ', p_eos,  ' ener: ', e_eos
  print *, 'h:    ', h_eos
  print *, 'c_v:  ', cv_eos, ' c_p : ', cp_eos
  print *, 'dpdT: ', dpdt_eos, ' dpdr: ', dpdr_eos
  print *, 'dedT: ', dedt_eos, ' dedr: ', dedr_eos
  print *, 'dpdX: ', dpdX_eos(:)
  print *, 'dhdX: ', dhdX_eos(:)


  ! test eos_input_ps
  p_eos = pres
  s_eos = entr
  call eos(eos_input_ps, den_eos, temp_eos, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           do_diag)

  print *, ' '
  print *, 'input: eos_input_ps'
  print *, 'pres: ', pres, 'entr: ', entr
  print *, 'dens: ', den_eos, ' temp: ', temp_eos
  print *, 'X: ', Xin
  print *, 'pres_eos: ', p_eos,  'entr_eos: ', s_eos
  print *, 'h:    ', h_eos, ' ener: ', e_eos
  print *, 'c_v:  ', cv_eos, ' c_p : ', cp_eos
  print *, 'dpdT: ', dpdt_eos, ' dpdr: ', dpdr_eos
  print *, 'dedT: ', dedt_eos, ' dedr: ', dedr_eos
  print *, 'dpdX: ', dpdX_eos(:)
  print *, 'dhdX: ', dhdX_eos(:)

end program testburn
