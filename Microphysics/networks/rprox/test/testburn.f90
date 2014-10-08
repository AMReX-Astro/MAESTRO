program test_burn

  use bl_types
  use bl_constants_module
  use bl_error_module
  use network
  use burner_module
  use eos_module

  implicit none

  real(kind=dp_t) :: temp, dens, dt
  real(kind=dp_t) :: Xin(nspec), Xout(nspec), rhoHnuc, rho_omegadot(nspec)


  call eos_init()
  call network_init()

  temp = 7.5e8_dp_t
  dens = 1.0e2_dp_t
  dt = 1.0e-6

  Xin = ZERO
  Xout = ZERO

  Xin(ih1) = 0.7_dp_t
  Xin(ihe4) = 0.3_dp_t

  call burner(dens,temp,Xin,dt,Xout,rho_omegadot,rhoHnuc)

  print *, Xin
  print *,'........'
  print *, Xout
  print *, '.......'
  print *, rho_omegadot
  print *, '.......'
  print *, log10(rhoHnuc/dens)

end program test_burn
  
