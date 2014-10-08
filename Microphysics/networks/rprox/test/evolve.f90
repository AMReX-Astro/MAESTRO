program evolve

  use network
  use burner_module
  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module

  implicit none

  real(kind=dp_t), parameter :: tstart = ZERO, &
       tstop = 1.0e5_dp_t, &
       dt = 1.0e-1_dp_t, &
       temp = 7.5e8_dp_t, &
       dens = 1.0e2_dp_t

  real(kind=dp_t) :: Xin(nspec), Xout(nspec), rho_omegadot(nspec), rho_Hnuc, t
  
  call eos_init()
  call network_init()

  Xin = ZERO
  Xout = ZERO
  t = tstart

  Xin(ih1) = 0.7_dp_t
  Xin(ihe4) = 0.3_dp_t

  do while(t <= tstop)

!     write(*,*) '...', Xin

     call burner(dens,temp, Xin, dt, Xout, rho_omegadot, rho_Hnuc)

     write(*,'(12g16.8)') log10(t), rho_omegadot, log10(rho_Hnuc/dens)
     
     t = t + dt
     Xin = Xout
     
  end do

end program evolve
