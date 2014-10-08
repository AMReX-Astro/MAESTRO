program evolve

  use network
  use burner_module
  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use eos_type_module

  implicit none

  real(kind=dp_t), parameter :: tstart = ZERO

  real(kind=dp_t), parameter :: dX_tol = 1.e-10
  real(kind=dp_t) :: dX_max

  real(kind=dp_t) :: t, tstop, dt
  real(kind=dp_t) :: temp, dens

  real(kind=dp_t) :: Xin(nspec), Xout(nspec), rho_omegadot(nspec), rho_Hnuc
  real(kind=dp_t) :: hin, hout
  real(kind=dp_t) :: p_fixed

  type (eos_t) :: eos_state

  ! initialize
  call eos_init()
  call network_init()

  t = tstart
  tstop = 100.0d0

  ! initial conditions
  dens = 2.d6
  temp = 1.e9

  Xin = ZERO

  ! approximately solar
  Xin(ih1) = 0.72_dp_t
  Xin(ihe4) = 0.24_dp_t
  Xin(io14) = 0.004_dp_t
  Xin(io15) = 0.03_dp_t
  Xin(img22) = 0.001_dp_t
!  Xin(ic12) = 0.001_dp_t
  Xin(is30) = 0.001_dp_t
  Xin(ini56) = 1.0_dp_t - sum(Xin(:))

  eos_state%rho = dens
  eos_state%T   = temp
  eos_state%xn  = Xin


  call eos(eos_input_rt, eos_state, .false.)

  ! the fixed pressure for this isobaric burn
  p_fixed = eos_state%p


  Xout = ZERO

  ! pick an initial dt -- we will adjust dt to allow dX to change by at
  ! most dX_tol from step to step, looking only at H and He
  dt = 1.e-12_dp_t

  do while(t <= tstop)

     ! get the initial enthalpy
     eos_state%rho = dens
     eos_state%T   = temp
     eos_state%xn  = Xin
     
     call eos(eos_input_rt, eos_state, .false.)

     hin = eos_state%h

     call burner(dens, temp, Xin, dt, Xout, rho_omegadot, rho_Hnuc)

     hout = hin + rho_Hnuc*dt/dens

     ! in the burner, we added energy but kept both the density and     
     ! pressure fixed -- this cannot be.  But our temperature
     ! evolution equation acted as if the pressure were constant.
     ! Here we reconcile that.

     ! get the new temperature in the zone, from the old density
     ! and the hout = hin + rho_Hnuc/rho
     eos_state%rho = dens
     eos_state%h   = hout
     eos_state%xn  = xout
     eos_state%p = p_fixed
        
     call eos(eos_input_ph, eos_state, .false.)
        
     temp = eos_state%T
     dens = eos_state%rho


     ! get the new dt -- we want to constrain the change in dX over
     ! the timestep.  so dX/dt = dX_tol/dt_want -- this would give the
     ! dt_want such that the change in dX should have been dX_tol
     dX_max = max(abs(Xout(ih1) - Xin(ih1))/Xin(ih1), abs(Xout(ihe4) - Xin(ihe4))/Xin(ihe4))
     dt = max(5.0*dt,(dX_tol/dX_max)*dt)
     dt = min(dt, 2.0d0)

     if (t == tstart) then
        print *, "time, dens, temp, X(H1), X(He4), X(C12), Hnuc (erg/g/s)"
     endif

     t = t + dt
     Xin = Xout

     ! output

     print *, t, real(dens), real(temp), real(Xout(ih1)), real(Xout(ihe4)), real(Xout(ic12)), real(rho_Hnuc/dens)


  end do

end program evolve
