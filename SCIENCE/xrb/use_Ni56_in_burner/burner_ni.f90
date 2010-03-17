! This module contains the triple-alpha reaction network burner.  
!
! Given the initial state of the system and the right-hand-side (f_rhs.f90) of
! a linear system of ODEs, this routine calls VODE to get the updated 
! composition, enthalpy and rho_omegadot.
! The temperature is evolved (in an isobaric, self-heating formalism) 
! concurrently with the species to prevent instabilities - see 
! Muller A&A 162, 103-108 (1986) for a disussion on the formation of 
! instabilities.  However, at the end of the call, the changes to temperature 
! are discarded.
!
! This burner provides an explicit Jacobian matrix to the DVODE solver.
!

module burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use network

  private
  public :: burner
  
contains

  subroutine burner(dens, temp, Xin, hin, dt, Xout, hout, rho_omegadot)

    implicit none

    real(kind=dp_t), intent(in   ) :: dens, temp, Xin(nspec), hin, dt
    real(kind=dp_t), intent(  out) :: Xout(nspec), hout, rho_omegadot(nspec)
  
    integer :: n
    real(kind=dp_t) :: enuc, dX

    logical, parameter :: verbose = .false.

    ! set the number of independent variables -- this should be temperature
    ! + the number of species which participate in the evolution equations
    ! 
    ! NOTE: fe56 does NOT participate in reactions BUT does come into play when
    !       determining abar and zbar for screening purposes and from partial
    !       derivs. in the EOS
    ! 
    integer, parameter :: NEQ = 1 + nspec

    ! allocate storage for the input state
    real(kind=dp_t), dimension(NEQ) :: y

    ! density, specific heat at constant pressure, c_p, and dhdX are needed
    ! in the righthand side routine, so we will pass these in through common
    ! Since evaluating the EOS is expensive, we don't call it for every RHS
    ! call -- T_eos and dT_crit control the frequency (see below)
    real(kind=dp_t) :: dens_zone, c_p_pass, dhdx_pass(nspec), T_eos, dT_crit
    common /zone_state/ dens_zone, c_p_pass, dhdx_pass, T_eos, dT_crit


    ! our problem is stiff, tell ODEPACK that. 21 means stiff, jacobian 
    ! function is supplied, 22 means stiff, figure out my jacobian through 
    ! differencing
    integer, parameter :: MF_ANALYTIC_JAC = 21, MF_NUMERICAL_JAC = 22


    ! tolerance parameters:
    !
    !  itol specifies whether to use an single absolute tolerance for
    !  all variables (1), or to pass an array of absolute tolerances, one
    !  for each variable with a scalar relative tol (2), a scalar absolute
    !  and array of relative tolerances (3), or arrays for both (4)
    !  
    !  The error is determined as e(i) = rtol*abs(y(i)) + atol, and must
    !  be > 0.  Since we have some compositions that may be 0 initially,
    !  we will specify both an absolute and a relative tolerance.
    !
    ! We will use arrays for both the absolute and relative tolerances, 
    ! since we want to be easier on the temperature than the species
    integer, parameter :: ITOL = 4
    real(kind=dp_t), dimension(NEQ) :: atol, rtol


    real(kind=dp_t) :: time
    

    ! we want to do a normal computation, and get the output values of y(t)
    ! after stepping though dt
    integer, PARAMETER :: ITASK = 1
  

    ! istate determines the state of the calculation.  A value of 1 meeans
    ! this is the first call to the problem -- this is what we will want.
    ! Note, istate is changed over the course of the calculation, so it
    ! cannot be a parameter
    integer :: istate


    ! we will override the maximum number of steps, so turn on the 
    ! optional arguments flag
    integer, parameter :: IOPT = 1
    
    ! declare a real work array of size 22 + 9*NEQ + 2*NEQ**2 and an
    ! integer work array of since 30 + NEQ
    integer, parameter :: LRW = 22 + 9*NEQ + 2*NEQ**2
    real(kind=dp_t), dimension(LRW) :: rwork
    
    integer, parameter :: LIW = 30 + NEQ
    integer, dimension(LIW) :: iwork
    
    ! indices of the species in the reaction network passed through 
    ! common block to other routines
    integer, save :: ihe4, ic12, ife56
    integer :: ihe4_pass, ic12_pass, ife56_pass
    common /species_index/ ihe4_pass, ic12_pass, ife56_pass

    ! indices of the reactions in the network passed through common
    integer, save :: ir3a, irg3a
    integer :: ir3a_pass, irg3a_pass
    common /reaction_index/ ir3a_pass, irg3a_pass

    real(kind=dp_t) :: rpar
    integer :: ipar

    EXTERNAL jac, f_rhs
    
    integer :: i

    logical, save :: firstCall = .true.

    if (firstCall) then

       if (.NOT. network_initialized) then
          call bl_error("ERROR in burner: must initialize network first")
       endif

       ihe4  = network_species_index("helium-4")
       ic12  = network_species_index("carbon-12")
       ife56 = network_species_index("nickel-56")

       ir3a  = network_reaction_index("forward")
       irg3a = network_reaction_index("backward")

       if (ihe4 < 0 .or. ic12 < 0 .or. ife56 < 0 .or. &
            ir3a < 0 .or. irg3a < 0) then
          call bl_error("ERROR in burner: species undefined")
       endif
       
       firstCall = .false.
    endif

    ihe4_pass  = ihe4
    ic12_pass  = ic12
    ife56_pass = ife56

    ir3a_pass  = ir3a
    irg3a_pass = irg3a

    
    ! set the parameters regarding how often to re-evaluate the 
    ! thermodynamics.  T_eos will always store the temperature
    ! that was used for the last EOS call.  dT_crit is the 
    ! relative temperature change beyond which we need to re-evaluate
    ! the thermodynamics
    !
    ! **NOTE** if the burner is not converging (and the temperatures
    ! are shooting up to unrealistically high values), you likely
    ! need to reduce dT_crit to ensure more frequent EOS calls.
    T_eos = temp
    dT_crit = 1.d20


    ! set the tolerances.  We will be more relaxed on the temperature
    ! since it is only used in evaluating the rates.  
    !
    ! **NOTE** if you reduce these tolerances, you probably will need
    ! to (a) decrease dT_crit, (b) increase the maximum number of 
    ! steps allowed.
    atol(1:nspec) = 1.d-12    ! mass fractions
    atol(nspec+1) = 1.d-8     ! temperature

    rtol(1:nspec) = 1.d-12    ! mass fractions
    rtol(nspec+1) = 1.d-5     ! temperature


    ! we want VODE to re-initialize each time we call it
    istate = 1

    rwork(:) = ZERO
    iwork(:) = 0


    ! set the maximum number of steps allowed (the VODE default is 500)
    iwork(6) = 15000


    ! initialize the integration time
    time = ZERO

    ! abundances are the first nspec-1 values and temperature is the last
    y(ihe4)    = Xin(ihe4)
    y(ic12)    = Xin(ic12)
    y(ife56)   = Xin(ife56)
    y(nspec+1) = temp

    ! set the thermodynamics that are passed via common to the RHS routine--
    ! these will be updated in f_rhs if the relative temperature change 
    ! exceeds dT_crit
    dens_zone = dens

    ! we need the specific heat at constant pressure and dhdX |_p.  Take
    ! T, rho, Xin as input
    den_eos(1) = dens
    temp_eos(1) = temp
    xn_eos(1,:) = Xin(:)

    call eos(eos_input_rt, den_eos, temp_eos, &
             npts, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             .false.)

    c_p_pass = cp_eos(1)
    dhdx_pass(:) = dhdX_eos(1,:)


    ! call the integration routine
    call dvode(f_rhs, NEQ, y, time, dt, ITOL, rtol, atol, ITASK, &
         istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_ANALYTIC_JAC,&
         rpar, ipar)

    if (istate < 0) then
       print *, 'ERROR: integration failed in net'
       print *, 'istate = ', istate
       print *, 'time = ', time
       call bl_error("ERROR in burner: integration failed")
    endif


    ! store the new mass fractions -- note, we discard the temperature
    ! here and instead compute the energy release from the binding
    ! energy -- make sure that they are positive
    Xout(ihe4)  = max(y(ihe4), ZERO)
    Xout(ic12)  = max(y(ic12), ZERO)
    Xout(ife56) = max(y(ife56), ZERO)


    ! compute the energy release and update the enthalpy.  Our convention
    ! is that the binding energies are negative, so the energy release is
    ! - sum_k { (Xout(k) - Xin(k)) ebin(k) }
    enuc = 0.0_dp_t
    do n = 1, nspec
       dX = Xout(n) - Xin(n) 

       enuc = enuc - ebin(n) * dX

       rho_omegadot(n) = dens * dX / dt
    enddo


    hout = hin + enuc

    if (verbose) then

       ! print out some integration statistics, if desired
       print *, 'integration summary: '
       print *, 'dens: ', dens, ' temp: ', temp
       print *, 'number of steps taken: ', iwork(11)
       print *, 'number of f evaluations: ', iwork(12)
    endif


  end subroutine burner

end module burner_module
