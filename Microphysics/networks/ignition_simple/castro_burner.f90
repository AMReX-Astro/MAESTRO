module castro_burner_module

  use bl_types
  use bl_constants_module
  use bl_error_module
  use eos_module
  use network

  private
  public :: burner

contains

  subroutine burner(dens, temp, Xin, ein, dt, time_in, Xout, eout)

    implicit none

    real(kind=dp_t), intent(in   ) :: dens, temp, Xin(nspec), ein, dt, time_in
    real(kind=dp_t), intent(  out) :: Xout(nspec), eout
  
    integer :: n
    real(kind=dp_t) :: enuc, dX

    logical, parameter :: verbose = .false.

    ! set the number of independent variables -- this should be temperature
    ! + the number of species
    integer, parameter :: NEQ = 1 + nspec

    ! allocate storage for the input state
    real(kind=dp_t), dimension(NEQ) :: y

    ! density, specific heat at constant pressure, c_p, and dhdX are needed
    ! in the righthand side routine, so we will pass these in through common
    ! Since evaluating the EOS is expensive, we don't call it for every RHS
    ! call -- T_eos and dT_crit control the frequency (see below)
    real(kind=dp_t) :: dens_zone, c_p_pass, dhdx_pass(nspec), T_eos, dT_crit
    common /zone_state/ dens_zone, c_p_pass, dhdx_pass, T_eos, dT_crit

    ! we will always refer to the species by integer indices that come from
    ! the network module -- this makes things robust to a shuffling of the 
    ! species ordering
    integer, save :: ic12, io16, img24
    integer :: ic12_pass, io16_pass, img24_pass
    common /species_index/ ic12_pass, io16_pass, img24_pass

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

    real(kind=dp_t) :: integration_time

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
    

    real(kind=dp_t) :: rpar
    integer :: ipar

    EXTERNAL jac, f_rhs
    
    logical, save :: firstCall = .true.

    if (firstCall) then

       if (.NOT. network_initialized) then
          call bl_error("ERROR in burner: must initialize network first")
       endif
     
       ic12 = network_species_index("carbon-12")
       io16 = network_species_index("oxygen-16")
       img24 = network_species_index("magnesium-24")
       
       if (ic12 < 0 .OR. io16 < 0 .OR. img24 < 0) then
          call bl_error("ERROR in burner: species undefined")
       endif
       
       firstCall = .false.
    endif

    ic12_pass  = ic12
    io16_pass  = io16
    img24_pass = img24

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
    integration_time = ZERO
    
    ! abundances are the first nspec values and temperature is the last
    y(ic12) = Xin(ic12)
    y(io16) = Xin(io16)
    y(img24) = Xin(img24)
    y(nspec+1) = temp

    ! set the thermodynamics that are passed via common to the RHS routine --
    ! these will be updated in f_rhs if the relative temperature change 
    ! exceeds dT_crit
    dens_zone = dens
    
    ! we need the specific heat at constant pressure and dhdX |_p.  Take
    ! T, rho, Xin as input
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
    
    c_p_pass = cp_eos(1)
    dhdx_pass(:) = dhdX_eos(1,:)

    ! call the integration routine
    call dvode(f_rhs, NEQ, y, integration_time, dt, ITOL, rtol, atol, ITASK, &
               istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_ANALYTIC_JAC, &
               rpar, ipar)

    if (istate < 0) then
       print *, 'ERROR: integration failed in net'
       print *, 'istate = ', istate
       print *, 'time = ', integration_time
       call bl_error("ERROR in burner: integration failed")
    endif


    ! store the new mass fractions -- note, we discard the temperature
    ! here and instead compute the energy release from the binding
    ! energy -- make sure that they are positive
    Xout(ic12)  = max(y(ic12), ZERO)
    Xout(io16)  = max(y(io16), ZERO)
    Xout(img24) = max(y(img24), ZERO)

    ! compute the energy release and update the enthalpy.  Our convention
    ! is that the binding energies are negative, so the energy release is
    ! - sum_k { (Xout(k) - Xin(k)) ebin(k) }
    enuc = 0.0_dp_t
    do n = 1, nspec
       dX = Xout(n) - Xin(n) 
       enuc = enuc - ebin(n) * dX
    enddo
    
    eout = ein + enuc

    if (verbose) then
       
       ! print out some integration statistics, if desired
       print *, 'integration summary: '
       print *, 'dens: ', dens, ' temp: ', temp
       print *, 'number of steps taken: ', iwork(11)
       print *, 'number of f evaluations: ', iwork(12)
    endif

  end subroutine burner

end module castro_burner_module
