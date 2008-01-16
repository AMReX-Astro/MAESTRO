! This module provides the burner for the dvode test problem.  This burner 
! should not be used in any real MAESTRO run.
!
! More information in the README file
!
module burner_module

  use bl_types
  use bl_constants_module
  use network
  use bl_error_module

contains

  subroutine burner(Xin, dt, tol, Xout)

    implicit none

    real(kind=dp_t), intent(in   ) :: Xin(nspec), dt, tol(2)
    real(kind=dp_t), intent(  out) :: Xout(nspec)
  
    ! allocate storage for the input state
    real(kind=dp_t), dimension(nspec) :: y, ydot


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

    integer, parameter :: ITOL = 4
    real(kind=dp_t), dimension(nspec) :: atol, rtol

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
    
    ! declare a real work array of size 22 + 9*nspec + 2*nspec**2 and an
    ! integer work array of since 30 + nspec
    integer, parameter :: LRW = 22 + 9*nspec + 2*nspec**2
    real(kind=dp_t), dimension(LRW) :: rwork
    
    integer, parameter :: LIW = 30 + nspec
    integer, dimension(LIW) :: iwork
    
    real(kind=dp_t) :: rpar
    integer :: ipar

    EXTERNAL jac, f_rhs, dvode
    
    integer :: i


    ! set the tolerances.  

    atol(:) = tol(1) 
    
    rtol(:) = tol(2) 


    ! we want VODE to re-initialize each time we call it
    istate = 1
    
    rwork(:) = ZERO
    iwork(:) = 0
    
    
    ! set the maximum number of steps allowed (the VODE default is 500)
    iwork(6) = 15000
    
    
    ! initialize the integration time
    time = ZERO
    
    ! abundances 
    y = Xin


    ! call the integration routine
    call dvode(f_rhs, nspec, y, time, dt, ITOL, rtol, atol, ITASK, &
               istate, IOPT, rwork, LRW, iwork, LIW, jac, MF_NUMERICAL_JAC,&
               rpar, ipar)

    if (istate < 0) then
       print *, 'ERROR: integration failed in net'
       print *, 'istate = ', istate
       print *, 'time = ', time
       call bl_error("ERROR in burner: integration failed")
    endif


    ! store the new mass fractions 
    Xout(:) = max(min(y(:), ONE), ZERO)

    
  end subroutine burner

end module burner_module
