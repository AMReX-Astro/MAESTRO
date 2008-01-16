! This program gives a test driver for the DVODE test problem.  
!
! Changing the value of max(alpha, beta) / min(alpha, beta) changes the 
! stiffness of the problem.  A very stiff problem, such as reaction networks,
! really gives DVODE a nice workout.  Adjusting the atol and rtol parameters
! can give some insight as to what these values should be to give a reasonable
! result for some given initial conditions.  More information about this 
! reaction network and the equations involved can be found in ../README.
!
! The function xbxa gives the analytic solution to the system given the 
! initial condition xb(t=0) = 0 and of course the condition sum_i x_i = 1.
!
program testburn

  use burner_module
  use network
  use bl_types
  use bl_constants_module

  implicit none

  real(kind=dp_t) :: x_dvode(nspec), x_new(nspec), tol(2)
  real(kind=dp_t) :: analytic_soln, dvode_soln

  real(kind=dp_t) :: tstep, time

  integer :: i, j, k

  ! number of burn steps to take
  integer, parameter :: imax = 1000

  ! rate parameters
  real(kind=dp_t), parameter :: alpha = 8.0d9, &
                                 beta = 4.0d-6

  ! tolerance parameters
  real(kind=dp_t), parameter :: atol = 1.0d-15, &
                                rtol = 1.0d-12

  ! species information
  integer :: ia, ib


  if (.not. network_initialized) then
     call network_init

     ia = network_species_index("A")
     ib = network_species_index("B")

  endif
     
  
! initial conditions
  x_dvode(ia) = ONE
  x_dvode(ib) = ZERO

  rates(ia,ia) = -alpha
  rates(ib,ia) =  alpha
  rates(ib,ib) =  -beta
  rates(ia,ib) =   beta

  time = ZERO

  ! some value of timestep...
  ! this doesn't seem to affect the solution too much as DVODE will take 
  ! a corrected timestep if this does not work
  tstep = 2.d-11

  tol(1) = atol
  tol(2) = rtol


  print *, ' alpha = ', alpha
  print *, '  beta = ', beta
  print *, ' tstep = ', tstep
  print *, '  atol = ', atol
  print *, '  rtol = ', rtol
  print *, ''
  print *, '  time           analytic               DVODE     ', &
           '        % difference'
  

! open a file for output
  open(unit=1,file='testburn.out',status='unknown',action='write')
  write(1,*)'# alpha = ', alpha
  write(1,*)'#  beta = ', beta
  write(1,*)'# tstep = ', tstep
  write(1,*)'#  atol = ', atol
  write(1,*)'#  rtol = ', rtol
  write(1,*)'# time           analytic               DVODE     ', &
           '        % difference'

  do i = 1, imax
     time = i * tstep

     call burner(x_dvode, tstep, tol, x_new)

     x_dvode(:) = x_new(:)

     dvode_soln = x_dvode(ib)/x_dvode(ia)

     analytic_soln = xbxa(time, alpha, beta)


     write(6,100) time, analytic_soln, dvode_soln,                          &
              (abs(analytic_soln - dvode_soln)/analytic_soln) * 1.d2
     ! dump out to a file
     write(1,100) time, analytic_soln, dvode_soln,                          &
              (abs(analytic_soln - dvode_soln)/analytic_soln) * 1.d2

  enddo

  close(1)

100 format (ES10.3, 3(1x, ES20.12))

contains
  
  function xbxa(t, a, b)

    real(kind=dp_t) :: t, a, b, xbxa


    xbxa = (1 - dexp(-(a + b)*t)) / ((b/a)  +  dexp(-(a + b)*t))

    return

  end function xbxa


end program testburn




