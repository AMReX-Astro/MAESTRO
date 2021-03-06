! This example is based on
!
!   http://www.radford.edu/~thompson/vodef90web/vodef90source/Double/Prologue_Examples/example1.f90
!
! Some notes from example1.f90:
!
!   The problem is from chemical kinetics, and consists of the following
!   three rate equations:
!
!     dy1/dt = -.04d0*y1 + 1.d4*y2*y3
!     dy2/dt = .04d0*y1 - 1.d4*y2*y3 - 3.d7*y2**2
!     dy3/dt = 3.d7*y2**2
!
!   on the interval from t = 0.0d0 to t = 4.d10, with initial
!   conditions y1 = 1.0d0, y2 = y3 = 0.0d0. The problem is stiff.
!
! Here we're going to solve the same problem, except:
!   * we'll evolve two solutions / initial conditions at the same time
!   * the Jacobian will be computed using the first solution only
!


module feval
  use bdf
  use bl_types
  implicit none
  integer, parameter :: neq = 3
  integer, parameter :: npt = 2
contains
  subroutine f(neq, y, t, ydot, upar)
    !$acc routine seq
    integer,  intent(in   ) :: neq
    real(dp_t), intent(in   ) :: y(neq), t
    real(dp_t), intent(  out) :: ydot(neq)
    real(dp_t), intent(inout), optional :: upar(:)
    
    ydot(1) = -.04d0*y(1) + 1.d4*y(2)*y(3)
    ydot(3) = 3.e7*y(2)*y(2)
    ydot(2) = -ydot(1) - ydot(3)
  end subroutine f
  subroutine J(neq, y, t, pd, upar)
    !$acc routine seq
    integer,  intent(in   ) :: neq
    real(dp_t), intent(in   ) :: y(neq), t
    real(dp_t), intent(  out) :: pd(neq,neq)
    real(dp_t), intent(inout), optional :: upar(:)
    
    pd(1,1) = -.04d0
    pd(1,2) = 1.d4*y(3)
    pd(1,3) = 1.d4*y(2)
    pd(2,1) = .04d0
    pd(2,3) = -pd(1,3)
    pd(3,2) = 6.e7*y(2)
    pd(2,2) = -pd(1,2) - pd(3,2)
  end subroutine J
end module feval

program test
  use bdf
  use feval
  implicit none

  type(bdf_ts)  :: ts
  double precision :: rtol(neq), atol(neq), dt
  double precision :: y0(neq,npt), t0, y1(neq,npt), t1

  integer :: i, ierr

  y0(:,1) = [ 1.d0, 0.d0, 0.d0 ]
  y0(:,2) = [ 0.98516927747181138d0, 3.3863452485889568d-5, 1.4796859075703273d-2 ]

  t0 = 0.d0
  t1 = 0.4d0

  rtol = 1.d-4
  atol = [ 1.d-8, 1.d-14, 1.d-6 ]
  dt   = 1.d-8

  call bdf_ts_build(ts, neq, npt, rtol, atol, max_order=3)

  do i = 1, 11
     call bdf_advance(ts, f, J, neq, npt, y0, t0, y1, t1, dt, .true., .false., ierr)
     print *, t1, y1(:,1), ierr, errors(ierr)
     print *, '                         ', y1(:,2), ts%dt
     if (ierr /= BDF_ERR_SUCCESS) exit
     y0 = y1
     t0 = t1
     t1 = 10*t1
     dt = 2*maxval(ts%dt)
  end do

  print *, ''
  print *, 'max stats for last interval'
  print *, 'number of steps taken      ', ts%n
  print *, 'number of function evals   ', maxval(ts%nfe)
  print *, 'number of jacobian evals   ', maxval(ts%nje)
  print *, 'number of lu decomps       ', maxval(ts%nlu)
  print *, 'number of solver iterations', maxval(ts%nit)
  print *, 'number of solver errors    ', maxval(ts%nse)

  call bdf_ts_destroy(ts)

end program test
