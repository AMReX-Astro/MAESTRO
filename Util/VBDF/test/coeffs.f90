module feval
  use bdf
  use bl_types
  implicit none
  integer, parameter :: neq = 1
contains
  subroutine f(neq, npt, y, t, ydot, upar)
    integer,    intent(in)    :: neq, npt
    real(dp_t), intent(in)    :: y(neq,npt), t(npt)
    real(dp_t), intent(  out) :: ydot(neq,npt)
    real(dp_t), intent(inout), optional :: upar(:,:)
  end subroutine f
  subroutine J(neq, npt, y, t, pd, upar)
    integer,    intent(in)    :: neq, npt
    real(dp_t), intent(in)    :: y(neq,npt), t(npt)
    real(dp_t), intent(  out) :: pd(neq,neq,npt)
    real(dp_t), intent(inout), optional :: upar(:,:)
  end subroutine J
end module feval

program test
  use bdf
  use feval

  double precision, parameter :: tol = 1.d-12
  integer,          parameter :: NPT = 1

  type(bdf_ts)  :: ts
  double precision :: rtol(neq), atol(neq)
  double precision :: v(0:6)

  rtol = 1
  atol = 1

  call bdf_ts_build(ts, neq, 1, rtol, atol, max_order=6)

  print *, "====> order 1"
  ts%k = 1
  call random_number(ts%h(:,NPT))
  call bdf_update(ts)
  print *, 'l', ts%l(0:1,NPT)
  call assert(all(ts%l(0:1,NPT) == [ 1.d0, 1.d0 ]), "error in l")

  ts%h = 1
  call bdf_update(ts)
  print *, 'tq', ts%tq
  call assert(all(abs(ts%tq(:,NPT) - [ 1.d0, 0.5d0, 0.2222222222222d0, 2.0d0 ]) < tol), "error in tq")

  print *, "====> order 2"
  ts%k = 2
  call random_number(ts%h(:,NPT))
  call bdf_update(ts)
  print *, 'l', ts%l(0:2,NPT)
  call assert(all(abs(ts%l(0:2,NPT) - [ 1.d0, 1.5d0, 0.5d0 ]) < tol), "error in l")

  ts%h = 1
  call bdf_update(ts)
  print *, 'tq', ts%tq
  call assert(all(abs(ts%tq(:,NPT) - [ 1.0d0, 0.222222222222222d0, 0.13636363636363638d0, 6.0d0 ]) < tol), "error in tq")

  print *, "====> order 3"
  ts%k = 3
  ts%h = 1
  call bdf_update(ts)
  print *, 'l', ts%l(0:3,NPT)
  v(0:3) = [ 1.d0, 1.8333333333333d0, 1.d0, 0.1666666666666d0 ]
  call assert(all(abs(ts%l(0:3,NPT) - v(0:3)) < tol), "error in l")

  ts%h(1,NPT) = 4
  call bdf_update(ts)
  print *, 'l', ts%l(0:3,NPT)
  v(0:3) = [ 1.d0, 1.8333333333333d0, 0.96d0, 0.12666666666666d0 ]
  call assert(all(abs(ts%l(0:3,NPT) - v(0:3)) < tol), "error in l")

  ts%h = 1
  call bdf_update(ts)
  print *, 'tq', ts%tq
  call assert(all(abs(ts%tq(:,NPT) - [ 1.3333333333333d0, 0.13636363636363638d0, 9.6d-2, 24.0d0 ]) < tol), "error in tq")
  
  call bdf_ts_destroy(ts)

contains

  subroutine assert(cond, message)
    logical, intent(in) :: cond
    character(*), intent(in) :: message

    if (.not. cond) then
       print *, "ASSERTION FAILED: ", message
       call abort()
    end if
  end subroutine assert

end program test
