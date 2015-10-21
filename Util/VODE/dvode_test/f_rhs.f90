! These routines provide the right-hand-side and a dummy jacobian for the 
! DVODE solver.  More information about how ODE's can be found in the README 
! file.
!
subroutine f_rhs(n, t, y, ydot, rpar, ipar)

  use bl_types

  implicit none

  integer,         intent(IN   ) :: n, ipar
  real(kind=dp_t), intent(IN   ) :: t, y(n), rpar
  real(kind=dp_t), intent(  OUT) :: ydot(n)

  call dydt(n, y, ydot)

  return

end subroutine f_rhs
  



subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  use bl_types

  integer :: neq
  integer :: ml, mu, nrpd

  real(kind=dp_t) :: y(neq)
  real(kind=dp_t) :: pd(neq,neq)

  real(kind=dp_t) :: rpar
  integer :: ipar

  real(kind=dp_t) :: t


  return
end subroutine jac




subroutine dydt(n, y,ydot)

  use network
  use bl_types

  integer,            intent(IN   ) :: n
  real(kind=dp_t),    intent(IN   ) :: y(n)
  real(kind=dp_t),    intent(  OUT) :: ydot(n)

  integer :: i, j


  do i = 1, n

     ydot(i) = 0.0d0

     do j = 1, n

        ydot(i) = ydot(i) + y(j) * rates(i,j)

     enddo

  enddo

  return
end subroutine dydt
