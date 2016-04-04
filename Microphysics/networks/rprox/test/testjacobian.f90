! a simple code to check the analytic Jacobian via numerical 
! differencing

program testjacobian

  use bl_types
  use bl_constants_module
  use bl_error_module
  use network
  use eos_module
  use eos_type_module
  use burner_module
  use rpar_indices

  implicit none

  real(kind=dp_t) :: dens, temp
  integer, parameter :: neq = nspec+1
  real(kind=dp_t), dimension(nspec) :: Xin
  real(kind=dp_t), dimension(neq) :: y, ydot
  real(kind=dp_t), dimension(neq) :: yp, ym
  real(kind=dp_t), dimension(neq) :: ydotp, ydotm
  real(kind=dp_t), dimension(neq,neq) :: pd
  real(kind=dp_t) :: enucdot

  real(kind=dp_t), allocatable :: rpar(:)
  integer :: ipar

  type(eos_t) :: eos_state

  integer :: i, j, n

  real(kind=dp_t), parameter :: delta = 0.001d0
  real(kind=dp_t), parameter :: SMALL = 1.d-12
  real(kind=dp_t) :: num_jac

  character(len=16) :: namei,namej

  call network_init()
  call eos_init()

  allocate(rpar(n_rpar_comps))

  dens = 2.6e9_dp_t
  temp = 7.e8_dp_t

  Xin = 0.005d0
  Xin(ih1) = 0.70d0
  Xin(ihe4) = 0.25d0

  rpar(irp_dens) = dens
  rpar(irp_T9_eos) = temp/1e9

  eos_state%rho = dens
  eos_state%T   = temp
  eos_state%xn  = Xin

  call eos(eos_input_rt, eos_state)

  rpar(irp_cp) = eos_state%cp
  rpar(irp_dhdX:irp_dhdX+nspec-1) = eos_state%dhdX

  print *, 'evaluating the RHS...'

  ! load the state
  y(1:nspec) = Xin / aion
  y(neq) = temp/1e9

  ! load the common block

  call f_rhs(neq, ZERO, y, ydot, rpar, ipar)

  call jac(neq, ZERO, y, 0, 0, pd, neq, rpar, ipar)

888 format(a,"-derivatives that don't match:")
999 format(5x, "df(",a,")/dy(",a,")", g18.10, g18.10, g18.10)

  do j = 1, neq

     yp(:) = y(:)
     ym(:) = y(:)

     yp(j) = (1.d0 + delta)*y(j) 
     call f_rhs(neq, ZERO, yp, ydotp, rpar, ipar)
     
     ym(j) = (1.d0 - delta)*y(j) 
     call f_rhs(neq, ZERO, ym, ydotm, rpar, ipar)        
     
     namej = short_spec_names(j)
     if (j==neq) namej = "T"

     write(*,888) trim(namej)

     do i = 1, neq
        
        num_jac = (ydotp(i) - ydotm(i))/(yp(j) - ym(j) + SMALL)

        namei = short_spec_names(i)
        if (i==neq) namei = "T"

        ! only dump the ones that don't match
        if (num_jac /= pd(i,j)) then
           if (num_jac /= ZERO) then
              write (*,999) trim(namei), &
                   trim(namej), num_jac, pd(i,j), (num_jac-pd(i,j))/num_jac
           else
              write (*,999) trim(namei), &
                   trim(namej), num_jac, pd(i,j)
           endif
        endif

     enddo
  enddo


end program testjacobian
