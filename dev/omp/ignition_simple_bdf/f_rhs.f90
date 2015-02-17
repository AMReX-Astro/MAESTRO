! My understanding of GPU/OpenACC constructs (subject to ignorance): 
!   vector --> vector of threads, typically a CUDA block.  A block of threads
!   executes in lock-step on data.  gang --> CUDA grid.  Grid is a group of blocks,
!   with each block being independent.  It's possible to get threads within a
!   block to share/communicate, but not between blocks.  As of now we have no
!   intention of making use of the ability to communicate within block

! VBDF-style vectorized RHS with OpenACC acceleration
subroutine f_rhs_vec(neq, npt, y, t, yd_out, upar)

  use bl_types, only: dp_t
  use bl_constants_module, only: ONE
  use network
  use network_indices
  use rpar_indices

  implicit none
  !$acc routine(screenz) seq

  ! our convention is that y(1:nspec) are the species (in the same
  ! order as defined in network.f90, and y(nspec+1) is the temperature
  
  integer,         intent(in   ) :: neq, npt
  real(kind=dp_t), intent(in   ) :: y(neq,npt), t
  real(kind=dp_t), intent(  out) :: yd_out(neq,npt)
  real(kind=dp_t), intent(inout), optional :: upar(:,:)

  integer :: k, n
  real(kind=dp_t) :: ymass(nspec,npt), yd(neq,npt)
  real(kind=dp_t) :: dens(npt), c_p(npt), dhdX(nspec,npt), X_O16(npt)
  real(kind=dp_t) :: temp(npt), T9, T9a, dT9dt, dT9adt

  real(kind=dp_t) :: rate, dratedt
  real(kind=dp_t) :: sc1212, dsc1212dt
  real(kind=dp_t) :: xc12tmp

  real(kind=dp_t), PARAMETER :: &
                     one_twelvth = 1.0d0/12.0d0, &
                     five_sixths = 5.0d0/ 6.0d0, &
                       one_third = 1.0d0/ 3.0d0, &
                      two_thirds = 2.0d0/ 3.0d0

  real(kind=dp_t) :: scratch, dscratchdt
  
  real(kind=dp_t) :: a, b, dadt, dbdt

  ! These are only needed to placate the restrictions of OpenACC
  real(kind=dp_t) :: otw_cpy, fs_cpy, oth_cpy, tt_cpy, o_cpy
  real(kind=dp_t), allocatable :: upar_acc(:,:)
  integer         :: ic12_cpy, img24_cpy, nspec_cpy, nspec_adv_cpy
  integer         :: irp_rate_cpy, irp_dratedt_cpy, irp_sc1212_cpy
  integer         :: irp_dsc1212dt_cpy, irp_xc12tmp_cpy

  otw_cpy   = one_twelvth
  fs_cpy    = five_sixths
  oth_cpy   = one_third
  tt_cpy    = two_thirds
  o_cpy     = ONE
  ic12_cpy  = ic12_
  img24_cpy = img24_
  irp_rate_cpy = irp_rate
  irp_dratedt_cpy = irp_dratedt
  irp_sc1212_cpy = irp_sc1212
  irp_dsc1212dt_cpy = irp_dsc1212dt
  irp_xc12tmp_cpy = irp_xc12tmp
  nspec_cpy = nspec
  nspec_adv_cpy = nspec_advance

  allocate(upar_acc(size(upar,1), size(upar,2)))

  dens(:) = upar(irp_dens,:)
  temp(:) = y(nspec_advance+1,:)

  c_p(:)  = upar(irp_cp,:)
  dhdX(:,:) = upar(irp_dhdX:irp_dhdX-1+nspec,:)
  X_O16(:)   = upar(irp_o16,:)

  ! compute the molar fractions -- needed for the screening
  ymass(ic12_,:) = y(1,:)/aion(ic12_)
  ymass(io16_,:) = X_O16(:)/aion(io16_)
  ymass(img24_,:) = (ONE - y(1,:) - X_O16)/aion(img24_)

 
  !OpenACC guidelines: 
  !  You should copyin read-only variables that will be shared by all
  !  kernels/gangs on the accelerator.  copyout shared write-only variables.
  !  For variables that will be needed for each thread executing a subset of
  !  loop iterations, create() them and declare them private() in the loop
  !  directive.  This will create a copy for each executing thread.
  !$acc data copyin(temp, dens, y, ymass, c_p, ebin, dhdx)                  &
  !$acc   copyin(aion, zion, nspec_cpy, nspec_adv_cpy)                      &
  !$acc   copyin(oth_cpy, fs_cpy, otw_cpy, tt_cpy, o_cpy)                   &
  !$acc   copyin(ic12_cpy,img24_cpy)                                        &
  !$acc   copyin(irp_rate_cpy, irp_dratedt_cpy, irp_sc1212_cpy)             &
  !$acc   copyin(irp_dsc1212dt_cpy, irp_xc12tmp_cpy, npt)                   &
  !$acc   copyout(yd, upar_acc)                                             &
  !$acc   create(sc1212, dsc1212dt, T9, dT9dt, T9a, dT9adt, scratch)    &
  !$acc   create(dscratchdt, a, dadt, b, dbdt, rate, dratedt, xc12tmp)    
 
   
  !$acc kernels default (none)  !default(none) is highly recommended so that
  !  all data must be given a clear scope by the programmer.  
  !  OpenACC/Cray's automatic scoping often leads
  !  to undesired, confusing results or poor performance.
  !$acc loop independent                                                    &
  !$acc   private(sc1212, dsc1212dt, T9, dT9dt, T9a, dT9adt)            &
  !$acc   private(scratch, dscratchdt, a, dadt, b, dbdt)                    &
  !$acc   private(rate, dratedt, xc12tmp, n)    
  do n = 1, npt
    ! call the screening routine
    call screenz(temp(n),dens(n),6.0d0,6.0d0,12.0d0,12.0d0,ymass(:,n),aion,zion,nspec,sc1212,dsc1212dt)
    
    ! compute some often used temperature constants
    T9     = temp(n)/1.d9
    !dT9dt  = ONE/1.d9
    dT9dt  = o_cpy/1.d9
    T9a    = T9/(1.0d0 + 0.0396d0*T9)
    dT9adt = (T9a / T9 - (T9a / (1.0d0 + 0.0396d0*T9)) * 0.0396d0) * dT9dt

    ! compute the CF88 rate
    !scratch    = T9a**one_third
    scratch    = T9a**oth_cpy
    !dscratchdt = one_third * T9a**(-2.0d0 * one_third) * dT9adt
    dscratchdt = oth_cpy * T9a**(-2.0d0 * oth_cpy) * dT9adt

    !a       = 4.27d26*T9a**five_sixths*T9**(-1.5d0)
    a       = 4.27d26*T9a**fs_cpy*T9**(-1.5d0)
    !dadt    = five_sixths * (a/T9a) * dT9adt - 1.5d0 * (a/T9) * dT9dt
    dadt    = fs_cpy * (a/T9a) * dT9adt - 1.5d0 * (a/T9) * dT9dt

    b       = dexp(-84.165d0/scratch - 2.12d-3*T9*T9*T9)
    dbdt    = (84.165d0 * dscratchdt/ scratch**2.0d0                            &
      - 3.0d0 * 2.12d-3 * T9 * T9 * dT9dt) * b

    rate    = a *  b
    dratedt = dadt * b + a * dbdt

    ! The change in number density of C12 is
    ! d(n12)/dt = - 2 * 1/2 (n12)**2 <sigma v>
    !
    ! where <sigma v> is the average of the relative velocity times the cross
    ! section for the reaction, and the factor accounting for the total number
    ! of particle pairs has a 1/2 because we are considering a reaction involving 
    ! identical particles (see Clayton p. 293).  Finally, the -2 means that for
    ! each reaction, we lose 2 carbon nuclei.
    !
    ! The corresponding Mg24 change is
    ! d(n24)/dt = + 1/2 (n12)**2 <sigma v>
    !
    ! note that no factor of 2 appears here, because we create only 1 Mg nuclei.
    !
    ! Switching over to mass fractions, using n = rho X N_A/A, where N_A is
    ! Avagadro's number, and A is the mass number of the nucleon, we get
    !
    ! d(X12)/dt = -2 *1/2 (X12)**2 rho N_A <sigma v> / A12
    !
    ! d(X24)/dt = + 1/2 (X12)**2 rho N_A <sigma v> (A24/A12**2)
    !
    ! these are equal and opposite.
    !
    ! The quantity [N_A <sigma v>] is what is tabulated in Caughlin and Fowler.

    ! we will always refer to the species by integer indices that come from
    ! the network module -- this makes things robust to a shuffling of the 
    ! species ordering

    xc12tmp = max(y(ic12_cpy,n),0.d0)
    yd(ic12_cpy,n) = -otw_cpy*dens(n)*sc1212*rate*xc12tmp**2

    ! now compute the change in temperature, using the evolution equation
    ! dT/dt = -(1/c_p) sum_k (xi_k + q_k) omega_k
    ! 
    ! we make use of the fact that omega(Mg24) = - omega(C12), and that
    ! omega(O16) = 0 in our simplified burner
    !yd(nspec_advance+1,n) =  ( (dhdx(img24_,n) - dhdx(ic12_,n)) + &
    !  (ebin(img24_) - ebin(ic12_)) )*yd(ic12_,n)/c_p(n)
    yd(nspec_adv_cpy+1,n) =  ( (dhdx(img24_cpy,n) - dhdx(ic12_cpy,n)) + &
      (ebin(img24_cpy) - ebin(ic12_cpy)) )*yd(ic12_cpy,n)/c_p(n)

    ! for Mr. Jacobian
    upar_acc(irp_rate_cpy,n)      = rate
    upar_acc(irp_dratedt_cpy,n)   = dratedt
    upar_acc(irp_sc1212_cpy,n)    = sc1212
    upar_acc(irp_dsc1212dt_cpy,n) = dsc1212dt
    upar_acc(irp_xc12tmp_cpy,n)   = xc12tmp
  enddo
  !$acc end kernels
  
  !$acc end data

  
  yd_out(ic12_,:) = yd(ic12_,:)
  yd_out(nspec_advance+1,:) = yd(nspec_advance+1,:)
  upar(irp_rate,:)      = upar_acc(irp_rate,:)     
  upar(irp_dratedt,:)   = upar_acc(irp_dratedt,:)  
  upar(irp_sc1212,:)    = upar_acc(irp_sc1212,:)   
  upar(irp_dsc1212dt,:) = upar_acc(irp_dsc1212dt,:)
  upar(irp_xc12tmp,:)   = upar_acc(irp_xc12tmp,:)  

  !do n = 1, npt
  !  print *, 'yd:  ', yd_out(:,n)
  !  print *, 'up:  ', upar(:,n)
  !enddo
  !stop 'debug'
  return

end subroutine f_rhs_vec

! DVODE-style RHS
subroutine f_rhs(n, t, y, ydot, rpar, ipar)

  use bl_types, only: dp_t
  use bl_constants_module, only: ONE
  use network
  use network_indices
  use rpar_indices

  implicit none

  ! our convention is that y(1:nspec) are the species (in the same
  ! order as defined in network.f90, and y(nspec+1) is the temperature
  integer,         intent(in   ) :: n, ipar(:)
  real(kind=dp_t), intent(in   ) :: y(n), t
  real(kind=dp_t), intent(  out) :: ydot(n)
  real(kind=dp_t), intent(inout) :: rpar(*)  !Works with DVODE
  !real(kind=dp_t), intent(inout) :: rpar(:)   !Works with VBDF's bdf_wrap

  integer :: k
  real(kind=dp_t) :: ymass(nspec)

  real(kind=dp_t) :: dens, c_p, dhdX(nspec), X_O16
  real(kind=dp_t) :: temp, T9, T9a, dT9dt, dT9adt

  real(kind=dp_t) :: rate, dratedt
  real(kind=dp_t) :: sc1212, dsc1212dt
  real(kind=dp_t) :: xc12tmp

  real(kind=dp_t), PARAMETER :: &
                     one_twelvth = 1.0d0/12.0d0, &
                     five_sixths = 5.0d0/ 6.0d0, &
                       one_third = 1.0d0/ 3.0d0, &
                      two_thirds = 2.0d0/ 3.0d0

  real(kind=dp_t) :: scratch, dscratchdt
  
  real(kind=dp_t) :: a, b, dadt, dbdt

  dens = rpar(irp_dens)
  temp = y(nspec_advance+1)

  c_p     = rpar(irp_cp)
  dhdX(:) = rpar(irp_dhdX:irp_dhdX-1+nspec)
  X_O16   = rpar(irp_o16)

  ! compute the molar fractions -- needed for the screening
  ymass(ic12_) = y(1)/aion(ic12_)
  ymass(io16_) = X_O16/aion(io16_)
  ymass(img24_) = (ONE - y(1) - X_O16)/aion(img24_)

  ! call the screening routine
  call screenz(temp,dens,6.0d0,6.0d0,12.0d0,12.0d0,ymass,aion,zion,nspec,     &
    sc1212, dsc1212dt)

  ! compute some often used temperature constants
  T9     = temp/1.d9
  dT9dt  = ONE/1.d9
  T9a    = T9/(1.0d0 + 0.0396d0*T9)
  dT9adt = (T9a / T9 - (T9a / (1.0d0 + 0.0396d0*T9)) * 0.0396d0) * dT9dt

  ! compute the CF88 rate
  scratch    = T9a**one_third
  dscratchdt = one_third * T9a**(-2.0d0 * one_third) * dT9adt

  a       = 4.27d26*T9a**five_sixths*T9**(-1.5d0)
  dadt    = five_sixths * (a/T9a) * dT9adt - 1.5d0 * (a/T9) * dT9dt

  b       = dexp(-84.165d0/scratch - 2.12d-3*T9*T9*T9)
  dbdt    = (84.165d0 * dscratchdt/ scratch**2.0d0                            &
    - 3.0d0 * 2.12d-3 * T9 * T9 * dT9dt) * b

  rate    = a *  b
  dratedt = dadt * b + a * dbdt

  ! The change in number density of C12 is
  ! d(n12)/dt = - 2 * 1/2 (n12)**2 <sigma v>
  !
  ! where <sigma v> is the average of the relative velocity times the cross
  ! section for the reaction, and the factor accounting for the total number
  ! of particle pairs has a 1/2 because we are considering a reaction involving 
  ! identical particles (see Clayton p. 293).  Finally, the -2 means that for
  ! each reaction, we lose 2 carbon nuclei.
  !
  ! The corresponding Mg24 change is
  ! d(n24)/dt = + 1/2 (n12)**2 <sigma v>
  !
  ! note that no factor of 2 appears here, because we create only 1 Mg nuclei.
  !
  ! Switching over to mass fractions, using n = rho X N_A/A, where N_A is
  ! Avagadro's number, and A is the mass number of the nucleon, we get
  !
  ! d(X12)/dt = -2 *1/2 (X12)**2 rho N_A <sigma v> / A12
  !
  ! d(X24)/dt = + 1/2 (X12)**2 rho N_A <sigma v> (A24/A12**2)
  !
  ! these are equal and opposite.
  !
  ! The quantity [N_A <sigma v>] is what is tabulated in Caughlin and Fowler.

  ! we will always refer to the species by integer indices that come from
  ! the network module -- this makes things robust to a shuffling of the 
  ! species ordering

  xc12tmp = max(y(ic12_),0.d0)
  ydot(ic12_) = -one_twelvth*dens*sc1212*rate*xc12tmp**2

  ! now compute the change in temperature, using the evolution equation
  ! dT/dt = -(1/c_p) sum_k (xi_k + q_k) omega_k
  ! 
  ! we make use of the fact that omega(Mg24) = - omega(C12), and that
  ! omega(O16) = 0 in our simplified burner
  ydot(nspec_advance+1) =  ( (dhdx(img24_) - dhdx(ic12_)) + &
    (ebin(img24_) - ebin(ic12_)) )*ydot(ic12_)/c_p


  ! for Mr. Jacobian
  rpar(irp_rate)      = rate
  rpar(irp_dratedt)   = dratedt
  rpar(irp_sc1212)    = sc1212
  rpar(irp_dsc1212dt) = dsc1212dt
  rpar(irp_xc12tmp)   = xc12tmp

  return

end subroutine f_rhs

! VBDF-style vectorized Jacobian
subroutine jac_vec(neq, npt, y, t, pd, upar)
  use bl_types
  use bl_constants_module
  use network
  use network_indices
  use rpar_indices

  implicit none

  integer        , intent(IN   ) :: neq, npt
  real(kind=dp_t), intent(IN   ) :: y(neq,npt), upar(:,:), t 
  real(kind=dp_t), intent(  OUT) :: pd(neq,neq,npt)

  real(kind=dp_t) :: dens(npt), c_p(npt), dhdX(nspec,npt), X_O16(npt)
  real(kind=dp_t) :: rate, dratedt, scorr, dscorrdt, xc12tmp

  integer :: itemp, n

  dens(:)    = upar(irp_dens,:)
  c_p(:)     = upar(irp_cp,:)
  dhdX(:,:)  = upar(irp_dhdX:irp_dhdX-1+nspec,:)
  X_O16(:)   = upar(irp_o16,:)

  do n=1, npt
    rate     = upar(irp_rate,n)     
    dratedt  = upar(irp_dratedt,n)  
    scorr    = upar(irp_sc1212,n)   
    dscorrdt = upar(irp_dsc1212dt,n)
    xc12tmp  = upar(irp_xc12tmp,n)  

    ! initialize
    pd(:,:,n)  = ZERO

    itemp = nspec_advance + 1


    ! carbon jacobian elements
    pd(ic12_, ic12_,n) = -(1.0d0/6.0d0)*dens(n)*scorr*rate*xc12tmp



    ! add the temperature derivatives: df(y_i) / dT
    pd(ic12_,itemp,n) = -(1.0d0/12.0d0)*(dens(n)*rate*xc12tmp**2*dscorrdt    &
                                     + dens(n)*scorr*xc12tmp**2*dratedt)  

    ! add the temperature jacobian elements df(T) / df(y)
    pd(itemp,ic12_,n) =  ( (dhdx(img24_,n) - dhdx(ic12_,n)) + &
                         (ebin(img24_) - ebin(ic12_)) )*pd(ic12_,ic12_,n)/c_p(n)


    ! add df(T) / dT
    pd(itemp,itemp,n) = ( (dhdx(img24_,n) - dhdx(ic12_,n)) + &
                        (ebin(img24_) - ebin(ic12_)) )*pd(ic12_,itemp,n)/c_p(n)
  end do

  return
end subroutine jac_vec

!
!
! DVODE-style Jacobian
subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  use bl_types
  use bl_constants_module
  use network
  use network_indices
  use rpar_indices

  implicit none

  integer        , intent(IN   ) :: neq, ml, mu, nrpd, ipar(:)
  real(kind=dp_t), intent(IN   ) :: y(neq), rpar(*), t      !Works with DVODE
  !real(kind=dp_t), intent(IN   ) :: y(neq), rpar(:), t     !Works with VBDF's bdf_wrap
  real(kind=dp_t), intent(  OUT) :: pd(neq,neq)

  real(kind=dp_t) :: dens, c_p, dhdX(nspec), X_O16
  real(kind=dp_t) :: rate, dratedt, scorr, dscorrdt, xc12tmp

  integer :: itemp

  dens    = rpar(irp_dens)
  c_p     = rpar(irp_cp)
  dhdX(:) = rpar(irp_dhdX:irp_dhdX-1+nspec)
  X_O16   = rpar(irp_o16)

  rate     = rpar(irp_rate)     
  dratedt  = rpar(irp_dratedt)  
  scorr    = rpar(irp_sc1212)   
  dscorrdt = rpar(irp_dsc1212dt)
  xc12tmp  = rpar(irp_xc12tmp)  

  ! initialize
  pd(:,:)  = ZERO

  itemp = nspec_advance + 1


  ! carbon jacobian elements
  pd(ic12_, ic12_) = -(1.0d0/6.0d0)*dens*scorr*rate*xc12tmp



  ! add the temperature derivatives: df(y_i) / dT
  pd(ic12_,itemp) = -(1.0d0/12.0d0)*(dens*rate*xc12tmp**2*dscorrdt    &
                                   + dens*scorr*xc12tmp**2*dratedt)  

  ! add the temperature jacobian elements df(T) / df(y)
  pd(itemp,ic12_) =  ( (dhdx(img24_) - dhdx(ic12_)) + &
                       (ebin(img24_) - ebin(ic12_)) )*pd(ic12_,ic12_)/c_p


  ! add df(T) / dT
  pd(itemp,itemp) = ( (dhdx(img24_) - dhdx(ic12_)) + &
                      (ebin(img24_) - ebin(ic12_)) )*pd(ic12_,itemp)/c_p


  return
end subroutine jac

