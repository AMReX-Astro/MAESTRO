module feval
   use bl_types, only: dp_t
   use bl_constants_module, only: ZERO, ONE, TWELFTH, FIVE6TH, THIRD
   use network, only: nspec, nspec_advance, aion, zion, ebin
   use network_indices, only: ic12_, io16_, img24_
   use rpar_indices, only: irp_dens, irp_cp, irp_dhdx, irp_o16, irp_rate, &
                           irp_dratedt, irp_sc1212, irp_dsc1212dt, irp_xc12tmp
   implicit none
contains
   subroutine f_rhs_vec(neq, npt, y, t, yd, upar)
      !$acc routine seq
      !$acc routine(screenz) seq
   
      ! our convention is that y(1:nspec) are the species (in the same
      ! order as defined in network.f90, and y(nspec+1) is the temperature
      
      integer,         intent(in   ) :: neq, npt
      real(kind=dp_t), intent(in   ) :: y(neq,npt), t
      !real(kind=dp_t), intent(  out) :: yd_out(neq,npt)
      real(kind=dp_t), intent(  out) :: yd(neq,npt)
      real(kind=dp_t), intent(inout) :: upar(:,:)
   
      integer :: k, n
      real(kind=dp_t) :: ymass(nspec,npt)!, yd(neq,npt)
      real(kind=dp_t) :: dens(npt), c_p(npt), dhdX(nspec,npt), X_O16(npt)
      real(kind=dp_t) :: temp(npt), T9, T9a, dT9dt, dT9adt
   
      real(kind=dp_t) :: rate, dratedt
      real(kind=dp_t) :: sc1212, dsc1212dt
      real(kind=dp_t) :: xc12tmp
   
      !real(kind=dp_t), PARAMETER :: &
      !                   one_twelvth = 1.0d0/12.0d0, &
      !                   five_sixths = 5.0d0/ 6.0d0, &
      !                     one_third = 1.0d0/ 3.0d0, &
      !                    two_thirds = 2.0d0/ 3.0d0
   
      !                   one_twelvth = 1.0d0/12.0d0, & --> TWELFTH
      !                   five_sixths = 5.0d0/ 6.0d0, & --> FIVE6TH
      !                     one_third = 1.0d0/ 3.0d0, & --> THIRD
      !                    two_thirds = 2.0d0/ 3.0d0    --> TWO3RD
      real(kind=dp_t) :: scratch, dscratchdt
      
      real(kind=dp_t) :: a, b, dadt, dbdt
  
       
    
      do n = 1, npt
         !initialize data
         dens(n) = upar(irp_dens,n)
         temp(n) = y(nspec_advance+1,n)
   
         c_p(n)  = upar(irp_cp,n)
         do k = 1, nspec
            dhdX(k,n) = upar(irp_dhdX - 1 + k, n)
         end do
         X_O16(n)   = upar(irp_o16,n)
   
         ! compute the molar fractions -- needed for the screening
         ymass(ic12_,n) = y(1,n)/aion(ic12_)
         ymass(io16_,n) = X_O16(n)/aion(io16_)
         ymass(img24_,n) = (ONE - y(1,n) - X_O16(n))/aion(img24_)

        ! call the screening routine
        call screenz(temp(n),dens(n),6.0d0,6.0d0,12.0d0,12.0d0,ymass(:,n),aion,zion,nspec,sc1212,dsc1212dt)
        
        ! compute some often used temperature constants
        T9     = temp(n)/1.d9
        dT9dt  = ONE/1.d9
        T9a    = T9/(ONE + 0.0396d0*T9)
        dT9adt = (T9a / T9 - (T9a / (ONE + 0.0396d0*T9)) * 0.0396d0) * dT9dt
   
        ! compute the CF88 rate
        scratch    = T9a**THIRD
        dscratchdt = THIRD * T9a**(-2.0d0 * THIRD) * dT9adt
   
        a       = 4.27d26*T9a**FIVE6TH*T9**(-1.5d0)
        dadt    = FIVE6TH * (a/T9a) * dT9adt - 1.5d0 * (a/T9) * dT9dt
   
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
   
        xc12tmp = max(y(ic12_,n),0.d0)
        yd(ic12_,n) = -TWELFTH*dens(n)*sc1212*rate*xc12tmp**2
   
        ! now compute the change in temperature, using the evolution equation
        ! dT/dt = -(1/c_p) sum_k (xi_k + q_k) omega_k
        ! 
        ! we make use of the fact that omega(Mg24) = - omega(C12), and that
        ! omega(O16) = 0 in our simplified burner
        yd(nspec_advance+1,n) =  ( (dhdx(img24_,n) - dhdx(ic12_,n)) + &
          (ebin(img24_) - ebin(ic12_)) )*yd(ic12_,n)/c_p(n)
   
        ! for Mr. Jacobian
        upar(irp_rate,n)      = rate
        upar(irp_dratedt,n)   = dratedt
        upar(irp_sc1212,n)    = sc1212
        upar(irp_dsc1212dt,n) = dsc1212dt
        upar(irp_xc12tmp,n)   = xc12tmp
      enddo
      
      !yd_out(ic12_,:) = yd(ic12_,:)
      !yd_out(nspec_advance+1,:) = yd(nspec_advance+1,:)
   
   end subroutine f_rhs_vec

   subroutine jac_vec(neq, npt, y, t, pd, upar)
      !$acc routine seq
   
      integer        , intent(in   ) :: neq, npt
      real(kind=dp_t), intent(in   ) :: y(neq,npt), t 
      real(kind=dp_t), intent(  out) :: pd(neq,neq,npt)
      real(kind=dp_t), intent(inout) :: upar(:,:)
   
      real(kind=dp_t) :: dens(npt), c_p(npt), dhdX(nspec,npt), X_O16(npt)
      real(kind=dp_t) :: rate, dratedt, scorr, dscorrdt, xc12tmp
   
      integer :: itemp, n, k, nn, kk
   
   
      do n=1, npt

        dens(n)    = upar(irp_dens,n)
        c_p(n)     = upar(irp_cp,n)
        do k = 1, nspec
           dhdX(k,n)  = upar(irp_dhdX - 1 + k, n)
        end do
        X_O16(n)   = upar(irp_o16,n)

        rate     = upar(irp_rate,n)     
        dratedt  = upar(irp_dratedt,n)  
        scorr    = upar(irp_sc1212,n)   
        dscorrdt = upar(irp_dsc1212dt,n)
        xc12tmp  = upar(irp_xc12tmp,n)  
   
        ! initialize
        do nn = 1, neq
           do kk = 1, neq
              pd(nn,kk,n) = ZERO
           end do
        end do
   
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
   end subroutine jac_vec
end module feval


!TODO: The below is old, delete it after done using it as development reference
!! VBDF-style vectorized RHS
!subroutine f_rhs_vec(neq, npt, y, t, yd_out, upar)
!  !$acc routine seq
!
!  use bl_types, only: dp_t
!  use bl_constants_module, only: ONE
!  use network
!  use network_indices
!  use rpar_indices
!
!  implicit none
!  !$acc routine(screenz) seq
!
!  ! our convention is that y(1:nspec) are the species (in the same
!  ! order as defined in network.f90, and y(nspec+1) is the temperature
!  
!  integer,         intent(in   ) :: neq, npt
!  real(kind=dp_t), intent(in   ) :: y(neq,npt), t
!  real(kind=dp_t), intent(  out) :: yd_out(neq,npt)
!  real(kind=dp_t), intent(inout) :: upar(:,:)
!
!  integer :: k, n
!  real(kind=dp_t) :: ymass(nspec,npt), yd(neq,npt)
!  real(kind=dp_t) :: dens(npt), c_p(npt), dhdX(nspec,npt), X_O16(npt)
!  real(kind=dp_t) :: temp(npt), T9, T9a, dT9dt, dT9adt
!
!  real(kind=dp_t) :: rate, dratedt
!  real(kind=dp_t) :: sc1212, dsc1212dt
!  real(kind=dp_t) :: xc12tmp
!
!  real(kind=dp_t), PARAMETER :: &
!                     one_twelvth = 1.0d0/12.0d0, &
!                     five_sixths = 5.0d0/ 6.0d0, &
!                       one_third = 1.0d0/ 3.0d0, &
!                      two_thirds = 2.0d0/ 3.0d0
!
!  real(kind=dp_t) :: scratch, dscratchdt
!  
!  real(kind=dp_t) :: a, b, dadt, dbdt
!
!  dens(:) = upar(irp_dens,:)
!  temp(:) = y(nspec_advance+1,:)
!
!  c_p(:)  = upar(irp_cp,:)
!  dhdX(:,:) = upar(irp_dhdX:irp_dhdX-1+nspec,:)
!  X_O16(:)   = upar(irp_o16,:)
!
!  ! compute the molar fractions -- needed for the screening
!  ymass(ic12_,:) = y(1,:)/aion(ic12_)
!  ymass(io16_,:) = X_O16(:)/aion(io16_)
!  ymass(img24_,:) = (ONE - y(1,:) - X_O16)/aion(img24_)
! 
!  do n = 1, npt
!    ! call the screening routine
!    call screenz(temp(n),dens(n),6.0d0,6.0d0,12.0d0,12.0d0,ymass(:,n),aion,zion,nspec,sc1212,dsc1212dt)
!    
!    ! compute some often used temperature constants
!    T9     = temp(n)/1.d9
!    dT9dt  = ONE/1.d9
!    T9a    = T9/(1.0d0 + 0.0396d0*T9)
!    dT9adt = (T9a / T9 - (T9a / (1.0d0 + 0.0396d0*T9)) * 0.0396d0) * dT9dt
!
!    ! compute the CF88 rate
!    scratch    = T9a**one_third
!    dscratchdt = one_third * T9a**(-2.0d0 * one_third) * dT9adt
!
!    a       = 4.27d26*T9a**five_sixths*T9**(-1.5d0)
!    dadt    = five_sixths * (a/T9a) * dT9adt - 1.5d0 * (a/T9) * dT9dt
!
!    b       = dexp(-84.165d0/scratch - 2.12d-3*T9*T9*T9)
!    dbdt    = (84.165d0 * dscratchdt/ scratch**2.0d0                            &
!      - 3.0d0 * 2.12d-3 * T9 * T9 * dT9dt) * b
!
!    rate    = a *  b
!    dratedt = dadt * b + a * dbdt
!
!    ! The change in number density of C12 is
!    ! d(n12)/dt = - 2 * 1/2 (n12)**2 <sigma v>
!    !
!    ! where <sigma v> is the average of the relative velocity times the cross
!    ! section for the reaction, and the factor accounting for the total number
!    ! of particle pairs has a 1/2 because we are considering a reaction involving 
!    ! identical particles (see Clayton p. 293).  Finally, the -2 means that for
!    ! each reaction, we lose 2 carbon nuclei.
!    !
!    ! The corresponding Mg24 change is
!    ! d(n24)/dt = + 1/2 (n12)**2 <sigma v>
!    !
!    ! note that no factor of 2 appears here, because we create only 1 Mg nuclei.
!    !
!    ! Switching over to mass fractions, using n = rho X N_A/A, where N_A is
!    ! Avagadro's number, and A is the mass number of the nucleon, we get
!    !
!    ! d(X12)/dt = -2 *1/2 (X12)**2 rho N_A <sigma v> / A12
!    !
!    ! d(X24)/dt = + 1/2 (X12)**2 rho N_A <sigma v> (A24/A12**2)
!    !
!    ! these are equal and opposite.
!    !
!    ! The quantity [N_A <sigma v>] is what is tabulated in Caughlin and Fowler.
!
!    ! we will always refer to the species by integer indices that come from
!    ! the network module -- this makes things robust to a shuffling of the 
!    ! species ordering
!
!    xc12tmp = max(y(ic12_,n),0.d0)
!    yd(ic12_,n) = -one_twelvth*dens(n)*sc1212*rate*xc12tmp**2
!
!    ! now compute the change in temperature, using the evolution equation
!    ! dT/dt = -(1/c_p) sum_k (xi_k + q_k) omega_k
!    ! 
!    ! we make use of the fact that omega(Mg24) = - omega(C12), and that
!    ! omega(O16) = 0 in our simplified burner
!    yd(nspec_advance+1,n) =  ( (dhdx(img24_,n) - dhdx(ic12_,n)) + &
!      (ebin(img24_) - ebin(ic12_)) )*yd(ic12_,n)/c_p(n)
!
!    ! for Mr. Jacobian
!    upar(irp_rate,n)      = rate
!    upar(irp_dratedt,n)   = dratedt
!    upar(irp_sc1212,n)    = sc1212
!    upar(irp_dsc1212dt,n) = dsc1212dt
!    upar(irp_xc12tmp,n)   = xc12tmp
!  enddo
!
!  
!  yd_out(ic12_,:) = yd(ic12_,:)
!  yd_out(nspec_advance+1,:) = yd(nspec_advance+1,:)
!
!  return
!
!end subroutine f_rhs_vec
!
!! DVODE-style RHS
!subroutine f_rhs(n, t, y, ydot, rpar, ipar)
!
!  use bl_types, only: dp_t
!  use bl_constants_module, only: ONE
!  use network
!  use network_indices
!  use rpar_indices
!
!  implicit none
!
!  ! our convention is that y(1:nspec) are the species (in the same
!  ! order as defined in network.f90, and y(nspec+1) is the temperature
!  integer,         intent(in   ) :: n, ipar(:)
!  real(kind=dp_t), intent(in   ) :: y(n), t
!  real(kind=dp_t), intent(  out) :: ydot(n)
!  real(kind=dp_t), intent(inout) :: rpar(*)  !For DVODE compatibility, use (*) instead of (:)
!
!  integer :: k
!  real(kind=dp_t) :: ymass(nspec)
!
!  real(kind=dp_t) :: dens, c_p, dhdX(nspec), X_O16
!  real(kind=dp_t) :: temp, T9, T9a, dT9dt, dT9adt
!
!  real(kind=dp_t) :: rate, dratedt
!  real(kind=dp_t) :: sc1212, dsc1212dt
!  real(kind=dp_t) :: xc12tmp
!
!  real(kind=dp_t), PARAMETER :: &
!                     one_twelvth = 1.0d0/12.0d0, &
!                     five_sixths = 5.0d0/ 6.0d0, &
!                       one_third = 1.0d0/ 3.0d0, &
!                      two_thirds = 2.0d0/ 3.0d0
!
!  real(kind=dp_t) :: scratch, dscratchdt
!  
!  real(kind=dp_t) :: a, b, dadt, dbdt
!
!  dens = rpar(irp_dens)
!  temp = y(nspec_advance+1)
!
!  c_p     = rpar(irp_cp)
!  dhdX(:) = rpar(irp_dhdX:irp_dhdX-1+nspec)
!  X_O16   = rpar(irp_o16)
!
!  ! compute the molar fractions -- needed for the screening
!  ymass(ic12_) = y(1)/aion(ic12_)
!  ymass(io16_) = X_O16/aion(io16_)
!  ymass(img24_) = (ONE - y(1) - X_O16)/aion(img24_)
!
!  ! call the screening routine
!  call screenz(temp,dens,6.0d0,6.0d0,12.0d0,12.0d0,ymass,aion,zion,nspec,     &
!    sc1212, dsc1212dt)
!
!  ! compute some often used temperature constants
!  T9     = temp/1.d9
!  dT9dt  = ONE/1.d9
!  T9a    = T9/(1.0d0 + 0.0396d0*T9)
!  dT9adt = (T9a / T9 - (T9a / (1.0d0 + 0.0396d0*T9)) * 0.0396d0) * dT9dt
!
!  ! compute the CF88 rate
!  scratch    = T9a**one_third
!  dscratchdt = one_third * T9a**(-2.0d0 * one_third) * dT9adt
!
!  a       = 4.27d26*T9a**five_sixths*T9**(-1.5d0)
!  dadt    = five_sixths * (a/T9a) * dT9adt - 1.5d0 * (a/T9) * dT9dt
!
!  b       = dexp(-84.165d0/scratch - 2.12d-3*T9*T9*T9)
!  dbdt    = (84.165d0 * dscratchdt/ scratch**2.0d0                            &
!    - 3.0d0 * 2.12d-3 * T9 * T9 * dT9dt) * b
!
!  rate    = a *  b
!  dratedt = dadt * b + a * dbdt
!
!  ! The change in number density of C12 is
!  ! d(n12)/dt = - 2 * 1/2 (n12)**2 <sigma v>
!  !
!  ! where <sigma v> is the average of the relative velocity times the cross
!  ! section for the reaction, and the factor accounting for the total number
!  ! of particle pairs has a 1/2 because we are considering a reaction involving 
!  ! identical particles (see Clayton p. 293).  Finally, the -2 means that for
!  ! each reaction, we lose 2 carbon nuclei.
!  !
!  ! The corresponding Mg24 change is
!  ! d(n24)/dt = + 1/2 (n12)**2 <sigma v>
!  !
!  ! note that no factor of 2 appears here, because we create only 1 Mg nuclei.
!  !
!  ! Switching over to mass fractions, using n = rho X N_A/A, where N_A is
!  ! Avagadro's number, and A is the mass number of the nucleon, we get
!  !
!  ! d(X12)/dt = -2 *1/2 (X12)**2 rho N_A <sigma v> / A12
!  !
!  ! d(X24)/dt = + 1/2 (X12)**2 rho N_A <sigma v> (A24/A12**2)
!  !
!  ! these are equal and opposite.
!  !
!  ! The quantity [N_A <sigma v>] is what is tabulated in Caughlin and Fowler.
!
!  ! we will always refer to the species by integer indices that come from
!  ! the network module -- this makes things robust to a shuffling of the 
!  ! species ordering
!
!  xc12tmp = max(y(ic12_),0.d0)
!  ydot(ic12_) = -one_twelvth*dens*sc1212*rate*xc12tmp**2
!
!  ! now compute the change in temperature, using the evolution equation
!  ! dT/dt = -(1/c_p) sum_k (xi_k + q_k) omega_k
!  ! 
!  ! we make use of the fact that omega(Mg24) = - omega(C12), and that
!  ! omega(O16) = 0 in our simplified burner
!  ydot(nspec_advance+1) =  ( (dhdx(img24_) - dhdx(ic12_)) + &
!    (ebin(img24_) - ebin(ic12_)) )*ydot(ic12_)/c_p
!
!
!  ! for Mr. Jacobian
!  rpar(irp_rate)      = rate
!  rpar(irp_dratedt)   = dratedt
!  rpar(irp_sc1212)    = sc1212
!  rpar(irp_dsc1212dt) = dsc1212dt
!  rpar(irp_xc12tmp)   = xc12tmp
!
!  return
!
!end subroutine f_rhs
!
!! VBDF-style vectorized Jacobian
!subroutine jac_vec(neq, npt, y, t, pd, upar)
!  !$acc routine seq
!  use bl_types
!  use bl_constants_module
!  use network
!  use network_indices
!  use rpar_indices
!
!  implicit none
!
!  integer        , intent(IN   ) :: neq, npt
!  real(kind=dp_t), intent(IN   ) :: y(neq,npt), t 
!  real(kind=dp_t), intent(  OUT) :: pd(neq,neq,npt)
!  real(dp_t),      intent(inout) :: upar(:,:)
!
!  real(kind=dp_t) :: dens(npt), c_p(npt), dhdX(nspec,npt), X_O16(npt)
!  real(kind=dp_t) :: rate, dratedt, scorr, dscorrdt, xc12tmp
!
!  integer :: itemp, n
!
!  dens(:)    = upar(irp_dens,:)
!  c_p(:)     = upar(irp_cp,:)
!  dhdX(:,:)  = upar(irp_dhdX:irp_dhdX-1+nspec,:)
!  X_O16(:)   = upar(irp_o16,:)
!
!  do n=1, npt
!    rate     = upar(irp_rate,n)     
!    dratedt  = upar(irp_dratedt,n)  
!    scorr    = upar(irp_sc1212,n)   
!    dscorrdt = upar(irp_dsc1212dt,n)
!    xc12tmp  = upar(irp_xc12tmp,n)  
!
!    ! initialize
!    pd(:,:,n)  = ZERO
!
!    itemp = nspec_advance + 1
!
!
!    ! carbon jacobian elements
!    pd(ic12_, ic12_,n) = -(1.0d0/6.0d0)*dens(n)*scorr*rate*xc12tmp
!
!
!    ! add the temperature derivatives: df(y_i) / dT
!    pd(ic12_,itemp,n) = -(1.0d0/12.0d0)*(dens(n)*rate*xc12tmp**2*dscorrdt    &
!                                     + dens(n)*scorr*xc12tmp**2*dratedt)  
!
!    ! add the temperature jacobian elements df(T) / df(y)
!    pd(itemp,ic12_,n) =  ( (dhdx(img24_,n) - dhdx(ic12_,n)) + &
!                         (ebin(img24_) - ebin(ic12_)) )*pd(ic12_,ic12_,n)/c_p(n)
!
!
!    ! add df(T) / dT
!    pd(itemp,itemp,n) = ( (dhdx(img24_,n) - dhdx(ic12_,n)) + &
!                        (ebin(img24_) - ebin(ic12_)) )*pd(ic12_,itemp,n)/c_p(n)
!  end do
!
!  return
!end subroutine jac_vec
!
!!
!!
!! DVODE-style Jacobian
!subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
!
!  use bl_types
!  use bl_constants_module
!  use network
!  use network_indices
!  use rpar_indices
!
!  implicit none
!
!  integer        , intent(IN   ) :: neq, ml, mu, nrpd, ipar(:)
!  real(kind=dp_t), intent(IN   ) :: y(neq), rpar(*), t      !For DVODE compatibility, use (*) instead of (:)
!  real(kind=dp_t), intent(  OUT) :: pd(neq,neq)
!
!  real(kind=dp_t) :: dens, c_p, dhdX(nspec), X_O16
!  real(kind=dp_t) :: rate, dratedt, scorr, dscorrdt, xc12tmp
!
!  integer :: itemp
!
!  dens    = rpar(irp_dens)
!  c_p     = rpar(irp_cp)
!  dhdX(:) = rpar(irp_dhdX:irp_dhdX-1+nspec)
!  X_O16   = rpar(irp_o16)
!
!  rate     = rpar(irp_rate)     
!  dratedt  = rpar(irp_dratedt)  
!  scorr    = rpar(irp_sc1212)   
!  dscorrdt = rpar(irp_dsc1212dt)
!  xc12tmp  = rpar(irp_xc12tmp)  
!
!  ! initialize
!  pd(:,:)  = ZERO
!
!  itemp = nspec_advance + 1
!
!
!  ! carbon jacobian elements
!  pd(ic12_, ic12_) = -(1.0d0/6.0d0)*dens*scorr*rate*xc12tmp
!
!
!
!  ! add the temperature derivatives: df(y_i) / dT
!  pd(ic12_,itemp) = -(1.0d0/12.0d0)*(dens*rate*xc12tmp**2*dscorrdt    &
!                                   + dens*scorr*xc12tmp**2*dratedt)  
!
!  ! add the temperature jacobian elements df(T) / df(y)
!  pd(itemp,ic12_) =  ( (dhdx(img24_) - dhdx(ic12_)) + &
!                       (ebin(img24_) - ebin(ic12_)) )*pd(ic12_,ic12_)/c_p
!
!
!  ! add df(T) / dT
!  pd(itemp,itemp) = ( (dhdx(img24_) - dhdx(ic12_)) + &
!                      (ebin(img24_) - ebin(ic12_)) )*pd(ic12_,itemp)/c_p
!
!
!  return
!end subroutine jac
!
