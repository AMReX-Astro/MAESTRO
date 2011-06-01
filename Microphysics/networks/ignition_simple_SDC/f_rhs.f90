subroutine f_rhs(n, t, y, ydot, rpar, ipar)

  use bl_types
  use bl_constants_module
  use network
  use eos_module
  
  use burner_aux_module, only : sdc_rho_pass, sdc_X_pass, p0_pass

  implicit none

  ! our convention is that y(1:nspec) are the species (in the same
  ! order as defined in network.f90, and y(nspec+1) is the density
  integer :: n
  real(kind=dp_t) :: y(n), ydot(n)

  integer :: k
  real(kind=dp_t) :: ymass(nspec)

  real(kind=dp_t) :: rpar
  integer :: ipar

  real(kind=dp_t) :: t

  real(kind=dp_t) :: dens, temp, T9, T9a, dT9dt, dT9adt

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

  integer, save :: ic12, io16, img24

  logical, save :: firstCall = .true.

  if (firstCall) then
     ic12 = network_species_index("carbon-12")
     io16 = network_species_index("oxygen-16")
     img24 = network_species_index("magnesium-24")

     firstCall = .false.
  end if



  ! compute the molar fractions -- needed for the screening
  ymass(ic12) = y(1)/aion(ic12)
  ymass(io16) = y(2)/aion(io16)
  ymass(img24) = y(3)/aion(img24)

  dens = y(nspec+1)

  ! compute temp with EOS -- note, here we are assuming use_tfromp = T
  den_eos(1) = dens
  xn_eos(1,:) = y(1:nspec)
  p_eos(1) = p0_pass

  call eos(eos_input_rp, den_eos, temp_eos, &
           npts, &
           xn_eos, &
           p_eos, h_eos, e_eos, &
           cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
           dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
           dpdX_eos, dhdX_eos, &
           gam1_eos, cs_eos, s_eos, &
           dsdt_eos, dsdr_eos, &
           .false., &
           pt_index_eos)

  temp = temp_eos(1)

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

  xc12tmp = max(y(ic12),0.d0)
  ydot(1) = -one_twelvth*dens*sc1212*rate*xc12tmp**2 + sdc_X_pass(1)
  ydot(2) = sdc_X_pass(2)
  ydot(3) = one_twelvth*dens*sc1212*rate*xc12tmp**2 + sdc_X_pass(3)

  ydot(nspec+1) =  sdc_rho_pass

  return

end subroutine f_rhs


subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)

  use bl_types
  use bl_constants_module
  use network

  ! we get the thermodynamic state through the burner_aux module -- we freeze
  ! these to the values are the top of the timestep to avoid costly
  ! EOS calls

  implicit none

  integer        , intent(IN   ) :: neq, ml, mu, nrpd, ipar
  real(kind=dp_t), intent(IN   ) :: y(neq), rpar, t
  real(kind=dp_t), intent(  OUT) :: pd(neq,neq)

  ! initialize
  pd(:,:)  = ZERO

  return
end subroutine jac

