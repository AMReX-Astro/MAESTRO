module feval
   use bl_types, only: dp_t
   use bl_constants_module, only: ZERO, ONE, TWELFTH, FIVE6TH, THIRD
   use network
   !use network_indices, only: ic12_, io16_, img24_
   use rpar_indices

   implicit none
  
   private
   public :: f_rhs_vec, jac_vec
   integer :: PT_INDEX = 1
contains
 


   ! The f_rhs routine provides the right-hand-side for the DVODE solver.
   ! It deals with molar abundances throughout (we expect that the input
   ! vector y has molar abundances, Y = X/A) for make_rates, and
   ! dydt routines.  It also checks to see if the temperature has changed
   ! much since the last call - if so, it updates the temperature to get
   ! a better estimate of the reaction rates.
   !
   ! The jac routine provides an explicit Jacobian to the DVODE solver.
   !
   
   !subroutine f_rhs(n, t, y, ydot, rpar, ipar)
   subroutine f_rhs_vec(y, t, ydot, rpar)
   !$acc routine seq
     use bl_types
     use bl_constants_module
     use network
     use eos_module
     use eos_type_module
     use rates_module
     use rpar_indices
   
     implicit none
   
     !integer,         intent(IN   ) :: n, ipar
     !real(kind=dp_t), intent(IN   ) :: t, y(n)
     real(kind=dp_t), intent(IN   ) :: t, y(neqs,burn_npts)
     real(kind=dp_t), intent(INOUT) :: rpar(n_rpar_comps,burn_npts)
     !real(kind=dp_t), intent(  OUT) :: ydot(n)
     real(kind=dp_t), intent(  OUT) :: ydot(neqs,burn_npts)
   
     integer :: k, n
   
     real(kind=dp_t), parameter :: T2T9 = 1.0e-9_dp_t
   
     real(kind=dp_t) :: dens, cp, dhdx(nspec), T9_eos, dT_crit, t9
     real(kind=dp_t) :: ymol(nspec), yd_sum
   
     !type(eos_t) :: eos_state
     
     ydot = ZERO
   
     ! several thermodynamic quantities come in via rpar
     dens = rpar(irp_dens,PT_INDEX)
     T9_eos = rpar(irp_T9_eos,PT_INDEX)
     dT_crit = rpar(irp_dTcrit,PT_INDEX)
  
     do n = 1, nspec
        ymol(n) = y(n,PT_INDEX)
     enddo
     t9 = y(nspec+1,PT_INDEX)
   
     !if (abs(t9 - T9_eos) > dT_crit*T9_eos) then
     !   T9_eos = t9
   
     !   eos_state%T = T9_eos/T2T9
     !   eos_state%rho = dens
     !   eos_state%xn = ymol*aion
   
     !   call eos(eos_input_rt, eos_state)
   
     !   rpar(irp_T9_eos) = T9_eos
     !   rpar(irp_cp) = eos_state%cp
     !   rpar(irp_dhdX:irp_dhdX+nspec-1) = eos_state%dhdX
     !endif
   
     ! more thermodyanmics that possibly were updated
     cp   = rpar(irp_cp,PT_INDEX)
     
     !dhdX = rpar(irp_dhdX:irp_dhdX+nspec-1)
     k=1
     do n=irp_dhdX, irp_dhdX+nspec-1
        dhdX(k) = rpar(n,PT_INDEX)
        k=k+1
     enddo
       
     ! build the rates; weak rates are the wk* variables
     call make_rates(t9, dens, ymol, rpar)
   
     ! set up the ODEs
     call make_ydots(ymol,t9,rpar,ydot(1:nspec,PT_INDEX),.false.)
   
   !.... t9
   ! shouldn't an aion be here?  putting it in breaks VODE convergence...
     !ydot(n) = -sum((dhdX+ebin)*ydot(1:nspec))/cp
     !ydot(n) = ydot(n) * T2T9
     yd_sum = ZERO
     do n=1, nspec
        yd_sum  = yd_sum + (dhdX(n)+ebin(n))*ydot(n,PT_INDEX)
     enddo
     ydot(neqs,PT_INDEX) = -yd_sum/cp * T2T9
   
     return
   
   !end subroutine f_rhs
   end subroutine f_rhs_vec
   
   subroutine make_rates(t9,dens,y,rpar)
   !$acc routine seq
   
     use bl_types
     use bl_constants_module
     use rates_module
     use network
     use rpar_indices
     
     implicit none
   
     real(kind=dp_t), intent(in   ) :: t9, dens, y(nspec,burn_npts)
     real(kind=dp_t), intent(inout) :: rpar(n_rpar_comps,burn_npts)
   
     ! locally used rates
     real(kind=dp_t) :: rate,dratedt,wk18ne,wk19ne
     real(kind=dp_t) :: r56pg,dr56pgdt,cutoni,dcutonidt,r57decay,r56eff,dr56effdt
     real(kind=dp_t) :: t9i32
   
     ! some numbers from appendix C in WW81; these should probably be
     ! updated with current rates
     real(kind=dp_t), parameter :: Lweak = 1.05d0, & ! this is for NS
   !                                Lweak = 0.107d0, & ! this is for lower 
                                                       ! densities
                                   la2 = ONE/FIFTEEN ! mean rate from 30s to 56ni
                                                     ! from p-capture and beta 
                                                     ! decays
                                   
     type (temp_t) :: tfactors
     integer :: n
    
     !rpar(irp_rates:irp_drtdt-1+nrat) = ZERO ! rates and dratesdt
     do n=irp_rates, irp_drtdt-1+nrat 
        rpar(n,PT_INDEX) = ZERO        ! rates and dratesdt
     enddo

     !rpar(irp_dlambCNOdh1:n_rpar_comps) = ZERO ! other constants
     do n=irp_dlambCNOdh1, n_rpar_comps
        rpar(n,PT_INDEX) = ZERO
     enddo
   
     tfactors = calc_tfactors(t9)
   
     ! some common parameters
     rpar(irp_rates-1+irLweak,PT_INDEX) = Lweak
     rpar(irp_rates-1+irla2,PT_INDEX)   = la2
   
     ! weak rates first
     !
     ! 14o(beta nu)14n
     call rate_o14_to_n14(tfactors,rate,dratedt)
     rpar(irp_rates-1+irwk14o,PT_INDEX) = rate
     ! 15o(beta nu)15n
     call rate_o15_to_n15(tfactors,rate,dratedt)
     rpar(irp_rates-1+irwk15o,PT_INDEX) = rate
     ! 17f(beta nu)17o
     call rate_f17_to_o17(tfactors,rate,dratedt)
     rpar(irp_rates-1+irwk17f,PT_INDEX) = rate
     ! these weak rates aren't needed outside of this routine
     ! 18ne(beta nu)18f
     call rate_ne18_to_f18(tfactors,wk18ne,dratedt) 
     ! 19ne(beta nu)19f
     call rate_ne19_to_f19(tfactors,wk19ne,dratedt)
   
     ! 12c(p,g)13n
     call rate_p_c12_to_n13(tfactors,rate,dratedt)
     rpar(irp_rates-1+irpg12c,PT_INDEX) = dens*rate
     rpar(irp_drtdt-1+irpg12c,PT_INDEX) = dens*dratedt
   
     ! triple alpha
     call rate_he4_he4_he4_to_c12(tfactors,rate,dratedt)
     rpar(irp_rates-1+ir3a,PT_INDEX) = dens*dens*rate
     rpar(irp_drtdt-1+ir3a,PT_INDEX) = dens*dens*dratedt
   
     ! 17f(p,g)18ne
     call rate_p_f17_to_ne18(tfactors,rate,dratedt)
     rpar(irp_rates-1+irpg17f,PT_INDEX) = dens*rate
     rpar(irp_drtdt-1+irpg17f,PT_INDEX) = dens*dratedt
   
     ! 17f(g,p)16o
     call rate_f17_to_p_o16(tfactors,rate,dratedt)
     rpar(irp_rates-1+irgp17f,PT_INDEX) = rate
     rpar(irp_drtdt-1+irgp17f,PT_INDEX) = dratedt
   
     ! 15o(a,g)19ne
     call rate_he4_o15_to_ne19(tfactors,rate,dratedt)
     rpar(irp_rates-1+irag15o,PT_INDEX) = dens*rate
     rpar(irp_drtdt-1+irag15o,PT_INDEX) = dens*dratedt
   
     ! 16o(a,g)20ne
     call rate_he4_o16_to_ne20(tfactors,rate,dratedt)
     rpar(irp_rates-1+irag16o,PT_INDEX) = dens*rate
     rpar(irp_drtdt-1+irag16o,PT_INDEX) = dens*dratedt
   
     ! 16o(p,g)17f
     call rate_p_o16_to_f17(tfactors,rate,dratedt)
     rpar(irp_rates-1+irpg16o,PT_INDEX) = dens*rate
     rpar(irp_drtdt-1+irpg16o,PT_INDEX) = dens*dratedt
   
     ! 14o(a,p)17f
     call rate_he4_o14_to_p_f17(tfactors,rate,dratedt)
     rpar(irp_rates-1+irap14o,PT_INDEX) = dens*rate
     rpar(irp_drtdt-1+irap14o,PT_INDEX) = dens*dratedt
   
     ! limit CNO as minimum between 14n(p,g)15o and 15o(beta nu)15n
     ! we store the limited rate in irlambCNO; this is lambda_CNO in WW81
     call rate_p_n14_to_o15(tfactors,rate,dratedt)
     rpar(irp_rates-1+irlambCNO,PT_INDEX) = min(rpar(irp_rates-1+irwk15o,PT_INDEX),rate*dens*y(ih1,PT_INDEX))
     if (rpar(irp_rates-1+irlambCNO,PT_INDEX) < rpar(irp_rates-1+irwk15o,PT_INDEX)) then
        rpar(irp_drtdt-1+irlambCNO,PT_INDEX) = dens*y(ih1,PT_INDEX)*dratedt
        rpar(irp_dlambCNOdh1,PT_INDEX) = rate*dens
     endif
   
     ! 22mg(...)30s
     ! check if this proceeds via p-captures or (a,p) reactions
     ! the Lweak is from WW81, eqn C15
     ! we store the rate in irlambda1; this is the lambda1 in WW81
     call rate_he4_si26_to_p_p29(tfactors,rate,dratedt)
     rpar(irp_rates-1+irlambda1,PT_INDEX) = max(rpar(irp_rates-1+irLweak,PT_INDEX),dens*y(ihe4,PT_INDEX)*rate)
     if (rpar(irp_rates-1+irlambda1,PT_INDEX) > rpar(irp_rates-1+irLweak,PT_INDEX)) then
          rpar(irp_drtdt-1+irlambda1,PT_INDEX) = dens*y(ihe4,PT_INDEX)*dratedt
          rpar(irp_dlambda1dhe4,PT_INDEX) = dens*rate
     ! use the sign of rpar(irp_rates-1+irlambda1) to indicate the value of delta1 in WW81
     ! if delta1 = 1, then we multiply the rate by -1
          rpar(irp_rates-1+irlambda1,PT_INDEX) = -ONE*rpar(irp_rates-1+irlambda1,PT_INDEX)
     endif
   
     ! 30s(...) 56ni 
     ! check if this proceeds via p-captures or (a,p) reactions
     ! use 44ti(a,p)v47 as a typical limiting rate for the (a,p) process
     ! store this in irlambda2; this is lambda2 in WW81
     call rate_he4_ti44_to_p_v47(tfactors,rate,dratedt)
     rpar(irp_rates-1+irlambda2,PT_INDEX) = max(rpar(irp_rates-1+irla2,PT_INDEX),dens*y(ihe4,PT_INDEX)*rate)
     if (rpar(irp_rates-1+irlambda2,PT_INDEX) > rpar(irp_rates-1+irla2,PT_INDEX)) then
          rpar(irp_drtdt-1+irlambda2,PT_INDEX) = dens*y(ihe4,PT_INDEX)*dratedt
          rpar(irp_dlambda2dhe4,PT_INDEX) = dens*rate
     ! use the sign of rpar(irp_rates-1+irlambda2) to indicate th value of delta2
     ! if delta2 = 1, then we multiply the rate by -1
          rpar(irp_rates-1+irlambda2,PT_INDEX) = -ONE*rpar(irp_rates-1+irlambda2,PT_INDEX)
     endif
   
     ! form s1 from WW81; branching ratio for 18ne beta decay (wk18ne) vs (a,p)
     ! store result in irs1
     ! 18ne(a,p)21na
     call rate_he4_ne18_to_p_na21(tfactors,rate,dratedt)
     rpar(irp_rates-1+irs1,PT_INDEX) = wk18ne / (wk18ne + dens*y(ihe4,PT_INDEX)*rate)
     rpar(irp_drtdt-1+irs1,PT_INDEX) = -rpar(irp_rates-1+irs1,PT_INDEX)*dens*y(ihe4,PT_INDEX)*dratedt &
          / (wk18ne + dens*y(ihe4,PT_INDEX)*rate)
     rpar(irp_drs1dhe4,PT_INDEX) = -rpar(irp_rates-1+irs1,PT_INDEX)*dens*rate &
          / (wk18ne + dens*y(ihe4,PT_INDEX)*rate)
   
     ! form r1 from WW81; ranching ratio for 19ne beta decay (wk19ne) vs (p,g)
     ! store result in irr1
     ! 19ne(p,g)20na
     call rate_p_ne19_to_na20(tfactors,rate,dratedt)
     rpar(irp_rates-1+irr1,PT_INDEX) = wk19ne / (wk19ne + dens*y(ih1,PT_INDEX)*rate)
     rpar(irp_drtdt-1+irr1,PT_INDEX) = -rpar(irp_rates-1+irr1,PT_INDEX)*dens*y(ih1,PT_INDEX)*dratedt &
          / (wk19ne + dens*y(ih1,PT_INDEX)*rate)
     rpar(irp_drr1dh1,PT_INDEX) = -rpar(irp_rates-1+irr1,PT_INDEX)*dens*rate &
          / (wk19ne + dens*y(ih1,PT_INDEX)*rate)
   
   
     !....
     !....  additional coding for proton capture on 56ni to heavier elements
     !....   kludge    56ni+56p -> 2 (56ni) at a rate given by min
     !....   of 56ni(pg) and 57cu decay rate
     !....
     !....  use 56ni rate from wallace and woosley 1981
     t9i32=tfactors%t9i*sqrt(tfactors%t9i)
     r56pg=dens*(1.29e-02_dp_t*exp(-4.897_dp_t*tfactors%t9i) &
          +7.065e+03_dp_t*exp(-20.33_dp_t*tfactors%t9i))*t9i32
     dr56pgdt = -THREE*HALF*r56pg*tfactors%t9i + &
          dens*t9i32*tfactors%t9i*tfactors%t9i* &
          (4.897_dp_t*1.29e-2_dp_t*exp(-4.897_dp_t*tfactors%t9i) &
          +20.33_dp_t*7.065e3_dp_t*exp(-20.33_dp_t*tfactors%t9i))
     !....  use generic proton separation energy of 400 kev
     !....  8.02 -> 4.64
     !      cutoni=2.08d-10*dens*exp(8.02*t9m1)/t932
     cutoni=2.08e-10_dp_t*dens*exp(4.642_dp_t*tfactors%t9i)*t9i32
     dcutonidt = cutoni*tfactors%t9i*(-THREE*HALF - 4.642_dp_t*tfactors%t9i)
     r57decay=3.54_dp_t
     r56eff=min(r56pg,cutoni*r57decay)
   !   rpar(irp_r56eff) = r56eff
   !   if (r56eff < r56pg) rpar(irp_dr56effdt) = r57decay*dcutonidt
     rpar(irp_r56eff,PT_INDEX) = ZERO
     rpar(irp_dr56effdt,PT_INDEX) = ZERO
   
   end subroutine make_rates
   
   subroutine make_ydots(ymol,t9,rpar,dydt,doing_dratesdt)
   !$acc routine
   
     use bl_types
     use bl_constants_module
     use network
     use rpar_indices
     
     implicit none
   
     real(kind=dp_t), intent(IN   ) :: ymol(nspec), t9
     logical,         intent(IN   ) :: doing_dratesdt
     real(kind=dp_t), intent(INOUT) :: rpar(n_rpar_comps,burn_npts)
     real(kind=dp_t), intent(  OUT) :: dydt(nspec,burn_npts)
     
     integer :: irp_start
     real(kind=dp_t) :: delta1, delta2
     real(kind=dp_t) :: dens
   
     ! initialize
     dydt = ZERO
     dens = rpar(irp_dens,PT_INDEX)
   
     ! check to see if we are doing this with the t-derivatives
     ! if so, offset our starting index in rpar
     irp_start = irp_rates
     if (doing_dratesdt) irp_start = irp_drtdt
   
     if (.not. doing_dratesdt) then
        delta1 = ZERO; delta2 = ZERO
        ! figure out the delta's; we used negative rates to indicate delta=1
        if (rpar(irp_rates-1+irlambda1,PT_INDEX) < ZERO) then
           delta1 = ONE
           rpar(irp_rates-1+irlambda1,PT_INDEX) = -ONE*rpar(irp_rates-1+irlambda1,PT_INDEX)
        endif
        if (rpar(irp_rates-1+irlambda2,PT_INDEX) < ZERO) then
           delta2 = ONE
           rpar(irp_rates-1+irlambda2,PT_INDEX) = -ONE*rpar(irp_rates-1+irlambda2,PT_INDEX)
        endif
        rpar(irp_delta1,PT_INDEX) = delta1
        rpar(irp_delta2,PT_INDEX) = delta2
     endif
   
   ! setup ODEs
   !
   !....
   !.... 12c = 1
   !....
         dydt(ic12,PT_INDEX)=-ymol(ic12)*ymol(ih1)*rpar(irp_start-1+irpg12c,PT_INDEX) &
              +ymol(ihe4)**3*rpar(irp_start-1+ir3a,PT_INDEX)/SIX &
              +ymol(io15)*rpar(irp_start-1+irlambCNO,PT_INDEX)
   !....
   !.... 14o = 2
   !....
         dydt(io14,PT_INDEX)=-ymol(io14)*ymol(ihe4)*rpar(irp_start-1+irap14o,PT_INDEX) &
              -ymol(io14)*rpar(irp_start-1+irwk14o,PT_INDEX) &
              +ymol(ic12)*ymol(ih1)*rpar(irp_start-1+irpg12c,PT_INDEX)
   !....
   !.... 15o = 3
   !....
         dydt(io15,PT_INDEX)=ymol(io14)*rpar(irp_start-1+irwk14o,PT_INDEX) &
              -ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irag15o,PT_INDEX) &
              -ymol(io15)*rpar(irp_start-1+irlambCNO,PT_INDEX) &
              +ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f,PT_INDEX)*rpar(irp_start-1+irs1,PT_INDEX) &
              +ymol(if17)*rpar(irp_start-1+irwk17f,PT_INDEX)
   !....
   !.... 16o = 4
   !....
         dydt(io16,PT_INDEX) = ymol(if17)*rpar(irp_start-1+irgp17f,PT_INDEX) &
              -ymol(io16)*ymol(ih1)*rpar(irp_start-1+irpg16o,PT_INDEX) &
              +ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irr1,PT_INDEX)*rpar(irp_start-1+irag15o,PT_INDEX) &
              -ymol(io16)*ymol(ihe4)*rpar(irp_start-1+irag16o,PT_INDEX)
   !....
   !.... 17f = 5
   !....
         dydt(if17,PT_INDEX)=ymol(io14)*ymol(ihe4)*rpar(irp_start-1+irap14o,PT_INDEX) &
              +ymol(io16)*ymol(ih1)*rpar(irp_start-1+irpg16o,PT_INDEX) &
              -ymol(if17)*rpar(irp_start-1+irgp17f,PT_INDEX) &
              -ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f,PT_INDEX) &
              -ymol(if17)*rpar(irp_start-1+irwk17f,PT_INDEX)
   !....
   !.... 22mg = 6
   !....
         dydt(img22,PT_INDEX)=ymol(io16)*ymol(ihe4)*rpar(irp_start-1+irag16o,PT_INDEX) &
              +ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f,PT_INDEX)*(ONE-rpar(irp_start-1+irs1,PT_INDEX)) &
              +ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irag15o,PT_INDEX)*(ONE-rpar(irp_start-1+irr1,PT_INDEX)) &
              -ymol(img22)*rpar(irp_start-1+irlambda1,PT_INDEX)
   !....
   !.... 30s = 7
   !....
         dydt(is30,PT_INDEX)=ymol(img22)*rpar(irp_start-1+irlambda1,PT_INDEX) &
              -ymol(is30)*rpar(irp_start-1+irlambda2,PT_INDEX)
   !....
   !.... amax (56ni) = 8  (note that WW81 have a typo -- they write lambda1 here)
   !....
         dydt(ini56,PT_INDEX)=ymol(is30)*rpar(irp_start-1+irlambda2,PT_INDEX)
   !....
   !.... 4he (alpha) = 9
   !....
         dydt(ihe4,PT_INDEX)=-ymol(ihe4)**3*HALF*rpar(irp_start-1+ir3a,PT_INDEX) &
              +ymol(io15)*rpar(irp_start-1+irlambCNO,PT_INDEX) &
              -ymol(io14)*ymol(ihe4)*rpar(irp_start-1+irap14o,PT_INDEX) &
              +ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f,PT_INDEX)*rpar(irp_start-1+irs1,PT_INDEX) &
              -ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irag15o,PT_INDEX)*(ONE-rpar(irp_start-1+irr1,PT_INDEX)) &
              -ymol(io16)*ymol(ihe4)*rpar(irp_start-1+irag16o,PT_INDEX) &
              -ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f,PT_INDEX)*(ONE-rpar(irp_start-1+irs1,PT_INDEX)) &
              +ymol(if17)*rpar(irp_start-1+irwk17f,PT_INDEX) &
              -TWO*ymol(img22)*rpar(irp_start-1+irlambda1,PT_INDEX)*rpar(irp_delta1,PT_INDEX) &
              -6.5d0*ymol(is30)*rpar(irp_start-1+irlambda2,PT_INDEX)*rpar(irp_delta2,PT_INDEX)
   !....
   !.... 1h (p) = 10
   !....
         dydt(ih1,PT_INDEX)=-ymol(io14)*rpar(irp_start-1+irwk14o,PT_INDEX) &
              -ymol(io15)*rpar(irp_start-1+irlambCNO,PT_INDEX) &
              -TWO*ymol(ic12)*ymol(ih1)*rpar(irp_start-1+irpg12c,PT_INDEX) &
              +ymol(io14)*ymol(ihe4)*rpar(irp_start-1+irap14o,PT_INDEX) &
              -TWO*ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f,PT_INDEX)*rpar(irp_start-1+irs1,PT_INDEX) &
              +ymol(if17)*rpar(irp_start-1+irgp17f,PT_INDEX) &
              -ymol(io16)*ymol(ih1)*rpar(irp_start-1+irpg16o,PT_INDEX) &
              -ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irag15o,PT_INDEX)*rpar(irp_start-1+irr1,PT_INDEX) &
              -TWO*ymol(io16)*ymol(ihe4)*rpar(irp_start-1+irag16o,PT_INDEX) &
              -THREE*ymol(io15)*ymol(ihe4)*rpar(irp_start-1+irag15o,PT_INDEX)*(ONE-rpar(irp_start-1+irr1,PT_INDEX)) &
              -ymol(if17)*ymol(ih1)*rpar(irp_start-1+irpg17f,PT_INDEX)*(ONE-rpar(irp_start-1+irs1,PT_INDEX)) &
              -TWO*ymol(if17)*rpar(irp_start-1+irwk17f,PT_INDEX) &
              -EIGHT*ymol(img22)*rpar(irp_start-1+irlambda1,PT_INDEX)*(ONE-rpar(irp_delta1,PT_INDEX)) &
              -26.e0_dp_t*ymol(is30)*rpar(irp_start-1+irlambda2,PT_INDEX)*(ONE-rpar(irp_delta2,PT_INDEX))
   
   
         if (.not. doing_dratesdt) then
            dydt(ini56,PT_INDEX) = dydt(ini56,PT_INDEX)+ymol(ini56)*ymol(ih1)*rpar(irp_r56eff,PT_INDEX)
            dydt(ih1,PT_INDEX) = dydt(ih1,PT_INDEX)-56.0d0*ymol(ini56)*ymol(ih1)*rpar(irp_r56eff,PT_INDEX)
         else
            dydt(ini56,PT_INDEX) = dydt(ini56,PT_INDEX)+ymol(ini56)*ymol(ih1)*rpar(irp_dr56effdt,PT_INDEX)
            dydt(ih1,PT_INDEX) = dydt(ih1,PT_INDEX)-56.0d0*ymol(ini56)*ymol(ih1)*rpar(irp_dr56effdt,PT_INDEX)
         endif
   
   end subroutine make_ydots
   
   !subroutine jac(neq, t, y, ml, mu, pd, nrpd, rpar, ipar)
   subroutine jac_vec(y, t, pd, rpar)
   !$acc routine seq
     use bl_types
     use bl_constants_module
     use network
     use rpar_indices
   
     implicit none
   
     !integer        , intent(IN   ) :: neq, ml, mu, nrpd, ipar
     real(kind=dp_t), intent(IN   ) :: y(neqs,burn_npts), t
     real(kind=dp_t), intent(INOUT) :: rpar(n_rpar_comps,burn_npts)
     real(kind=dp_t), intent(  OUT) :: pd(neqs,neqs,burn_npts)
   
     real(kind=dp_t), parameter :: T2T9 = 1.0e-9_dp_t
     real(kind=dp_t) :: ymol(nspec), t9
     real(kind=dp_t) :: cp, dhdX(nspec)
     real(kind=dp_t) :: psum
     integer :: it9
     integer :: i, j
   
   
     ! initialize
     pd(:,:,:)  = ZERO
     ymol = y(1:nspec,PT_INDEX)
     it9 = neqs
     t9 = y(it9,PT_INDEX)
     cp = rpar(irp_cp,PT_INDEX)
     dhdX = rpar(irp_dhdX:irp_dhdX+nspec-1,PT_INDEX)
   
     ! carbon-12
     pd(ic12,ic12,PT_INDEX) = -ymol(ih1)*rpar(irp_rates-1+irpg12c,PT_INDEX)
     pd(ic12,io15,PT_INDEX) = rpar(irp_rates-1+irlambCNO,PT_INDEX)
     pd(ic12,ihe4,PT_INDEX) = HALF*ymol(ihe4)*ymol(ihe4)*rpar(irp_rates-1+ir3a,PT_INDEX)
     pd(ic12,ih1,PT_INDEX)  = -ymol(ic12)*rpar(irp_rates-1+irpg12c,PT_INDEX) &
          +ymol(io15)*rpar(irp_dlambCNOdh1,PT_INDEX)
   
     ! oxygen-14
     pd(io14,ic12,PT_INDEX) = ymol(ih1)*rpar(irp_rates-1+irpg12c,PT_INDEX)
     pd(io14,io14,PT_INDEX) = -ymol(ihe4)*rpar(irp_rates-1+irap14o,PT_INDEX) &
          -rpar(irp_rates-1+irwk14o,PT_INDEX)
     pd(io14,ihe4,PT_INDEX) = -ymol(io14)*rpar(irp_rates-1+irap14o,PT_INDEX)
     pd(io14,ih1,PT_INDEX)  = ymol(ic12)*rpar(irp_rates-1+irpg12c,PT_INDEX)
   
     ! oxygen-15
     pd(io15,io14,PT_INDEX) = rpar(irp_rates-1+irwk14o,PT_INDEX)
     pd(io15,io15,PT_INDEX) = -ymol(ihe4)*rpar(irp_rates-1+irag15o,PT_INDEX) &
          -rpar(irp_rates-1+irlambCNO,PT_INDEX)
     pd(io15,if17,PT_INDEX) = ymol(ih1)*rpar(irp_rates-1+irpg17f,PT_INDEX)*rpar(irp_rates-1+irs1,PT_INDEX) &
          +rpar(irp_rates-1+irwk17f,PT_INDEX)
     pd(io15,ihe4,PT_INDEX) = -ymol(io15)*rpar(irp_rates-1+irag15o,PT_INDEX) &
          +ymol(if17)*ymol(ih1)*rpar(irp_rates-1+irpg17f,PT_INDEX)*rpar(irp_drs1dhe4,PT_INDEX)
     pd(io15,ih1,PT_INDEX)  = ymol(if17)*rpar(irp_rates-1+irpg17f,PT_INDEX)*rpar(irp_rates-1+irs1,PT_INDEX) &
          -ymol(io15)*rpar(irp_dlambCNOdh1,PT_INDEX)
   
     ! oxygen-16
     pd(io16,io15,PT_INDEX) = ymol(ihe4)*rpar(irp_rates-1+irr1,PT_INDEX)*rpar(irp_rates-1+irag15o,PT_INDEX)
     pd(io16,io16,PT_INDEX) = -ymol(ih1)*rpar(irp_rates-1+irpg16o,PT_INDEX) &
          -ymol(ihe4)*rpar(irp_rates-1+irag16o,PT_INDEX)
     pd(io16,if17,PT_INDEX) = rpar(irp_rates-1+irgp17f,PT_INDEX)
     pd(io16,ihe4,PT_INDEX) = ymol(io15)*rpar(irp_rates-1+irr1,PT_INDEX)*rpar(irp_rates-1+irag15o,PT_INDEX) &
          -ymol(io16)*rpar(irp_rates-1+irag16o,PT_INDEX)
     pd(io16,ih1,PT_INDEX)  = -ymol(io16)*rpar(irp_rates-1+irpg16o,PT_INDEX) &
          +ymol(io15)*ymol(ihe4)*rpar(irp_drr1dh1,PT_INDEX)*rpar(irp_rates-1+irag15o,PT_INDEX)
   
     ! flourine-17
     pd(if17,io14,PT_INDEX) = ymol(ihe4)*rpar(irp_rates-1+irap14o,PT_INDEX)
     pd(if17,io16,PT_INDEX) = ymol(ih1)*rpar(irp_rates-1+irpg16o,PT_INDEX)
     pd(if17,if17,PT_INDEX) = -rpar(irp_rates-1+irgp17f,PT_INDEX) &
          -ymol(ih1)*rpar(irp_rates-1+irpg17f,PT_INDEX) &
          -rpar(irp_rates-1+irwk17f,PT_INDEX)
     pd(if17,ihe4,PT_INDEX) = ymol(io14)*rpar(irp_rates-1+irap14o,PT_INDEX)
     pd(if17,ih1,PT_INDEX)  = ymol(io16)*rpar(irp_rates-1+irpg16o,PT_INDEX) &
          -ymol(if17)*rpar(irp_rates-1+irpg17f,PT_INDEX)
   
     ! magnesium-22
     pd(img22,io15,PT_INDEX) = ymol(ihe4)*rpar(irp_rates-1+irag15o,PT_INDEX)*(ONE-rpar(irp_rates-1+irr1,PT_INDEX))
     pd(img22,io16,PT_INDEX) = ymol(ihe4)*rpar(irp_rates-1+irag16o,PT_INDEX)
     pd(img22,if17,PT_INDEX) = ymol(ih1)*rpar(irp_rates-1+irpg17f,PT_INDEX)*(ONE-rpar(irp_rates-1+irs1,PT_INDEX))
     pd(img22,img22,PT_INDEX) = -rpar(irp_rates-1+irlambda1,PT_INDEX)
     pd(img22,ihe4,PT_INDEX) = ymol(io16)*rpar(irp_rates-1+irag16o,PT_INDEX) &
          +ymol(io15)*rpar(irp_rates-1+irag15o,PT_INDEX)*(ONE-rpar(irp_rates-1+irr1,PT_INDEX)) &
          -ymol(if17)*ymol(ih1)*rpar(irp_rates-1+irpg17f,PT_INDEX)*rpar(irp_drs1dhe4,PT_INDEX) &
          -ymol(img22)*rpar(irp_dlambda1dhe4,PT_INDEX)
     pd(img22,ih1,PT_INDEX)  = ymol(if17)*rpar(irp_rates-1+irpg17f,PT_INDEX)*(ONE-rpar(irp_rates-1+irs1,PT_INDEX)) &
          -ymol(io15)*ymol(ihe4)*rpar(irp_rates-1+irag15o,PT_INDEX)*rpar(irp_drr1dh1,PT_INDEX)
   
     ! sulfur-30
     pd(is30,img22,PT_INDEX) = rpar(irp_rates-1+irlambda1,PT_INDEX)
     pd(is30,is30,PT_INDEX)  = -rpar(irp_rates-1+irlambda2,PT_INDEX)
     pd(is30,ihe4,PT_INDEX)  = ymol(img22)*rpar(irp_dlambda1dhe4,PT_INDEX) &
          -ymol(is30)*rpar(irp_dlambda2dhe4,PT_INDEX)
   
     ! nickel-56
     pd(ini56,is30,PT_INDEX) = rpar(irp_rates-1+irlambda2,PT_INDEX)
     pd(ini56,ini56,PT_INDEX) = ymol(ih1)*rpar(irp_r56eff,PT_INDEX)
     pd(ini56,ihe4,PT_INDEX) = ymol(is30)*rpar(irp_dlambda2dhe4,PT_INDEX)
     pd(ini56,ih1,PT_INDEX) = ymol(ini56)*rpar(irp_r56eff,PT_INDEX)
   
     ! helium-4
     pd(ihe4,io14,PT_INDEX) = -ymol(ihe4)*rpar(irp_rates-1+irap14o,PT_INDEX)
     pd(ihe4,io15,PT_INDEX) = rpar(irp_rates-1+irlambCNO,PT_INDEX) &
          -ymol(ihe4)*rpar(irp_rates-1+irag15o,PT_INDEX)*(ONE-rpar(irp_rates-1+irr1,PT_INDEX))
     pd(ihe4,io16,PT_INDEX) = -ymol(ihe4)*rpar(irp_rates-1+irag16o,PT_INDEX)
     pd(ihe4,if17,PT_INDEX) = ymol(ih1)*rpar(irp_rates-1+irpg17f,PT_INDEX)*rpar(irp_rates-1+irs1,PT_INDEX) &
          -ymol(ih1)*rpar(irp_rates-1+irpg17f,PT_INDEX)*(ONE-rpar(irp_rates-1+irs1,PT_INDEX)) &
          +rpar(irp_rates-1+irwk17f,PT_INDEX)
     pd(ihe4,img22,PT_INDEX) = -TWO*rpar(irp_rates-1+irlambda1,PT_INDEX)*rpar(irp_delta1,PT_INDEX)
     pd(ihe4,is30,PT_INDEX) = -6.5d0*rpar(irp_rates-1+irlambda2,PT_INDEX)*rpar(irp_delta2,PT_INDEX)
     pd(ihe4,ihe4,PT_INDEX) = -THREE*ymol(ihe4)*ymol(ihe4)*HALF*rpar(irp_rates-1+ir3a,PT_INDEX) &
          -ymol(io14)*rpar(irp_rates-1+irap14o,PT_INDEX) &
          -ymol(io16)*rpar(irp_rates-1+irag16o,PT_INDEX) &
          -ymol(io15)*rpar(irp_rates-1+irag15o,PT_INDEX)*(ONE-rpar(irp_rates-1+irr1,PT_INDEX)) &
          +ymol(if17)*ymol(ih1)*rpar(irp_rates-1+irpg17f,PT_INDEX)*rpar(irp_drs1dhe4,PT_INDEX) &
          +ymol(if17)*ymol(ih1)*rpar(irp_rates-1+irpg17f,PT_INDEX)*rpar(irp_drs1dhe4,PT_INDEX) &
          -TWO*ymol(img22)*rpar(irp_dlambda1dhe4,PT_INDEX)*rpar(irp_delta1,PT_INDEX) &
          -6.5d0*ymol(is30)*rpar(irp_dlambda2dhe4,PT_INDEX)*rpar(irp_delta2,PT_INDEX)
     pd(ihe4,ih1,PT_INDEX)  = ymol(if17)*rpar(irp_rates-1+irpg17f,PT_INDEX)*rpar(irp_rates-1+irs1,PT_INDEX) &
          -ymol(if17)*rpar(irp_rates-1+irpg17f,PT_INDEX)*(ONE-rpar(irp_rates-1+irs1,PT_INDEX)) &
          +ymol(io15)*rpar(irp_dlambCNOdh1,PT_INDEX) &
          +ymol(io15)*ymol(ihe4)*rpar(irp_rates-1+irag15o,PT_INDEX)*rpar(irp_drr1dh1,PT_INDEX)
   
     ! hydrogen-1
     pd(ih1,ic12,PT_INDEX) = -TWO*ymol(ih1)*rpar(irp_rates-1+irpg12c,PT_INDEX)
     pd(ih1,io14,PT_INDEX) = ymol(ihe4)*rpar(irp_rates-1+irap14o,PT_INDEX) &
          -rpar(irp_rates-1+irwk14o,PT_INDEX)
     pd(ih1,io15,PT_INDEX) = -rpar(irp_rates-1+irlambCNO,PT_INDEX) &
          -ymol(ihe4)*rpar(irp_rates-1+irag15o,PT_INDEX)*rpar(irp_rates-1+irr1,PT_INDEX) &
          -THREE*ymol(ihe4)*rpar(irp_rates-1+irag15o,PT_INDEX)*(ONE-rpar(irp_rates-1+irr1,PT_INDEX))
     pd(ih1,io16,PT_INDEX) = -ymol(ih1)*rpar(irp_rates-1+irpg16o,PT_INDEX) &
          -TWO*ymol(ihe4)*rpar(irp_rates-1+irag16o,PT_INDEX)
     pd(ih1,if17,PT_INDEX) = -TWO*ymol(ih1)*rpar(irp_rates-1+irpg17f,PT_INDEX)*rpar(irp_rates-1+irs1,PT_INDEX) &
          +rpar(irp_rates-1+irgp17f,PT_INDEX) &
          -ymol(ih1)*rpar(irp_rates-1+irpg17f,PT_INDEX)*(ONE-rpar(irp_rates-1+irs1,PT_INDEX)) &
          -TWO*rpar(irp_rates-1+irwk17f,PT_INDEX)
     pd(ih1,img22,PT_INDEX) = -EIGHT*rpar(irp_rates-1+irlambda1,PT_INDEX)*(ONE-rpar(irp_delta1,PT_INDEX))
     pd(ih1,is30,PT_INDEX)  = -26.d0*rpar(irp_rates-1+irlambda2,PT_INDEX)*(ONE-rpar(irp_delta2,PT_INDEX))
     pd(ih1,ini56,PT_INDEX) = -56.0d0*ymol(ih1)*rpar(irp_r56eff,PT_INDEX)
     pd(ih1,ihe4,PT_INDEX) = ymol(io14)*rpar(irp_rates-1+irap14o,PT_INDEX) &
          -ymol(io15)*rpar(irp_rates-1+irag15o,PT_INDEX)*rpar(irp_rates-1+irr1,PT_INDEX) &
          -TWO*ymol(io16)*rpar(irp_rates-1+irag16o,PT_INDEX) &
          -THREE*ymol(io15)*rpar(irp_rates-1+irag15o,PT_INDEX)*(ONE-rpar(irp_rates-1+irr1,PT_INDEX)) &
          -ymol(if17)*ymol(ih1)*rpar(irp_rates-1+irpg17f,PT_INDEX)*rpar(irp_drs1dhe4,PT_INDEX) &
          -EIGHT*ymol(img22)*rpar(irp_dlambda1dhe4,PT_INDEX)*(ONE-rpar(irp_delta1,PT_INDEX)) &
          -26.d0*ymol(is30)*rpar(irp_dlambda2dhe4,PT_INDEX)*(ONE-rpar(irp_delta2,PT_INDEX))
     pd(ih1,ih1,PT_INDEX)  = -TWO*ymol(ic12)*rpar(irp_rates-1+irpg12c,PT_INDEX) &
          -TWO*ymol(if17)*rpar(irp_rates-1+irpg17f,PT_INDEX)*rpar(irp_rates-1+irs1,PT_INDEX) &
          -ymol(io16)*rpar(irp_rates-1+irpg16o,PT_INDEX) &
          -ymol(if17)*rpar(irp_rates-1+irpg17f,PT_INDEX)*(ONE-rpar(irp_rates-1+irs1,PT_INDEX)) &
          -ymol(io15)*rpar(irp_dlambCNOdh1,PT_INDEX) &
          +TWO*ymol(io15)*ymol(ihe4)*rpar(irp_rates-1+irag15o,PT_INDEX)*rpar(irp_drr1dh1,PT_INDEX) &
          -56.0d0*ymol(ini56)*rpar(irp_r56eff,PT_INDEX)
   
     ! temperature derivatives df(Y)/df(T)
     call make_ydots(ymol,t9,rpar,pd(1:nspec,it9,PT_INDEX),.true.)
   
     ! temperature jacobian elements
     do i = 1, neqs
        psum = 0.0d0
        do j = 1, nspec
           psum = psum + (dhdX(j) + ebin(j))*pd(j,i,PT_INDEX)
        enddo
        pd(it9,i,PT_INDEX) = -psum
     enddo
   
     pd(it9,:,PT_INDEX) = pd(it9,:,PT_INDEX) * T2T9 / cp
   
     return
   !end subroutine jac
   end subroutine jac_vec

end module feval
