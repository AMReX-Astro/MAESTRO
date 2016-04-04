module burner_module

   use bl_types
   use bl_constants_module
   use bl_error_module
   use eos_module, only: eos_input_rt, eos
   use eos_type_module
   use network
   use bdf
   use parallel
   
   private
   public :: burner_vec

   integer, public, save :: nst = 0
   integer, public, save :: nfe = 0
   integer, public, save :: nje = 0
   integer, public, save :: nlu = 0
   integer, public, save :: nit = 0

contains

   subroutine burner_vec(dens, temp, Xin, dt, Xout, rho_omegadot, rho_Hnuc)

      ! outputs:
      !   Xout are the mass fractions after burning through timestep dt
      !   rho_omegadot = rho dX/dt
      !   rho_Hnuc = - sum_k q_k rho_omegadot_k  [erg / cm^3 / s]

      use rpar_indices
      use bdf

      implicit none

      !!Arguments
      real(kind=dp_t), intent(in   ) :: dens(:), temp(:), Xin(:,:), dt
      real(kind=dp_t), intent(  out) :: Xout(:,:), rho_omegadot(:,:), rho_Hnuc(:)
  
      !!Local constants
      integer, parameter :: NEQ = 1 + nspec_advance   ! set the number of independent variables -- this should be temperature
                                                      ! + the number of species
      integer, parameter :: MAX_ORDER = 5   !This is arbitrary, should investigate other values
      logical, parameter :: RESET = .true.  !.true. means we want to initialize the bdf_ts object
      logical, parameter :: REUSE = .false. !.false. means don't reuse the Jacobian
      real(kind=dp_t), parameter :: DT0 = 1.0d-9 !Initial dt to be used in getting from 
                                                 !t to tout.  Also arbitrary,
                                                 !multiple values should be
                                                 !explored.
      !!Local variables 
      integer :: n, i, j, npt, bdf_npt, ierr, ierr_tot  ! npt is the size of the vectors
                                           !passed in, basically the number of
                                           !hydro cells. bdf_npt refers to bdf's
                                           !own notion of vectors, which we're
                                           !currently not exploiting.
      integer, save :: ic12, io16, img24   ! we will always refer to the species by integer indices that come from
                                           ! the network module -- this makes things robust to a shuffling of the 
                                           ! species ordering
      logical, save :: firstCall = .true.
      real(kind=dp_t), dimension(NEQ) :: y, atol, rtol   ! input state, abs and rel tolerances
      real(kind=dp_t), allocatable :: eos_cp(:), eos_dhdX(:,:)
      real(kind=dp_t), allocatable :: upar(:,:), y0(:,:), y1(:,:)
      real(kind=dp_t) :: t0, t1, enuc, dX
      type (eos_t) :: eos_state
      type(bdf_ts), allocatable :: ts(:)

      !!Execution
      if (firstCall) then

         if (.NOT. network_initialized) then
            call bl_error("ERROR in burner: must initialize network first")
         endif
       
         ic12 = network_species_index("carbon-12")
         io16 = network_species_index("oxygen-16")
         img24 = network_species_index("magnesium-24")
         
         if (ic12 < 0 .OR. io16 < 0 .OR. img24 < 0) then
            call bl_error("ERROR in burner: species undefined")
         endif
         
         firstCall = .false.
      endif

      ! configure data that's the same for all pts
      
      ! set the tolerances.  We will be more relaxed on the temperature
      ! since it is only used in evaluating the rates.  
      atol(1:nspec_advance) = 1.d-12    ! mass fractions
      atol(nspec_advance+1) = 1.d-8     ! temperature
      rtol(1:nspec_advance) = 1.d-12    ! mass fractions
      rtol(nspec_advance+1) = 1.d-5     ! temperature

      ! Call EoS on all points.  Maestro's being redesigned to have already done
      ! this, so this is just a temporary hack for the purposes of rapid GPU development.
      ! Also build array of vbdf time-stepper derived types.
      npt = size(dens)
      bdf_npt = 1
      ierr_tot = 0
      allocate(eos_cp(npt), eos_dhdX(nspec, npt))
      allocate(ts(npt))
      allocate(upar(n_rpar_comps, bdf_npt))
      allocate(y0(NEQ, bdf_npt), y1(NEQ, bdf_npt))

      !$acc enter data create(ts)
      !NOTE: We create ts and then update all members below (even scalar)
      !      instead of, say, doing a copyin followed by updates.  This is
      !      necessary because of how PGI implements these directives, and to
      !      some extent it makes sense regardless of compiler.  Doing a copyin
      !      and then updating serves to overwrite the device (GPU) pointer with
      !      the host (CPU) pointer.  Once deep copy is widely implemented we
      !      won't have to worry so much about this.
      do i = 1, npt
         ! we need the specific heat at constant pressure and dhdX |_p.  Take
         ! T, rho, Xin as input
         eos_state%rho   = dens(i)
         eos_state%T     = temp(i)
         eos_state%xn(:) = Xin(:,i)
            
         call eos(eos_input_rt, eos_state, .false.)

         eos_cp(i) = eos_state%cp
         eos_dhdX(:,i) = eos_state%dhdX(:)

         !We need to build and allocate bdf_ts objects before the OpenACC region,
         !because you cannot allocate within OpenACC regions.

         ! Build the bdf_ts time-stepper object
         call bdf_ts_build(ts(i), NEQ, bdf_npt, rtol, atol, MAX_ORDER, upar)

         !Now we update all non-dynamic data members, meaning those that aren't
         !allocatables, pointers, etc.  They can be arrays as long as they're
         !static.  To make my previous note more concrete, if we did something
         !like `update device(ts(i))` it would overwrite the device's pointer
         !address with that of the host, according to PGI/NVIDIA consults.
         !Updating individual members avoids this.
          
         !$acc update device(       &
         !$acc    ts(i)%neq,        &
         !$acc    ts(i)%npt,        &
         !$acc    ts(i)%max_order,  &
         !$acc    ts(i)%max_steps,  &
         !$acc    ts(i)%max_iters,  &
         !$acc    ts(i)%verbose,    &
         !$acc    ts(i)%dt_min,     &
         !$acc    ts(i)%eta_min,    &
         !$acc    ts(i)%eta_max,    &
         !$acc    ts(i)%eta_thresh, &
         !$acc    ts(i)%max_j_age,  &
         !$acc    ts(i)%max_p_age,  &
         !$acc    ts(i)%debug,      &
         !$acc    ts(i)%dump_unit,  &
         !$acc    ts(i)%t,          &
         !$acc    ts(i)%t1,         &
         !$acc    ts(i)%dt,         &
         !$acc    ts(i)%dt_nwt,     &
         !$acc    ts(i)%k,          &
         !$acc    ts(i)%n,          &
         !$acc    ts(i)%j_age,      &
         !$acc    ts(i)%p_age,      &
         !$acc    ts(i)%k_age,      &
         !$acc    ts(i)%tq,         &
         !$acc    ts(i)%tq2save,    &
         !$acc    ts(i)%temp_data,  &
         !$acc    ts(i)%refactor,   &
         !$acc    ts(i)%nfe,        &
         !$acc    ts(i)%nje,        &
         !$acc    ts(i)%nlu,        &
         !$acc    ts(i)%nit,        &
         !$acc    ts(i)%nse,        &
         !$acc    ts(i)%ncse,       &
         !$acc    ts(i)%ncit,       &
         !$acc    ts(i)%ncdtmin)

         !Now it's time to deal with dynamic data.  At the moment, they only
         !exist as pointers on the device.  For PGI at least, doing a copyin on
         !dynamic data serves to create, allocate, and then attach each dynamic
         !data member to the corresponding pointer in ts(i).
         !$acc enter data copyin(   &
         !$acc    ts(i)%rtol,       &
         !$acc    ts(i)%atol,       &
         !$acc    ts(i)%J,          &
         !$acc    ts(i)%P,          &
         !$acc    ts(i)%z(:,:,0:),  &
         !$acc    ts(i)%z0(:,:,0:), &
         !$acc    ts(i)%h(0:),      &
         !$acc    ts(i)%l(0:),      &
         !$acc    ts(i)%shift(0:),  &
         !$acc    ts(i)%upar,       &
         !$acc    ts(i)%y,          &
         !$acc    ts(i)%yd,         &
         !$acc    ts(i)%rhs,        &
         !$acc    ts(i)%e,          &
         !$acc    ts(i)%e1,         &
         !$acc    ts(i)%ewt,        &
         !$acc    ts(i)%b,          &
         !$acc    ts(i)%ipvt,       &
         !$acc    ts(i)%A(0:,0:))
      end do

      !NOTE: The (:)'s are not necessary.  I'm putting them here just for 
      !      clarity about the shape.
      !$acc data                                                               &
      !$acc copyin(dens(:), temp(:), eos_cp(:), eos_dhdX(:,:), Xin(:,:))       &
      !$acc copyout(Xout(:,:), rho_omegadot(:,:), rho_Hnuc(:)) 

      !$acc parallel loop gang vector present(dens, temp, eos_cp, eos_dhdX,    &
      !$acc    Xin, Xout, rho_omegadot, rho_Hnuc, ebin, ts)         &
      !$acc    private(ierr, y, y0, y1) reduction(+:ierr_tot)
      do i = 1, npt
         ! abundances are the first nspec_advance values and temperature is the last
         y(ic12) = Xin(ic12,i)
         y(nspec_advance+1) = temp(i)
         if (i==32) then
            ts(i)%temp_data(1,1) = Xin(ic12,i)
            ts(i)%temp_data(2,1) = temp(i)
         endif
         
         ! density, specific heat at constant pressure, c_p, and dhdX are needed
         ! in the righthand side routine, so we will pass these in through the
         ! burner_aux module.
         !
         ! Since evaluating the EOS is expensive, we don't call it for every RHS
         ! call -- instead we freeze these values over the timestep.
         ! Since we are only integrating C12, we will need the O16 mass fraction
         ! in the RHS routine to compute the screening (and we know that the
         ! Mg24 abundance is constraint so things add to 1).
         ts(i)%upar(irp_dens,1) = dens(i)
         ts(i)%upar(irp_cp,1)   = eos_cp(i)
         !ts(i)%upar(irp_dhdX:irp_dhdX-1+nspec,1) = eos_dhdX(:,i)
         j=1
         do n = irp_dhdX, irp_dhdX-1+nspec
            !We replace the array notation assignment commented out above
            !because such operations often cause errors on GPU.
            ts(i)%upar(n,1) = eos_dhdX(j,i)
            j = j + 1
         end do
         ts(i)%upar(irp_o16,1)  = Xin(io16,i)

         !y0(:,1) = y
         do n = 1, NEQ
            y0(n,1) = y(n)
         end do
         t0 = ZERO
         t1 = dt
         call bdf_advance(ts(i), NEQ, bdf_npt, y0, t0, y1, t1, &
                          DT0, RESET, REUSE, ierr, .true.)
         !y = y1(:,1)
         do n = 1, NEQ
            y(n) = y1(n,1)
         end do
         ierr_tot = ierr_tot + ierr

         ! store the new mass fractions -- note, we discard the temperature
         ! here and instead compute the energy release from the binding
         ! energy -- make sure that they are positive
         Xout(ic12,i)  = max(y(ic12), ZERO)
         Xout(io16,i)  = Xin(io16,i)
         Xout(img24,i) = ONE - Xout(ic12,i) - Xout(io16,i)

         ! compute the energy release.  Our convention is that the binding 
         ! energies are negative, so the energy release is
         ! - sum_k { (Xout(k) - Xin(k)) ebin(k) }
         !
         ! since this version of the network only evolves C12, we can
         ! compute the energy release easily
         enuc = (ebin(img24) - ebin(ic12))*(Xout(ic12,i) - Xin(ic12,i))

         ! also compute the density-weighted creation rates, rho_omegadot
         do n = 1, nspec
            dX = Xout(n,i) - Xin(n,i) 
            rho_omegadot(n,i) = dens(i) * dX / dt
         enddo

         rho_Hnuc(i) = dens(i)*enuc/dt
         !rho_Hnuc(i) = 3.25 
      end do

      !$acc update host(ts(32)%temp_data)
      ! Cleanup
      do i=1, npt
         !$acc exit data delete(     &
         !$acc    ts(i)%rtol(:),     &
         !$acc    ts(i)%atol(:),     &
         !$acc    ts(i)%J(:,:,:),    &
         !$acc    ts(i)%P(:,:,:),    &
         !$acc    ts(i)%z(:,:,:),    &
         !$acc    ts(i)%z0(:,:,:),   &
         !$acc    ts(i)%h(:),        &
         !$acc    ts(i)%l(:),        &
         !$acc    ts(i)%shift(:),    &
         !$acc    ts(i)%upar(:,:),   &
         !$acc    ts(i)%y(:,:),      &
         !$acc    ts(i)%yd(:,:),     &
         !$acc    ts(i)%rhs(:,:),    &
         !$acc    ts(i)%e(:,:),      &
         !$acc    ts(i)%e1(:,:),     &
         !$acc    ts(i)%ewt(:,:),    &
         !$acc    ts(i)%b(:,:),      &
         !$acc    ts(i)%ipvt(:,:),   &
         !$acc    ts(i)%A(:,:))
 
         call bdf_ts_destroy(ts(i))
      end do

      !WARNING! Do *not* do copyout, it'll break
      !$acc exit data delete(ts(:))
     
      !$acc end data 
      print *, 'Xout:         ', Xout(io16,32)
      print *, 'rho_omegadot: ', rho_omegadot(io16,32)
      print *, 'rho_Hnuc:     ', rho_Hnuc(32)
      print *, 'ierr_tot:     ', ierr_tot
      print *, 'temp_data:    ', ts(32)%temp_data
      !TODO: Here I'm using fact I know success is 0, need to update this since
      !      we're looping over cells and have ierr_tot now instead of single ierr
      if (ierr_tot /= BDF_ERR_SUCCESS) then
         print *, 'ERROR: integration failed'
         print *, 'sum(ierr) for all GPU threads: ', ierr_tot
         call  bl_error("ERROR in burner: integration failed")
      endif
    
      !TODO: This worked for dvode -- implement in bdf
      !nst = nst + iwork(11)
      !nfe = nfe + iwork(12)
      !nje = nje + iwork(13)
      !nlu = nlu + iwork(19)
      !nit = nit + iwork(20)
   end subroutine burner_vec
end module burner_module
