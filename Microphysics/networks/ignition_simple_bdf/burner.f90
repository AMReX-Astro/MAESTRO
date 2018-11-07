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
   public :: burner, burner_vec

   !integer, public, save :: nst = 0
   !integer, public, save :: nfe = 0
   !integer, public, save :: nje = 0
   !integer, public, save :: nlu = 0
   !integer, public, save :: nit = 0

contains
   
   !subroutine burner(dens, temp, Xin, dt, Xout, rho_omegadot, rho_Hnuc)

   !   ! A wrapper for compatibility purposes
   !   ! outputs:
   !   !   Xout are the mass fractions after burning through timestep dt
   !   !   rho_omegadot = rho dX/dt
   !   !   rho_Hnuc = - sum_k q_k rho_omegadot_k  [erg / cm^3 / s]

   !   use burner_data
   !   use bdf

   !   implicit none

   !   !!Arguments
   !   real(kind=dp_t), intent(in   ) :: dens, temp, Xin(:), dt
   !   real(kind=dp_t), intent(  out) :: Xout(:), rho_omegadot(:), rho_Hnuc

   !   real(kind=dp_t) :: dens_arr(1), temp_arr(1), Xin_arr(nspec,1)
   !   real(kind=dp_t) :: Xout_arr(nspec,1), rho_omegadot_arr(nspec,1), rho_Hnuc_arr(1)

   !   dens_arr(1) = dens
   !   temp_arr(1) = temp
   !   Xin_arr(:,1) = Xin(:)

   !   call burner_vec(dens_arr, temp_arr, Xin_arr, dt, Xout_arr, &
   !                   rho_omegadot_arr, rho_Hnuc_arr)
   !   
   !   Xout(:) = Xout_arr(:,1)
   !   rho_omegadot(:) = rho_omegadot_arr(:,1)
   !   rho_Hnuc = rho_Hnuc_arr(1)

   !end subroutine burner

   subroutine burner(dens, temp, Xin, dt, Xout, rho_omegadot, rho_Hnuc)
   !$acc routine seq

      ! outputs:
      !   Xout are the mass fractions after burning through timestep dt
      !   rho_omegadot = rho dX/dt
      !   rho_Hnuc = - sum_k q_k rho_omegadot_k  [erg / cm^3 / s]

      use burner_data, only: neqs, burn_max_order, burn_npts, reset, reuse,    &
                             burn_npts, irp_dens, irp_cp, irp_dhdx, irp_o16,   &
                             irp_rate, irp_dratedt, irp_sc1212, irp_dsc1212dt, &
                             irp_xc12tmp, n_rpar_comps
      use network, only: nspec, nspec_advance, ic12_, io16_, img24_
      use bdf

      implicit none

      !!Arguments
      real(kind=dp_t), intent(in   ) :: dens, temp, Xin(nspec), dt
      real(kind=dp_t), intent(  out) :: Xout(nspec), rho_omegadot(nspec), rho_Hnuc
  
      !!Local constants
      real(kind=dp_t), parameter :: DT0 = 1.0d-9 !Initial dt to be used in getting from 
                                                 !t to tout.  This is arbitrary,
                                                 !multiple values should be
                                                 !explored, or an optimal value
                                                 !should be calculated.
      !!Local variables 
      integer :: n, i, j, npt, ierr, ierr_tot  ! npt is the size of the vectors
                                           !passed in, basically the number of
                                           !hydro cells. 
      !logical, save :: firstCall = .true.
      real(kind=dp_t), dimension(neqs) :: y, atol, rtol   ! input state, abs and rel tolerances
      real(kind=dp_t) :: eos_cp, eos_dhdX(nspec)
      real(kind=dp_t) :: upar(n_rpar_comps, burn_npts), y0(neqs, burn_npts), y1(neqs, burn_npts)
      real(kind=dp_t) :: t0, t1, enuc, dX
      type(eos_t)  :: eos_state
      type(bdf_ts) :: ts

      !!Execution
      !if (firstCall) then

      !   if (.NOT. network_initialized) then
      !      call bl_error("ERROR in burner: must initialize network first")
      !   endif
      ! 
      !   ic12 = network_species_index("carbon-12")
      !   io16 = network_species_index("oxygen-16")
      !   img24 = network_species_index("magnesium-24")
      !   
      !   if (ic12 < 0 .OR. io16 < 0 .OR. img24 < 0) then
      !      call bl_error("ERROR in burner: species undefined")
      !   endif
      !   
      !   firstCall = .false.
      !endif

      ! configure data that's the same for all pts
      
      ! set the tolerances.  We will be more relaxed on the temperature
      ! since it is only used in evaluating the rates.  
      atol(1:nspec_advance) = 1.d-12    ! mass fractions
      atol(nspec_advance+1) = 1.d-8     ! temperature
      rtol(1:nspec_advance) = 1.d-12    ! mass fractions
      rtol(nspec_advance+1) = 1.d-5     ! temperature

      ierr_tot = 0
      ! Call EoS.  Maestro's being redesigned to have already done
      ! this, so this is just a temporary hack for the purposes of rapid GPU development.
      ! we need the specific heat at constant pressure and dhdX |_p.  Take
      ! T, rho, Xin as input
      eos_state%rho   = dens
      eos_state%T     = temp
      do i = 1, nspec
         eos_state%xn(i) = Xin(i)
      enddo
         
      !call eos(eos_input_rt, eos_state)

      eos_cp = eos_state%cp
      do i = 1, nspec
         eos_dhdX(i) = eos_state%dhdX(i)
      enddo

      ! Build the bdf_ts time-stepper object
      call bdf_ts_build(ts, rtol, atol, upar)

      ! abundances are the first nspec_advance values and temperature is the last
      y(ic12_) = Xin(ic12_)
      y(nspec_advance+1) = temp
      
      ! density, specific heat at constant pressure, c_p, and dhdX are needed
      ! in the righthand side routine, so we will pass these in through the
      ! burner_aux module.
      !
      ! Since evaluating the EOS is expensive, we don't call it for every RHS
      ! call -- instead we freeze these values over the timestep.
      ! Since we are only integrating C12, we will need the O16 mass fraction
      ! in the RHS routine to compute the screening (and we know that the
      ! Mg24 abundance is constraint so things add to 1).
      ts%upar(irp_dens,1) = dens
      !ts%upar(irp_cp,1)   = eos_cp
      ts%upar(irp_cp,1)   = 19451875.637384996
      !ts(i)%upar(irp_dhdX:irp_dhdX-1+nspec,1) = eos_dhdX(i)
      j=1
      do i = irp_dhdX, irp_dhdX-1+nspec
         !We replace the array notation assignment commented out above
         !because such operations often cause errors on GPU.
         ts%upar(i,1) = eos_dhdX(j)
         j = j + 1
      end do
      ts%upar(irp_o16,1)  = Xin(io16_)

      !y0(:,1) = y
      do i = 1, neqs
         y0(i,1) = y(i)
      end do
      t0 = ZERO
      t1 = dt
      call bdf_advance(ts, y0, t0, y1, t1, &
                       DT0, reset, reuse, ierr, .true.)
      !y = y1(:,1)
      do i = 1, neqs
         y(i) = y1(i,1)
      end do
      ierr_tot = ierr_tot + ierr

      ! store the new mass fractions -- note, we discard the temperature
      ! here and instead compute the energy release from the binding
      ! energy -- make sure that they are positive
      Xout(ic12_)  = max(y(ic12_), ZERO)
      Xout(io16_)  = Xin(io16_)
      Xout(img24_) = ONE - Xout(ic12_) - Xout(io16_)

      ! compute the energy release.  Our convention is that the binding 
      ! energies are negative, so the energy release is
      ! - sum_k { (Xout(k) - Xin(k)) ebin(k) }
      !
      ! since this version of the network only evolves C12, we can
      ! compute the energy release easily
      enuc = (ebin(img24_) - ebin(ic12_))*(Xout(ic12_) - Xin(ic12_))

      ! also compute the density-weighted creation rates, rho_omegadot
      do i = 1, nspec
         dX = Xout(i) - Xin(i) 
         rho_omegadot(i) = dens * dX / dt
      enddo

      rho_Hnuc = dens*enuc/dt

      call bdf_ts_destroy(ts)

      !TODO: Here I'm using fact I know success is 0, need to update this since
      !      we're looping over cells and have ierr_tot now instead of single ierr
      !if (ierr_tot /= BDF_ERR_SUCCESS) then
      !   print *, 'ERROR: integration failed'
      !   print *, 'sum(ierr) for all GPU threads: ', ierr_tot
      !   call  bl_error("ERROR in burner: integration failed")
      !endif
    
      !TODO: This worked for dvode -- implement in bdf
      !nst = nst + iwork(11)
      !nfe = nfe + iwork(12)
      !nje = nje + iwork(13)
      !nlu = nlu + iwork(19)
      !nit = nit + iwork(20)
   end subroutine burner
end module burner_module
