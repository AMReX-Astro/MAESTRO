module advance_timestep_module

  use bl_types           , only: dp_t
  use bl_constants_module, only: ZERO, HALF, TWO
  use multifab_module
  use ml_layout_module   , only: ml_layout
  use define_bc_module   , only: bc_tower
  use parallel           , only: parallel_IOProcessor, parallel_IOProcessorNode, &
                                 parallel_wtime, parallel_reduce, parallel_barrier, &
                                 MPI_MAX
  use particle_module

  implicit none

  private

  public :: advance_timestep

contains
    
  subroutine advance_timestep(init_mode,mla,uold,sold,unew,snew, &
                              gpi,pi,normal,rho0_old,rhoh0_old, &
                              rho0_new,rhoh0_new,p0_old,p0_new,tempbar,gamma1bar,w0, &
                              rho_omegadot2,rho_Hnuc2,rho_Hext,diff_new,intra,&
                              div_coeff_old,div_coeff_new, &
                              grav_cell_old,dx,dt,dtold,the_bc_tower, &
                              dSdt,Source_old,Source_new,etarho_ec,etarho_cc, &
                              psi,sponge,hgrhs,tempbar_init,particles)

    use bl_prof_module              , only : bl_prof_timer, build, destroy
    use      pre_advance_module     , only : advance_premac
    use velocity_advance_module     , only : velocity_advance
    use  density_advance_module     , only : density_advance
    use enthalpy_advance_module     , only : enthalpy_advance
    use make_div_coeff_module       , only : make_div_coeff
    use make_w0_module              , only : make_w0
    use advect_base_module          , only : advect_base_dens, advect_base_enthalpy
    use react_state_module          , only : react_state, instantaneous_reaction_rates
    use make_S_module               , only : make_S
    use average_module              , only : average
    use phihalf_module              , only : make_S_at_halftime, make_at_halftime
    use extraphalf_module           , only : extrap_to_halftime
    use thermal_conduct_module      , only : thermal_conduct_predictor, &
                                             thermal_conduct_corrector, &
                                             make_explicit_thermal_hterm
    use make_explicit_thermal_module, only : make_explicit_thermal, make_thermal_coeffs 
    use make_grav_module            , only : make_grav_cell
    use make_eta_module             , only : make_etarho_planar, make_etarho_spherical
    use make_psi_module             , only : make_psi_planar, make_psi_spherical
    use fill_3d_module              , only : put_1d_array_on_cart, make_w0mac, make_s0mac
    use cell_to_edge_module         , only : cell_to_edge
    use make_gamma_module           , only : make_gamma
    use rhoh_vs_t_module            , only : makePfromRhoH, makeTfromRhoP, makeTfromRhoH
    use diag_module                 , only : diag
    use sanity_module               , only : sanity_check
    use enforce_HSE_module          , only : enforce_HSE

    use macrhs_module               , only : make_macrhs
    use macproject_module           , only : macproject

    use hgrhs_module                , only : make_hgrhs, correct_hgrhs
    use hgproject_module            , only : hgproject
    use proj_parameters             , only : pressure_iters_comp, regular_timestep_comp

    use variables                   , only : nscal, ntrac, temp_comp, rho_comp, rhoh_comp, &
                                             foextrap_comp, spec_comp
    use geometry                    , only : nlevs_radial, spherical, nr_fine, compute_cutoff_coords
    use network                     , only : nspec
    use probin_module               , only : barrier_timers, evolve_base_state, fix_base_state, &
                                             use_etarho, dpdt_factor, verbose, &
                                             use_tfromp, use_thermal_diffusion, &
                                             use_delta_gamma1_term, nodal, mach_max_abort, &
                                             prob_lo, prob_hi, use_particles, sdc_iters, &
                                             enthalpy_pred_type, species_pred_type
    use time_module                 , only : time
    use addw0_module                , only : addw0
    use pred_parameters
    use make_intra_coeffs_module    , only: make_intra_coeffs

    use fabio_module

    logical,         intent(in   ) :: init_mode
    type(ml_layout), intent(inout) :: mla
    type(multifab),  intent(in   ) ::   uold(:)
    type(multifab),  intent(inout) ::   sold(:)
    type(multifab),  intent(inout) ::   unew(:)
    type(multifab),  intent(inout) ::   snew(:)
    type(multifab),  intent(inout) ::  gpi(:)
    type(multifab),  intent(inout) ::   pi(:)
    type(multifab),  intent(in   ) :: normal(:)
    real(dp_t)    ,  intent(inout) ::  rho0_old(:,0:)
    real(dp_t)    ,  intent(inout) :: rhoh0_old(:,0:)
    real(dp_t)    ,  intent(inout) ::  rho0_new(:,0:)
    real(dp_t)    ,  intent(inout) :: rhoh0_new(:,0:)
    real(dp_t)    ,  intent(inout) ::    p0_old(:,0:)
    real(dp_t)    ,  intent(inout) ::    p0_new(:,0:)
    real(dp_t)    ,  intent(inout) ::   tempbar(:,0:)
    real(dp_t)    ,  intent(inout) ::   tempbar_init(:,0:)
    real(dp_t)    ,  intent(inout) :: gamma1bar(:,0:)
    real(dp_t)    ,  intent(inout) ::        w0(:,0:)
    type(multifab),  intent(inout) :: rho_omegadot2(:)
    type(multifab),  intent(inout) :: rho_Hnuc2(:)
    type(multifab),  intent(inout) :: rho_Hext(:)
    type(multifab),  intent(inout) ::  diff_new(:)
    type(multifab),  intent(inout) ::     intra(:)
    real(dp_t)    ,  intent(inout) :: div_coeff_old(:,0:)
    real(dp_t)    ,  intent(inout) :: div_coeff_new(:,0:)
    real(dp_t)    ,  intent(inout) :: grav_cell_old(:,0:)
    real(dp_t)    ,  intent(in   ) :: dx(:,:),dt,dtold
    type(bc_tower),  intent(in   ) :: the_bc_tower
    type(multifab),  intent(inout) ::       dSdt(:)
    type(multifab),  intent(inout) :: Source_old(:)
    type(multifab),  intent(inout) :: Source_new(:)
    real(dp_t)    ,  intent(inout) ::  etarho_ec(:,0:)
    real(dp_t)    ,  intent(inout) ::  etarho_cc(:,0:)
    real(dp_t)    ,  intent(inout) ::        psi(:,0:)
    type(multifab),  intent(in   ) :: sponge(:)
    type(multifab),  intent(inout) ::  hgrhs(:)
    type(particle_container), intent(inout) :: particles

    ! local variables
    type(multifab) ::                shat(mla%nlevel)
    type(multifab) ::             rhohalf(mla%nlevel)
    type(multifab) ::             cphalf(mla%nlevel)
    type(multifab) ::             xihalf(mla%nlevel)
    type(multifab) ::       w0_force_cart(mla%nlevel)
    type(multifab) ::              macrhs(mla%nlevel)
    type(multifab) ::              macphi(mla%nlevel)
    type(multifab) ::           hgrhs_old(mla%nlevel)
    type(multifab) ::          Source_nph(mla%nlevel)
    type(multifab) ::            diff_old(mla%nlevel)
    type(multifab) ::            diff_hat(mla%nlevel)
    type(multifab) ::      diff_hterm_new(mla%nlevel)
    type(multifab) ::      diff_hterm_hat(mla%nlevel)
    type(multifab) ::        div_coeff_3d(mla%nlevel)
    type(multifab) :: div_coeff_cart_edge(mla%nlevel,mla%dim)
    type(multifab) ::              gamma1(mla%nlevel)
    type(multifab) ::          etarhoflux(mla%nlevel)

    ! coefficients for thermal conduction stuff
    type(multifab) ::         Tcoeff_old(mla%nlevel)
    type(multifab) ::         hcoeff_old(mla%nlevel)
    type(multifab) ::        Xkcoeff_old(mla%nlevel)
    type(multifab) ::         pcoeff_old(mla%nlevel)
    type(multifab) ::         Tcoeff_new(mla%nlevel)
    type(multifab) ::         hcoeff_new(mla%nlevel)
    type(multifab) ::        Xkcoeff_new(mla%nlevel)
    type(multifab) ::         pcoeff_new(mla%nlevel)

    ! used for dpdt volume discrepancy
    type(multifab) ::            peos_old(mla%nlevel)
    type(multifab) ::            peos_nph(mla%nlevel)
    type(multifab) ::            peos_new(mla%nlevel)
    type(multifab) ::        peosbar_cart(mla%nlevel)
    type(multifab) ::        delta_p_term(mla%nlevel)

    ! only used if delta_gamma1_term = T
    type(multifab) ::   delta_gamma1_term(mla%nlevel)
    type(multifab) ::        delta_gamma1(mla%nlevel)

    type(multifab) ::          scal_force(mla%nlevel)
    type(multifab) ::               w0mac(mla%nlevel,mla%dim)
    type(multifab) ::                umac(mla%nlevel,mla%dim)
    type(multifab) ::               sedge(mla%nlevel,mla%dim)
    type(multifab) ::               sflux(mla%nlevel,mla%dim)
    type(multifab) ::          sdc_source(mla%nlevel)
    type(multifab) ::                aofs(mla%nlevel)

    real(kind=dp_t), allocatable ::        grav_cell_nph(:,:)
    real(kind=dp_t), allocatable ::        grav_cell_new(:,:)
    real(kind=dp_t), allocatable ::             rho0_nph(:,:)
    real(kind=dp_t), allocatable ::               p0_nph(:,:)
    real(kind=dp_t), allocatable ::     p0_minus_peosbar(:,:)
    real(kind=dp_t), allocatable ::              peosbar(:,:)
    real(kind=dp_t), allocatable ::             w0_force(:,:)
    real(kind=dp_t), allocatable ::                 Sbar(:,:)
    real(kind=dp_t), allocatable ::        div_coeff_nph(:,:)
    real(kind=dp_t), allocatable ::        gamma1bar_old(:,:)
    real(kind=dp_t), allocatable ::      gamma1bar_temp1(:,:)
    real(kind=dp_t), allocatable ::      gamma1bar_temp2(:,:)
    real(kind=dp_t), allocatable :: delta_gamma1_termbar(:,:)
    real(kind=dp_t), allocatable ::               w0_old(:,:)
    real(kind=dp_t), allocatable ::       div_coeff_edge(:,:)
    real(kind=dp_t), allocatable ::  rho0_predicted_edge(:,:)

    integer    :: i,n,comp,proj_type,nlevs,dm,misdc
    real(dp_t) :: halfdt

    ! need long int to store numbers greater than 2^31
    integer(kind=ll_t) :: numcell

    real(kind=dp_t) :: advect_time , advect_time_start , advect_time_max
    real(kind=dp_t) :: macproj_time, macproj_time_start, macproj_time_max
    real(kind=dp_t) :: ndproj_time , ndproj_time_start , ndproj_time_max
    real(kind=dp_t) :: thermal_time, thermal_time_start, thermal_time_max
    real(kind=dp_t) :: react_time  , react_time_start  , react_time_max
    real(kind=dp_t) :: misc_time   , misc_time_start   , misc_time_max

    type(bl_prof_timer), save :: bpt

    call build(bpt, "advance_timestep")

    misc_time_start = parallel_wtime()

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< CALLING advance_timestep with dt =',dt 
    end if

    if (evolve_base_state) then
       call bl_error("evolve_base_state not supported with SDC")
    end if
    if (enthalpy_pred_type == 1) then
       call bl_error("enthalpy_pred_type == 1 not supported with SDC")
    end if

    allocate(       grav_cell_nph(nlevs_radial,0:nr_fine-1))
    allocate(       grav_cell_new(nlevs_radial,0:nr_fine-1))
    allocate(            rho0_nph(nlevs_radial,0:nr_fine-1))
    allocate(              p0_nph(nlevs_radial,0:nr_fine-1))
    allocate(    p0_minus_peosbar(nlevs_radial,0:nr_fine-1))
    allocate(             peosbar(nlevs_radial,0:nr_fine-1))
    allocate(            w0_force(nlevs_radial,0:nr_fine-1))
    allocate(                Sbar(nlevs_radial,0:nr_fine-1))
    allocate(       div_coeff_nph(nlevs_radial,0:nr_fine-1))
    allocate(       gamma1bar_old(nlevs_radial,0:nr_fine-1))
    allocate(     gamma1bar_temp1(nlevs_radial,0:nr_fine-1))
    allocate(     gamma1bar_temp2(nlevs_radial,0:nr_fine-1))
    allocate(delta_gamma1_termbar(nlevs_radial,0:nr_fine-1))
    allocate(              w0_old(nlevs_radial,0:nr_fine))
    allocate(      div_coeff_edge(nlevs_radial,0:nr_fine))
    allocate( rho0_predicted_edge(nlevs_radial,0:nr_fine))

    advect_time  = 0.d0
    macproj_time = 0.d0
    ndproj_time  = 0.d0
    thermal_time = 0.d0
    react_time   = 0.d0
    misc_time    = 0.d0

    nlevs = mla%nlevel
    dm = mla%dim

    halfdt = half*dt

    if (verbose .ge. 1) then

       if (parallel_IOProcessor()) then
          do n = 1, nlevs
             write(6,*) 'level: ', n
             write(6,*) '   number of boxes = ', nboxes(pi(n))
             write(6,*) '   maximum zones   = ', (extent(mla%mba%pd(n),i),i=1,dm)
          end do
       end if

       do n=1,nlevs
          numcell = multifab_volume(pi(n),.false.)
          if (parallel_IOProcessor()) then
             write(6,*) 'Number of valid cells at level        ',n,numcell
          end if
          numcell = multifab_volume(pi(n),.true.)
          if (parallel_IOProcessor()) then
             write(6,*) 'Number of valid + ghost cells at level',n,numcell
          end if
       end do
       if (parallel_IOProcessor()) then
          write(6,*) ''
       end if
    end if

    ! Initialize these to previous values
    w0_old        = w0
    gamma1bar_old = gamma1bar

    if (barrier_timers) call parallel_barrier()
    misc_time = misc_time + parallel_wtime() - misc_time_start
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 1 -- Compute advection velocities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    advect_time_start = parallel_wtime()

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP 1: Compute advection velocities'
    end if

    do n=1,nlevs
       call multifab_build(Source_nph(n), mla%la(n), 1, 1)
    end do

    if (time .eq. ZERO) then
       call make_S_at_halftime(mla,Source_nph,Source_old,Source_new, &
                               the_bc_tower%bc_tower_array)
    else
       call extrap_to_halftime(mla,Source_nph,dSdt,Source_old,dt, &
                               the_bc_tower%bc_tower_array)
    end if

    do n=1,nlevs
       call multifab_build(delta_p_term(n), mla%la(n), 1, 0)
       call setval(delta_p_term(n),ZERO,all=.true.)
    end do

    ! compute p0_minus_peosbar = p0_old - peosbar (for making w0) and
    ! compute delta_p_term = peos_old - peosbar_cart (for RHS of projections)
    if (dpdt_factor .gt. ZERO ) then
    
       do n=1,nlevs
          call multifab_build(peos_old(n), mla%la(n), 1, 0)
       end do

       ! peos_old now holds the thermodynamic p computed from sold(rho,h,X)
       call makePfromRhoH(sold,sold,peos_old,mla,the_bc_tower%bc_tower_array)

       ! compute peosbar = Avg(peos_old)
       call average(mla,peos_old,peosbar,dx,1)

       ! compute p0_minus_peosbar = p0_old - peosbar
       p0_minus_peosbar = p0_old - peosbar

       do n=1,nlevs
          call multifab_build(peosbar_cart(n), mla%la(n), 1, 0)
       end do
       
       ! compute peosbar_cart from peosbar
       call put_1d_array_on_cart(peosbar,peosbar_cart,foextrap_comp, &
                                 .false.,.false.,dx,the_bc_tower%bc_tower_array,mla)

       ! compute delta_p_term = peos_old - peosbar_cart
       do n=1,nlevs
          call multifab_copy(delta_p_term(n), peos_old(n))
          call multifab_sub_sub(delta_p_term(n), peosbar_cart(n))
       end do
       
       do n=1,nlevs
          call destroy(peosbar_cart(n))
       end do

    else

       ! this should have no effect if dpdt_factor .le. 0
       p0_minus_peosbar = ZERO

    end if

    if (dm .eq. 3) then
       do n=1,nlevs
          call multifab_build(w0_force_cart(n),mla%la(n),dm,1)
          call setval(w0_force_cart(n),ZERO,all=.true.)
          do comp=1,dm
             call multifab_build_edge(w0mac(n,comp),mla%la(n),1,1,comp)
             call setval(w0mac(n,comp),ZERO,all=.true.)
          end do
       end do
    end if

!    if (evolve_base_state) then
!
!       call average(mla,Source_nph,Sbar,dx,1)
!
!       call make_w0(w0,w0_old,w0_force,Sbar,rho0_old,rho0_old,p0_old,p0_old,gamma1bar_old, &
!                    gamma1bar_old,p0_minus_peosbar,psi,etarho_ec,etarho_cc,dt,dtold)
!
!       if (spherical .eq. 1) then
!          call make_w0mac(mla,w0,w0mac,dx,the_bc_tower%bc_tower_array)
!       end if
!
!       if (dm .eq. 3) then
!          call put_1d_array_on_cart(w0_force,w0_force_cart,foextrap_comp,.false., &
!                                    .true.,dx,the_bc_tower%bc_tower_array,mla)
!       end if
!
!    else

       ! these should have no effect if evolve_base_state = F
       w0_force = ZERO
       Sbar = ZERO

!    end if

    do n=1,nlevs
       do comp=1,dm
          call multifab_build_edge(umac(n,comp), mla%la(n),1,1,comp)
       end do
    end do
    
    call advance_premac(uold,sold,umac,gpi,normal,w0,w0mac,w0_force,w0_force_cart, &
                        rho0_old,grav_cell_old,dx,dt,the_bc_tower%bc_tower_array,mla)

    if (dm .eq. 3) then
       do n=1,nlevs
          call destroy(w0_force_cart(n))
       end do
    end if

    do n=1,nlevs
       call multifab_build(macrhs(n),            mla%la(n), 1, 0)
       call multifab_build(delta_gamma1_term(n), mla%la(n), 1, 0)
       call setval(delta_gamma1_term(n), ZERO, all=.true.)
    end do

    if (barrier_timers) call parallel_barrier()
    advect_time = advect_time + parallel_wtime() - advect_time_start

    macproj_time_start = parallel_wtime()

    call make_macrhs(macrhs,rho0_old,Source_nph,delta_gamma1_term,Sbar,div_coeff_old,dx, &
                     gamma1bar_old,gamma1bar_old,p0_old,p0_old,delta_p_term,dt)

    do n=1,nlevs
       call destroy(delta_gamma1_term(n))
       call destroy(Source_nph(n))
       call destroy(delta_p_term(n))
    end do

    do n=1,nlevs
       call multifab_build(macphi(n), mla%la(n), 1, 1)
       call setval(macphi(n), ZERO, all=.true.)
    end do

    ! MAC projection
    if (spherical .eq. 1) then
       do n=1,nlevs
          do comp=1,dm
             call multifab_build_edge(div_coeff_cart_edge(n,comp), mla%la(n),1,1,comp)
          end do
       end do

       call make_s0mac(mla,div_coeff_old,div_coeff_cart_edge,dx,foextrap_comp, &
                       the_bc_tower%bc_tower_array)

       call macproject(mla,umac,macphi,sold,dx,the_bc_tower,macrhs, &
                       div_coeff_cart_edge=div_coeff_cart_edge)

       do n=1,nlevs
          do comp=1,dm
             call destroy(div_coeff_cart_edge(n,comp))
          end do
       end do
    else
       call cell_to_edge(div_coeff_old,div_coeff_edge)
       call macproject(mla,umac,macphi,sold,dx,the_bc_tower, &
                       macrhs,div_coeff_1d=div_coeff_old,div_coeff_1d_edge=div_coeff_edge)
    end if

    do n=1,nlevs
       call destroy(macrhs(n))
    end do

    if (barrier_timers) call parallel_barrier()
    macproj_time = macproj_time + (parallel_wtime() - macproj_time_start)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 2: Predictor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP 2: Predictor'
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 2A: Compute advective flux divergences
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    advect_time_start = parallel_wtime()

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP 2A: Compute advective flux divergences'
    end if
    
!    if (evolve_base_state) then
!       call advect_base_dens(w0,rho0_old,rho0_new,rho0_predicted_edge,dt)
!       call compute_cutoff_coords(rho0_new)
!    else
       rho0_new = rho0_old
!    end if

    do n=1,nlevs
       call multifab_build(diff_old(n), mla%la(n), 1, 0)
    end do
    
    ! diff_old is the forcing for rhoh or temperature
    if(use_thermal_diffusion) then

       do n=1,nlevs
          call multifab_build(Tcoeff_old(n),  mla%la(n), 1,     1)
          call multifab_build(hcoeff_old(n),  mla%la(n), 1,     1)
          call multifab_build(Xkcoeff_old(n), mla%la(n), nspec, 1)
          call multifab_build(pcoeff_old(n),  mla%la(n), 1,     1)
       end do

       ! compute transport coefficients
       call make_thermal_coeffs(sold,Tcoeff_old,hcoeff_old,Xkcoeff_old,pcoeff_old)

       ! compute diff_old
       call make_explicit_thermal(mla,dx,diff_old,sold,Tcoeff_old,hcoeff_old, &
                                  Xkcoeff_old,pcoeff_old,p0_old,the_bc_tower)

    else

       do n=1,nlevs
          call setval(diff_old(n),ZERO,all=.true.)
       end do

    end if

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '            :  density and tracer_advance '
    end if

    do n=1,nlevs
       call multifab_build(shat(n),mla%la(n),nscal,nghost(sold(n)))
       do comp = 1,dm
          call multifab_build_edge(sedge(n,comp),mla%la(n),nscal,0,comp)
          call multifab_build_edge(sflux(n,comp),mla%la(n),nscal,0,comp)
       end do
       call multifab_build(scal_force(n), mla%la(n), nscal, 1)
       call multifab_build_edge(etarhoflux(n), mla%la(n), 1, 0, dm)
       call setval(etarhoflux(n),ZERO,all=.true.)
    end do

    call density_advance(mla,1,sold,shat,sedge,sflux,scal_force,intra,umac, &
                         w0,w0mac,etarhoflux, &
                         rho0_old,rho0_new,p0_new,rho0_predicted_edge, &
                         dx,dt,the_bc_tower%bc_tower_array)

    ! Now compute the new etarho
!    if (evolve_base_state) then
!       if (use_etarho) then
!
!          if (spherical .eq. 0) then
!             call make_etarho_planar(etarho_ec,etarho_cc,etarhoflux,mla)
!          else
!             call make_etarho_spherical(s1,s2,umac,w0mac,rho0_old,rho0_new,dx,normal, &
!                                        etarho_ec,etarho_cc,mla,the_bc_tower%bc_tower_array)
!          endif
!
!       endif
!    end if

    ! Correct the base state by "averaging"
!    if (use_etarho .and. evolve_base_state) then
!       call average(mla,s2,rho0_new,dx,rho_comp)
!       call compute_cutoff_coords(rho0_new)
!    end if

!    if (evolve_base_state) then
!       call make_grav_cell(grav_cell_new,rho0_new)
!    else
       grav_cell_new = grav_cell_old
       grav_cell_nph = grav_cell_old
       rho0_nph = rho0_old
!    end if

!    if (evolve_base_state) then
!
!       ! set new p0 through HSE
!       p0_new = p0_old
!       call enforce_HSE(rho0_new,p0_new,grav_cell_new)
!
!       ! make psi
!       if (spherical .eq. 0) then
!          call make_psi_planar(etarho_cc,psi)
!       else
!          ! compute p0_nph
!          p0_nph = HALF*(p0_old+p0_new)
!
!          do n=1,nlevs
!             call multifab_build(gamma1(n), mla%la(n), 1, 0)
!          end do
!
!          ! compute gamma1bar^{(1)} and store it in gamma1bar_temp1
!          call make_gamma(mla,gamma1,s1,p0_old,dx)
!          call average(mla,gamma1,gamma1bar_temp1,dx,1)
!
!          ! compute gamma1bar^{(2),*} and store it in gamma1bar_temp2
!          call make_gamma(mla,gamma1,s2,p0_new,dx)
!          call average(mla,gamma1,gamma1bar_temp2,dx,1)
!
!          do n=1,nlevs
!             call destroy(gamma1(n))
!          end do
!
!          ! compute gamma1bar^{nph,*} and store it in gamma1bar_temp2
!          gamma1bar_temp2 = HALF*(gamma1bar_temp1+gamma1bar_temp2)
!
!          ! make time-centered psi
!          call make_psi_spherical(psi,w0,gamma1bar_temp2,p0_nph,Sbar)
!       end if
!
!    else
       p0_new = p0_old
       p0_nph = p0_old
!    end if

!    if (evolve_base_state) then
!
!       ! compute rhoh0_old by "averaging"
!       call average(mla,s1,rhoh0_old,dx,rhoh_comp)
!
!       call advect_base_enthalpy(w0,rho0_old,rhoh0_old,rhoh0_new,rho0_predicted_edge,psi,dt)
!    else
       rhoh0_new = rhoh0_old
!    end if

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '            : enthalpy_advance '
    end if

    call enthalpy_advance(mla,1,uold,sold,shat,sedge,sflux,scal_force,intra, &
                          diff_old,umac,w0,w0mac, &
                          rho0_old,rhoh0_old,rho0_new,rhoh0_new,p0_old,p0_new, &
                          tempbar,psi,dx,dt,the_bc_tower%bc_tower_array)

    do n = 1, nlevs
       do comp = 1,dm
          call destroy(sedge(n,comp))
          call destroy(sflux(n,comp))
       end do
       call destroy(scal_force(n))
    end do

    ! extract aofs = (shat - sold) / dt
    do n=1,nlevs
       call multifab_build(aofs(n), mla%la(n), nscal, 0)
       call multifab_copy_c(aofs(n), 1, shat(n), 1, nscal, 0)
       call multifab_sub_sub_c(aofs(n), 1, sold(n), 1, nscal, 0)
       call multifab_div_div_s_c(aofs(n), 1, dt, nscal, 0)
    end do

    if (barrier_timers) call parallel_barrier()
    advect_time = advect_time + parallel_wtime() - advect_time_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 2B: Compute diffusive flux divergence
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    thermal_time_start = parallel_wtime()

    do n=1,nlevs
       call multifab_build(diff_hat(n), mla%la(n), 1, 0)
       call setval(diff_hat(n),ZERO,all=.true.)
    end do

    if (use_thermal_diffusion) then
       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP 2B: Compute diffusive flux divergence'
       end if

       call thermal_conduct_predictor(mla,dx,dt,sold,shat,p0_old,p0_new, &
                                      hcoeff_old,Xkcoeff_old,pcoeff_old, &
                                      aofs,intra,the_bc_tower)

       ! compute diff_hat using shat, p0_new, and old coefficients
       call make_explicit_thermal(mla,dx,diff_hat,shat, &
                                  Tcoeff_old,hcoeff_old,Xkcoeff_old,pcoeff_old, &
                                  p0_new,the_bc_tower)
    end if

    if (barrier_timers) call parallel_barrier()
    thermal_time = thermal_time + (parallel_wtime() - thermal_time_start)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 2C -- Advance thermodynamic variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    react_time_start = parallel_wtime()

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP 2C : Advance thermodynamic variables '
    end if

    ! build sdc_source
    do n=1,nlevs
       call multifab_build(sdc_source(n), mla%la(n), nscal, 0)
       call multifab_setval_c(sdc_source(n), 0.d0, 1, nscal, all=.true.)
       call multifab_plus_plus_c(sdc_source(n), rhoh_comp, diff_old(n), 1, 1, 0)
       call multifab_plus_plus_c(sdc_source(n), rhoh_comp, diff_hat(n), 1, 1, 0)
       call multifab_mult_mult_s_c(sdc_source(n), rhoh_comp, 0.5d0, 1, 0)
       call multifab_plus_plus_c(sdc_source(n), 1, aofs(n), 1, nscal, 0)
    end do

    call react_state(mla,tempbar_init,sold,snew,rho_Hext,p0_new, &
                     dt,dx,sdc_source,the_bc_tower%bc_tower_array)

    ! extract IR = [ (snew - sold)/dt - sdc_source ] 
    do n=1,nlevs
       call multifab_setval       (intra(n), 0.d0, all=.true.)
       call multifab_copy_c       (intra(n), 1, snew(n), 1,       nscal, 0)
       call multifab_sub_sub_c    (intra(n), 1, sold(n), 1,       nscal, 0)
       call multifab_div_div_s_c  (intra(n), 1, dt,               nscal, 0)
       call multifab_sub_sub_c    (intra(n), 1, sdc_source(n), 1, nscal, 0)
    end do

    ! massage the rhoh intra term into the proper form, depending on
    ! what we are predicting.  Note: we do this before we deal with
    ! the species terms, since some enthalpy types need this default
    ! species intra.

    ! first create rhohalf -- a lot of forms need this.
    do n=1,nlevs
       call multifab_build(rhohalf(n), mla%la(n), 1, 1)
    end do
    call make_at_halftime(rhohalf,sold,snew,rho_comp,1, &
                          the_bc_tower%bc_tower_array,mla)

    if (enthalpy_pred_type == predict_rhohprime) then

       call bl_error("enthalpy_pred_type = predict_rhohprime not supported")

    else if (enthalpy_pred_type == predict_h) then

       ! we want this in terms of h, not (rho h)
       do n=1,nlevs
          call multifab_div_div_c(intra(n),rhoh_comp,rhohalf(n),1,1,1)
       end do

    else if ((enthalpy_pred_type == predict_T_then_rhohprime) .or. &
             (enthalpy_pred_type == predict_T_then_h)) then

       ! for predict_T_*, the intra force needs to be in the temp_comp
       ! slot, since temperature is what is predicted.

       ! first make the thermodynamic coefficients at the half-time
       do n=1,nlevs
          call multifab_build(cphalf(n), mla%la(n), 1, 1)
          call multifab_build(xihalf(n), mla%la(n), nspec, 1)
       enddo

       call make_intra_coeffs(sold,snew,cphalf,xihalf)

       ! overwrite intra(temp_comp).  We want to create
       ! I_T = (1 / (rho c_p)) [ (rhoh_new - rhoh_old)/dt - A_rhoh -
       !     sum_k xi_k ( (rhoX_new - rhoX_old)/dt - A_rhoX ) ]
       do n=1,nlevs
          call multifab_copy_c(intra(n), temp_comp, intra(n), rhoh_comp, 1, 1)
          do comp=1, nspec
             ! multiple xi by intra and store in xi
             call multifab_mult_mult_c(xihalf(n), comp, &
                                       intra(n),  spec_comp+comp-1, 1, 1)

             ! subtract from intra temp
             call multifab_sub_sub_c(intra(n), temp_comp, xihalf(n), comp, 1, 1)
             
          enddo

          call multifab_div_div_c(intra(n), temp_comp, rhohalf(n), 1, 1, 1)
          call multifab_div_div_c(intra(n), temp_comp, cphalf(n),  1, 1, 1)

       end do

       ! clean-up
       do n=1,nlevs
          call destroy(cphalf(n))
          call destroy(xihalf(n))
       enddo

    endif     

    ! for some species_pred_types, we need to make intra in terms of
    ! X, NOT rhoX
    if ( (species_pred_type == predict_rhoprime_and_X) .or. &
         (species_pred_type == predict_rho_and_X) ) then

       do n=1,nlevs
          do comp=spec_comp,spec_comp+nspec-1
             call multifab_div_div_c(intra(n),comp,rhohalf(n),1,1,1)
          end do
       end do

    endif

    do n=1,nlevs
       call destroy(rhohalf(n))
    end do

    if (barrier_timers) call parallel_barrier()
    react_time = react_time + parallel_wtime() - react_time_start
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 3 -- Update advection velocities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (use_thermal_diffusion) then

       do n=1,nlevs
          call multifab_build(Tcoeff_new(n),  mla%la(n), 1,     1)
          call multifab_build(hcoeff_new(n),  mla%la(n), 1,     1)
          call multifab_build(Xkcoeff_new(n), mla%la(n), nspec, 1)
          call multifab_build(pcoeff_new(n),  mla%la(n), 1,     1)
          call multifab_build(diff_hterm_new(n), mla%la(n), 1, 0)
       end do
    end if

    misc_time_start = parallel_wtime()

    ! compute gamma1bar
!    if (evolve_base_state) then
!
!       do n=1,nlevs
!          call multifab_build(gamma1(n), mla%la(n), 1, 0)
!       end do
!       
!       call make_gamma(mla,gamma1,snew,p0_new,dx)
!       call average(mla,gamma1,gamma1bar,dx,1)
!
!       do n=1,nlevs
!          call destroy(gamma1(n))
!       end do
!
!       call make_div_coeff(div_coeff_new,rho0_new,p0_new,gamma1bar,grav_cell_new)
!
!    else
        
    ! Just copy div_coeff_new from div_coeff_old if not evolving the base state
    div_coeff_new = div_coeff_old

!    end if

    div_coeff_nph = HALF*(div_coeff_old + div_coeff_new)

    if (barrier_timers) call parallel_barrier()
    misc_time = misc_time + parallel_wtime() - misc_time_start

    if (sdc_iters .ge. 1) then
    
       advect_time_start = parallel_wtime()

       ! reset cutoff coordinates to old time value
       call compute_cutoff_coords(rho0_old)

       if (use_thermal_diffusion) then

          call make_thermal_coeffs(snew,Tcoeff_new,hcoeff_new,Xkcoeff_new,pcoeff_new)

          ! compute diff_new using snew, p0_new, and new coefficients
          call make_explicit_thermal(mla,dx,diff_new,snew, &
                                     Tcoeff_new,hcoeff_new,Xkcoeff_new,pcoeff_new, &
                                     p0_new,the_bc_tower)

          ! compute only the h term in diff_new
          call make_explicit_thermal_hterm(mla,dx,diff_hterm_new,snew,hcoeff_new,the_bc_tower)
          
       else
          
          do n=1,nlevs
             call setval(diff_new(n),ZERO,all=.true.)
          end do
          
       end if
       
       do n=1,nlevs
          call multifab_build(delta_gamma1_term(n), mla%la(n), 1, 0)
          call multifab_build(delta_gamma1(n), mla%la(n), 1, 0)
       end do
       
       call instantaneous_reaction_rates(mla,snew,rho_omegadot2,rho_Hnuc2)
       
       call make_S(Source_new,delta_gamma1_term,delta_gamma1,snew,uold,rho_omegadot2, &
                   rho_Hnuc2,rho_Hext,diff_new,p0_old,gamma1bar,delta_gamma1_termbar,psi,dx, &
                   mla,the_bc_tower%bc_tower_array)
       
       do n=1,nlevs
          call destroy(delta_gamma1(n))
       end do
       
       do n=1,nlevs
          call multifab_build(Source_nph(n), mla%la(n), 1, 1)
       end do
       
       call make_S_at_halftime(mla,Source_nph,Source_old,Source_new,the_bc_tower%bc_tower_array)
       
       do n=1,nlevs
          call multifab_build(delta_p_term(n), mla%la(n), 1, 0)
          call setval(delta_p_term(n),ZERO,all=.true.)
       end do
       
       ! compute p0_minus_peosbar = p0_nph - peosbar (for making w0)
       ! and delta_p_term = peos_nph - peosbar_cart (for RHS of projection)
       if (dpdt_factor .gt. ZERO) then

          do n=1,nlevs
             call multifab_build(peos_new(n), mla%la(n), 1, 0)
          end do

          ! peos_new now holds the thermodynamic p computed from snew(rho h X)
          call makePfromRhoH(snew,snew,peos_new,mla,the_bc_tower%bc_tower_array)

          do n=1,nlevs
             call multifab_build(peos_nph(n), mla%la(n), 1, 0)
          end do

          ! compute peos_nph = (1/2)*(peos_old+peos_new)
          do n=1,nlevs
             call multifab_copy(peos_nph(n), peos_old(n))
             call multifab_plus_plus(peos_nph(n), peos_new(n))
             call multifab_div_div_s(peos_nph(n), TWO)
          enddo

          do n=1,nlevs
             call destroy(peos_old(n))
             call destroy(peos_new(n))
          end do

          ! compute peosbar = Avg(peos_nph)
          call average(mla,peos_nph,peosbar,dx,1)

          ! compute p0_nph = (1/2)*(p0_old+p0_new)
          p0_nph = HALF*(p0_old + p0_new)

          ! compute p0_minus_peosbar = p0_nph - peosbar
          p0_minus_peosbar = p0_nph - peosbar

          do n=1,nlevs
             call multifab_build(peosbar_cart(n), mla%la(n), 1, 0)
          end do

          ! compute peosbar_cart from peosbar
          call put_1d_array_on_cart(peosbar,peosbar_cart,foextrap_comp, &
               .false.,.false.,dx,the_bc_tower%bc_tower_array,mla)

          ! compute delta_p_term = peos_nph - peosbar_cart
          do n=1,nlevs
             call multifab_copy(delta_p_term(n), peos_nph(n))
             call multifab_sub_sub(delta_p_term(n), peosbar_cart(n))
          end do

          do n=1,nlevs
             call destroy(peos_nph(n))
             call destroy(peosbar_cart(n))
          end do

       end if

       if (dm .eq. 3) then
          do n=1,nlevs
             call multifab_build(w0_force_cart(n), mla%la(n), dm, 1)
             call setval(w0_force_cart(n),ZERO,all=.true.)
          end do
       end if

       !    if (evolve_base_state) then
       !
       !       call average(mla,Source_nph,Sbar,dx,1)
       !
       !       if(use_delta_gamma1_term) then
       !          ! add delta_gamma1_termbar to Sbar
       !          Sbar = Sbar + delta_gamma1_termbar
       !       end if
       !
       !       call make_w0(w0,w0_old,w0_force,Sbar,rho0_old,rho0_new,p0_old,p0_new, &
       !                    gamma1bar_old,gamma1bar,p0_minus_peosbar, &
       !                    psi,etarho_ec,etarho_cc,dt,dtold)
       !
       !       if (spherical .eq. 1) then
       !          call make_w0mac(mla,w0,w0mac,dx,the_bc_tower%bc_tower_array)
       !       end if
       !
       !       if (dm .eq. 3) then
       !          call put_1d_array_on_cart(w0_force,w0_force_cart,foextrap_comp,.false., &
       !                                    .true.,dx,the_bc_tower%bc_tower_array,mla)
       !       end if
       !    end if

       call advance_premac(uold,sold,umac,gpi,normal,w0,w0mac,w0_force,w0_force_cart, &
                           rho0_old,grav_cell_old,dx,dt,the_bc_tower%bc_tower_array,mla)

       if (barrier_timers) call parallel_barrier()
       advect_time = advect_time + parallel_wtime() - advect_time_start

       macproj_time_start = parallel_wtime()

       do n=1,nlevs
          call multifab_build(macrhs(n), mla%la(n), 1, 0)
       end do

       ! note delta_gamma1_term here is not time-centered
       call make_macrhs(macrhs,rho0_old,Source_nph,delta_gamma1_term,Sbar,div_coeff_nph,dx, &
                        gamma1bar_old,gamma1bar,p0_old,p0_new,delta_p_term,dt)

       do n=1,nlevs
          call destroy(delta_gamma1_term(n))
          call destroy(Source_nph(n))
          call destroy(delta_p_term(n))
       end do

       do n=1,nlevs
          call multifab_build(rhohalf(n), mla%la(n), 1, 1)
       end do

       call make_at_halftime(rhohalf,sold,snew,rho_comp,1,the_bc_tower%bc_tower_array,mla)

       ! MAC projection !
       if (spherical .eq. 1) then
          do n=1,nlevs
             do comp=1,dm
                call multifab_build_edge(div_coeff_cart_edge(n,comp), mla%la(n),1,1,comp)
             end do
          end do

          call make_s0mac(mla,div_coeff_nph,div_coeff_cart_edge,dx,foextrap_comp, &
                          the_bc_tower%bc_tower_array)

          call macproject(mla,umac,macphi,rhohalf,dx,the_bc_tower,macrhs, &
                          div_coeff_cart_edge=div_coeff_cart_edge)

          do n=1,nlevs
             do comp=1,dm
                call destroy(div_coeff_cart_edge(n,comp))
             end do
          end do
       else
          call cell_to_edge(div_coeff_nph,div_coeff_edge)
          call macproject(mla,umac,macphi,rhohalf,dx,the_bc_tower,macrhs, &
                          div_coeff_1d=div_coeff_nph,div_coeff_1d_edge=div_coeff_edge)
       end if

       ! advect the particles through dt using umac.  Then redistribute them
       if (.not. init_mode .and. use_particles) then

          call addw0(umac,the_bc_tower%bc_tower_array,mla,w0,w0mac,mult=1.d0)

          call move_advect(particles,mla,umac,dx,dt,prob_lo,prob_hi)

          call addw0(umac,the_bc_tower%bc_tower_array,mla,w0,w0mac,mult=-1.d0)

       end if

       do n=1,nlevs
          call destroy(rhohalf(n))
          call destroy(macrhs(n))
          call destroy(macphi(n))
       end do

       if (barrier_timers) call parallel_barrier()
       macproj_time = macproj_time + (parallel_wtime() - macproj_time_start)

    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 4: Corrector loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do misdc=1,sdc_iters

       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP 4: Corrector loop (MISDC iter = ', misdc, ')   '
       end if

       advect_time_start = parallel_wtime()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 4A: Compute advective flux divergences
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP 4A: Compute advective flux divergences (MISDC iter = ', misdc, ')   '
       end if

!    if (evolve_base_state) then
!       call advect_base_dens(w0,rho0_old,rho0_new,rho0_predicted_edge,dt)
!       call compute_cutoff_coords(rho0_new)
!    else
       rho0_new = rho0_old
!    end if

       do n=1,nlevs
          call setval(etarhoflux(n),ZERO,all=.true.)
       end do
       
       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '            :  density and tracer advance'
       end if

       ! Build the sedge array.
       do n=1,nlevs
          do comp = 1,dm
             call multifab_build_edge(sedge(n,comp),mla%la(n),nscal,0,comp)
             call multifab_build_edge(sflux(n,comp),mla%la(n),nscal,0,comp)
          end do
          call multifab_build(scal_force(n), mla%la(n), nscal, 1)
       end do
       
       call density_advance(mla,2,sold,shat,sedge,sflux,scal_force,intra,umac, &
                            w0,w0mac,etarhoflux, &
                            rho0_old,rho0_new,p0_new,rho0_predicted_edge,dx,dt, &
                            the_bc_tower%bc_tower_array)

    ! Now compute the new etarho
!    if (evolve_base_state) then
!       if (use_etarho) then
!
!          if (spherical .eq. 0) then
!             call make_etarho_planar(etarho_ec,etarho_cc,etarhoflux,mla)
!          else
!             call make_etarho_spherical(s1,s2,umac,w0mac,rho0_old,rho0_new,dx,normal, &
!                                        etarho_ec,etarho_cc,mla,the_bc_tower%bc_tower_array)
!          endif
!
!       endif
!    end if

    ! Correct the base state using "averaging"
!    if (use_etarho .and. evolve_base_state) then
!       call average(mla,s2,rho0_new,dx,rho_comp)
!       call compute_cutoff_coords(rho0_new)
!    end if

!    if (evolve_base_state) then
!       call make_grav_cell(grav_cell_new,rho0_new)
!       rho0_nph = HALF*(rho0_old+rho0_new)
!       call make_grav_cell(grav_cell_nph,rho0_nph)
!    else
       grav_cell_new = grav_cell_old
       grav_cell_nph = grav_cell_old
       rho0_nph = rho0_old
!    end if

!    if (evolve_base_state) then
!       
!       ! set new p0 through HSE
!       p0_new = p0_old
!       call enforce_HSE(rho0_new,p0_new,grav_cell_new)
!       p0_nph = HALF*(p0_old+p0_new)
!
!       ! make psi
!       if (spherical .eq. 0) then
!          call make_psi_planar(etarho_cc,psi)
!       else
!          p0_nph = HALF*(p0_old+p0_new)
!
!          do n=1,nlevs
!             call multifab_build(gamma1(n), mla%la(n), 1, 0)
!          end do
!
!          ! compute gamma1bar^{(2)} and store it in gamma1bar_temp2
!          call make_gamma(mla,gamma1,s2,p0_new,dx)
!          call average(mla,gamma1,gamma1bar_temp2,dx,1)
!
!          do n=1,nlevs
!             call destroy(gamma1(n))
!          end do
!
!          ! compute gamma1bar^{nph} and store it in gamma1bar_temp2
!          gamma1bar_temp2 = HALF*(gamma1bar_temp1+gamma1bar_temp2)
!
!          call make_psi_spherical(psi,w0,gamma1bar_temp2,p0_nph,Sbar)
!       end if
!
!    else
       p0_new = p0_old
       p0_nph = p0_old
 !   end if

!    if (evolve_base_state) then
!       call advect_base_enthalpy(w0,rho0_old,rhoh0_old,rhoh0_new, &
!                                 rho0_predicted_edge,psi,dt)
!    else
       rhoh0_new = rhoh0_old
!    end if

       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '            : enthalpy_advance'
       end if

       call enthalpy_advance(mla,2,uold,sold,shat,sedge,sflux,scal_force, &
                             intra,diff_old,umac,w0,w0mac, &
                             rho0_old,rhoh0_old,rho0_new,rhoh0_new,p0_old,p0_new, &
                             tempbar,psi,dx,dt,the_bc_tower%bc_tower_array)

       do n=1,nlevs
          do comp = 1,dm
             call destroy(sedge(n,comp))
             call destroy(sflux(n,comp))
          end do
          call destroy(scal_force(n))
       end do

       ! extract aofs = (shat - sold) / dt
       do n=1,nlevs
          call multifab_copy_c(aofs(n), 1, shat(n), 1, nscal, 0)
          call multifab_sub_sub_c(aofs(n), 1, sold(n), 1, nscal, 0)
          call multifab_div_div_s_c(aofs(n), 1, dt, nscal, 0)
       end do

       if (barrier_timers) call parallel_barrier()
       advect_time = advect_time + parallel_wtime() - advect_time_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 4B: Compute diffusive flux divergences
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP 4B: Compute diffusive flux divergence (MISDC iter = ', misdc, ')   '
       end if

       thermal_time_start = parallel_wtime()

       if (use_thermal_diffusion) then

          call thermal_conduct_corrector(mla,dx,dt,sold,shat,snew,p0_old,p0_new, &
                                         hcoeff_old,Xkcoeff_old,pcoeff_old, &
                                         hcoeff_new,Xkcoeff_new,pcoeff_new, &
                                         intra,the_bc_tower)
          
          ! compute diff_hat using shat, p0_new, and new coefficients from previous iteration
          call make_explicit_thermal(mla,dx,diff_hat,shat, &
                                     Tcoeff_new,hcoeff_new,Xkcoeff_new,pcoeff_new, &
                                     p0_new,the_bc_tower)


          ! compute only the h term in diff_hat
          do n=1,nlevs
             call multifab_build(diff_hterm_hat(n), mla%la(n), 1, 0)
             call setval(diff_hterm_hat(n),ZERO,all=.true.)
          end do
          call make_explicit_thermal_hterm(mla,dx,diff_hterm_hat,shat,hcoeff_new,the_bc_tower)
          
       end if

       if (barrier_timers) call parallel_barrier()
       thermal_time = thermal_time + (parallel_wtime() - thermal_time_start)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 4C: Advance thermodynamic variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       react_time_start = parallel_wtime()
       
       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP 4C: Advance thermodynamic variables (MISDC iter = ', misdc, ')   '
       end if

       ! build sdc_source
       do n=1,nlevs
          call multifab_setval_c(sdc_source(n), 0.d0, 1, nscal, all=.true.)
          call multifab_plus_plus_c(sdc_source(n), rhoh_comp, diff_old(n), 1, 1, 0)
          call multifab_plus_plus_c(sdc_source(n), rhoh_comp, diff_hat(n), 1, 1, 0)
          call multifab_plus_plus_c(sdc_source(n), rhoh_comp, diff_hterm_hat(n), 1, 1, 0)
          call multifab_sub_sub_c(sdc_source(n), rhoh_comp, diff_hterm_new(n), 1, 1, 0)
          call multifab_mult_mult_s_c(sdc_source(n), rhoh_comp, 0.5d0, 1, 0)
          call multifab_plus_plus_c(sdc_source(n), 1, aofs(n), 1, nscal, 0)
       end do

       do n=1,nlevs
          call multifab_destroy(diff_hterm_hat(n))
       end do

       call react_state(mla,tempbar_init,sold,snew,rho_Hext,p0_new, &
                        dt,dx,sdc_source,the_bc_tower%bc_tower_array)

       ! extract IR = [ (snew - sold)/dt - sdc_source ] 
       do n=1,nlevs
          call multifab_setval       (intra(n), 0.d0, all=.true.)
          call multifab_copy_c       (intra(n), 1, snew(n), 1,       nscal, 0)
          call multifab_sub_sub_c    (intra(n), 1, sold(n), 1,       nscal, 0)
          call multifab_div_div_s_c  (intra(n), 1, dt,               nscal, 0)
          call multifab_sub_sub_c    (intra(n), 1, sdc_source(n), 1, nscal, 0)
       end do

       ! massage the rhoh intra term into the proper form, depending on
       ! what we are predicting.  Note: we do this before we deal with
       ! the species terms, since some enthalpy types need this default
       ! species intra.
       
       ! first create rhohalf -- a lot of forms need this.
       do n=1,nlevs
          call multifab_build(rhohalf(n), mla%la(n), 1, 1)
       end do

       call make_at_halftime(rhohalf,sold,snew,rho_comp,1, &
                             the_bc_tower%bc_tower_array,mla)

       if (enthalpy_pred_type == predict_rhohprime) then
          
          call bl_error("enthalpy_pred_type = predict_rhohprime not supported")

       else if (enthalpy_pred_type == predict_h) then

          ! we want this in terms of h, not (rho h)
          do n=1,nlevs
             call multifab_div_div_c(intra(n),rhoh_comp,rhohalf(n),1,1,1)
          end do

       else if ((enthalpy_pred_type == predict_T_then_rhohprime) .or. &
                (enthalpy_pred_type == predict_T_then_h)) then

          ! for predict_T_*, the intra force needs to be in the temp_comp
          ! slot, since temperature is what is predicted.

          ! first make the thermodynamic coefficients at the half-time
          do n=1,nlevs
             call multifab_build(cphalf(n), mla%la(n), 1, 1)
             call multifab_build(xihalf(n), mla%la(n), nspec, 1)
          enddo

          call make_intra_coeffs(sold,snew,cphalf,xihalf)

          ! overwrite intra(temp_comp).  We want to create
          ! I_T = (1 / (rho c_p)) [ (rhoh_new - rhoh_old)/dt - A_rhoh -
          !     sum_k xi_k ( (rhoX_new - rhoX_old)/dt - A_rhoX ) ]
          do n=1,nlevs
             call multifab_copy_c(intra(n), temp_comp, intra(n), rhoh_comp, 1, 1)

             do comp=1,nspec

                ! multiple xi by intra and store in xi
                call multifab_mult_mult_c(xihalf(n), comp, &
                                          intra(n),  spec_comp+comp-1, 1, 1)

                ! subtract from intra temp
                call multifab_sub_sub_c(intra(n), temp_comp, xihalf(n), comp, 1, 1)
             
             enddo

             call multifab_div_div_c(intra(n), temp_comp, rhohalf(n), 1, 1, 1)
             call multifab_div_div_c(intra(n), temp_comp, cphalf(n),  1, 1, 1)

          end do
          
          ! clean-up
          do n=1,nlevs
             call destroy(cphalf(n))
             call destroy(xihalf(n))
          enddo
          
       endif

       ! for some species_pred_types, we need to make intra in terms of
       ! X, NOT rhoX
       if ( (species_pred_type == predict_rhoprime_and_X) .or. &
            (species_pred_type == predict_rho_and_X) ) then
          
          do n=1,nlevs
             do comp=spec_comp,spec_comp+nspec-1
                call multifab_div_div_c(intra(n),comp,rhohalf(n),1,1,1)
             end do
          end do
       endif
       
       do n=1,nlevs
          call destroy(rhohalf(n))
       enddo
       
       if (barrier_timers) call parallel_barrier()
       react_time = react_time + parallel_wtime() - react_time_start
       
       misc_time_start = parallel_wtime()

       if (use_thermal_diffusion) then

          call make_thermal_coeffs(snew,Tcoeff_new,hcoeff_new,Xkcoeff_new,pcoeff_new)

          ! compute diff_new using snew, p0_new, and new coefficients
          call make_explicit_thermal(mla,dx,diff_new,snew, &
                                     Tcoeff_new,hcoeff_new,Xkcoeff_new,pcoeff_new, &
                                     p0_new,the_bc_tower)

       end if
       
    ! compute gamma1bar
!    if (evolve_base_state) then
!
!       do n=1,nlevs
!          call multifab_build(gamma1(n), mla%la(n), 1, 0)
!       end do
!
!       call make_gamma(mla,gamma1,snew,p0_new,dx)
!       call average(mla,gamma1,gamma1bar,dx,1)
!
!       do n=1,nlevs
!          call destroy(gamma1(n))
!       end do
!
!       !  We used to call this even if evolve_base was false,but we don't need to
!       call make_div_coeff(div_coeff_new,rho0_new,p0_new,gamma1bar,grav_cell_new)
!
!    else
        
       ! Just copy div_coeff_new from div_coeff_old if not evolving the base state
       div_coeff_new = div_coeff_old

!    end if

       div_coeff_nph = HALF*(div_coeff_old+div_coeff_new)

       if (barrier_timers) call parallel_barrier()
       misc_time = misc_time + parallel_wtime() - misc_time_start
       
    end do ! end loop over misdc iterations

    if (use_thermal_diffusion) then
       do n=1,nlevs
          call destroy(shat(n))
          call destroy(Tcoeff_old(n))
          call destroy(hcoeff_old(n))
          call destroy(Xkcoeff_old(n))
          call destroy(pcoeff_old(n))
          call destroy(Tcoeff_new(n))
          call destroy(hcoeff_new(n))
          call destroy(Xkcoeff_new(n))
          call destroy(pcoeff_new(n))
          call destroy(diff_hterm_new(n))
       end do
    end if

    do n=1,nlevs
       call destroy(etarhoflux(n))
       call destroy(diff_old(n))
       call destroy(diff_hat(n))
       call destroy(aofs(n))
       call destroy(sdc_source(n))
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 5 -- Advance velocity and dynamic pressure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    thermal_time_start = parallel_wtime()
    
    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP 5 : Advance velocity and dynamic pressure'
    end if
    
    do n=1,nlevs
       call multifab_build(delta_gamma1_term(n), mla%la(n), 1, 0)
       call multifab_build(delta_gamma1(n), mla%la(n), 1, 0)
    end do

    thermal_time = thermal_time + (parallel_wtime() - thermal_time_start)

    ndproj_time_start = parallel_wtime()

    call instantaneous_reaction_rates(mla,snew,rho_omegadot2,rho_Hnuc2)

    call make_S(Source_new,delta_gamma1_term,delta_gamma1,snew,uold,rho_omegadot2, &
                rho_Hnuc2,rho_Hext,diff_new,p0_new,gamma1bar,delta_gamma1_termbar,psi,dx, &
                mla,the_bc_tower%bc_tower_array)

    do n=1,nlevs
       call destroy(delta_gamma1(n))
    end do

!    if (evolve_base_state) then
!       call average(mla,Source_new,Sbar,dx,1)
!
!       if(use_delta_gamma1_term) then
!          ! add delta_gamma1_termbar to Sbar
!          Sbar = Sbar + delta_gamma1_termbar
!       end if
!
!    end if
    
    ! define dSdt = (Source_new - Source_old) / dt
    do n=1,nlevs
       call multifab_copy(dSdt(n),Source_new(n))
       call multifab_sub_sub(dSdt(n),Source_old(n))
       call multifab_div_div_s(dSdt(n),dt)
    end do

    if (barrier_timers) call parallel_barrier()
    ndproj_time = ndproj_time + parallel_wtime() - ndproj_time_start
    
    advect_time_start = parallel_wtime()

    ! Define rho at half time using the new rho from Step 8
    do n=1,nlevs
       call multifab_build(rhohalf(n), mla%la(n), 1, 1)
    end do

    call make_at_halftime(rhohalf,sold,snew,rho_comp,1,the_bc_tower%bc_tower_array,mla)

    call velocity_advance(mla,uold,unew,sold,rhohalf,umac,gpi,normal,w0,w0mac,w0_force, &
                          w0_force_cart,rho0_old,rho0_nph,grav_cell_old,grav_cell_nph, &
                          dx,dt,the_bc_tower%bc_tower_array,sponge)

    if (init_mode) then
       ! throw away w0 by setting w0 = w0_old
       w0 = w0_old
    end if

    do n=1,nlevs
       do comp=1,dm
          call destroy(umac(n,comp))
       end do
    end do

    if (dm .eq. 3) then
       do n=1,nlevs
          do comp=1,dm
             call destroy(w0mac(n,comp))
          end do
          call destroy(w0_force_cart(n))
       end do
    end if

    if (barrier_timers) call parallel_barrier()
    advect_time = advect_time + parallel_wtime() - advect_time_start
       
    ndproj_time_start = parallel_wtime()

    ! Project the new velocity field.
    if (init_mode) then

       proj_type = pressure_iters_comp

       do n=1,nlevs
          call multifab_build(hgrhs_old(n), mla%la(n), 1, 0, nodal)
          call multifab_copy(hgrhs_old(n),hgrhs(n))
       end do
       call make_hgrhs(the_bc_tower,mla,hgrhs,Source_new,delta_gamma1_term, &
                       Sbar,div_coeff_nph,dx)
       do n=1,nlevs
          call multifab_sub_sub(hgrhs(n),hgrhs_old(n))
          call multifab_div_div_s(hgrhs(n),dt)
       end do

    else

       proj_type = regular_timestep_comp

       call make_hgrhs(the_bc_tower,mla,hgrhs,Source_new,delta_gamma1_term, &
                       Sbar,div_coeff_nph,dx)

       ! compute delta_p_term = peos_new - peosbar_cart (for RHS of projection)
       if (dpdt_factor .gt. ZERO) then

          do n=1,nlevs
             call multifab_build(peos_new(n), mla%la(n), 1, 0)
          enddo

          ! peos_new now holds the thermodynamic p computed from snew(rho h X)
          call makePfromRhoH(snew,snew,peos_new,mla,the_bc_tower%bc_tower_array)

          ! compute peosbar = Avg(peos_new)
          call average(mla,peos_new,peosbar,dx,1)

          ! no need to compute p0_minus_peosbar since make_w0 is not called after here

          do n=1,nlevs
             call multifab_build(peosbar_cart(n), mla%la(n), 1, 0)
          end do

          ! compute peosbar_cart from peosbar
          call put_1d_array_on_cart(peosbar,peosbar_cart,foextrap_comp, &
                                    .false.,.false.,dx,the_bc_tower%bc_tower_array,mla)

          do n=1,nlevs
             call multifab_build(delta_p_term(n), mla%la(n), 1, 0)
          end do

          ! compute delta_p_term = peos_new - peosbar_cart
          do n=1,nlevs
             call multifab_copy(delta_p_term(n), peos_new(n))
             call multifab_sub_sub(delta_p_term(n), peosbar_cart(n))
          end do

          do n=1,nlevs
             call destroy(peos_new(n))
             call destroy(peosbar_cart(n))
          end do
          
          call correct_hgrhs(the_bc_tower,mla,rho0_new,hgrhs,div_coeff_nph,dx,dt, &
                             gamma1bar,p0_new,delta_p_term)
          
          do n=1,nlevs
             call destroy(delta_p_term(n))
          enddo
          
       end if

    end if

    do n=1,nlevs
       call destroy(delta_gamma1_term(n))
    end do

    do n=1,nlevs
       call multifab_build(div_coeff_3d(n), mla%la(n), 1, 1)
    end do
       
    call put_1d_array_on_cart(div_coeff_nph,div_coeff_3d,foextrap_comp,.false., &
                              .false.,dx,the_bc_tower%bc_tower_array,mla)

    call hgproject(proj_type,mla,unew,uold,rhohalf,pi,gpi,dx,dt,the_bc_tower,div_coeff_3d,hgrhs)

    do n=1,nlevs
       call destroy(div_coeff_3d(n))
       call destroy(rhohalf(n))
    end do
    
    ! If doing pressure iterations then put hgrhs_old into hgrhs to be returned to varden.
    if (init_mode) then
       do n=1,nlevs
          call multifab_copy(hgrhs(n),hgrhs_old(n))
          call destroy(hgrhs_old(n))
       end do
    end if

    if (barrier_timers) call parallel_barrier()
    ndproj_time = ndproj_time + (parallel_wtime() - ndproj_time_start)

    misc_time_start = parallel_wtime()

    if (.not. init_mode) then
       
       grav_cell_old = grav_cell_new

       if (.not. fix_base_state) then
          ! compute tempbar by "averaging"
          call average(mla,snew,tempbar,dx,temp_comp)
       end if

       ! output any runtime diagnostics
       ! pass in the new time value, time+dt
       call diag(time+dt,dt,dx,snew,rho_Hnuc2,rho_Hext,diff_new,rho_omegadot2,&
                 rho0_new,rhoh0_new,p0_new,tempbar, &
                 gamma1bar,div_coeff_new, &
                 unew,w0,normal, &
                 mla,the_bc_tower)


       ! perform sanity checks, if desired
       if (mach_max_abort > ZERO) then
          call sanity_check(time+dt,dx,snew, &
                 rho0_new,rhoh0_new,p0_new,tempbar, &
                 gamma1bar,div_coeff_new, &
                 unew,w0,normal, &
                 mla,the_bc_tower)
       endif

    end if

    if (barrier_timers) call parallel_barrier()
    misc_time = misc_time + parallel_wtime() - misc_time_start

    call destroy(bpt)

    call parallel_reduce(advect_time_max, advect_time, MPI_MAX, &
                         proc=parallel_IOProcessorNode())

    call parallel_reduce(macproj_time_max,   macproj_time, MPI_MAX, &
                         proc=parallel_IOProcessorNode())

    call parallel_reduce(ndproj_time_max,   ndproj_time, MPI_MAX, &
                         proc=parallel_IOProcessorNode())

    call parallel_reduce(react_time_max,  react_time, MPI_MAX, &
                         proc=parallel_IOProcessorNode())

    call parallel_reduce(misc_time_max,  misc_time, MPI_MAX, &
                         proc=parallel_IOProcessorNode())

    if(use_thermal_diffusion) then
      call parallel_reduce(thermal_time_max,  thermal_time, MPI_MAX, &
                           proc=parallel_IOProcessorNode())
    end if

    if (parallel_IOProcessor()) then
       write(6,*) 'Timing summary:'
       write(6,*) '   Advection       : ', advect_time_max , ' seconds'
       write(6,*) '   MAC   Projection: ', macproj_time_max, ' seconds'
       write(6,*) '   Nodal Projection: ', ndproj_time_max , ' seconds'
       if (use_thermal_diffusion) &
          write(6,*) '   Thermal         : ', thermal_time_max, ' seconds'
       write(6,*) '   Reactions       : ', react_time_max  , ' seconds'
       write(6,*) '   Misc            : ', misc_time_max   , ' seconds'
       write(6,*) ' '
    endif
    
  end subroutine advance_timestep

end module advance_timestep_module
