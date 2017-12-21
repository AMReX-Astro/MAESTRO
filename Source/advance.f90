! advance_timestep advances the solution through one timestep,
! proceeding through the 12 steps described in the multilevel paper
! (Nonaka et al. 2010).

module advance_timestep_module

  use bl_types           , only: dp_t
  use bl_constants_module, only: ZERO, HALF, TWO
  use multifab_module !   , only: multifab, multifab_build, multifab_build_edge, &
                      !           destroy, setval, nghost, &
                      !           extent, multifab_volume, nboxes, &
                      !           multifab_copy, multifab_copy_c, &
                      !           multifab_sub_sub, multifab_div_div_s, multifab_plus_plus
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
                              rho0_new,rhoh0_new,p0_old,p0_new,tempbar,gamma1bar_old, &
                              gamma1bar_new,w0, &
                              rho_omegadot2,rho_Hnuc2,rho_Hext,thermal2,&
                              beta0_old,beta0_new, &
                              grav_cell_old,grav_cell_new,dx,dt,dtold,the_bc_tower, &
                              dSdt,S_cc_old,S_cc_new,etarho_ec,etarho_cc, &
                              psi,sponge,nodalrhs,tempbar_init,particles)

    use bl_prof_module              , only : bl_prof_timer, build, destroy
    use      pre_advance_module     , only : advance_premac
    use velocity_advance_module     , only : velocity_advance
    use  density_advance_module     , only : density_advance
    use enthalpy_advance_module     , only : enthalpy_advance
    use make_beta0_module           , only : make_beta0
    use make_w0_module              , only : make_w0
    use advect_base_module          , only : advect_base_dens, advect_base_enthalpy
    use react_state_module          , only : react_state
    use make_S_cc_module            , only : make_S_cc
    use average_module              , only : average
    use phihalf_module              , only : make_S_at_halftime, make_at_halftime
    use extraphalf_module           , only : extrap_to_halftime
    use thermal_conduct_module      , only : thermal_conduct
    use make_explicit_thermal_module, only : make_explicit_thermal, make_thermal_coeffs 
    use make_grav_module            , only : make_grav_cell
    use make_eta_module             , only : make_etarho_planar, make_etarho_spherical
    use make_psi_module             , only : make_psi_planar, make_psi_spherical
    use fill_3d_module              , only : put_1d_array_on_cart, make_w0mac, make_s0mac
    use cell_to_edge_module         , only : cell_to_edge
    use make_gamma_module           , only : make_gamma1bar
    use rhoh_vs_t_module            , only : makePfromRhoH, makeTfromRhoP, makeTfromRhoH
    use diag_module                 , only : diag
    use sanity_module               , only : sanity_check
    use enforce_HSE_module          , only : enforce_HSE

    use macrhs_module               , only : make_macrhs
    use macproject_module           , only : macproject

    use make_nodalrhs_module        , only : make_nodalrhs, correct_nodalrhs
    use hgproject_module            , only : hgproject
    use proj_parameters             , only : pressure_iters_comp, regular_timestep_comp

    use variables                   , only : nscal, temp_comp, rho_comp, rhoh_comp, pi_comp, &
                                             foextrap_comp
    use geometry                    , only : nlevs_radial, spherical, nr_fine, compute_cutoff_coords
    use network                     , only : nspec
    use probin_module               , only : barrier_timers, evolve_base_state, fix_base_state, &
                                             use_etarho, dpdt_factor, verbose, &
                                             use_tfromp, use_thermal_diffusion, &
                                             use_delta_gamma1_term, nodal, mach_max_abort, &
                                             prob_lo, prob_hi, use_particles, ppm_trace_forces
    use time_module                 , only : time
    use addw0_module                , only : addw0
    use make_pi_cc_module           , only : make_pi_cc
    
    logical,         intent(in   ) :: init_mode
    type(ml_layout), intent(inout) :: mla
    type(multifab),  intent(in   ) ::   uold(:)
    type(multifab),  intent(in   ) ::   sold(:)
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
    real(dp_t)    ,  intent(inout) :: gamma1bar_old(:,0:)
    real(dp_t)    ,  intent(inout) :: gamma1bar_new(:,0:)
    real(dp_t)    ,  intent(inout) ::        w0(:,0:)
    type(multifab),  intent(inout) :: rho_omegadot2(:)
    type(multifab),  intent(inout) :: rho_Hnuc2(:)
    type(multifab),  intent(inout) :: rho_Hext(:)
    type(multifab),  intent(inout) ::  thermal2(:)
    real(dp_t)    ,  intent(inout) :: beta0_old(:,0:)
    real(dp_t)    ,  intent(inout) :: beta0_new(:,0:)
    real(dp_t)    ,  intent(inout) :: grav_cell_old(:,0:)
    real(dp_t)    ,  intent(inout) :: grav_cell_new(:,0:)
    real(dp_t)    ,  intent(in   ) :: dx(:,:),dt,dtold
    type(bc_tower),  intent(in   ) :: the_bc_tower
    type(multifab),  intent(inout) ::       dSdt(:)
    type(multifab),  intent(inout) :: S_cc_old(:)
    type(multifab),  intent(inout) :: S_cc_new(:)
    real(dp_t)    ,  intent(inout) ::  etarho_ec(:,0:)
    real(dp_t)    ,  intent(inout) ::  etarho_cc(:,0:)
    real(dp_t)    ,  intent(inout) ::        psi(:,0:)
    type(multifab),  intent(in   ) :: sponge(:)
    type(multifab),  intent(inout) ::  nodalrhs(:)
    type(particle_container), intent(inout) :: particles

    ! local
    type(multifab) ::             rhohalf(mla%nlevel)
    type(multifab) ::       w0_force_cart(mla%nlevel)
    type(multifab) ::              macrhs(mla%nlevel)
    type(multifab) ::              macphi(mla%nlevel)
    type(multifab) ::         nodalrhs_old(mla%nlevel)
    type(multifab) ::            S_cc_nph(mla%nlevel)
    type(multifab) ::            thermal1(mla%nlevel)
    type(multifab) ::                  s1(mla%nlevel)
    type(multifab) ::                  s2(mla%nlevel)
    type(multifab) ::              s2star(mla%nlevel)
    type(multifab) ::   delta_gamma1_term(mla%nlevel)
    type(multifab) ::        delta_gamma1(mla%nlevel)
    type(multifab) ::      beta0_cart(mla%nlevel)
    type(multifab) ::          etarhoflux(mla%nlevel)
    type(multifab) ::        peosbar_cart(mla%nlevel)
    type(multifab) ::        delta_p_term(mla%nlevel)
    type(multifab) ::              Tcoeff(mla%nlevel)
    type(multifab) ::             hcoeff1(mla%nlevel)
    type(multifab) ::            Xkcoeff1(mla%nlevel)
    type(multifab) ::             pcoeff1(mla%nlevel)
    type(multifab) ::             hcoeff2(mla%nlevel)
    type(multifab) ::            Xkcoeff2(mla%nlevel)
    type(multifab) ::             pcoeff2(mla%nlevel)
    type(multifab) ::          scal_force(mla%nlevel)
    type(multifab) ::           delta_chi(mla%nlevel)

    type(multifab) ::               w0mac(mla%nlevel,mla%dim)
    type(multifab) ::                umac(mla%nlevel,mla%dim)
    type(multifab) ::               sedge(mla%nlevel,mla%dim)
    type(multifab) ::               sflux(mla%nlevel,mla%dim)
    type(multifab) :: beta0_cart_edge(mla%nlevel,mla%dim)

    real(kind=dp_t), allocatable ::        grav_cell_nph(:,:)
    real(kind=dp_t), allocatable ::             rho0_nph(:,:)
    real(kind=dp_t), allocatable ::               p0_nph(:,:)
    real(kind=dp_t), allocatable ::     p0_minus_peosbar(:,:)
    real(kind=dp_t), allocatable ::              peosbar(:,:)
    real(kind=dp_t), allocatable ::             w0_force(:,:)
    real(kind=dp_t), allocatable ::                 Sbar(:,:)
    real(kind=dp_t), allocatable ::        beta0_nph(:,:)
    real(kind=dp_t), allocatable ::      gamma1bar_temp1(:,:)
    real(kind=dp_t), allocatable ::      gamma1bar_temp2(:,:)
    real(kind=dp_t), allocatable :: delta_gamma1_termbar(:,:)
    real(kind=dp_t), allocatable ::               w0_old(:,:)
    real(kind=dp_t), allocatable ::       beta0_edge(:,:)
    real(kind=dp_t), allocatable ::  rho0_predicted_edge(:,:)
    real(kind=dp_t), allocatable ::         delta_chi_w0(:,:)

    integer    :: i,n,comp,proj_type,nlevs,dm

    ! need long int to store numbers greater than 2^31
    integer(kind=ll_t) :: numcell

    ! keep track of wallclock time of various parts of the code
    real(kind=dp_t) :: advect_time , advect_time_start , advect_time_max
    real(kind=dp_t) :: macproj_time, macproj_time_start, macproj_time_max
    real(kind=dp_t) :: ndproj_time , ndproj_time_start , ndproj_time_max
    real(kind=dp_t) :: thermal_time, thermal_time_start, thermal_time_max
    real(kind=dp_t) :: react_time  , react_time_start  , react_time_max
    real(kind=dp_t) :: misc_time   , misc_time_start   , misc_time_max

    integer :: nreduce
    real(kind=dp_t), allocatable :: times_local(:), times_global(:)
    
    type(bl_prof_timer), save :: bpt

    call build(bpt, "advance_timestep")

    nlevs = mla%nlevel
    dm = mla%dim

    ! keep track of wallclock time of various parts of the code
    advect_time  = 0.d0
    macproj_time = 0.d0
    ndproj_time  = 0.d0
    thermal_time = 0.d0
    react_time   = 0.d0
    misc_time    = 0.d0

    misc_time_start = parallel_wtime()

    allocate(       grav_cell_nph(nlevs_radial,0:nr_fine-1))
    allocate(            rho0_nph(nlevs_radial,0:nr_fine-1))
    allocate(              p0_nph(nlevs_radial,0:nr_fine-1))
    allocate(    p0_minus_peosbar(nlevs_radial,0:nr_fine-1))
    allocate(             peosbar(nlevs_radial,0:nr_fine-1))
    allocate(            w0_force(nlevs_radial,0:nr_fine-1))
    allocate(                Sbar(nlevs_radial,0:nr_fine-1))
    allocate(           beta0_nph(nlevs_radial,0:nr_fine-1))
    allocate(     gamma1bar_temp1(nlevs_radial,0:nr_fine-1))
    allocate(     gamma1bar_temp2(nlevs_radial,0:nr_fine-1))
    allocate(delta_gamma1_termbar(nlevs_radial,0:nr_fine-1))
    allocate(              w0_old(nlevs_radial,0:nr_fine))
    allocate(          beta0_edge(nlevs_radial,0:nr_fine))
    allocate( rho0_predicted_edge(nlevs_radial,0:nr_fine))
    allocate(        delta_chi_w0(nlevs_radial,0:nr_fine-1))

    if (verbose .ge. 1) then

       if (parallel_IOProcessor()) then
          do n = 1, nlevs
             print *, 'level: ', n
             print *, '   number of boxes = ', nboxes(pi(n)%la)
             print *, '   maximum zones   = ', (extent(mla%mba%pd(n),i),i=1,dm)
          end do
       end if

       do n=1,nlevs
          numcell = multifab_volume(pi(n),.false.)
          if (parallel_IOProcessor()) then
             print*,'Number of valid cells at level        ',n,numcell
          end if
          numcell = multifab_volume(pi(n),.true.)
          if (parallel_IOProcessor()) then
             print*,'Number of valid + ghost cells at level',n,numcell
          end if
       end do
       if (parallel_IOProcessor()) then
          print*,''
       end if
    end if

    if (barrier_timers) call parallel_barrier()
    misc_time = misc_time + parallel_wtime() - misc_time_start
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 1 -- react the full state and then base state through dt/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    react_time_start = parallel_wtime()

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  1 : react state     '
    end if

    do n=1,nlevs
       call multifab_build(s1(n),            mla%la(n), nscal, nghost(sold(n)))
    end do

    ! note: rho_omegadot2 and rho_Hnuc2 are just temporaries, and are overwritten in the
    ! next call to react_state before they are used
    call react_state(mla,tempbar_init,sold,s1,rho_omegadot2,rho_Hnuc2,rho_Hext,p0_old, &
                     half*dt,dx,the_bc_tower%bc_tower_array)

    do n=1,nlevs
       call destroy(rho_Hext(n))
    end do

    if (barrier_timers) call parallel_barrier()
    react_time = react_time + parallel_wtime() - react_time_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 2 -- define average expansion at time n+1/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    advect_time_start = parallel_wtime()
    
    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< CALLING advance_timestep with dt =',dt 
       write(6,*) '<<< STEP  2 : make w0 >>> '
    end if
    
    do n=1,nlevs
       call multifab_build(S_cc_nph(n), mla%la(n), 1, 0)
    end do

    if (time .eq. ZERO) then
        ! this is either a pressure iteration or the first time step
        ! set S_cc_nph = (1/2) (S_cc_old + S_cc_new)
       call make_S_at_halftime(mla,S_cc_nph,S_cc_old,S_cc_new, &
                               the_bc_tower%bc_tower_array)
    else
       ! set S_cc_nph = S_cc_old + (dt/2) * dSdt
       call extrap_to_halftime(mla,S_cc_nph,dSdt,S_cc_old,dt, &
                               the_bc_tower%bc_tower_array)
    end if

    do n=1,nlevs
       call multifab_build(delta_p_term(n), mla%la(n), 1, 0)
    end do

    ! compute p0_minus_peosbar = p0_old - peosbar (for making w0) and
    ! compute delta_p_term = peos_old - peosbar_cart (for RHS of projections)
    if (dpdt_factor .gt. ZERO ) then

       ! peos_old now holds the thermodynamic p computed from sold(rho,h,X)
       call makePfromRhoH(sold,sold,delta_p_term,mla,the_bc_tower%bc_tower_array)

       ! compute peosbar = Avg(peos_old)
       call average(mla,delta_p_term,peosbar,dx,1)

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
          call multifab_sub_sub(delta_p_term(n), peosbar_cart(n))
       end do
       
       do n=1,nlevs
          call destroy(peosbar_cart(n))
       end do

    else

       ! these should have no effect if dpdt_factor .le. 0
       p0_minus_peosbar = ZERO
       do n=1,nlevs
          call setval(delta_p_term(n),ZERO,all=.true.)
       end do

    end if

    if (spherical .eq. 1) then
       do n=1,nlevs
          call multifab_build(w0_force_cart(n),mla%la(n),dm,1)
          call setval(w0_force_cart(n),ZERO,all=.true.)
       end do
    end if

    if (dm .eq. 3) then
       do n=1,nlevs
          do comp=1,dm
             call multifab_build_edge(w0mac(n,comp),mla%la(n),1,1,comp)
             call setval(w0mac(n,comp),ZERO,all=.true.)
          end do
       end do
    end if

    if (evolve_base_state) then

       ! compute Sbar = average(S_cc_nph)
       call average(mla,S_cc_nph,Sbar,dx,1)

       ! save old-time value
       w0_old = w0

       ! compute w0, w0_force, and delta_chi_w0
       call make_w0(w0,w0_old,w0_force,Sbar,rho0_old,rho0_old,p0_old,p0_old,gamma1bar_old, &
                    gamma1bar_old,p0_minus_peosbar,psi,etarho_ec,etarho_cc,dt,dtold, &
                    delta_chi_w0,.true.)

       if (spherical .eq. 1) then
          ! put w0 on Cartesian edges
          call make_w0mac(mla,w0,w0mac,dx,the_bc_tower%bc_tower_array)
          ! put w0_force on Cartesian cells
          call put_1d_array_on_cart(w0_force,w0_force_cart,foextrap_comp,.false., &
                                    .true.,dx,the_bc_tower%bc_tower_array,mla)
       end if

    else

       ! these should have no effect if evolve_base_state = F
       w0_force = ZERO
       Sbar = ZERO

    end if
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 3 -- construct the advective velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  3 : create MAC velocities>>> '
    end if

    do n=1,nlevs
       do comp=1,dm
          call multifab_build_edge(umac(n,comp), mla%la(n),1,1,comp)
       end do
    end do
    
    ! compute unprojected MAC velocities
    call advance_premac(uold,sold,umac,gpi,normal,w0,w0mac,w0_force,w0_force_cart, &
                        rho0_old,grav_cell_old,dx,dt,the_bc_tower%bc_tower_array,mla)

    if (spherical .eq. 1) then
       do n=1,nlevs
          call destroy(w0_force_cart(n))
       end do
    end if

    do n=1,nlevs
       call multifab_build(delta_gamma1_term(n), mla%la(n), 1, 0)
       call multifab_build(macrhs(n),            mla%la(n), 1, 0)
       call multifab_build(delta_chi(n),         mla%la(n), 1, 0)

       call setval(delta_gamma1_term(n), ZERO, all=.true.)
       call setval(delta_chi(n),         ZERO, all=.true.)
    end do

    if (barrier_timers) call parallel_barrier()
    advect_time = advect_time + parallel_wtime() - advect_time_start

    macproj_time_start = parallel_wtime()

    ! compute RHS for MAC projection
    call make_macrhs(macrhs,rho0_old,S_cc_nph,delta_gamma1_term,Sbar,beta0_old,dx, &
                     gamma1bar_old,p0_old,delta_p_term,dt,delta_chi,.true.)

    do n=1,nlevs
       call destroy(delta_gamma1_term(n))
       call destroy(S_cc_nph(n))
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
             call multifab_build_edge(beta0_cart_edge(n,comp), mla%la(n),1,1,comp)
          end do
       end do

       call make_s0mac(mla,beta0_old,beta0_cart_edge,dx,foextrap_comp, &
                       the_bc_tower%bc_tower_array)

       call macproject(mla,umac,macphi,sold,dx,the_bc_tower,macrhs, &
                       beta0_cart_edge=beta0_cart_edge)

       do n=1,nlevs
          do comp=1,dm
             call destroy(beta0_cart_edge(n,comp))
          end do
       end do
    else
       call cell_to_edge(beta0_old,beta0_edge)
       call macproject(mla,umac,macphi,sold,dx,the_bc_tower, &
                       macrhs,beta0_1d=beta0_old,beta0_1d_edge=beta0_edge)
    end if

    do n=1,nlevs
       call destroy(macrhs(n))
    end do

    if (barrier_timers) call parallel_barrier()
    macproj_time = macproj_time + (parallel_wtime() - macproj_time_start)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 4 -- advect the base state and full state through dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    advect_time_start = parallel_wtime()

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  4 : advect base        '
    end if
    
    ! advect the base state density
    if (evolve_base_state) then
       call advect_base_dens(w0,rho0_old,rho0_new,rho0_predicted_edge,dt)
       call compute_cutoff_coords(rho0_new)
    else
       rho0_new = rho0_old
    end if

    do n=1,nlevs
       call multifab_build(thermal1(n), mla%la(n), 1, 0)
    end do
    
    ! thermal is the forcing for rhoh or temperature
    if(use_thermal_diffusion) then
       do n=1,nlevs
          call multifab_build( Tcoeff(n),  mla%la(n), 1,     1)
          call multifab_build(hcoeff1(n),  mla%la(n), 1,     1)
          call multifab_build(Xkcoeff1(n), mla%la(n), nspec, 1)
          call multifab_build(pcoeff1(n),  mla%la(n), 1,     1)
       end do

       call make_thermal_coeffs(s1,Tcoeff,hcoeff1,Xkcoeff1,pcoeff1)

       call make_explicit_thermal(mla,dx,thermal1,s1,Tcoeff,hcoeff1,Xkcoeff1,pcoeff1, &
                                  p0_old,the_bc_tower)

       do n=1,nlevs
          call destroy(Tcoeff(n))
       end do
    else
       do n=1,nlevs
          call setval(thermal1(n),ZERO,all=.true.)
       end do
    end if

    ! copy temperature into s2 for seeding eos calls only
    ! temperature will be overwritten later after enthalpy advance
    do n=1,nlevs
       call multifab_build(s2(n), mla%la(n), nscal, nghost(sold(n)))
       call multifab_copy_c(s2(n), temp_comp, s1(n), temp_comp, 1, nghost(sold(n)))
    end do

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '            :  density_advance >>> '
       write(6,*) '            :   tracer_advance >>> '
    end if

    do n=1,nlevs
       do comp = 1,dm
          call multifab_build_edge(sedge(n,comp),mla%la(n),nscal,0,comp)
          call multifab_build_edge(sflux(n,comp),mla%la(n),nscal,0,comp)
       end do
       if (ppm_trace_forces == 0) then
          call multifab_build(scal_force(n), mla%la(n), nscal, 1)
       else
          ! we need more ghostcells if we are tracing the forces
          call multifab_build(scal_force(n), mla%la(n), nscal, nghost(sold(n)))
       endif
       call multifab_build_edge(etarhoflux(n), mla%la(n), 1, 0, dm)
       ! set etarhoflux to zero
       call setval(etarhoflux(n),ZERO,all=.true.)
    end do

    ! advect rhoX, rho, and tracers
    call density_advance(mla,1,s1,s2,sedge,sflux,scal_force,umac,w0,w0mac,etarhoflux, &
                         rho0_old,rho0_new,p0_new,rho0_predicted_edge, &
                         dx,dt,the_bc_tower%bc_tower_array)

    ! compute the new etarho
    if (evolve_base_state .and. use_etarho) then
       if (spherical .eq. 0) then
          call make_etarho_planar(etarho_ec,etarho_cc,etarhoflux,mla)
       else
          call make_etarho_spherical(s1,s2,umac,w0mac,rho0_old,rho0_new,dx,normal, &
                                     etarho_ec,etarho_cc,mla,the_bc_tower%bc_tower_array)
       endif
    end if

    ! correct the base state density by "averaging"
    if (evolve_base_state .and. use_etarho) then
       call average(mla,s2,rho0_new,dx,rho_comp)
       call compute_cutoff_coords(rho0_new)
    end if

    if (evolve_base_state) then
       call make_grav_cell(grav_cell_new,rho0_new)
    else
       grav_cell_new = grav_cell_old
    end if

    ! base state pressure update
    if (evolve_base_state) then

       ! set new p0 through HSE
       p0_new = p0_old
       call enforce_HSE(rho0_new,p0_new,grav_cell_new)

       ! make psi
       if (spherical .eq. 0) then
          call make_psi_planar(etarho_cc,psi)
       else
          ! compute p0_nph
          p0_nph = HALF*(p0_old+p0_new)

          ! compute gamma1bar^{(1)} and store it in gamma1bar_temp1
          call make_gamma1bar(mla,s1,gamma1bar_temp1,p0_old,dx)

          ! compute gamma1bar^{(2),*} and store it in gamma1bar_temp2
          call make_gamma1bar(mla,s2,gamma1bar_temp2,p0_new,dx)

          ! compute gamma1bar^{nph,*} and store it in gamma1bar_temp2
          gamma1bar_temp2 = HALF*(gamma1bar_temp1+gamma1bar_temp2)

          ! make time-centered psi
          call make_psi_spherical(psi,w0,gamma1bar_temp2,p0_nph,Sbar)
       end if

    else

       p0_new = p0_old

    end if

    ! base state enthalpy update
    if (evolve_base_state) then

       ! compute rhoh0_old by "averaging"
       call average(mla,s1,rhoh0_old,dx,rhoh_comp)

       call advect_base_enthalpy(w0,rho0_old,rhoh0_old,rhoh0_new,rho0_predicted_edge,psi,dt)
    else
       rhoh0_new = rhoh0_old
    end if

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '            : enthalpy_advance >>> '
    end if

    call enthalpy_advance(mla,1,uold,s1,s2,sedge,sflux,scal_force,thermal1,umac,w0,w0mac, &
                          rho0_old,rhoh0_old,rho0_new,rhoh0_new,p0_old,p0_new, &
                          tempbar,psi,dx,dt,the_bc_tower%bc_tower_array)

    do n = 1, nlevs
       do comp = 1,dm
          call destroy(sedge(n,comp))
          call destroy(sflux(n,comp))
             call destroy(umac(n,comp))
       end do
       call destroy(scal_force(n))
    end do

    if (barrier_timers) call parallel_barrier()
    advect_time = advect_time + parallel_wtime() - advect_time_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 4a (Option I) -- Add thermal conduction (only enthalpy terms)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    thermal_time_start = parallel_wtime()

    if (use_thermal_diffusion) then
       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  4a: thermal conduct >>>'
       end if

       call thermal_conduct(mla,dx,dt,s1,hcoeff1,Xkcoeff1,pcoeff1,hcoeff1,Xkcoeff1,pcoeff1, &
                            s2,p0_old,p0_new,the_bc_tower)
    end if

    if (barrier_timers) call parallel_barrier()
    thermal_time = thermal_time + (parallel_wtime() - thermal_time_start)

    misc_time_start = parallel_wtime()

    ! pass temperature through for seeding the temperature update eos call
    ! pi goes along for the ride
    do n=1,nlevs
       call multifab_copy_c(s2(n),temp_comp,s1(n),temp_comp,1,nghost(sold(n)))
       call multifab_copy_c(s2(n),pi_comp,  s1(n),pi_comp,  1,nghost(sold(n)))
    end do

    ! now update temperature
    if (use_tfromp) then
       call makeTfromRhoP(s2,p0_new,mla,the_bc_tower%bc_tower_array,dx)
    else
       call makeTfromRhoH(s2,p0_new,mla,the_bc_tower%bc_tower_array,dx)
    end if

    if (use_thermal_diffusion) then
       ! make a copy of s2star since these are needed to compute
       ! coefficients in the call to thermal_conduct_full_alg
       do n=1,nlevs
          call multifab_build(s2star(n), mla%la(n), nscal, nghost(sold(n)))
          call multifab_copy_c(s2star(n), 1, s2(n), 1, nscal, nghost(sold(n)))
       end do

    end if

    if (barrier_timers) call parallel_barrier()
    misc_time = misc_time + parallel_wtime() - misc_time_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 5 -- react the full state and then base state through dt/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    react_time_start = parallel_wtime()

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  5 : react state     '
    end if

    do n=1,nlevs
       call multifab_build(rho_Hext(n), mla%la(n), 1, 0)
    end do
    
    call react_state(mla,tempbar_init,s2,snew,rho_omegadot2,rho_Hnuc2,rho_Hext,p0_new, &
                     half*dt,dx,the_bc_tower%bc_tower_array)

    do n=1,nlevs
       call destroy(s2(n))
    end do

    if (barrier_timers) call parallel_barrier()
    react_time = react_time + parallel_wtime() - react_time_start
    
    misc_time_start = parallel_wtime()

    if (evolve_base_state) then

       ! compute beta0 and gamma1bar
       call make_gamma1bar(mla,snew,gamma1bar_new,p0_new,dx)
       call make_beta0(beta0_new,rho0_new,p0_new,gamma1bar_new,grav_cell_new)

    else
        
       ! Just pass beta0 and gamma1bar through if not evolving base state
       beta0_new = beta0_old
       gamma1bar_new = gamma1bar_old

    end if

    beta0_nph = HALF*(beta0_old + beta0_new)

    if (barrier_timers) call parallel_barrier()
    misc_time = misc_time + parallel_wtime() - misc_time_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 6 -- define a new average expansion rate at n+1/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    advect_time_start = parallel_wtime()
    
    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  6 : make new S and new w0 >>> '
    end if

    ! reset cutoff coordinates to old time value
    if (evolve_base_state) then
       call compute_cutoff_coords(rho0_old)
    end if

    if(use_thermal_diffusion) then

       do n=1,nlevs
          call multifab_build(Tcoeff(n),  mla%la(n), 1,     1)
          call multifab_build(hcoeff2(n),  mla%la(n), 1,     1)
          call multifab_build(Xkcoeff2(n), mla%la(n), nspec, 1)
          call multifab_build(pcoeff2(n),  mla%la(n), 1,     1)
       end do

       call make_thermal_coeffs(snew,Tcoeff,hcoeff2,Xkcoeff2,pcoeff2)

       call make_explicit_thermal(mla,dx,thermal2,snew,Tcoeff,hcoeff2,Xkcoeff2,pcoeff2, &
                                  p0_new,the_bc_tower)

       do n=1,nlevs
          call destroy(hcoeff2(n))
          call destroy(Xkcoeff2(n))
          call destroy(pcoeff2(n))
       end do

    else
       do n=1,nlevs
          call setval(thermal2(n),ZERO,all=.true.)
       end do
    end if

    do n=1,nlevs
       call multifab_build(delta_gamma1_term(n), mla%la(n), 1, 0)
       call multifab_build(delta_gamma1(n), mla%la(n), 1, 0)
    end do

    ! p0 is only used for the delta_gamma1_term
    call make_S_cc(S_cc_new,delta_gamma1_term,delta_gamma1, &
                   snew,uold, &
                   normal, &
                   rho_omegadot2,rho_Hnuc2,rho_Hext,thermal2, &
                   p0_old,gamma1bar_new,delta_gamma1_termbar,psi, &
                   dx,mla,the_bc_tower%bc_tower_array)

    do n=1,nlevs
       call destroy(rho_Hext(n))
       call destroy(delta_gamma1(n))
    end do

    do n=1,nlevs
       call multifab_build(S_cc_nph(n), mla%la(n), 1, 0)
    end do

    call make_S_at_halftime(mla,S_cc_nph,S_cc_old,S_cc_new,the_bc_tower%bc_tower_array)

    do n=1,nlevs
       call multifab_build(delta_p_term(n), mla%la(n), 1, 0)
    end do

    ! compute p0_minus_peosbar = p0_new - peosbar (for making w0)
    ! and delta_p_term = peos_new - peosbar_cart (for RHS of projection)
    if (dpdt_factor .gt. ZERO) then

       ! peos_new now holds the thermodynamic p computed from snew(rho h X)
       call makePfromRhoH(snew,snew,delta_p_term,mla,the_bc_tower%bc_tower_array)

       ! compute peosbar = Avg(peos_new)
       call average(mla,delta_p_term,peosbar,dx,1)

       ! compute p0_minus_peosbar = p0_new - peosbar
       p0_minus_peosbar = p0_new - peosbar

       do n=1,nlevs
          call multifab_build(peosbar_cart(n), mla%la(n), 1, 0)
       end do

       ! compute peosbar_cart from peosbar
       call put_1d_array_on_cart(peosbar,peosbar_cart,foextrap_comp, &
                                 .false.,.false.,dx,the_bc_tower%bc_tower_array,mla)

       ! compute delta_p_term = peos_new - peosbar_cart
       do n=1,nlevs
          call multifab_sub_sub(delta_p_term(n), peosbar_cart(n))
       end do

       do n=1,nlevs
          call destroy(peosbar_cart(n))
       end do

    else

       ! this should have no effect if dpdt_factor .le. 0
       do n=1,nlevs
          call setval(delta_p_term(n),ZERO,all=.true.)
       end do

    end if

    if (spherical .eq. 1) then
       do n=1,nlevs
          call multifab_build(w0_force_cart(n), mla%la(n), dm, 1)
          call setval(w0_force_cart(n),ZERO,all=.true.)
       end do
    end if

    if (evolve_base_state) then

       call average(mla,S_cc_nph,Sbar,dx,1)

       if(use_delta_gamma1_term) then
          ! add delta_gamma1_termbar to Sbar
          Sbar = Sbar + delta_gamma1_termbar
       end if

       ! compute w0, w0_force, and delta_chi_w0
       call make_w0(w0,w0_old,w0_force,Sbar,rho0_old,rho0_new,p0_old,p0_new, &
                    gamma1bar_old,gamma1bar_new,p0_minus_peosbar, &
                    psi,etarho_ec,etarho_cc,dt,dtold,delta_chi_w0,.false.)

       if (spherical .eq. 1) then
          ! put w0 on Cartesian edges
          call make_w0mac(mla,w0,w0mac,dx,the_bc_tower%bc_tower_array)
          ! put w0_force on Cartesian cells
          call put_1d_array_on_cart(w0_force,w0_force_cart,foextrap_comp,.false., &
                                    .true.,dx,the_bc_tower%bc_tower_array,mla)
       end if
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 7 -- redo the construction of the advective velocity using the 
!! current w0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  7 : create MAC velocities >>> '
    end if

    do n=1,nlevs
       do comp=1,dm
          call multifab_build_edge(umac(n,comp),mla%la(n),1,1,comp)
       end do
    end do

    call advance_premac(uold,sold,umac,gpi,normal,w0,w0mac,w0_force,w0_force_cart, &
                        rho0_old,grav_cell_old,dx,dt,the_bc_tower%bc_tower_array,mla)

    if (barrier_timers) call parallel_barrier()
    advect_time = advect_time + parallel_wtime() - advect_time_start

    macproj_time_start = parallel_wtime()

    do n=1,nlevs
       call multifab_build(macrhs(n), mla%la(n), 1, 0)
    end do

    ! note delta_gamma1_term here is not time-centered
    call make_macrhs(macrhs,rho0_new,S_cc_nph,delta_gamma1_term,Sbar,beta0_nph,dx, &
                     gamma1bar_new,p0_new,delta_p_term,dt,delta_chi,.false.)

    do n=1,nlevs
       call destroy(delta_gamma1_term(n))
       call destroy(S_cc_nph(n))
       call destroy(delta_p_term(n))
       call destroy(delta_chi(n))
    end do

    do n=1,nlevs
       call multifab_build(rhohalf(n), mla%la(n), 1, 1)
    end do

    call make_at_halftime(rhohalf,sold,snew,rho_comp,1,the_bc_tower%bc_tower_array,mla)

    ! MAC projection !
    if (spherical .eq. 1) then
       do n=1,nlevs
          do comp=1,dm
             call multifab_build_edge(beta0_cart_edge(n,comp), mla%la(n),1,1,comp)
          end do
       end do

       call make_s0mac(mla,beta0_nph,beta0_cart_edge,dx,foextrap_comp, &
                       the_bc_tower%bc_tower_array)

       call macproject(mla,umac,macphi,rhohalf,dx,the_bc_tower,macrhs, &
                       beta0_cart_edge=beta0_cart_edge)

       do n=1,nlevs
          do comp=1,dm
             call destroy(beta0_cart_edge(n,comp))
          end do
       end do
    else
       call cell_to_edge(beta0_nph,beta0_edge)
       call macproject(mla,umac,macphi,rhohalf,dx,the_bc_tower,macrhs, &
                       beta0_1d=beta0_nph,beta0_1d_edge=beta0_edge)
    end if


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! advect the particles through dt using umac.  Then redistribute
    !! them
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 8 -- advect the base state and full state through dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    advect_time_start = parallel_wtime()

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  8 : advect base   '
    end if

    if (evolve_base_state) then
       call advect_base_dens(w0,rho0_old,rho0_new,rho0_predicted_edge,dt)
       call compute_cutoff_coords(rho0_new)
    else
       rho0_new = rho0_old
    end if

    do n=1,nlevs
       call multifab_build(s2(n), mla%la(n), nscal, nghost(sold(n)))
       ! copy temperature into s2 for seeding eos calls only
       ! temperature will be overwritten later after enthalpy advance
       call multifab_copy_c(s2(n), temp_comp, s1(n), temp_comp, 1, nghost(sold(n)))

       call setval(etarhoflux(n),ZERO,all=.true.)
    end do

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '            :  density_advance >>>'
       write(6,*) '            :   tracer_advance >>>'
    end if

    ! Build the sedge array.
    do n=1,nlevs
       do comp = 1,dm
          call multifab_build_edge(sedge(n,comp),mla%la(n),nscal,0,comp)
          call multifab_build_edge(sflux(n,comp),mla%la(n),nscal,0,comp)
       end do
       if (ppm_trace_forces == 0) then
          call multifab_build(scal_force(n), mla%la(n), nscal, 1)
       else
          call multifab_build(scal_force(n), mla%la(n), nscal, nghost(sold(n)))
       endif
    end do

    call density_advance(mla,2,s1,s2,sedge,sflux,scal_force,umac,w0,w0mac,etarhoflux, &
                         rho0_old,rho0_new,p0_new,rho0_predicted_edge,dx,dt, &
                         the_bc_tower%bc_tower_array)

    ! Now compute the new etarho
    if (evolve_base_state .and. use_etarho) then
       if (spherical .eq. 0) then
          call make_etarho_planar(etarho_ec,etarho_cc,etarhoflux,mla)
       else
          call make_etarho_spherical(s1,s2,umac,w0mac,rho0_old,rho0_new,dx,normal, &
                                     etarho_ec,etarho_cc,mla,the_bc_tower%bc_tower_array)
       endif
    end if

    do n=1,nlevs
       call destroy(etarhoflux(n))
    end do

    ! Correct the base state density using "averaging"
    if (evolve_base_state .and. use_etarho) then
       call average(mla,s2,rho0_new,dx,rho_comp)
       call compute_cutoff_coords(rho0_new)
    end if

    if (evolve_base_state) then
       call make_grav_cell(grav_cell_new,rho0_new)
       rho0_nph = HALF*(rho0_old+rho0_new)
       call make_grav_cell(grav_cell_nph,rho0_nph)
    else
       grav_cell_new = grav_cell_old
       rho0_nph = rho0_old
       grav_cell_nph = grav_cell_old
    end if

    ! base state pressure update
    if (evolve_base_state) then
       
       ! set new p0 through HSE
       p0_new = p0_old
       call enforce_HSE(rho0_new,p0_new,grav_cell_new)
       p0_nph = HALF*(p0_old+p0_new)

       ! make psi
       if (spherical .eq. 0) then
          call make_psi_planar(etarho_cc,psi)
       else
          p0_nph = HALF*(p0_old+p0_new)

          ! compute gamma1bar^{(2)} and store it in gamma1bar_temp2
          call make_gamma1bar(mla,s2,gamma1bar_temp2,p0_new,dx)

          ! compute gamma1bar^{nph} and store it in gamma1bar_temp2
          gamma1bar_temp2 = HALF*(gamma1bar_temp1+gamma1bar_temp2)

          call make_psi_spherical(psi,w0,gamma1bar_temp2,p0_nph,Sbar)
       end if

    else

       p0_new = p0_old

    end if

    ! base state enthalpy update
    if (evolve_base_state) then
       call advect_base_enthalpy(w0,rho0_old,rhoh0_old,rhoh0_new, &
                                 rho0_predicted_edge,psi,dt)
    else
       rhoh0_new = rhoh0_old
    end if

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '            : enthalpy_advance >>>'
    end if

    call enthalpy_advance(mla,2,uold,s1,s2,sedge,sflux,scal_force,thermal1,umac,w0,w0mac, &
                          rho0_old,rhoh0_old,rho0_new,rhoh0_new,p0_old,p0_new, &
                          tempbar,psi,dx,dt,the_bc_tower%bc_tower_array)

    do n=1,nlevs
       do comp = 1,dm
          call destroy(sedge(n,comp))
          call destroy(sflux(n,comp))
       end do
       call destroy(scal_force(n))
       call destroy(thermal1(n))
    end do

    if (barrier_timers) call parallel_barrier()
    advect_time = advect_time + parallel_wtime() - advect_time_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 8a (Option I) -- Add thermal conduction (only enthalpy terms)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    thermal_time_start = parallel_wtime()

    if (use_thermal_diffusion) then
       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  8a: thermal conduct >>>'
       end if

       do n=1,nlevs
          call multifab_build(Tcoeff(n),  mla%la(n), 1,     1)
          call multifab_build(hcoeff2(n),  mla%la(n), 1,     1)
          call multifab_build(Xkcoeff2(n), mla%la(n), nspec, 1)
          call multifab_build(pcoeff2(n),  mla%la(n), 1,     1)
       end do

       call make_thermal_coeffs(s2star,Tcoeff,hcoeff2,Xkcoeff2,pcoeff2)

       do n=1,nlevs
          call destroy(s2star(n))
          call destroy(Tcoeff(n))
       end do

       call thermal_conduct(mla,dx,dt,s1,hcoeff1,Xkcoeff1,pcoeff1,hcoeff2,Xkcoeff2,pcoeff2, &
                            s2,p0_old,p0_new,the_bc_tower)

       do n=1,nlevs
          call destroy(hcoeff1(n))
          call destroy(Xkcoeff1(n))
          call destroy(pcoeff1(n))
          call destroy(hcoeff2(n))
          call destroy(Xkcoeff2(n))
          call destroy(pcoeff2(n))
       end do
    end if

    if (barrier_timers) call parallel_barrier()
    thermal_time = thermal_time + (parallel_wtime() - thermal_time_start)

    misc_time_start = parallel_wtime()

    ! pass temperature through for seeding the temperature update eos call
    ! pi just goes along for the ride too
    do n=1,nlevs
       call multifab_copy_c(s2(n),temp_comp,s1(n),temp_comp,1,nghost(sold(n)))
       call multifab_copy_c(s2(n),pi_comp,  s1(n),pi_comp,  1,nghost(sold(n)))
    end do

    ! now update temperature
    if (use_tfromp) then
       call makeTfromRhoP(s2,p0_new,mla,the_bc_tower%bc_tower_array,dx)
    else
       call makeTfromRhoH(s2,p0_new,mla,the_bc_tower%bc_tower_array,dx)
    end if

    do n=1,nlevs
       call destroy(s1(n))
    end do

    if (barrier_timers) call parallel_barrier()
    misc_time = misc_time + parallel_wtime() - misc_time_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 9 -- react the full state and then base state through dt/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    react_time_start = parallel_wtime()

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  9 : react state '
    end if

    do n=1,nlevs
       call multifab_build(rho_Hext(n), mla%la(n), 1, 0)
    end do

    call react_state(mla,tempbar_init,s2,snew,rho_omegadot2,rho_Hnuc2,rho_Hext,p0_new, &
                     half*dt,dx,the_bc_tower%bc_tower_array)

    do n=1,nlevs
       call destroy(s2(n))
    end do

    if (barrier_timers) call parallel_barrier()
    react_time = react_time + parallel_wtime() - react_time_start

    misc_time_start = parallel_wtime()

    if (evolve_base_state) then

       ! compute beta0 and gamma1bar
       call make_gamma1bar(mla,snew,gamma1bar_new,p0_new,dx)
       call make_beta0(beta0_new,rho0_new,p0_new,gamma1bar_new,grav_cell_new)

    end if

    beta0_nph = HALF*(beta0_old+beta0_new)

    if (barrier_timers) call parallel_barrier()
    misc_time = misc_time + parallel_wtime() - misc_time_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 10 -- compute S^{n+1} for the final projection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ndproj_time_start = parallel_wtime()
    
    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP 10 : make new S >>>'
    end if
          
    if(use_thermal_diffusion) then
       do n=1,nlevs
          call multifab_build(Tcoeff(n),  mla%la(n), 1,     1)
          call multifab_build(hcoeff2(n),  mla%la(n), 1,     1)
          call multifab_build(Xkcoeff2(n), mla%la(n), nspec, 1)
          call multifab_build(pcoeff2(n),  mla%la(n), 1,     1)
       end do

       call make_thermal_coeffs(snew,Tcoeff,hcoeff2,Xkcoeff2,pcoeff2)

       call make_explicit_thermal(mla,dx,thermal2,snew,Tcoeff,hcoeff2,Xkcoeff2,pcoeff2, &
                                  p0_new,the_bc_tower)

       do n=1,nlevs
          call destroy(Tcoeff(n))
          call destroy(hcoeff2(n))
          call destroy(Xkcoeff2(n))
          call destroy(pcoeff2(n))
       end do
    end if
    
    do n=1,nlevs
       call multifab_build(delta_gamma1_term(n), mla%la(n), 1, 0)
       call multifab_build(delta_gamma1(n), mla%la(n), 1, 0)
    end do

    ! p0 is only used for the delta_gamma1_term
    call make_S_cc(S_cc_new,delta_gamma1_term,delta_gamma1, &
                   snew,uold, &
                   normal, &
                   rho_omegadot2,rho_Hnuc2,rho_Hext,thermal2, &
                   p0_new,gamma1bar_new,delta_gamma1_termbar,psi, &
                   dx,mla,the_bc_tower%bc_tower_array)

    do n=1,nlevs
       call destroy(delta_gamma1(n))
    end do

    if (evolve_base_state) then
       call average(mla,S_cc_new,Sbar,dx,1)

       if(use_delta_gamma1_term) then
          ! add delta_gamma1_termbar to Sbar
          Sbar = Sbar + delta_gamma1_termbar
       end if

    end if
    
    ! define dSdt = (S_cc_new - S_cc_old) / dt
    do n=1,nlevs
       call multifab_copy(dSdt(n),S_cc_new(n))
       call multifab_sub_sub(dSdt(n),S_cc_old(n))
       call multifab_div_div_s(dSdt(n),dt)
    end do

    if (barrier_timers) call parallel_barrier()
    ndproj_time = ndproj_time + parallel_wtime() - ndproj_time_start
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 11 -- update the velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    advect_time_start = parallel_wtime()

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP 11 : update and project new velocity >>>'
    end if
    
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
       end do
    end if

    if (spherical .eq. 1) then
       do n=1,nlevs
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
          call multifab_build(nodalrhs_old(n), mla%la(n), 1, 0, nodal)
          call multifab_copy(nodalrhs_old(n),nodalrhs(n))
       end do
       call make_nodalrhs(the_bc_tower,mla,nodalrhs,S_cc_new,delta_gamma1_term, &
                         Sbar,beta0_nph,dx)
       do n=1,nlevs
          call multifab_sub_sub(nodalrhs(n),nodalrhs_old(n))
          call multifab_div_div_s(nodalrhs(n),dt)
       end do

    else

       proj_type = regular_timestep_comp

       call make_nodalrhs(the_bc_tower,mla,nodalrhs,S_cc_new,delta_gamma1_term, &
                         Sbar,beta0_nph,dx)

       ! compute delta_p_term = peos_new - peosbar_cart (for RHS of projection)
       if (dpdt_factor .gt. ZERO) then

          do n=1,nlevs
             call multifab_build(delta_p_term(n), mla%la(n), 1, 0)
             call multifab_build(peosbar_cart(n), mla%la(n), 1, 0)
          end do

          ! peos_new now holds the thermodynamic p computed from snew(rho h X)
          call makePfromRhoH(snew,snew,delta_p_term,mla,the_bc_tower%bc_tower_array)

          ! compute peosbar = Avg(peos_new)
          call average(mla,delta_p_term,peosbar,dx,1)

          ! no need to compute p0_minus_peosbar since make_w0 is not called after here

          ! compute peosbar_cart from peosbar
          call put_1d_array_on_cart(peosbar,peosbar_cart,foextrap_comp, &
                                    .false.,.false.,dx,the_bc_tower%bc_tower_array,mla)

          ! compute delta_p_term = peos_new - peosbar_cart
          do n=1,nlevs
             call multifab_sub_sub(delta_p_term(n), peosbar_cart(n))
          end do
          
          call correct_nodalrhs(the_bc_tower,mla,rho0_new,nodalrhs,beta0_nph,dx,dt, &
                               gamma1bar_new,p0_new,delta_p_term)
          
          do n=1,nlevs
             call destroy(delta_p_term(n))
             call destroy(peosbar_cart(n))
          enddo
          
       end if

    end if

    do n=1,nlevs
       call destroy(delta_gamma1_term(n))
    end do

    do n=1,nlevs
       call multifab_build(beta0_cart(n), mla%la(n), 1, 1)
    end do
       
    call put_1d_array_on_cart(beta0_nph,beta0_cart,foextrap_comp,.false., &
                              .false.,dx,the_bc_tower%bc_tower_array,mla)

    call hgproject(proj_type,mla,unew,uold,rhohalf,pi,gpi,dx,dt,the_bc_tower, &
                   beta0_cart,nodalrhs)

    call make_pi_cc(mla,pi,snew,pi_comp,the_bc_tower%bc_tower_array,beta0_cart)

    do n=1,nlevs
       call destroy(beta0_cart(n))
       call destroy(rhohalf(n))
    end do
    
    ! If doing pressure iterations then put nodalrhs_old into nodalrhs to be returned to varden.
    if (init_mode) then
       do n=1,nlevs
          call multifab_copy(nodalrhs(n),nodalrhs_old(n))
          call destroy(nodalrhs_old(n))
       end do
    end if


    if (barrier_timers) call parallel_barrier()
    ndproj_time = ndproj_time + (parallel_wtime() - ndproj_time_start)

    misc_time_start = parallel_wtime()

    if (.not. init_mode) then
       
       if (.not. fix_base_state) then
          ! compute tempbar by "averaging"
          call average(mla,snew,tempbar,dx,temp_comp)
       end if

       ! output any runtime diagnostics
       ! pass in the new time value, time+dt
       call diag(time+dt,dt,dx,snew,rho_Hnuc2,rho_Hext,thermal2,rho_omegadot2,&
                 rho0_new,rhoh0_new,p0_new,tempbar, &
                 gamma1bar_new,beta0_new, &
                 unew,w0,normal, &
                 mla,the_bc_tower)


       ! perform sanity checks, if desired
       if (mach_max_abort > ZERO) then
          call sanity_check(time+dt,dx,snew, &
                 rho0_new,rhoh0_new,p0_new,tempbar, &
                 gamma1bar_new,beta0_new, &
                 unew,w0,normal, &
                 mla,the_bc_tower)
       endif

    end if

    if (barrier_timers) call parallel_barrier()
    misc_time = misc_time + parallel_wtime() - misc_time_start

    call destroy(bpt)

    nreduce = 5
    if (use_thermal_diffusion) nreduce = nreduce + 1
    allocate (times_local(nreduce))
    allocate (times_global(nreduce))

    ! pack
    times_local(1) = advect_time
    times_local(2) = macproj_time
    times_local(3) = ndproj_time
    times_local(4) = react_time
    times_local(5) = misc_time
    if (use_thermal_diffusion) times_local(6) = thermal_time
    
    call parallel_reduce(times_global, times_local, MPI_MAX, &
                         proc=parallel_IOProcessorNode())
 
    if (parallel_IOProcessor()) then
       print *, 'Timing summary:'
       print *, '   Advection       : ', times_global(1), ' seconds'
       print *, '   MAC   Projection: ', times_global(2), ' seconds'
       print *, '   Nodal Projection: ', times_global(3), ' seconds'
       if (use_thermal_diffusion) &
          print *, '   Thermal         : ', times_global(6), ' seconds'
       print *, '   Reactions       : ', times_global(4), ' seconds'
       print *, '   Misc            : ', times_global(5), ' seconds'
       print *, ' '
    endif

    deallocate(times_local, times_global)

  end subroutine advance_timestep

end module advance_timestep_module
