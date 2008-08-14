module advance_timestep_module

  implicit none

  private

  public :: advance_timestep

contains
    
  subroutine advance_timestep(init_mode,mla,uold,sold,unew,snew, &
                              gpres,pres,normal,rho0_old,rhoh0_old, &
                              rho0_new,rhoh0_new,p0_old,p0_new,tempbar,gamma1bar,w0, &
                              rho_omegadot2,div_coeff_old,div_coeff_new, &
                              grav_cell_old,dx,time,dt,dtold,the_bc_tower, &
                              dSdt,Source_old,Source_new,etarho,etarho_cc,psi,sponge,hgrhs, &
                              istep)

    use bl_prof_module
    use ml_layout_module
    use bl_constants_module
    use multifab_module
    use pre_advance_module
    use velocity_advance_module
    use scalar_advance_module
    use density_advance_module
    use macrhs_module
    use macproject_module
    use hgrhs_module
    use hgproject_module
    use proj_parameters
    use bc_module
    use box_util_module
    use make_div_coeff_module
    use make_w0_module
    use advect_base_module
    use correct_base_module
    use react_base_module
    use react_state_module
    use make_S_module
    use average_module
    use phihalf_module
    use extraphalf_module
    use thermal_conduct_module
    use make_explicit_thermal_module
    use add_react_to_thermal_module
    use variables, only: nscal, press_comp, temp_comp, rho_comp, rhoh_comp, foextrap_comp
    use geometry, only: spherical, nr_fine, r_end_coord, anelastic_cutoff_coord, &
         base_cutoff_density_coord, burning_cutoff_density_coord
    use network, only: nspec
    use make_grav_module
    use make_eta_module
    use make_psi_module
    use fill_3d_module
    use cell_to_edge_module
    use define_bc_module
    use make_gamma_module
    use rhoh_vs_t_module
    use probin_module
    use bl_prof_module
    use ml_restriction_module
    use multifab_physbc_module
    use multifab_fill_ghost_module

    
    logical,         intent(in   ) :: init_mode
    type(ml_layout), intent(inout) :: mla
    type(multifab),  intent(in   ) :: uold(:)
    type(multifab),  intent(in   ) :: sold(:)
    type(multifab),  intent(inout) :: unew(:)
    type(multifab),  intent(inout) :: snew(:)
    type(multifab),  intent(inout) :: gpres(:)
    type(multifab),  intent(inout) :: pres(:)
    type(multifab),  intent(in   ) :: normal(:)
    real(dp_t)    ,  intent(inout) :: rho0_old(:,0:)
    real(dp_t)    ,  intent(inout) :: rhoh0_old(:,0:)
    real(dp_t)    ,  intent(inout) :: rho0_new(:,0:)
    real(dp_t)    ,  intent(inout) :: rhoh0_new(:,0:)
    real(dp_t)    ,  intent(inout) :: p0_old(:,0:)
    real(dp_t)    ,  intent(inout) :: p0_new(:,0:)
    real(dp_t)    ,  intent(inout) :: tempbar(:,0:)
    real(dp_t)    ,  intent(inout) :: gamma1bar(:,0:)
    real(dp_t)    ,  intent(inout) :: w0(:,0:)
    type(multifab),  intent(inout) :: rho_omegadot2(:)
    real(dp_t)    ,  intent(inout) :: div_coeff_old(:,0:)
    real(dp_t)    ,  intent(inout) :: div_coeff_new(:,0:)
    real(dp_t)    ,  intent(in   ) :: grav_cell_old(:,0:)
    real(dp_t)    ,  intent(in   ) :: dx(:,:),time,dt,dtold
    type(bc_tower),  intent(in   ) :: the_bc_tower
    type(multifab),  intent(inout) :: dSdt(:)
    type(multifab),  intent(inout) :: Source_old(:)
    type(multifab),  intent(inout) :: Source_new(:)
    real(dp_t)    ,  intent(inout) :: etarho(:,0:)
    real(dp_t)    ,  intent(inout) :: etarho_cc(:,0:)
    real(dp_t)    ,  intent(inout) :: psi(:,0:)
    type(multifab),  intent(in   ) :: sponge(:)
    type(multifab),  intent(inout) :: hgrhs(:)
    integer       ,  intent(in   ) :: istep

    ! local
    type(multifab) :: rhohalf(mla%nlevel)
    type(multifab) :: w0mac(mla%nlevel,mla%dim)
    type(multifab) :: w0_force_cart_vec(mla%nlevel)
    type(multifab) :: macrhs(mla%nlevel)
    type(multifab) :: macphi(mla%nlevel)
    type(multifab) :: hgrhs_old(mla%nlevel)
    type(multifab) :: Source_nph(mla%nlevel)
    type(multifab) :: thermal(mla%nlevel)
    type(multifab) :: s2star(mla%nlevel)
    type(multifab) :: rho_omegadot2_hold(mla%nlevel)
    type(multifab) :: s1(mla%nlevel)
    type(multifab) :: s2(mla%nlevel)
    type(multifab) :: delta_gamma1_term(mla%nlevel)
    type(multifab) :: delta_gamma1(mla%nlevel)
    type(multifab) :: rho_omegadot1(mla%nlevel)
    type(multifab) :: rho_Hext(mla%nlevel)
    type(multifab) :: div_coeff_3d(mla%nlevel)
    type(multifab) :: gamma1(mla%nlevel)
    type(multifab) :: umac(mla%nlevel,mla%dim)
    type(multifab) :: utrans(mla%nlevel,mla%dim)
    type(multifab) :: etarhoflux(mla%nlevel)
    type(multifab) :: ptherm_old(mla%nlevel)
    type(multifab) :: ptherm_nph(mla%nlevel)
    type(multifab) :: ptherm_new(mla%nlevel)
    type(multifab) :: p0_cart(mla%nlevel)
    type(multifab) :: delta_p_term(mla%nlevel)

    real(dp_t), allocatable :: grav_cell_nph(:,:)
    real(dp_t), allocatable :: grav_cell_new(:,:)
    real(dp_t), allocatable :: rho0_nph(:,:)
    real(dp_t), allocatable :: p0_nph(:,:)
    real(dp_t), allocatable :: p0_minus_pthermbar(:,:)
    real(dp_t), allocatable :: w0_force(:,:)
    real(dp_t), allocatable :: w0_old(:,:)
    real(dp_t), allocatable :: Sbar(:,:)
    real(dp_t), allocatable :: div_coeff_nph(:,:)
    real(dp_t), allocatable :: div_coeff_edge(:,:)
    real(dp_t), allocatable :: rho_omegadotbar(:,:,:)
    real(dp_t), allocatable :: rho_Hextbar(:,:)
    real(dp_t), allocatable :: rhoh0_1(:,:)
    real(dp_t), allocatable :: rhoh0_2(:,:)
    real(dp_t), allocatable :: rho0_predicted_edge(:,:)
    real(dp_t), allocatable :: gamma1bar_old(:,:)
    real(dp_t), allocatable :: delta_gamma1_termbar(:,:)

    integer    :: r,n,dm,comp,nlevs,ng_s,proj_type
    real(dp_t) :: halfdt,eps_in
    logical    :: nodal(mla%dim),umac_nodal_flag(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "advance_timestep")

    dm = mla%dim
    nlevs = mla%nlevel

    allocate(       grav_cell_nph(nlevs,0:nr_fine-1))
    allocate(       grav_cell_new(nlevs,0:nr_fine-1))
    allocate(            rho0_nph(nlevs,0:nr_fine-1))
    allocate(              p0_nph(nlevs,0:nr_fine-1))
    allocate(  p0_minus_pthermbar(nlevs,0:nr_fine-1))
    allocate(            w0_force(nlevs,0:nr_fine-1))
    allocate(              w0_old(nlevs,0:nr_fine  ))
    allocate(                Sbar(nlevs,0:nr_fine-1))
    allocate(       div_coeff_nph(nlevs,0:nr_fine-1))
    allocate(      div_coeff_edge(nlevs,0:nr_fine  ))
    allocate(     rho_omegadotbar(nlevs,0:nr_fine-1,nspec))
    allocate(         rho_Hextbar(nlevs,0:nr_fine-1))
    allocate(             rhoh0_1(nlevs,0:nr_fine-1))
    allocate(             rhoh0_2(nlevs,0:nr_fine-1))
    allocate( rho0_predicted_edge(nlevs,0:nr_fine  ))
    allocate(       gamma1bar_old(nlevs,0:nr_fine-1))
    allocate(delta_gamma1_termbar(nlevs,0:nr_fine-1))

    ! Set this to zero so if evolve_base_state = F there is no effect in update_vel
    w0_force = ZERO

    ! Set Sbar to zero so if evolve_base_state = F then it doesn't affect rhs of projections
    Sbar = ZERO

    ! Set these to results from last
    w0_old = w0
    gamma1bar_old = gamma1bar

    nodal = .true.
    ng_s = sold(1)%ng
    halfdt = half*dt

    ! compute the coordinates of the anelastic cutoff
    anelastic_cutoff_coord(1) = r_end_coord(1,1)+1
    do r=0,r_end_coord(1,1)
       if (rho0_old(1,r) .lt. anelastic_cutoff .and. &
            anelastic_cutoff_coord(1) .eq. r_end_coord(1,1)+1) then
          anelastic_cutoff_coord(1) = r
          exit
       end if
    end do
    do n=2,nlevs
       anelastic_cutoff_coord(n) = 2*anelastic_cutoff_coord(n-1)
    end do

    ! compute the coordinates of the base cutoff density
    base_cutoff_density_coord(1) = r_end_coord(1,1)+1
    do r=0,r_end_coord(1,1)
       if (rho0_old(1,r) .lt. base_cutoff_density .and. &
            base_cutoff_density_coord(1) .eq. r_end_coord(1,1)+1) then
          base_cutoff_density_coord(1) = r
          exit
       end if
    end do
    do n=2,nlevs
       base_cutoff_density_coord(n) = 2*base_cutoff_density_coord(n-1)
    end do

    ! compute the coordinates of the burning cutoff density
    burning_cutoff_density_coord(1) = r_end_coord(1,1)+1
    do r=0,r_end_coord(1,1)
       if (rho0_old(1,r) .lt. burning_cutoff_density .and. &
            burning_cutoff_density_coord(1) .eq. r_end_coord(1,1)+1) then
          burning_cutoff_density_coord(1) = r
          exit
       end if
    end do
    do n=2,nlevs
       burning_cutoff_density_coord(n) = 2*burning_cutoff_density_coord(n-1)
    end do
    
    ! tempbar is only used as an initial guess for eos calls
    call average(mla,sold,tempbar,dx,temp_comp)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 1 -- define average expansion at time n+1/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< CALLING advance_timestep with dt =',dt 
       write(6,*) '<<< STEP  1 : make w0 >>> '
    end if
    
    do n=1,nlevs
       call multifab_build(Source_nph(n), mla%la(n), 1, 1)
    end do

    if (init_mode) then
       call make_S_at_halftime(nlevs,mla,Source_nph,Source_old,Source_new, &
                               the_bc_tower%bc_tower_array)
    else
       call extrap_to_halftime(nlevs,mla,Source_nph,dSdt,Source_old,dt, &
                               the_bc_tower%bc_tower_array)
    end if

    if (spherical .eq. 1) then
       do n=1,nlevs
          do comp=1,dm
             umac_nodal_flag = .false.
             umac_nodal_flag(comp) = .true.
             call multifab_build(w0mac(n,comp), mla%la(n),  1, 1, nodal = umac_nodal_flag)
             call setval(w0mac(n,comp), ZERO, all=.true.)
          end do
       end do
    end if

    do n=1,nlevs
       call multifab_build(ptherm_old(n), mla%la(n), 1, 0)
       call multifab_build(p0_cart(n), mla%la(n), 1, 1)
    end do
    
    ! ptherm_old now holds the thermodynamic p computed from sold(rho,h,X)
    call makePfromRhoH(nlevs,sold,ptherm_old,tempbar,mla,the_bc_tower%bc_tower_array,dx)

    ! p0_cart now holds p0_old
    call put_1d_array_on_cart(nlevs,p0_old,p0_cart,foextrap_comp,.false.,.false.,dx, &
                              the_bc_tower%bc_tower_array,mla)

    ! p0_cart now holds (p0_old - ptherm_old)
    do n=1,nlevs
       call multifab_sub_sub(p0_cart(n),ptherm_old(n))
    enddo

    ! fill the ghostcells of p0_cart
    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(p0_cart(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(p0_cart(nlevs),1,foextrap_comp,1, &
                            the_bc_tower%bc_tower_array(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(p0_cart(n-1),p0_cart(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(p0_cart(n),p0_cart(n-1), &
                                         p0_cart(n)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n), &
                                         1,foextrap_comp,1)
       enddo

    end if
    
    ! p0_minus_pthermbar now holds 
    ! Avg(p0_cart) = Avg(p0_old - ptherm_old) = p0_old - pthermbar_old
    call average(mla,p0_cart,p0_minus_pthermbar,dx,1)

    do n=1,nlevs
       call multifab_build(delta_p_term(n), mla%la(n), 1, 1)
    enddo

    ! delta_p_term now holds p0_minus_pthermbar = p0_old - pthermbar_old
    call put_1d_array_on_cart(nlevs,p0_minus_pthermbar,delta_p_term,foextrap_comp, &
                              .false.,.false.,dx,the_bc_tower%bc_tower_array,mla)

    ! fill the ghostcells of delta_p_term
    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(delta_p_term(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(delta_p_term(nlevs),1,foextrap_comp,1, &
                            the_bc_tower%bc_tower_array(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(delta_p_term(n-1),delta_p_term(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(delta_p_term(n),delta_p_term(n-1), &
                                         delta_p_term(n)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n), &
                                         1,foextrap_comp,1)
       enddo

    end if

    ! finish computing delta_p_term = (p0_old - pthermbar_old) - (p0_old - ptherm_old)
    ! = ptherm_old - pthermbar_old
    do n = 1,nlevs
       call multifab_sub_sub(delta_p_term(n),p0_cart(n),1)
    enddo

    do n=1,nlevs
       call destroy(p0_cart(n))
    enddo

    if (dm .eq. 3) then
       do n=1,nlevs
          call multifab_build(w0_force_cart_vec(n), mla%la(n), dm, 1)
          call setval(w0_force_cart_vec(n),ZERO,all=.true.)
       end do
    end if

    if (evolve_base_state) then

       call average(mla,Source_nph,Sbar,dx,1)

       call make_w0(nlevs,w0,w0_old,w0_force,Sbar,rho0_old,rho0_old,p0_old,p0_old, &
                    gamma1bar,gamma1bar,p0_minus_pthermbar,psi,etarho,etarho_cc,dt,dtold)

       if (spherical .eq. 1) then
          call put_w0_on_edges(mla,w0,w0mac,dx,div_coeff_old,the_bc_tower)
       end if

       if (dm .eq. 3) then
          call put_1d_array_on_cart(nlevs,w0_force,w0_force_cart_vec,foextrap_comp, &
                                    .false.,.true.,dx,the_bc_tower%bc_tower_array,mla, &
                                    normal)
       end if

    end if
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 2 -- construct the advective velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  2 : create MAC velocities>>> '
    end if

    do n=1,nlevs
       do comp=1,dm
          umac_nodal_flag = .false.
          umac_nodal_flag(comp) = .true.
          call multifab_build(  umac(n,comp), mla%la(n),  1, 1, nodal = umac_nodal_flag)
          call multifab_build(utrans(n,comp), mla%la(n),  1, 1, nodal = umac_nodal_flag)
       end do
    end do
    
    call advance_premac(nlevs,uold,sold,umac,utrans,gpres,normal,w0,w0mac, &
                        w0_force,w0_force_cart_vec,rho0_old,grav_cell_old,dx,dt, &
                        the_bc_tower%bc_tower_array,mla)

    if (dm .eq. 3) then
       do n=1,nlevs
          call destroy(w0_force_cart_vec(n))
       end do
    end if

    do n=1,nlevs
       call multifab_build(delta_gamma1_term(n), mla%la(n), 1, 0)
       call multifab_build(macrhs(n),            mla%la(n), 1, 0)
       call setval(delta_gamma1_term(n), ZERO, all=.true.)
    end do

    call make_macrhs(nlevs,macrhs,rho0_old,Source_nph,delta_gamma1_term,Sbar, &
                     div_coeff_old,dx, &
                     gamma1bar,gamma1bar,p0_old,p0_old,delta_p_term,dt)

    do n=1,nlevs
       call destroy(delta_gamma1_term(n))
       call destroy(Source_nph(n))
       call destroy(delta_p_term(n))
    end do

    do n=1,nlevs
       call multifab_build(macphi(n), mla%la(n), 1, 1)
       call setval(macphi(n), ZERO, all=.true.)
    end do

    ! MAC projection !
    if (spherical .eq. 1) then
       do n=1,nlevs
          call multifab_build(div_coeff_3d(n), mla%la(nlevs), 1, 1)
       end do

       call put_1d_array_on_cart(nlevs,div_coeff_old,div_coeff_3d,foextrap_comp,.false., &
                                 .false.,dx,the_bc_tower%bc_tower_array,mla)

       call macproject(mla,umac,macphi,sold,dx,the_bc_tower, &
                       press_comp, macrhs,div_coeff_3d=div_coeff_3d)

       do n=1,nlevs
          call destroy(div_coeff_3d(n))
       end do
    else
       do n=1,nlevs
          call cell_to_edge(n,div_coeff_old(n,:),div_coeff_edge(n,:))
       end do
       call macproject(mla,umac,macphi,sold,dx,the_bc_tower, &
                       press_comp, macrhs,div_coeff_1d=div_coeff_old, &
                       div_coeff_half_1d=div_coeff_edge)
    end if

    if(do_half_alg) then
       do n=1,nlevs
          call destroy(macphi(n))
       end do
    end if

    do n=1,nlevs
       call destroy(macrhs(n))
    end do
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 3 -- react the full state and then base state through dt/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  3 : react state     '
    end if

    do n=1,nlevs
       call multifab_build(s1(n),            mla%la(n), nscal, ng_s)
       call multifab_build(rho_omegadot1(n), mla%la(n), nspec, 1)
       call multifab_build(rho_Hext(n),      mla%la(n), 1,     1)
    end do

    call react_state(nlevs,mla,sold,s1,rho_omegadot1,rho_Hext,halfdt,dx, &
                     the_bc_tower%bc_tower_array,time)
    
    if (full_rhoh0_evolution) then
       if (evolve_base_state) then
          if (parallel_IOProcessor() .and. verbose .ge. 1) then
             write(6,*) '            : react  base >>> '
          end if
          do comp=1,nspec
             call average(mla,rho_omegadot1,rho_omegadotbar(:,:,comp),dx,comp)
          end do
          call average(mla,rho_Hext,rho_Hextbar,dx,1)
          call react_base(nlevs,rhoh0_old,rho_omegadotbar,rho_Hextbar,halfdt,rhoh0_1)
       else
          rhoh0_1 = rhoh0_old
       end if
    else
       call average(mla,s1,rhoh0_1,dx,rhoh_comp)
    end if

    if (evolve_base_state) then
       do n=1,nlevs
          call multifab_build(gamma1(n), mla%la(n), 1, 1)
       end do
       
       call make_gamma(nlevs,mla,gamma1,s1,p0_old,tempbar,dx, &
                       the_bc_tower%bc_tower_array)
       call average(mla,gamma1,gamma1bar,dx,1)

       do n=1,nlevs
          call destroy(gamma1(n))
       end do
    end if

    do n=1,nlevs
       call destroy(rho_Hext(n))
    end do

    do n=1,nlevs
       call make_grav_cell(n,grav_cell_new(n,:),rho0_old(n,:))
    end do

    call make_div_coeff(nlevs,div_coeff_new,rho0_old,p0_old,gamma1bar,grav_cell_new)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 4 -- advect the base state and full state through dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  4 : advect base        '
    end if
    
    if (evolve_base_state) then
       call advect_base(nlevs,w0,Sbar,p0_old,p0_new,rho0_old,rho0_new,rhoh0_1,rhoh0_2, &
                        gamma1bar,div_coeff_new,rho0_predicted_edge,psi,dx(:,dm),dt)
    else
       rho0_new = rho0_old
       rhoh0_2 = rhoh0_1
       p0_new = p0_old
    end if

    do n=1,nlevs
       call multifab_build(thermal(n), mla%la(n), 1, 1)
    end do
    
    ! thermal is the forcing for rhoh or temperature
    if(use_thermal_diffusion) then
       call make_explicit_thermal(mla,dx,thermal,s1,p0_old, &
                                  the_bc_tower,temp_diffusion_formulation)
    else
       do n=1,nlevs
          call setval(thermal(n),ZERO,all=.true.)
       end do
    end if
    
    ! if we are predicting temperature at edges, we need to add the reaction 
    ! terms to thermal
    if ( (enthalpy_pred_type .eq. predict_T_then_rhohprime) .or. &
         (enthalpy_pred_type .eq. predict_T_then_h        ) ) then

       do n=1,nlevs
          call multifab_build(rho_Hext(n), mla%la(n), 1, 1)
       end do

       if(istep .le. 1) then
          call add_react_to_thermal(nlevs,thermal,rho_omegadot1,s1,rho_Hext, &
                                    the_bc_tower%bc_tower_array,mla,dx,time)
       else
          call add_react_to_thermal(nlevs,thermal,rho_omegadot2,s1,rho_Hext, &
                                    the_bc_tower%bc_tower_array,mla,dx,time)
          
          if(.not. do_half_alg) then
             do n=1,nlevs
                call multifab_build(rho_omegadot2_hold(n), mla%la(n), nspec, 1)
                call multifab_copy_c(rho_omegadot2_hold(n),1,rho_omegadot2(n),1,nspec,1)
             end do
          end if
       end if

       do n=1,nlevs
          call destroy(rho_Hext(n))
       end do

    end if
            
    if(do_half_alg) then
       do n=1,nlevs
          call destroy(rho_omegadot1(n))
       end do
    end if

    do n=1,nlevs
       call multifab_build(s2(n), mla%la(n), nscal, ng_s)
    end do

    ! Build etarhoflux here so that we can call correct_base before make_etarho.
    umac_nodal_flag = .false.
    umac_nodal_flag(dm) = .true.
    do n=1,nlevs
       call multifab_build(etarhoflux(n), mla%la(n), 1, nodal = umac_nodal_flag)
       call setval(etarhoflux(n),ZERO,all=.true.)
    end do

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '            : density_advance >>> '
       write(6,*) '            : scalar_advance >>> '
    end if

    call density_advance(nlevs,mla,1,s1,s2,&
                         umac,w0,w0mac,etarhoflux,normal, &
                         rho0_old,rho0_new,&
                         p0_new,rho0_predicted_edge, &
                         dx,dt,the_bc_tower%bc_tower_array)
    call scalar_advance(nlevs,mla,1,uold,s1,s2,thermal, &
                        umac,w0,w0mac,etarhoflux,utrans,normal, &
                        rho0_old,rhoh0_1, &
                        rho0_new,rhoh0_2, &
                        p0_old,p0_new,tempbar,psi,rho0_predicted_edge, &
                        dx,dt,the_bc_tower%bc_tower_array)

    ! Correct the base state using the lagged etarho and psi
    if (use_etarho .and. evolve_base_state) then
       call correct_base(nlevs,rho0_new,etarho,dx(:,dm),dt)
    end if

    ! Now compute the new etarho and psi
    if (evolve_base_state) then
       if (use_etarho) then

          if (spherical .eq. 0) then
             call make_etarho_planar(nlevs,etarho,etarho_cc,etarhoflux,mla)
          else
             call make_etarho_spherical(nlevs,s1,s2,umac,rho0_old,rho0_new,dx,normal, &
                                        etarho,etarho_cc,mla,the_bc_tower%bc_tower_array)
          endif

       endif

       call make_psi(nlevs,etarho,psi,w0,gamma1bar,p0_old,p0_new,Sbar)
    end if


    do n=1,nlevs
       call destroy(thermal(n))
    end do

    if(.not. do_half_alg) then
       do n=1,nlevs
          do comp=1,dm
             call destroy(umac(n,comp))
             call destroy(utrans(n,comp))
          end do
       end do
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 4a (Option I) -- Add thermal conduction (only enthalpy terms)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if (use_thermal_diffusion) then
       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  4a: thermal conduct >>>'
       end if
       
       if(do_half_alg) then
          call thermal_conduct_half_alg(mla,dx,dt,s1,s2,p0_old,p0_new, &
                                        tempbar,the_bc_tower)
       else
          call thermal_conduct_full_alg(mla,dx,dt,s1,s1,s2,p0_old,p0_new, &
                                        tempbar,the_bc_tower)
          
          ! make a copy of s2star since these are needed to compute
          ! coefficients in the call to thermal_conduct_full_alg
          do n=1,nlevs
             call multifab_build(s2star(n), mla%la(n), nscal, ng_s)
             call multifab_copy_c(s2star(n), 1, s2(n), 1, nscal, ng_s)
          end do
       end if
    end if

    if(do_half_alg) then
       do n=1,nlevs
          call destroy(s1(n))
       end do
    end if
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 5 -- react the full state and then base state through dt/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP  5 : react state     '
    end if

    do n=1,nlevs
       call multifab_build(rho_Hext(n), mla%la(n), 1, 1)
    end do
    
    call react_state(nlevs,mla,s2,snew,rho_omegadot2,rho_Hext,halfdt,dx, &
                     the_bc_tower%bc_tower_array,time)

    do n=1,nlevs
       call destroy(s2(n))
    end do

    if (evolve_base_state .and. full_rhoh0_evolution) then
       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '            : react  base >>> '
       end if
       do comp=1,nspec
          call average(mla,rho_omegadot2,rho_omegadotbar(:,:,comp),dx,comp)
       end do
       call average(mla,rho_Hext,rho_Hextbar,dx,1)
       call react_base(nlevs,rhoh0_2,rho_omegadotbar,rho_Hextbar,halfdt,rhoh0_new)
    else
       rhoh0_new = rhoh0_2
    end if

    if (evolve_base_state) then
       do n=1,nlevs
          call multifab_build(gamma1(n), mla%la(n), 1, 1)
       end do
       
       call make_gamma(nlevs,mla,gamma1,snew,p0_new,tempbar,dx, &
                       the_bc_tower%bc_tower_array)
       call average(mla,gamma1,gamma1bar,dx,1)

       do n=1,nlevs
          call destroy(gamma1(n))
       end do
    end if

    do n=1,nlevs
       call make_grav_cell(n,grav_cell_new(n,:),rho0_new(n,:))
    end do

    call make_div_coeff(nlevs,div_coeff_new,rho0_new,p0_new,gamma1bar,grav_cell_new)
    
    div_coeff_nph = HALF*(div_coeff_old + div_coeff_new)

    rho0_nph = HALF*(rho0_old+rho0_new)

    ! Define base state at half time for use in velocity advance
    do n=1,nlevs
       call make_grav_cell(n,grav_cell_nph(n,:),rho0_nph(n,:))
    end do
    
    if(.not. do_half_alg) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 6 -- define a new average expansion rate at n+1/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  6 : make new S and new w0 >>> '
       end if
       
       do n=1,nlevs
          call multifab_build(thermal(n), mla%la(n), 1, 1)
       end do

       if(use_thermal_diffusion) then
          call make_explicit_thermal(mla,dx,thermal,snew,p0_new, &
                                     the_bc_tower,temp_diffusion_formulation)
       else
          do n=1,nlevs
             call setval(thermal(n),ZERO,all=.true.)
          end do
       end if
       
       do n=1,nlevs
          call multifab_build(delta_gamma1_term(n), mla%la(n), 1, 0)
          call multifab_build(delta_gamma1(n), mla%la(n), 1, 0)
       end do

       ! p0 is only used for the delta_gamma1_term
       call make_S(nlevs,Source_new,delta_gamma1_term,delta_gamma1, &
                   snew,uold,rho_omegadot2,rho_Hext,thermal, &
                   p0_old,gamma1bar,delta_gamma1_termbar,psi,dx, &
                   mla,the_bc_tower%bc_tower_array)
       
       do n=1,nlevs
          call destroy(rho_Hext(n))
          call destroy(thermal(n))
          call destroy(delta_gamma1(n))
       end do

       do n=1,nlevs
          call multifab_build(Source_nph(n), mla%la(n), 1, 1)
       end do

       call make_S_at_halftime(nlevs,mla,Source_nph,Source_old,Source_new, &
                               the_bc_tower%bc_tower_array)

       do n=1,nlevs
          call multifab_build(ptherm_new(n), mla%la(n), 1, 0)
          call multifab_build(ptherm_nph(n), mla%la(n), 1, 0)
          call multifab_build(p0_cart(n), mla%la(n), 1, 1)
       end do
       
       ! ptherm_new now holds the thermodynamic p computed from snew(rho h X)
       call makePfromRhoH(nlevs,snew,ptherm_new,tempbar,mla,the_bc_tower%bc_tower_array,dx)

       p0_nph = HALF*(p0_old + p0_new)

       ! p0_cart now holds p0_nph
       call put_1d_array_on_cart(nlevs,p0_nph,p0_cart,foextrap_comp,.false.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla)

       do n=1,nlevs
          call multifab_copy(ptherm_nph(n), ptherm_old(n))
          call multifab_plus_plus(ptherm_nph(n), ptherm_new(n))
          call multifab_div_div_s(ptherm_nph(n), TWO)
       enddo

       ! p0_cart now holds (p0_nph - ptherm_nph)
       do n=1, nlevs
          call multifab_sub_sub(p0_cart(n),ptherm_nph(n))
       enddo

       do n=1,nlevs
          call destroy(ptherm_new(n))
          call destroy(ptherm_nph(n))
       end do

       ! fill the ghostcells of p0_cart
       if (nlevs .eq. 1) then

          ! fill ghost cells for two adjacent grids at the same level
          ! this includes periodic domain boundary ghost cells
          call multifab_fill_boundary(p0_cart(nlevs))

          ! fill non-periodic domain boundary ghost cells
          call multifab_physbc(p0_cart(nlevs),1,foextrap_comp,1, &
                               the_bc_tower%bc_tower_array(nlevs))
       
       else

          ! the loop over nlevs must count backwards to make sure the finer grids are 
          ! done first
          do n=nlevs,2,-1

             ! set level n-1 data to be the average of the level n data covering it
             call ml_cc_restriction(p0_cart(n-1)    ,p0_cart(n)    ,mla%mba%rr(n-1,:))

             ! fill level n ghost cells using interpolation from level n-1 data
             ! note that multifab_fill_boundary and multifab_physbc are called for
             ! both levels n-1 and n
             call multifab_fill_ghost_cells(p0_cart(n),p0_cart(n-1), &
                                            p0_cart(n)%ng,mla%mba%rr(n-1,:), &
                                            the_bc_tower%bc_tower_array(n-1), &
                                            the_bc_tower%bc_tower_array(n), &
                                            1,foextrap_comp,1)
          enddo

       end if

       ! p0_minus_pthermbar now holds
       ! Avg(p0_cart) = Avg(p0_nph - ptherm_nph) = p0_nph - pthermbar_nph
       call average(mla,p0_cart,p0_minus_pthermbar,dx,1)

       do n=1,nlevs
          call multifab_build(delta_p_term(n), mla%la(n), 1, 1)
       enddo

       ! delta_p_term now holds p0_minus_pthermbar = p0_nph - pthermbar_nph
       call put_1d_array_on_cart(nlevs,p0_minus_pthermbar,delta_p_term,foextrap_comp, &
                                 .false.,.false.,dx,the_bc_tower%bc_tower_array,mla)

       ! now fill the ghostcells of delta_p_term
       if (nlevs .eq. 1) then

          ! fill ghost cells for two adjacent grids at the same level
          ! this includes periodic domain boundary ghost cells
          call multifab_fill_boundary(delta_p_term(nlevs))

          ! fill non-periodic domain boundary ghost cells
          call multifab_physbc(delta_p_term(nlevs),1,foextrap_comp,1, &
                               the_bc_tower%bc_tower_array(nlevs))
       
       else

          ! the loop over nlevs must count backwards to make sure the finer grids are 
          ! done first
          do n=nlevs,2,-1

             ! set level n-1 data to be the average of the level n data covering it
             call ml_cc_restriction(delta_p_term(n-1),delta_p_term(n),mla%mba%rr(n-1,:))

             ! fill level n ghost cells using interpolation from level n-1 data
             ! note that multifab_fill_boundary and multifab_physbc are called for
             ! both levels n-1 and n
             call multifab_fill_ghost_cells(delta_p_term(n),delta_p_term(n-1), &
                                            delta_p_term(n)%ng,mla%mba%rr(n-1,:), &
                                            the_bc_tower%bc_tower_array(n-1), &
                                            the_bc_tower%bc_tower_array(n), &
                                            1,foextrap_comp,1)
          enddo

       end if

       ! finish computing delta_p_term = (p0_nph - pthermbar_nph) - (p0_nph - ptherm_nph)
       ! = ptherm_nph - pthermbar_nph
       do n = 1,nlevs
          call multifab_sub_sub(delta_p_term(n),p0_cart(n),1)
       enddo

       do n=1,nlevs
          call multifab_destroy(p0_cart(n))
       enddo

       if (dm .eq. 3) then
          do n=1,nlevs
             call multifab_build(w0_force_cart_vec(n), mla%la(n), dm, 1)
             call setval(w0_force_cart_vec(n),ZERO,all=.true.)
          end do
       end if

       if (evolve_base_state) then
       
          call average(mla,Source_nph,Sbar,dx,1)

          if(use_delta_gamma1_term) then
             ! add delta_gamma1_termbar to Sbar
             Sbar = Sbar + delta_gamma1_termbar
          end if

          call make_w0(nlevs,w0,w0_old,w0_force,Sbar,rho0_old,rho0_new,p0_old,p0_new, &
                       gamma1bar_old,gamma1bar,p0_minus_pthermbar,psi,etarho,etarho_cc,dt, &
                       dtold)

          if (spherical .eq. 1) then
             call put_w0_on_edges(mla,w0,w0mac,dx,div_coeff_nph,the_bc_tower)
          end if

          if (dm .eq. 3) then
             call put_1d_array_on_cart(nlevs,w0_force,w0_force_cart_vec,foextrap_comp, &
                                       .false.,.true.,dx,the_bc_tower%bc_tower_array,mla, &
                                       normal)
          end if
       end if
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 7 -- redo the construction of the advective velocity using the current w0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  7 : create MAC velocities >>> '
       end if

       do n=1,nlevs
          do comp=1,dm
             umac_nodal_flag = .false.
             umac_nodal_flag(comp) = .true.
             call multifab_build(  umac(n,comp), mla%la(n),  1, 1, nodal = umac_nodal_flag)
             call multifab_build(utrans(n,comp), mla%la(n),  1, 1, nodal = umac_nodal_flag)
          end do
       end do

       call advance_premac(nlevs,uold,sold,umac,utrans,gpres,normal,w0,w0mac, &
                           w0_force,w0_force_cart_vec,rho0_old,grav_cell_old,dx,dt, &
                           the_bc_tower%bc_tower_array,mla)

       do n=1,nlevs
          call multifab_build(macrhs(n), mla%la(n), 1, 0)
       end do

       ! note delta_gamma1_term here is not time-centered
       call make_macrhs(nlevs,macrhs,rho0_old,Source_nph,delta_gamma1_term,Sbar, &
                        div_coeff_nph,dx, &
                        gamma1bar_old,gamma1bar,p0_old,p0_new,delta_p_term,dt)
    
       do n=1,nlevs
          call destroy(delta_gamma1_term(n))
          call destroy(Source_nph(n))
          call destroy(delta_p_term(n))
       end do

       do n=1,nlevs
          call multifab_build(rhohalf(n), mla%la(n), 1, 1)
       end do

       call make_at_halftime(nlevs,rhohalf,sold,snew,rho_comp,1, &
                             the_bc_tower%bc_tower_array,mla)
       
       ! MAC projection !
       if (spherical .eq. 1) then
          do n=1,nlevs
             call multifab_build(div_coeff_3d(n), mla%la(nlevs), 1, 1)
          end do

          call put_1d_array_on_cart(nlevs,div_coeff_nph,div_coeff_3d,foextrap_comp,.false., &
                                    .false.,dx,the_bc_tower%bc_tower_array,mla)

          call macproject(mla,umac,macphi,rhohalf,dx,the_bc_tower, &
                          press_comp,macrhs,div_coeff_3d=div_coeff_3d)

          do n=1,nlevs
             call destroy(div_coeff_3d(n))
          end do
       else
          do n=1,nlevs
             call cell_to_edge(n,div_coeff_nph(n,:),div_coeff_edge(n,:))
          end do
          call macproject(mla,umac,macphi,rhohalf,dx,the_bc_tower, &
                          press_comp,macrhs,div_coeff_1d=div_coeff_nph, &
                          div_coeff_half_1d=div_coeff_edge)
       end if

       do n=1,nlevs
          call destroy(rhohalf(n))
          call destroy(macrhs(n))
          call destroy(macphi(n))
       end do
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 8 -- advect the base state and full state through dt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  8 : advect base   '
       end if

       if (evolve_base_state) then
          call advect_base(nlevs,w0,Sbar,p0_old,p0_new,rho0_old,rho0_new,rhoh0_1,rhoh0_2, &
                           gamma1bar,div_coeff_nph,rho0_predicted_edge,psi,dx(:,dm),dt)
       else
          rho0_new = rho0_old
          rhoh0_2 = rhoh0_1
          p0_new = p0_old
       end if
              
       do n=1,nlevs
          call multifab_build(thermal(n), mla%la(n), 1, 1)
       end do

       ! thermal is the forcing for rhoh or temperature
       if(use_thermal_diffusion) then
          call make_explicit_thermal(mla,dx,thermal,s1,p0_old, &
                                     the_bc_tower,temp_diffusion_formulation)
       else
          do n=1,nlevs
             call setval(thermal(n),ZERO,all=.true.)
          end do
       end if
       
       ! if we are predicting temperature at edges, we need to add the reaction 
       ! terms to thermal
       if ( (enthalpy_pred_type .eq. predict_T_then_rhohprime) .or. &
            (enthalpy_pred_type .eq. predict_T_then_h        ) ) then

          do n=1,nlevs
             call multifab_build(rho_Hext(n), mla%la(n), 1, 1)
          end do

          if(istep .le. 1) then
             call add_react_to_thermal(nlevs,thermal,rho_omegadot1,s1,rho_Hext, &
                                       the_bc_tower%bc_tower_array,mla,dx,time)
          else
             call add_react_to_thermal(nlevs,thermal,rho_omegadot2_hold,s1,rho_Hext, &
                                       the_bc_tower%bc_tower_array,mla,dx,time)

             do n=1,nlevs
                call destroy(rho_omegadot2_hold(n))
             end do
          end if

          do n=1,nlevs
             call destroy(rho_Hext(n))
          end do


       end if

       do n=1,nlevs
          call destroy(rho_omegadot1(n))
       end do
       
       do n=1,nlevs
          call multifab_build(s2(n), mla%la(n), nscal, ng_s)
          call setval(etarhoflux(n),ZERO,all=.true.)
       end do

       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '            : density_advance >>>'
          write(6,*) '            : scalar_advance >>>'
       end if

       call density_advance(nlevs,mla,2,s1,s2,&
                            umac,w0,w0mac,etarhoflux,normal, &
                            rho0_old,rho0_new,&
                            p0_new,rho0_predicted_edge, &
                            dx,dt,the_bc_tower%bc_tower_array)
       call scalar_advance(nlevs,mla,2,uold,s1,s2,thermal, &
                           umac,w0,w0mac,etarhoflux,utrans,normal, &
                           rho0_old,rhoh0_1, &
                           rho0_new,rhoh0_2, &
                           p0_old,p0_new,tempbar,psi,rho0_predicted_edge, &
                           dx,dt,the_bc_tower%bc_tower_array)

       ! Correct the base state using the lagged etarho and psi
       if (use_etarho .and. evolve_base_state) then
          call correct_base(nlevs,rho0_new,etarho,dx(:,dm),dt)
       end if

       ! Now compute the new etarho and psi
       if (evolve_base_state) then
          if (use_etarho) then

             if (spherical .eq. 0) then
                call make_etarho_planar(nlevs,etarho,etarho_cc,etarhoflux,mla)
             else
                call make_etarho_spherical(nlevs,s1,s2,umac,rho0_old,rho0_new,dx,normal, &
                                           etarho,etarho_cc,mla,the_bc_tower%bc_tower_array)
             endif

          endif

          call make_psi(nlevs,etarho,psi,w0,gamma1bar,p0_old,p0_new,Sbar)
       end if

       do n=1,nlevs
          call destroy(thermal(n))
          call destroy(etarhoflux(n))
       end do
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 8a (Option I) -- Add thermal conduction (only enthalpy terms)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       if (use_thermal_diffusion) then
          if (parallel_IOProcessor() .and. verbose .ge. 1) then
             write(6,*) '<<< STEP  8a: thermal conduct >>>'
          end if
          
          call thermal_conduct_full_alg(mla,dx,dt,s1,s2star,s2,p0_old,p0_new, &
                                        tempbar,the_bc_tower)

          do n=1,nlevs
             call destroy(s2star(n))
          end do
       end if

       do n=1,nlevs
          call destroy(s1(n))
       end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 9 -- react the full state and then base state through dt/2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  9 : react state '
       end if

       do n=1,nlevs
          call multifab_build(rho_Hext(n), mla%la(n), 1, 1)
       end do
       
       call react_state(nlevs,mla,s2,snew,rho_omegadot2,rho_Hext,halfdt,dx,&
                        the_bc_tower%bc_tower_array,time)

       do n=1,nlevs
          call destroy(s2(n))
       end do

       if (evolve_base_state .and. full_rhoh0_evolution) then
          if (parallel_IOProcessor() .and. verbose .ge. 1) then
             write(6,*) '            : react  base >>>'
          end if
          do comp=1,nspec
             call average(mla,rho_omegadot2,rho_omegadotbar(:,:,comp),dx,comp)
          end do
          call average(mla,rho_Hext,rho_Hextbar,dx,1)
          call react_base(nlevs,rhoh0_2,rho_omegadotbar,rho_Hextbar,halfdt,rhoh0_new)
       else
          rhoh0_new = rhoh0_2
       end if

       if (evolve_base_state) then
          do n=1,nlevs
             call multifab_build(gamma1(n), mla%la(n), 1, 1)
          end do
          
          call make_gamma(nlevs,mla,gamma1,snew,p0_new,tempbar,dx, &
                          the_bc_tower%bc_tower_array)
          call average(mla,gamma1,gamma1bar,dx,1)
          
          do n=1,nlevs
             call destroy(gamma1(n))
          end do
       end if
       
       do n=1,nlevs
          call make_grav_cell(n,grav_cell_new(n,:),rho0_new(n,:))
       end do

       call make_div_coeff(nlevs,div_coeff_new,rho0_new,p0_new,gamma1bar,grav_cell_new)
              
    end if ! end if corresponding to .not. do_half_alg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 10 -- compute S^{n+1} for the final projection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP 10 : make new S >>>'
    end if
          
    do n=1,nlevs
       call multifab_destroy(ptherm_old(n))
       call multifab_build(thermal(n), mla%la(n), 1, 1)
    end do

    if(use_thermal_diffusion) then
       call make_explicit_thermal(mla,dx,thermal,snew,p0_new, &
                                  the_bc_tower,temp_diffusion_formulation)
    else
       do n=1,nlevs
          call setval(thermal(n),ZERO,all=.true.)
       end do
    end if
    
    do n=1,nlevs
       call multifab_build(delta_gamma1_term(n), mla%la(n), 1, 0)
       call multifab_build(delta_gamma1(n), mla%la(n), 1, 0)
    end do

    ! p0 is only used for the delta_gamma1_term
    call make_S(nlevs,Source_new,delta_gamma1_term,delta_gamma1, &
                snew,uold,rho_omegadot2,rho_Hext,thermal, &
                p0_new,gamma1bar,delta_gamma1_termbar,psi,dx, &
                mla,the_bc_tower%bc_tower_array)

    do n=1,nlevs
       call destroy(rho_Hext(n))
       call destroy(thermal(n))
       call destroy(delta_gamma1(n))
    end do

    if (evolve_base_state) then
       call average(mla,Source_new,Sbar,dx,1)

       if(use_delta_gamma1_term) then
          ! add delta_gamma1_termbar to Sbar
          Sbar = Sbar + delta_gamma1_termbar
       end if

    end if
    
    ! define dSdt = (Source_new - Source_old) / dt
    do n=1,nlevs
       call multifab_copy(dSdt(n),Source_new(n))
       call multifab_sub_sub(dSdt(n),Source_old(n))
       call multifab_div_div_s(dSdt(n),dt)
    end do
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! STEP 11 -- update the velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if (parallel_IOProcessor() .and. verbose .ge. 1) then
       write(6,*) '<<< STEP 11 : update and project new velocity >>>'
    end if
    
    ! Define rho at half time using the new rho from Step 8!
    do n=1,nlevs
       call multifab_build(rhohalf(n), mla%la(n), 1, 1)
    end do

    call make_at_halftime(nlevs,rhohalf,sold,snew,rho_comp,1,the_bc_tower%bc_tower_array,mla)
    
    call velocity_advance(nlevs,mla,uold,unew,sold,rhohalf,umac,utrans,gpres, &
                          normal,w0,w0mac,w0_force,w0_force_cart_vec, &
                          rho0_old,rho0_nph, &
                          grav_cell_old,grav_cell_nph,dx,dt, &
                          the_bc_tower%bc_tower_array,sponge)

    if (init_mode) then
       ! throw away w0 by setting w0 = w0_old
       w0 = w0_old
    end if

    do n=1,nlevs
       do comp=1,dm
          call destroy(umac(n,comp))
          call destroy(utrans(n,comp))
       end do
    end do

    if (spherical .eq. 1) then
       do n=1,nlevs
          do comp=1,dm
             call destroy(w0mac(n,comp))
          end do
       end do
    end if

    if (dm .eq. 3) then
       do n=1,nlevs
          call destroy(w0_force_cart_vec(n))
       end do
    end if

    ! Define beta at half time using the div_coeff_new from step 9!
    div_coeff_nph = HALF*(div_coeff_old+div_coeff_new)
       
    ! Project the new velocity field.
    if (init_mode) then

       proj_type = pressure_iters_comp

       do n=1,nlevs
          call multifab_build(hgrhs_old(n), mla%la(n), 1, 0, nodal)
          call multifab_copy(hgrhs_old(n),hgrhs(n))
       end do
       call make_hgrhs(nlevs,the_bc_tower,mla,hgrhs,Source_new,delta_gamma1_term, &
                       Sbar,div_coeff_new,dx)
       do n=1,nlevs
          call multifab_sub_sub(hgrhs(n),hgrhs_old(n))
          call multifab_div_div_s(hgrhs(n),dt)
       end do

    else

       proj_type = regular_timestep_comp
       call make_hgrhs(nlevs,the_bc_tower,mla,hgrhs,Source_new,delta_gamma1_term, &
                       Sbar,div_coeff_new,dx)

       if (dpdt_factor .gt. ZERO) then

          do n=1,nlevs
             call multifab_build(ptherm_new(n), mla%la(n), 1, 0)
             call multifab_build(p0_cart(n), mla%la(n), 1, 1)
          enddo

          ! ptherm_new now holds the thermodynamic p computed from snew(rho h X)
          call makePfromRhoH(nlevs,snew,ptherm_new,tempbar,mla, &
                             the_bc_tower%bc_tower_array,dx)
          
          ! p0_cart now holds p0_new
          call put_1d_array_on_cart(nlevs,p0_new,p0_cart,foextrap_comp,.false.,.false.,dx, &
                                    the_bc_tower%bc_tower_array,mla)
          
          ! p0_cart now holds (p0_new - ptherm_new)
          do n=1,nlevs
             call multifab_sub_sub(p0_cart(n),ptherm_new(n))
          enddo

          do n=1,nlevs
             call destroy(ptherm_new(n))
          end do
          
          ! fill ghostcells for p0_cart
          if (nlevs .eq. 1) then

             ! fill ghost cells for two adjacent grids at the same level
             ! this includes periodic domain boundary ghost cells
             call multifab_fill_boundary(p0_cart(nlevs))
       
             ! fill non-periodic domain boundary ghost cells
             call multifab_physbc(p0_cart(nlevs),1,foextrap_comp,1, &
                                  the_bc_tower%bc_tower_array(nlevs))

          else

             ! the loop over nlevs must count backwards to make sure the finer grids are 
             ! done first
             do n=nlevs,2,-1

                ! set level n-1 data to be the average of the level n data covering it
                call ml_cc_restriction(p0_cart(n-1)    ,p0_cart(n)    ,mla%mba%rr(n-1,:))

                ! fill level n ghost cells using interpolation from level n-1 data
                ! note that multifab_fill_boundary and multifab_physbc are called for
                ! both levels n-1 and n
                call multifab_fill_ghost_cells(p0_cart(n),p0_cart(n-1), &
                                               p0_cart(n)%ng,mla%mba%rr(n-1,:), &
                                               the_bc_tower%bc_tower_array(n-1), &
                                               the_bc_tower%bc_tower_array(n), &
                                               1,foextrap_comp,1)
             enddo
             
          end if

          ! p0_minus_pthermbar now holds
          ! Avg(p0_cart) = Avg(p0_new - ptherm_new) = p0_new - pthermbar_new
          call average(mla,p0_cart,p0_minus_pthermbar,dx,1)

          do n=1,nlevs
             call multifab_build(delta_p_term(n), mla%la(n), 1, 1)
          enddo

          ! delta_p_term now holds p0_minus_pthermbar = p0_new - pthermbar_new
          call put_1d_array_on_cart(nlevs,p0_minus_pthermbar,delta_p_term,foextrap_comp, &
                                    .false.,.false.,dx,the_bc_tower%bc_tower_array,mla)
          
          ! fill the ghostcells of delta_p_term
          if (nlevs .eq. 1) then

             ! fill ghost cells for two adjacent grids at the same level
             ! this includes periodic domain boundary ghost cells
             call multifab_fill_boundary(delta_p_term(nlevs))

             ! fill non-periodic domain boundary ghost cells
             call multifab_physbc(delta_p_term(nlevs),1,foextrap_comp,1, &
                                  the_bc_tower%bc_tower_array(nlevs))

          else

             ! the loop over nlevs must count backwards to make sure the finer grids are 
             ! done first
             do n=nlevs,2,-1

                ! set level n-1 data to be the average of the level n data covering it
                call ml_cc_restriction(delta_p_term(n-1),delta_p_term(n),mla%mba%rr(n-1,:))

                ! fill level n ghost cells using interpolation from level n-1 data
                ! note that multifab_fill_boundary and multifab_physbc are called for
                ! both levels n-1 and n
                call multifab_fill_ghost_cells(delta_p_term(n),delta_p_term(n-1), &
                                               delta_p_term(n)%ng,mla%mba%rr(n-1,:), &
                                               the_bc_tower%bc_tower_array(n-1), &
                                               the_bc_tower%bc_tower_array(n), &
                                               1,foextrap_comp,1)
             enddo

          end if

          ! finish computing delta_p_term = (p0_new - pthermbar_new) - (p0_new - ptherm_new)
          ! = ptherm_new - pthermbar_new
          do n = 1,nlevs
             call multifab_sub_sub(delta_p_term(n),p0_cart(n),1)
          enddo
          
          do n=1,nlevs
             call destroy(p0_cart(n))
          enddo
          
          call correct_hgrhs(nlevs,the_bc_tower,mla,rho0_new,hgrhs,div_coeff_new,dx,dt, &
                             gamma1bar,p0_new,delta_p_term)
          
          do n=1,nlevs
             call destroy(delta_p_term(n))
          enddo
          
       end if

    end if

    do n=1,nlevs
       call destroy(delta_gamma1_term(n))
    end do

    if (spherical .eq. 1) then
       do n=1,nlevs
          call multifab_build(div_coeff_3d(n), mla%la(nlevs), 1, 1)
       end do
       
       call put_1d_array_on_cart(nlevs,div_coeff_nph,div_coeff_3d,foextrap_comp,.false., &
                                 .false.,dx,the_bc_tower%bc_tower_array,mla)

       eps_in = 1.d-12
       call hgproject(proj_type, mla, unew, uold, rhohalf, pres, gpres, dx, dt, &
                      the_bc_tower, press_comp, &
                      hgrhs, div_coeff_3d=div_coeff_3d, eps_in = eps_in)

       do n=1,nlevs
          call destroy(div_coeff_3d(n))
       end do
    else
       call hgproject(proj_type, mla, unew, uold, rhohalf, pres, gpres, dx, dt, &
                      the_bc_tower, press_comp, &
                      hgrhs, div_coeff_1d=div_coeff_nph)
    end if

    do n=1,nlevs
       call destroy(rhohalf(n))
    end do
    
    ! If doing pressure iterations then put hgrhs_old into hgrhs to be returned to varden.
    if (init_mode) then
       do n=1,nlevs
          call multifab_copy(hgrhs(n),hgrhs_old(n))
          call destroy(hgrhs_old(n))
       end do
    end if

    call destroy(bpt)
    
  end subroutine advance_timestep

end module advance_timestep_module
