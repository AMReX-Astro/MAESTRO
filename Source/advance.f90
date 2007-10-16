module advance_timestep_module

  use BoxLib
  use omp_module
  use f2kcli
  use list_box_module
  use ml_boxarray_module
  use layout_module
  use multifab_module
  use pre_advance_module
  use velocity_advance_module
  use scalar_advance_module
  use macrhs_module
  use macproject_module
  use hgrhs_module
  use hgproject_module
  use proj_parameters
  use bc_module
  use bl_mem_stat_module
  use bl_timer_module
  use box_util_module
  use bl_IO_module
  use fabio_module
  use make_div_coeff_module
  use make_w0_module
  use advect_base_module
  use react_base_module
  use react_state_module
  use make_S_module
  use average_module
  use phihalf_module
  use extraphalf_module
  use thermal_conduct_module
  use make_explicit_thermal_module
  use add_react_to_thermal_module
  use variables
  use network
  use probin_module

  contains

    subroutine advance_timestep(init_mode, &
                                mla,uold,sold,s1,s2,unew,snew,umac,uedge,sedge,utrans,gp,p, &
                                scal_force,normal, &
                                s0_old,s0_1,s0_2,s0_new,p0_old,p0_1,p0_2,p0_new,gam1,w0, &
                                rho_omegadot1, rho_omegadot2, rho_Hext, &
                                div_coeff_old,div_coeff_new,&
                                grav_cell_old, &
                                dx,time,dt,dtold,the_bc_tower, &
                                anelastic_cutoff,verbose,mg_verbose,cg_verbose,&
                                dSdt,Source_old,Source_new,gamma1_term, &
                                sponge,do_sponge,hgrhs,istep)

    implicit none

    logical, intent(in) :: init_mode

    type(ml_layout),intent(inout) :: mla
    type(multifab), intent(inout) :: uold(:)
    type(multifab), intent(inout) :: sold(:)
    type(multifab), intent(inout) ::   s1(:)
    type(multifab), intent(inout) ::   s2(:)
    type(multifab), intent(inout) :: unew(:)
    type(multifab), intent(inout) :: snew(:)
    type(multifab), intent(inout) :: umac(:,:)
    type(multifab), intent(inout) :: uedge(:,:)
    type(multifab), intent(inout) :: sedge(:,:)
    type(multifab), intent(inout) :: utrans(:,:)
    type(multifab), intent(inout) :: gp(:)
    type(multifab), intent(inout) :: p(:)
    type(multifab), intent(inout) :: scal_force(:)
    type(multifab), intent(in   ) :: normal(:)
    type(multifab), intent(inout) :: rho_omegadot1(:)
    type(multifab), intent(inout) :: rho_omegadot2(:)
    type(multifab), intent(inout) :: rho_Hext(:)
    type(multifab), intent(inout) :: dSdt(:)
    type(multifab), intent(inout) :: Source_old(:)
    type(multifab), intent(inout) :: Source_new(:)
    type(multifab), intent(inout) :: gamma1_term(:)
    type(multifab), intent(inout) :: hgrhs(:)
    real(dp_t)    , intent(inout) :: s0_old(0:,:)
    real(dp_t)    , intent(inout) :: s0_1(0:,:)
    real(dp_t)    , intent(inout) :: s0_2(0:,:)
    real(dp_t)    , intent(inout) :: s0_new(0:,:)
    real(dp_t)    , intent(inout) :: p0_old(0:)
    real(dp_t)    , intent(inout) :: p0_1(0:)
    real(dp_t)    , intent(inout) :: p0_2(0:)
    real(dp_t)    , intent(inout) :: p0_new(0:)
    real(dp_t)    , intent(inout) :: gam1(0:)
    real(dp_t)    , intent(inout) :: w0(0:)
    real(dp_t)    , intent(in   ) :: div_coeff_old(0:)
    real(dp_t)    , intent(inout) :: div_coeff_new(0:)
    real(dp_t)    , intent(in   ) :: grav_cell_old(0:)
    real(dp_t)    , intent(in   ) :: dx(:,:), time, dt, dtold
    type(bc_tower), intent(in   ) :: the_bc_tower
    real(dp_t)    , intent(in   ) :: anelastic_cutoff
    integer       , intent(in   ) :: verbose,mg_verbose,cg_verbose,istep

    type(multifab), intent(in   ) :: sponge(:)
    logical       , intent(in   ) :: do_sponge

    type(multifab), allocatable :: rhohalf(:)
    type(multifab), allocatable :: w0_cart_vec(:)
    type(multifab), allocatable :: w0_force_cart_vec(:)
    type(multifab), allocatable :: macrhs(:)
    type(multifab), allocatable :: macphi(:)
    type(multifab), allocatable ::  hgrhs_old(:)
    type(multifab), allocatable :: Source_nph(:)
    type(multifab), allocatable :: thermal(:)

    ! Only needed for spherical.eq.1 
    type(multifab) , allocatable :: div_coeff_3d(:)
    real(kind=dp_t), pointer     :: dp(:,:,:,:)

    real (dp_t), allocatable :: grav_cell_nph(:)
    real (dp_t), allocatable :: grav_cell_new(:)

    real(dp_t)    , allocatable ::        s0_nph(:,:)
    real(dp_t)    , allocatable ::      w0_force(:)
    real(dp_t)    , allocatable ::        w0_old(:)
    real(dp_t)    , allocatable ::          Sbar(:,:)
    real(dp_t)    , allocatable :: div_coeff_nph(:)
    real(dp_t)    , allocatable :: div_coeff_edge(:)
    real(dp_t)    , allocatable :: rho_omegadotbar1(:,:)
    real(dp_t)    , allocatable :: rho_omegadotbar2(:,:)
    real(dp_t)    , allocatable :: rho_Hextbar(:,:)
    type(box)      ::  fine_domain
    real(dp_t)     :: halfdt, eps_in
    integer :: i,j,n,dm,nlevs
    integer :: nr,ng_cell,proj_type
    integer, allocatable :: lo(:),hi(:)
    logical :: nodal(mla%dim)

    nlevs = size(uold)
    dm    = mla%dim

    allocate(lo(dm),hi(dm))

    ng_cell = uold(1)%ng

    halfdt = half * dt

    allocate(Source_nph(nlevs))

    allocate(rhohalf(nlevs))
    allocate(w0_cart_vec(nlevs))
    allocate(w0_force_cart_vec(nlevs))
    allocate(macrhs(nlevs))
    allocate(macphi(nlevs))
    allocate( hgrhs_old(nlevs))
    allocate(thermal(nlevs))
    
    ! nr is the number of zones in a cell-centered basestate quantity
    nr    = size(s0_old,dim=1)

    allocate(          s0_nph(0:nr-1,nscal))
    allocate(            Sbar(0:nr-1,1))

    allocate(   div_coeff_nph(0:nr-1))
    allocate(  div_coeff_edge(0:nr))

    allocate(   grav_cell_nph(0:nr-1))
    allocate(   grav_cell_new(0:nr-1))

    allocate(rho_omegadotbar1(0:nr-1,nspec))
    allocate(rho_omegadotbar2(0:nr-1,nspec))
    allocate(     rho_Hextbar(0:nr-1,1))

    if (spherical.eq.1) &
      allocate(div_coeff_3d(nlevs))

    allocate(w0_force(0:nr-1))

    allocate(w0_old(0:nr))

    ! Set w0_old to w0 from last time step.
    w0_old(:) = w0(:)

    nodal = .true.
    do n = 1,nlevs
       call multifab_build(   rhohalf(n), mla%la(n),     1, 1)
       call multifab_build(Source_nph(n), mla%la(n),     1, 0)
       call multifab_build(    macrhs(n), mla%la(n),     1, 0)
       call multifab_build(    macphi(n), mla%la(n),     1, 1)
       call multifab_build( hgrhs_old(n), mla%la(n),     1, 0, nodal)
       call multifab_build(   thermal(n), mla%la(n),     1, 0)
       
       call setval(rhohalf(n),ZERO,all=.true.)
       call setval(Source_nph(n),ZERO,all=.true.)
       call setval(macrhs(n),ZERO,all=.true.)
       call setval(macphi(n),ZERO,all=.true.)
       call setval(hgrhs_old(n),ZERO,all=.true.)
       call setval(thermal(n),ZERO,all=.true.)
       
       if (dm.eq.3) then
          call multifab_build(      w0_cart_vec(n), mla%la(n),dm,1)
          call multifab_build(w0_force_cart_vec(n), mla%la(n),dm,1)
          call setval(w0_cart_vec(n),ZERO,all=.true.)
          call setval(w0_force_cart_vec(n),ZERO,all=.true.)
       end if

       if (spherical.eq.1) then
          call multifab_build(div_coeff_3d(n),mla%la(nlevs),1,1)
          call setval(div_coeff_3d(n),ZERO,all=.true.)
       endif
       
    end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 1 -- define average expansion at time n+1/2
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< CALLING advance_timestep with dt =',dt 
          write(6,*) '<<< STEP  1 : make w0 >>> '
        end if

        do n = 1, nlevs

           if (init_mode) then
              call make_S_at_halftime(Source_nph(n),Source_old(n),Source_new(n))
           else
              call extrap_to_halftime(Source_nph(n),dSdt(n),Source_old(n),dt)
           endif

        end do

        call average(Source_nph,Sbar,dx,1,1)
    
        call make_w0(w0,w0_old,w0_force,Sbar(:,1),p0_old,s0_old(:,rho_comp),gam1,dt,dtold,verbose)

        if (dm .eq. 3) then
          do n = 1, nlevs
             call make_w0_cart(w0      ,w0_cart_vec(n),normal(n),dx(n,:)) 
             call make_w0_cart(w0_force,w0_force_cart_vec(n),normal(n),dx(n,:)) 
          end do
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 2 -- construct the advective velocity
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  2 : create MAC velocities>>> '
        end if

        do n = 1,nlevs

           call advance_premac(uold(n), sold(n),&
                               umac(n,:), uedge(n,:), utrans(n,:),&
                               gp(n), normal(n), w0, w0_cart_vec(n), &
                               s0_old, grav_cell_old, &
                               dx(n,:),dt,the_bc_tower%bc_tower_array(n))
        end do

        do n = 1, nlevs
           call make_macrhs(macrhs(n),Source_nph(n),gamma1_term(n),Sbar(:,1),div_coeff_old,dx(n,:))
        end do

        ! MAC projection !
        if (spherical .eq. 1) then
          do n = 1, nlevs
            do i = 1,div_coeff_3d(n)%nboxes
              if (multifab_remote(div_coeff_3d(n),i)) cycle
              dp => dataptr(div_coeff_3d(n), i)
              lo =  lwb(get_box(div_coeff_3d(n), i))
              hi =  upb(get_box(div_coeff_3d(n), i))
              call fill_3d_data(dp(:,:,:,1),div_coeff_old,lo,hi,dx(nlevs,:),1)
            end do
            call multifab_fill_boundary(div_coeff_3d(n))
          end do
          call macproject(mla,umac,macphi,sold,dx,the_bc_tower,verbose,mg_verbose,cg_verbose,press_comp,&
                          macrhs,div_coeff_3d=div_coeff_3d)
        else
          call cell_to_edge(div_coeff_old,div_coeff_edge)
          call macproject(mla,umac,macphi,sold,dx,the_bc_tower,verbose,mg_verbose,cg_verbose,press_comp,&
                          macrhs,div_coeff_1d=div_coeff_old,div_coeff_half_1d=div_coeff_edge)
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 3 -- react the full state and then base state through dt/2
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  3 : react state     '
          write(6,*) '            : react  base >>> '
        end if

        do n = 1,nlevs
           call react_state(sold(n),s1(n),rho_omegadot1(n),rho_Hext(n), &
                            halfdt,dx(n,:), &
                            the_bc_tower%bc_tower_array(n),time)
        end do

        call average(rho_omegadot1,rho_omegadotbar1,dx,1,nspec)
        call average(rho_Hext,rho_Hextbar,dx,1,1)
        if (evolve_base_state) then
          call react_base(p0_old,s0_old,rho_omegadotbar1,rho_Hextbar(:,1),halfdt,p0_1,s0_1,gam1)
        else
          p0_1(:  ) = p0_old(:  )
          s0_1(:,:) = s0_old(:,:)
        end if
        call make_grav_cell(grav_cell_new,s0_1(:,rho_comp))
        call make_div_coeff(div_coeff_new,s0_1(:,rho_comp),p0_1, &
                            gam1,grav_cell_new,anelastic_cutoff)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 4 -- advect the base state and full state through dt
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  4 : advect base        '
          write(6,*) '            : scalar_advance >>> '
        end if

        if (evolve_base_state) then
          call advect_base(w0,Sbar(:,1),p0_1,p0_2,s0_1,s0_2,gam1, &
                           div_coeff_new,dx(1,dm),dt,anelastic_cutoff)
        else
          p0_2(:  ) = p0_1(:  )
          s0_2(:,:) = s0_1(:,:)
        end if
        
        if(use_thermal_diffusion) then
           call make_explicit_thermal(mla,dx,thermal,s1,p0_1, &
                                      mg_verbose,cg_verbose,the_bc_tower, &
                                      temp_diffusion_formulation)
        else
          do n = 1,nlevs
             call setval(thermal(n),ZERO)
          end do
        endif

        do n=1,nlevs
           if(istep .le. 1) then
              call add_react_to_thermal(thermal(n),rho_omegadot1(n),s1(n))
           else
              call add_react_to_thermal(thermal(n),rho_omegadot2(n),s1(n))
           endif
        enddo

        do n = 1,nlevs
           call scalar_advance (1, uold(n), s1(n), s2(n), thermal(n),&
                                umac(n,:), w0, w0_cart_vec(n), sedge(n,:), utrans(n,:),&
                                scal_force(n), normal(n), s0_1, s0_2, p0_1, p0_2, &
                                dx(n,:),dt, &
                                the_bc_tower%bc_tower_array(n), &
                                verbose)
        end do

        do n = 2, nlevs
           fine_domain = layout_get_pd(mla%la(n))
           call multifab_fill_ghost_cells(s2(n),s2(n-1),fine_domain, &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,rho_comp,nscal)
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 4a (Option I) -- Add thermal conduction (only enthalpy terms)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (use_thermal_diffusion) then
           if (parallel_IOProcessor() .and. verbose .ge. 1) then
              write(6,*) '<<< STEP  4a: thermal conduct >>>'
           end if

           if(do_half_alg) then
              call thermal_conduct_half_alg(mla,dx,dt,s1,s2,p0_1,p0_2, &
                                            s0_1(:,temp_comp), s0_2(:,temp_comp), &
                                            mg_verbose,cg_verbose,the_bc_tower)
           else
              
           endif
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 5 -- react the full state and then base state through dt/2
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  5 : react state     '
          write(6,*) '            : react  base >>> '
        end if

        do n = 1,nlevs
          call react_state(s2(n),snew(n),rho_omegadot2(n),rho_Hext(n), &
                           halfdt,dx(n,:), &
                           the_bc_tower%bc_tower_array(n),time)
        end do
        call average(rho_omegadot2,rho_omegadotbar2,dx,1,nspec)
        call average(rho_Hext,rho_Hextbar,dx,1,1)
        if (evolve_base_state) then
          call react_base(p0_2,s0_2,rho_omegadotbar2,rho_Hextbar(:,1),halfdt,p0_new,s0_new,gam1)
        else
          p0_new(:  ) = p0_2(:  )
          s0_new(:,:) = s0_2(:,:)
        end if
        call make_grav_cell(grav_cell_new,s0_new(:,rho_comp))
        call make_div_coeff(div_coeff_new,s0_new(:,rho_comp),p0_new, &
                            gam1,grav_cell_new,anelastic_cutoff)

        ! Define rho at half time !
        do n = 1,nlevs
           call make_at_halftime(rhohalf(n),sold(n),snew(n),rho_comp,1,dx(n,:), &
                                 the_bc_tower%bc_tower_array(n))
        end do

        ! Define base state at half time for use in velocity advance!
        do j = 0, nr-1
           s0_nph(j,:) = HALF * (s0_old(j,:) + s0_new(j,:))
        enddo

        call make_grav_cell(grav_cell_nph,s0_nph(:,rho_comp))

        ! Define beta at half time !
        do j = 0, nr-1
           div_coeff_nph(j) = HALF * (div_coeff_old(j) + div_coeff_new(j))
        enddo

        if(.not. do_half_alg) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 6 -- define a new average expansion rate at n+1/2
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  6 : make new S and new w0 >>> '
        end if

        if(use_thermal_diffusion) then
           call make_explicit_thermal(mla,dx,thermal,snew,p0_new, &
                                      mg_verbose,cg_verbose,the_bc_tower, &
                                      temp_diffusion_formulation)
        else
          do n = 1,nlevs
             call setval(thermal(n),ZERO)
          end do
        endif

        do n = 1, nlevs

           call make_S(Source_new(n),gamma1_term(n),snew(n), &
                       rho_omegadot2(n),rho_Hext(n),thermal(n), &
                       s0_old(:,temp_comp),gam1,dx(n,:))
           call make_S_at_halftime(Source_nph(n),Source_old(n),Source_new(n))
           call average(Source_nph,Sbar,dx,1,1)

        end do

        call make_w0(w0,w0_old,w0_force,Sbar(:,1),p0_new,s0_new(:,rho_comp),gam1,dt,dtold,verbose)

        if (dm .eq. 3) then
           do n = 1, nlevs
             call make_w0_cart(w0,w0_cart_vec(n),normal(n),dx(n,:)) 
             call make_w0_cart(w0_force,w0_force_cart_vec(n),normal(n),dx(n,:)) 
           end do
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 7 -- redo the construction of the advective velocity using the current w0
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  7 : create MAC velocities >>> '
        end if

        do n = 1,nlevs
           call advance_premac(uold(n), sold(n),&
                               umac(n,:), uedge(n,:), utrans(n,:),&
                               gp(n),  normal(n), w0, w0_cart_vec(n), &
                               s0_old, grav_cell_old, &
                               dx(n,:),dt,the_bc_tower%bc_tower_array(n))
        end do

        do n = 1, nlevs
           call make_macrhs(macrhs(n),Source_nph(n),gamma1_term(n),Sbar(:,1),div_coeff_nph,dx(n,:))
        end do

        ! MAC projection !
        if (spherical .eq. 1) then
          do n = 1, nlevs
            do i = 1,div_coeff_3d(n)%nboxes
              if (multifab_remote(div_coeff_3d(n),i)) cycle
              dp => dataptr(div_coeff_3d(n), i)
              lo =  lwb(get_box(div_coeff_3d(n), i))
              hi =  upb(get_box(div_coeff_3d(n), i))
              call fill_3d_data(dp(:,:,:,1),div_coeff_nph,lo,hi,dx(nlevs,:),1)
            end do
            call multifab_fill_boundary(div_coeff_3d(n))
          end do
          call macproject(mla,umac,macphi,rhohalf,dx,the_bc_tower,verbose,mg_verbose,cg_verbose,&
                          press_comp,macrhs,div_coeff_3d=div_coeff_3d)
        else
          call cell_to_edge(div_coeff_nph,div_coeff_edge)
          call macproject(mla,umac,macphi,rhohalf,dx,the_bc_tower,verbose,mg_verbose,cg_verbose,&
                          press_comp,macrhs,div_coeff_1d=div_coeff_nph,div_coeff_half_1d=div_coeff_edge)
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 8 -- advect the base state and full state through dt
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  8 : advect base   '
          write(6,*) '            : scalar_advance >>>'
        end if
        if (evolve_base_state) then
          call advect_base(w0,Sbar(:,1),p0_1,p0_2,s0_1,s0_2,gam1, &
                           div_coeff_nph,dx(1,dm),dt,anelastic_cutoff)
        else
          p0_2(:  ) = p0_1(:  )
          s0_2(:,:) = s0_1(:,:)
        end if

        if(use_thermal_diffusion) then
           call make_explicit_thermal(mla,dx,thermal,s1,p0_1, &
                                      mg_verbose,cg_verbose,the_bc_tower, &
                                      temp_diffusion_formulation)
        else
          do n = 1,nlevs
             call setval(thermal(n),ZERO)
          end do
        endif

         do n=1,nlevs
            if(istep .le. 1) then
               call add_react_to_thermal(thermal(n),rho_omegadot1(n),s1(n))
            else
               call add_react_to_thermal(thermal(n),rho_omegadot2(n),s1(n))
            endif
         enddo

        do n = 1,nlevs
           call scalar_advance (2, uold(n), s1(n), s2(n), thermal(n), &
                                umac(n,:), w0, w0_cart_vec(n), sedge(n,:), utrans(n,:),&
                                scal_force(n), normal(n), s0_1, s0_2, p0_1, p0_2, &
                                dx(n,:),dt, &
                                the_bc_tower%bc_tower_array(n), &
                                verbose)
        end do

        do n = 2, nlevs
           fine_domain = layout_get_pd(mla%la(n))
           call multifab_fill_ghost_cells(s2(n),s2(n-1),fine_domain, &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,rho_comp,nscal)
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 8a (Option I) -- Add thermal conduction (only enthalpy terms)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!        if (use_thermal_diffusion) then
!           if (parallel_IOProcessor() .and. verbose .ge. 1) then
!              write(6,*) '<<< STEP  8a: thermal conduct >>>'
!           end if
!
!        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 9 -- react the full state and then base state through dt/2
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  9 : react state '
          write(6,*) '            : react  base >>>'
        end if
        do n = 1,nlevs
          call react_state(s2(n),snew(n),rho_omegadot2(n),rho_Hext(n), &
                           halfdt,dx(n,:),&
                           the_bc_tower%bc_tower_array(n),time)
        end do
        call average(rho_omegadot2,rho_omegadotbar2,dx,1,nspec)
        call average(rho_Hext,rho_Hextbar,dx,1,1)
        if (evolve_base_state) then
          call react_base(p0_2,s0_2,rho_omegadotbar2,rho_Hextbar(:,1),halfdt,p0_new,s0_new,gam1)
        else
          p0_new(:  ) = p0_2(:  )
          s0_new(:,:) = s0_2(:,:)
        end if
        call make_grav_cell(grav_cell_new,s0_new(:,rho_comp))
        call make_div_coeff(div_coeff_new,s0_new(:,rho_comp),p0_new, &
                            gam1,grav_cell_new,anelastic_cutoff)


        ! endif corresponding to .not. do_half_alg
        endif

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 10 -- compute S^{n+1} for the final projection
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP 10 : make new S >>>'
        end if

        if(use_thermal_diffusion) then
           call make_explicit_thermal(mla,dx,thermal,snew,p0_new, &
                                      mg_verbose,cg_verbose,the_bc_tower, &
                                      temp_diffusion_formulation)
        else
          do n = 1,nlevs
             call setval(thermal(n),ZERO)
          end do
        endif

        do n = 1, nlevs
           call make_S(Source_new(n),gamma1_term(n),snew(n), &
                       rho_omegadot2(n),rho_Hext(n),thermal(n), &
                       s0_new(:,temp_comp),gam1,dx(n,:))
        end do
        call average(Source_new,Sbar,dx,1,1)

        ! define dSdt = (Source_new - Source_old) / dt
        do n = 1,nlevs
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
        do n = 1,nlevs
           call make_at_halftime(rhohalf(n),sold(n),snew(n),rho_comp,1,dx(n,:), &
                                 the_bc_tower%bc_tower_array(n))
        end do

        do n = 1,nlevs
           call velocity_advance(uold(n),unew(n),sold(n),rhohalf(n),&
                                 umac(n,:),uedge(n,:), &
                                 utrans(n,:),gp(n),p(n), &
                                 normal(n), w0, w0_cart_vec(n), &
                                 w0_force, w0_force_cart_vec(n), &
                                 s0_old, grav_cell_old, s0_nph, grav_cell_nph, &
                                 dx(n,:),dt, &
                                 the_bc_tower%bc_tower_array(n), &
                                 sponge(n),do_sponge,verbose)
        end do

        do n = 2, nlevs
           fine_domain = layout_get_pd(mla%la(n))
           call multifab_fill_ghost_cells(unew(n),unew(n-1),fine_domain, &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,1,dm)
        end do

        ! Define beta at half time using the div_coeff_new from step 9!
        do j = 0, nr-1
          div_coeff_nph(j) = HALF * (div_coeff_old(j) + div_coeff_new(j))
        end do

!       Project the new velocity field.
        if (init_mode) then
          proj_type = pressure_iters
          do n = 1,nlevs
             call multifab_copy(hgrhs_old(n),hgrhs(n))
             call make_hgrhs(hgrhs(n),Source_new(n),gamma1_term(n),Sbar(:,1),div_coeff_new,dx(n,:))
             call multifab_sub_sub(hgrhs(n),hgrhs_old(n))
             call multifab_div_div_s(hgrhs(n),dt)
          end do
        else
          proj_type = regular_timestep
          do n = 1,nlevs
             call make_hgrhs(hgrhs(n),Source_new(n),gamma1_term(n),Sbar(:,1),div_coeff_new,dx(n,:))
          end do
        end if

        if (spherical .eq. 1) then
          do n = 1,nlevs
            do i = 1,div_coeff_3d(n)%nboxes
              if (multifab_remote(div_coeff_3d(n),i)) cycle
              dp => dataptr(div_coeff_3d(n), i)
              lo =  lwb(get_box(div_coeff_3d(n), i))
              hi =  upb(get_box(div_coeff_3d(n), i))
              call fill_3d_data(dp(:,:,:,1),div_coeff_nph,lo,hi,dx(nlevs,:),1)
            end do
          end do
          eps_in = 1.d-12
          call hgproject(proj_type, mla, unew, uold, rhohalf, p, gp, dx, dt, &
                         the_bc_tower, verbose, mg_verbose, cg_verbose, press_comp, &
                         hgrhs, div_coeff_3d=div_coeff_3d, eps_in = eps_in)
        else
          call hgproject(proj_type, mla, unew, uold, rhohalf, p, gp, dx, dt, &
                         the_bc_tower, verbose, mg_verbose, cg_verbose, press_comp, &
                         hgrhs, div_coeff_1d=div_coeff_nph)
        end if

!       If doing pressure iterations then put hgrhs_old into hgrhs to be returned to varden.
        if (init_mode) then
          do n = 1,nlevs
            call multifab_copy(hgrhs(n),hgrhs_old(n))
          end do
        end if

        do n = 2, nlevs
           fine_domain = layout_get_pd(mla%la(n))
           call multifab_fill_ghost_cells(unew(n),unew(n-1),fine_domain, &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,1,dm)
        end do

        do n = 1, nlevs
          call destroy(Source_nph(n))
          call destroy(macrhs(n))
          call destroy(macphi(n))
          call destroy( hgrhs_old(n))
          call destroy(thermal(n))
          call destroy(rhohalf(n))
          if (spherical .eq. 1) &
            call destroy(div_coeff_3d(n))
        end do
        deallocate(Source_nph)
        deallocate(macrhs)
        deallocate(macphi)
        deallocate( hgrhs_old)
        deallocate(thermal)
        deallocate(rhohalf)

        if (spherical .eq. 1) &
          deallocate(div_coeff_3d)
       
        if (dm .eq. 3) then
          do n = 1, nlevs
            call destroy(w0_cart_vec(n))
            call destroy(w0_force_cart_vec(n))
          end do
          deallocate(w0_cart_vec)
          deallocate(w0_force_cart_vec)
        end if

        deallocate(Sbar)
        deallocate(s0_nph)
        deallocate(div_coeff_nph)
        deallocate(div_coeff_edge)
        deallocate(grav_cell_nph)
        deallocate(grav_cell_new)

        deallocate(w0_old,w0_force)

        deallocate(lo,hi)

    end subroutine advance_timestep

end module advance_timestep_module
