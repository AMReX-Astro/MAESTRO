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
  use variables
  use network

  contains

    subroutine advance_timestep(init_mode, &
                                mla,uold,sold,s1,s2,unew,snew,umac,uedge,sedge,utrans,gp,p, &
                                scal_force,normal, &
                                s0_old,s0_1,s0_2,s0_new,p0_old,p0_1,p0_2,p0_new,temp0,gam1,w0, &
                                rho_omegadot1, rho_omegadot2, rho_Hext, &
                                div_coeff_old,div_coeff_new,&
                                dx,time,dt,dtold,the_bc_tower, &
                                anelastic_cutoff,verbose,mg_verbose,cg_verbose,&
                                Source_nm1,Source_old,Source_new,gamma1_term)

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
    type(multifab), intent(inout) :: Source_nm1(:)
    type(multifab), intent(inout) :: Source_old(:)
    type(multifab), intent(inout) :: Source_new(:)
    type(multifab), intent(inout) :: gamma1_term(:)
    real(dp_t)    , intent(inout) :: s0_old(:,:)
    real(dp_t)    , intent(inout) :: s0_1(:,:)
    real(dp_t)    , intent(inout) :: s0_2(:,:)
    real(dp_t)    , intent(inout) :: s0_new(:,:)
    real(dp_t)    , intent(inout) :: p0_old(:)
    real(dp_t)    , intent(inout) :: p0_1(:)
    real(dp_t)    , intent(inout) :: p0_2(:)
    real(dp_t)    , intent(inout) :: p0_new(:)
    real(dp_t)    , intent(inout) :: temp0(:)
    real(dp_t)    , intent(inout) :: gam1(:)
    real(dp_t)    , intent(inout) :: w0(:)
    real(dp_t)    , intent(in   ) :: div_coeff_old(:)
    real(dp_t)    , intent(inout) :: div_coeff_new(:)
    real(dp_t)    , intent(in   ) :: dx(:,:), time, dt, dtold
    type(bc_tower), intent(in   ) :: the_bc_tower
    real(dp_t)    , intent(in   ) :: anelastic_cutoff
    integer       , intent(in   ) :: verbose,mg_verbose,cg_verbose

    type(multifab), allocatable :: rhohalf(:)
    type(multifab), allocatable :: w0_cart_vec(:)
    type(multifab), allocatable :: macrhs(:)
    type(multifab), allocatable ::  hgrhs(:)
    type(multifab), allocatable :: Source_nph(:)

    ! Only needed for spherical.eq.1 
    type(multifab) , allocatable :: div_coeff_3d(:)
    real(kind=dp_t), pointer     :: dp(:,:,:,:)

    real(dp_t)    , allocatable ::        s0_nph(:,:)
    real(dp_t)    , allocatable ::          Sbar(:,:)
    real(dp_t)    , allocatable :: div_coeff_nph(:)
    real(dp_t)    , allocatable :: div_coeff_edge(:)
    real(dp_t)    , allocatable ::      grav_cell(:)
    real(dp_t)    , allocatable :: rho_omegadotbar1(:,:)
    real(dp_t)    , allocatable :: rho_omegadotbar2(:,:)
    real(dp_t)    , allocatable :: rho_Hextbar(:,:)
    type(bc_level) ::  bc
    type(box)      ::  fine_domain
    real(dp_t)     :: halfdt, half_time, new_time, eps_in
    integer :: i,j,n,dm,nscal,nlevs,comp
    integer :: nr,ng_cell
    integer, allocatable :: lo(:),hi(:)
    logical :: nodal(mla%dim)

    nlevs = size(uold)
    dm    = mla%dim

    allocate(lo(dm),hi(dm))

    ng_cell = uold(1)%ng

    halfdt = half * dt
    half_time = dt + halfdt
     new_time = dt +      dt

    allocate(Source_nph(nlevs))

    allocate(rhohalf(nlevs))
    allocate(w0_cart_vec(nlevs))
    allocate(macrhs(nlevs))
    allocate( hgrhs(nlevs))

    nscal = multifab_ncomp(sold(1))
    nr    = size(s0_old,dim=1)

    allocate(          s0_nph(nr,nscal))
    allocate(            Sbar(nr,1))
    allocate(   div_coeff_nph(nr))
    allocate(  div_coeff_edge(nr+1))
    allocate(       grav_cell(nr))
    allocate(rho_omegadotbar1(nr,nspec))
    allocate(rho_omegadotbar2(nr,nspec))
    allocate(rho_Hextbar(nr,1))

    if (spherical.eq.1) &
      allocate(div_coeff_3d(nlevs))

    nodal = .true.
    do n = 1,nlevs
     call multifab_build(   rhohalf(n), mla%la(n),     1, 1)
     call multifab_build(Source_nph(n), mla%la(n),     1, 0)
     call multifab_build(    macrhs(n), mla%la(n),     1, 0)
     call multifab_build(     hgrhs(n), mla%la(n),     1,       0, nodal)
     if (dm.eq.3) &
         call multifab_build(w0_cart_vec(n), mla%la(n),dm,1)
     if (spherical.eq.1) &
         call multifab_build(div_coeff_3d(n),mla%la(nlevs),1,0)
    end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 1 -- define average expansion at time n+1/2
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  1 : make w0 >>> '
        end if

        do n = 1, nlevs

           if (init_mode) then
              call make_S_at_halftime(Source_nph(n),Source_old(n),Source_new(n))
           else
              call extrap_to_halftime(Source_nph(n),Source_nm1(n),Source_old(n),dtold,dt)
           endif

           call average(Source_nph(n),Sbar,dx(n,:))
    
           call make_w0(w0,Sbar(:,1),p0_old,s0_old(:,rho_comp),temp0,gam1,dx(n,dm),dt,verbose)
           if (dm .eq. 3) then
             call make_w0_cart(w0,w0_cart_vec(n),normal(n),dx(n,:)) 
           end if
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 2 -- construct the advective velocity
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  2 : create MAC velocities>>> '
        end if

        do n = 1,nlevs

           call advance_premac(uold(n), sold(n),&
                               umac(n,:), uedge(n,:), utrans(n,:),&
                               gp(n), p(n), normal(n), w0, w0_cart_vec(n), s0_old, &
                               dx(n,:),time,dt, &
                               the_bc_tower%bc_tower_array(n), &
                               verbose)
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
              call fill_3d_data(dp(:,:,:,1),div_coeff_old,lo,hi,dx(nlevs,:),0)
            end do
          end do
          call macproject(mla,umac,sold,dx,the_bc_tower,verbose,mg_verbose,cg_verbose,press_comp,&
                          macrhs,div_coeff_3d=div_coeff_3d)
        else
          call put_1d_beta_on_edges(div_coeff_old,div_coeff_edge)
          call macproject(mla,umac,sold,dx,the_bc_tower,verbose,mg_verbose,cg_verbose,press_comp,&
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
           call react_state(sold(n),s1(n),rho_omegadot1(n),rho_Hext(n),temp0, &
                            halfdt,dx(n,:), &
                            the_bc_tower%bc_tower_array(n),time)
           call multifab_fill_boundary(s1(n))
        end do

        call average(rho_omegadot1(1),rho_omegadotbar1,dx(1,:))
        call average(rho_Hext(1),rho_Hextbar,dx(1,:))
        call react_base(p0_old,s0_old,temp0,rho_omegadotbar1,rho_Hextbar(:,1),dx(1,dm),halfdt,p0_1,s0_1,gam1)
        call make_grav_cell(grav_cell,s0_1(:,rho_comp))
        call make_div_coeff(div_coeff_new,s0_1(:,rho_comp),p0_1, &
                            gam1,grav_cell,dx(1,dm),anelastic_cutoff)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 4 -- advect the base state and full state through dt
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  4 : advect base        '
          write(6,*) '            : scalar_advance >>> '
        end if

        call put_1d_beta_on_edges(div_coeff_new,div_coeff_edge)
        call advect_base(w0,Sbar(:,1),p0_1,p0_2,s0_1,s0_2,temp0,gam1, &
                         div_coeff_edge,&
                         dx(1,dm),dt,anelastic_cutoff)
        call make_grav_cell(grav_cell,s0_2(:,rho_comp))
        call make_div_coeff(div_coeff_new,s0_2(:,rho_comp),p0_2, &
                            gam1,grav_cell,dx(1,dm),anelastic_cutoff)

        do n = 1,nlevs
           call scalar_advance (uold(n), s1(n), s2(n), &
                                umac(n,:), w0, w0_cart_vec(n), sedge(n,:), utrans(n,:),&
                                scal_force(n), normal(n), s0_1, s0_2, p0_1, p0_2, temp0, &
                                dx(n,:),time,dt, &
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

        do n = 1,nlevs
           call multifab_fill_boundary(s2(n))
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 5 -- react the full state and then base state through dt/2
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  5 : react state     '
          write(6,*) '            : react  base >>> '
        end if

        do n = 1,nlevs
          call react_state(s2(n),snew(n),rho_omegadot2(n),rho_Hext(n),temp0, &
                           halfdt,dx(n,:), &
                           the_bc_tower%bc_tower_array(n),time)
          call multifab_fill_boundary(snew(n))
        end do
        call average(rho_omegadot2(1),rho_omegadotbar2,dx(1,:))
        call average(rho_Hext(1),rho_Hextbar,dx(1,:))
        call react_base(p0_2,s0_2,temp0,rho_omegadotbar2,rho_Hextbar(:,1),dx(1,dm),halfdt,p0_new,s0_new,gam1)
        call make_grav_cell(grav_cell,s0_new(:,rho_comp))
        call make_div_coeff(div_coeff_new,s0_new(:,rho_comp),p0_new, &
                            gam1,grav_cell,dx(1,dm),anelastic_cutoff)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 6 -- define a new average expansion rate at n+1/2
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  6 : make new S and new w0 >>> '
        end if

        do n = 1, nlevs
           call make_S(Source_new(n),gamma1_term(n),snew(n),uold(n), &
                       rho_omegadot2(n),rho_Hext(n),p0_new,temp0,gam1,dx(n,:),time)
           call make_S_at_halftime(Source_nph(n),Source_old(n),Source_new(n))
           call average(Source_nph(n),Sbar,dx(n,:))

           call make_w0(w0,Sbar(:,1),p0_new,s0_new(:,rho_comp),temp0,gam1,dx(n,dm),dt,verbose)
           if (dm .eq. 3) then
             call make_w0_cart(w0,w0_cart_vec(n),normal(n),dx(n,:)) 
           end if
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 7 -- redo the construction of the advective velocity using the current w0
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  7 : create MAC velocities >>> '
        end if

        do n = 1,nlevs
           call advance_premac(uold(n), sold(n),&
                               umac(n,:), uedge(n,:), utrans(n,:),&
                               gp(n), p(n), normal(n), w0, w0_cart_vec(n), s0_old, &
                               dx(n,:),time,dt, &
                               the_bc_tower%bc_tower_array(n), &
                               verbose)
        end do

        ! Define rho at half time !
        do n = 1,nlevs
           call make_at_halftime(rhohalf(n),sold(n),snew(n),rho_comp,1,dx(n,:), &
                                 the_bc_tower%bc_tower_array(n))
           call multifab_fill_boundary(rhohalf(n))
        end do

        ! Define base state at half time for use in velocity advance!
   
        do j = 1,size(s0_nph,dim=1)
          s0_nph(j,:) = HALF * (s0_old(j,:) + s0_new(j,:))
        end do

        ! Define beta at half time !
        do j = 1,size(div_coeff_old,dim=1)
          div_coeff_nph(j) = HALF * (div_coeff_old(j) + div_coeff_new(j))
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
              call fill_3d_data(dp(:,:,:,1),div_coeff_nph,lo,hi,dx(nlevs,:),0)
            end do
          end do
          call macproject(mla,umac,rhohalf,dx,the_bc_tower,verbose,mg_verbose,cg_verbose,&
                          press_comp,macrhs,div_coeff_3d=div_coeff_3d)
        else
          call put_1d_beta_on_edges(div_coeff_nph,div_coeff_edge)
          call macproject(mla,umac,rhohalf,dx,the_bc_tower,verbose,mg_verbose,cg_verbose,&
                          press_comp,macrhs,div_coeff_1d=div_coeff_nph,div_coeff_half_1d=div_coeff_edge)
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 8 -- advect the base state and full state through dt
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  8 : advect base   '
          write(6,*) '            : scalar_advance >>>'
        end if
        call advect_base(w0,Sbar(:,1),p0_1,p0_2,s0_1,s0_2,temp0,gam1, &
                         div_coeff_edge,&
                         dx(1,dm),dt,anelastic_cutoff)
        call make_grav_cell(grav_cell,s0_2(:,rho_comp))
        call make_div_coeff(div_coeff_new,s0_2(:,rho_comp),p0_2, &
                            gam1,grav_cell,dx(1,dm),anelastic_cutoff)

        do n = 1,nlevs
           call scalar_advance (uold(n), s1(n), s2(n), &
                                umac(n,:), w0, w0_cart_vec(n), sedge(n,:), utrans(n,:),&
                                scal_force(n), normal(n), s0_1, s0_2, p0_1, p0_2, temp0, &
                                dx(n,:),time,dt, &
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

        do n = 1,nlevs
           call multifab_fill_boundary(s2(n))
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 9 -- react the full state and then base state through dt/2
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP  9 : react state '
          write(6,*) '            : react  base >>>'
        end if
        do n = 1,nlevs
          call react_state(s2(n),snew(n),rho_omegadot2(n),rho_Hext(n),temp0, &
                           halfdt,dx(n,:),&
                           the_bc_tower%bc_tower_array(n),time)
          call multifab_fill_boundary(snew(n))
        end do
        call average(rho_omegadot2(1),rho_omegadotbar2,dx(1,:))
        call average(rho_Hext(1),rho_Hextbar,dx(1,:))
        call react_base(p0_2,s0_2,temp0,rho_omegadotbar2,rho_Hextbar(:,1),dx(1,dm),halfdt,p0_new,s0_new,gam1)
        call make_grav_cell(grav_cell,s0_new(:,rho_comp))
        call make_div_coeff(div_coeff_new,s0_new(:,rho_comp),p0_new, &
                            gam1,grav_cell,dx(1,dm),anelastic_cutoff)


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 10 -- compute S^{n+1} for the final projection
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) '<<< STEP 10 : make new S >>>'
        end if

        do n = 1, nlevs
           call make_S(Source_new(n),gamma1_term(n),snew(n),uold(n), &
                       rho_omegadot2(n),rho_Hext(n),p0_new,temp0,gam1,dx(n,:),time)
           call average(Source_new(n),Sbar,dx(n,:))
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
           call multifab_fill_boundary(rhohalf(n))
        end do

        do n = 1,nlevs
           call velocity_advance(uold(n),unew(n),sold(n),rhohalf(n),&
                                 umac(n,:),uedge(n,:), &
                                 utrans(n,:),gp(n),p(n), &
                                 normal(n), w0, w0_cart_vec(n), s0_old, s0_nph, &
                                 dx(n,:),time,dt, &
                                 the_bc_tower%bc_tower_array(n), &
                                 verbose)
        end do

        do n = 2, nlevs
           fine_domain = layout_get_pd(mla%la(n))
           call multifab_fill_ghost_cells(unew(n),unew(n-1),fine_domain, &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,1,dm)
        end do

        do n = 1,nlevs
           call multifab_fill_boundary(unew(n))
        end do

        do n = 1,nlevs
           call make_hgrhs(hgrhs(n),Source_new(n),gamma1_term(n),Sbar(:,1),div_coeff_new,dx(n,:), &
                           layout_get_pd(mla%la(n)))
        end do

        ! Define beta at half time using the div_coeff_new from step 9!
        do j = 1,size(div_coeff_old,dim=1)
          div_coeff_nph(j) = HALF * (div_coeff_old(j) + div_coeff_new(j))
        end do

!       Project the new velocity field.
        if (spherical .eq. 1) then
          do n = 1,nlevs
            do i = 1,div_coeff_3d(n)%nboxes
              if (multifab_remote(div_coeff_3d(n),i)) cycle
              dp => dataptr(div_coeff_3d(n), i)
              lo =  lwb(get_box(div_coeff_3d(n), i))
              hi =  upb(get_box(div_coeff_3d(n), i))
              call fill_3d_data(dp(:,:,:,1),div_coeff_nph,lo,hi,dx(nlevs,:),0)
            end do
          end do
          eps_in = 1.d-10
          call hgproject(mla, unew, rhohalf, p, gp, dx, dt, &
                         the_bc_tower, verbose, mg_verbose, cg_verbose, press_comp, &
                         hgrhs, div_coeff_3d=div_coeff_3d, eps_in = eps_in)
        else
          call hgproject(mla, unew, rhohalf, p, gp, dx, dt, &
                         the_bc_tower, verbose, mg_verbose, cg_verbose, press_comp, &
                         hgrhs, div_coeff_1d=div_coeff_nph)
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
          call destroy( hgrhs(n))
          call destroy(rhohalf(n))
          if (spherical .eq. 1) &
            call destroy(div_coeff_3d(n))
        end do
        deallocate(Source_nph)
        deallocate(macrhs)
        deallocate( hgrhs)
        deallocate(rhohalf)

        if (spherical .eq. 1) &
          deallocate(div_coeff_3d)
       
        if (dm .eq. 3) then
          do n = 1, nlevs
            call destroy(w0_cart_vec(n))
          end do
          deallocate(w0_cart_vec)
        end if

        deallocate(  Sbar)
        deallocate(s0_nph)
        deallocate(div_coeff_nph)
        deallocate(div_coeff_edge)
        deallocate(grav_cell)

        deallocate(lo,hi)

    end subroutine advance_timestep

end module advance_timestep_module
