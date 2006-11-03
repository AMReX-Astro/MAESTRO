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
  use make_w0_module
  use advect_base_module
  use react_base_module
  use react_state_module
  use make_S_module
  use average_module
  use rhohalf_module
  use variables
  use network

  contains

    subroutine advance_timestep(mla,uold,sold,s1,s2,unew,snew,umac,uedge,sedge,utrans,gp,p, &
                                force,scal_force,&
                                s0_old,s0_1,s0_2,s0_new,s0_nph,p0_old,p0_1,p0_2,p0_new,temp0,gam1,w0, &
                                rho_omegadot1, rho_omegadot2, &
                                div_coeff_old,div_coeff_new,&
                                dx,time,dt,the_bc_tower, &
                                anelastic_cutoff,verbose,mg_verbose,cg_verbose,&
                                Source_old,Source_new,spherical)


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
    type(multifab), intent(inout) :: force(:)
    type(multifab), intent(inout) :: scal_force(:)
    type(multifab), intent(inout) :: rho_omegadot1(:)
    type(multifab), intent(inout) :: rho_omegadot2(:)
    type(multifab), intent(inout) :: Source_old(:)
    type(multifab), intent(inout) :: Source_new(:)
    real(dp_t)    , intent(inout) :: s0_old(:,:)
    real(dp_t)    , intent(inout) :: s0_1(:,:)
    real(dp_t)    , intent(inout) :: s0_2(:,:)
    real(dp_t)    , intent(inout) :: s0_new(:,:)
    real(dp_t)    , intent(inout) :: s0_nph(:,:)
    real(dp_t)    , intent(inout) :: p0_old(:)
    real(dp_t)    , intent(inout) :: p0_1(:)
    real(dp_t)    , intent(inout) :: p0_2(:)
    real(dp_t)    , intent(inout) :: p0_new(:)
    real(dp_t)    , intent(inout) :: temp0(:)
    real(dp_t)    , intent(inout) :: gam1(:)
    real(dp_t)    , intent(inout) :: w0(:)
    real(dp_t)    , intent(in   ) :: div_coeff_old(:)
    real(dp_t)    , intent(inout) :: div_coeff_new(:)
    real(dp_t)    , intent(in   ) :: dx(:,:), time, dt
    type(bc_tower), intent(in   ) :: the_bc_tower
    real(dp_t)    , intent(in   ) :: anelastic_cutoff
    integer       , intent(in   ) :: verbose,mg_verbose,cg_verbose
    integer       , intent(in   ) :: spherical

    type(multifab), allocatable :: rhohalf(:)
    type(multifab), allocatable :: macrhs(:)
    type(multifab), allocatable ::  hgrhs(:)
    type(multifab), allocatable :: Source_nph(:)
    real(dp_t)    , allocatable :: Sbar(:,:)
    real(dp_t)    , allocatable :: div_coeff_nph(:)
    real(dp_t)    , allocatable :: div_coeff_edge(:)
    real(dp_t)    , allocatable ::      grav_edge(:)
    real(dp_t)    , allocatable :: rho_omegadotbar1(:,:)
    real(dp_t)    , allocatable :: rho_omegadotbar2(:,:)
    type(bc_level) ::  bc
    type(box)      ::  fine_domain
    real(dp_t)     :: halfdt, half_time, new_time
    integer :: n,dm,nlevs,comp,bc_comp
    logical :: nodal(mla%dim)

    nlevs = size(uold)
    dm    = mla%dim

    halfdt = half * dt
    half_time = dt + halfdt
     new_time = dt +      dt

    allocate(Source_nph(nlevs))

    allocate(rhohalf(nlevs))
    allocate(macrhs(nlevs))
    allocate( hgrhs(nlevs))

    allocate(            Sbar(extent(mla%mba%pd(1),dm),1))
    allocate(   div_coeff_nph(extent(mla%mba%pd(1),dm)))
    allocate(  div_coeff_edge(extent(mla%mba%pd(1),dm)+1))
    allocate(       grav_edge(extent(mla%mba%pd(1),dm)+1))
    allocate(rho_omegadotbar1(extent(mla%mba%pd(1),dm),nspec))
    allocate(rho_omegadotbar2(extent(mla%mba%pd(1),dm),nspec))

    nodal = .true.
    do n = 1,nlevs
     call multifab_build(   rhohalf(n), mla%la(n),     1, 1)
     call multifab_build(Source_nph(n), mla%la(n),     1, 0)
     call multifab_build(    macrhs(n), mla%la(n),     1, 0)
     call multifab_build(     hgrhs(n), mla%la(n),     1,       0, nodal)
    end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 1 !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        print *,'<<< STEP 1 >>>'
        do n = 1, nlevs
           call make_S(Source_nph(n),sold(n),p0_old,temp0,gam1,dx(n,:),time)
           call average(Source_nph(n),Sbar)
           call make_w0(w0,Sbar(:,1),p0_old,s0_old(:,rho_comp),temp0,gam1,dx(n,dm),dt,spherical)
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 2 !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        print *,'<<< STEP 2 >>>'

        do n = 1,nlevs
           call advance_premac(uold(n), sold(n),&
                               umac(n,:), uedge(n,:), utrans(n,:),&
                               gp(n), p(n), force(n), &
                               s0_old, &
                               dx(n,:),time,dt, &
                               the_bc_tower%bc_tower_array(n), &
                               verbose)
        end do

        do n = 1, nlevs
           call make_macrhs(macrhs(n),Source_nph(n),Sbar(:,1),div_coeff_old)
        end do

        ! MAC projection !
        call put_beta_on_edges(div_coeff_old,div_coeff_edge)
        call macproject(mla,umac,sold,dx,the_bc_tower,verbose,mg_verbose,cg_verbose,press_comp,&
                        macrhs,div_coeff_old,div_coeff_edge)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 3 !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        print *,'<<< STEP 3 >>>'

        do n = 1,nlevs
          call react_state(sold(n),s1(n),temp0,rho_omegadot1(n),halfdt)
          call multifab_fill_boundary(s1(n))
        end do

        call average(rho_omegadot1(1),rho_omegadotbar1)
        call react_base(p0_old,s0_old,temp0,rho_omegadotbar1,dx(1,dm),halfdt,p0_1,s0_1,gam1)
        call make_grav_edge(grav_edge,s0_1(:,rho_comp),dx(1,dm),spherical)
        call make_div_coeff(div_coeff_new,s0_1(:,rho_comp),p0_1, &
                            gam1,grav_edge,dx(1,dm),anelastic_cutoff)
        call put_beta_on_edges(div_coeff_new,div_coeff_edge)


!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 4 !!
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        print *,'<<< STEP 4 >>>'

        call advect_base(w0,Sbar(:,1),p0_1,p0_2,s0_1,s0_2,temp0,gam1, &
                         div_coeff_edge,&
                         dx(1,dm),dt,anelastic_cutoff,spherical)
        call make_grav_edge(grav_edge,s0_2(:,rho_comp),dx(1,dm),spherical)
        call make_div_coeff(div_coeff_new,s0_2(:,rho_comp),p0_2, &
                            gam1,grav_edge,dx(1,dm),anelastic_cutoff)
        call put_beta_on_edges(div_coeff_new,div_coeff_edge)

        do n = 1,nlevs
           call scalar_advance (uold(n), s1(n), s2(n), &
                                umac(n,:), w0, sedge(n,:), utrans(n,:),&
                                scal_force(n), s0_1, s0_2, p0_1, p0_2, &
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

!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       !! STEP 5 !!
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        print *,'<<< STEP 5 >>>'
        do n = 1,nlevs
          call react_state(s2(n),snew(n),temp0,rho_omegadot2(n),halfdt)
          call multifab_fill_boundary(s2(n))
        end do
        call average(rho_omegadot2(1),rho_omegadotbar2)
        call react_base(p0_2,s0_2,temp0,rho_omegadotbar2,dx(1,dm),halfdt,p0_new,s0_new,gam1)
        call make_grav_edge(grav_edge,s0_new(:,rho_comp),dx(1,dm),spherical)
        call make_div_coeff(div_coeff_new,s0_new(:,rho_comp),p0_new, &
                            gam1,grav_edge,dx(1,dm),anelastic_cutoff)
        call put_beta_on_edges(div_coeff_old,div_coeff_edge)

!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 6 !!
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        print *,'<<< STEP 6 >>>'
        do n = 1, nlevs
           call make_S(Source_nph(n),snew(n),p0_new,temp0,gam1,dx(n,:),time)
           call make_S(Source_new(n),snew(n),p0_new,temp0,gam1,dx(n,:),time)
           call average(Source_nph(n),Sbar)
           call make_w0(w0,Sbar(:,1),p0_new,s0_new(:,rho_comp),temp0,gam1,dx(n,dm),dt,spherical)
        end do

!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 7 !!
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        print *,'<<< STEP 7 >>>'

        do n = 1,nlevs
           call advance_premac(uold(n), sold(n),&
                               umac(n,:), uedge(n,:), utrans(n,:),&
                               gp(n), p(n), force(n), &
                               s0_old, &
                               dx(n,:),time,dt, &
                               the_bc_tower%bc_tower_array(n), &
                               verbose)
        end do

        ! Define rho at half time !
        do n = 1,nlevs
          call make_rhohalf(rhohalf(n),sold(n),snew(n))
          call multifab_fill_boundary(rhohalf(n))
        end do

        ! Define beta at half time !
        do j = 1,size(div_coeff_old,dim=1)
          div_coeff_nph(j) = HALF * (div_coeff_old(j) + div_coeff_new(j))
        end do
        call put_beta_on_edges(div_coeff_nph,div_coeff_edge)

        do n = 1, nlevs
           call make_macrhs(macrhs(n),Source_nph(n),Sbar(:,1),div_coeff_nph)
        end do

        ! MAC projection !
        call macproject(mla,umac,rhohalf,dx,the_bc_tower,verbose,mg_verbose,cg_verbose,&
                        press_comp,macrhs,div_coeff_nph,div_coeff_edge)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 8 !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        print *,'<<< STEP 8 >>>'
        call advect_base(w0,Sbar(:,1),p0_1,p0_2,s0_1,s0_2,temp0,gam1, &
                         div_coeff_edge,&
                         dx(1,dm),dt,anelastic_cutoff,spherical)
        call make_grav_edge(grav_edge,s0_2(:,rho_comp),dx(1,dm),spherical)
        call make_div_coeff(div_coeff_new,s0_2(:,rho_comp),p0_2, &
                            gam1,grav_edge,dx(1,dm),anelastic_cutoff)
        call put_beta_on_edges(div_coeff_new,div_coeff_edge)

        do n = 1,nlevs
           call scalar_advance (uold(n), s1(n), s2(n), &
                                umac(n,:), w0, sedge(n,:), utrans(n,:),&
                                scal_force(n), s0_1, s0_2, p0_1, p0_2, &
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

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 9 !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        print *,'<<< STEP 9 >>>'
        do n = 1,nlevs
          call react_state(s2(n),snew(n),temp0,rho_omegadot1(n),halfdt)
          call multifab_fill_boundary(snew(n))
        end do
        call average(rho_omegadot2(1),rho_omegadotbar2)
        call react_base(p0_2,s0_2,temp0,rho_omegadotbar2,dx(1,dm),halfdt,p0_new,s0_new,gam1)
        call make_grav_edge(grav_edge,s0_new(:,rho_comp),dx(1,dm),spherical)
        call make_div_coeff(div_coeff_new,s0_new(:,rho_comp),p0_new, &
                            gam1,grav_edge,dx(1,dm),anelastic_cutoff)
        call put_beta_on_edges(div_coeff_new,div_coeff_edge)



        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! STEP 10 !!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        print *,'<<< STEP 10 >>>'
        do n = 1,nlevs
           call velocity_advance(uold(n),unew(n),sold(n),snew(n),rhohalf(n),&
                                 umac(n,:),uedge(n,:), &
                                 utrans(n,:),gp(n),p(n), &
                                 force(n), w0, s0_old, s0_nph, &
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
           call average(Source_new(n),Sbar)
           call make_hgrhs(hgrhs(n),Source_new(n),Sbar(:,1),div_coeff_new)
        end do

!       Project the new velocity field.
        call hgproject(mla, unew, rhohalf, p, gp, dx, dt, &
                       the_bc_tower, verbose, mg_verbose, cg_verbose, press_comp, &
                       hgrhs, div_coeff_nph)

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
        end do
        deallocate(Source_nph)
        deallocate(macrhs)
        deallocate( hgrhs)
        deallocate(rhohalf)

        deallocate(  Sbar)
        deallocate(div_coeff_nph)
        deallocate(div_coeff_edge)
        deallocate(grav_edge)

    end subroutine advance_timestep

end module advance_timestep_module
