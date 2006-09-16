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
  use base_state_module
  use variables
  use network

  contains

    subroutine advance_timestep(mla,uold,sold,shalf,unew,snew,umac,uedge,sedge,utrans,gp,p,force,scal_force,&
                                s0_old,s0_new,s0_nph,p0_old,p0_new,temp0,gam1,w0, &
                                div_coeff_n,div_coeff_nph,div_coeff_half,div_coef_type, &
                                dx,time,dt,the_bc_tower, &
                                visc_coef,diff_coef,anelastic_cutoff,grav,verbose,mg_verbose,cg_verbose)


    type(ml_layout),intent(inout) :: mla
    type(multifab), intent(inout) :: uold(:)
    type(multifab), intent(inout) :: sold(:)
    type(multifab), intent(inout) :: shalf(:)
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
    real(dp_t)    , intent(inout) :: s0_old(:,:)
    real(dp_t)    , intent(inout) :: s0_new(:,:)
    real(dp_t)    , intent(inout) :: s0_nph(:,:)
    real(dp_t)    , intent(inout) :: p0_old(:)
    real(dp_t)    , intent(inout) :: p0_new(:)
    real(dp_t)    , intent(inout) :: temp0(:)
    real(dp_t)    , intent(inout) :: gam1(:)
    real(dp_t)    , intent(inout) :: w0(:)
    real(dp_t)    , intent(inout) :: div_coeff_n(:)
    real(dp_t)    , intent(inout) :: div_coeff_nph(:)
    real(dp_t)    , intent(inout) :: div_coeff_half(:)
    real(dp_t)    , intent(in   ) :: dx(:,:), time, dt
    type(bc_tower), intent(in   ) :: the_bc_tower
    real(dp_t)    , intent(in   ) :: visc_coef,diff_coef,grav
    real(dp_t)    , intent(in   ) :: anelastic_cutoff
    integer       , intent(in   ) :: div_coef_type
    integer       , intent(in   ) :: verbose,mg_verbose,cg_verbose

    type(multifab), allocatable :: macrhs(:)
    type(multifab), allocatable ::  hgrhs(:)
    type(bc_level) ::  bc
    type(box)      ::  fine_domain
    real(dp_t)     :: half_time, new_time
    real(dp_t)     :: visc_mu
    integer :: n,dm,nlevs,comp,bc_comp
    integer :: pred_vs_corr
    logical :: nodal(mla%dim)

    nlevs = size(uold)
    dm    = mla%dim

    half_time = dt + HALF*dt
     new_time = dt +      dt

    allocate(macrhs(nlevs))
    allocate( hgrhs(nlevs))

    nodal = .true.
    do n = 1,nlevs
     call multifab_build(    macrhs(n), mla%la(n),     1, 0)
     call multifab_build(  hgrhs(n), mla%la(n),     1,       0, nodal)
    end do

        do n = 1,nlevs
           call advance_premac(uold(n), sold(n),&
                               umac(n,:), uedge(n,:), utrans(n,:),&
                               gp(n), p(n), force(n), &
                               s0_old, temp0, &
                               dx(n,:),time,dt, &
                               the_bc_tower%bc_tower_array(n), &
                               visc_coef,verbose)
        end do

        do n = 1, nlevs
           call make_macrhs(macrhs(n),sold(n),uold(n),div_coeff_n,p0_old,temp0,gam1,dx(n,:),half_time)
        end do

!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       MAC projection 
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call macproject(mla,umac,sold,dx,the_bc_tower,verbose,mg_verbose,cg_verbose,press_comp,&
                        macrhs,div_coeff_n,div_coeff_half)

!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Update rho, (rho h), and (rho X)_i
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        w0 = ZERO
        pred_vs_corr = 1
        do n = 1,nlevs
           call scalar_advance (uold(n), sold(n), snew(n), shalf(n), &
                                umac(n,:), w0, sedge(n,:), utrans(n,:),&
                                scal_force(n), s0_old, s0_old, s0_old, &
                                p0_old, p0_old, temp0, &
                                dx(n,:),time,dt, &
                                the_bc_tower%bc_tower_array(n), &
                                diff_coef,verbose,div_coeff_n,div_coeff_half, &
                                pred_vs_corr)
        end do

        do n = 2, nlevs
           fine_domain = layout_get_pd(mla%la(n))
           call multifab_fill_ghost_cells(snew(n),snew(n-1),fine_domain, &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,rho_comp,nscal)
        end do

        if (diff_coef > ZERO) then
          do comp = rhoh_comp, spec_comp+nspec-1
            bc_comp = dm+comp
            visc_mu = HALF*dt*diff_coef
            call diff_scalar_solve(mla,snew,dx,visc_mu,the_bc_tower,comp,bc_comp,&
                                   mg_verbose,cg_verbose)
          end do
        end if

        do n = 1,nlevs
           call multifab_fill_boundary(snew(n))
        end do

!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Update the base state
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call eval_base_state(w0,p0_old,p0_new, &
                             s0_old,s0_nph,s0_new,temp0, &
                             gam1,div_coeff_n,div_coeff_nph,div_coeff_half, &
                             shalf(1),grav,dx(1,:), &
                             dt,time,div_coef_type,anelastic_cutoff)

!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Update rho and (rho h) again.
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        pred_vs_corr = 2
        do n = 1,nlevs
           call scalar_advance (uold(n), sold(n), snew(n), shalf(n), &
                                umac(n,:), w0, sedge(n,:), utrans(n,:),&
                                scal_force(n), s0_old, s0_new, s0_nph, &
                                p0_old, p0_new, temp0, &
                                dx(n,:),time,dt, &
                                the_bc_tower%bc_tower_array(n), &
                                diff_coef,verbose,div_coeff_n,div_coeff_half, &
                                pred_vs_corr)
        end do

        do n = 2, nlevs
           fine_domain = layout_get_pd(mla%la(n))
           call multifab_fill_ghost_cells(snew(n),snew(n-1),fine_domain, &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,rho_comp,nscal)
        end do

        if (diff_coef > ZERO) then
          comp = 2
          bc_comp = dm+comp
          visc_mu = HALF*dt*diff_coef
          call diff_scalar_solve(mla,snew,dx,visc_mu,the_bc_tower,comp,bc_comp,&
                                 mg_verbose,cg_verbose)
        end if

        do n = 1,nlevs
           call multifab_fill_boundary(snew(n))
        end do

!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Update velocity
!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do n = 1,nlevs
           call velocity_advance(uold(n),unew(n),sold(n),snew(n),shalf(n),&
                                 umac(n,:),uedge(n,:), &
                                 utrans(n,:),gp(n),p(n), &
                                 force(n), w0, s0_old, s0_nph, &
                                 dx(n,:),time,dt, &
                                 the_bc_tower%bc_tower_array(n), &
                                 visc_coef,verbose)
        end do

        do n = 2, nlevs
           fine_domain = layout_get_pd(mla%la(n))
           call multifab_fill_ghost_cells(unew(n),unew(n-1),fine_domain, &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,1,dm)
        end do

        if (visc_coef > ZERO) then
           visc_mu = HALF*dt*visc_coef
           call visc_solve(mla,unew,shalf,dx,visc_mu,the_bc_tower,mg_verbose,cg_verbose)
        end if

        do n = 1,nlevs
           call multifab_fill_boundary(unew(n))
        end do

        do n = 1,nlevs
           call make_hgrhs(hgrhs(n),macrhs(n),snew(n),unew(n),div_coeff_nph,p0_new,temp0,gam1,dx(n,:),new_time)
        end do

!       Project the new velocity field.
        call hgproject(mla, unew, shalf, p, gp, dx, dt, &
                       the_bc_tower, verbose, mg_verbose, cg_verbose, press_comp, hgrhs, div_coeff_nph)

        do n = 2, nlevs
           fine_domain = layout_get_pd(mla%la(n))
           call multifab_fill_ghost_cells(unew(n),unew(n-1),fine_domain, &
                                          ng_cell,mla%mba%rr(n-1,:), &
                                          the_bc_tower%bc_tower_array(n-1)%adv_bc_level_array(0,:,:,:), &
                                          1,1,dm)
        end do

        do n = 1, nlevs
          call destroy(macrhs(n))
          call destroy( hgrhs(n))
        end do
        deallocate(macrhs)
        deallocate( hgrhs)

    end subroutine advance_timestep

end module advance_timestep_module
