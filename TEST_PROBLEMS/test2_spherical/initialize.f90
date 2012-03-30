module initialize_module

  use define_bc_module
  use ml_layout_module
  use multifab_module
  use bc_module
  use probin_module
  use variables, only: nscal, rho_comp, rhoh_comp, temp_comp
  use geometry
  use network, only: nspec
  use bl_constants_module
  use base_state_module
  use base_io_module

  implicit none

  private

  public :: initialize_from_restart, initialize_with_fixed_grids, &
            initialize_with_adaptive_grids, initialize_bc, initialize_dx

contains
    
  subroutine initialize_from_restart(mla,restart,time,dt,pmask,dx,uold,sold,gpres,pres, &
                                     dSdt,Source_old,Source_new, &
                                     rho_omegadot2,rho_Hnuc2,rho_Hext,thermal2,the_bc_tower, &
                                     div_coeff_old,div_coeff_new,gamma1bar,gamma1bar_hold, &
                                     s0_init,rho0_old,rhoh0_old,rho0_new,rhoh0_new,p0_init, &
                                     p0_old,p0_new,w0,etarho_ec,etarho_cc,psi, &
                                     tempbar,grav_cell)

    use restart_module
    use ml_restriction_module
    use multifab_fill_ghost_module
    use multifab_physbc_module
    use probin_module, only : drdxfac

    type(ml_layout),intent(out)   :: mla
    integer       , intent(in   ) :: restart
    real(dp_t)    , intent(  out) :: time,dt
    logical       , intent(in   ) :: pmask(:)
    real(dp_t)    , pointer       :: dx(:,:)
    type(multifab), pointer       :: uold(:),sold(:),gpres(:),pres(:),dSdt(:)
    type(multifab), pointer       :: Source_old(:),Source_new(:)
    type(multifab), pointer       :: rho_omegadot2(:),rho_Hnuc2(:),rho_Hext(:),thermal2(:)
    type(bc_tower), intent(  out) :: the_bc_tower
    real(dp_t)    , pointer       :: div_coeff_old(:,:),div_coeff_new(:,:),gamma1bar(:,:)
    real(dp_t)    , pointer       :: gamma1bar_hold(:,:),s0_init(:,:,:),rho0_old(:,:)
    real(dp_t)    , pointer       :: rhoh0_old(:,:),rho0_new(:,:),rhoh0_new(:,:),p0_init(:,:)
    real(dp_t)    , pointer       :: p0_old(:,:),p0_new(:,:),w0(:,:),etarho_ec(:,:)
    real(dp_t)    , pointer       :: etarho_cc(:,:),psi(:,:),tempbar(:,:),grav_cell(:,:)

    ! local
    type(ml_boxarray) :: mba

    real(dp_t) :: lenx,leny,lenz,max_dist

    integer :: n,ng_s

    type(multifab), pointer :: chkdata(:)
    type(multifab), pointer :: chk_p(:)
    type(multifab), pointer :: chk_dsdt(:)
    type(multifab), pointer :: chk_src_old(:)
    type(multifab), pointer :: chk_src_new(:)
    type(multifab), pointer :: chk_rho_omegadot2(:)
    type(multifab), pointer :: chk_rho_Hnuc2(:)
    type(multifab), pointer :: chk_rho_Hext(:)
    type(multifab), pointer :: chk_thermal2(:)

    character(len=5)   :: check_index
    character(len=6)   :: check_index6
    character(len=256) :: check_file_name

    ! create mba, chk stuff, time, and dt
    call fill_restart_data(restart, mba, chkdata, chk_p, chk_dsdt, chk_src_old, &
                           chk_src_new, chk_rho_omegadot2, chk_rho_Hnuc2, &
                           chk_rho_Hext,chk_thermal2, time, dt)

    ! create mla
    call ml_layout_build(mla,mba,pmask)

    ! check for proper nesting
    if (.not. ml_boxarray_properly_nested(mla%mba, 3, pmask)) then
       call bl_error('restart_grids not properly nested')
    end if

    ! initialize nlevs
    nlevs = mla%nlevel
    nlevs_radial = merge(1, nlevs, spherical .eq. 1)

    ! initialize boundary conditions
    call initialize_bc(the_bc_tower,nlevs,pmask)
    do n = 1,nlevs
       call bc_tower_level_build(the_bc_tower,n,mla%la(n))
    end do

    ! allocate states
    allocate(uold(nlevs),sold(nlevs),gpres(nlevs),pres(nlevs))
    allocate(dSdt(nlevs),Source_old(nlevs),Source_new(nlevs))
    allocate(rho_omegadot2(nlevs),rho_Hnuc2(nlevs),rho_Hext(nlevs),thermal2(nlevs))

    if (ppm_type .eq. 2) then
       ng_s = 4
    else
       ng_s = 3
    end if

    ! build and fill states
    do n = 1,nlevs
       call multifab_build(         uold(n), mla%la(n),    dm, ng_s)
       call multifab_build(         sold(n), mla%la(n), nscal, ng_s)
       call multifab_build(        gpres(n), mla%la(n),    dm, 1)
       call multifab_build(         pres(n), mla%la(n),     1, 1, nodal)
       call multifab_build(         dSdt(n), mla%la(n),     1, 0)
       call multifab_build(   Source_old(n), mla%la(n),     1, 1)
       call multifab_build(   Source_new(n), mla%la(n),     1, 1)
       call multifab_build(rho_omegadot2(n), mla%la(n), nspec, 0)
       call multifab_build(    rho_Hnuc2(n), mla%la(n),     1, 0)
       call multifab_build(     rho_Hext(n), mla%la(n),     1, 0)
       call multifab_build(     thermal2(n), mla%la(n),     1, 0)
    end do

    do n=1,nlevs
       call multifab_copy_c( uold(n),1,chkdata(n),1                ,dm)
       call multifab_copy_c( sold(n),1,chkdata(n),rho_comp+dm      ,nscal)
       call multifab_copy_c(gpres(n),1,chkdata(n),rho_comp+dm+nscal,dm)
       call destroy(chkdata(n)%la)
       call destroy(chkdata(n))
    end do
    
    do n=1,nlevs
       call multifab_copy_c(pres(n),1,chk_p(n),1,1)       
       call destroy(chk_p(n)%la)
       call destroy(chk_p(n))
    end do
    
    do n=1,nlevs
       call multifab_copy_c(dSdt(n),1,chk_dsdt(n),1,1)
       call destroy(chk_dsdt(n)%la)
       call destroy(chk_dsdt(n))
    end do
    
    do n=1,nlevs
       call multifab_copy_c(Source_old(n),1,chk_src_old(n),1,1)
       call destroy(chk_src_old(n)%la)
       call destroy(chk_src_old(n))
    end do

    do n=1,nlevs
       call multifab_copy_c(Source_new(n),1,chk_src_new(n),1,1)
       call destroy(chk_src_new(n)%la)
       call destroy(chk_src_new(n))
    end do
    
    ! Note: rho_omegadot2, rho_Hnuc2, rho_Hext and thermal2 are not actually needed other
    ! than to have them available when we print a plotfile immediately after
    ! restart.  They are recomputed before they are used.

    do n=1,nlevs
       call multifab_copy_c(rho_omegadot2(n),1,chk_rho_omegadot2(n),1,nspec)
       call destroy(chk_rho_omegadot2(n)%la)
       call destroy(chk_rho_omegadot2(n))
    end do

    do n=1,nlevs
       call multifab_copy_c(rho_Hnuc2(n),1,chk_rho_Hnuc2(n),1,1)
       call destroy(chk_rho_Hnuc2(n)%la)
       call destroy(chk_rho_Hnuc2(n))
    end do

    if (plot_Hext) then
       do n=1,nlevs
          call multifab_copy_c(rho_Hext(n),1,chk_rho_Hext(n),1,1)
          call destroy(chk_rho_Hext(n)%la)
          call destroy(chk_rho_Hext(n))
       end do
       deallocate(chk_rho_Hext)
    else
       do n=1,nlevs
          call setval(rho_Hext(n),ZERO,all=.true.)
       end do
    end if

    if (use_thermal_diffusion) then
       do n=1,nlevs
          call multifab_copy_c(thermal2(n),1,chk_thermal2(n),1,0)
          call destroy(chk_thermal2(n)%la)
          call destroy(chk_thermal2(n))
       end do
       deallocate(chk_thermal2)
    else
       do n=1,nlevs
          call setval(thermal2(n),ZERO,all=.true.)
       end do
    end if
    
    deallocate(chkdata, chk_p, chk_dsdt, chk_src_old, chk_src_new)
    deallocate(chk_rho_omegadot2, chk_rho_Hnuc2)

    ! initialize dx
    call initialize_dx(dx,mba,nlevs)

    ! initialize cutoff arrays
    call init_cutoff(nlevs)

    ! now that we have dx we can initialize nr_fine and dr_fine
    if (spherical .eq. 1) then

       ! for spherical, we will now require that dr_fine = dx
       dr_fine = dx(nlevs,1) / dble(drdxfac)
       
       lenx = HALF * (prob_hi(1) - prob_lo(1))
       leny = HALF * (prob_hi(2) - prob_lo(2))
       lenz = HALF * (prob_hi(3) - prob_lo(3))
       
       max_dist = sqrt(lenx**2 + leny**2 + lenz**2)
       nr_fine = int(max_dist / dr_fine) + 1
       
    else
       
       nr_fine = extent(mla%mba%pd(nlevs),dm)
       dr_fine = (prob_hi(dm)-prob_lo(dm)) / dble(nr_fine)
       
    end if

    ! create numdisjointchunks, r_start_coord, r_end_coord
    call init_multilevel(sold)

    ! now that we have nr_fine and dr_fine we can create nr, dr, r_cc_loc, r_edge_loc
    call init_radial(nlevs,mba)

    ! now that we have nr_fine we can allocate 1d arrays
    call initialize_1d_arrays(nlevs,div_coeff_old,div_coeff_new,gamma1bar,gamma1bar_hold, &
                              s0_init,rho0_old,rhoh0_old,rho0_new,rhoh0_new,p0_init, &
                              p0_old,p0_new,w0,etarho_ec,etarho_cc,psi,tempbar, &
                              grav_cell)

    ! now that we have nr we can read in the base state
    if (restart <= 99999) then
       write(unit=check_index,fmt='(i5.5)') restart
       check_file_name = trim(check_base_name) // check_index
    else
       write(unit=check_index6,fmt='(i6.6)') restart
       check_file_name = trim(check_base_name) // check_index6
    endif

    ! note: still need to load/store tempbar
    call read_base_state(restart, check_file_name, &
                         rho0_old, rhoh0_old, p0_old, gamma1bar, w0, &
                         etarho_ec, etarho_cc, div_coeff_old, psi, tempbar)

    ! fill ghost cells
    ! this need to be done after read_base_state since in some problems, the inflow
    ! boundary conditions are set in read_base_state
    if (nlevs .eq. 1) then
        
       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(sold(nlevs))
       call multifab_fill_boundary(uold(nlevs))
       
       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(sold(nlevs),rho_comp,dm+rho_comp,nscal, &
                            the_bc_tower%bc_tower_array(nlevs))
       call multifab_physbc(uold(nlevs),       1,          1,   dm, &
                            the_bc_tower%bc_tower_array(nlevs))
       
    else

       ! the loop over nlevs must count backwards to make sure the finer grids are 
       ! done first
       do n = nlevs,2,-1
          
          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(uold(n-1),uold(n),mla%mba%rr(n-1,:))
          call ml_cc_restriction(sold(n-1),sold(n),mla%mba%rr(n-1,:))
          
          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(uold(n),uold(n-1), &
                                         uold(n)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n  ), &
                                         1,1,dm,fill_crse_input=.false.)
          call multifab_fill_ghost_cells(sold(n),sold(n-1), &
                                         sold(n)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n  ), &
                                         rho_comp,dm+rho_comp,nscal,fill_crse_input=.false.)
       end do
       
    end if

    call destroy(mba)

  end subroutine initialize_from_restart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initialize_with_fixed_grids(mla,time,dt,pmask,dx,uold,sold,gpres,pres, &
                                         dSdt,Source_old,Source_new, &
                                         rho_omegadot2,rho_Hnuc2,rho_Hext,thermal2, &
                                         the_bc_tower,div_coeff_old,div_coeff_new, &
                                         gamma1bar,gamma1bar_hold,s0_init,rho0_old, &
                                         rhoh0_old,rho0_new,rhoh0_new,p0_init, &
                                         p0_old,p0_new,w0,etarho_ec,etarho_cc, &
                                         psi,tempbar,grav_cell)

    use box_util_module
    use init_module
    use average_module
    use restrict_base_module
    use probin_module, only : drdxfac
    
    type(ml_layout),intent(out  ) :: mla
    real(dp_t)    , intent(inout) :: time,dt
    logical       , intent(in   ) :: pmask(:)
    real(dp_t)    , pointer       :: dx(:,:)
    type(multifab), pointer       :: uold(:),sold(:),gpres(:),pres(:),dSdt(:)
    type(multifab), pointer       :: Source_old(:),Source_new(:)
    type(multifab), pointer       :: rho_omegadot2(:),rho_Hnuc2(:),rho_Hext(:),thermal2(:)
    type(bc_tower), intent(  out) :: the_bc_tower
    real(dp_t)    , pointer       :: div_coeff_old(:,:),div_coeff_new(:,:),gamma1bar(:,:)
    real(dp_t)    , pointer       :: gamma1bar_hold(:,:),s0_init(:,:,:),rho0_old(:,:)
    real(dp_t)    , pointer       :: rhoh0_old(:,:),rho0_new(:,:),rhoh0_new(:,:),p0_init(:,:)
    real(dp_t)    , pointer       :: p0_old(:,:),p0_new(:,:),w0(:,:),etarho_ec(:,:)
    real(dp_t)    , pointer       :: etarho_cc(:,:),psi(:,:),tempbar(:,:),grav_cell(:,:)

    ! local
    type(ml_boxarray) :: mba

    real(dp_t) :: lenx,leny,lenz,max_dist

    integer :: n,ng_s
    
    ! set time and dt
    time = ZERO
    dt = 1.d20

    ! create mba
    call read_a_hgproj_grid(mba,test_set)

    ! create mla
    call ml_layout_build(mla,mba,pmask)
    
    ! check for proper nesting
    if (.not. ml_boxarray_properly_nested(mla%mba, 3, pmask)) then
       call bl_error('fixed_grids not properly nested')
    end if
    
    ! initialize nlevs
    nlevs = mla%nlevel
    nlevs_radial = merge(1, nlevs, spherical .eq. 1)

    ! initialize boundary conditions
    call initialize_bc(the_bc_tower,nlevs,pmask)
    do n = 1,nlevs
       call bc_tower_level_build(the_bc_tower,n,mla%la(n))
    end do

    ! allocate states
    allocate(uold(nlevs),sold(nlevs),gpres(nlevs),pres(nlevs))
    allocate(dSdt(nlevs),Source_old(nlevs),Source_new(nlevs))
    allocate(rho_omegadot2(nlevs),rho_Hnuc2(nlevs),rho_Hext(nlevs),thermal2(nlevs))

    if (ppm_type .eq. 2) then
       ng_s = 4
    else
       ng_s = 3
    end if

    ! build states
    do n = 1,nlevs
       call multifab_build(         uold(n), mla%la(n),    dm, ng_s)
       call multifab_build(         sold(n), mla%la(n), nscal, ng_s)
       call multifab_build(        gpres(n), mla%la(n),    dm, 1)
       call multifab_build(         pres(n), mla%la(n),     1, 1, nodal)
       call multifab_build(         dSdt(n), mla%la(n),     1, 0)
       call multifab_build(   Source_old(n), mla%la(n),     1, 1)
       call multifab_build(   Source_new(n), mla%la(n),     1, 1)
       call multifab_build(rho_omegadot2(n), mla%la(n), nspec, 0)
       call multifab_build(    rho_Hnuc2(n), mla%la(n),     1, 0)
       call multifab_build(     rho_Hext(n), mla%la(n),     1, 0)
       call multifab_build(     thermal2(n), mla%la(n),     1, 0)

       call setval(         uold(n), ZERO, all=.true.)
       call setval(         sold(n), ZERO, all=.true.)
       call setval(        gpres(n), ZERO, all=.true.)
       call setval(         pres(n), ZERO, all=.true.)
       call setval(   Source_old(n), ZERO, all=.true.)
       call setval(   Source_new(n), ZERO, all=.true.)
       call setval(         dSdt(n), ZERO, all=.true.)
       call setval(rho_omegadot2(n), ZERO, all=.true.)
       call setval(    rho_Hnuc2(n), ZERO, all=.true.)
       call setval(    rho_Hext(n), ZERO, all=.true.)
       call setval(     thermal2(n), ZERO, all=.true.)
    end do
    ! initialize dx
    call initialize_dx(dx,mba,nlevs)

    ! initialize cutoff arrays
    call init_cutoff(nlevs)

    ! now that we have dx we can initialize nr_fine and dr_fine
    if (spherical .eq. 1) then

       ! for spherical, we will now require that dr_fine = dx
       dr_fine = dx(nlevs,1) / dble(drdxfac)
       
       lenx = HALF * (prob_hi(1) - prob_lo(1))
       leny = HALF * (prob_hi(2) - prob_lo(2))
       lenz = HALF * (prob_hi(3) - prob_lo(3))
       
       max_dist = sqrt(lenx**2 + leny**2 + lenz**2)
       nr_fine = int(max_dist / dr_fine) + 1
       
    else
       
       nr_fine = extent(mla%mba%pd(nlevs),dm)
       dr_fine = (prob_hi(dm)-prob_lo(dm)) / dble(nr_fine)
       
    end if

    ! create numdisjointchunks, r_start_coord, r_end_coord
    call init_multilevel(sold)

    ! now that we have nr_fine and dr_fine we can create nr, dr, r_cc_loc, r_edge_loc
    call init_radial(nlevs,mba)

    ! now that we have nr_fine we can allocate 1d arrays
    call initialize_1d_arrays(nlevs,div_coeff_old,div_coeff_new,gamma1bar,gamma1bar_hold, &
                              s0_init,rho0_old,rhoh0_old,rho0_new,rhoh0_new,p0_init, &
                              p0_old,p0_new,w0,etarho_ec,etarho_cc,psi,tempbar, &
                              grav_cell)

    ! now that we have dr and nr we can fill initial state
    if (spherical .eq. 1) then
       call init_base_state(1,model_file,s0_init(1,:,:),p0_init(1,:),dx(nlevs,:))
    else
       ! init_base_state requires loop backwards over levels
       do n=nlevs,1,-1
          call init_base_state(n,model_file,s0_init(n,:,:),p0_init(n,:),dx(n,:))
       end do
    end if

    call initveldata(uold,s0_init,p0_init,dx,the_bc_tower%bc_tower_array,mla)
    call initscalardata(sold,s0_init,p0_init,dx,the_bc_tower%bc_tower_array,mla)

    if (evolve_base_state) then

       if (perturb_model) then
          ! force rho0 to be the average density
          call average(mla,sold,rho0_old,dx,rho_comp)
       else
          ! s0_init already contains the average
          rho0_old = s0_init(:,:,rho_comp)
          call restrict_base(rho0_old,.true.)
          call fill_ghost_base(rho0_old,.true.)
       end if

       ! this will be overwritten if we are multilevel or if perturb_model = T
       ! but we copy it anyway so s0_init can remain local
       rhoh0_old = s0_init(:,:,rhoh_comp)
       call restrict_base(rhoh0_old,.true.)
       call fill_ghost_base(rhoh0_old,.true.)

       ! this will be overwritten if we are multilevel or if perturb_model = T
       ! but we copy it anyway for the initial condition
       p0_old = p0_init
       call restrict_base(p0_old,.true.)
       call fill_ghost_base(p0_old,.true.)

    else

       if (do_smallscale) then
          ! leave rho0_old = rhoh0_old ZERO
          ! but we still need p0
          p0_old = p0_init
       else
          rho0_old = s0_init(:,:,rho_comp)
          rhoh0_old = s0_init(:,:,rhoh_comp)
          p0_old = p0_init
       end if

    end if

    ! this will be overwritten if we are multilevel or if perturb_model = T
    ! but we copy it anyway for the initial condition
    tempbar = s0_init(:,:,temp_comp)
    call restrict_base(tempbar,.true.)
    call fill_ghost_base(tempbar,.true.)
    
    call destroy(mba)

  end subroutine initialize_with_fixed_grids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initialize_with_adaptive_grids(mla,time,dt,pmask,dx,uold,sold,gpres,pres, &
                                            dSdt,Source_old,Source_new, &
                                            rho_omegadot2,rho_Hnuc2,rho_Hext,thermal2, &
                                            the_bc_tower,div_coeff_old,div_coeff_new, &
                                            gamma1bar,gamma1bar_hold,s0_init,rho0_old, &
                                            rhoh0_old,rho0_new,rhoh0_new,p0_init, &
                                            p0_old,p0_new,w0,etarho_ec,etarho_cc, &
                                            psi,tempbar,grav_cell)

    use probin_module, only: n_cellx, n_celly, n_cellz, regrid_int, max_grid_size, &
         max_grid_size_base, ref_ratio, max_levs
    use init_module
    use average_module
    use restrict_base_module
    use make_new_grids_module
    use probin_module, only : drdxfac
    use multifab_physbc_module
    use ml_restriction_module
    use multifab_fill_ghost_module


    type(ml_layout),intent(out  ) :: mla
    real(dp_t)    , intent(inout) :: time,dt
    logical       , intent(in   ) :: pmask(:)
    real(dp_t)    , pointer       :: dx(:,:)
    type(multifab), pointer       :: uold(:),sold(:),gpres(:),pres(:),dSdt(:)
    type(multifab), pointer       :: Source_old(:),Source_new(:)
    type(multifab), pointer       :: rho_omegadot2(:),rho_Hnuc2(:),rho_Hext(:),thermal2(:)
    type(bc_tower), intent(  out) :: the_bc_tower
    real(dp_t)    , pointer       :: div_coeff_old(:,:),div_coeff_new(:,:),gamma1bar(:,:)
    real(dp_t)    , pointer       :: gamma1bar_hold(:,:),s0_init(:,:,:),rho0_old(:,:)
    real(dp_t)    , pointer       :: rhoh0_old(:,:),rho0_new(:,:),rhoh0_new(:,:),p0_init(:,:)
    real(dp_t)    , pointer       :: p0_old(:,:),p0_new(:,:),w0(:,:),etarho_ec(:,:)
    real(dp_t)    , pointer       :: etarho_cc(:,:),psi(:,:),tempbar(:,:),grav_cell(:,:)

    ! local
    type(ml_boxarray) :: mba

    type(layout) :: la_array(max_levs)
    type(box)    :: bxs

    real(dp_t) :: lenx,leny,lenz,max_dist
    integer    :: buf_wid,n,ng_s,nl
    integer    :: lo(dm), hi(dm)
    logical    :: new_grid

    ! set time and dt
    time = ZERO
    dt = 1.d20

    buf_wid = regrid_int

    ! set up hi & lo to carry indexing info
    lo = 0
    hi(1) = n_cellx-1
    if (dm > 1) then   
       hi(2) = n_celly - 1        
       if (dm > 2)  then
          hi(3) = n_cellz -1
       endif
    endif

    ! mba is big enough to hold max_levs levels
    call ml_boxarray_build_n(mba,max_levs,dm)
    do n = 1, max_levs-1
       mba%rr(n,:) = ref_ratio
    enddo

    ! allocate states
    allocate(uold(max_levs),sold(max_levs),gpres(max_levs),pres(max_levs))
    allocate(dSdt(max_levs),Source_old(max_levs),Source_new(max_levs))
    allocate(rho_omegadot2(max_levs),rho_Hnuc2(max_levs),rho_Hext(max_levs),thermal2(max_levs))

    ! Build the level 1 boxarray
    call box_build_2(bxs,lo,hi)
    call boxarray_build_bx(mba%bas(1),bxs)
    call boxarray_maxsize(mba%bas(1),max_grid_size_base)

    ! build pd(:)
    mba%pd(1) = bxs
    do n = 2, max_levs
       mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
    enddo

    ! initialize dx
    call initialize_dx(dx,mba,max_levs)

    ! initialize cutoff arrays
    call init_cutoff(max_levs)

    if (ppm_type .eq. 2) then
       ng_s = 4
    else
       ng_s = 3
    end if

    ! now that we have dx we can initialize nr_fine and dr_fine
    if (spherical .eq. 1) then

       ! for spherical, we will now require that dr_fine = dx
       dr_fine = dx(max_levs,1) / dble(drdxfac)

       lenx = HALF * (prob_hi(1) - prob_lo(1))
       leny = HALF * (prob_hi(2) - prob_lo(2))
       lenz = HALF * (prob_hi(3) - prob_lo(3))

       max_dist = sqrt(lenx**2 + leny**2 + lenz**2)
       nr_fine = int(max_dist / dr_fine) + 1

    else

       nr_fine = extent(mba%pd(max_levs),dm)
       dr_fine = (prob_hi(dm)-prob_lo(dm)) / dble(nr_fine)

    end if

    ! now that we have nr_fine and dr_fine we can create nr, dr, r_cc_loc, r_edge_loc
    call init_radial(max_levs,mba)

    ! now that we have nr_fine we can allocate 1d arrays
    call initialize_1d_arrays(max_levs,div_coeff_old,div_coeff_new,gamma1bar, &
                              gamma1bar_hold,s0_init,rho0_old,rhoh0_old,rho0_new, &
                              rhoh0_new,p0_init,p0_old,p0_new,w0,etarho_ec,etarho_cc, &
                              psi,tempbar,grav_cell)

    ! now that we have dr and nr we can fill initial state
    if (spherical .eq. 1) then
       call init_base_state(1,model_file,s0_init(1,:,:),p0_init(1,:),dx(max_levs,:))
    else
       ! init_base_state requires loop backwards over levels
       do n=max_levs,1,-1
          call init_base_state(n,model_file,s0_init(n,:,:),p0_init(n,:),dx(n,:))
       end do
    end if

    tempbar = s0_init(:,:,temp_comp)

    ! Initialize bc's
    call initialize_bc(the_bc_tower,max_levs,pmask)

    ! Build the level 1 layout.
    call layout_build_ba(la_array(1),mba%bas(1),mba%pd(1),pmask)

    ! Build the level 1 data only
    call multifab_build(sold(1), la_array(1), nscal, ng_s)

    ! Define bc_tower at level 1.
    call bc_tower_level_build(the_bc_tower,1,la_array(1))

    nlevs = 1
    nlevs_radial = 1

    if (max_levs > 1) then

       ! Initialize the level 1 data only.
       call initscalardata_on_level(1,sold(1),s0_init(1,:,:),p0_init(1,:), &
                                    dx(1,:),the_bc_tower%bc_tower_array(1))

       new_grid = .true.
       nl = 1
       
       do while ( (nl .lt. max_levs) .and. (new_grid) )
          
          ! Do we need finer grids?
          if (nl .eq. 1) then
             call multifab_fill_boundary(sold(1))
             call multifab_physbc(sold(1),rho_comp,dm+rho_comp,nscal, &
                                  the_bc_tower%bc_tower_array(1))
          else
             do n=nl,2,-1
                call ml_cc_restriction(sold(n-1),sold(n),mba%rr(n-1,:))
                call multifab_fill_ghost_cells(sold(n),sold(n-1), &
                                               sold(n)%ng,mba%rr(n-1,:), &
                                               the_bc_tower%bc_tower_array(n-1), &
                                               the_bc_tower%bc_tower_array(n), &
                                               rho_comp,dm+rho_comp,nscal, &
                                               fill_crse_input=.false.)
             enddo
          endif

          call make_new_grids(new_grid,la_array(nl),la_array(nl+1),sold(nl),dx(nl,1), &
                              buf_wid,ref_ratio,nl,max_grid_size,tempbar)
          
          if (new_grid) then
              
             call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))
             
             ! Build the level nl+1 data only.
             call multifab_build(sold(nl+1),la_array(nl+1),nscal,ng_s)
             
             ! Define bc_tower at level nl+1.
             call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))
             
             if (spherical .eq. 1) then
                call initscalardata_on_level(nl+1,sold(nl+1),s0_init(1,:,:), &
                                             p0_init(1,:),dx(nl+1,:), &
                                             the_bc_tower%bc_tower_array(nl+1))
             else
                ! fills the physical region of each level with problem data
                call initscalardata_on_level(nl+1,sold(nl+1),s0_init(nl+1,:,:), &
                                             p0_init(nl+1,:),dx(nl+1,:), &
                                             the_bc_tower%bc_tower_array(nl+1))
             end if

             nlevs = nl+1
             nlevs_radial = merge(1, nlevs, spherical .eq. 1)
             nl = nl + 1
             
          endif ! if (new_grid) 
          
       enddo
       
       do n = 1,nlevs
          call destroy(sold(n))
       end do

       nlevs = nl
       nlevs_radial = merge(1, nlevs, spherical .eq. 1)

       ! check for proper nesting
       if (nlevs .ge. 3) then

          call enforce_proper_nesting(mba,la_array,max_grid_size)

          ! enforce_proper_nesting can create new grids at coarser levels
          ! this makes sure the boundary conditions are properly defined everywhere
          do n=2,nlevs
             call bc_tower_level_build(the_bc_tower,n,la_array(n))
          end do

       end if
       
    else

       call destroy(sold(1))
       
    end if  ! end if (maxlev > 1)

    call ml_layout_restricted_build(mla,mba,nlevs,pmask)
    
    nlevs = mla%nlevel

    if (nlevs .ne. max_levs) then
       call bl_error('initialize_with_adaptive_grids: nlevs .ne. max_levs not supported yet')
    end if

    nlevs_radial = merge(1, nlevs, spherical .eq. 1)
    
    do n = 1, nlevs
       call destroy(la_array(n))
    end do
    
    ! build states
    do n = 1,nlevs
       call multifab_build(         uold(n), mla%la(n),    dm, ng_s)
       call multifab_build(         sold(n), mla%la(n), nscal, ng_s)
       call multifab_build(        gpres(n), mla%la(n),    dm, 1)
       call multifab_build(         pres(n), mla%la(n),     1, 1, nodal)
       call multifab_build(         dSdt(n), mla%la(n),     1, 0)
       call multifab_build(   Source_old(n), mla%la(n),     1, 1)
       call multifab_build(   Source_new(n), mla%la(n),     1, 1)
       call multifab_build(rho_omegadot2(n), mla%la(n), nspec, 0)
       call multifab_build(    rho_Hnuc2(n), mla%la(n),     1, 0)
       call multifab_build(     rho_Hext(n), mla%la(n),     1, 0)
       call multifab_build(     thermal2(n), mla%la(n),     1, 0)

       call setval(         uold(n), ZERO, all=.true.)
       call setval(         sold(n), ZERO, all=.true.)
       call setval(        gpres(n), ZERO, all=.true.)
       call setval(         pres(n), ZERO, all=.true.)
       call setval(   Source_old(n), ZERO, all=.true.)
       call setval(   Source_new(n), ZERO, all=.true.)
       call setval(         dSdt(n), ZERO, all=.true.)
       call setval(rho_omegadot2(n), ZERO, all=.true.)
       call setval(    rho_Hnuc2(n), ZERO, all=.true.)
       call setval(    rho_Hext(n), ZERO, all=.true.)
       call setval(     thermal2(n), ZERO, all=.true.)
    end do

    ! create numdisjointchunks, r_start_coord, r_end_coord
    call init_multilevel(sold)

    call initveldata(uold,s0_init,p0_init,dx,the_bc_tower%bc_tower_array,mla)
    call initscalardata(sold,s0_init,p0_init,dx,the_bc_tower%bc_tower_array,mla)

    if (evolve_base_state) then

       if (perturb_model) then
          ! force rho0 to be the average density
          call average(mla,sold,rho0_old,dx,rho_comp)
       else
          ! s0_init already contains the average
          rho0_old = s0_init(:,:,rho_comp)
          call restrict_base(rho0_old,.true.)
          call fill_ghost_base(rho0_old,.true.)
       end if

       ! this will be overwritten if we are multilevel or if perturb_model = T
       ! but we copy it anyway so s0_init can remain local
       rhoh0_old = s0_init(:,:,rhoh_comp)
       call restrict_base(rhoh0_old,.true.)
       call fill_ghost_base(rhoh0_old,.true.)

       ! this will be overwritten if we are multilevel or if perturb_model = T
       ! but we copy it anyway for the initial condition
       p0_old = p0_init
       call restrict_base(p0_old,.true.)
       call fill_ghost_base(p0_old,.true.)

    else

       if (do_smallscale) then
          ! leave rho0_old = rhoh0_old ZERO
          ! but we still need p0
          p0_old = p0_init
       else
          rho0_old = s0_init(:,:,rho_comp)
          rhoh0_old = s0_init(:,:,rhoh_comp)
          p0_old = p0_init
       end if

    end if

    ! this will be overwritten if we are multilevel or if perturb_model = T
    ! but we copy it anyway for the initial condition
    tempbar = s0_init(:,:,temp_comp)
    call restrict_base(tempbar,.true.)
    call fill_ghost_base(tempbar,.true.)

    call destroy(mba)

  end subroutine initialize_with_adaptive_grids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine initialize_bc(the_bc_tower,num_levs,pmask)

    use bc_module
    use probin_module, only : bcx_lo, bcx_hi, bcy_lo, bcy_hi, bcz_lo, bcz_hi

    type(bc_tower), intent(  out) :: the_bc_tower
    integer       , intent(in   ) :: num_levs
    logical       , intent(in   ) :: pmask(:)
    
    integer :: domain_phys_bc(dm,2)

    ! Define the physical boundary conditions on the domain
    ! Put the bc values from the inputs file into domain_phys_bc
    domain_phys_bc(1,1) = bcx_lo
    domain_phys_bc(1,2) = bcx_hi
    if (pmask(1)) then
       domain_phys_bc(1,:) = BC_PER
       if (bcx_lo .ne. -1 .or. bcx_hi .ne. -1) &
            call bl_error('MUST HAVE BCX = -1 if PMASK = T')
    end if
    if (dm > 1) then
       domain_phys_bc(2,1) = bcy_lo
       domain_phys_bc(2,2) = bcy_hi
       if (pmask(2)) then
          domain_phys_bc(2,:) = BC_PER
          if (bcy_lo .ne. -1 .or. bcy_hi .ne. -1) &
               call bl_error('MUST HAVE BCY = -1 if PMASK = T') 
       end if
    end if
    if (dm > 2) then
       domain_phys_bc(3,1) = bcz_lo
       domain_phys_bc(3,2) = bcz_hi
       if (pmask(3)) then
          domain_phys_bc(3,:) = BC_PER
          if (bcz_lo .ne. -1 .or. bcz_hi .ne. -1) &
               call bl_error('MUST HAVE BCZ = -1 if PMASK = T')
       end if
    end if
    
    ! Initialize the_bc_tower object.
    call bc_tower_init(the_bc_tower,num_levs,dm,domain_phys_bc)
    
  end subroutine initialize_bc

  subroutine initialize_dx(dx,mba,num_levs)

    real(dp_t)       , pointer     :: dx(:,:)
    type(ml_boxarray), intent(in ) :: mba
    integer          , intent(in ) :: num_levs
    
    integer :: n,d
    
    allocate(dx(num_levs,dm))
    
    do d=1,dm
       dx(1,d) = (prob_hi(d)-prob_lo(d)) / real(extent(mba%pd(1),d),kind=dp_t)
    end do
    do n=2,num_levs
       dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
    end do

  end subroutine initialize_dx

  subroutine initialize_1d_arrays(num_levs,div_coeff_old,div_coeff_new,gamma1bar, &
                                  gamma1bar_hold,s0_init,rho0_old,rhoh0_old,rho0_new, &
                                  rhoh0_new,p0_init,p0_old,p0_new,w0,etarho_ec,etarho_cc, &
                                  psi,tempbar,grav_cell)

    integer    , intent(in) :: num_levs    
    real(dp_t) , pointer    :: div_coeff_old(:,:),div_coeff_new(:,:),gamma1bar(:,:)
    real(dp_t) , pointer    :: gamma1bar_hold(:,:),s0_init(:,:,:),rho0_old(:,:)
    real(dp_t) , pointer    :: rhoh0_old(:,:),rho0_new(:,:),rhoh0_new(:,:),p0_init(:,:)
    real(dp_t) , pointer    :: p0_old(:,:),p0_new(:,:),w0(:,:),etarho_ec(:,:),etarho_cc(:,:)
    real(dp_t) , pointer    :: psi(:,:),tempbar(:,:),grav_cell(:,:)
    
    if (spherical .eq. 0) then
       allocate(div_coeff_old (num_levs,0:nr_fine-1))
       allocate(div_coeff_new (num_levs,0:nr_fine-1))
       allocate(gamma1bar     (num_levs,0:nr_fine-1))
       allocate(gamma1bar_hold(num_levs,0:nr_fine-1))
       allocate(s0_init       (num_levs,0:nr_fine-1,nscal))
       allocate(rho0_old      (num_levs,0:nr_fine-1))
       allocate(rhoh0_old     (num_levs,0:nr_fine-1))
       allocate(rho0_new      (num_levs,0:nr_fine-1))
       allocate(rhoh0_new     (num_levs,0:nr_fine-1))
       allocate(p0_init       (num_levs,0:nr_fine-1))
       allocate(p0_old        (num_levs,0:nr_fine-1))
       allocate(p0_new        (num_levs,0:nr_fine-1))
       allocate(w0            (num_levs,0:nr_fine))
       allocate(etarho_ec     (num_levs,0:nr_fine))
       allocate(etarho_cc     (num_levs,0:nr_fine-1))
       allocate(psi           (num_levs,0:nr_fine-1))
       allocate(tempbar       (num_levs,0:nr_fine-1))
       allocate(grav_cell     (num_levs,0:nr_fine-1))
    else
       allocate(div_coeff_old (1,0:nr_fine-1))
       allocate(div_coeff_new (1,0:nr_fine-1))
       allocate(gamma1bar     (1,0:nr_fine-1))
       allocate(gamma1bar_hold(1,0:nr_fine-1))
       allocate(s0_init       (1,0:nr_fine-1,nscal))
       allocate(rho0_old      (1,0:nr_fine-1))
       allocate(rhoh0_old     (1,0:nr_fine-1))
       allocate(rho0_new      (1,0:nr_fine-1))
       allocate(rhoh0_new     (1,0:nr_fine-1))
       allocate(p0_init       (1,0:nr_fine-1))
       allocate(p0_old        (1,0:nr_fine-1))
       allocate(p0_new        (1,0:nr_fine-1))
       allocate(w0            (1,0:nr_fine))
       allocate(etarho_ec     (1,0:nr_fine))
       allocate(etarho_cc     (1,0:nr_fine-1))
       allocate(psi           (1,0:nr_fine-1))
       allocate(tempbar       (1,0:nr_fine-1))
       allocate(grav_cell     (1,0:nr_fine-1))
    end if

    div_coeff_old = ZERO
    div_coeff_new = ZERO
    gamma1bar = ZERO
    gamma1bar_hold = ZERO
    s0_init = ZERO
    rho0_old = ZERO
    rhoh0_old = ZERO
    rho0_new = ZERO
    rhoh0_new = ZERO
    p0_init = ZERO
    p0_old = ZERO
    p0_new = ZERO
    w0 = ZERO
    etarho_ec = ZERO
    etarho_cc = ZERO
    psi = ZERO
    tempbar = ZERO
    grav_cell = ZERO

  end subroutine initialize_1d_arrays
  
end module initialize_module
