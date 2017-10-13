module initialize_module

  use define_bc_module
  use ml_layout_module
  use multifab_module
  use bc_module
  use variables, only: nscal, rho_comp, rhoh_comp, temp_comp, foextrap_comp
  use geometry, only: spherical, nr_irreg, nr_fine, dr_fine, nlevs_radial, &
       init_cutoff, init_multilevel, init_radial, destroy_geometry, &
       compute_cutoff_coords, initialize_dx
  use network, only: nspec
  use bl_constants_module
  use base_state_module

  implicit none

  private

  public :: initialize_from_restart, initialize_with_fixed_grids, &
            initialize_with_adaptive_grids

contains
    
  subroutine initialize_from_restart(mla,restart,dt,pmask,dx,uold,sold,gpi,pi, &
                                     dSdt,S_cc_old,S_cc_new, &
                                     rho_omegadot2,rho_Hnuc2,rho_Hext,thermal2,the_bc_tower, &
                                     div_coeff_old,div_coeff_new,gamma1bar,gamma1bar_hold, &
                                     s0_init,rho0_old,rhoh0_old,rho0_new,rhoh0_new,p0_init, &
                                     p0_old,p0_new,w0,etarho_ec,etarho_cc,psi, &
                                     tempbar,tempbar_init,grav_cell)

    use restart_module
    use multifab_fill_ghost_module
    use probin_module, only : drdxfac, restart_into_finer, octant, max_levs, &
         ppm_type, bds_type, plot_Hext, use_thermal_diffusion, &
         prob_lo, prob_hi, nodal, &
         check_base_name, use_tfromp, cflfac, dm_in, restart_with_vel_field, &
         model_file, do_smallscale, fix_base_state, max_grid_size_1, &
         change_max_grid_size_1, drive_initial_convection, use_tpert_in_tagging
    use average_module
    use make_grav_module
    use enforce_HSE_module
    use rhoh_vs_t_module
    use make_gamma_module
    use make_div_coeff_module
    use fill_3d_module
    use estdt_module
    use regrid_module
    use init_scalar_module
    use time_module, only: time
    use base_io_module
    use aux_data_module
    use ml_restrict_fill_module
    use pert_form_module, only: put_in_pert_form

    type(ml_layout),intent(out)   :: mla
    integer       , intent(inout) :: restart
    real(dp_t)    , intent(  out) :: dt
    logical       , intent(in   ) :: pmask(:)
    real(dp_t)    , pointer       :: dx(:,:)
    type(multifab), pointer       :: uold(:),sold(:),gpi(:),pi(:),dSdt(:)
    type(multifab), pointer       :: S_cc_old(:),S_cc_new(:)
    type(multifab), pointer       :: rho_omegadot2(:),rho_Hnuc2(:),rho_Hext(:),thermal2(:)
    type(bc_tower), intent(  out) :: the_bc_tower
    real(dp_t)    , pointer       :: div_coeff_old(:,:),div_coeff_new(:,:),gamma1bar(:,:)
    real(dp_t)    , pointer       :: gamma1bar_hold(:,:),s0_init(:,:,:),rho0_old(:,:)
    real(dp_t)    , pointer       :: rhoh0_old(:,:),rho0_new(:,:),rhoh0_new(:,:),p0_init(:,:)
    real(dp_t)    , pointer       :: p0_old(:,:),p0_new(:,:),w0(:,:),etarho_ec(:,:)
    real(dp_t)    , pointer       :: etarho_cc(:,:),psi(:,:),tempbar(:,:),tempbar_init(:,:),grav_cell(:,:)

    ! local
    type(ml_boxarray) :: mba, mba_old
    type(box)         :: domain
    integer           :: domhi(dm_in)

    real(dp_t) :: lenx,leny,lenz,max_dist,p0_temp

    integer :: n,ng_s,nr_fine_old,r,dm,nlevs

    type(multifab), pointer :: chkdata(:)
    type(multifab), pointer :: chk_p(:)
    type(multifab), pointer :: chk_dsdt(:)
    type(multifab), pointer :: chk_src_old(:)
    type(multifab), pointer :: chk_src_new(:)
    type(multifab), pointer :: chk_rho_omegadot2(:)
    type(multifab), pointer :: chk_rho_Hnuc2(:)
    type(multifab), pointer :: chk_rho_Hext(:)
    type(multifab), pointer :: chk_thermal2(:)

    type(multifab), allocatable :: gamma1(:)
    type(multifab), pointer :: tpert_mf(:)

    type(layout) :: la

    type(box)    :: bxs

    real(dp_t), allocatable :: psi_temp(:,:)
    real(dp_t), allocatable :: etarho_cc_temp(:,:)
    real(dp_t), allocatable :: w0_temp(:,:)
  
    character(len=5)   :: check_index
    character(len=6)   :: check_index6
    character(len=7)   :: check_index7
    character(len=256) :: check_file_name

    ! sanity check
    if (drive_initial_convection .and. restart_with_vel_field) then
       call bl_error("restart_with_vel_field and drive_initial_convection both T")
    endif

    ! create mba, chk stuff, time, and dt
    call fill_restart_data(restart, mba_old, chkdata, chk_p, chk_dsdt, chk_src_old, &
                           chk_src_new, chk_rho_omegadot2, chk_rho_Hnuc2, &
                           chk_rho_Hext,chk_thermal2, dt)

    if (change_max_grid_size_1) then
       ! Change max grid size if inputs file calls for it
       ! rebuild level 1 box array
       ! for now just copy higher level boxarrays
       ! regrid will take care of the higher levels when it's called
       call ml_boxarray_copy(mba, mba_old)
       call box_build_2(bxs,lwb(mba_old%pd(1)),upb(mba_old%pd(1)))
       call destroy(mba%bas(1))
       call boxarray_build_bx(mba%bas(1),bxs)
       call boxarray_maxsize(mba%bas(1),max_grid_size_1)
    else
       call ml_boxarray_copy(mba, mba_old)
    endif

    ! create mla
    call ml_layout_build(mla,mba,pmask)

    dm = mla%dim

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
    allocate(uold(nlevs),sold(nlevs),gpi(nlevs),pi(nlevs))
    allocate(dSdt(nlevs),S_cc_old(nlevs),S_cc_new(nlevs))
    allocate(rho_omegadot2(nlevs),rho_Hnuc2(nlevs),rho_Hext(nlevs))
    allocate(thermal2(nlevs))

    if (ppm_type .eq. 2 .or. bds_type .eq. 1) then
       ng_s = 4
    else
       ng_s = 3
    end if

    ! build and fill states
    do n = 1,nlevs
       call multifab_build(         uold(n), mla%la(n),    dm, ng_s)
       call multifab_build(         sold(n), mla%la(n), nscal, ng_s)
       call multifab_build(          gpi(n), mla%la(n),    dm, 0)
       call multifab_build(           pi(n), mla%la(n),     1, 0, nodal)
       call multifab_build(         dSdt(n), mla%la(n),     1, 0)
       call multifab_build(     S_cc_old(n), mla%la(n),     1, 0)
       call multifab_build(     S_cc_new(n), mla%la(n),     1, 0)
       call multifab_build(rho_omegadot2(n), mla%la(n), nspec, 0)
       call multifab_build(    rho_Hnuc2(n), mla%la(n),     1, 0)
       call multifab_build(     rho_Hext(n), mla%la(n),     1, 0)
       call multifab_build(     thermal2(n), mla%la(n),     1, 0)
    end do

    do n=1,nlevs
       call multifab_copy_c( uold(n),1,chkdata(n),1                ,dm)
       call multifab_copy_c( sold(n),1,chkdata(n),rho_comp+dm      ,nscal)
       call multifab_copy_c(  gpi(n),1,chkdata(n),rho_comp+dm+nscal,dm)
       la = get_layout(chkdata(n))
       call destroy(chkdata(n))
       call destroy(la)
    end do
    
    do n=1,nlevs
       call multifab_copy_c(pi(n),1,chk_p(n),1,1)
       la = get_layout(chk_p(n))
       call destroy(chk_p(n))
       call destroy(la)
    end do
    
    do n=1,nlevs
       call multifab_copy_c(dSdt(n),1,chk_dsdt(n),1,1)
       la = get_layout(chk_dsdt(n))
       call destroy(chk_dsdt(n))
       call destroy(la)
    end do
    
    do n=1,nlevs
       call multifab_copy_c(S_cc_old(n),1,chk_src_old(n),1,1)
       la = get_layout(chk_src_old(n))
       call destroy(chk_src_old(n))
       call destroy(la)
    end do

    do n=1,nlevs
       call multifab_copy_c(S_cc_new(n),1,chk_src_new(n),1,1)
       la = get_layout(chk_src_new(n))
       call destroy(chk_src_new(n))
       call destroy(la)
    end do
    
    ! Note: rho_omegadot2, rho_Hnuc2, rho_Hext, and thermal2 are not
    ! actually needed other than to have them available when we print
    ! a plotfile immediately after restart.  They are recomputed
    ! before they are used.

    do n=1,nlevs
       call multifab_copy_c(rho_omegadot2(n),1,chk_rho_omegadot2(n),1,nspec)
       la = get_layout(chk_rho_omegadot2(n))
       call destroy(chk_rho_omegadot2(n))
       call destroy(la)
    end do

    do n=1,nlevs
       call multifab_copy_c(rho_Hnuc2(n),1,chk_rho_Hnuc2(n),1,1)
       la = get_layout(chk_rho_Hnuc2(n))
       call destroy(chk_rho_Hnuc2(n))
       call destroy(la)
    end do

    if (plot_Hext) then
       do n=1,nlevs
          call multifab_copy_c(rho_Hext(n),1,chk_rho_Hext(n),1,1)
          la = get_layout(chk_rho_Hext(n))
          call destroy(chk_rho_Hext(n))
          call destroy(la)
       end do
       deallocate(chk_rho_Hext)
    end if

    if (use_thermal_diffusion) then
       do n=1,nlevs
          call multifab_copy_c(thermal2(n),1,chk_thermal2(n),1,1)
          la = get_layout(chk_thermal2(n))
          call destroy(chk_thermal2(n))
          call destroy(la)
       end do
       deallocate(chk_thermal2)
    else
       do n=1,nlevs
          call setval(thermal2(n),ZERO,all=.true.)
       end do
    end if

    ! reset the state data if restarting with vel information
    if (restart_with_vel_field) then
       do n = 1, nlevs
          call setval(         sold(n), ZERO, all=.true.)
          call setval(          gpi(n), ZERO, all=.true.)
          call setval(           pi(n), ZERO, all=.true.)
          call setval(     S_cc_old(n), ZERO, all=.true.)
          call setval(     S_cc_new(n), ZERO, all=.true.)
          call setval(         dSdt(n), ZERO, all=.true.)
          call setval(rho_omegadot2(n), ZERO, all=.true.)
          call setval(    rho_Hnuc2(n), ZERO, all=.true.)
          call setval(     rho_Hext(n), ZERO, all=.true.)
          call setval(     thermal2(n), ZERO, all=.true.)
       end do

     endif
    
    deallocate(chkdata, chk_p, chk_dsdt, chk_src_old, chk_src_new)
    deallocate(chk_rho_omegadot2, chk_rho_Hnuc2)

    ! initialize dx
    call initialize_dx(dx,mla%mba,nlevs)

    ! initialize cutoff arrays
    call init_cutoff(nlevs)

    ! now that we have dx we can initialize nr_fine and dr_fine
    if (spherical .eq. 1) then

       dr_fine = dx(nlevs,1) / dble(drdxfac)
       
       if (.not. octant) then
          lenx = HALF * (prob_hi(1) - prob_lo(1))
          leny = HALF * (prob_hi(2) - prob_lo(2))
          lenz = HALF * (prob_hi(3) - prob_lo(3))
       else
          lenx = prob_hi(1) - prob_lo(1)
          leny = prob_hi(2) - prob_lo(2)
          lenz = prob_hi(3) - prob_lo(3)
       end if
       
       max_dist = sqrt(lenx**2 + leny**2 + lenz**2)
       nr_fine = int(max_dist / dr_fine) + 1

       ! compute nr_irreg
       domain = get_pd(get_layout(sold(nlevs)))
       domhi  = upb(domain)+1
       if (.not. octant) then
          nr_irreg = (3*(domhi(1)/2-0.5d0)**2-0.75d0)/2.d0
       else
          nr_irreg = (3*(domhi(1)-0.5d0)**2-0.75d0)/2.d0
       endif
    else
       
       nr_fine = extent(mla%mba%pd(nlevs),dm)
       dr_fine = (prob_hi(dm)-prob_lo(dm)) / dble(nr_fine)
       
    end if

    ! create numdisjointchunks, r_start_coord, r_end_coord
    call init_multilevel(sold)

    ! now that we have nr_fine and dr_fine we can create nr, dr, r_cc_loc, r_edge_loc
    call init_radial(nlevs,mla%mba)

    ! now that we have nr_fine we can allocate 1d arrays
    call initialize_1d_arrays(nlevs,div_coeff_old,div_coeff_new,gamma1bar,gamma1bar_hold, &
                              s0_init,rho0_old,rhoh0_old,rho0_new,rhoh0_new,p0_init, &
                              p0_old,p0_new,w0,etarho_ec,etarho_cc,psi,tempbar,tempbar_init, &
                              grav_cell)

    if (restart_with_vel_field) then

       if (spherical .eq. 1) then
          call init_base_state(1,model_file,s0_init(1,:,:),p0_init(1,:),dx(max_levs,:))
       else
          ! init_base_state requires loop backwards over levels
          do n=max_levs,1,-1
             call init_base_state(n,model_file,s0_init(n,:,:),p0_init(n,:),dx(n,:))
          end do
       end if

       p0_old       = p0_init
       rho0_old     = s0_init(:,:,rho_comp)
       rhoh0_old    = s0_init(:,:,rhoh_comp)
       tempbar      = s0_init(:,:,temp_comp)
       tempbar_init = s0_init(:,:,temp_comp)

       call initscalardata(sold,s0_init,p0_init,dx,the_bc_tower%bc_tower_array,mla)
       
       if (fix_base_state) then
          call compute_cutoff_coords(rho0_old)
          call make_grav_cell(grav_cell,rho0_old)
          call destroy(mba)
          return
       end if

       if (do_smallscale) then
          ! first compute cutoff coordinates using initial density profile
          call compute_cutoff_coords(rho0_old)
          ! set rho0_old = rhoh0_old = ZERO
          rho0_old  = ZERO
          rhoh0_old = ZERO
       else
          ! set rho0 to be the average
          call average(mla,sold,rho0_old,dx,rho_comp)
          call compute_cutoff_coords(rho0_old)

          ! compute p0 with HSE
          call make_grav_cell(grav_cell,rho0_old)
          call enforce_HSE(rho0_old,p0_old,grav_cell)

          ! call eos with r,p as input to recompute T,h
          call makeTHfromRhoP(sold,p0_old,the_bc_tower%bc_tower_array,mla,dx)

          ! set rhoh0 to be the average
          call average(mla,sold,rhoh0_old,dx,rhoh_comp)
       end if

       ! set tempbar to be the average
       call average(mla,sold,tempbar,dx,temp_comp)
       tempbar_init = tempbar

       ! reset the time, timestep size, and restart integer
       time = ZERO
       dt = 1.d20
       restart = -1

       ! end of restart_with_vel_field = T

    else

       ! now that we have nr we can read in the base state
       if (restart <= 99999) then
          write(unit=check_index,fmt='(i5.5)') restart
          check_file_name = trim(check_base_name) // check_index
       else if (restart <= 999999) then
          write(unit=check_index6,fmt='(i6.6)') restart
          check_file_name = trim(check_base_name) // check_index6
       else
          write(unit=check_index7,fmt='(i7.7)') restart
          check_file_name = trim(check_base_name) // check_index7
       endif

       ! note: still need to load/store tempbar
       call read_base_state(restart, check_file_name, &
            rho0_old, rhoh0_old, p0_old, gamma1bar, w0, &
            etarho_ec, etarho_cc, div_coeff_old, psi, tempbar, tempbar_init)

       if (do_smallscale) then
          call average(mla,sold,rho0_old,dx,rho_comp)
          call compute_cutoff_coords(rho0_old)
          rho0_old = ZERO
       endif

       ! read in any auxillary data
       call read_aux_data(restart, check_file_name)

    endif

    ! restrict data and fill all ghost cells
    ! this need to be done after read_base_state since in some problems, the inflow
    ! boundary conditions are set in read_base_state
    call ml_restrict_and_fill(nlevs,sold,mba%rr,the_bc_tower%bc_tower_array, &
                              icomp=rho_comp, &
                              bcomp=dm+rho_comp, &
                              nc=nscal, &
                              ng=sold(1)%ng)
    call ml_restrict_and_fill(nlevs,uold,mba%rr,the_bc_tower%bc_tower_array, &
                              icomp=1, &
                              bcomp=1, &
                              nc=dm, &
                              ng=uold(1)%ng)

    if (restart_into_finer .and. spherical .eq. 0) then
       call bl_error('restart_into_finer only currently supported for spherical')
    end if

    if (restart_into_finer .and. time .eq. 0) then
       call bl_error('restart_into_finer does not work if time = 0')
    end if

    if (restart_into_finer .and. spherical .eq. 1) then

       nr_fine_old = nr_fine
       
       allocate(psi_temp      (1,0:nr_fine_old-1))
       allocate(etarho_cc_temp(1,0:nr_fine_old-1))
       allocate(w0_temp       (1,0:nr_fine_old))
       
       ! deallocate the following:
       ! bct%bc_tower_array(i)%phys_bc_level_array
       ! bct%bc_tower_array(i)%adv_bc_level_array
       ! bct%bc_tower_array(i)%ell_bc_level_array
       ! bct%bc_tower_array
       ! bct%domain_bc
       call bc_tower_destroy(the_bc_tower)

       ! this calls bc_tower_init sets bct%nlevels, bct%dim, allocates
       ! bcg%bc_tower_array and bct%domain_bc.  Then bct%bc_tower_array(n)%ngrids
       ! is set to a null value, and bct%domain_bc(:,:) = phys_bc_in(:,:)
       call initialize_bc(the_bc_tower,max_levs,pmask)

       ! build the bc_tower for level 1 only
       call bc_tower_level_build(the_bc_tower,1,mla%la(1))

       ! destroy these before we reset nlevs
       do n=1,nlevs
          call multifab_destroy(S_cc_new(n))
          call multifab_destroy(rho_omegadot2(n))
          call multifab_destroy(rho_Hext(n))
          call multifab_destroy(thermal2(n))
       end do

       ! regrid
       ! this also rebuilds mla and the_bc_tower
       if (use_tpert_in_tagging) then

          ! destroy this before we reset nlevs
          do n=1,nlevs
             call multifab_destroy(rho_Hnuc2(n))
          end do

          ! create tpert
          allocate(tpert_mf(nlevs))
          do n = 1,nlevs
             call build(tpert_mf(n), mla%la(n), 1, 0)
             call multifab_copy_c(tpert_mf(n),1,sold(n),temp_comp,1)
          enddo
          
          call put_in_pert_form(mla,tpert_mf,tempbar,dx,1, &
                                foextrap_comp,.true., &
                                the_bc_tower%bc_tower_array)

          call regrid(restart,mla,uold,sold,gpi,pi,dSdt,S_cc_old, &
                      dx,the_bc_tower, &
                      rho0_old,rhoh0_old,.true.,tpert_mf)

          do n = 1,nlevs
             call destroy(tpert_mf(n))
          enddo
          deallocate(tpert_mf)

       else
          call regrid(restart,mla,uold,sold,gpi,pi,dSdt,S_cc_old, &
                      dx,the_bc_tower, &
                      rho0_old,rhoh0_old,.true.,rho_Hnuc2)

          do n = 1,nlevs
             call destroy(rho_Hnuc2(n))
          enddo
       endif

       ! nlevs is local so we need to reset it
       nlevs = mla%nlevel

       if (nlevs .ne. max_levs) then
          call bl_error('restart_into_finer: nlevs .ne. max_levs not supported yet')
       end if

       ! rebuild these with the new ml_layout
       do n=1,nlevs
          call multifab_build(     S_cc_new(n), mla%la(n),     1, 1)
          call multifab_build(rho_omegadot2(n), mla%la(n), nspec, 0)
          call multifab_build(    rho_Hnuc2(n), mla%la(n),     1, 0)
          call multifab_build(     rho_Hext(n), mla%la(n),     1, 0)
          call multifab_build(     thermal2(n), mla%la(n),     1, 0)
       end do

       ! we set these to zero because they won't affect the solution
       ! they only affect an immediately generated plotfile
       do n=1,nlevs
          call setval(     S_cc_new(n),ZERO,all=.true.)
          call setval(rho_omegadot2(n),ZERO,all=.true.)
          call setval(    rho_Hnuc2(n),ZERO,all=.true.)
          call setval(     rho_Hext(n),ZERO,all=.true.)
          call setval(     thermal2(n),ZERO,all=.true.)
       end do

       ! compute dx at the new level
       deallocate(dx)
       call initialize_dx(dx,mla%mba,nlevs)

       ! compute dr_fine and nr_fine assuming a finer dx
       dr_fine = dx(nlevs,1) / dble(drdxfac)

       if (.not. octant) then
          lenx = HALF * (prob_hi(1) - prob_lo(1))
          leny = HALF * (prob_hi(2) - prob_lo(2))
          lenz = HALF * (prob_hi(3) - prob_lo(3))
       else
          lenx = prob_hi(1) - prob_lo(1)
          leny = prob_hi(2) - prob_lo(2)
          lenz = prob_hi(3) - prob_lo(3)
       end if

       max_dist = sqrt(lenx**2 + leny**2 + lenz**2)
       nr_fine = int(max_dist / dr_fine) + 1
       
       ! compute nr_irreg
       domain = get_pd(get_layout(sold(nlevs)))
       domhi  = upb(domain)+1
       if (.not. octant) then
          nr_irreg = (3*(domhi(1)/2-0.5d0)**2-0.75d0)/2.d0
       else
          nr_irreg = (3*(domhi(1)-0.5d0)**2-0.75d0)/2.d0
       endif


       ! deallocate arrays in geometry.f90 including:
       ! dr,r_cc_loc,r_edge_loc,r_start_coord,r_end_coord,nr,numdisjointchunks,
       ! anelastic_cutoff_coord,base_cutoff_density_coord,burning_cutoff_density_coord
       call destroy_geometry()

       ! set numdisjointchunks = 1
       ! set r_start_coord(1,1) = 0
       ! set r_end_coord(1,1) = nr_fine-1
       call init_multilevel(sold)

       ! initialize arrays in geometry.f90
       call init_radial(nlevs,mla%mba)
       call init_cutoff(nlevs)

       ! make temporary copy of old psi, etarho_cc, and w0
       psi_temp       = psi
       etarho_cc_temp = etarho_cc
       w0_temp        = w0

       ! copy outer pressure for reference
       p0_temp = p0_old(1,nr_fine_old-1)

       ! deallocate 1D arrays
       deallocate(div_coeff_old,div_coeff_new,gamma1bar,gamma1bar_hold,s0_init,rho0_old)
       deallocate(rhoh0_old,rho0_new,rhoh0_new,p0_init,p0_old,p0_new,w0,etarho_ec)
       deallocate(etarho_cc,psi,tempbar,tempbar_init,grav_cell)

       ! reallocate 1D arrays
       call initialize_1d_arrays(nlevs,div_coeff_old,div_coeff_new,gamma1bar, &
                                 gamma1bar_hold,s0_init,rho0_old,rhoh0_old,rho0_new, &
                                 rhoh0_new,p0_init,p0_old,p0_new,w0,etarho_ec,etarho_cc, &
                                 psi,tempbar,tempbar_init,grav_cell)

       ! copy outer pressure for reference
       p0_old(1,nr_fine-1) = p0_temp

       ! fill psi and etarho_cc using linear interpolation
       do r=0,nr_fine-1
          if (r .eq. 0) then
             psi      (1,0) = psi_temp      (1,0)
             etarho_cc(1,0) = etarho_cc_temp(1,0)
          else if ((r+1)/2 .ge. nr_fine_old) then
             psi      (1,r) = psi_temp      (1,nr_fine_old-1)
             etarho_cc(1,r) = etarho_cc_temp(1,nr_fine_old-1)
          else
             if (mod(r,2) .eq. 0) then
                psi      (1,r) = 0.75d0*psi_temp      (1,r/2)+0.25d0*psi_temp      (1,r/2-1)
                etarho_cc(1,r) = 0.75d0*etarho_cc_temp(1,r/2)+0.25d0*etarho_cc_temp(1,r/2-1)
             else
                psi      (1,r) = 0.75d0*psi_temp      (1,r/2)+0.25d0*psi_temp      (1,r/2+1)
                etarho_cc(1,r) = 0.75d0*etarho_cc_temp(1,r/2)+0.25d0*etarho_cc_temp(1,r/2+1)
             end if
          end if
       end do

       ! fill w0 using linear interpolation
       do r=0,nr_fine
          if (r .gt. 2*nr_fine_old) then
             w0(1,r) = w0_temp(1,nr_fine_old)
          else
             if (mod(r,2) .eq. 0) then
                w0(1,r) = w0_temp(1,r/2)
             else
                w0(1,r) = 0.5d0*w0_temp(1,r/2)+0.5d0*w0_temp(1,r/2+1)
             end if
          end if
       end do

       ! put eta on base state edges
       ! note that in spherical the base state has no refinement
       ! the 0th value of etarho = 0, since U dot . e_r must be 
       ! zero at the center (since e_r is not defined there)
       etarho_ec(1,0) = ZERO
       do r=1,nr_fine-1
          etarho_ec(1,r) = HALF*(etarho_cc(1,r) + etarho_cc(1,r-1))
       enddo
       ! probably should do some better extrapolation here eventually
       etarho_ec(1,nr_fine) = etarho_cc(1,nr_fine-1)

       ! compute rho0 by calling average
       call average(mla,sold,rho0_old,dx,rho_comp)
       call compute_cutoff_coords(rho0_old)

       ! compute gravity
       call make_grav_cell(grav_cell,rho0_old)

       ! compute p0 by HSE
       call enforce_HSE(rho0_old,p0_old,grav_cell)

       ! compute temperature with EOS
       if (use_tfromp) then
          ! compute full state T = T(rho,p0,X)
          call makeTfromRhoP(sold,p0_old,mla,the_bc_tower%bc_tower_array,dx)
       else
          ! compute full state T = T(rho,h,X)
          call makeTfromRhoH(sold,p0_old,mla,the_bc_tower%bc_tower_array,dx)
       end if

       ! force tempbar to be the average of temp
       call average(mla,sold,tempbar,dx,temp_comp)
       
       ! for restarting into finer, we lose the drive_initial_convection 
       ! functionality, so we just set tempbar_init = tempbar
       tempbar_init = tempbar

       ! compute gamma1 (just for use in computing gamma1bar)
       allocate(gamma1(nlevs))
       do n=1,nlevs
          call multifab_build(gamma1(n), mla%la(n), 1, 0)
       end do
       call make_gamma(mla,gamma1,sold,p0_old,dx)

       ! compute gamma1bar
       call average(mla,gamma1,gamma1bar,dx,1)

       ! deallocate gamma1
       do n=1,nlevs
          call destroy(gamma1(n))
       end do
       deallocate(gamma1)
     
       ! compute div_coeff_old
       call make_div_coeff(div_coeff_old,rho0_old,p0_old,gamma1bar,grav_cell)

       ! recompute time step
       dt = 1.d20
       call estdt(mla,the_bc_tower,uold,sold,gpi,S_cc_old,dSdt, &
                  w0,rho0_old,p0_old,gamma1bar,grav_cell,dx,cflfac,dt)

    end if ! end spherical restart_into_finer initialization

    call destroy(mba_old)
    call destroy(mba)

  end subroutine initialize_from_restart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initialize_with_fixed_grids(mla,dt,pmask,dx,uold,sold,gpi,pi, &
                                         dSdt,S_cc_old,S_cc_new, &
                                         rho_omegadot2,rho_Hnuc2,rho_Hext, &
                                         thermal2, &
                                         the_bc_tower,div_coeff_old,div_coeff_new, &
                                         gamma1bar,gamma1bar_hold,s0_init,rho0_old, &
                                         rhoh0_old,rho0_new,rhoh0_new,p0_init, &
                                         p0_old,p0_new,w0,etarho_ec,etarho_cc, &
                                         psi,tempbar,tempbar_init,grav_cell)

    use box_util_module
    use init_scalar_module
    use init_vel_module
    use average_module
    use restrict_base_module
    use probin_module, only : drdxfac, octant, test_set, ppm_type, bds_type, nodal, &
         prob_lo, prob_hi, model_file, do_smallscale, dm_in, fix_base_state
    use make_grav_module
    use enforce_HSE_module
    use rhoh_vs_t_module
    use time_module, only: time
    
    type(ml_layout),intent(out  ) :: mla
    real(dp_t)    , intent(inout) :: dt
    logical       , intent(in   ) :: pmask(:)
    real(dp_t)    , pointer       :: dx(:,:)
    type(multifab), pointer       :: uold(:),sold(:),gpi(:),pi(:),dSdt(:)
    type(multifab), pointer       :: S_cc_old(:),S_cc_new(:)
    type(multifab), pointer       :: rho_omegadot2(:),rho_Hnuc2(:),rho_Hext(:),thermal2(:)
    type(bc_tower), intent(  out) :: the_bc_tower
    real(dp_t)    , pointer       :: div_coeff_old(:,:),div_coeff_new(:,:),gamma1bar(:,:)
    real(dp_t)    , pointer       :: gamma1bar_hold(:,:),s0_init(:,:,:),rho0_old(:,:)
    real(dp_t)    , pointer       :: rhoh0_old(:,:),rho0_new(:,:),rhoh0_new(:,:),p0_init(:,:)
    real(dp_t)    , pointer       :: p0_old(:,:),p0_new(:,:),w0(:,:),etarho_ec(:,:)
    real(dp_t)    , pointer       :: etarho_cc(:,:),psi(:,:),tempbar(:,:),tempbar_init(:,:),grav_cell(:,:)

    ! local
    type(ml_boxarray) :: mba
    type(box)         :: domain
    integer           :: domhi(dm_in)

    real(dp_t) :: lenx,leny,lenz,max_dist

    integer :: n,ng_s,dm,nlevs
    
    ! set time and dt
    time = ZERO
    dt = 1.d20

    ! create mba
    call read_a_hgproj_grid(mba,test_set)

    ! create mla
    call ml_layout_build(mla,mba,pmask)
    
    dm = mla%dim

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
    allocate(uold(nlevs),sold(nlevs),gpi(nlevs),pi(nlevs))
    allocate(dSdt(nlevs),S_cc_old(nlevs),S_cc_new(nlevs))
    allocate(rho_omegadot2(nlevs),rho_Hnuc2(nlevs),rho_Hext(nlevs),thermal2(nlevs))

    if (ppm_type .eq. 2 .or. bds_type .eq. 1) then
       ng_s = 4
    else
       ng_s = 3
    end if

    ! build states
    do n = 1,nlevs
       call multifab_build(         uold(n), mla%la(n),    dm, ng_s)
       call multifab_build(         sold(n), mla%la(n), nscal, ng_s)
       call multifab_build(          gpi(n), mla%la(n),    dm, 0)
       call multifab_build(           pi(n), mla%la(n),     1, 0, nodal)
       call multifab_build(         dSdt(n), mla%la(n),     1, 0)
       call multifab_build(     S_cc_old(n), mla%la(n),     1, 0)
       call multifab_build(     S_cc_new(n), mla%la(n),     1, 0)
       call multifab_build(rho_omegadot2(n), mla%la(n), nspec, 0)
       call multifab_build(    rho_Hnuc2(n), mla%la(n),     1, 0)
       call multifab_build(     rho_Hext(n), mla%la(n),     1, 0)
       call multifab_build(     thermal2(n), mla%la(n),     1, 0)

       call setval(         uold(n), ZERO, all=.true.)
       call setval(         sold(n), ZERO, all=.true.)
       call setval(          gpi(n), ZERO, all=.true.)
       call setval(           pi(n), ZERO, all=.true.)
       call setval(     S_cc_old(n), ZERO, all=.true.)
       call setval(     S_cc_new(n), ZERO, all=.true.)
       call setval(         dSdt(n), ZERO, all=.true.)
       call setval(rho_omegadot2(n), ZERO, all=.true.)
       call setval(    rho_Hnuc2(n), ZERO, all=.true.)
       call setval(     rho_Hext(n), ZERO, all=.true.)
       call setval(     thermal2(n), ZERO, all=.true.)
    end do
    ! initialize dx
    call initialize_dx(dx,mba,nlevs)

    ! initialize cutoff arrays
    call init_cutoff(nlevs)

    ! now that we have dx we can initialize nr_fine and dr_fine
    if (spherical .eq. 1) then

       dr_fine = dx(nlevs,1) / dble(drdxfac)
       
       if (.not. octant) then
          lenx = HALF * (prob_hi(1) - prob_lo(1))
          leny = HALF * (prob_hi(2) - prob_lo(2))
          lenz = HALF * (prob_hi(3) - prob_lo(3))
       else
          lenx = prob_hi(1) - prob_lo(1)
          leny = prob_hi(2) - prob_lo(2)
          lenz = prob_hi(3) - prob_lo(3)
       end if
       
       max_dist = sqrt(lenx**2 + leny**2 + lenz**2)
       nr_fine = int(max_dist / dr_fine) + 1

       ! compute nr_irreg
       domain = get_pd(get_layout(sold(nlevs)))
       domhi  = upb(domain)+1
       if (.not.octant) then
          nr_irreg = (3*(domhi(1)/2-0.5d0)**2-0.75d0)/2.d0
       else
          nr_irreg = (3*(domhi(1)-0.5d0)**2-0.75d0)/2.d0
       endif

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
                              p0_old,p0_new,w0,etarho_ec,etarho_cc,psi,tempbar,tempbar_init, &
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

    p0_old       = p0_init
    rho0_old     = s0_init(:,:,rho_comp)
    rhoh0_old    = s0_init(:,:,rhoh_comp)
    tempbar      = s0_init(:,:,temp_comp)
    tempbar_init = s0_init(:,:,temp_comp)

    if (fix_base_state) then
       call compute_cutoff_coords(rho0_old)
       call make_grav_cell(grav_cell,rho0_old)
       call destroy(mba)
       return
    end if

    if (do_smallscale) then
       ! first compute cutoff coordinates using initial density profile
       call compute_cutoff_coords(rho0_old)
       ! set rho0_old = rhoh0_old = ZERO
       rho0_old  = ZERO
       rhoh0_old = ZERO
    else
       ! set rho0 to be the average
       call average(mla,sold,rho0_old,dx,rho_comp)
       call compute_cutoff_coords(rho0_old)

       ! compute p0 with HSE
       call make_grav_cell(grav_cell,rho0_old)
       call enforce_HSE(rho0_old,p0_old,grav_cell)

       ! call eos with r,p as input to recompute T,h
       call makeTHfromRhoP(sold,p0_old,the_bc_tower%bc_tower_array,mla,dx)

       ! set rhoh0 to be the average
       call average(mla,sold,rhoh0_old,dx,rhoh_comp)
    end if

    ! set tempbar to be the average
    call average(mla,sold,tempbar,dx,temp_comp)
    tempbar_init = tempbar

    call destroy(mba)

  end subroutine initialize_with_fixed_grids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initialize_with_adaptive_grids(mla,dt,pmask,dx,uold,sold,gpi,pi, &
                                            dSdt,S_cc_old,S_cc_new, &
                                            rho_omegadot2,rho_Hnuc2,rho_Hext, &
                                            thermal2, &
                                            the_bc_tower,div_coeff_old,div_coeff_new, &
                                            gamma1bar,gamma1bar_hold,s0_init,rho0_old, &
                                            rhoh0_old,rho0_new,rhoh0_new,p0_init, &
                                            p0_old,p0_new,w0,etarho_ec,etarho_cc, &
                                            psi,tempbar,tempbar_init,grav_cell)

    use probin_module, only: n_cellx, n_celly, n_cellz, &
         amr_buf_width, max_grid_size_1, max_grid_size_2, max_grid_size_3, &
         ref_ratio, max_levs, octant
    use init_scalar_module
    use init_vel_module
    use average_module
    use restrict_base_module
    use make_new_grids_module
    use probin_module, only : drdxfac, ppm_type, bds_type, prob_lo, prob_hi, do_smallscale, &
         model_file, nodal, dm_in, fix_base_state, use_tpert_in_tagging
    use ml_restrict_fill_module
    use make_grav_module
    use enforce_HSE_module
    use rhoh_vs_t_module
    use time_module, only: time
    use pert_form_module, only: put_in_pert_form_one_level

    type(ml_layout),intent(out  ) :: mla
    real(dp_t)    , intent(inout) :: dt
    logical       , intent(in   ) :: pmask(:)
    real(dp_t)    , pointer       :: dx(:,:)
    type(multifab), pointer       :: uold(:),sold(:),gpi(:),pi(:),dSdt(:)
    type(multifab), pointer       :: S_cc_old(:),S_cc_new(:)
    type(multifab), pointer       :: rho_omegadot2(:),rho_Hnuc2(:),rho_Hext(:),thermal2(:)
    type(bc_tower), intent(  out) :: the_bc_tower
    real(dp_t)    , pointer       :: div_coeff_old(:,:),div_coeff_new(:,:),gamma1bar(:,:)
    real(dp_t)    , pointer       :: gamma1bar_hold(:,:),s0_init(:,:,:),rho0_old(:,:)
    real(dp_t)    , pointer       :: rhoh0_old(:,:),rho0_new(:,:),rhoh0_new(:,:),p0_init(:,:)
    real(dp_t)    , pointer       :: p0_old(:,:),p0_new(:,:),w0(:,:),etarho_ec(:,:)
    real(dp_t)    , pointer       :: etarho_cc(:,:),psi(:,:),tempbar(:,:),tempbar_init(:,:),grav_cell(:,:)

    ! local
    type(ml_boxarray) :: mba
    type(box)         :: domain
    integer           :: domhi(dm_in)

    type(multifab), pointer :: tpert_mf(:)

    type(layout) :: la_array(max_levs)
    type(box)    :: bxs

    real(dp_t) :: lenx,leny,lenz,max_dist
    integer    :: n,ng_s,nl
    integer    :: lo(dm_in), hi(dm_in), dm, nlevs
    logical    :: new_grid

    ! set time and dt
    time = ZERO
    dt = 1.d20

    dm = dm_in

    ! set up hi & lo to carry indexing info
    lo = 0
    hi(1) = n_cellx - 1
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
    allocate(uold(max_levs),sold(max_levs),gpi(max_levs),pi(max_levs))
    allocate(dSdt(max_levs),S_cc_old(max_levs),S_cc_new(max_levs))
    allocate(rho_omegadot2(max_levs),rho_Hnuc2(max_levs),rho_Hext(max_levs),thermal2(max_levs))

    ! Build the level 1 boxarray
    call box_build_2(bxs,lo,hi)
    call boxarray_build_bx(mba%bas(1),bxs)
    call boxarray_maxsize(mba%bas(1),max_grid_size_1)

    ! build pd(:)
    mba%pd(1) = bxs
    do n = 2, max_levs
       mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
    enddo

    ! initialize dx
    call initialize_dx(dx,mba,max_levs)

    ! initialize cutoff arrays
    call init_cutoff(max_levs)

    ! set number of ghost cells required by hydrodynamics stencil
    if (ppm_type .eq. 2 .or. bds_type .eq. 1) then
       ng_s = 4
    else
       ng_s = 3
    end if

    ! now that we have dx we can initialize nr_fine and dr_fine
    if (spherical .eq. 1) then

       dr_fine = dx(max_levs,1) / dble(drdxfac)

       if (.not. octant) then
          lenx = HALF * (prob_hi(1) - prob_lo(1))
          leny = HALF * (prob_hi(2) - prob_lo(2))
          lenz = HALF * (prob_hi(3) - prob_lo(3))
       else
          lenx = prob_hi(1) - prob_lo(1)
          leny = prob_hi(2) - prob_lo(2)
          lenz = prob_hi(3) - prob_lo(3)
       end if

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
                              psi,tempbar,tempbar_init,grav_cell)

    ! now that we have dr and nr we can fill initial state
    if (spherical .eq. 1) then
       call init_base_state(1,model_file,s0_init(1,:,:),p0_init(1,:),dx(max_levs,:))
    else
          ! init_base_state requires loop backwards over levels
       do n=max_levs,1,-1
          call init_base_state(n,model_file,s0_init(n,:,:),p0_init(n,:),dx(n,:))
       end do
    end if

    ! Initialize bc's
    call initialize_bc(the_bc_tower,max_levs,pmask)

    ! Build the level 1 layout.
    call layout_build_ba(la_array(1),mba%bas(1),mba%pd(1),pmask)

    ! Define bc_tower at level 1.
    call bc_tower_level_build(the_bc_tower,1,la_array(1))

    nlevs = 1
    nlevs_radial = 1

    if (max_levs > 1) then
       ! Build and initialize the level 1 data 
       ! so we can use it to build finer levels
       call multifab_build(sold(1), la_array(1), nscal, ng_s)
       call initscalardata_on_level(1,sold(1),s0_init(1,:,:),p0_init(1,:), &
                                    dx(1,:),the_bc_tower%bc_tower_array(1))

       ! Fill level 1's ghost cells and apply boundary conditions 
       ! so we can use them in tagging boxes for refinement
       call ml_restrict_and_fill(1,sold,mba%rr,the_bc_tower%bc_tower_array, &
                                 icomp=rho_comp, &
                                 bcomp=dm+rho_comp, &
                                 nc=nscal, &
                                 ng=sold(1)%ng)

       new_grid = .true.
       nl = 1

       ! Loop over all levels, building fine grids as we go.  This loop's
       ! ultimate purpose is to add properly refined grids to 'mba' so we
       ! can use it to build 'mla'
       do while ( (nl .lt. max_levs) .and. (new_grid) )
          if (nl >  1) then
             ! We made a new fine level that needs to be filled.
             !   Note: We fill the data here instead of at the bottom of the
             !   loop so that we only spend resources building and filling
             !   temporary data structures if we need to do another refinement. 

             ! Build the level nl data only.
             call multifab_build(sold(nl),la_array(nl),nscal,ng_s)
             
             ! Define bc_tower at level nl.
             call bc_tower_level_build(the_bc_tower,nl,la_array(nl))
            
             ! Fill the nl level valid region 
             if (spherical .eq. 1) then
                call initscalardata_on_level(nl,sold(nl),s0_init(1,:,:), &
                                             p0_init(1,:),dx(nl,:), &
                                             the_bc_tower%bc_tower_array(nl))
             else
                call initscalardata_on_level(nl,sold(nl),s0_init(nl,:,:), &
                                             p0_init(nl,:),dx(nl,:), &
                                             the_bc_tower%bc_tower_array(nl))
             end if

             ! restrict data and fill all ghost cells
             call ml_restrict_and_fill(nl,sold,mba%rr,the_bc_tower%bc_tower_array, &
                                       icomp=rho_comp, &
                                       bcomp=dm+rho_comp, &
                                       nc=nscal, &
                                       ng=sold(1)%ng)
          endif

          if (use_tpert_in_tagging) then

             call average_one_level(nl,sold,tempbar,temp_comp)

             allocate(tpert_mf(max_levs))
             do n = 1, nl
                call build(tpert_mf(n), la_array(n), 1, 0)
                call multifab_copy_c(tpert_mf(n),1,sold(n),temp_comp,1)
             enddo

             call put_in_pert_form_one_level(nl,tpert_mf(nl),tempbar,dx(nl,:),1,.true.)

             ! Do we need finer grids?  
             !   Note: currently the answer to this question is always yes.  Maestro 
             !   doesn't currently support having fewer levels than max_levs
             if (nl .eq. 1) then
                call make_new_grids(new_grid,la_array(nl),la_array(nl+1),sold(nl),dx(nl,1), &
                                    amr_buf_width,ref_ratio,nl,max_grid_size_2,tpert_mf(nl))
             else
                call make_new_grids(new_grid,la_array(nl),la_array(nl+1),sold(nl),dx(nl,1), &
                                 amr_buf_width,ref_ratio,nl,max_grid_size_3,tpert_mf(nl))
             end if

             do n=1,nl
                call destroy(tpert_mf(n))
             enddo
             deallocate(tpert_mf)
          
          else
             ! Do we need finer grids?  
             !   Note: currently the answer to this question is always yes.  Maestro 
             !   doesn't currently support having fewer levels than max_levs
             if (nl .eq. 1) then
                call make_new_grids(new_grid,la_array(nl),la_array(nl+1),sold(nl),dx(nl,1), &
                                    amr_buf_width,ref_ratio,nl,max_grid_size_2)
             else
                call make_new_grids(new_grid,la_array(nl),la_array(nl+1),sold(nl),dx(nl,1), &
                                    amr_buf_width,ref_ratio,nl,max_grid_size_3)
             end if

          endif
         
          ! If we tagged boxes for refinement, then build the new grids! 
          if (new_grid) then
             ! Add the new level's boxes to mba's boxarray list
             call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))
            
             ! We must enforce proper nesting as we build the grids so that we
             ! can fillpatch() them
             !   Note: Proper nesting only needs to be enforced for nl >= 2
             !   because Maestro requires that level 1 cover the entire domain
             !   and thus level 2 will always be properly nested.
             if(nl .ge. 2) then
               ! Test if grids are already properly nested
               if (.not. ml_boxarray_properly_nested(mba, ng_s, pmask, 2, nl+1)) then
                 do n = 2,nl  
                   ! Delete old multifabs so that we can rebuild them.
                   call destroy(  sold(n))
                 enddo

                 ! This changes mba and la_array such that all levels are properly nested
                 call enforce_proper_nesting(mba,la_array,max_grid_size_2,max_grid_size_3)

                 ! Loop over all the lower levels which we might have changed when we enforced proper nesting.
                 do n = 2,nl
                   ! This makes sure the boundary conditions are properly defined everywhere
                   call bc_tower_level_build(the_bc_tower,n,la_array(n)) 
   
                   ! Rebuild and fill the lower level data.
                   ! TODO: It would be more efficient to only destroy/rebuild
                   !   the data when we know for sure enforce_proper_nesting()
                   !   changed the layout.
                   call multifab_build(sold(n),la_array(n),nscal,ng_s)
                   ! fills the physical (valid) region of each level with problem data
                   if (spherical .eq. 1) then
                     call initscalardata_on_level(n,sold(n),s0_init(1,:,:), &
                                                  p0_init(1,:),dx(n,:), &
                                                  the_bc_tower%bc_tower_array(n))
                   else
                     call initscalardata_on_level(n,sold(n),s0_init(n,:,:), &
                                                  p0_init(n,:),dx(n,:), &
                                                  the_bc_tower%bc_tower_array(n))
                   end if

                 end do

                 ! restrict data and fill all ghost cells
                 call ml_restrict_and_fill(nl,sold,mba%rr,the_bc_tower%bc_tower_array, &
                                           icomp=rho_comp, &
                                           bcomp=dm+rho_comp, &
                                           nc=nscal, &
                                           ng=sold(1)%ng)

               endif
             endif !end proper nesting enforcement

             nlevs = nl+1
             nlevs_radial = merge(1, nlevs, spherical .eq. 1)
             nl = nl + 1
             
          endif ! if (new_grid) 
          
       enddo
       
       do n = 1,nlevs-1
          call destroy(sold(n))
       end do

       nlevs = nl
       nlevs_radial = merge(1, nlevs, spherical .eq. 1)

       ! A final enforcement of proper nesting
       if (nlevs .ge. 3) then
          call enforce_proper_nesting(mba,la_array,max_grid_size_2,max_grid_size_3)
       end if
       
    end if  ! end if (maxlev > 1)
    
    do n = 1, nlevs
       call destroy(la_array(n))
    end do

    call ml_layout_restricted_build(mla,mba,nlevs,pmask)
    
    nlevs = mla%nlevel
    
    if (nlevs .ne. max_levs) then
       call bl_error('initialize_with_adaptive_grids: nlevs .ne. max_levs not supported yet')
    end if

    ! this makes sure the boundary conditions are properly defined everywhere
    do n = 1,nlevs
       call bc_tower_level_build(the_bc_tower,n,mla%la(n))
    end do

    nlevs_radial = merge(1, nlevs, spherical .eq. 1)
    
    ! build states
    do n = 1,nlevs
       call multifab_build(         uold(n), mla%la(n),    dm, ng_s)
       call multifab_build(         sold(n), mla%la(n), nscal, ng_s)
       call multifab_build(          gpi(n), mla%la(n),    dm, 0)
       call multifab_build(           pi(n), mla%la(n),     1, 0, nodal)
       call multifab_build(         dSdt(n), mla%la(n),     1, 0)
       call multifab_build(     S_cc_old(n), mla%la(n),     1, 0)
       call multifab_build(     S_cc_new(n), mla%la(n),     1, 0)
       call multifab_build(rho_omegadot2(n), mla%la(n), nspec, 0)
       call multifab_build(    rho_Hnuc2(n), mla%la(n),     1, 0)
       call multifab_build(     rho_Hext(n), mla%la(n),     1, 0)
       call multifab_build(     thermal2(n), mla%la(n),     1, 0)

       call setval(         uold(n), ZERO, all=.true.)
       call setval(         sold(n), ZERO, all=.true.)
       call setval(          gpi(n), ZERO, all=.true.)
       call setval(           pi(n), ZERO, all=.true.)
       call setval(     S_cc_old(n), ZERO, all=.true.)
       call setval(     S_cc_new(n), ZERO, all=.true.)
       call setval(         dSdt(n), ZERO, all=.true.)
       call setval(rho_omegadot2(n), ZERO, all=.true.)
       call setval(    rho_Hnuc2(n), ZERO, all=.true.)
       call setval(     rho_Hext(n), ZERO, all=.true.)
       call setval(     thermal2(n), ZERO, all=.true.)
    end do

    ! compute nr_irreg
    domain = get_pd(get_layout(sold(nlevs)))
    domhi  = upb(domain)+1
    if (.not. octant) then
       nr_irreg = (3*(domhi(1)/2-0.5d0)**2-0.75d0)/2.d0
    else
       nr_irreg = (3*(domhi(1)-0.5d0)**2-0.75d0)/2.d0
    endif

    ! create numdisjointchunks, r_start_coord, r_end_coord
    call init_multilevel(sold)

    call initveldata(uold,s0_init,p0_init,dx,the_bc_tower%bc_tower_array,mla)
    call initscalardata(sold,s0_init,p0_init,dx,the_bc_tower%bc_tower_array,mla)

    p0_old       = p0_init
    rho0_old     = s0_init(:,:,rho_comp)
    rhoh0_old    = s0_init(:,:,rhoh_comp)
    tempbar      = s0_init(:,:,temp_comp)
    tempbar_init = s0_init(:,:,temp_comp)

    if (fix_base_state) then
       call compute_cutoff_coords(rho0_old)
       call make_grav_cell(grav_cell,rho0_old)
       call destroy(mba)
       return
    end if

    if (do_smallscale) then
       ! first compute cutoff coordinates using initial density profile
       call compute_cutoff_coords(rho0_old)
       ! set rho0_old = rhoh0_old = ZERO
       rho0_old  = ZERO
       rhoh0_old = ZERO
    else
       ! set rho0 to be the average
       call average(mla,sold,rho0_old,dx,rho_comp)
       call compute_cutoff_coords(rho0_old)

       ! compute p0 with HSE
       call make_grav_cell(grav_cell,rho0_old)
       call enforce_HSE(rho0_old,p0_old,grav_cell)

       ! call eos with r,p as input to recompute T,h
       call makeTHfromRhoP(sold,p0_old,the_bc_tower%bc_tower_array,mla,dx)

       ! set rhoh0 to be the average
       call average(mla,sold,rhoh0_old,dx,rhoh_comp)
    end if

    ! set tempbar to be the average
    call average(mla,sold,tempbar,dx,temp_comp)
    tempbar_init = tempbar

    call destroy(mba)

  end subroutine initialize_with_adaptive_grids


  subroutine initialize_1d_arrays(num_levs,div_coeff_old,div_coeff_new,gamma1bar, &
                                  gamma1bar_hold,s0_init,rho0_old,rhoh0_old,rho0_new, &
                                  rhoh0_new,p0_init,p0_old,p0_new,w0,etarho_ec,etarho_cc, &
                                  psi,tempbar,tempbar_init,grav_cell)

    integer    , intent(in) :: num_levs    
    real(dp_t) , pointer    :: div_coeff_old(:,:),div_coeff_new(:,:),gamma1bar(:,:)
    real(dp_t) , pointer    :: gamma1bar_hold(:,:),s0_init(:,:,:),rho0_old(:,:)
    real(dp_t) , pointer    :: rhoh0_old(:,:),rho0_new(:,:),rhoh0_new(:,:),p0_init(:,:)
    real(dp_t) , pointer    :: p0_old(:,:),p0_new(:,:),w0(:,:),etarho_ec(:,:),etarho_cc(:,:)
    real(dp_t) , pointer    :: psi(:,:),tempbar(:,:),tempbar_init(:,:),grav_cell(:,:)
    
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
       allocate(tempbar_init  (num_levs,0:nr_fine-1))
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
       allocate(tempbar_init  (1,0:nr_fine-1))
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
    tempbar_init = ZERO
    grav_cell = ZERO

  end subroutine initialize_1d_arrays
  
end module initialize_module
