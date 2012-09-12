module regrid_module

  use BoxLib
  use ml_boxarray_module
  use ml_layout_module
  use define_bc_module
  use restart_module
  use init_module
  use box_util_module
  use make_new_grids_module
  use initialize_module
  use fillpatch_module
  use ml_prolongation_module
  use multifab_physbc_module
  use multifab_fill_ghost_module
  use ml_restriction_module
  use pert_form_module
  use convert_rhoX_to_X_module

  implicit none

  private

  public :: regrid

contains

  subroutine regrid(mla,uold,sold,gpres,pres,dSdt,src,dx,the_bc_tower,rho0,rhoh0)

    use probin_module, only : nodal, pmask, regrid_int, max_grid_size, ref_ratio, max_levs, &
         ppm_type
    use geometry, only: dm, nlevs, nlevs_radial, spherical, nr_fine
    use variables, only: nscal, rho_comp, rhoh_comp, foextrap_comp
    use network, only: nspec
    use average_module

    type(ml_layout), intent(inout) :: mla
    type(multifab),  pointer       :: uold(:),sold(:),gpres(:),pres(:)
    type(multifab),  pointer       :: dSdt(:),src(:)
    real(dp_t)    ,  pointer       :: dx(:,:)
    type(bc_tower),  intent(inout) :: the_bc_tower
    real(kind=dp_t), intent(in   ) :: rho0(:,0:),rhoh0(:,0:)

    ! local
    logical           :: new_grid
    integer           :: n, nl, buf_wid, d, ng_s
    type(layout)      :: la_array(max_levs)
    type(ml_layout)   :: mla_old
    type(ml_boxarray) :: mba
    real(kind=dp_t)   :: tempbar(max_levs,0:nr_fine-1)

    ! These are copies to hold the old data.
    type(multifab) :: uold_temp(nlevs), sold_temp(nlevs), gpres_temp(nlevs), pres_temp(nlevs)
    type(multifab) :: dSdt_temp(nlevs), src_temp(nlevs)

    if (ppm_type .eq. 2) then
       ng_s = 4
    else
       ng_s = 3
    end if

    if (max_levs < 2) then
       call bl_error('Dont call regrid with max_levs < 2')
    end if

    call ml_layout_build(mla_old,mla%mba,mla%pmask)

    do n = 1,nlevs

       ! Create copies of the old data.
       call multifab_build(  uold_temp(n),mla_old%la(n),   dm, ng_s)
       call multifab_build(  sold_temp(n),mla_old%la(n),nscal, ng_s)
       call multifab_build( gpres_temp(n),mla_old%la(n),   dm, 1)
       call multifab_build(  pres_temp(n),mla_old%la(n),    1, 1, nodal)
       call multifab_build(  dSdt_temp(n),mla_old%la(n),    1, 0)
       call multifab_build(   src_temp(n),mla_old%la(n),    1, 1)

       call multifab_copy_c(  uold_temp(n),1,  uold(n),1,   dm)
       call multifab_copy_c(  sold_temp(n),1,  sold(n),1,nscal)
       call multifab_copy_c( gpres_temp(n),1, gpres(n),1,   dm)
       call multifab_copy_c(  pres_temp(n),1,  pres(n),1,    1)
       call multifab_copy_c(  dSdt_temp(n),1,  dSdt(n),1,    1)
       call multifab_copy_c(   src_temp(n),1,   src(n),1,    1)

       ! Get rid of the old data structures so we can create new ones 
       ! with the same names.
       call multifab_destroy(  uold(n))
       call multifab_destroy(  sold(n))
       call multifab_destroy( gpres(n))
       call multifab_destroy(  pres(n))
       call multifab_destroy(  dSdt(n))
       call multifab_destroy(   src(n))

    end do

    call destroy(mla)

    buf_wid = regrid_int

    ! mba is big enough to hold max_levs levels
    ! even though we know we had nlevs last time, we might 
    ! want more or fewer levels after regrid (if nlevs < max_levs)
    call ml_boxarray_build_n(mba,max_levs,dm)

    do n = 1, max_levs-1
       mba%rr(n,:) = ref_ratio
    enddo

    if (associated(uold)) then
       deallocate(uold,sold,pres,gpres,dSdt,src)
    end if

    allocate(uold(max_levs),sold(max_levs),pres(max_levs),gpres(max_levs))
    allocate(dSdt(max_levs),src(max_levs))

    ! Copy the level 1 boxarray
    call copy(mba%bas(1),mla_old%mba%bas(1))

    ! Copy the pd(:)
    mba%pd(1) = mla_old%mba%pd(1)
    do n = 2, max_levs
       mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
    enddo

    ! Build the level 1 layout.
    call layout_build_ba(la_array(1),mba%bas(1),mba%pd(1),pmask)

    ! Build the level 1 data only.
    call multifab_build(  uold(1), la_array(1),    dm, ng_s)
    call multifab_build(  sold(1), la_array(1), nscal, ng_s)
    call multifab_build( gpres(1), la_array(1),    dm, 1)
    call multifab_build(  pres(1), la_array(1),     1, 1, nodal)
    call multifab_build(  dSdt(1), la_array(1),     1, 0)
    call multifab_build(   src(1), la_array(1),     1, 1)

    ! Copy the level 1 data from the "old" temporaries.
    call multifab_copy_c(  uold(1),1,  uold_temp(1) ,1,   dm)
    call multifab_copy_c(  sold(1),1,  sold_temp(1) ,1,nscal)
    call multifab_copy_c( gpres(1),1, gpres_temp(1), 1,   dm)
    call multifab_copy_c(  pres(1),1,  pres_temp(1) ,1,    1)
    call multifab_copy_c(  dSdt(1),1,  dSdt_temp(1), 1,    1)
    call multifab_copy_c(   src(1),1,   src_temp(1), 1,    1)

    nl       = 1
    new_grid = .true.

    do while ( (nl .lt. max_levs) .and. (new_grid) )

       ! Do we need finer grids?

       ! Need to fill ghost cells here in case we use them in tagging
       call multifab_fill_boundary(sold(nl))
       call multifab_physbc(sold(nl),rho_comp,dm+rho_comp,nscal, &
                            the_bc_tower%bc_tower_array(nl))

       call average_one_level(nl,sold,tempbar,temp_comp)

       call make_new_grids(new_grid,la_array(nl),la_array(nl+1),sold(nl),dx(nl,1),buf_wid,&
                           ref_ratio,nl,max_grid_size,tempbar)

       if (new_grid) then

          call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))

          ! Build the level nl+1 data only.
          call multifab_build(  uold(nl+1), la_array(nl+1),    dm, ng_s)
          call multifab_build(  sold(nl+1), la_array(nl+1), nscal, ng_s)
          call multifab_build( gpres(nl+1), la_array(nl+1),    dm, 1)
          call multifab_build(  pres(nl+1), la_array(nl+1),     1, 1, nodal)
          call multifab_build(  dSdt(nl+1), la_array(nl+1),     1, 0)
          call multifab_build(   src(nl+1), la_array(nl+1),     1, 1)

          ! Define bc_tower at level nl+1.
          call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))

          ! Fill the data in the new level nl+1 state -- first from the coarser data.
          call fillpatch(uold(nl+1),uold(nl), &
                         ng_s,mba%rr(nl,:), &
                         the_bc_tower%bc_tower_array(nl  ), &
                         the_bc_tower%bc_tower_array(nl+1), &
                         1,1,1,dm)
          call fillpatch(sold(nl+1),sold(nl), &
                         ng_s,mba%rr(nl,:), &
                         the_bc_tower%bc_tower_array(nl  ), &
                         the_bc_tower%bc_tower_array(nl+1), &
                         1,1,dm+rho_comp,nscal)
          do d=1,dm
             call fillpatch(gpres(nl+1),gpres(nl), &
                            1,mba%rr(nl,:), &
                            the_bc_tower%bc_tower_array(nl  ), &
                            the_bc_tower%bc_tower_array(nl+1), &
                            d,d,foextrap_comp,1)
          end do
          call fillpatch(dSdt(nl+1),dSdt(nl), &
                         0,mba%rr(nl,:), &
                         the_bc_tower%bc_tower_array(nl  ), &
                         the_bc_tower%bc_tower_array(nl+1), &
                         1,1,foextrap_comp,1)
          call fillpatch(src(nl+1),src(nl), &
                         1,mba%rr(nl,:), &
                         the_bc_tower%bc_tower_array(nl  ), &
                         the_bc_tower%bc_tower_array(nl+1), &
                         1,1,foextrap_comp,1)

          ! We interpolate p differently because it is nodal, not cell-centered
          call ml_prolongation(pres(nl+1),pres(nl),layout_get_pd(la_array(nl+1)), &
                               mba%rr(nl,:))

          ! Copy from old data at current level, if it exists
          if (mla_old%nlevel .ge. nl+1) then
             call multifab_copy_c(  uold(nl+1),1,  uold_temp(nl+1),1,   dm)
             call multifab_copy_c(  sold(nl+1),1,  sold_temp(nl+1),1,nscal)
             call multifab_copy_c( gpres(nl+1),1, gpres_temp(nl+1),1,   dm)
             call multifab_copy_c(  pres(nl+1),1,  pres_temp(nl+1),1,    1)
             call multifab_copy_c(  dSdt(nl+1),1,  dSdt_temp(nl+1),1,    1)
             call multifab_copy_c(   src(nl+1),1,   src_temp(nl+1),1,    1)
          end if

          nlevs = nl+1
          nlevs_radial = merge(1, nlevs, spherical .eq. 1)
          nl = nl + 1

       endif

    enddo

    if (spherical .eq. 1) then

       ! convert (rho X) --> X in sold 
       call convert_rhoX_to_X(sold_temp,.true.,mla,the_bc_tower%bc_tower_array)

       ! convert rho -> rho' in sold_temp
       call put_in_pert_form(mla_old,sold_temp,rho0,dx,rho_comp,foextrap_comp,.true., &
                             the_bc_tower%bc_tower_array)

       ! convert (rho h) -> (rho h)' in sold_temp
       call put_in_pert_form(mla_old,sold_temp,rhoh0,dx,rhoh_comp,foextrap_comp,.true., &
                             the_bc_tower%bc_tower_array)

    end if

    do n = 1,nl
       call destroy(  sold(n))
       call destroy(  uold(n))
       call destroy( gpres(n))
       call destroy(  pres(n))
       call destroy(  dSdt(n))
       call destroy(   src(n))
    end do

    nlevs = nl
    nlevs_radial = merge(1, nlevs, spherical .eq. 1)

    call ml_layout_restricted_build(mla,mba,nlevs,pmask)

    ! check for proper nesting
    if (nlevs .ge. 3) then
       call enforce_proper_nesting(mba,la_array,max_grid_size)
    end if

    do n = 1,nl
       call destroy(la_array(n))
    end do

    call destroy(mla)

    call ml_layout_restricted_build(mla,mba,nlevs,pmask)

    ! this makes sure the boundary conditions are properly defined everywhere
    do n = 1,nlevs
       call bc_tower_level_build(the_bc_tower,n,mla%la(n))
    end do

    ! Build the level 1 data again.
    call multifab_build(  uold(1), mla%la(1),    dm, ng_s)
    call multifab_build(  sold(1), mla%la(1), nscal, ng_s)
    call multifab_build( gpres(1), mla%la(1),    dm, 1)
    call multifab_build(  pres(1), mla%la(1),     1, 1, nodal)
    call multifab_build(  dSdt(1), mla%la(1),     1, 0)
    call multifab_build(   src(1), mla%la(1),     1, 1)

    ! Copy the level 1 data from the "old" temporaries again.
    call multifab_copy_c(  uold(1),1,  uold_temp(1) ,1,   dm)
    call multifab_copy_c(  sold(1),1,  sold_temp(1) ,1,nscal)
    call multifab_copy_c( gpres(1),1, gpres_temp(1),1,    dm)
    call multifab_copy_c(  pres(1),1,  pres_temp(1) ,1,    1)
    call multifab_copy_c(  dSdt(1),1,  dSdt_temp(1), 1,    1)
    call multifab_copy_c(   src(1),1,   src_temp(1), 1,    1)

    nlevs = mla%nlevel
    nlevs_radial = merge(1, nlevs, spherical .eq. 1)

    ! Now make the data for the final time.
    do nl = 1,nlevs-1

       call multifab_build(  uold(nl+1), mla%la(nl+1),    dm, ng_s)
       call multifab_build(  sold(nl+1), mla%la(nl+1), nscal, ng_s)
       call multifab_build( gpres(nl+1), mla%la(nl+1),    dm, 1)
       call multifab_build(  pres(nl+1), mla%la(nl+1),     1, 1, nodal)
       call multifab_build(  dSdt(nl+1), mla%la(nl+1),     1, 0)
       call multifab_build(   src(nl+1), mla%la(nl+1),     1, 1)
       
       ! Define bc_tower at level nl+1.
       call bc_tower_level_build(the_bc_tower,nl+1,mla%la(nl+1))

       ! Fill the data in the new level nl+1 state -- first from the coarser data.
       call fillpatch(uold(nl+1),uold(nl), &
                      ng_s,mba%rr(nl,:), &
                      the_bc_tower%bc_tower_array(nl  ), &
                      the_bc_tower%bc_tower_array(nl+1), &
                      1,1,1,dm)
       call fillpatch(sold(nl+1),sold(nl), &
                      ng_s,mba%rr(nl,:), &
                      the_bc_tower%bc_tower_array(nl  ), &
                      the_bc_tower%bc_tower_array(nl+1), &
                      1,1,dm+rho_comp,nscal)
       do d=1,dm
          call fillpatch(gpres(nl+1),gpres(nl), &
                         1,mba%rr(nl,:), &
                         the_bc_tower%bc_tower_array(nl  ), &
                         the_bc_tower%bc_tower_array(nl+1), &
                         d,d,foextrap_comp,1)
       end do
       call fillpatch(dSdt(nl+1),dSdt(nl), &
                      0,mba%rr(nl,:), &
                      the_bc_tower%bc_tower_array(nl  ), &
                      the_bc_tower%bc_tower_array(nl+1), &
                      1,1,foextrap_comp,1)
       call fillpatch(src(nl+1),src(nl), &
                      1,mba%rr(nl,:), &
                      the_bc_tower%bc_tower_array(nl  ), &
                      the_bc_tower%bc_tower_array(nl+1), &
                      1,1,foextrap_comp,1)

       ! We interpolate p differently because it is nodal, not cell-centered
       call ml_prolongation(pres(nl+1),pres(nl),layout_get_pd(mla%la(nl+1)),mba%rr(nl,:))

       ! Copy from old data at current level, if it exists
       if (mla_old%nlevel .ge. nl+1) then
          call multifab_copy_c(  uold(nl+1),1,  uold_temp(nl+1),1,   dm)
          call multifab_copy_c(  sold(nl+1),1,  sold_temp(nl+1),1,nscal)
          call multifab_copy_c( gpres(nl+1),1, gpres_temp(nl+1),1,   dm)
          call multifab_copy_c(  pres(nl+1),1,  pres_temp(nl+1),1,    1)
          call multifab_copy_c(  dSdt(nl+1),1,  dSdt_temp(nl+1),1,    1)
          call multifab_copy_c(   src(nl+1),1,   src_temp(nl+1),1,    1)
       end if

       call destroy(  uold_temp(nl+1))
       call destroy(  sold_temp(nl+1))
       call destroy( gpres_temp(nl+1))
       call destroy(  pres_temp(nl+1))
       call destroy(  dSdt_temp(nl+1))
       call destroy(   src_temp(nl+1))

    end do

    if (spherical .eq. 1) then

       ! convert rho' -> rho in sold
       call put_in_pert_form(mla,sold,rho0,dx,rho_comp,dm+rho_comp,.false., &
                             the_bc_tower%bc_tower_array)

       ! convert X --> (rho X) in sold 
       call convert_rhoX_to_X(sold,.false.,mla,the_bc_tower%bc_tower_array)

       ! convert (rho h)' -> (rho h) in sold
       call put_in_pert_form(mla,sold,rhoh0,dx,rhoh_comp,dm+rhoh_comp,.false., &
                             the_bc_tower%bc_tower_array)

    end if

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(uold(nlevs))
       call multifab_fill_boundary(sold(nlevs))
       call multifab_fill_boundary(gpres(nlevs))
       call multifab_fill_boundary(pres(nlevs))
       call multifab_fill_boundary(src(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(uold(nlevs),1,1,dm,the_bc_tower%bc_tower_array(nlevs))
       call multifab_physbc(sold(nlevs),1,dm+rho_comp,nscal, &
                            the_bc_tower%bc_tower_array(nlevs))
       call multifab_physbc(pres(nlevs),1,foextrap_comp,1,the_bc_tower%bc_tower_array(nlevs))
       do d=1,dm
          call multifab_physbc(gpres(nlevs),d,foextrap_comp,1, &
                               the_bc_tower%bc_tower_array(nlevs))
       end do
       call multifab_physbc(src(nlevs),1,foextrap_comp,1,the_bc_tower%bc_tower_array(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(uold(n-1),uold(n),mla%mba%rr(n-1,:))
          call ml_cc_restriction(sold(n-1),sold(n),mla%mba%rr(n-1,:))
          call ml_cc_restriction(gpres(n-1),gpres(n),mla%mba%rr(n-1,:))
          call ml_cc_restriction(src(n-1),src(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(uold(n),uold(n-1),ng_s,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n),1,1,dm, &
                                         fill_crse_input=.false.)
          call multifab_fill_ghost_cells(sold(n),sold(n-1),ng_s,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n), &
                                         rho_comp,dm+rho_comp,nscal,fill_crse_input=.false.)
          do d=1,dm
             call multifab_fill_ghost_cells(gpres(n),gpres(n-1),1,mla%mba%rr(n-1,:), &
                                            the_bc_tower%bc_tower_array(n-1), &
                                            the_bc_tower%bc_tower_array(n), &
                                            d,foextrap_comp,1,fill_crse_input=.false.)
          end do
          call multifab_fill_ghost_cells(src(n),src(n-1),1,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n),1,foextrap_comp,1, &
                                         fill_crse_input=.false.)

       enddo

    end if

    call destroy(mba)

    call destroy(  uold_temp(1))
    call destroy(  sold_temp(1))
    call destroy( gpres_temp(1))
    call destroy(  pres_temp(1))
    call destroy(  dSdt_temp(1))
    call destroy(   src_temp(1))

    call destroy(mla_old)

  end subroutine regrid

end module regrid_module
