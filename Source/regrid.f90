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

  implicit none

  private

  public :: regrid

contains

  subroutine regrid(mla,uold,sold,gpres,pres,dSdt,Source_old,rho_omegadot2,dx,the_bc_tower)

    use probin_module, only : nlevs, nodal, pmask, regrid_int, max_grid_size, ref_ratio, &
         max_levs
    use geometry, only: dm
    use variables, only: nscal

    type(ml_layout),intent(inout) :: mla
    type(multifab), pointer       :: uold(:),sold(:),gpres(:),pres(:)
    type(multifab), pointer       :: dSdt,Source_old(:),rho_omegadot2(:)
    real(dp_t)    , pointer       :: dx(:,:)
    type(bc_tower), intent(inout) :: the_bc_tower

    logical           :: new_grid
    integer           :: n, nl, buf_wid
    type(layout)      :: la_array(max_levs)
    type(ml_layout)   :: mla_old
    type(ml_boxarray) :: mba

    ! These are copies to hold the old data.
    type(multifab) :: uold_temp(nlevs), sold_temp(nlevs), gpres_temp(nlevs), pres_temp(nlevs)
    type(multifab) :: dSdt_temp(nlevs), srcold_temp(nlevs), rw2_temp(nlevs)

    if (max_levs < 2) &
         call bl_error('Dont call regrid with max_levs < 2')

    call ml_layout_build(mla_old,mla%mba,mla%pmask)

    do n = 1,nlevs

       ! Create copies of the old data.
       call multifab_build( uold_temp(n),mla_old%la(n),   dm, 3)
       call multifab_build( sold_temp(n),mla_old%la(n),nscal, 3)
       call multifab_build(gpres_temp(n),mla_old%la(n),   dm, 1)
       call multifab_build( pres_temp(n),mla_old%la(n),    1, 1, nodal)

       call multifab_copy_c( uold_temp(n),1, uold(n),1,   dm)
       call multifab_copy_c( sold_temp(n),1, sold(n),1,nscal)
       call multifab_copy_c(gpres_temp(n),1,gpres(n),1,   dm)
       call multifab_copy_c( pres_temp(n),1, pres(n),1,    1)

       ! Get rid of the old data structures so we can create new ones 
       ! with the same names.
       call multifab_destroy( uold(n))
       call multifab_destroy( sold(n))
       call multifab_destroy(gpres(n))
       call multifab_destroy( pres(n))

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
       deallocate(uold,sold,pres,gpres)
    end if

    allocate(uold(max_levs),sold(max_levs),pres(max_levs),gpres(max_levs))

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
    call multifab_build( uold(1), la_array(1),    dm, 3)
    call multifab_build( sold(1), la_array(1), nscal, 3)
    call multifab_build(gpres(1), la_array(1),    dm, 1)
    call multifab_build( pres(1), la_array(1),     1, 1, nodal)

    ! Copy the level 1 data from the "old" temporaries.
    call multifab_copy_c( uold(1),1, uold_temp(1) ,1,   dm)
    call multifab_copy_c( sold(1),1, sold_temp(1) ,1,nscal)
    call multifab_copy_c(gpres(1),1,gpres_temp(1), 1,   dm)
    call multifab_copy_c( pres(1),1, pres_temp(1) ,1,    1)

    nl       = 1
    new_grid = .true.

    do while ( (nl .lt. max_levs) .and. (new_grid) )

       ! Do we need finer grids?

       call make_new_grids(new_grid,la_array(nl),la_array(nl+1),sold(nl),dx(nl,1),buf_wid,&
                           ref_ratio,nl,max_grid_size)

       if (new_grid) then

          call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))

          ! Build the level nl+1 data only.
          call multifab_build( uold(nl+1), la_array(nl+1),    dm, 3)
          call multifab_build( sold(nl+1), la_array(nl+1), nscal, 3)
          call multifab_build(gpres(nl+1), la_array(nl+1),    dm, 1)
          call multifab_build( pres(nl+1), la_array(nl+1),     1, 1, nodal)

          ! Define bc_tower at level nl+1.
          call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))

          ! Fill the data in the new level nl+1 state -- first from the coarser data.
          call fillpatch(uold(nl+1),uold(nl), &
                         3,mba%rr(nl,:), &
                         the_bc_tower%bc_tower_array(nl  ), &
                         the_bc_tower%bc_tower_array(nl+1), &
                         1,1,1,dm)
          call fillpatch(sold(nl+1),sold(nl), &
                         3,mba%rr(nl,:), &
                         the_bc_tower%bc_tower_array(nl  ), &
                         the_bc_tower%bc_tower_array(nl+1), &
                         1,1,dm+1,nscal)
          call fillpatch(gpres(nl+1),gpres(nl), &
                         1,mba%rr(nl,:), &
                         the_bc_tower%bc_tower_array(nl  ), &
                         the_bc_tower%bc_tower_array(nl+1), &
                         1,1,1,dm)

          ! We interpolate p differently because it is nodal, not cell-centered
          call ml_prolongation(pres(nl+1),pres(nl),layout_get_pd(la_array(nl+1)), &
                               mba%rr(nl,:))

          ! Copy from old data at current level, if it exists
          if (mla_old%nlevel .ge. nl+1) then
             call multifab_copy_c( uold(nl+1),1, uold_temp(nl+1),1,   dm)
             call multifab_copy_c( sold(nl+1),1, sold_temp(nl+1),1,nscal)
             call multifab_copy_c(gpres(nl+1),1,gpres_temp(nl+1),1,   dm)
             call multifab_copy_c( pres(nl+1),1, pres_temp(nl+1),1,    1)
          end if

          nlevs = nl+1
          nl = nl + 1

       endif

    enddo

    do n = 1,nl
       call destroy( sold(n))
       call destroy( uold(n))
       call destroy(gpres(n))
       call destroy( pres(n))
    end do

    nlevs = nl

    call ml_layout_restricted_build(mla,mba,nlevs,pmask)

    ! check for proper nesting
    if (nlevs .ge. 3) &
         call enforce_proper_nesting(mba,la_array,max_grid_size)

    do n = 1,nl
       call destroy(la_array(n))
    end do

    call destroy(mla)

    call ml_layout_restricted_build(mla,mba,nlevs,pmask)

    ! Build the level 1 data again.
    call multifab_build( uold(1), mla%la(1),    dm, 3)
    call multifab_build( sold(1), mla%la(1), nscal, 3)
    call multifab_build(gpres(1), mla%la(1),    dm, 1)
    call multifab_build( pres(1), mla%la(1),     1, 1, nodal)

    ! Copy the level 1 data from the "old" temporaries again.
    call multifab_copy_c( uold(1),1, uold_temp(1) ,1,   dm)
    call multifab_copy_c( sold(1),1, sold_temp(1) ,1,nscal)
    call multifab_copy_c(gpres(1),1,gpres_temp(1),1,    dm)
    call multifab_copy_c( pres(1),1, pres_temp(1) ,1,    1)

    nlevs = mla%nlevel

    ! Now make the data for the final time.
    do nl = 1,nlevs-1

       call multifab_build( uold(nl+1), mla%la(nl+1),    dm, 3)
       call multifab_build( sold(nl+1), mla%la(nl+1), nscal, 3)
       call multifab_build(gpres(nl+1), mla%la(nl+1),    dm, 1)
       call multifab_build( pres(nl+1), mla%la(nl+1),     1, 1, nodal)
       
       ! Define bc_tower at level nl+1.
       call bc_tower_level_build(the_bc_tower,nl+1,mla%la(nl+1))

       ! Fill the data in the new level nl+1 state -- first from the coarser data.
       call fillpatch(uold(nl+1),uold(nl), &
                      3,mba%rr(nl,:), &
                      the_bc_tower%bc_tower_array(nl  ), &
                      the_bc_tower%bc_tower_array(nl+1), &
            1,1,1,dm)
       call fillpatch(sold(nl+1),sold(nl), &
                      3,mba%rr(nl,:), &
                      the_bc_tower%bc_tower_array(nl  ), &
                      the_bc_tower%bc_tower_array(nl+1), &
                      1,1,dm+1,nscal)
       call fillpatch(gpres(nl+1),gpres(nl), &
                      1,mba%rr(nl,:), &
                      the_bc_tower%bc_tower_array(nl  ), &
                      the_bc_tower%bc_tower_array(nl+1), &
                      1,1,1,dm)

       ! We interpolate p differently because it is nodal, not cell-centered
       call ml_prolongation(pres(nl+1),pres(nl),layout_get_pd(mla%la(nl+1)),mba%rr(nl,:))

       ! Copy from old data at current level, if it exists
       if (mla_old%nlevel .ge. nl+1) then
          call multifab_copy_c( uold(nl+1),1, uold_temp(nl+1),1,   dm)
          call multifab_copy_c( sold(nl+1),1, sold_temp(nl+1),1,nscal)
          call multifab_copy_c(gpres(nl+1),1,gpres_temp(nl+1),1,   dm)
          call multifab_copy_c( pres(nl+1),1, pres_temp(nl+1),1,    1)
       end if

       call destroy( uold_temp(nl+1))
       call destroy( sold_temp(nl+1))
       call destroy(gpres_temp(nl+1))
       call destroy( pres_temp(nl+1))

    end do

    call destroy(mba)

    call destroy( uold_temp(1))
    call destroy( sold_temp(1))
    call destroy(gpres_temp(1))
    call destroy( pres_temp(1))

    call destroy(mla_old)

  end subroutine regrid

end module regrid_module
