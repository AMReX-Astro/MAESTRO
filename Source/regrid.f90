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

  subroutine regrid(mla,u,s,gp,p,dx,the_bc_tower)

    use probin_module, only : nlevs, nodal, pmask, regrid_int, max_grid_size, ref_ratio, &
         max_levs
    use geometry, only: dm
    use variables, only: nscal

    type(ml_layout),intent(inout) :: mla
    type(multifab), pointer       :: u(:),s(:),gp(:),p(:)
    real(dp_t)    , pointer       :: dx(:,:)
    type(bc_tower), intent(inout) :: the_bc_tower

    logical           :: new_grid
    integer           :: n, nl, buf_wid
    type(layout)      :: la_array(max_levs)
    type(ml_layout)   :: mla_old
    type(ml_boxarray) :: mba

    ! These are copies to hold the old data.
    type(multifab) :: uold(nlevs), sold(nlevs), gpold(nlevs), pold(nlevs)

    if (max_levs < 2) &
         call bl_error('Dont call regrid with max_levs < 2')

    call ml_layout_build(mla_old,mla%mba,mla%pmask)

    do n = 1,nlevs
       !
       ! Create copies of the old data.
       call multifab_build( uold(n),mla_old%la(n),   dm, 3)
       call multifab_build( sold(n),mla_old%la(n),nscal, 3)
       call multifab_build(gpold(n),mla_old%la(n),   dm, 1)
       call multifab_build( pold(n),mla_old%la(n),    1, 1, nodal)
       call multifab_copy_c( uold(n),1, u(n),1,   dm)
       call multifab_copy_c( sold(n),1, s(n),1,nscal)
       call multifab_copy_c(gpold(n),1,gp(n),1,   dm)
       call multifab_copy_c( pold(n),1, p(n),1,    1)
       !
       ! Get rid of the old data structures so we can create new ones 
       ! with the same names.
       call multifab_destroy( u(n))
       call multifab_destroy( s(n))
       call multifab_destroy(gp(n))
       call multifab_destroy( p(n))

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

    if (associated(u)) then
       deallocate(u,s,p,gp)
    end if

    allocate(u(max_levs),s(max_levs),p(max_levs),gp(max_levs))

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
    ! call make_new_state(la_array(1),u(1),s(1),gp(1),p(1)) 

    ! Copy the level 1 data from the "old" temporaries.
    call multifab_copy_c( u(1),1, uold(1) ,1,   dm)
    call multifab_copy_c( s(1),1, sold(1) ,1,nscal)
    call multifab_copy_c(gp(1),1,gpold(1),1,    dm)
    call multifab_copy_c( p(1),1, pold(1) ,1,    1)

    nl       = 1
    new_grid = .true.

    do while ( (nl .lt. max_levs) .and. (new_grid) )

       ! Do we need finer grids?

       call make_new_grids(new_grid,la_array(nl),la_array(nl+1),s(nl),dx(nl,1),buf_wid,&
            ref_ratio,nl,max_grid_size)

       if (new_grid) then

          call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))

          ! Build the level nl+1 data only.
          ! call make_new_state(la_array(nl+1),u(nl+1),s(nl+1),gp(nl+1),p(nl+1)) 

          ! Define bc_tower at level nl+1.
          call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))

          ! Fill the data in the new level nl+1 state -- first from the coarser data.
          call fillpatch(u(nl+1),u(nl), &
               3,mba%rr(nl,:), &
               the_bc_tower%bc_tower_array(nl  ), &
               the_bc_tower%bc_tower_array(nl+1), &
               1,1,1,dm)
          call fillpatch(s(nl+1),s(nl), &
               3,mba%rr(nl,:), &
               the_bc_tower%bc_tower_array(nl  ), &
               the_bc_tower%bc_tower_array(nl+1), &
               1,1,dm+1,nscal)
          call fillpatch(gp(nl+1),gp(nl), &
               1,mba%rr(nl,:), &
               the_bc_tower%bc_tower_array(nl  ), &
               the_bc_tower%bc_tower_array(nl+1), &
               1,1,1,dm)

          ! We interpolate p differently because it is nodal, not cell-centered
          call ml_prolongation(p(nl+1),p(nl),layout_get_pd(la_array(nl+1)),mba%rr(nl,:))

          ! Copy from old data at current level, if it exists
          if (mla_old%nlevel .ge. nl+1) then
             call multifab_copy_c( u(nl+1),1, uold(nl+1),1,   dm)
             call multifab_copy_c( s(nl+1),1, sold(nl+1),1,nscal)
             call multifab_copy_c(gp(nl+1),1,gpold(nl+1),1,   dm)
             call multifab_copy_c( p(nl+1),1, pold(nl+1),1,    1)
          end if

          nlevs = nl+1
          nl = nl + 1

       endif

    enddo

    do n = 1,nl
       call destroy( s(n))
       call destroy( u(n))
       call destroy(gp(n))
       call destroy( p(n))
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
    ! call make_new_state(mla%la(1),u(1),s(1),gp(1),p(1)) 

    ! Copy the level 1 data from the "old" temporaries again.
    call multifab_copy_c( u(1),1, uold(1) ,1,   dm)
    call multifab_copy_c( s(1),1, sold(1) ,1,nscal)
    call multifab_copy_c(gp(1),1,gpold(1),1,    dm)
    call multifab_copy_c( p(1),1, pold(1) ,1,    1)

    nlevs = mla%nlevel

    ! Now make the data for the final time.
    do nl = 1,nlevs-1

       ! call make_new_state(mla%la(nl+1),u(nl+1),s(nl+1),gp(nl+1),p(nl+1))

       ! Define bc_tower at level nl+1.
       call bc_tower_level_build(the_bc_tower,nl+1,mla%la(nl+1))

       ! Fill the data in the new level nl+1 state -- first from the coarser data.
       call fillpatch(u(nl+1),u(nl), &
            3,mba%rr(nl,:), &
            the_bc_tower%bc_tower_array(nl  ), &
            the_bc_tower%bc_tower_array(nl+1), &
            1,1,1,dm)
       call fillpatch(s(nl+1),s(nl), &
            3,mba%rr(nl,:), &
            the_bc_tower%bc_tower_array(nl  ), &
            the_bc_tower%bc_tower_array(nl+1), &
            1,1,dm+1,nscal)
       call fillpatch(gp(nl+1),gp(nl), &
            1,mba%rr(nl,:), &
            the_bc_tower%bc_tower_array(nl  ), &
            the_bc_tower%bc_tower_array(nl+1), &
            1,1,1,dm)

       ! We interpolate p differently because it is nodal, not cell-centered
       call ml_prolongation(p(nl+1),p(nl),layout_get_pd(mla%la(nl+1)),mba%rr(nl,:))

       ! Copy from old data at current level, if it exists
       if (mla_old%nlevel .ge. nl+1) then
          call multifab_copy_c( u(nl+1),1, uold(nl+1),1,   dm)
          call multifab_copy_c( s(nl+1),1, sold(nl+1),1,nscal)
          call multifab_copy_c(gp(nl+1),1,gpold(nl+1),1,   dm)
          call multifab_copy_c( p(nl+1),1, pold(nl+1),1,    1)
       end if

       call destroy( uold(nl+1))
       call destroy( sold(nl+1))
       call destroy(gpold(nl+1))
       call destroy( pold(nl+1))

    end do

    call destroy(mba)

    call destroy( uold(1))
    call destroy( sold(1))
    call destroy(gpold(1))
    call destroy( pold(1))

    call destroy(mla_old)

  end subroutine regrid

end module regrid_module
