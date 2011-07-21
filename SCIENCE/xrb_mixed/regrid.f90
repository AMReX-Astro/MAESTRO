module regrid_module

  use BoxLib
  use ml_boxarray_module
  use ml_layout_module
  use define_bc_module
  use restart_module
  use box_util_module

  implicit none

  private

  public :: regrid

contains

  subroutine regrid(nstep,mla,uold,sold,gpi,pi,dSdt,src,dx,the_bc_tower,rho0,rhoh0,is_restart,rhoHdot)

    use fillpatch_module
    use ml_prolongation_module
    use multifab_physbc_module
    use multifab_fill_ghost_module
    use ml_restriction_module
    use make_new_grids_module
    use convert_rhoX_to_X_module
    use pert_form_module

    use probin_module, only : verbose, nodal, pmask, &
         regrid_int, amr_buf_width, &
         max_grid_size_2, max_grid_size_3, ref_ratio, max_levs, &
         ppm_type, bds_type, dump_grid_file
    use geometry, only: nlevs_radial, spherical
    use variables, only: nscal, rho_comp, rhoh_comp, foextrap_comp
    use network, only: nspec

    integer       ,  intent(in   ) :: nstep
    type(ml_layout), intent(inout) :: mla
    type(multifab),  pointer       :: uold(:),sold(:),gpi(:),pi(:)
    type(multifab),  pointer       :: dSdt(:),src(:)
    type(multifab),  pointer       :: rhoHdot(:)
    real(dp_t)    ,  pointer       :: dx(:,:)
    type(bc_tower),  intent(inout) :: the_bc_tower
    real(kind=dp_t), intent(in   ) :: rho0(:,0:),rhoh0(:,0:)
    logical        , intent(in   ) :: is_restart

    ! local
    logical           :: new_grid
    integer           :: n, nl, d, ng_s, dm, nlevs, ng_buffer
    type(layout)      :: la_array(max_levs)
    type(ml_layout)   :: mla_old
    type(ml_boxarray) :: mba
    integer           :: un
    logical           :: lexist

    ! These are copies to hold the old data.
    type(multifab) :: uold_temp(max_levs), sold_temp(max_levs), gpi_temp(max_levs)
    type(multifab) :: pi_temp(max_levs), dSdt_temp(max_levs), src_temp(max_levs)
    type(multifab) :: rhoHdot_temp(max_levs)

    dm    = mla%dim
    nlevs = mla%nlevel

    if (verbose .ge. 1) then
       if (parallel_IOProcessor()) print*,'Calling regrid'
    end if

    if (ppm_type .eq. 2 .or. bds_type .eq. 1) then
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
       call multifab_build(   gpi_temp(n),mla_old%la(n),   dm, 1)
       call multifab_build(    pi_temp(n),mla_old%la(n),    1, 1, nodal)
       call multifab_build(  dSdt_temp(n),mla_old%la(n),    1, 0)
       call multifab_build(   src_temp(n),mla_old%la(n),    1, 1)
       call multifab_build(rhoHdot_temp(n),mla_old%la(n),   1, 0)

       call multifab_copy_c(  uold_temp(n),1,  uold(n),1,   dm)
       call multifab_copy_c(  sold_temp(n),1,  sold(n),1,nscal)
       call multifab_copy_c(   gpi_temp(n),1,   gpi(n),1,   dm)
       call multifab_copy_c(    pi_temp(n),1,    pi(n),1,    1)
       call multifab_copy_c(  dSdt_temp(n),1,  dSdt(n),1,    1)
       call multifab_copy_c(   src_temp(n),1,   src(n),1,    1)
       call multifab_copy_c(rhoHdot_temp(n),1,rhoHdot(n),1,  1)

       ! Get rid of the old data structures so we can create new ones 
       ! with the same names.
       call multifab_destroy(  uold(n))
       call multifab_destroy(  sold(n))
       call multifab_destroy(   gpi(n))
       call multifab_destroy(    pi(n))
       call multifab_destroy(  dSdt(n))
       call multifab_destroy(   src(n))
       call multifab_destroy(rhoHdot(n))

    end do

    call destroy(mla)

    ! mba is big enough to hold max_levs levels
    ! even though we know we had nlevs last time, we might 
    ! want more or fewer levels after regrid (if nlevs < max_levs)
    call ml_boxarray_build_n(mba,max_levs,dm)

    do n = 1, max_levs-1
       mba%rr(n,:) = ref_ratio
    enddo

    if (associated(uold)) then
       deallocate(uold,sold,pi,gpi,dSdt,src)
       deallocate(rhoHdot)
    end if

    allocate(uold(max_levs),sold(max_levs),pi(max_levs),gpi(max_levs))
    allocate(dSdt(max_levs),src(max_levs))
    allocate(rhoHdot(max_levs))

    ! Copy the level 1 boxarray
    call copy(mba%bas(1),mla_old%mba%bas(1))

    ! Copy the pd(:)
    mba%pd(1) = mla_old%mba%pd(1)
    do n = 2, max_levs
       mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
    enddo

    ! Build the level 1 layout.
    call layout_build_ba(la_array(1),mba%bas(1),mba%pd(1),pmask)

    ! Build and fill the level 1 data only.
    call build_and_fill_data(1,la_array(1),mla_old, &
                             uold     ,sold     ,gpi     ,pi     ,dSdt     ,src,      rhoHdot, &
                             uold_temp,sold_temp,gpi_temp,pi_temp,dSdt_temp,src_temp, rhoHdot_temp, &
                             the_bc_tower,dm,ng_s,mba%rr(1,:))

    nl       = 1
    new_grid = .true.

    ng_buffer = 2

    do while ( (nl .lt. max_levs) .and. (new_grid) )

       ! Do we need finer grids?

       ! Need to fill ghost cells here in case we use them in tagging
       call multifab_fill_boundary(sold(nl))
       call multifab_physbc(sold(nl),rho_comp,dm+rho_comp,nscal, &
                            the_bc_tower%bc_tower_array(nl))

       if (nl .eq. 1) then
          call make_new_grids(new_grid,la_array(nl),la_array(nl+1),sold(nl),dx(nl,1), &
                              amr_buf_width,ref_ratio,nl,max_grid_size_2,rhoHdot(nl))
       else
          call make_new_grids(new_grid,la_array(nl),la_array(nl+1),sold(nl),dx(nl,1), &
                              amr_buf_width,ref_ratio,nl,max_grid_size_3,rhoHdot(nl))
       end if

       if (new_grid) then

          call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))

          ! Need to enforce proper nesting within the grid creation procedure so that we can
          !  fillpatch the new levels.
          if (nl .ge. 2) then

             ! Test on whether grids are already properly nested
             if (.not. ml_boxarray_properly_nested(mba, ng_buffer, pmask, 2, nl+1)) then

                call enforce_proper_nesting(mba,la_array,max_grid_size_2,max_grid_size_3)

                ! Loop over all the lower levels which we might have changed when we enforced proper nesting.
                do n = 2,nl
   
                   ! This makes sure the boundary conditions are properly defined everywhere
                   call bc_tower_level_build(the_bc_tower,n,la_array(n))
   
                   ! Delete old multifabs so that we can rebuild them.
                   call destroy(  sold(n))
                   call destroy(  uold(n))
                   call destroy(   gpi(n))
                   call destroy(    pi(n))
                   call destroy(  dSdt(n))
                   call destroy(   src(n))
                   call destroy(rhoHdot(n))
   
                   ! Rebuild the lower level data again if it changed.
                   call build_and_fill_data(n,la_array(n),mla_old, &
                                            uold     ,sold     ,gpi     ,pi     ,dSdt     ,src,      rhoHdot, &
                                            uold_temp,sold_temp,gpi_temp,pi_temp,dSdt_temp,src_temp, rhoHdot_temp, &
                                            the_bc_tower,dm,ng_s,mba%rr(n-1,:))
                end do
             end if
          end if

          ! Define bc_tower at level nl+1.
          call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))

          ! Build the level nl+1 data only.
          call build_and_fill_data(nl+1,la_array(nl+1),mla_old, &
                                   uold     ,sold     ,gpi     ,pi     ,dSdt     ,src,      rhoHdot, &
                                   uold_temp,sold_temp,gpi_temp,pi_temp,dSdt_temp,src_temp, rhoHdot_temp, &
                                   the_bc_tower,dm,ng_s,mba%rr(nl,:))

          nlevs = nl+1
          nlevs_radial = merge(1, nlevs, spherical .eq. 1)
          nl = nl+1

       endif

    enddo

    nlevs = nl
    nlevs_radial = merge(1, nlevs, spherical .eq. 1)

    ! Note: this build actually sets mla%la(n) = la_array(n) so we mustn't delete la_array(n)
    call build(mla,mba,la_array,pmask,nlevs)

    ! This makes sure the boundary conditions are properly defined everywhere
    do n = 1, nlevs
       call bc_tower_level_build(the_bc_tower,n,la_array(n))
    end do

    if (spherical .eq. 1) then

       if (is_restart) nlevs = nlevs-1

       ! convert (rho X) --> X in sold_temp 
       call convert_rhoX_to_X(sold_temp,.true.,mla_old,the_bc_tower%bc_tower_array)

       ! convert rho -> rho' in sold_temp
       call put_in_pert_form(mla_old,sold_temp,rho0,dx,rho_comp,foextrap_comp,.true., &
                             the_bc_tower%bc_tower_array)

       ! convert (rho h) -> (rho h)' in sold_temp
       call put_in_pert_form(mla_old,sold_temp,rhoh0,dx,rhoh_comp,foextrap_comp,.true., &
                             the_bc_tower%bc_tower_array)

       ! Delete old multifabs so that we can rebuild them.
       do n = 1, nlevs
          call destroy(  sold(n))
          call destroy(  uold(n))
          call destroy(   gpi(n))
          call destroy(    pi(n))
          call destroy(  dSdt(n))
          call destroy(   src(n))
          call destroy(rhoHdot(n))
       end do

       ! Rebuild using the perturbational form
       do n = 1, nlevs
          call build_and_fill_data(n,la_array(n),mla_old, &
                                   uold     ,sold     ,gpi     ,pi     ,dSdt     ,src,      rhoHdot, &
                                   uold_temp,sold_temp,gpi_temp,pi_temp,dSdt_temp,src_temp, rhoHdot_temp, &
                                   the_bc_tower,dm,ng_s,mba%rr(n-1,:))
       end do

       if (is_restart) nlevs = nlevs+1

       ! convert rho' -> rho in sold
       call put_in_pert_form(mla,sold,rho0,dx,rho_comp,dm+rho_comp,.false., &
                             the_bc_tower%bc_tower_array)

       ! convert (rho h)' -> (rho h) in sold
       call put_in_pert_form(mla,sold,rhoh0,dx,rhoh_comp,dm+rhoh_comp,.false., &
                             the_bc_tower%bc_tower_array)

       ! convert X --> (rho X) in sold 
       call convert_rhoX_to_X(sold,.false.,mla,the_bc_tower%bc_tower_array)

    end if

    do nl = 1,nlevs
       if (mla_old%nlevel .ge. nl) then
          call destroy(  uold_temp(nl))
          call destroy(  sold_temp(nl))
          call destroy(   gpi_temp(nl))
          call destroy(    pi_temp(nl))
          call destroy(  dSdt_temp(nl))
          call destroy(   src_temp(nl))
          call destroy(rhoHdot_temp(nl))
       end if
    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(uold(nlevs))
       call multifab_fill_boundary(sold(nlevs))
       call multifab_fill_boundary(gpi(nlevs))
       call multifab_fill_boundary(pi(nlevs))
       call multifab_fill_boundary(src(nlevs))
       call multifab_fill_boundary(rhoHdot(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(uold(nlevs),1,1,dm,the_bc_tower%bc_tower_array(nlevs))
       call multifab_physbc(sold(nlevs),1,dm+rho_comp,nscal, &
                            the_bc_tower%bc_tower_array(nlevs))
       call multifab_physbc(pi(nlevs),1,foextrap_comp,1,the_bc_tower%bc_tower_array(nlevs))
       do d=1,dm
          call multifab_physbc(gpi(nlevs),d,foextrap_comp,1, &
                               the_bc_tower%bc_tower_array(nlevs))
       end do
       call multifab_physbc(src(nlevs),1,foextrap_comp,1,the_bc_tower%bc_tower_array(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(uold(n-1),uold(n),mla%mba%rr(n-1,:))
          call ml_cc_restriction(sold(n-1),sold(n),mla%mba%rr(n-1,:))
          call ml_cc_restriction(gpi(n-1),gpi(n),mla%mba%rr(n-1,:))
          call ml_cc_restriction(src(n-1),src(n),mla%mba%rr(n-1,:))
          call ml_cc_restriction(rhoHdot(n-1),rhoHdot(n),mla%mba%rr(n-1,:))

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
             call multifab_fill_ghost_cells(gpi(n),gpi(n-1),1,mla%mba%rr(n-1,:), &
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

    ! optionally output details of the grid structure
    if (dump_grid_file .and. parallel_IOProcessor()) then

       un = unit_new()
       inquire(file="grids.out", exist=lexist)

       if (lexist .and. (nstep /= 1)) then
          open(unit=un, file="grids.out", status="old", position="append")
       else
          open(unit=un, file="grids.out", status="unknown")
       endif

998    format(1x,"step = ", i7)
999    format(1x,"level",4x,"        # valid cells",4x,"# valid + ghost cells",4x,"       # nboxes")
1000   format(1x,i5     ,4x, 11x, i10,              4x, 11x, i10, 4x, 5x, i10)

       
       write (un, 998) nstep
       write (un, 999) 
       do n = 1, nlevs
          write (un,1000) n, multifab_volume(sold(n),.false.), multifab_volume(sold(n),.true.), nboxes(sold(n))
       enddo
       write (un, *) " "
       close (un)

    endif

!   if (verbose .and. parallel_IOProcessor()) then
!       print *,'New grids after regridding:'
!       do n = 2, nlevs
!          print *,'...at level ',n
!          call print(mba%bas(n))
!       end do
!   end if

    call destroy(mba)

    call destroy(mla_old)

  end subroutine regrid

  subroutine build_and_fill_data(lev,la,mla_old, &
                                 uold     ,sold     ,gpi     ,pi     ,dSdt     ,src,      rhoHdot, &
                                 uold_temp,sold_temp,gpi_temp,pi_temp,dSdt_temp,src_temp, rhoHdot_temp, &
                                 the_bc_tower,dm,ng_s,rr)

    use fillpatch_module
    use ml_prolongation_module
    use probin_module, only : nodal
    use variables    , only : nscal, rho_comp, rhoh_comp, foextrap_comp

    integer                    , intent(in   ) :: lev, dm, ng_s, rr(:)
    type(layout)               , intent(in   ) :: la
    type(ml_layout)            , intent(in   ) :: mla_old
    type(bc_tower)             , intent(inout) :: the_bc_tower
    type(multifab)   ,  pointer :: uold(:),sold(:),gpi(:),pi(:),dSdt(:),src(:)
    type(multifab)   ,  pointer :: rhoHdot(:)
    type(multifab)             , intent(in   ) :: uold_temp(:),sold_temp(:),gpi_temp(:),pi_temp(:)
    type(multifab)             , intent(in   ) :: dSdt_temp(:),src_temp(:), rhoHdot_temp(:)
 
    integer :: d 

    ! Build the level lev data only.
    call multifab_build(  uold(lev), la,    dm, ng_s)
    call multifab_build(  sold(lev), la, nscal, ng_s)
    call multifab_build(   gpi(lev), la,    dm, 1)
    call multifab_build(    pi(lev), la,     1, 1, nodal)
    call multifab_build(  dSdt(lev), la,     1, 0)
    call multifab_build(   src(lev), la,     1, 1)
    call multifab_build(rhoHdot(lev),la,     1, 0)

    ! Fill the data in the new level lev state -- first from the coarser data if lev > 1.

    if (lev .gt. 1) then

       call fillpatch(uold(lev),uold(lev-1), &
                      ng_s,rr, &
                      the_bc_tower%bc_tower_array(lev-1), &
                      the_bc_tower%bc_tower_array(lev  ), &
                      1,1,1,dm)
       call fillpatch(sold(lev),sold(lev-1), &
                      ng_s,rr, &
                      the_bc_tower%bc_tower_array(lev-1), &
                      the_bc_tower%bc_tower_array(lev  ), &
                      1,1,dm+rho_comp,nscal)
       do d=1,dm
          call fillpatch(gpi(lev),gpi(lev-1), &
                         1,rr, &
                         the_bc_tower%bc_tower_array(lev-1), &
                         the_bc_tower%bc_tower_array(lev  ), &
                         d,d,foextrap_comp,1)
       end do
       call fillpatch(dSdt(lev),dSdt(lev-1), &
                      0,rr, &
                      the_bc_tower%bc_tower_array(lev-1), &
                      the_bc_tower%bc_tower_array(lev  ), &
                      1,1,foextrap_comp,1)
       call fillpatch(src(lev),src(lev-1), &
                      1,rr, &
                      the_bc_tower%bc_tower_array(lev-1), &
                      the_bc_tower%bc_tower_array(lev  ), &
                      1,1,foextrap_comp,1) 
       call fillpatch(rhoHdot(lev),rhoHdot(lev-1), &
                      0,rr, &
                      the_bc_tower%bc_tower_array(lev-1), &
                      the_bc_tower%bc_tower_array(lev  ), &
                      1,1,foextrap_comp,1) 
       ! We interpolate p differently because it is nodal, not cell-centered
       call ml_nodal_prolongation(pi(lev), pi(lev-1), rr)

    end if

    ! Copy from old data at current level, if it exists
    if (mla_old%nlevel .ge. lev) then
       call multifab_copy_c(  uold(lev),1,  uold_temp(lev),1,   dm)
       call multifab_copy_c(  sold(lev),1,  sold_temp(lev),1,nscal)
       call multifab_copy_c(   gpi(lev),1,   gpi_temp(lev),1,   dm)
       call multifab_copy_c(    pi(lev),1,    pi_temp(lev),1,    1)
       call multifab_copy_c(  dSdt(lev),1,  dSdt_temp(lev),1,    1)
       call multifab_copy_c(   src(lev),1,   src_temp(lev),1,    1)
       call multifab_copy_c(rhoHdot(lev),1,rhoHdot_temp(lev),1,  1)
    end if

  end subroutine build_and_fill_data

end module regrid_module
