module regrid_module

  use BoxLib
  use bl_prof_module
  use ml_boxarray_module
  use ml_layout_module
  use define_bc_module
  use box_util_module
  use tag_boxes_module, only : tagging_needs_ghost_cells

  implicit none

  private

  public :: regrid

contains

  subroutine regrid(nstep,mla,uold,sold,gpi,dSdt,S_cc,dx,the_bc_tower, &
                    rho0,rhoh0,init_into_finer,tag_mf)

    use fillpatch_module
    use ml_prolongation_module
    use make_new_grids_module
    use convert_rhoX_to_X_module
    use pert_form_module
    use ml_restrict_fill_module
    use probin_module, only : verbose, nodal, pmask, &
         amr_buf_width, &
         max_grid_size_2, max_grid_size_3, ref_ratio, max_levs, &
         ppm_type, bds_type, dump_grid_file
    use geometry, only: nlevs_radial, spherical
    use variables, only: nscal, rho_comp, rhoh_comp, foextrap_comp
    use network, only: nspec

    integer       ,  intent(in   ) :: nstep
    type(ml_layout), intent(inout) :: mla
    type(multifab),  pointer       :: uold(:),sold(:),gpi(:)
    type(multifab),  pointer       :: dSdt(:),S_cc(:)
    type(multifab),  pointer       :: tag_mf(:)
    real(dp_t)    ,  pointer       :: dx(:,:)
    type(bc_tower),  intent(inout) :: the_bc_tower
    real(kind=dp_t), intent(in   ) :: rho0(:,0:),rhoh0(:,0:)
    logical        , intent(in   ) :: init_into_finer

    ! local
    logical           :: new_grid
    integer           :: n, nl, ng_s, dm, nlevs, ng_buffer
    type(layout)      :: la_array(max_levs)
    type(ml_boxarray) :: mba
    integer           :: un
    logical           :: lexist

    ! These are copies to hold the old data.
    integer         :: nlevs_temp
    type(ml_layout) :: mla_temp
    type(bc_tower)  :: the_bc_tower_temp
    type(multifab)  :: uold_temp(max_levs), sold_temp(max_levs), gpi_temp(max_levs)
    type(multifab)  :: dSdt_temp(max_levs), S_cc_temp(max_levs)
    type(multifab)  :: tag_mf_temp(max_levs)
    type(multifab), allocatable :: uold_opt(:), sold_opt(:), gpi_opt(:), &
         dSdt_opt(:), S_cc_opt(:), tag_mf_opt(:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "regrid")

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

    ! create a copy of the original mla and bc_tower
    ! we need this for putting the original state into perturbational form
    ! for filling new grids
    call ml_layout_build_mla(mla_temp,mla)
    call initialize_bc(the_bc_tower_temp,nlevs,mla_temp%pmask)
    do n = 1,nlevs
       call bc_tower_level_build(the_bc_tower_temp,n,mla_temp%la(n))
    end do    

    nlevs_temp = nlevs

    do n = 1,nlevs

       ! Create copies of the old data.
       call multifab_build(   uold_temp(n),mla_temp%la(n),   dm, ng_s)
       call multifab_build(   sold_temp(n),mla_temp%la(n),nscal, ng_s)
       call multifab_build(    gpi_temp(n),mla_temp%la(n),   dm, 0)
       call multifab_build(   dSdt_temp(n),mla_temp%la(n),    1, 0)
       call multifab_build(   S_cc_temp(n),mla_temp%la(n),    1, 0)
       call multifab_build( tag_mf_temp(n),mla_temp%la(n),    1, 0)

       call multifab_copy_c(   uold_temp(n),1,   uold(n),1,   dm)
       call multifab_copy_c(   sold_temp(n),1,   sold(n),1,nscal)
       call multifab_copy_c(    gpi_temp(n),1,    gpi(n),1,   dm)
       call multifab_copy_c(   dSdt_temp(n),1,   dSdt(n),1,    1)
       call multifab_copy_c(   S_cc_temp(n),1,   S_cc(n),1,    1)
       call multifab_copy_c( tag_mf_temp(n),1, tag_mf(n),1,    1)

       ! Get rid of the old data structures so we can create new ones 
       ! with the same names.
       call multifab_destroy(  uold(n))
       call multifab_destroy(  sold(n))
       call multifab_destroy(   gpi(n))
       call multifab_destroy(  dSdt(n))
       call multifab_destroy(  S_cc(n))
       call multifab_destroy(tag_mf(n))

    end do

    ! mba is big enough to hold max_levs levels
    ! even though we know we had nlevs last time, we might 
    ! want more or fewer levels after regrid (if nlevs < max_levs)
    call ml_boxarray_build_n(mba,max_levs,dm)

    do n = 1, max_levs-1
       mba%rr(n,:) = ref_ratio
    enddo

    if (associated(uold)) then
       deallocate(uold,sold,gpi,dSdt,S_cc,tag_mf)
    end if

    allocate(uold(max_levs),sold(max_levs),gpi(max_levs))
    allocate(dSdt(max_levs),S_cc(max_levs))
    allocate(tag_mf(max_levs))

    ! Copy the level 1 boxarray
    call copy(mba%bas(1),mla_temp%mba%bas(1))

    ! Copy the pd(:)
    mba%pd(1) = mla_temp%mba%pd(1)
    do n = 2, max_levs
       mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
    enddo

    ! Build the level 1 layout.
    call layout_build_ba(la_array(1),mba%bas(1),mba%pd(1),pmask, &
                         mapping=LA_EXPLICIT, explicit_mapping=get_proc(mla%la(1)))

    call destroy(mla)

    ! This makes sure the boundary conditions are properly defined everywhere
    call bc_tower_level_build(the_bc_tower,1,la_array(1))

    ! Build and fill the level 1 data only.
    call build_and_fill_data(1,la_array(1),mla_temp, &
                             uold     ,sold     ,gpi     ,dSdt     ,S_cc,     tag_mf, &
                             uold_temp,sold_temp,gpi_temp,dSdt_temp,S_cc_temp,tag_mf_temp, &
                             the_bc_tower,dm,ng_s,mba%rr(1,:))

    nl       = 1
    new_grid = .true.

    ng_buffer = 4

    do while ( (nl .lt. max_levs) .and. (new_grid) )

       ! Do we need finer grids?

       if (tagging_needs_ghost_cells) then
          ! Need to fill ghost cells here in case we use them in tagging
          if (nl .eq. 1) then
             call multifab_fill_boundary(sold(nl))
             call multifab_physbc(sold(nl),rho_comp,dm+rho_comp,nscal, &
                                  the_bc_tower%bc_tower_array(nl))
          else
             call multifab_fill_ghost_cells(sold(nl),sold(nl-1),sold(nl)%ng,mba%rr((nl-1),:), &
                                            the_bc_tower%bc_tower_array(nl-1), &
                                            the_bc_tower%bc_tower_array(nl), &
                                            rho_comp,dm+rho_comp,nscal)
          end if
       end if

       if (nl .eq. 1) then
          call make_new_grids(new_grid,la_array(nl),la_array(nl+1),sold(nl),dx(nl,1), &
                              amr_buf_width,ref_ratio,nl,max_grid_size_2,tag_mf(nl))
       else
          call make_new_grids(new_grid,la_array(nl),la_array(nl+1),sold(nl),dx(nl,1), &
                              amr_buf_width,ref_ratio,nl,max_grid_size_3,tag_mf(nl))
       end if

       if (new_grid) then

          call copy(mba%bas(nl+1),get_boxarray(la_array(nl+1)))

          ! Need to enforce proper nesting within the grid creation procedure so that we can
          !  fillpatch the new levels.
          if (nl .ge. 2) then

             ! Test on whether grids are already properly nested
             if (.not. ml_boxarray_properly_nested(mba, ng_buffer, pmask, 2, nl+1)) then

                do n = 2,nl  
                   ! Delete old multifabs so that we can rebuild them.
                   call destroy(  sold(n))
                   call destroy(  uold(n))
                   call destroy(   gpi(n))
                   call destroy(  dSdt(n))
                   call destroy(  S_cc(n))
                   call destroy(tag_mf(n))
                enddo

                call enforce_proper_nesting(mba,la_array,max_grid_size_2,max_grid_size_3)

                ! Loop over all the lower levels which we might have changed when we enforced proper nesting.
                do n = 2,nl

                   ! This makes sure the boundary conditions are properly defined everywhere
                   call bc_tower_level_build(the_bc_tower,n,la_array(n)) 
   
                   ! Rebuild the lower level data again if it changed.
                   call build_and_fill_data(n,la_array(n),mla_temp, &
                                            uold     ,sold     ,gpi     ,dSdt     ,S_cc,     tag_mf, &
                                            uold_temp,sold_temp,gpi_temp,dSdt_temp,S_cc_temp,tag_mf_temp, &
                                            the_bc_tower,dm,ng_s,mba%rr(n-1,:))
                end do

             end if
          end if

          ! Define bc_tower at level nl+1.
          call bc_tower_level_build(the_bc_tower,nl+1,la_array(nl+1))

          ! Build the level nl+1 data only.
          call build_and_fill_data(nl+1,la_array(nl+1),mla_temp, &
                                   uold     ,sold     ,gpi     ,dSdt     ,S_cc,     tag_mf, &
                                   uold_temp,sold_temp,gpi_temp,dSdt_temp,S_cc_temp,tag_mf_temp, &
                                   the_bc_tower,dm,ng_s,mba%rr(nl,:))

          nlevs = nl+1
          nlevs_radial = merge(1, nlevs, spherical .eq. 1)
          nl = nl+1

       endif

    enddo

    nlevs = nl
    nlevs_radial = merge(1, nlevs, spherical .eq. 1)

    call ml_layout_build_la_array(mla,la_array,mba,pmask,nlevs)

    ! We need to move data if a layout in la_array is not used in mla.
    ! We also need to destroy any unused layouts.
    allocate(uold_opt(nlevs),sold_opt(nlevs),gpi_opt(nlevs))
    allocate(dSdt_opt(nlevs),S_cc_opt(nlevs),tag_mf_opt(nlevs))
    do n=1,nlevs
       if (mla%la(n) .ne. la_array(n)) then
          call multifab_build(  uold_opt(n),mla%la(n),   dm, ng_s)
          call multifab_build(  sold_opt(n),mla%la(n),nscal, ng_s)
          call multifab_build(   gpi_opt(n),mla%la(n),   dm, 0)
          call multifab_build(  dSdt_opt(n),mla%la(n),    1, 0)
          call multifab_build(  S_cc_opt(n),mla%la(n),    1, 0)
          call multifab_build(tag_mf_opt(n),mla%la(n),    1, 0)

          call multifab_copy_c(  uold_opt(n),1,  uold(n),1,   dm)
          call multifab_copy_c(  sold_opt(n),1,  sold(n),1,nscal)
          call multifab_copy_c(   gpi_opt(n),1,   gpi(n),1,   dm)
          call multifab_copy_c(  dSdt_opt(n),1,  dSdt(n),1,    1)
          call multifab_copy_c(  S_cc_opt(n),1,  S_cc(n),1,    1)
          call multifab_copy_c(tag_mf_opt(n),1,tag_mf(n),1,  1)

          call multifab_destroy(  uold(n))
          call multifab_destroy(  sold(n))
          call multifab_destroy(   gpi(n))
          call multifab_destroy(  dSdt(n))
          call multifab_destroy(  S_cc(n))
          call multifab_destroy(tag_mf(n))
          
          call destroy(la_array(n))

          uold  (n) =   uold_opt(n)
          sold  (n) =   sold_opt(n)
          gpi   (n) =    gpi_opt(n)
          dSdt  (n) =   dSdt_opt(n)
          S_cc  (n) =   S_cc_opt(n)
          tag_mf(n) = tag_mf_opt(n)
       end if
    end do
    deallocate(uold_opt,sold_opt,gpi_opt,dSdt_opt,S_cc_opt,tag_mf_opt)

    ! This makes sure the boundary conditions are properly defined everywhere
    do n = 1, nlevs
       call bc_tower_level_build(the_bc_tower,n,mla%la(n))
    end do

    if (spherical .eq. 1) then

       ! convert (rho X) --> X in sold_temp 
       call convert_rhoX_to_X(sold_temp,.true.,mla_temp,the_bc_tower_temp%bc_tower_array)

       ! convert rho -> rho' in sold_temp
       call put_in_pert_form(mla_temp,sold_temp,rho0,dx,rho_comp,foextrap_comp,.true., &
                             the_bc_tower_temp%bc_tower_array)

       ! convert (rho h) -> (rho h)' in sold_temp
       call put_in_pert_form(mla_temp,sold_temp,rhoh0,dx,rhoh_comp,foextrap_comp,.true., &
                             the_bc_tower_temp%bc_tower_array)

       ! Delete old multifabs so that we can rebuild them.
       do n = 1, nlevs
          call destroy(  sold(n))
          call destroy(  uold(n))
          call destroy(   gpi(n))
          call destroy(  dSdt(n))
          call destroy(  S_cc(n))
          call destroy(tag_mf(n))
       end do

       ! Rebuild using the perturbational form
       do n = 1, nlevs
          call build_and_fill_data(n,mla%la(n),mla_temp, &
                                   uold     ,sold     ,gpi     ,dSdt     ,S_cc,     tag_mf, &
                                   uold_temp,sold_temp,gpi_temp,dSdt_temp,S_cc_temp,tag_mf_temp, &
                                   the_bc_tower,dm,ng_s,mba%rr(n-1,:))
       end do

       ! convert rho' -> rho in sold
       call put_in_pert_form(mla,sold,rho0,dx,rho_comp,dm+rho_comp,.false., &
                             the_bc_tower%bc_tower_array)

       ! convert (rho h)' -> (rho h) in sold
       call put_in_pert_form(mla,sold,rhoh0,dx,rhoh_comp,dm+rhoh_comp,.false., &
                             the_bc_tower%bc_tower_array)

       ! convert X --> (rho X) in sold 
       call convert_rhoX_to_X(sold,.false.,mla,the_bc_tower%bc_tower_array)

    end if

    do nl = 1,nlevs_temp
       call destroy(  uold_temp(nl))
       call destroy(  sold_temp(nl))
       call destroy(   gpi_temp(nl))
       call destroy(  dSdt_temp(nl))
       call destroy(  S_cc_temp(nl))
       call destroy(tag_mf_temp(nl))
    end do

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,sold,mla%mba%rr,the_bc_tower%bc_tower_array, &
                              icomp=1, &
                              bcomp=dm+rho_comp, &
                              nc=nscal, &
                              ng=sold(1)%ng)

    ! restrict data and fill all ghost cells
    call ml_restrict_and_fill(nlevs,uold,mla%mba%rr,the_bc_tower%bc_tower_array, &
                              icomp=1, &
                              bcomp=1, &
                              nc=dm, &
                              ng=uold(1)%ng)

    ! restrict data (no ghost cells)
    call ml_restrict_and_fill(nlevs,gpi,mla%mba%rr,the_bc_tower%bc_tower_array, &
                              icomp=1, &
                              bcomp=foextrap_comp, &
                              nc=dm, &
                              ng=gpi(1)%ng, &
                              same_boundary=.true.)

    ! restrict data (no ghost cells)
    call ml_restrict_and_fill(nlevs,S_cc,mla%mba%rr,the_bc_tower%bc_tower_array, &
                              icomp=1, &
                              bcomp=foextrap_comp, &
                              nc=1, &
                              ng=S_cc(1)%ng)

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
          write (un,1000) n, multifab_volume(sold(n),.false.), multifab_volume(sold(n),.true.), nboxes(sold(n)%la)
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
    call destroy(mla_temp)
    call bc_tower_destroy(the_bc_tower_temp)

    call destroy(bpt)

  end subroutine regrid

  subroutine build_and_fill_data(lev,la,mla_temp, &
                                 uold     ,sold     ,gpi     ,dSdt     ,S_cc,     tag_mf, &
                                 uold_temp,sold_temp,gpi_temp,dSdt_temp,S_cc_temp,tag_mf_temp, &
                                 the_bc_tower,dm,ng_s,rr)

    use fillpatch_module
    use ml_prolongation_module
    use probin_module, only : nodal
    use variables    , only : nscal, rho_comp, foextrap_comp

    use fabio_module

    integer                    , intent(in   ) :: lev, dm, ng_s, rr(:)
    type(layout)               , intent(in   ) :: la
    type(ml_layout)            , intent(in   ) :: mla_temp
    type(bc_tower)             , intent(inout) :: the_bc_tower
    type(multifab)   ,  pointer :: uold(:),sold(:),gpi(:),dSdt(:),S_cc(:)
    type(multifab)   ,  pointer :: tag_mf(:)
    type(multifab)             , intent(in   ) :: uold_temp(:),sold_temp(:),gpi_temp(:)
    type(multifab)             , intent(in   ) :: dSdt_temp(:),S_cc_temp(:), tag_mf_temp(:)
 
    integer :: d 

    ! Build the level lev data only.
    call multifab_build(  uold(lev), la,    dm, ng_s)
    call multifab_build(  sold(lev), la, nscal, ng_s)
    call multifab_build(   gpi(lev), la,    dm, 0)
    call multifab_build(  dSdt(lev), la,     1, 0)
    call multifab_build(  S_cc(lev), la,     1, 0)
    call multifab_build(tag_mf(lev), la,     1, 0)

    ! Fill the data in the new level lev state -- first from the coarser data if lev > 1.

    if (lev .gt. 1) then

       call fillpatch(uold(lev),uold(lev-1), &
                      0,rr, &
                      the_bc_tower%bc_tower_array(lev-1), &
                      the_bc_tower%bc_tower_array(lev  ), &
                      1,1,1,dm,no_final_physbc_input=.true.)
       call fillpatch(sold(lev),sold(lev-1), &
                      0,rr, &
                      the_bc_tower%bc_tower_array(lev-1), &
                      the_bc_tower%bc_tower_array(lev  ), &
                      1,1,dm+rho_comp,nscal,no_final_physbc_input=.true.)
       do d=1,dm
          call fillpatch(gpi(lev),gpi(lev-1), &
                         0,rr, &
                         the_bc_tower%bc_tower_array(lev-1), &
                         the_bc_tower%bc_tower_array(lev  ), &
                         d,d,foextrap_comp,1,no_final_physbc_input=.true.)
       end do
       call fillpatch(dSdt(lev),dSdt(lev-1), &
                      0,rr, &
                      the_bc_tower%bc_tower_array(lev-1), &
                      the_bc_tower%bc_tower_array(lev  ), &
                      1,1,foextrap_comp,1,no_final_physbc_input=.true.)
       call fillpatch(S_cc(lev),S_cc(lev-1), &
                      0,rr, &
                      the_bc_tower%bc_tower_array(lev-1), &
                      the_bc_tower%bc_tower_array(lev  ), &
                      1,1,foextrap_comp,1,no_final_physbc_input=.true.) 
       call fillpatch(tag_mf(lev),tag_mf(lev-1), &
                      0,rr, &
                      the_bc_tower%bc_tower_array(lev-1), &
                      the_bc_tower%bc_tower_array(lev  ), &
                      1,1,foextrap_comp,1,no_final_physbc_input=.true.) 

    end if

    ! Copy from old data at current level, if it exists
    if (mla_temp%nlevel .ge. lev) then
       call multifab_copy_c(  uold(lev),1,  uold_temp(lev),1,   dm)
       call multifab_copy_c(  sold(lev),1,  sold_temp(lev),1,nscal)
       call multifab_copy_c(   gpi(lev),1,   gpi_temp(lev),1,   dm)
       call multifab_copy_c(  dSdt(lev),1,  dSdt_temp(lev),1,    1)
       call multifab_copy_c(  S_cc(lev),1,  S_cc_temp(lev),1,    1)
       call multifab_copy_c(tag_mf(lev),1,tag_mf_temp(lev),1,    1)
    end if

  end subroutine build_and_fill_data

end module regrid_module
