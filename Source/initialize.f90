module initialize_module

  use define_bc_module
  use ml_layout_module
  use multifab_module
  use bc_module
  use probin_module, only: nlevs, nodal, test_set, prob_lo, prob_hi, bcx_lo, bcx_hi, &
       bcy_lo, bcy_hi, bcz_lo, bcz_hi
  use variables, only: nscal, rho_comp
  use geometry, only: dm
  use network, only: nspec

  implicit none

  private

  public :: initialize_from_restart, initialize_with_fixed_grids, &
       initialize_with_adaptive_grids, initialize_bc

contains
    
  subroutine initialize_from_restart(mla,restart,time,dt,dx,pmask,uold,sold,gpres,pres, &
                                     dSdt,Source_old,rho_omegadot2,the_bc_tower)

    use restart_module
    use ml_restriction_module
    use multifab_fill_ghost_module
    use multifab_physbc_module

    type(ml_layout),intent(out)   :: mla
    integer       , intent(in   ) :: restart
    real(dp_t)    , intent(  out) :: time,dt
    real(dp_t)    , pointer       :: dx(:,:)
    logical       , intent(in   ) :: pmask(:)
    type(multifab), pointer       :: uold(:),sold(:),gpres(:),pres(:),dSdt(:)
    type(multifab), pointer       :: Source_old(:),rho_omegadot2(:)
    type(bc_tower), intent(  out) :: the_bc_tower

    ! local
    type(multifab), pointer :: chkdata(:)
    type(multifab), pointer :: chk_p(:)
    type(multifab), pointer :: chk_dsdt(:)
    type(multifab), pointer :: chk_src_old(:)
    type(multifab), pointer :: chk_rho_omegadot2(:)

    integer     , allocatable :: domain_phys_bc(:,:)
    type(box)   , allocatable :: domain_boxes(:)

    type(ml_boxarray) :: mba

    integer :: n,d

    call fill_restart_data(restart, mba, chkdata, chk_p, chk_dsdt, chk_src_old, &
                           chk_rho_omegadot2, time, dt)

    call ml_layout_build(mla,mba,pmask)

    nlevs = mla%nlevel

    allocate(uold(nlevs),sold(nlevs),gpres(nlevs),pres(nlevs))
    allocate(dSdt(nlevs),Source_old(nlevs),rho_omegadot2(nlevs))

    do n = 1,nlevs
       call multifab_build(         uold(n), mla%la(n),    dm, 3)
       call multifab_build(         sold(n), mla%la(n), nscal, 3)
       call multifab_build(        gpres(n), mla%la(n),    dm, 1)
       call multifab_build(         pres(n), mla%la(n),     1, 1, nodal)
       call multifab_build(         dSdt(n), mla%la(n),     1, 0)
       call multifab_build(   Source_old(n), mla%la(n),     1, 1)
       call multifab_build(rho_omegadot2(n), mla%la(n), nspec, 1)
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
       call multifab_copy_c(rho_omegadot2(n),1,chk_rho_omegadot2(n),1,nspec)
       call destroy(chk_rho_omegadot2(n)%la)
       call destroy(chk_rho_omegadot2(n))
    end do
    
    deallocate(chkdata, chk_p, chk_dsdt, chk_src_old, chk_rho_omegadot2)

    allocate(dx(nlevs,dm))
    
    do d=1,dm
       dx(1,d) = (prob_hi(d)-prob_lo(d)) / real(extent(mla%mba%pd(1),d),kind=dp_t)
    end do
    do n = 2,nlevs
       dx(n,:) = dx(n-1,:) / mla%mba%rr(n-1,:)
    end do

    allocate(domain_phys_bc(dm,2))
    allocate(domain_boxes(nlevs))
    
    do n = 1,nlevs
       domain_boxes(n) = layout_get_pd(mla%la(n))
    end do

    ! Put the bc values from the inputs file into domain_phys_bc
    domain_phys_bc(1,1) = bcx_lo
    domain_phys_bc(1,2) = bcx_hi
    if (pmask(1)) then
       if (bcx_lo .ne. -1 .or. bcx_hi .ne. -1) then
          call bl_error('MUST HAVE BCX = -1 if PMASK = T')
       endif
    end if
    if (dm > 1) then
       domain_phys_bc(2,1) = bcy_lo
       domain_phys_bc(2,2) = bcy_hi
       if (pmask(2)) then
          if (bcy_lo .ne. -1 .or. bcy_hi .ne. -1) then
             call bl_error('MUST HAVE BCY = -1 if PMASK = T') 
          end if
       end if
    end if
    if (dm > 2) then
       domain_phys_bc(3,1) = bcz_lo
       domain_phys_bc(3,2) = bcz_hi
       if (pmask(3)) then
          if (bcz_lo .ne. -1 .or. bcz_hi .ne. -1) then
             call bl_error('MUST HAVE BCZ = -1 if PMASK = T')
          end if
       end if
    end if

    do d=1,dm
       if ( pmask(d) ) domain_phys_bc(d,:) = BC_PER
    end do

    ! Build the arrays for each grid from the domain_bc arrays.
    call bc_tower_build(the_bc_tower,mla,domain_phys_bc,domain_boxes)

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
                                         3,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n  ), &
                                         1,1,dm)
          call multifab_fill_ghost_cells(sold(n),sold(n-1), &
                                         3,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n  ), &
                                         1,dm+rho_comp,nscal)
       end do
       
    end if

    call destroy(mba)

  end subroutine initialize_from_restart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initialize_with_fixed_grids(mla,pmask,dx,uold,sold,gpres,pres, &
                                         dSdt,Source_old,rho_omegadot2,the_bc_tower)

    use box_util_module
    
    type(ml_layout),intent(out)   :: mla
    logical       , intent(in   ) :: pmask(:)
    real(dp_t)    , pointer       :: dx(:,:)
    type(multifab), pointer       :: uold(:),sold(:),gpres(:),pres(:),dSdt(:)
    type(multifab), pointer       :: Source_old(:),rho_omegadot2(:)
    type(bc_tower), intent(  out) :: the_bc_tower
    
    integer     , allocatable :: domain_phys_bc(:,:)
    type(box)   , allocatable :: domain_boxes(:)

    type(ml_boxarray)         :: mba

    integer :: n,d
    
    call read_a_hgproj_grid(mba,test_set)
    call ml_layout_build(mla,mba,pmask)
    
    ! check for proper nesting
    if (.not. ml_boxarray_properly_nested(mla%mba, 3, pmask)) then
       call bl_error('fixed_grids not properly nested')
    end if
    
    nlevs = mla%nlevel
    
    allocate(uold(nlevs),sold(nlevs),gpres(nlevs),pres(nlevs))
    allocate(dSdt(nlevs),Source_old(nlevs),rho_omegadot2(nlevs))
    
    do n = 1,nlevs
       call multifab_build(         uold(n), mla%la(n),    dm, 3)
       call multifab_build(         sold(n), mla%la(n), nscal, 3)
       call multifab_build(        gpres(n), mla%la(n),    dm, 1)
       call multifab_build(         pres(n), mla%la(n),     1, 1, nodal)
       call multifab_build(         dSdt(n), mla%la(n),     1, 0)
       call multifab_build(   Source_old(n), mla%la(n),     1, 1)
       call multifab_build(rho_omegadot2(n), mla%la(n), nspec, 1)

       call setval(       uold(n), 0.0_dp_t, all=.true.)
       call setval(       sold(n), 0.0_dp_t, all=.true.)
       call setval(      gpres(n), 0.0_dp_t, all=.true.)
       call setval(       pres(n), 0.0_dp_t, all=.true.)
       call setval( Source_old(n), 0.0_dp_t, all=.true.)
       call setval(       dSdt(n), 0.0_dp_t, all=.true.)
       call setval(rho_omegadot2(n),0.0_dp_t,all=.true.)
    end do

    allocate(dx(nlevs,dm))

    do d=1,dm
       dx(1,d) = (prob_hi(d)-prob_lo(d)) / real(extent(mla%mba%pd(1),d),kind=dp_t)
    end do
    do n = 2,nlevs
       dx(n,:) = dx(n-1,:) / mla%mba%rr(n-1,:)
    end do

    allocate(domain_phys_bc(dm,2))
    allocate(domain_boxes(nlevs))
    
    do n = 1,nlevs
       domain_boxes(n) = layout_get_pd(mla%la(n))
    end do

    ! Put the bc values from the inputs file into domain_phys_bc
    domain_phys_bc(1,1) = bcx_lo
    domain_phys_bc(1,2) = bcx_hi
    if (pmask(1)) then
       if (bcx_lo .ne. -1 .or. bcx_hi .ne. -1) then
          call bl_error('MUST HAVE BCX = -1 if PMASK = T')
       endif
    end if
    if (dm > 1) then
       domain_phys_bc(2,1) = bcy_lo
       domain_phys_bc(2,2) = bcy_hi
       if (pmask(2)) then
          if (bcy_lo .ne. -1 .or. bcy_hi .ne. -1) then
             call bl_error('MUST HAVE BCY = -1 if PMASK = T') 
          end if
       end if
    end if
    if (dm > 2) then
       domain_phys_bc(3,1) = bcz_lo
       domain_phys_bc(3,2) = bcz_hi
       if (pmask(3)) then
          if (bcz_lo .ne. -1 .or. bcz_hi .ne. -1) then
             call bl_error('MUST HAVE BCZ = -1 if PMASK = T')
          end if
       end if
    end if

    do d=1,dm
       if ( pmask(d) ) domain_phys_bc(d,:) = BC_PER
    end do

    ! Build the arrays for each grid from the domain_bc arrays.
    call bc_tower_build(the_bc_tower,mla,domain_phys_bc,domain_boxes)

    call destroy(mba)

  end subroutine initialize_with_fixed_grids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initialize_with_adaptive_grids(mla,pmask,dx,uold,sold,gpres,pres, &
                                            dSdt,Source_old,rho_omegadot2,the_bc_tower)

    use probin_module, only: n_cellx, n_celly, n_cellz, regrid_int, max_grid_size, &
         ref_ratio, max_levs

    type(ml_layout),intent(out)   :: mla
    logical       , intent(in   ) :: pmask(:)
    real(dp_t)    , pointer       :: dx(:,:)
    type(multifab), pointer       :: uold(:),sold(:),gpres(:),pres(:),dSdt(:)
    type(multifab), pointer       :: Source_old(:),rho_omegadot2(:)
    type(bc_tower), intent(  out) :: the_bc_tower

    integer     , allocatable :: domain_phys_bc(:,:)
    type(box)   , allocatable :: domain_boxes(:)

    integer           :: buf_wid
    type(layout)      :: la_array(max_levs)
    type(box)         :: bxs
    type(ml_boxarray) :: mba

    logical :: new_grid
    integer :: lo(dm), hi(dm)
    integer :: n, d

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

     allocate(uold(max_levs),sold(max_levs),gpres(max_levs),pres(max_levs))
     allocate(dSdt(max_levs),Source_old(max_levs),rho_omegadot2(max_levs))

       ! Build the level 1 boxarray
     call box_build_2(bxs,lo,hi)
     call boxarray_build_bx(mba%bas(1),bxs)
     call boxarray_maxsize(mba%bas(1),max_grid_size)

     ! build pd(:)
     mba%pd(1) = bxs
     do n = 2, max_levs
        mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
     enddo

    allocate(dx(max_levs,dm))

    do d=1,dm
       dx(1,d) = (prob_hi(d)-prob_lo(d)) / real(extent(mla%mba%pd(1),d),kind=dp_t)
    end do
    do n = 2,max_levs
       dx(n,:) = dx(n-1,:) / mla%mba%rr(n-1,:)
    end do

    allocate(domain_phys_bc(dm,2))
    allocate(domain_boxes(max_levs))
    
    do n = 1,max_levs
       domain_boxes(n) = layout_get_pd(mla%la(n))
    end do

    ! Put the bc values from the inputs file into domain_phys_bc
    domain_phys_bc(1,1) = bcx_lo
    domain_phys_bc(1,2) = bcx_hi
    if (pmask(1)) then
       if (bcx_lo .ne. -1 .or. bcx_hi .ne. -1) then
          call bl_error('MUST HAVE BCX = -1 if PMASK = T')
       endif
    end if
    if (dm > 1) then
       domain_phys_bc(2,1) = bcy_lo
       domain_phys_bc(2,2) = bcy_hi
       if (pmask(2)) then
          if (bcy_lo .ne. -1 .or. bcy_hi .ne. -1) then
             call bl_error('MUST HAVE BCY = -1 if PMASK = T') 
          end if
       end if
    end if
    if (dm > 2) then
       domain_phys_bc(3,1) = bcz_lo
       domain_phys_bc(3,2) = bcz_hi
       if (pmask(3)) then
          if (bcz_lo .ne. -1 .or. bcz_hi .ne. -1) then
             call bl_error('MUST HAVE BCZ = -1 if PMASK = T')
          end if
       end if
    end if

    do d=1,dm
       if ( pmask(d) ) domain_phys_bc(d,:) = BC_PER
    end do

    ! Build the arrays for each grid from the domain_bc arrays.
    call bc_tower_build(the_bc_tower,mla,domain_phys_bc,domain_boxes)

    ! this may be modified later since we don't know how many levels we actually
    ! need until we start initializing the data
    nlevs = max_levs

  end subroutine initialize_with_adaptive_grids

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine initialize_bc()

  end subroutine initialize_bc

end module initialize_module
