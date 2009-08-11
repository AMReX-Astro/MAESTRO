module mac_multigrid_module

  use bl_types
  use ml_layout_module
  use define_bc_module
  use multifab_module
  use bndry_reg_module
  use bl_constants_module
  use sparse_solve_module
  use create_umac_grown_module
  use impose_phys_bcs_on_edges_module

  implicit none

  private

  public :: mac_multigrid

contains 
  subroutine mac_multigrid(mla,rh,phi,fine_flx,alpha,beta,dx,the_bc_tower,bc_comp, &
                           stencil_order,ref_ratio,umac_norm)
    use mg_module
    use coeffs_module
    use ml_solve_module
    use probin_module, only : mg_bottom_solver, max_mg_bottom_nlevels, verbose, mg_verbose, cg_verbose
    use geometry, only: dm, nlevs

    type(ml_layout), intent(in   )        :: mla
    integer        , intent(in   )        :: stencil_order
    integer        , intent(in   )        :: ref_ratio(:,:)
    real(dp_t)     , intent(in)           :: dx(:,:)
    type(bc_tower) , intent(in)           :: the_bc_tower
    integer        , intent(in   )        :: bc_comp
    type(multifab) , intent(in   )        :: alpha(:), beta(:)
    type(multifab) , intent(inout)        ::    rh(:),  phi(:)
    type(bndry_reg), intent(inout)        :: fine_flx(2:)
    real(dp_t)     , intent(in), optional :: umac_norm(:)

    type(layout  ) :: la
    type(boxarray) :: pdv
    type(box     ) :: pd

    type(multifab), allocatable :: coeffs(:)

    type(sparse)    :: sparse_object
    type(mg_tower)  :: mgt(mla%nlevel)
    integer         :: i, ns

    ! Bottom MGT stuff
    type(mg_tower)  :: bottom_mgt
    real(dp_t)      ::  coarse_xa(dm),  coarse_xb(dm)
    real(dp_t)      :: coarse_pxa(dm), coarse_pxb(dm)
    type(layout)    :: old_coarse_la,new_coarse_la
    type(layout)    :: old_la_grown, new_la_grown
    type(box)       :: coarse_pd,bxs
    type(boxarray)  :: ba_cc,new_coarse_ba
    type(multifab)  :: stored_coeffs, stored_coeffs_grown
    type(multifab)  :: new_coeffs_grown
    type(multifab), allocatable :: coarse_coeffs(:)
    integer         :: j,nx,mglev,bottom_box_size
    real(dp_t), pointer :: sc_orig(:,:,:,:), sc_grown(:,:,:,:)

    ! MG solver defaults
    integer    :: bottom_solver, bottom_max_iter
    integer    :: max_iter, min_width, max_nlevel
    integer    :: n, nu1, nu2, gamma, cycle, smoother
    integer    :: max_nlevel_in,do_diagnostics
    real(dp_t) :: eps,abs_eps,omega,bottom_solver_eps
    real(dp_t) ::  xa(dm),  xb(dm)
    real(dp_t) :: pxa(dm), pxb(dm)
    real(dp_t) :: coarse_dx(dm)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mac_multigrid")

    !! Defaults:

    max_nlevel        = mgt(nlevs)%max_nlevel
    max_iter          = mgt(nlevs)%max_iter
    eps               = mgt(nlevs)%eps
    abs_eps           = mgt(nlevs)%abs_eps
    smoother          = mgt(nlevs)%smoother
    nu1               = mgt(nlevs)%nu1
    nu2               = mgt(nlevs)%nu2
    gamma             = mgt(nlevs)%gamma
    omega             = mgt(nlevs)%omega
    cycle             = mgt(nlevs)%cycle
    bottom_solver     = mgt(nlevs)%bottom_solver
    bottom_solver_eps = mgt(nlevs)%bottom_solver_eps
    bottom_max_iter   = mgt(nlevs)%bottom_max_iter
    min_width         = mgt(nlevs)%min_width

    ! Note: put this here to minimize asymmetries - ASA
    if (nlevs .eq. 1) then
       eps = 1.d-12
    else if (nlevs .eq. 2) then
       eps = 1.d-11
    else
       eps = 1.d-10
    end if

    abs_eps = -1.0_dp_t
    if (present(umac_norm)) then
       do n = 1,nlevs
          abs_eps = max(abs_eps, umac_norm(n) / dx(n,1))
       end do
       abs_eps = eps * abs_eps
    end if

    if ( mg_bottom_solver >= 0 ) then
        if (mg_bottom_solver == 4 .and. phi(1)%nboxes == 1) then
           if (parallel_IOProcessor()) then
              print *,'Dont use mg_bottom_solver == 4 with only one grid -- '
              print *,'  Reverting to default bottom solver ',bottom_solver
           end if
        else if (mg_bottom_solver == 4 .and. max_mg_bottom_nlevels < 2) then
           if (parallel_IOProcessor()) then
              print *,'Dont use mg_bottom_solver == 4 with max_mg_bottom_nlevels < 2'
              print *,'  Reverting to default bottom solver ',bottom_solver
           end if
        else
           bottom_solver = mg_bottom_solver
        end if
    end if

    bottom_solver_eps = 1.d-3

    ! Note: put this here for robustness
    max_iter = 100

    if ( bottom_solver /= 0 .AND. max_iter == mgt(nlevs)%max_iter ) &
         max_iter = 1000

    ns = 1 + dm*3

    do n = nlevs, 1, -1

       if (n == 1) then
          max_nlevel_in = max_nlevel
       else
          if ( all(ref_ratio(n-1,:) == 2) ) then
             max_nlevel_in = 1
          else if ( all(ref_ratio(n-1,:) == 4) ) then
             max_nlevel_in = 2
          else
             call bl_error("MAC_MULTIGRID: confused about ref_ratio")
          end if
       end if

       pd = layout_get_pd(mla%la(n))

       call mg_tower_build(mgt(n), mla%la(n), pd, &
                           the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,bc_comp),&
                           dh = dx(n,:), &
                           ns = ns, &
                           smoother = smoother, &
                           nu1 = nu1, &
                           nu2 = nu2, &
                           gamma = gamma, &
                           cycle = cycle, &
                           omega = omega, &
                           bottom_solver = bottom_solver, &
                           bottom_max_iter = bottom_max_iter, &
                           bottom_solver_eps = bottom_solver_eps, &
                           max_iter = max_iter, &
                           max_nlevel = max_nlevel_in, &
                           min_width = min_width, &
                           eps = eps, &
                           abs_eps = abs_eps, &
                           verbose = mg_verbose, &
                           cg_verbose = cg_verbose, &
                           nodal = rh(nlevs)%nodal)

    end do

    !! Fill coefficient array
    do n = nlevs,1,-1

       allocate(coeffs(mgt(n)%nlevels))

       la = mla%la(n)
       pd = layout_get_pd(la)

       call multifab_build(coeffs(mgt(n)%nlevels), la, 1+dm, 1)
       call multifab_copy_c(coeffs(mgt(n)%nlevels),1,alpha(n),1, 1,ng=alpha(n)%ng)
       call multifab_copy_c(coeffs(mgt(n)%nlevels),2, beta(n),1,dm,ng= beta(n)%ng)

       do i = mgt(n)%nlevels-1, 1, -1
          call multifab_build(coeffs(i), mgt(n)%ss(i)%la, 1+dm, 1)
          call setval(coeffs(i), ZERO, 1, dm+1, all=.true.)
          call coarsen_coeffs(coeffs(i+1),coeffs(i))
       end do

       if (n > 1) then
          xa = HALF*ref_ratio(n-1,:)*mgt(n)%dh(:,mgt(n)%nlevels)
          xb = HALF*ref_ratio(n-1,:)*mgt(n)%dh(:,mgt(n)%nlevels)
       else
          xa = ZERO
          xb = ZERO
       end if

       pxa = ZERO
       pxb = ZERO
       do i = mgt(n)%nlevels, 1, -1
          pdv = layout_boxarray(mgt(n)%ss(i)%la)
          call stencil_fill_cc(mgt(n)%ss(i), coeffs(i), mgt(n)%dh(:,i), &
                               pdv, mgt(n)%mm(i), xa, xb, pxa, pxb, pd, stencil_order, &
                               the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,bc_comp))
       end do

       if ( n == 1 .and. bottom_solver == 3 ) then
          call sparse_build(mgt(n)%sparse_object, mgt(n)%ss(1), &
                            mgt(n)%mm(1), mgt(n)%ss(1)%la, stencil_order, mgt(nlevs)%verbose)
       end if

       ! Need to hang on to these coeffs at the bottom level of this mg
       if ( n == 1 .and. bottom_solver == 4 ) then
          call multifab_build(stored_coeffs, mgt(n)%ss(1)%la, nc=1+dm, ng=1)
          call multifab_copy_c(stored_coeffs,1,coeffs(n),1,dm+1,ng=coeffs(1)%ng)
       end if

       do i = mgt(n)%nlevels, 1, -1
          call destroy(coeffs(i))
       end do
       deallocate(coeffs)

    end do

    ! START OF BOTTOM_SOLVER == 4
    if (bottom_solver == 4) then

       ! Get the old/new coarse problem domain
       old_coarse_la = mgt(1)%ss(1)%la
       coarse_pd = layout_get_pd(old_coarse_la)

       ! Get the new coarse boxarray and layout
       call box_build_2(bxs,coarse_pd%lo(1:dm),coarse_pd%hi(1:dm))
       call boxarray_build_bx(new_coarse_ba,bxs)

       ! This is how many levels could be built if we made just one grid
       n = max_mg_levels(new_coarse_ba,min_width)

       ! This is the user-imposed limit
       n = min(n,max_mg_bottom_nlevels)

       if ( n .eq. 1) then
          call bl_error("DONT USE MG_BOTTOM_SOLVER == 4 WHEN BOTTOM GRID NOT PROPERLY DIVISIBLE : n = 1 ")
       end if

       bottom_box_size = 2**n

       do j = 1,dm
          nx = extent(bxs,j)
          if ( (bottom_box_size * (nx/bottom_box_size)) .ne. nx ) then 
             call bl_error("DONT USE MG_BOTTOM_SOLVER == 4 WHEN BOTTOM GRID NOT PROPERLY DIVISIBLE ")
          end if
       end do

       call boxarray_maxsize(new_coarse_ba,bottom_box_size)
       call layout_build_ba(new_coarse_la,new_coarse_ba,coarse_pd,old_coarse_la%lap%pmask)

       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          call print(layout_get_pd(old_coarse_la),"COARSE PD")
          print *,'ORIG MG NBOXES ',old_coarse_la%lap%nboxes
          print *,'NEW  MG NBOXES ',new_coarse_la%lap%nboxes
       end if

       coarse_dx(:) = dx(1,:) * 2**(mgt(1)%nlevels-1)

       call mg_tower_build(bottom_mgt, new_coarse_la, coarse_pd, &
                           the_bc_tower%bc_tower_array(1)%ell_bc_level_array(0,:,:,bc_comp),&
                           dh = coarse_dx, &
                           ns = ns, &
                           smoother = smoother, &
                           nu1 = nu1, &
                           nu2 = nu2, &
                           gamma = gamma, &
                           cycle = cycle, &
                           omega = omega, &
                           bottom_solver = 1, &
                           bottom_max_iter = bottom_max_iter, &
                           bottom_solver_eps = bottom_solver_eps, &
                           max_iter = max_iter, &
                           max_nlevel = max_nlevel, &
                           min_width = min_width, &
                           eps = eps, &
                           abs_eps = abs_eps, &
                           verbose = mg_verbose, &
                           cg_verbose = cg_verbose, &
                           nodal = rh(1)%nodal)

       ! START SPECIAL COPY
       ! Here we do special stuff to be able to copy the ghost cells of stored_coeffs into
       !   the ghost cells of coarse_coeffs(bottom)

       ! Make sure to do this before the copy so we get all the data
       call multifab_fill_boundary(stored_coeffs)

       mglev = bottom_mgt%nlevels

       allocate(coarse_coeffs(mglev))
       call multifab_build(coarse_coeffs(mglev), new_coarse_la, 1+dm, 1)
       call setval(coarse_coeffs(mglev),ZERO,all=.true.)

       do j = 1, dm
          call boxarray_build_copy(ba_cc,get_boxarray(stored_coeffs))
          call boxarray_grow(ba_cc,1,j, 1)
          call layout_build_ba(old_la_grown,ba_cc,pmask=old_coarse_la%lap%pmask, &
                               explicit_mapping=get_proc(old_coarse_la))
          call destroy(ba_cc)
          call multifab_build(stored_coeffs_grown,old_la_grown,1,ng=0)

          ! Note that stored_coeffs_grown only has one component at a time
          do i = 1, stored_coeffs_grown%nboxes
             if (remote(stored_coeffs_grown,i)) cycle 
             sc_orig  => dataptr(stored_coeffs      ,i,get_pbox(stored_coeffs_grown,i),j+1,1)
             sc_grown => dataptr(stored_coeffs_grown,i,get_pbox(stored_coeffs_grown,i),  1,1)
             sc_grown = sc_orig
          end do

          ! Note that new_coeffs_grown only has one component at a time
          call boxarray_build_copy(ba_cc,new_coarse_ba)
          call boxarray_grow(ba_cc,1,j, 1)
          call layout_build_ba(new_la_grown,ba_cc,pmask=old_coarse_la%lap%pmask, &
                               explicit_mapping=get_proc(new_coarse_la))
          call destroy(ba_cc)
          call multifab_build(new_coeffs_grown,new_la_grown,1,ng=0)
          call multifab_copy_c(new_coeffs_grown,1,stored_coeffs_grown,1,1)
          call destroy(stored_coeffs_grown)
          call destroy(old_la_grown)

          do i = 1, new_coeffs_grown%nboxes
             if (remote(new_coeffs_grown,i)) cycle 
             sc_orig  => dataptr(coarse_coeffs(mglev),i,get_pbox(new_coeffs_grown,i),j+1,1)
             sc_grown => dataptr(new_coeffs_grown    ,i,get_pbox(new_coeffs_grown,i),1  ,1)
             sc_orig = sc_grown
          end do

          call destroy(new_coeffs_grown)
          call destroy(new_la_grown)

       end do
       call destroy(new_coarse_ba)
       !   END SPECIAL COPY

       do i = mglev-1, 1, -1
          call multifab_build(coarse_coeffs(i), bottom_mgt%ss(i)%la, 1+dm, 1)
          call setval(coarse_coeffs(i), ZERO, 1, dm+1, all=.true.)
          call coarsen_coeffs(coarse_coeffs(i+1),coarse_coeffs(i))
       end do

       coarse_xa = ZERO
       coarse_xb = ZERO
       coarse_pxa = ZERO
       coarse_pxb = ZERO

       do i = mglev, 1, -1
          pdv = layout_boxarray(bottom_mgt%ss(i)%la)
          call stencil_fill_cc(bottom_mgt%ss(i), coarse_coeffs(i), bottom_mgt%dh(:,i), &
                               pdv, bottom_mgt%mm(i), coarse_xa, coarse_xb, coarse_pxa, coarse_pxb, &
                               coarse_pd, stencil_order, &
                               the_bc_tower%bc_tower_array(1)%ell_bc_level_array(0,:,:,bc_comp))
       end do

       do i = mglev, 1, -1
          call destroy(coarse_coeffs(i))
       end do
       deallocate(coarse_coeffs)
       call destroy(stored_coeffs)
    end if
    ! END   OF BOTTOM_SOLVER == 4

    if (mg_verbose >= 3) then
       do_diagnostics = 1
    else
       do_diagnostics = 0
    end if

    if (bottom_solver == 4) then
       call ml_cc_solve(mla, mgt, rh, phi, fine_flx, ref_ratio, do_diagnostics, bottom_mgt)
    else
       call ml_cc_solve(mla, mgt, rh, phi, fine_flx, ref_ratio, do_diagnostics)
    end if

    if ( bottom_solver == 3 ) then
       call sparse_destroy(sparse_object)
    else if (bottom_solver == 4) then
       call destroy(new_coarse_la)
    end if

    do n = 1, nlevs
       call mg_tower_destroy(mgt(n))
    end do

    if (bottom_solver == 4) then
       call mg_tower_destroy(bottom_mgt)
    end if

    call destroy(bpt)

  end subroutine mac_multigrid

end module mac_multigrid_module
