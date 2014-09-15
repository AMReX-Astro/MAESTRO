module mac_multigrid_module

  use bl_types
  use ml_layout_module
  use define_bc_module
  use multifab_module
  use bndry_reg_module
  use stencil_types_module
  use bl_constants_module

  implicit none

  private

  public :: mac_multigrid

contains 

  subroutine mac_multigrid(mla,rh,phi,fine_flx,alpha,beta,dx,&
                           the_bc_tower,bc_comp,stencil_order,ref_ratio,&
                           rel_solver_eps,abs_solver_eps)

    use cc_stencil_fill_module, only : stencil_fill_cc_all_mglevels
    use mg_module             , only : mg_tower, mg_tower_build, mg_tower_destroy
    use ml_solve_module       , only : ml_cc_solve
    use probin_module, only : mg_bottom_solver, max_mg_bottom_nlevels, &
                              mg_verbose, cg_verbose, mg_bottom_nu, mg_cycle_type
    use mg_eps_module, only: eps_mac_bottom

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rh(:),phi(:)
    type(bndry_reg), intent(inout) :: fine_flx(2:)
    type(multifab) , intent(in   ) :: alpha(:), beta(:,:)
    real(dp_t)     , intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: bc_comp
    integer        , intent(in   ) :: stencil_order
    integer        , intent(in   ) :: ref_ratio(:,:)
    real(dp_t)     , intent(in   ) :: rel_solver_eps
    real(dp_t)     , intent(in   ) :: abs_solver_eps

    type(layout  ) :: la
    type(box     ) :: pd

    type(multifab), allocatable :: cell_coeffs(:)
    type(multifab), allocatable :: edge_coeffs(:,:)

    type(mg_tower)  :: mgt(mla%nlevel)
    integer         :: dm,nlevs

    ! MG solver defaults
    integer :: bottom_solver, bottom_max_iter
    integer    :: max_iter
    integer    :: min_width
    integer    :: max_nlevel
    integer    :: d, n, nu1, nu2, smoother
    integer    :: max_nlevel_in
    integer    :: do_diagnostics
    integer    :: stencil_type
    real(dp_t) :: bottom_solver_eps
    real(dp_t) ::  xa(mla%dim),  xb(mla%dim)
    real(dp_t) :: pxa(mla%dim), pxb(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mac_multigrid")

    dm = mla%dim
    nlevs = mla%nlevel

    !! Defaults:
    max_nlevel        = mgt(nlevs)%max_nlevel
    max_iter          = mgt(nlevs)%max_iter
    smoother          = mgt(nlevs)%smoother
    nu1               = mgt(nlevs)%nu1
    nu2               = mgt(nlevs)%nu2
    bottom_solver     = mgt(nlevs)%bottom_solver
    bottom_solver_eps = mgt(nlevs)%bottom_solver_eps
    bottom_max_iter   = mgt(nlevs)%bottom_max_iter
    min_width         = mgt(nlevs)%min_width

    if ( mg_bottom_solver >= 0 ) then
        if (mg_bottom_solver == 4 .and. nboxes(phi(1)%la) == 1) then
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

    bottom_solver_eps = eps_mac_bottom

    ! Note: put this here for robustness
    max_iter = 100

    if ( bottom_solver /= 0 .AND. max_iter == mgt(nlevs)%max_iter ) &
         max_iter = 1000

    do n = nlevs, 1, -1

       if (dm == 1) then
          max_nlevel_in = 1
       else if (n == 1) then
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

       stencil_type = CC_CROSS_STENCIL

       call mg_tower_build(mgt(n), mla%la(n), pd, &
                           the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,bc_comp),&
                           stencil_type, &
                           dh = dx(n,:), &
                           smoother = smoother, &
                           nu1 = nu1, &
                           nu2 = nu2, &
                           nub = mg_bottom_nu, &
                           cycle_type = mg_cycle_type, &
                           bottom_solver = bottom_solver, &
                           bottom_max_iter = bottom_max_iter, &
                           bottom_solver_eps = bottom_solver_eps, &
                           max_iter = max_iter, &
                           max_nlevel = max_nlevel_in, &
                           max_bottom_nlevel = max_mg_bottom_nlevels, &
                           min_width  = min_width, &
                           eps        = rel_solver_eps, &
                           abs_eps    = abs_solver_eps, &
                           verbose    = mg_verbose, &
                           cg_verbose = cg_verbose, &
                           nodal = nodal_flags(rh(nlevs)))
    end do

    !! Fill coefficient array
    do n = nlevs,1,-1

       allocate(cell_coeffs(mgt(n)%nlevels))
       allocate(edge_coeffs(mgt(n)%nlevels,dm))

       la = mla%la(n)

       call multifab_build(cell_coeffs(mgt(n)%nlevels),la,1,1)
       call multifab_copy_c(cell_coeffs(mgt(n)%nlevels),1,alpha(n),1, 1,ng=nghost(alpha(n)))

       do d = 1, dm
          call multifab_build_edge(edge_coeffs(mgt(n)%nlevels,d),la,1,1,d)
          call multifab_copy_c(edge_coeffs(mgt(n)%nlevels,d),1,beta(n,d),1,1,ng=nghost(beta(n,d)))
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
       call stencil_fill_cc_all_mglevels(mgt(n), cell_coeffs, edge_coeffs, xa, xb, pxa, pxb, stencil_order, &
                                         the_bc_tower%bc_tower_array(n)%ell_bc_level_array(0,:,:,bc_comp))

       call destroy(cell_coeffs(mgt(n)%nlevels))
       deallocate(cell_coeffs)

       do d = 1, dm
          call destroy(edge_coeffs(mgt(n)%nlevels,d))
       end do
       deallocate(edge_coeffs)

    end do

    if (mg_verbose >= 3) then
       do_diagnostics = 1
    else
       do_diagnostics = 0
    end if

    call ml_cc_solve(mla, mgt, rh, phi, fine_flx, do_diagnostics)

    do n = 1, nlevs
       call mg_tower_destroy(mgt(n))
    end do

!   if (bottom_solver == 4) then
!      call mg_tower_destroy(mgt%bottom_mgt)
!   end if

    call destroy(bpt)

  end subroutine mac_multigrid

end module mac_multigrid_module
