module mac_multigrid_module

  use bl_types
  use ml_layout_module
  use define_bc_module
  use multifab_module
  use bndry_reg_module
  use bl_constants_module

  implicit none

  private

  public :: mac_multigrid

contains

  subroutine mac_multigrid(mla,rh,phi,fine_flx,alpha,beta,dx,the_bc_tower,bc_comp,&
                           stencil_order,rel_solver_eps,abs_solver_eps)

    use cc_stencil_fill_module, only : stencil_fill_cc_all_mglevels
    use mg_module             , only : mg_tower, mg_tower_build, mg_tower_destroy
    use ml_solve_module       , only : ml_cc_solve
    use probin_module         , only : mg_verbose, cg_verbose, mg_bottom_solver, max_mg_bottom_nlevels
    use stencil_types_module

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rh(:),phi(:)
    type(bndry_reg), intent(inout) :: fine_flx(2:)
    type(multifab) , intent(in   ) :: alpha(:), beta(:,:)
    real(dp_t)     , intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: bc_comp
    integer        , intent(in   ) :: stencil_order
    real(dp_t)     , intent(in   ) :: rel_solver_eps
    real(dp_t)     , intent(in   ) :: abs_solver_eps

    integer         :: nlevs

    ! MG solver defaults
    integer :: do_diagnostics

    type(bl_prof_timer), save :: bpt
 
    call build(bpt, "mac_multigrid")

    nlevs = mla%nlevel

    if (mg_verbose >= 3) then
       do_diagnostics = 1
    else
       do_diagnostics = 0
    end if

    call ml_cc_solve(mla,rh,phi,fine_flx,alpha,beta,dx,the_bc_tower,bc_comp, &
                     eps = rel_solver_eps, &
                     abs_eps = abs_solver_eps, &
                     bottom_solver_eps = 1.d-3, &
                     bottom_solver = mg_bottom_solver, &
                     max_bottom_nlevel = max_mg_bottom_nlevels, &
                     do_diagnostics = do_diagnostics, &
                     verbose = mg_verbose, &
                     cg_verbose = cg_verbose, &
                     stencil_order = stencil_order)

    call destroy(bpt)

  end subroutine mac_multigrid

end module mac_multigrid_module
