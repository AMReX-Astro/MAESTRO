module mac_hypre_module
 
  use bl_types
  use ml_layout_module
  use define_bc_module
  use multifab_module
  use bndry_reg_module
 
  implicit none
 
  private
 
  public :: mac_hypre
 
contains

  subroutine mac_hypre(mla,rh,phi,fine_flx,alpha,beta,dx,the_bc_tower,bc_comp, &
                       stencil_order,ref_ratio,rel_solver_eps,abs_solver_eps)

    use mac_multigrid_module, only : mac_multigrid

    type(ml_layout), intent(in   )        :: mla
    integer        , intent(in   )        :: stencil_order
    integer        , intent(in   )        :: ref_ratio(:,:)
    real(dp_t)     , intent(in)           :: dx(:,:)
    type(bc_tower) , intent(in)           :: the_bc_tower
    integer        , intent(in   )        :: bc_comp
    type(multifab) , intent(in   )        :: alpha(:), beta(:,:)
    type(multifab) , intent(inout)        ::    rh(:),  phi(:)
    type(bndry_reg), intent(inout)        :: fine_flx(2:)
    real(dp_t)     , intent(in)           :: rel_solver_eps 
    real(dp_t)     , intent(in)           :: abs_solver_eps 

    call mac_multigrid(mla,rh,phi,fine_flx,alpha,beta,dx,the_bc_tower,bc_comp,&
                       stencil_order,ref_ratio,rel_solver_eps,abs_solver_eps)

  end subroutine mac_hypre

end module mac_hypre_module
