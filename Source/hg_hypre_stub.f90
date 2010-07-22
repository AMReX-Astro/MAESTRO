module hg_hypre_module

  use bl_types
  use ml_layout_module
  use define_bc_module
  use multifab_module

  implicit none

  private

  public :: hg_hypre

contains 

  subroutine hg_hypre(mla,rh,unew,rhohalf,phi,dx,the_bc_tower, &
                      stencil_type,rel_solver_eps,abs_solver_eps,divu_rhs)
 
    use hg_multigrid_module, only : hg_multigrid

    type(ml_layout), intent(in   ) :: mla
    type(multifab ), intent(inout) :: rh(:)
    type(multifab ), intent(inout) :: unew(:)
    type(multifab ), intent(in   ) :: rhohalf(:)
    type(multifab ), intent(inout) :: phi(:)
    real(dp_t)     , intent(in)    :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: stencil_type
    real(dp_t)     , intent(in)    :: rel_solver_eps
    real(dp_t)     , intent(in)    :: abs_solver_eps

    type(multifab ), intent(inout), optional :: divu_rhs(:)

    call hg_multigrid(mla,rh,unew,rhohalf,phi,dx,the_bc_tower, &
                      stencil_type,rel_solver_eps,abs_solver_eps,divu_rhs)

  end subroutine hg_hypre

end module hg_hypre_module
