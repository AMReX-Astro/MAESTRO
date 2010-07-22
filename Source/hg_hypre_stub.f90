module hg_hypre_module

  use bl_types
  use mg_module
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: hg_hypre

contains 

  subroutine hg_hypre(mla,rh,unew,rhohalf,phi,dx,the_bc_tower,stencil_type,divu_rhs,eps_in)
 
    use hg_multigrid_module, only : hg_multigrid

    type(ml_layout), intent(inout) :: mla
    type(multifab ), intent(inout) ::   rh(:)
    type(multifab ), intent(inout) :: unew(:)
    type(multifab ), intent(in   ) :: rhohalf(:)
    type(multifab ), intent(inout) :: phi(:)
    real(dp_t)     , intent(in)    :: dx(:,:)
    type(bc_tower ), intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: stencil_type

    type(multifab ), intent(in   ), optional :: divu_rhs(:)
    real(dp_t)     , intent(in)   , optional :: eps_in 

    call hg_multigrid(mla,rh,unew,rhohalf,phi,dx,the_bc_tower, &
                      stencil_type,divu_rhs,eps_in)

  end subroutine hg_hypre

end module hg_hypre_module
