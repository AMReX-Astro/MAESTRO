module add_thermal_to_force_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: add_thermal_to_force

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_thermal_to_force(nlevs,force,thermal,the_bc_level,mla,dx)

    use multifab_physbc_module
    use multifab_fill_ghost_module
    use variables, only: rhoh_comp, foextrap_comp 
    use ml_restriction_module, only : ml_cc_restriction_c

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: force(:)
    type(multifab) , intent(in   ) :: thermal(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    ! local
    integer :: n

    do n=1,nlevs
       call multifab_plus_plus_c(force(n),rhoh_comp,thermal(n),1,1)

       call multifab_fill_boundary_c(force(n),rhoh_comp,1)
       call multifab_physbc(force(n),rhoh_comp,foextrap_comp,1,dx(n,:), &
                            the_bc_level(n))
    enddo
    
    do n=nlevs,2,-1
       call ml_cc_restriction_c(force(n-1),rhoh_comp,force(n),rhoh_comp, &
                                mla%mba%rr(n-1,:),1)

       call multifab_fill_ghost_cells(force(n),force(n-1), &
                                      force(n)%ng,mla%mba%rr(n-1,:), &
                                      the_bc_level(n-1),the_bc_level(n), &
                                      rhoh_comp,foextrap_comp,1)
    enddo
       
  end subroutine add_thermal_to_force
  
end module add_thermal_to_force_module
