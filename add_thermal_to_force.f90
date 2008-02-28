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

  subroutine add_thermal_to_force(nlevs,force,thermal,the_bc_level,mla)

    use bl_prof_module
    use multifab_physbc_module
    use multifab_fill_ghost_module
    use variables, only: rhoh_comp, foextrap_comp 
    use ml_restriction_module, only : ml_cc_restriction_c

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: force(:)
    type(multifab) , intent(in   ) :: thermal(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla

    ! local
    integer :: n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "add_thermal_to_force")

    do n=1,nlevs
       call multifab_plus_plus_c(force(n),rhoh_comp,thermal(n),1,1)
    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(force(nlevs),rhoh_comp,1)

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(force(nlevs),rhoh_comp,foextrap_comp,1,the_bc_level(nlevs))

    else
    
       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(force(n-1),rhoh_comp,force(n),rhoh_comp, &
                                   mla%mba%rr(n-1,:),1)
          
          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(force(n),force(n-1), &
                                         force(n)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         rhoh_comp,foextrap_comp,1)
       enddo
       
    end if

    call destroy(bpt)
       
  end subroutine add_thermal_to_force
  
end module add_thermal_to_force_module
