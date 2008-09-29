! If flag = .true. then convert the state array species from (rho X)
! to X.  If flag = .false. then convert X to (rho X).  Note, this only
! applies when we are coming in with the full state (not the
! perturbational state).  This does not touch the base state.

module convert_rhoX_to_X_module

  use multifab_module

  implicit none

  private

  public :: convert_rhoX_to_X, convert_rhoh_to_h
  
contains

  subroutine convert_rhoX_to_X(s,flag,mla,the_bc_level)

    use network, only: nspec
    use variables, only: spec_comp, foextrap_comp, rho_comp
    use ml_layout_module
    use define_bc_module
    use ml_restriction_module, only: ml_cc_restriction
    use multifab_fill_ghost_module
    use multifab_physbc_module
    use geometry, only: dm, nlevs

    type(multifab) , intent(inout) :: s(:)
    logical        , intent(in   ) :: flag
    type(ml_layout), intent(inout) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local variables
    integer :: n,comp,bc_comp

    if (flag) then
       do n=1,nlevs
          do comp=spec_comp,spec_comp+nspec-1
             call multifab_div_div_c(s(n),comp,s(n),rho_comp,1)
          end do
       end do
    else
       do n=1,nlevs
          do comp=spec_comp,spec_comp+nspec-1
             call multifab_mult_mult_c(s(n),comp,s(n),rho_comp,1)
          end do
       end do
    end if

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(s(nlevs),spec_comp,nspec)

       do comp = spec_comp,spec_comp+nspec-1

          if (flag) then
             bc_comp = foextrap_comp
          else
             bc_comp = dm+comp
          end if

          ! fill non-periodic domain boundary ghost cells
          call multifab_physbc(s(nlevs),comp,bc_comp,1,the_bc_level(nlevs))
       end do

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))
          
          do comp = spec_comp,spec_comp+nspec-1
             if (flag) then
                bc_comp = foextrap_comp
             else
                bc_comp = dm+comp
             end if
             
             ! fill level n ghost cells using interpolation from level n-1 data
             ! note that multifab_fill_boundary and multifab_physbc are called for
             ! both levels n-1 and n
             call multifab_fill_ghost_cells(s(n),s(n-1), &
                                            s(n)%ng,mla%mba%rr(n-1,:), &
                                            the_bc_level(n-1),the_bc_level(n), &
                                            comp,bc_comp,1)
          end do

       end do

    end if
    
  end subroutine convert_rhoX_to_X

  subroutine convert_rhoh_to_h(s,flag,mla,the_bc_level)

    use variables, only: rho_comp, rhoh_comp, foextrap_comp
    use ml_layout_module
    use define_bc_module
    use ml_restriction_module, only: ml_cc_restriction
    use multifab_fill_ghost_module
    use multifab_physbc_module
    use geometry, only: dm, nlevs

    type(multifab) , intent(inout) :: s(:)
    logical        , intent(in   ) :: flag
    type(ml_layout), intent(inout) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local variables
    integer :: n,bc_comp

    if (flag) then
       do n=1,nlevs
          call multifab_div_div_c(s(n),rhoh_comp,s(n),rho_comp,1)
       end do
    else
       do n=1,nlevs
          call multifab_mult_mult_c(s(n),rhoh_comp,s(n),rho_comp,1)
       end do
    end if
    
    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(s(nlevs),rhoh_comp,1)

       if (flag) then
          bc_comp = foextrap_comp
       else
          bc_comp = dm+rhoh_comp
       end if

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s(nlevs),rhoh_comp,bc_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))
          
          if (flag) then
             bc_comp = foextrap_comp
          else
             bc_comp = dm+rhoh_comp
          end if
             
          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s(n),s(n-1), &
                                         s(n)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         rhoh_comp,bc_comp,1)
       end do

    end if
    
  end subroutine convert_rhoh_to_h
  
end module convert_rhoX_to_X_module
