! If flag = .true. then convert the state array species from (rho X)
! to X.  If flag = .false. then convert X to (rho X).  Note, this only
! applies when we are coming in with the full state (not the
! perturbational state).  This does not touch the base state.

module convert_rhoX_to_X_module

  use multifab_module, only: multifab, nghost, multifab_div_div_c, multifab_mult_mult_c
  use bl_prof_module, only: bl_prof_timer, build, destroy
  use ml_restrict_fill_module

  implicit none

  private

  public :: convert_rhoX_to_X, convert_rhoh_to_h
  
contains

  subroutine convert_rhoX_to_X(s,flag,mla,the_bc_level)

    use network, only: nspec
    use variables, only: spec_comp, foextrap_comp, rho_comp
    use ml_layout_module, only: ml_layout
    use define_bc_module, only: bc_level

    type(multifab) , intent(inout) :: s(:)
    logical        , intent(in   ) :: flag
    type(ml_layout), intent(inout) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local variables
    integer :: n,comp,dm,nlevs

    type(bl_prof_timer), save :: bpt

    call build(bpt, "convert_rhoX_to_X")

    dm = mla%dim
    nlevs = mla%nlevel

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

    ! restrict data and fill all ghost cells
    if (flag) then

       call ml_restrict_and_fill(nlevs,s,mla%mba%rr,the_bc_level, &
                                 icomp=spec_comp, &
                                 bcomp=foextrap_comp, &
                                 nc=nspec, &
                                 ng=s(1)%ng, &
                                 same_boundary=.true.)

    else

       call ml_restrict_and_fill(nlevs,s,mla%mba%rr,the_bc_level, &
                                 icomp=spec_comp, &
                                 bcomp=dm+spec_comp, &
                                 nc=nspec, &
                                 ng=s(1)%ng)

    end if

    call destroy(bpt)
    
  end subroutine convert_rhoX_to_X

  subroutine convert_rhoh_to_h(s,flag,mla,the_bc_level)

    use variables, only: rho_comp, rhoh_comp, foextrap_comp
    use ml_layout_module, only: ml_layout
    use define_bc_module, only: bc_level

    type(multifab) , intent(inout) :: s(:)
    logical        , intent(in   ) :: flag
    type(ml_layout), intent(inout) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local variables
    integer :: n,dm,nlevs

    type(bl_prof_timer), save :: bpt

    call build(bpt, "convert_rhoh_to_h")

    dm = mla%dim
    nlevs = mla%nlevel

    if (flag) then
       do n=1,nlevs
          call multifab_div_div_c(s(n),rhoh_comp,s(n),rho_comp,1)
       end do
    else
       do n=1,nlevs
          call multifab_mult_mult_c(s(n),rhoh_comp,s(n),rho_comp,1)
       end do
    end if

    ! restrict data and fill all ghost cells
    if (flag) then

       call ml_restrict_and_fill(nlevs,s,mla%mba%rr,the_bc_level, &
                                 icomp=rhoh_comp, &
                                 bcomp=foextrap_comp, &
                                 nc=1, &
                                 ng=s(1)%ng)

    else

       call ml_restrict_and_fill(nlevs,s,mla%mba%rr,the_bc_level, &
                                 icomp=rhoh_comp, &
                                 bcomp=dm+rhoh_comp, &
                                 nc=1, &
                                 ng=s(1)%ng)

    end if

    call destroy(bpt)
    
  end subroutine convert_rhoh_to_h
  
end module convert_rhoX_to_X_module
