! If flag = .true. then convert the state array species from (rho X)
! to X.  If flag = .false. then convert X to (rho X).  Note, this only
! applies when we are coming in with the full state (not the
! perturbational state).  This does not touch the base state.

module convert_rhoX_to_X_module

  use multifab_module

  implicit none

  private

  public :: convert_rhoX_to_X
  
contains

  subroutine convert_rhoX_to_X(nlevs,s,dx,flag,mla,the_bc_level)

    use geometry, only: spherical
    use network, only: nspec
    use variables, only: spec_comp, foextrap_comp, nscal
    use ml_layout_module
    use define_bc_module
    use ml_restriction_module, only: ml_cc_restriction
    use multifab_fill_ghost_module

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    logical        , intent(in   ) :: flag
    type(ml_layout), intent(inout) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local variables
    real(kind=dp_t), pointer::  sp(:,:,:,:)
    integer :: i,lo(s(1)%dim),hi(s(1)%dim),ng,dm,n,comp

    ng = s(1)%ng
    dm = s(1)%dim

    do n=1,nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n),i) ) cycle
          sp => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call convert_rhoX_to_X_2d(n,sp(:,:,1,:),lo,hi,ng,flag)
          case (3)
             call convert_rhoX_to_X_3d(n,sp(:,:,:,:),lo,hi,ng,flag)
          end select
       end do

       call multifab_fill_boundary_c(s(n),spec_comp,nspec)

    end do

    do n=nlevs,2,-1
       call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))

       call multifab_fill_ghost_cells(s(n),s(n-1), &
                                      s(n)%ng,mla%mba%rr(n-1,:), &
                                      the_bc_level(n-1),the_bc_level(n), &
                                      spec_comp,dm+spec_comp,nspec)

    end do

  end subroutine convert_rhoX_to_X

  subroutine convert_rhoX_to_X_2d(n,s,lo,hi,ng,flag)

    use network, only: nspec
    use variables, only: spec_comp, rho_comp
    use geometry, only: nr

    integer        , intent(in   ) :: n
    integer        , intent(in   ) ::  lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,:)
    logical        , intent(in   ) :: flag

    ! Local variables
    integer         :: i,j,r,comp

    if (flag) then

       ! convert (rho X) -> X
       do comp = spec_comp, spec_comp+nspec-1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                s(i,j,comp) = s(i,j,comp)/s(i,j,rho_comp)
             end do
          end do
       end do

    else

       ! convert X -> (rho X)
       do comp = spec_comp, spec_comp+nspec-1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                s(i,j,comp) = s(i,j,rho_comp)*s(i,j,comp)
             end do
          end do
       end do

    end if

  end subroutine convert_rhoX_to_X_2d

  subroutine convert_rhoX_to_X_3d(n,s,lo,hi,ng,flag)

    use network, only: nspec
    use variables, only: spec_comp, rho_comp
    use geometry, only: nr

    integer        , intent(in   ) :: n
    integer        , intent(in   ) ::  lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    logical        , intent(in   ) :: flag

    ! Local variables
    integer         :: i,j,k,r,comp

    if (flag) then

       ! convert (rho X) -> X
       do comp = spec_comp, spec_comp+nspec-1
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   s(i,j,k,comp) = s(i,j,k,comp)/s(i,j,k,rho_comp)
                end do
             end do
          end do
       end do

    else

       ! convert X -> (rho X)
       do comp = spec_comp, spec_comp+nspec-1
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   s(i,j,k,comp) = s(i,j,k,rho_comp)*s(i,j,k,comp)
                end do
             end do
          end do
       end do

    end if

  end subroutine convert_rhoX_to_X_3d

end module convert_rhoX_to_X_module
