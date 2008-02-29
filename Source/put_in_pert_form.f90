! If flag = .true. then put the state into perturbational form.  That
! is, construct S' = S - S_0 for each of the desired state variables.
! Here, the base state, S_0, is left unchanged.
!
! If flag = .false. then convert from perturbational form to full 
! state form (S = S' + S_0)
!

module pert_form_module

  use multifab_module

  implicit none

  private

  public :: put_in_pert_form

  
contains

  subroutine put_in_pert_form(nlevs,s,base,dx,startcomp,numcomp,flag,mla,the_bc_level)

    use geometry, only: spherical
    use variables, only: foextrap_comp, nscal
    use ml_layout_module
    use define_bc_module
    use ml_restriction_module, only: ml_cc_restriction
    use multifab_fill_ghost_module
    use multifab_physbc_module

    integer        , intent(in   ) :: nlevs,startcomp,numcomp
    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: base(:,0:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    logical        , intent(in   ) :: flag
    type(ml_layout), intent(inout) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local variables
    real(kind=dp_t), pointer::  sp(:,:,:,:)
    integer :: lo(s(1)%dim),hi(s(1)%dim)
    integer :: i,ng,dm,n,comp,bc_comp

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
             call pert_form_2d(sp(:,:,1,:),base(n,:,:),lo,hi,ng,startcomp,numcomp,flag)
          case (3)
             if (spherical .eq. 1) then
                call pert_form_3d_sphr(n,sp(:,:,:,:),base(n,:,:),lo,hi,ng,dx(n,:), &
                                       startcomp,numcomp,flag)
             else
                call pert_form_3d_cart(sp(:,:,:,:),base(n,:,:),lo,hi,ng,startcomp,numcomp,flag)
             end if
          end select
       end do
    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(s(nlevs),startcomp,numcomp)

       do comp = startcomp,startcomp+numcomp-1
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
          
          do comp = startcomp,startcomp+numcomp-1
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

  end subroutine put_in_pert_form

  subroutine pert_form_2d(s,s0,lo,hi,ng,startcomp,numcomp,flag)

    integer        , intent(in   ) ::  lo(:),hi(:),ng,startcomp,numcomp
    real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(in   ) :: s0(0:,:)
    logical        , intent(in   ) :: flag

    ! Local variables
    integer         :: i,j,comp

    if (flag) then
       do comp = startcomp, startcomp+numcomp-1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                s(i,j,comp) = s(i,j,comp) - s0(j,comp)
             end do
          end do
       end do
    else
       do comp = startcomp, startcomp+numcomp-1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                s(i,j,comp) = s(i,j,comp) + s0(j,comp)
             end do
          end do
       end do
    end if

  end subroutine pert_form_2d

  subroutine pert_form_3d_cart(s,s0,lo,hi,ng,startcomp,numcomp,flag)

    integer        , intent(in   ) ::  lo(:),hi(:),ng,startcomp,numcomp
    real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: s0(0:,:)
    logical        , intent(in   ) :: flag

    ! Local variables
    integer         :: i,j,k,comp

    if (flag) then
       do comp = startcomp, startcomp+numcomp-1
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   s(i,j,k,comp) = s(i,j,k,comp) - s0(k,comp)
                end do
             end do
          end do
       end do
    else
       do comp = startcomp, startcomp+numcomp-1
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   s(i,j,k,comp) = s(i,j,k,comp) + s0(k,comp)
                end do
             end do
          end do
       end do
    end if

  end subroutine pert_form_3d_cart

  subroutine pert_form_3d_sphr(n,s,s0,lo,hi,ng,dx,startcomp,numcomp,flag)

    use fill_3d_module

    integer        , intent(in   ) :: n,lo(:),hi(:),ng,startcomp,numcomp
    real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: s0(0:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    logical        , intent(in   ) :: flag

    real(kind=dp_t), allocatable :: s0_cart(:,:,:)
    integer                      :: i,j,k,comp

    allocate(s0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

    if (flag) then
       do comp = startcomp,startcomp+numcomp-1
          call fill_3d_data(n,s0_cart,s0(0:,comp),lo,hi,dx,0)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   s(i,j,k,comp) = s(i,j,k,comp) - s0_cart(i,j,k)
                end do
             end do
          end do
       end do
    else
       do comp = startcomp,startcomp+numcomp-1
          call fill_3d_data(n,s0_cart,s0(0:,comp),lo,hi,dx,0)
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   s(i,j,k,comp) = s(i,j,k,comp) + s0_cart(i,j,k)
                end do
             end do
          end do
       end do
    end if

    deallocate(s0_cart)

  end subroutine pert_form_3d_sphr

end module pert_form_module
