module pert_form_module

  ! Put the state into perturbational form

  use multifab_module

  implicit none

  private

  public :: put_in_pert_form

contains

  subroutine put_in_pert_form(nlevs,s,base,dx,comp,ncomp,flag)

    use geometry, only: spherical

    integer        , intent(in   ) :: nlevs,comp,ncomp
    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: base(:,0:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    logical        , intent(in   ) :: flag

    ! Local variables
    real(kind=dp_t), pointer::  sp(:,:,:,:)
    integer :: i,lo(s(1)%dim),hi(s(1)%dim),ng,dm,n

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
             call pert_form_2d(sp(:,:,1,:),base(n,:,:),lo,hi,ng,comp,ncomp,flag)
          case (3)
             if (spherical .eq. 1) then
                call pert_form_3d_sphr(n,sp(:,:,:,:),base(n,:,:),lo,hi,ng,dx(n,:), &
                                       comp,ncomp,flag)
             else
                call pert_form_3d_cart(sp(:,:,:,:),base(n,:,:),lo,hi,ng,comp,ncomp,flag)
             end if
          end select
       end do

       call multifab_fill_boundary_c(s(n),comp,ncomp)

    end do

  end subroutine put_in_pert_form

  subroutine pert_form_2d(s,s0,lo,hi,ng,incomp,ncomp,flag)

    integer        , intent(in   ) ::  lo(:),hi(:),ng,incomp,ncomp
    real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(in   ) :: s0(0:,:)
    logical        , intent(in   ) :: flag

    ! Local variables
    integer         :: i,j,comp

    if (flag) then
       do comp = incomp, incomp+ncomp-1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                s(i,j,comp) = s(i,j,comp) - s0(j,comp)
             end do
          end do
       end do
    else
       do comp = incomp, incomp+ncomp-1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                s(i,j,comp) = s(i,j,comp) + s0(j,comp)
             end do
          end do
       end do
    end if

  end subroutine pert_form_2d

  subroutine pert_form_3d_cart(s,s0,lo,hi,ng,incomp,ncomp,flag)

    integer        , intent(in   ) ::  lo(:),hi(:),ng,incomp,ncomp
    real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: s0(0:,:)
    logical        , intent(in   ) :: flag

    ! Local variables
    integer         :: i,j,k,comp

    if (flag) then
       do comp = incomp, incomp+ncomp-1
          do k = lo(3),hi(3)
             do j = lo(2),hi(2)
                do i = lo(1),hi(1)
                   s(i,j,k,comp) = s(i,j,k,comp) - s0(k,comp)
                end do
             end do
          end do
       end do
    else
       do comp = incomp, incomp+ncomp-1
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

  subroutine pert_form_3d_sphr(n,s,s0,lo,hi,ng,dx,incomp,ncomp,flag)

    use fill_3d_module

    integer        , intent(in   ) :: n,lo(:),hi(:),ng,incomp,ncomp
    real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: s0(0:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    logical        , intent(in   ) :: flag

    real(kind=dp_t), allocatable :: s0_cart(:,:,:)
    integer                      :: i,j,k,comp

    allocate(s0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

    if (flag) then
       do comp = incomp,incomp+ncomp-1
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
       do comp = incomp,incomp+ncomp-1
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
