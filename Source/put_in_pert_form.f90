module pert_form_module

  ! Put the state into perturbational form

  use bl_constants_module
  use geometry
  use variables
  use fill_3d_module

  implicit none

contains

  subroutine put_in_pert_form(s,base,dx,comp,ncomp,flag)

      integer        , intent(in   ) :: comp,ncomp
      type(multifab) , intent(inout) :: s
      real(kind=dp_t), intent(in   ) :: base(0:,:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      logical        , intent(in   ) :: flag

      ! Local variables
      real(kind=dp_t), pointer::  sp(:,:,:,:)
      integer :: i,lo(s%dim),hi(s%dim),ng,dm

      ng = s%ng
      dm = s%dim

      do i = 1, s%nboxes
         if ( multifab_remote(s, i) ) cycle
          sp => dataptr(s     , i)
          lo =  lwb(get_box(s, i))
          hi =  upb(get_box(s, i))
         select case (dm)
            case (2)
              call pert_form_2d(sp(:,:,1,:),base,lo,hi,ng,comp,ncomp,flag)
            case (3)
              if (spherical .eq. 1) then
                call pert_form_3d_sphr(sp(:,:,:,:),base,lo,hi,ng,dx,comp,ncomp,flag)
              else
                call pert_form_3d_cart(sp(:,:,:,:),base,lo,hi,ng,comp,ncomp,flag)
              end if
         end select
      end do
      call multifab_fill_boundary_c(s,comp,ncomp)

end subroutine put_in_pert_form

  subroutine pert_form_2d(s,s0,lo,hi,ng,comp,ncomp,flag)

      integer        , intent(in   ) ::  lo(:),hi(:),ng,comp,ncomp
      real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,:)
      real(kind=dp_t), intent(in   ) :: s0(0:,:)
      logical        , intent(in   ) :: flag

      ! Local variables
      integer         :: i,j,n

      if (flag) then
        do n = comp, comp+ncomp-1
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          s(i,j,n) = s(i,j,n) - s0(j,n)
        end do
        end do
        end do
      else
        do n = comp, comp+ncomp-1
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          s(i,j,n) = s(i,j,n) + s0(j,n)
        end do
        end do
        end do
      end if

  end subroutine pert_form_2d

  subroutine pert_form_3d_cart(s,s0,lo,hi,ng,comp,ncomp,flag)

      integer        , intent(in   ) ::  lo(:),hi(:),ng,comp,ncomp
      real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real(kind=dp_t), intent(in   ) :: s0(0:,:)
      logical        , intent(in   ) :: flag

      ! Local variables
      integer         :: i,j,k,n

      if (flag) then
        do n = comp, comp+ncomp-1
        do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          s(i,j,k,n) = s(i,j,k,n) - s0(k,n)
        end do
        end do
        end do
        end do
      else
        do n = comp, comp+ncomp-1
        do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          s(i,j,k,n) = s(i,j,k,n) + s0(k,n)
        end do
        end do
        end do
        end do
      end if

  end subroutine pert_form_3d_cart

  subroutine pert_form_3d_sphr(s,s0,lo,hi,ng,dx,comp,ncomp,flag)

      integer        , intent(in   ) ::  lo(:),hi(:),ng,comp,ncomp
      real(kind=dp_t), intent(inout) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real(kind=dp_t), intent(in   ) :: s0(0:,:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      logical        , intent(in   ) :: flag

      real(kind=dp_t), allocatable :: s0_cart(:,:,:)
      integer                      :: i,j,k,n

      allocate(s0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

      if (flag) then
        do n = comp,comp+ncomp-1
          call fill_3d_data(s0_cart,s0(0:,n),lo,hi,dx,0)
          do k = lo(3),hi(3)
          do j = lo(2),hi(2)
          do i = lo(1),hi(1)
            s(i,j,k,n) = s(i,j,k,n) - s0_cart(i,j,k)
          end do
          end do
          end do
        end do
      else
        do n = comp,comp+ncomp-1
          call fill_3d_data(s0_cart,s0(0:,n),lo,hi,dx,0)
          do k = lo(3),hi(3)
          do j = lo(2),hi(2)
          do i = lo(1),hi(1)
            s(i,j,k,n) = s(i,j,k,n) + s0_cart(i,j,k)
          end do
          end do
          end do
        end do
      end if

      deallocate(s0_cart)

  end subroutine pert_form_3d_sphr

end module pert_form_module
