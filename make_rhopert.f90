module rhopert_module

  use bl_types
  use bl_constants_module
  use variables
  use multifab_module
  use setbc_module
  use define_bc_module

  implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine make_rhopert(rhopert,s,rho0)

      type(multifab)  , intent(inout) :: rhopert
      type(multifab)  , intent(in   ) :: s
      real (kind=dp_t), intent(in   ) :: rho0(:)

      real(kind=dp_t), pointer::  sp(:,:,:,:)
      real(kind=dp_t), pointer:: rpp(:,:,:,:)
      integer :: lo(s%dim),hi(s%dim)
      integer :: i,dm,ng

      dm = s%dim
      ng = s%ng

      do i = 1, s%nboxes
         if ( multifab_remote(s, i) ) cycle
         rpp => dataptr(rhopert, i)
          sp => dataptr(s, i)
         lo =  lwb(get_box(s, i))
         hi =  upb(get_box(s, i))
         select case (dm)
            case (2)
              call make_rhopert_2d(rpp(:,:,1,1),sp(:,:,1,rho_comp),rho0,lo,hi,ng)
            case (3)
              call make_rhopert_3d(rpp(:,:,:,1),sp(:,:,:,rho_comp),rho0,lo,hi,ng)
         end select
      end do

   end subroutine make_rhopert

   subroutine make_rhopert_2d (rhopert,rho,rho0,lo,hi,ng)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:), ng
      real (kind=dp_t), intent(  out) :: rhopert(lo(1)   :,lo(2)   :)
      real (kind=dp_t), intent(in   ) ::     rho(lo(1)-ng:,lo(2)-ng:)
      real (kind=dp_t), intent(in   ) ::    rho0(0:)

!     Local variables
      integer          :: i, j

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        rhopert(i,j) = rho(i,j) - rho0(j)
      end do
      end do
 
   end subroutine make_rhopert_2d

   subroutine make_rhopert_3d (rhopert,rho,rho0,lo,hi,ng)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:), ng
      real (kind=dp_t), intent(  out) :: rhopert(lo(1)   :,lo(2)   :,lo(3)   :)
      real (kind=dp_t), intent(in   ) ::     rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
      real (kind=dp_t), intent(in   ) ::    rho0(0:)

!     Local variables
      integer          :: i, j, k

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        rhopert(i,j,k) = rho(i,j,k) - rho0(k)
      end do
      end do
      end do
 
   end subroutine make_rhopert_3d

end module rhopert_module
