module enthalpy_module

  use bl_types
  use bc_module
  use multifab_module

  implicit none

contains

   subroutine make_enthalpy (enthalpy,comp,s)

      integer        , intent(in   ) :: comp
      type(multifab) , intent(inout) :: enthalpy
      type(multifab) , intent(in   ) :: s

      real(kind=dp_t), pointer:: sp(:,:,:,:)
      real(kind=dp_t), pointer:: rp(:,:,:,:)
      integer :: lo(s%dim),hi(s%dim),ng,dm
      integer :: i

      ng = s%ng
      dm = s%dim

      do i = 1, s%nboxes
         if ( multifab_remote(s, i) ) cycle
         sp => dataptr(s, i)
         rp => dataptr(enthalpy, i)
         lo =  lwb(get_box(s, i))
         hi =  upb(get_box(s, i))
         select case (dm)
            case (2)
              call make_enthalpy_2d(rp(:,:,1,comp),sp(:,:,1,:), lo, hi, ng)
            case (3)
              call make_enthalpy_3d(rp(:,:,:,comp),sp(:,:,:,:), lo, hi, ng)
         end select
      end do

   end subroutine make_enthalpy

   subroutine make_enthalpy_2d (enthalpy,s,lo,hi,ng)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) :: enthalpy(lo(1):,lo(2):)  
      real (kind = dp_t), intent(in   ) ::    s(lo(1)-ng:,lo(2)-ng:,:)

!     Local variables
      integer :: i, j

      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           enthalpy(i,j) = s(i,j,2) / s(i,j,1)
        enddo
      enddo

   end subroutine make_enthalpy_2d

   subroutine make_enthalpy_3d (enthalpy,s,lo,hi,ng)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) :: enthalpy(lo(1):,lo(2):, lo(3):)  
      real (kind = dp_t), intent(in   ) ::    s(lo(1)-ng:,lo(2)-ng:, lo(3)-ng:,:)

!     Local variables
      integer :: i, j, k

      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           enthalpy(i,j,k) = s(i,j,k,2) / s(i,j,k,1)
        enddo
      enddo
      end do

   end subroutine make_enthalpy_3d

end module enthalpy_module
