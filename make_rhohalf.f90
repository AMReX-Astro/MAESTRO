module rhohalf_module

  use bl_types
  use bl_constants_module
  use bc_module
  use variables
  use multifab_module

  implicit none

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine make_rhohalf (rhohalf,sold,snew)

      type(multifab) , intent(inout) :: rhohalf
      type(multifab) , intent(in   ) :: sold
      type(multifab) , intent(in   ) :: snew

      real(kind=dp_t), pointer:: rhp(:,:,:,:)
      real(kind=dp_t), pointer:: rop(:,:,:,:)
      real(kind=dp_t), pointer:: rnp(:,:,:,:)
      integer :: lo(rhohalf%dim),hi(rhohalf%dim),ng_h,ng_o,dm
      integer :: i

      dm = rhohalf%dim
      ng_h = rhohalf%ng
      ng_o = sold%ng

      do i = 1, rhohalf%nboxes
         if ( multifab_remote(rhohalf, i) ) cycle
         rhp => dataptr(rhohalf, i)
         rop => dataptr(sold, i)
         rnp => dataptr(snew, i)
         lo =  lwb(get_box(rhohalf, i))
         hi =  upb(get_box(rhohalf, i))
         select case (dm)
            case (2)
              call make_rhohalf_2d(rhp(:,:,1,1),rop(:,:,1,rho_comp),rnp(:,:,1,rho_comp),&
                                   lo,hi,ng_h,ng_o)
            case (3)
              call make_rhohalf_3d(rhp(:,:,:,1),rop(:,:,:,rho_comp),rnp(:,:,:,rho_comp),&
                                   lo,hi,ng_h,ng_o)
         end select
      end do

   end subroutine make_rhohalf

   subroutine make_rhohalf_2d (rhohalf,rhoold,rhonew,lo,hi,ng_half,ng_old)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:), ng_half, ng_old
      real (kind=dp_t), intent(  out) :: rhohalf(lo(1)-ng_half:,lo(2)-ng_half:)
      real (kind=dp_t), intent(in   ) ::  rhoold(lo(1)-ng_old :,lo(2)-ng_old :)
      real (kind=dp_t), intent(in   ) ::  rhonew(lo(1)-ng_old :,lo(2)-ng_old :)

!     Local variables
      integer          :: i, j

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        rhohalf(i,j) = HALF * (rhoold(i,j) + rhonew(i,j))
      end do
      end do
 
   end subroutine make_rhohalf_2d

   subroutine make_rhohalf_3d (rhohalf,rhoold,rhonew,lo,hi,ng_half,ng_old)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:), ng_half, ng_old
      real (kind=dp_t), intent(  out) :: rhohalf(lo(1)-ng_half:,lo(2)-ng_half:,lo(3)-ng_half:)
      real (kind=dp_t), intent(in   ) ::  rhoold(lo(1)-ng_old:,lo(2)-ng_old:,lo(3)-ng_old:)
      real (kind=dp_t), intent(in   ) ::  rhonew(lo(1)-ng_old:,lo(2)-ng_old:,lo(3)-ng_old:)

!     Local variables
      integer          :: i, j, k

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        rhohalf(i,j,k) = HALF * (rhoold(i,j,k) + rhonew(i,j,k))
      end do
      end do
      end do
 
   end subroutine make_rhohalf_3d

end module rhohalf_module
