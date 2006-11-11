module phihalf_module

  use bl_types
  use bl_constants_module
  use variables
  use multifab_module
  use setbc_module
  use define_bc_module

  implicit none

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine make_at_halftime (phihalf,sold,snew,in_comp,out_comp,dx,the_bc_level)

      type(multifab) , intent(inout) :: phihalf
      type(multifab) , intent(in   ) :: sold
      type(multifab) , intent(in   ) :: snew
      integer        , intent(in   ) :: in_comp,out_comp
      real(kind=dp_t), intent(in   ) :: dx(:)
      type(bc_level) , intent(in   ) :: the_bc_level

      real(kind=dp_t), pointer:: rhp(:,:,:,:)
      real(kind=dp_t), pointer:: rop(:,:,:,:)
      real(kind=dp_t), pointer:: rnp(:,:,:,:)
      integer :: lo(phihalf%dim),hi(phihalf%dim),ng_h,ng_o,dm
      integer :: i,bc_comp

      dm = phihalf%dim
      ng_h = phihalf%ng
      ng_o = sold%ng

      do i = 1, phihalf%nboxes
         if ( multifab_remote(phihalf, i) ) cycle
         rhp => dataptr(phihalf, i)
         rop => dataptr(sold, i)
         rnp => dataptr(snew, i)
         lo =  lwb(get_box(phihalf, i))
         hi =  upb(get_box(phihalf, i))
         select case (dm)
            case (2)
              call make_at_halftime_2d(rhp(:,:,1,out_comp),rop(:,:,1,in_comp),rnp(:,:,1,in_comp),&
                                   lo,hi,ng_h,ng_o)
              if (ng_h .gt. 0) then
                bc_comp = dm+in_comp
                call setbc_2d(rhp(:,:,1,out_comp), lo, ng_h, &
                              the_bc_level%adv_bc_level_array(i,:,:,bc_comp),dx,bc_comp)
              end if
            case (3)
              call make_at_halftime_3d(rhp(:,:,:,out_comp),rop(:,:,:,in_comp),rnp(:,:,:,in_comp),&
                                   lo,hi,ng_h,ng_o)
              if (ng_h .gt. 0) then
                bc_comp = dm+in_comp
                call setbc_3d(rhp(:,:,:,out_comp), lo, ng_h, &
                              the_bc_level%adv_bc_level_array(i,:,:,bc_comp),dx,bc_comp)
              end if
         end select
      end do

   end subroutine make_at_halftime

   subroutine make_at_halftime_2d (phihalf,phiold,phinew,lo,hi,ng_half,ng_old)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:), ng_half, ng_old
      real (kind=dp_t), intent(  out) :: phihalf(lo(1)-ng_half:,lo(2)-ng_half:)
      real (kind=dp_t), intent(in   ) ::  phiold(lo(1)-ng_old :,lo(2)-ng_old :)
      real (kind=dp_t), intent(in   ) ::  phinew(lo(1)-ng_old :,lo(2)-ng_old :)

!     Local variables
      integer          :: i, j

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        phihalf(i,j) = HALF * (phiold(i,j) + phinew(i,j))
      end do
      end do
 
   end subroutine make_at_halftime_2d

   subroutine make_at_halftime_3d (phihalf,phiold,phinew,lo,hi,ng_half,ng_old)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:), ng_half, ng_old
      real (kind=dp_t), intent(  out) :: phihalf(lo(1)-ng_half:,lo(2)-ng_half:,lo(3)-ng_half:)
      real (kind=dp_t), intent(in   ) ::  phiold(lo(1)-ng_old:,lo(2)-ng_old:,lo(3)-ng_old:)
      real (kind=dp_t), intent(in   ) ::  phinew(lo(1)-ng_old:,lo(2)-ng_old:,lo(3)-ng_old:)

!     Local variables
      integer          :: i, j, k

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        phihalf(i,j,k) = HALF * (phiold(i,j,k) + phinew(i,j,k))
      end do
      end do
      end do
 
   end subroutine make_at_halftime_3d

end module phihalf_module
