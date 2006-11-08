module macrhs_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use eos_module
  use network
  use variables

  implicit none

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine make_macrhs (macrhs,Source,Sbar,div_coeff)

      type(multifab) , intent(inout) :: macrhs
      type(multifab) , intent(in   ) :: Source
      real(kind=dp_t), intent(in   ) :: Sbar(:)
      real(kind=dp_t), intent(in   ) :: div_coeff(:)

      real(kind=dp_t), pointer:: mp(:,:,:,:),sp(:,:,:,:)
      integer :: lo(Source%dim),hi(Source%dim),dm
      integer :: i

      dm = Source%dim

      do i = 1, Source%nboxes
         if ( multifab_remote(Source, i) ) cycle
         mp => dataptr(macrhs, i)
         sp => dataptr(Source, i)
         lo =  lwb(get_box(Source, i))
         hi =  upb(get_box(Source, i))
         select case (dm)
            case (2)
              call make_macrhs_2d(lo,hi,mp(:,:,1,1),sp(:,:,1,1),Sbar,div_coeff)
            case (3)
              call make_macrhs_3d(lo,hi,mp(:,:,:,1),sp(:,:,:,1),Sbar,div_coeff)
         end select
      end do

   end subroutine make_macrhs

   subroutine make_macrhs_2d (lo,hi,rhs,Source,Sbar,div_coeff)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:)
      real (kind=dp_t), intent(  out) :: rhs(lo(1):,lo(2):)  
      real (kind=dp_t), intent(in   ) :: Source(lo(1):,lo(2):)  
      real (kind=dp_t), intent(in   ) :: Sbar(lo(2):)  
      real (kind=dp_t), intent(in   ) :: div_coeff(lo(2):)  

!     Local variables
      integer :: i, j

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        rhs(i,j) = div_coeff(j) * (Source(i,j) - Sbar(j))
      end do
      end do
 
   end subroutine make_macrhs_2d

   subroutine make_macrhs_3d (lo,hi,rhs,Source,Sbar,div_coeff)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:)
      real (kind=dp_t), intent(  out) :: rhs(lo(1):,lo(2):,lo(3):)  
      real (kind=dp_t), intent(in   ) :: Source(lo(1):,lo(2):,lo(3):)  
      real (kind=dp_t), intent(in   ) :: Sbar(lo(3):)  
      real (kind=dp_t), intent(in   ) :: div_coeff(lo(3):)  

!     Local variables
      integer :: i, j, k

      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        rhs(i,j,k) = div_coeff(k) * (Source(i,j,k) - Sbar(k))
      end do
      end do
      end do
 
   end subroutine make_macrhs_3d

end module macrhs_module
