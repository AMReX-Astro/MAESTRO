module macrhs_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use network
  use variables
  use geometry
  use fill_3d_module

  implicit none

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine make_macrhs (macrhs,Source,gamma1_term,Sbar,div_coeff,dx)

      type(multifab) , intent(inout) :: macrhs
      type(multifab) , intent(in   ) :: Source
      type(multifab) , intent(in   ) :: gamma1_term
      real(kind=dp_t), intent(in   ) :: Sbar(:)
      real(kind=dp_t), intent(in   ) :: div_coeff(:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      real(kind=dp_t), pointer:: mp(:,:,:,:),sp(:,:,:,:),gp(:,:,:,:)
      integer :: lo(Source%dim),hi(Source%dim),dm
      integer :: i

      dm = Source%dim

      do i = 1, Source%nboxes
         if ( multifab_remote(Source, i) ) cycle
         mp => dataptr(macrhs, i)
         sp => dataptr(Source, i)
         gp => dataptr(gamma1_term, i)
         lo =  lwb(get_box(Source, i))
         hi =  upb(get_box(Source, i))
         select case (dm)
            case (2)
              call make_macrhs_2d(lo,hi,mp(:,:,1,1),sp(:,:,1,1),gp(:,:,1,1),Sbar,div_coeff,dx)
            case (3)
              call make_macrhs_3d(lo,hi,mp(:,:,:,1),sp(:,:,:,1),gp(:,:,:,1),Sbar,div_coeff,dx)
         end select
      end do

   end subroutine make_macrhs

   subroutine make_macrhs_2d (lo,hi,rhs,Source,gamma1_term,Sbar,div_coeff,dx)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:)
      real (kind=dp_t), intent(  out) :: rhs(lo(1):,lo(2):)  
      real (kind=dp_t), intent(in   ) :: Source(lo(1):,lo(2):)  
      real (kind=dp_t), intent(in   ) :: gamma1_term(lo(1):,lo(2):)  
      real (kind=dp_t), intent(in   ) :: Sbar(lo(2):)  
      real (kind=dp_t), intent(in   ) :: div_coeff(lo(2):)  
      real (kind=dp_t), intent(in   ) :: dx(:)

!     Local variables
      integer :: i, j
      integer :: rr

      rr = int( dx(2) / dr + 1.d-12)

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        rhs(i,j) = div_coeff(j) * (Source(i,j) - Sbar(rr*j) + gamma1_term(i,j))
      end do
      end do
 
   end subroutine make_macrhs_2d

   subroutine make_macrhs_3d (lo,hi,rhs,Source,gamma1_term,Sbar,div_coeff,dx)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:)
      real (kind=dp_t), intent(  out) :: rhs(lo(1):,lo(2):,lo(3):)  
      real (kind=dp_t), intent(in   ) :: Source(lo(1):,lo(2):,lo(3):)  
      real (kind=dp_t), intent(in   ) :: gamma1_term(lo(1):,lo(2):,lo(3):)  
      real (kind=dp_t), intent(in   ) :: Sbar(lo(3):)  
      real (kind=dp_t), intent(in   ) :: div_coeff(lo(3):)  
      real (kind=dp_t), intent(in   ) :: dx(:)

!     Local variables
      integer :: i, j, k
      integer :: rr
      real (kind=dp_t), allocatable :: div_cart(:,:,:),Sbar_cart(:,:,:)

      if (spherical .eq. 1) then

        allocate(div_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call fill_3d_data(div_cart,div_coeff,lo,hi,dx,0)

        allocate(Sbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call fill_3d_data(Sbar_cart,Sbar,lo,hi,dx,0)

        do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          rhs(i,j,k) = div_cart(i,j,k) * (Source(i,j,k) - Sbar_cart(i,j,k))
        end do
        end do
        end do

        deallocate(Sbar_cart,div_cart)

      else

        rr = int( dx(2) / dr + 1.d-12)

        do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          rhs(i,j,k) = div_coeff(k) * (Source(i,j,k) - Sbar(rr*k))
        end do
        end do
        end do
      end if
 
   end subroutine make_macrhs_3d

end module macrhs_module
