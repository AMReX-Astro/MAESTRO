module XfromrhoX_module

  use bl_types
  use bc_module
  use multifab_module
  use variables
  use network

  implicit none

contains

   subroutine make_XfromrhoX (plotdata,comp,s)

      integer        , intent(in   ) :: comp
      type(multifab) , intent(in   ) :: s
      type(multifab) , intent(inout) :: plotdata

      real(kind=dp_t), pointer:: sp(:,:,:,:)
      real(kind=dp_t), pointer:: pp(:,:,:,:)
      integer :: lo(s%dim),hi(s%dim),ng,dm
      integer :: i

      ng = s%ng
      dm = s%dim

      do i = 1, s%nboxes
         if ( multifab_remote(s, i) ) cycle
         sp => dataptr(s, i)
         pp => dataptr(plotdata, i)
         lo =  lwb(get_box(s, i))
         hi =  upb(get_box(s, i))
         select case (dm)
            case (2)
              call makeXfromrhoX_2d(pp(:,:,1,comp:),sp(:,:,1,1),sp(:,:,1,spec_comp:), lo, hi, ng)
            case (3)
              call makeXfromrhoX_3d(pp(:,:,:,comp:),sp(:,:,:,1),sp(:,:,:,spec_comp:), lo, hi, ng)
         end select
      end do

   end subroutine make_XfromrhoX

   subroutine makeXfromrhoX_2d (X,rho,rhoX,lo,hi,ng)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) ::    X(lo(1)   :,lo(2)   :,:)  
      real (kind = dp_t), intent(in   ) :: rhoX(lo(1)-ng:,lo(2)-ng:,:)
      real (kind = dp_t), intent(in   ) ::  rho(lo(1)-ng:,lo(2)-ng:  )

!     Local variables
      integer :: i, j, n

      do n = 1, nspec
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
         X(i,j,n) = rhoX(i,j,n) / rho(i,j)
      enddo
      enddo
      enddo

   end subroutine makeXfromrhoX_2d

   subroutine makeXfromrhoX_3d (X,rho,rhoX,lo,hi,ng)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) ::    X(lo(1)   :,lo(2)   :,lo(3)   :,:)  
      real (kind = dp_t), intent(in   ) :: rhoX(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(in   ) ::  rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:  )

!     Local variables
      integer :: i, j, k, n

      do n = 1, nspec
      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
         X(i,j,k,n) = rhoX(i,j,k,n) / rho(i,j,k)
      enddo
      enddo
      enddo
      enddo

   end subroutine makeXfromrhoX_3d

end module XfromrhoX_module
