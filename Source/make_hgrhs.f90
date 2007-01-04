module hgrhs_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use geometry
  use fill_3d_module

  implicit none

contains


!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine make_hgrhs (hgrhs,Source,gamma1_term,Sbar,div_coeff,dx,domain)

      type(multifab) , intent(inout) :: hgrhs
      type(multifab) , intent(in   ) :: Source
      type(multifab) , intent(in   ) :: gamma1_term
      real(kind=dp_t), intent(in   ) :: Sbar(0:)
      real(kind=dp_t), intent(in   ) :: div_coeff(0:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      type(box      ), intent(in   ) :: domain

      type(multifab)  :: rhs_cc
      real(kind=dp_t), pointer:: hp(:,:,:,:),sp(:,:,:,:),gp(:,:,:,:),rp(:,:,:,:)
      integer :: lo(Source%dim),hi(Source%dim)
      integer :: i,dm,lo_z,hi_z

      dm = Source%dim
      lo_z = domain%lo(dm)
      hi_z = domain%hi(dm)

      call multifab_build(rhs_cc,Source%la,1,1)
      call setval(rhs_cc,ZERO,all=.true.)

      do i = 1, Source%nboxes
         if ( multifab_remote(Source, i) ) cycle
         rp => dataptr(rhs_cc, i)
         sp => dataptr(Source, i)
         gp => dataptr(gamma1_term, i)
         lo =  lwb(get_box(Source, i))
         hi =  upb(get_box(Source, i))
         select case (dm)
            case (2)
              call make_rhscc_2d(lo,hi,rp(:,:,1,1),sp(:,:,1,1),gp(:,:,1,1),Sbar,div_coeff)
            case (3)
              call make_rhscc_3d(lo,hi,rp(:,:,:,1),sp(:,:,:,1),gp(:,:,:,1),Sbar,div_coeff,dx)
         end select
      end do
      call multifab_fill_boundary(rhs_cc)

      call setval(hgrhs,ZERO,all=.true.)
      do i = 1, Source%nboxes
         if ( multifab_remote(Source, i) ) cycle
         hp => dataptr(hgrhs, i)
         rp => dataptr(rhs_cc, i)
         lo =  lwb(get_box(Source, i))
         hi =  upb(get_box(Source, i))
         select case (dm)
            case (2)
              call make_hgrhs_2d(lo,hi,lo_z,hi_z,hp(:,:,1,1),rp(:,:,1,1))
            case (3)
              call make_hgrhs_3d(lo,hi,lo_z,hi_z,hp(:,:,:,1),rp(:,:,:,1))
         end select
      end do

      call multifab_destroy(rhs_cc)

   end subroutine make_hgrhs

   subroutine make_rhscc_2d(lo,hi,rhs_cc,Source,gamma1_term,Sbar,div_coeff)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:)
      real (kind=dp_t), intent(  out) :: rhs_cc(lo(1)-1:,lo(2)-1:)
      real (kind=dp_t), intent(in   ) :: Source(lo(1):,lo(2):)
      real (kind=dp_t), intent(in   ) :: gamma1_term(lo(1):,lo(2):)  
      real (kind=dp_t), intent(in   ) ::      Sbar(0:)
      real (kind=dp_t), intent(in   ) :: div_coeff(0:)

!     Local variables
      integer :: i, j

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        rhs_cc(i,j) = div_coeff(j) * (Source(i,j) - Sbar(j) + gamma1_term(i,j))
      end do
      end do

   end subroutine make_rhscc_2d

   subroutine make_hgrhs_2d(lo,hi,lo_z,hi_z,rhs,rhs_cc)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:), lo_z, hi_z
      real (kind=dp_t), intent(  out) :: rhs(lo(1):,lo(2):)  
      real (kind=dp_t), intent(in   ) :: rhs_cc(lo(1)-1:,lo(2)-1:)

!     Local variables
      integer :: i, j

      do j = lo(2),hi(2)+1
        if (j.ne.lo_z .and. j.ne.(hi_z+1)) then
          do i = lo(1), hi(1)+1
            rhs(i,j) = FOURTH * ( rhs_cc(i,j  ) + rhs_cc(i-1,j  ) &
                                 +rhs_cc(i,j-1) + rhs_cc(i-1,j-1) )
          enddo
        end if
      enddo

   end subroutine make_hgrhs_2d

   subroutine make_rhscc_3d(lo,hi,rhs_cc,Source,gamma1_term,Sbar,div_coeff,dx)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:)
      real (kind=dp_t), intent(  out) :: rhs_cc(lo(1)-1:,lo(2)-1:,lo(3)-1:)  
      real (kind=dp_t), intent(in   ) :: Source(lo(1):,lo(2):,lo(3):)  
      real (kind=dp_t), intent(in   ) :: gamma1_term(lo(1):,lo(2):,lo(3):)  
      real (kind=dp_t), intent(in   ) :: Sbar(0:)
      real (kind=dp_t), intent(in   ) :: div_coeff(0:)
      real (kind=dp_t), intent(in   ) :: dx(:)

!     Local variables
      integer :: i, j,k
      real (kind=dp_t), allocatable :: Sbar_cart(:,:,:),div_cart(:,:,:)

     if (spherical .eq. 1) then
 
        allocate(div_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call fill_3d_data(div_cart,div_coeff,lo,hi,dx,0)
 
        allocate(Sbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call fill_3d_data(Sbar_cart,Sbar,lo,hi,dx,0)
 
        do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          rhs_cc(i,j,k) = div_cart(i,j,k) * (Source(i,j,k) - Sbar_cart(i,j,k))
        end do
        end do  
        end do
       
        deallocate(Sbar_cart,div_cart) 
       
      else  
       
        do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          rhs_cc(i,j,k) = div_coeff(k) * (Source(i,j,k) - Sbar(k))
        end do 
        end do 
        end do 

      end if

   end subroutine make_rhscc_3d

   subroutine make_hgrhs_3d(lo,hi,lo_z,hi_z,rhs,rhs_cc)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:), lo_z, hi_z
      real (kind=dp_t), intent(  out) ::    rhs(lo(1)  :,lo(2)  :,lo(3)  :)  
      real (kind=dp_t), intent(in   ) :: rhs_cc(lo(1)-1:,lo(2)-1:,lo(3)-1:)  

!     Local variables
      integer :: i, j,k
 
        do k = lo(3), hi(3)+1
          if (k.ne.lo_z .and. k.ne.hi_z+1) then
            do j = lo(2), hi(2)+1
            do i = lo(1), hi(1)+1
              rhs(i,j,k) = EIGHTH * ( rhs_cc(i,j  ,k-1) + rhs_cc(i-1,j  ,k-1) &
                                     +rhs_cc(i,j-1,k-1) + rhs_cc(i-1,j-1,k-1) &
                                     +rhs_cc(i,j-1,k  ) + rhs_cc(i-1,j-1,k  ) &
                                     +rhs_cc(i,j-1,k  ) + rhs_cc(i-1,j-1,k  ) )
            enddo
            enddo
          end if
        enddo

   end subroutine make_hgrhs_3d

end module hgrhs_module
