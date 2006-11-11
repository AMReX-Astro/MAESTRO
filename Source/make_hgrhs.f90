module hgrhs_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use geometry
  use fill_3d_module

  implicit none

contains


!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine make_hgrhs (hgrhs,Source,Sbar,div_coeff,dx)

      type(multifab) , intent(inout) :: hgrhs
      type(multifab) , intent(in   ) :: Source
      real(kind=dp_t), intent(in   ) :: Sbar(:)
      real(kind=dp_t), intent(in   ) :: div_coeff(:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      real(kind=dp_t), pointer:: hp(:,:,:,:),sp(:,:,:,:)
      integer :: lo(Source%dim),hi(Source%dim),dm
      integer :: i

      dm = Source%dim

      do i = 1, Source%nboxes
         if ( multifab_remote(Source, i) ) cycle
         hp => dataptr(hgrhs, i)
         sp => dataptr(Source, i)
         lo =  lwb(get_box(Source, i))
         hi =  upb(get_box(Source, i))
         select case (dm)
            case (2)
              call make_hgrhs_2d(lo,hi,hp(:,:,1,1),sp(:,:,1,1),Sbar,div_coeff)
            case (3)
              call make_hgrhs_3d(lo,hi,hp(:,:,:,1),sp(:,:,:,1),Sbar,div_coeff,dx)
         end select
      end do

   end subroutine make_hgrhs

   subroutine make_hgrhs_2d(lo,hi,rhs,Source,Sbar,div_coeff)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:)
      real (kind=dp_t), intent(  out) :: rhs(lo(1):,lo(2):)  
      real (kind=dp_t), intent(in   ) :: Source(lo(1):,lo(2):)  
      real (kind=dp_t), intent(in   ) :: Sbar(lo(2):)
      real (kind=dp_t), intent(in   ) :: div_coeff(lo(2):)

!     Local variables
      integer :: i, j
      real (kind=dp_t), allocatable :: rhs_cc(:,:)

      allocate(rhs_cc(lo(1):hi(1),lo(2):hi(2)))

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        rhs_cc(i,j) = div_coeff(j) * (Source(i,j) - Sbar(j))
      end do
      end do

!     HACK : THIS ASSUMES PERIODIC AT SIDES AND NO HEATING NEAR TOP OR BOTTOM!!!
      rhs(:,lo(2)  ) = ZERO
      rhs(:,hi(2)  ) = ZERO
      rhs(:,hi(2)+1) = ZERO
      do j = lo(2)+1, hi(2)-1
        rhs(lo(1)  ,j) = FOURTH * ( rhs_cc(lo(1),j) + rhs_cc(lo(1),j-1) + &
                                    rhs_cc(hi(1),j) + rhs_cc(hi(1),j-1) )
        rhs(hi(1)+1,j) = rhs(lo(1),j)
        do i = lo(1)+1, hi(1)
          rhs(i,j) = FOURTH * ( rhs_cc(i,j  ) + rhs_cc(i-1,j  ) &
                               +rhs_cc(i,j-1) + rhs_cc(i-1,j-1) )
        enddo
      enddo

      deallocate(rhs_cc)

   end subroutine make_hgrhs_2d

   subroutine make_hgrhs_3d(lo,hi,rhs,Source,Sbar,div_coeff,dx)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:)
      real (kind=dp_t), intent(  out) :: rhs(lo(1):,lo(2):,lo(3):)  
      real (kind=dp_t), intent(in   ) :: Source(lo(1):,lo(2):,lo(3):)  
      real (kind=dp_t), intent(in   ) :: Sbar(lo(3):)
      real (kind=dp_t), intent(in   ) :: div_coeff(lo(3):)
      real (kind=dp_t), intent(in   ) :: dx(:)

!     Local variables
      integer :: i, j,k
      real (kind=dp_t), allocatable :: rhs_cc(:,:,:),Sbar_cart(:,:,:),div_cart(:,:,:)

      allocate(rhs_cc(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

     if (spherical .eq. 1) then
 
        allocate(div_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call fill_3d_data(div_cart,div_coeff,dx,0)
 
        allocate(Sbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        call fill_3d_data(Sbar_cart,Sbar,dx,0)
 
        do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          rhs_cc(i,j,k) = div_cart(i,j,k) * (Source(i,j,k) - Sbar_cart(i,j,k))
        end do
        end do  
        end do
       
        deallocate(Sbar_cart,div_cart) 

        ! HACK : THIS ASSUMES OUTFLOW AT ALL SIDES AND NO HEATING NEAR ANY SIDES
        rhs = ZERO
        do k = lo(3)+1, hi(3)
        do j = lo(2)+1, hi(2)
        do i = lo(1)+1, hi(1)
          rhs(i,j,k) = EIGHTH * ( rhs_cc(i,j  ,k-1) + rhs_cc(i-1,j  ,k-1) &
                                 +rhs_cc(i,j-1,k-1) + rhs_cc(i-1,j-1,k-1) &
                                 +rhs_cc(i,j-1,k  ) + rhs_cc(i-1,j-1,k  ) &
                                 +rhs_cc(i,j-1,k  ) + rhs_cc(i-1,j-1,k  ) )
        enddo
        enddo
        enddo
       
      else  
       
        do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)
          rhs_cc(i,j,k) = div_coeff(k) * (Source(i,j,k) - Sbar(k))
        end do 
        end do 
        end do 

        ! HACK : THIS ASSUMES PERIODIC AT SIDES AND NO HEATING NEAR TOP OR BOTTOM!!!
        rhs(:,:,lo(3)  ) = ZERO
        rhs(:,:,hi(3)  ) = ZERO
        rhs(:,:,hi(3)+1) = ZERO
        do k = lo(3)+1, hi(3)-1
        do j = lo(2)+1, hi(2)
          rhs(lo(1)  ,j,k) = EIGHTH * ( rhs_cc(lo(1),j,k-1) + rhs_cc(lo(1),j-1,k-1) &
                                       +rhs_cc(lo(1),j,k  ) + rhs_cc(lo(1),j-1,k  ) &
                                       +rhs_cc(hi(1),j,k-1) + rhs_cc(hi(1),j-1,k-1) &
                                       +rhs_cc(hi(1),j,k  ) + rhs_cc(hi(1),j-1,k  ) )
          rhs(hi(1)+1,j,k) = rhs(lo(1),j,k)
                            
        enddo
        do i = lo(1)+1, hi(1)
          rhs(i,lo(2)  ,k) = EIGHTH * ( rhs_cc(i,lo(2),k-1) + rhs_cc(i-1,lo(2),k-1) &
                                       +rhs_cc(i,lo(2),k  ) + rhs_cc(i-1,lo(2),k  ) &
                                       +rhs_cc(i,hi(2),k-1) + rhs_cc(i-1,hi(2),k-1) &
                                       +rhs_cc(i,hi(2),k  ) + rhs_cc(i-1,hi(2),k  ) )
          rhs(i,hi(2)+1,k) =  rhs(i,lo(2),k)
                            
        enddo
        rhs(lo(1)  ,lo(2)  ,k) = EIGHTH * ( rhs_cc(lo(1),lo(2),k-1) + rhs_cc(lo(1),lo(2),k) &
                                           +rhs_cc(hi(1),lo(2),k-1) + rhs_cc(hi(1),lo(2),k) &
                                           +rhs_cc(lo(1),hi(2),k-1) + rhs_cc(lo(1),hi(2),k) &
                                           +rhs_cc(hi(1),hi(2),k-1) + rhs_cc(hi(1),hi(2),k) )
        rhs(hi(1)+1,lo(2)  ,k) = rhs(lo(1),lo(2),k)
        rhs(lo(1)  ,hi(2)+1,k) = rhs(lo(1),lo(2),k)
        rhs(hi(1)+1,hi(2)+1,k) = rhs(lo(1),lo(2),k)
        do j = lo(2)+1, hi(2)
        do i = lo(1)+1, hi(1)
          rhs(i,j,k) = EIGHTH * ( rhs_cc(i,j  ,k-1) + rhs_cc(i-1,j  ,k-1) &
                                 +rhs_cc(i,j-1,k-1) + rhs_cc(i-1,j-1,k-1) &
                                 +rhs_cc(i,j-1,k  ) + rhs_cc(i-1,j-1,k  ) &
                                 +rhs_cc(i,j-1,k  ) + rhs_cc(i-1,j-1,k  ) )
        enddo
        enddo
        enddo
      end if

   end subroutine make_hgrhs_3d

end module hgrhs_module
