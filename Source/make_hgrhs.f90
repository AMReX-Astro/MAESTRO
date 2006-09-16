module hgrhs_module

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_module
  use macrhs_module

  implicit none

contains


!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine make_hgrhs (hgrhs,macrhs,s,u,div_coeff,p0,t0,gam1,dx,time)

      type(multifab) , intent(inout) :: hgrhs
      type(multifab) , intent(inout) :: macrhs
      type(multifab) , intent(in   ) :: s,u
      real(kind=dp_t), intent(in   ) :: div_coeff(:),p0(:),t0(:),gam1(:)
      real(kind=dp_t), intent(in   ) :: dx(:), time

      real(kind=dp_t), pointer:: mp(:,:,:,:),cp(:,:,:,:)
      real(kind=dp_t), pointer:: sp(:,:,:,:),up(:,:,:,:)
      integer :: lo(s%dim),hi(s%dim),ng,dm
      integer :: i

      ng = s%ng
      dm = s%dim

      do i = 1, s%nboxes
         if ( multifab_remote(s, i) ) cycle
         mp => dataptr(hgrhs, i)
         cp => dataptr(macrhs, i)
         sp => dataptr(s, i)
         up => dataptr(u, i)
         lo =  lwb(get_box(s, i))
         hi =  upb(get_box(s, i))
         select case (dm)
            case (2)
              call make_hgrhs_2d(lo,hi,mp(:,:,1,1),cp(:,:,1,1),sp(:,:,1,:),up(:,:,1,:), &
                                 ng, div_coeff, p0, t0, gam1, dx, time)
            case (3)
              call make_hgrhs_3d(lo,hi,mp(:,:,:,1),cp(:,:,:,1),sp(:,:,:,:),up(:,:,:,:), &
                                 ng, div_coeff, p0, t0, gam1, dx, time)
         end select
      end do

   end subroutine make_hgrhs

   subroutine make_hgrhs_2d (lo,hi,rhs,rhs_cc,s,u,ng,div_coeff,p0,t0,gam1,dx,time)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:), ng
      real (kind=dp_t), intent(  out) ::    rhs(lo(1):,lo(2):)  
      real (kind=dp_t), intent(  out) :: rhs_cc(lo(1):,lo(2):)  
      real (kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,:)
      real (kind=dp_t), intent(in   ) :: u(lo(1)-ng:,lo(2)-ng:,:)
      real (kind=dp_t), intent(in   ) :: div_coeff(lo(2):)
      real (kind=dp_t), intent(in   ) ::        p0(lo(2):)
      real (kind=dp_t), intent(in   ) ::        t0(lo(2):)
      real (kind=dp_t), intent(in   ) ::      gam1(lo(2):)
      real (kind=dp_t), intent(in   ) :: dx(:),time

!     Local variables
      integer :: i, j

      call make_macrhs_2d(lo,hi,rhs_cc,s,u,ng,div_coeff,p0,t0,gam1,dx,time)

!     HACK : THIS ASSUMES EXTRAP AT SIDES AND NO HEATING NEAR TOP OR BOTTOM!!!
      rhs(:,lo(2)  ) = ZERO
      rhs(:,hi(2)  ) = ZERO
      rhs(:,hi(2)+1) = ZERO
      do j = lo(2)+1, hi(2)-1
        rhs(lo(1)  ,j) = HALF * ( rhs_cc(lo(1),j) + rhs_cc(lo(1),j-1) )
        rhs(hi(1)+1,j) = HALF * ( rhs_cc(hi(1),j) + rhs_cc(hi(1),j-1) )
        do i = lo(1)+1, hi(1)
          rhs(i,j) = FOURTH * ( rhs_cc(i,j  ) + rhs_cc(i-1,j  ) &
                               +rhs_cc(i,j-1) + rhs_cc(i-1,j-1) )
        enddo
      enddo

   end subroutine make_hgrhs_2d

   subroutine make_hgrhs_3d (lo,hi,rhs,rhs_cc,s,u,ng,div_coeff,p0,t0,gam1,dx,time)

      implicit none

      integer         , intent(in   ) :: lo(:), hi(:), ng
      real (kind=dp_t), intent(  out) ::    rhs(lo(1):,lo(2):,lo(3):)  
      real (kind=dp_t), intent(  out) :: rhs_cc(lo(1):,lo(2):,lo(3):)  
      real (kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind=dp_t), intent(in   ) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind=dp_t), intent(in   ) :: div_coeff(lo(3):)
      real (kind=dp_t), intent(in   ) ::        p0(lo(3):)
      real (kind=dp_t), intent(in   ) ::        t0(lo(3):)
      real (kind=dp_t), intent(in   ) ::      gam1(lo(3):)
      real (kind=dp_t), intent(in   ) :: dx(:),time

!     Local variables
      integer :: i, j, k

      call make_macrhs_3d(lo,hi,rhs_cc,s,u,ng,div_coeff,p0,t0,gam1,dx,time)

!     HACK : THIS ASSUMES EXTRAP AT SIDES AND NO HEATING NEAR TOP OR BOTTOM!!!
      rhs(:,:,lo(3)  ) = ZERO
      rhs(:,:,hi(3)  ) = ZERO
      rhs(:,:,hi(3)+1) = ZERO
      do k = lo(3)+1, hi(3)-1
        do j = lo(2)+1, hi(2)
          rhs(lo(1)  ,j,k) = FOURTH * ( rhs_cc(lo(1),j,k-1) + rhs_cc(lo(1),j-1,k-1) &
                                       +rhs_cc(lo(1),j,k  ) + rhs_cc(lo(1),j-1,k  ) )
          rhs(hi(1)+1,j,k) = FOURTH * ( rhs_cc(hi(1),j,k-1) + rhs_cc(hi(1),j-1,k-1) &
                                       +rhs_cc(hi(1),j,k  ) + rhs_cc(hi(1),j-1,k  ) )
        enddo
        do i = lo(1)+1, hi(1)
          rhs(i,lo(2)  ,k) = FOURTH * ( rhs_cc(i,lo(2),k-1) + rhs_cc(i-1,lo(2),k-1) &
                                       +rhs_cc(i,lo(2),k  ) + rhs_cc(i-1,lo(2),k  ) )
          rhs(i,hi(2)+1,k) = FOURTH * ( rhs_cc(i,hi(2),k-1) + rhs_cc(i-1,hi(2),k-1) &
                                       +rhs_cc(i,hi(2),k  ) + rhs_cc(i-1,hi(2),k  ) )
        enddo
        rhs(lo(1)  ,lo(2)  ,k) = HALF * (rhs_cc(lo(1),lo(2),k-1) + rhs_cc(lo(1),lo(2),k) )
        rhs(hi(1)+1,lo(2)  ,k) = HALF * (rhs_cc(hi(1),lo(2),k-1) + rhs_cc(hi(1),lo(2),k) )
        rhs(lo(1)  ,hi(2)+1,k) = HALF * (rhs_cc(lo(1),hi(2),k-1) + rhs_cc(lo(1),hi(2),k) )
        rhs(hi(1)+1,hi(2)+1,k) = HALF * (rhs_cc(hi(1),hi(2),k-1) + rhs_cc(hi(1),hi(2),k) )
        do j = lo(2)+1, hi(2)
        do i = lo(1)+1, hi(1)
          rhs(i,j,k) = EIGHTH * ( rhs_cc(i,j  ,k-1) + rhs_cc(i-1,j  ,k-1) &
                                 +rhs_cc(i,j-1,k-1) + rhs_cc(i-1,j-1,k-1) &
                                 +rhs_cc(i,j-1,k  ) + rhs_cc(i-1,j-1,k  ) &
                                 +rhs_cc(i,j-1,k  ) + rhs_cc(i-1,j-1,k  ) )
        enddo
        enddo
      enddo

   end subroutine make_hgrhs_3d

end module hgrhs_module
