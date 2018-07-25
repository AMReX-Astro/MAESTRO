module make_gpi_module

  use bl_types
  use mg_module
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private
  public :: make_gpi

contains 

    subroutine make_gpi(gpi,s,dx,mla)
      use variables                   , only : pi_comp
      
      type(multifab), intent(inout) :: gpi(:)
      type(multifab), intent(in   ) :: s(:)
      real(dp_t)    , intent(in   ) :: dx(:,:)
      type(ml_layout), intent(in   ) :: mla

      integer :: lo(gpi(1)%dim),hi(gpi(1)%dim)
      integer :: i,n,ng_s,ng_gp,dm,nlevs

      real(kind=dp_t), pointer :: gp(:,:,:,:) 
      real(kind=dp_t), pointer :: sp(:,:,:,:) 

      type(bl_prof_timer), save :: bpt

      call build(bpt, "make_gi")

      ng_s  = nghost(s(1))
      ng_gp = nghost(gpi(1))
      dm = mla%dim
      nlevs = mla%nlevel  

      do n = 1, nlevs

         do i = 1, nfabs(s(n))
            lo = lwb(get_box(gpi(n),i))
            hi = upb(get_box(gpi(n),i))
            gp => dataptr(gpi(n),i)
            sp  => dataptr(s(n),i)
            select case (dm)
            case (1)
               call mkgpi_1d(gp(:,1,1,1), ng_gp, sp(:,1,1,pi_comp), ng_s, &
                              lo, hi, dx(n,:))
            case (2)
               call mkgpi_2d(gp(:,:,1,:), ng_gp, sp(:,:,1,pi_comp), ng_s, &
                              lo, hi, dx(n,:))
            case (3)
               call mkgpi_3d(gp(:,:,:,:), ng_gp, sp(:,:,:,pi_comp), ng_s, &
                              lo, hi, dx(n,:))
            end select
         end do

      end do

      call destroy(bpt)

    end subroutine make_gpi

    !   ********************************************************************************* !

    subroutine mkgpi_1d(gpi,ng_gp,pi,ng_p,lo,hi,dx)

      integer        , intent(in   ) :: ng_gp,ng_p
      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(inout) ::  gpi(lo(1)-ng_gp:)
      real(kind=dp_t), intent(inout) ::   pi(lo(1)-ng_p :)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i

      do i = lo(1),hi(1)
         gpi(i) = ( pi(i+1) - pi(i-1) ) / (TWO*dx(1))
      end do

    end subroutine mkgpi_1d

    !   ********************************************************************************* !

    subroutine mkgpi_2d(gpi,ng_gp,pi,ng_p,lo,hi,dx)

      integer        , intent(in   ) :: ng_gp,ng_p
      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(inout) ::  gpi(lo(1)-ng_gp:,lo(2)-ng_gp:,:)
      real(kind=dp_t), intent(inout) ::   pi(lo(1)-ng_p :,lo(2)-ng_p :)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            gpi(i,j,1) = (pi(i+1,j) - pi(i-1,j)) /(TWO * dx(1))
            gpi(i,j,2) = (pi(i,j+1) - pi(i,j-1)) /(TWO * dx(2))
         end do
      end do

    end subroutine mkgpi_2d

    !   ******************************************************************************** !

    subroutine mkgpi_3d(gpi,ng_gp,pi,ng_p,lo,hi,dx)

      integer        , intent(in   ) :: ng_gp,ng_p
      integer        , intent(in   ) :: lo(:),hi(:)
      real(kind=dp_t), intent(inout) ::  gpi(lo(1)-ng_gp:,lo(2)-ng_gp:,lo(3)-ng_gp:,:)
      real(kind=dp_t), intent(inout) ::   pi(lo(1)-ng_p :,lo(2)-ng_p :,lo(3)-ng_p :)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j,k

      !$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               gpi(i,j,k,1) = (pi(i+1,j,k  ) - pi(i-1,j,k  )) /(TWO*dx(1))
               gpi(i,j,k,2) = (pi(i,j+1,k  ) - pi(i,j-1,k  )) /(TWO*dx(2))
               gpi(i,j,k,3) = (pi(i,j  ,k+1) - pi(i,j  ,k-1)) /(TWO*dx(3))
            end do
         end do
      end do
      !$OMP END PARALLEL DO

    end subroutine mkgpi_3d
    
end module make_gpi_module
