module rhopert_module

  use bl_types
  use bc_module
  use multifab_module
  use variables

  implicit none

contains

  subroutine make_rhopert (rhopert,comp,s,rho0)

    integer        , intent(in   ) :: comp
    type(multifab) , intent(inout) :: rhopert
    type(multifab) , intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: rho0(:)
    
    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: pp(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng,dm
    integer :: i
    
    ng = s%ng
    dm = s%dim
    
    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       sp => dataptr(s, i)
       pp => dataptr(rhopert, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (2)
          call makerhopert_2d(pp(:,:,1,comp), sp(:,:,1,rho_comp), &
                              lo, hi, ng, rho0)
       case (3)
          call makerhopert_3d(pp(:,:,:,comp), sp(:,:,:,rho_comp), &
                              lo, hi, ng, rho0)
       end select
    end do

  end subroutine make_rhopert

  subroutine makerhopert_2d (rhopert,s,lo,hi,ng,rho0)

    implicit none
    
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: rhopert(lo(1):,lo(2):)  
    real (kind = dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:)
    real (kind = dp_t), intent(in   ) :: rho0(lo(2):)

!     Local variables
    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          rhopert(i,j) = s(i,j) - rho0(j)
       enddo
    enddo
    
  end subroutine makerhopert_2d
  
  subroutine makerhopert_3d (rhopert,s,lo,hi,ng,rho0)

    implicit none
    
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: rhopert(lo(1):,lo(2):, lo(3):)  
    real (kind = dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:, lo(3)-ng:)
    real (kind = dp_t), intent(in   ) :: rho0(lo(3):)
    
    !     Local variables
    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhopert(i,j,k) = s(i,j,k) - rho0(k)
          enddo
       enddo
    end do
    
  end subroutine makerhopert_3d
  
end module rhopert_module
