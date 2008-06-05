! This routine returns the externally imposed (i.e. not reactions)
! heating source term to the enthalpy equation (actually rho * H
! is returned, where H has units of erg/g/s).

module heating_module

  use bl_types

  implicit none
  
  private
  public :: get_rho_Hext
  
contains

  subroutine get_rho_Hext(nlevs,mla,s,rho_Hext,dx,time)

    use multifab_module
    use ml_layout_module
    use ml_restriction_module

    integer, intent(in) :: nlevs
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),time

    ! local
    integer                  :: n,i,ng,dm
    integer                  :: lo(s(1)%dim),hi(s(1)%dim)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: hp(:,:,:,:)

    ng = s(1)%ng
    dm = s(1)%dim

    do n=1,nlevs

       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          sp => dataptr(s(n) , i)
          hp => dataptr(rho_Hext(n) , i)
          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))
          select case (dm)
          case (2)
             call get_rho_Hext_2d(hp(:,:,1,1), sp(:,:,1,:), lo, hi, ng, dx(n,:), time)
          case (3)
             call get_rho_Hext_3d(hp(:,:,:,1), sp(:,:,:,:), lo, hi, ng, dx(n,:), time)
          end select
       end do

    end do

    ! the loop over nlevs must count backwards to make sure the finer grids are done first
    do n=nlevs,2,-1
       ! set level n-1 data to be the average of the level n data covering it
       call ml_cc_restriction(rho_Hext(n-1), rho_Hext(n), mla%mba%rr(n-1,:))
    end do

  end subroutine get_rho_Hext
  
  subroutine get_rho_Hext_2d(rho_Hext,s,lo,hi,ng,dx,time)
    
    use bl_constants_module
    
    integer, intent(in) :: lo(:), hi(:), ng
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1):,lo(2):)
    real(kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    
    rho_Hext = 0.0_dp_t
    
  end subroutine get_rho_Hext_2d
  
  subroutine get_rho_Hext_3d(rho_Hext,s,lo,hi,ng,dx,time)
    
    use bl_constants_module
    use geometry, only: center
    use variables, only: rho_comp

    integer, intent(in) :: lo(:), hi(:), ng
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1):,lo(2):,lo(3):)
    real(kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time

    integer :: i, j, k
    real(kind=dp_t) :: x, y, z, r, width

    width = 5.d7

    do k = lo(3), hi(3)
       z = (dble(k) + HALF)*dx(3) - center(3)

       do j = lo(2), hi(2)
          y = (dble(j) + HALF)*dx(2) - center(2)

          do i = lo(1), hi(1)
             x = (dble(i) + HALF)*dx(1) - center(1)

             r = sqrt(x**2 + y**2 + z**2)

             rho_Hext(i,j,k) = 1.d13 * s(i,j,k,rho_comp) * exp(-r**2/width**2)
             
          enddo
       enddo
    enddo

    

    
  end subroutine get_rho_Hext_3d
  
end module heating_module
