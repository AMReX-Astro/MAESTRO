! This routine returns the externally imposed (i.e. not reactions)
! heating source term to the enthalpy equation (actually rho * H
! is returned, where H has units of erg/g/s).

module heating_module

  use bl_types

  implicit none
  
  private
  public :: get_rho_Hext
  
contains

  subroutine get_rho_Hext(mla,s,rho_Hext,dx,time)

    use geometry, only: nlevs
    use multifab_module
    use ml_layout_module
    use ml_restriction_module
    use geometry, only: dm

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),time

    ! local
    integer                  :: n,i,ng_s,ng_h
    integer                  :: lo(dm),hi(dm)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: hp(:,:,:,:)

    ng_s = s(1)%ng
    ng_h = rho_Hext(1)%ng

    do n=1,nlevs

       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          sp => dataptr(s(n) , i)
          hp => dataptr(rho_Hext(n) , i)
          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))
          select case (dm)
          case (2)
             call get_rho_Hext_2d(hp(:,:,1,1), ng_h, sp(:,:,1,:), ng_s, lo, hi, dx(n,:), time)
          case (3)
             call get_rho_Hext_3d(hp(:,:,:,1), ng_h, sp(:,:,:,:), ng_s, lo, hi, dx(n,:), time)
          end select
       end do

    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(rho_Hext(nlevs),1,1)

       ! fill non-periodic domain boundary ghost cells
!       call multifab_physbc(rho_Hext(nlevs),1,foextrap_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1
          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(rho_Hext(n-1), rho_Hext(n), mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
!          call multifab_fill_ghost_cells(rho_Hext(n),rho_Hext(n-1), &
!                                         ng_h,mla%mba%rr(n-1,:), &
!                                         the_bc_level(n-1),the_bc_level(n), &
!                                         1,foextrap_comp,1)      
       end do

    end if

  end subroutine get_rho_Hext
  
  subroutine get_rho_Hext_2d(rho_Hext,ng_h,s,ng_s,lo,hi,dx,time)
    
    use bl_constants_module
    use variables, only: rho_comp
    
    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    
    integer         :: i,j
    real(kind=dp_t) :: y,width

    width = 5.d7

    do j=lo(2),hi(2)

       y = (dble(j) + HALF)*dx(2)

       do i=lo(1),hi(1)

          rho_Hext(i,j) = 1.d16*s(i,j,rho_comp) * exp(-y**2/width**2)

       end do

    end do
   
  end subroutine get_rho_Hext_2d
  
  subroutine get_rho_Hext_3d(rho_Hext,ng_h,s,ng_s,lo,hi,dx,time)
    
    use bl_constants_module
    use geometry, only: center
    use variables, only: rho_comp

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:,lo(3)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    
    rho_Hext = 0.d0

  end subroutine get_rho_Hext_3d
  
end module heating_module
