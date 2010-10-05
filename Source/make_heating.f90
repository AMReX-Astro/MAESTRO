! This routine returns the externally imposed (i.e. not reactions)
! heating source term to the enthalpy equation (actually rho * H
! is returned, where H has units of erg/g/s).

module heating_module

  use bl_types

  implicit none
  
  private
  public :: get_rho_Hext
  
contains

  subroutine get_rho_Hext(mla,s,rho_Hext,dx,time,dt)

    use multifab_module
    use ml_layout_module
    use ml_restriction_module, only : ml_cc_restriction
    use variables, only: foextrap_comp

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),time,dt

    ! local
    integer                  :: n,i,ng_s,ng_h,dm,nlevs
    integer                  :: lo(mla%dim),hi(mla%dim)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: hp(:,:,:,:)
    type(bl_prof_timer), save :: bpt

    call build(bpt, "get_rho_Hext")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_s = nghost(s(1))
    ng_h = nghost(rho_Hext(1))

    do n=1,nlevs

       do i = 1, nboxes(s(n))
          if ( multifab_remote(s(n), i) ) cycle
          sp => dataptr(s(n) , i)
          hp => dataptr(rho_Hext(n) , i)
          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))
          select case (dm)
          case (1)
             call get_rho_Hext_1d(hp(:,1,1,1),ng_h,sp(:,1,1,:),ng_s,lo,hi,dx(n,:),time)
          case (2)
             call get_rho_Hext_2d(hp(:,:,1,1),ng_h,sp(:,:,1,:),ng_s,lo,hi,dx(n,:),time)
          case (3)
             call get_rho_Hext_3d(hp(:,:,:,1),ng_h,sp(:,:,:,:),ng_s,lo,hi,dx(n,:),time)
          end select
       end do

    end do

    ! the loop over nlevs must count backwards to make sure the finer grids are done first
    do n=nlevs,2,-1
       ! set level n-1 data to be the average of the level n data covering it
       call ml_cc_restriction(rho_Hext(n-1), rho_Hext(n), mla%mba%rr(n-1,:))
    end do

    call destroy(bpt)

  end subroutine get_rho_Hext
  
  subroutine get_rho_Hext_1d(rho_Hext,ng_h,s,ng_s,lo,hi,dx,time)
    
    use bl_constants_module
    
    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    
    rho_Hext = 0.0_dp_t
    
  end subroutine get_rho_Hext_1d
  
  subroutine get_rho_Hext_2d(rho_Hext,ng_h,s,ng_s,lo,hi,dx,time)
    
    use bl_constants_module
    
    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    
    rho_Hext = 0.0_dp_t
    
  end subroutine get_rho_Hext_2d
  
  subroutine get_rho_Hext_3d(rho_Hext,ng_h,s,ng_s,lo,hi,dx,time)
    
    use bl_constants_module
    
    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:,lo(3)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    
    rho_Hext = 0.0_dp_t
    
  end subroutine get_rho_Hext_3d
  
end module heating_module
