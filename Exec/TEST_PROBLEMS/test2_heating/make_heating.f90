! This routine returns the externally imposed (i.e. not reactions)
! heating source term to the enthalpy equation (actually rho * H
! is returned, where H has units of erg/g/s).

module heating_module

  use bl_types
  use define_bc_module

  implicit none
  
  private
  public :: get_rho_Hext
  
contains

  subroutine get_rho_Hext(mla,tempbar_init,s,rho_Hext,the_bc_level,dx,dt)

    use multifab_module
    use ml_layout_module
    use ml_restriction_module
    use variables, only: foextrap_comp

    type(ml_layout), intent(in   ) :: mla
    real(kind=dp_t), intent(in   ) :: tempbar_init(:,0:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt

    ! local
    integer                  :: n,i,ng_s,ng_h,dm, nlevs
    integer                  :: lo(mla%dim),hi(mla%dim)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: hp(:,:,:,:)
    type(bl_prof_timer), save :: bpt

    call build(bpt, "get_rho_Hext")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_s = s(1)%ng
    ng_h = rho_Hext(1)%ng

    do n=1,nlevs

       do i = 1, nfabs(s(n))
          sp => dataptr(s(n) , i)
          hp => dataptr(rho_Hext(n) , i)
          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))
          select case (dm)
          case (2)
             call get_rho_Hext_2d(hp(:,:,1,1),ng_h,sp(:,:,1,:),ng_s,lo,hi,dx(n,:))
          case (3)
             call get_rho_Hext_3d(hp(:,:,:,1),ng_h,sp(:,:,:,:),ng_s,lo,hi,dx(n,:))
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
  
  subroutine get_rho_Hext_2d(rho_Hext,ng_h,s,ng_s,lo,hi,dx)
    
    use bl_constants_module
    use variables, only: rho_comp
    
    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    integer         :: i,j
    real(kind=dp_t) :: y,width

    width = 2.5d6

    do j=lo(2),hi(2)

       y = (dble(j) + HALF)*dx(2)

       do i=lo(1),hi(1)

          rho_Hext(i,j) = 1.d17*s(i,j,rho_comp) * exp(-(y-6.5d7)**2/width**2)

       end do

    end do
   
  end subroutine get_rho_Hext_2d
  
  subroutine get_rho_Hext_3d(rho_Hext,ng_h,s,ng_s,lo,hi,dx)
    
    use bl_constants_module
    use geometry, only: center
    use variables, only: rho_comp

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:,lo(3)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    rho_Hext = 0.d0

  end subroutine get_rho_Hext_3d
  
end module heating_module
