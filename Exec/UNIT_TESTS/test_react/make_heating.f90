! This routine returns the externally imposed (i.e. not reactions)
! heating source term to the enthalpy equation (actually rho * H
! is returned, where H has units of erg/g/s).

! This is a customized heating module for use in the reactions unit test.
! A simple Gaussian heat profile is generated and used in the testing of
! the react_state() subroutine.

module heating_module

  use bl_types

  implicit none
  
  private
  public :: get_rho_Hext
  
contains

  subroutine get_rho_Hext(mla,tempbar_init,s,rho_Hext,the_bc_level,dx,dt)

    use multifab_module
    use define_bc_module
    use ml_layout_module
    use ml_cc_restriction_module, only : ml_cc_restriction
    use variables, only: foextrap_comp

    type(ml_layout), intent(in   ) :: mla
    real(kind=dp_t), intent(in   ) :: tempbar_init(:,0:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt

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

       do i = 1, nfabs(s(n))
          sp => dataptr(s(n) , i)
          hp => dataptr(rho_Hext(n) , i)
          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))
          select case (dm)
          case (2)
             call bl_error('ERROR: 2D not implemented in test_react.')
          case (3)
             call get_rho_Hext_3d(hp(:,:,:,1),ng_h,sp(:,:,:,:),ng_s, &
                                  lo,hi,dx(n,:))
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

  subroutine get_rho_Hext_3d(rho_Hext,ng_h,s,ng_s,lo,hi,dx)
    
    use bl_constants_module 
    use bl_error_module
    use variables,     only: rho_comp
    use probin_module, only: spherical_in

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:,lo(3)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
   
    integer         :: ii, jj, kk
    real(kind=dp_t) :: x, y, z, r, r0, sig, amp
 
    if(spherical_in /= 0) then
      call bl_error('ERROR: Spherical geometry not implemented in test_react.')
    endif
     
    amp = 1.d17 
    r0  = ((dble(hi(3)) + HALF)*dx(3) - (dble(lo(3)) + HALF)*dx(3))/2.0_dp_t
    sig = r0*HALF

    do kk = lo(3), hi(3)
      z = (dble(kk) + HALF)*dx(3)
      do jj = lo(2), hi(2)
        y = (dble(jj) + HALF)*dx(2)
        do ii = lo(1), hi(1)
          x = (dble(ii) + HALF)*dx(1)

          rho_Hext(ii,jj,kk) = s(ii,jj,kk,rho_comp) * amp*exp(-(x-r0)**2/sig**2)*exp(-(y-r0)**2/sig**2)*exp(-(z-r0)**2/sig**2)
        enddo
      enddo
    enddo

  end subroutine get_rho_Hext_3d
  
end module heating_module
