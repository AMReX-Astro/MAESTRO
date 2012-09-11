! This routine returns the externally imposed (i.e. not reactions)
! heating source term to the enthalpy equation (actually rho * H
! is returned, where H has units of erg/g/s).

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
    use ml_restriction_module, only : ml_cc_restriction
    use variables, only: foextrap_comp
    use geometry, only: spherical
    use bl_constants_module, only: ZERO
    use probin_module, only: drive_initial_convection
    use fill_3d_module, only: put_1d_array_on_cart


    type(ml_layout), intent(in   ) :: mla
    real(kind=dp_t), intent(in   ) :: tempbar_init(:,0:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: rho_Hext(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt

    ! local
    integer                  :: n,i,ng_s,ng_h,ng_tc,dm,nlevs
    integer                  :: lo(mla%dim),hi(mla%dim)
    real(kind=dp_t), pointer ::  sp(:,:,:,:)
    real(kind=dp_t), pointer ::  hp(:,:,:,:)
    real(kind=dp_t), pointer :: tcp(:,:,:,:)

    type(multifab) :: tempbar_init_cart(mla%nlevel)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "get_rho_Hext")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_s = nghost(s(1))
    ng_h = nghost(rho_Hext(1))

    ! put tempbar_init on Cart
    if (spherical == 1) then
       do n=1, nlevs
          ! tempbar_init_cart will hold the initial tempbar on a Cartesian
          ! grid to be used if drive_initial_convection is true
          call build(tempbar_init_cart(n),mla%la(n),1,0)
          call setval(tempbar_init_cart(n), ZERO, all=.true.)
       enddo

       if (drive_initial_convection) then
          ! fill all components
          call put_1d_array_on_cart(tempbar_init,tempbar_init_cart, &
                                    foextrap_comp,.false.,.false.,dx, &
                                    the_bc_level,mla)
       endif
    endif



    do n=1,nlevs

       do i = 1, nfabs(s(n))
          sp => dataptr(s(n) , i)
          hp => dataptr(rho_Hext(n) , i)
          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))
          select case (dm)
          case (1)
             call get_rho_Hext_1d(tempbar_init(n,:), &
                                  hp(:,1,1,1),ng_h,sp(:,1,1,:),ng_s, &
                                  lo,hi,dx(n,:))
          case (2)
             call get_rho_Hext_2d(tempbar_init(n,:), &
                                  hp(:,:,1,1),ng_h,sp(:,:,1,:),ng_s, &
                                  lo,hi,dx(n,:))
          case (3)
             if (spherical == 1) then
                tcp => dataptr(tempbar_init_cart(n), i)
                ng_tc = nghost(tempbar_init_cart(1))
                call get_rho_Hext_3d_sph(tcp(:,:,:,1),ng_tc, &
                                         hp(:,:,:,1),ng_h,sp(:,:,:,:),ng_s, &
                                         lo,hi,dx(n,:))
             else
                call get_rho_Hext_3d(tempbar_init(n,:), &
                                     hp(:,:,:,1),ng_h,sp(:,:,:,:),ng_s, &
                                     lo,hi,dx(n,:))
             endif
          end select
       end do

    end do

    ! the loop over nlevs must count backwards to make sure the finer grids are done first
    do n=nlevs,2,-1
       ! set level n-1 data to be the average of the level n data covering it
       call ml_cc_restriction(rho_Hext(n-1), rho_Hext(n), mla%mba%rr(n-1,:))
    end do

    call destroy(bpt)

    if (spherical == 1) then
       do n = 1, nlevs
          call destroy(tempbar_init_cart(n))
       enddo
    endif

  end subroutine get_rho_Hext
  
  subroutine get_rho_Hext_1d(tempbar_init,rho_Hext,ng_h,s,ng_s,lo,hi,dx)
    
    use bl_constants_module
    
    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: tempbar_init(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    rho_Hext = 0.0_dp_t
    
  end subroutine get_rho_Hext_1d
  
  subroutine get_rho_Hext_2d(tempbar_init,rho_Hext,ng_h,s,ng_s,lo,hi,dx)
    
    use bl_constants_module
    use variables, only: rho_comp
    use time_module, only: time
    use geometry, only: center

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: tempbar_init(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    integer :: i,j

    real(kind=dp_t) :: x,y,r

    rho_Hext = 0.0_dp_t

    if (time .le. 2.d-5) then

       do j=lo(2),hi(2)
          y = (dble(j)+0.5d0)*dx(2)
          do i=lo(1),hi(1)
             x = (dble(i)+0.5d0)*dx(1)

             r = sqrt( (x-center(1))**2 +(y-center(2))**2 ) / 0.025d0

             rho_Hext(i,j) = 1.d11*(1.d0 + tanh(2.d0-r))

          end do
       end do

    end if
    
  end subroutine get_rho_Hext_2d
  
  subroutine get_rho_Hext_3d(tempbar_init,rho_Hext,ng_h,s,ng_s,lo,hi,dx)
    
    use bl_constants_module
    
    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:,lo(3)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: tempbar_init(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    rho_Hext = 0.0_dp_t
    
  end subroutine get_rho_Hext_3d

  subroutine get_rho_Hext_3d_sph(tempbar_init_cart,ng_tc, &
                                 rho_Hext,ng_h,s,ng_s, &
                                 lo,hi,dx)
    
    use bl_constants_module
    
    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h, ng_tc
    real(kind=dp_t), intent(in   ) :: tempbar_init_cart(lo(1)-ng_tc:,lo(2)-ng_tc:,lo(3)-ng_tc:)    
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:,lo(3)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    rho_Hext = 0.0_dp_t
    
  end subroutine get_rho_Hext_3d_sph
  
end module heating_module
