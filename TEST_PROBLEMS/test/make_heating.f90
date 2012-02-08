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
    use ml_restriction_module
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
    use probin_module, only: prob_lo
    use time_module, only: time
    
    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer :: i,j
    real(kind=dp_t) :: x,x0,x1,x2
    real(kind=dp_t) :: y,y0,y1,y2,y_layer
    real(kind=dp_t) :: r0,r1,r2
    real(kind=dp_t) :: ey,Hmax

    Hmax = 0.0_dp_t

    if (time <= 2.0_dp_t) then

       ! First point at (0.5,.65)
       x0 = 5.0d7
       y0 = 6.5d7

       ! Second point at (1.2,..85)
       x1 = 1.2d8
       y1 = 8.5d7

       ! Third point at (2.0,.75)
       x2 = 2.d8
       y2 = 7.5d7

       y_layer = y2

       do j = lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+HALF)*dx(2)
          ey = exp(-(y-y_layer)*(y-y_layer)/1.e14)
          do i = lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+HALF)*dx(1)

             r0 = sqrt( (x-x0)**2 +(y-y0)**2 ) / 2.5e6
             r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / 2.5e6
             r2 = sqrt( (x-x2)**2 +(y-y2)**2 ) / 2.5e6

             ! rho_Hext(i,j) = (ey &
             !      + .00625_dp_t * exp(-((x-x0)**2 +(y-y0)**2)/0.25e14) &
             !      + .01875_dp_t * exp(-((x-x1)**2 +(y-y1)**2)/0.25e14) &
             !      + .01250_dp_t * exp(-((x-x2)**2 +(y-y2)**2)/0.25e14) ) * 1.d17

             ! rho_Hext(i,j) = (  .00625_dp_t * exp(-((x-x0)**2 +(y-y0)**2)/0.25e14) &
             !      + .01875_dp_t * exp(-((x-x1)**2 +(y-y1)**2)/0.25e14) &
             !      + .01250_dp_t * exp(-((x-x2)**2 +(y-y2)**2)/0.25e14) ) * 1.d17

             rho_Hext(i,j) = ( ey + &
                  .00625_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0_dp_t-r0))) &
                  + .01875_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0_dp_t-r1))) &
                  + .01250_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0_dp_t-r2))) ) * 1.d17

             ! rho_Hext(i,j) = ey * 1.d17

             ! HACK NO HEATING
             ! rho_Hext(i,j) = ZERO

             Hmax = max(Hmax,rho_Hext(i,j))

             rho_Hext(i,j) = rho_Hext(i,j) * s(i,j,rho_comp)
          end do
       end do

    end if

  end subroutine get_rho_Hext_2d
    
  subroutine get_rho_Hext_3d(rho_Hext,ng_h,s,ng_s,lo,hi,dx)

    use bl_constants_module
    use variables, only: rho_comp
    use probin_module, only: prob_lo
    use time_module, only: time

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_h
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1)-ng_h:,lo(2)-ng_h:,lo(3)-ng_h:)
    real(kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer :: i,j,k
    real(kind=dp_t) :: x,x0,x1,x2
    real(kind=dp_t) :: z,z0,z1,z2,z_layer
    real(kind=dp_t) :: r0,r1,r2
    real(kind=dp_t) :: ez,Hmax

    Hmax = 0.0_dp_t

    if (time <= 2.0_dp_t) then

       ! First point at (0.5,.65)
       x0 = 5.0d7
       z0 = 6.5d7

       ! Second point at (1.2,..85)
       x1 = 1.2d8
       z1 = 8.5d7

       ! Third point at (2.0,.75)
       x2 = 2.d8
       z2 = 7.5d7

       z_layer = z2

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             z = prob_lo(3) + (dble(k)+HALF)*dx(3)
             ez = exp(-(z-z_layer)*(z-z_layer)/1.e14)
             do i = lo(1),hi(1)
                ! y = prob_lo(2) + (dble(j)+HALF)*dx(2)
                x = prob_lo(1) + (dble(i)+HALF)*dx(1)

                r0 = sqrt( (x-x0)**2 +(z-z0)**2 ) / 2.5e6
                r1 = sqrt( (x-x1)**2 +(z-z1)**2 ) / 2.5e6
                r2 = sqrt( (x-x2)**2 +(z-z2)**2 ) / 2.5e6

                ! rho_Hext(i,j,k) = (ez &
                !      + .00625_dp_t * exp(-((x-x0)**2 +(z-z0)**2)/0.25e14) &
                !      + .01875_dp_t * exp(-((x-x1)**2 +(z-z1)**2)/0.25e14) &
                !      + .01250_dp_t * exp(-((x-x2)**2 +(z-z2)**2)/0.25e14) ) * 1.d17

                ! rho_Hext(i,j,k) = (  .00625_dp_t * exp(-((x-x0)**2 +(z-z0)**2)/0.25e14) &
                !      + .01875_dp_t * exp(-((x-x1)**2 +(z-z1)**2)/0.25e14) &
                !      + .01250_dp_t * exp(-((x-x2)**2 +(z-z2)**2)/0.25e14) ) * 1.d17

                rho_Hext(i,j,k) = (  .00625_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0_dp_t-r0))) &
                     + .01875_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0_dp_t-r1))) &
                     + .01250_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0_dp_t-r2))) ) * 1.d17

                ! HACK NO HEATING
                rho_Hext(i,j,k) = ZERO

                Hmax = max(Hmax,rho_Hext(i,j,k))

                rho_Hext(i,j,k) = rho_Hext(i,j,k) * s(i,j,k,rho_comp)
             end do
          end do
       end do

    end if

  end subroutine get_rho_Hext_3d

end module heating_module
