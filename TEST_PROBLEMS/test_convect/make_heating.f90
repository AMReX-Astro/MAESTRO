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
    use variables, only: rho_comp
    use probin_module, only: prob_lo_x, prob_lo_y

    integer, intent(in) :: lo(:), hi(:), ng
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1):,lo(2):)
    real(kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time

    integer :: i,j
    real(kind=dp_t) :: x
    real(kind=dp_t) :: y,y_layer
    real(kind=dp_t) :: ey,Hmax
    real(kind=dp_t) :: pi,L_x

    rho_Hext = 0.0_dp_t
    Hmax = 0.0_dp_t

    L_x = 2.5d8
    pi = 3.1415926535897932384626433d0

    if (time <= 200.0) then

       y_layer = 1.25d8

       do j = lo(2),hi(2)
          y = (dble(j)+HALF)*dx(2) + prob_lo_y
          ey = exp(-(y-y_layer)*(y-y_layer)/1.e14)
          do i = lo(1),hi(1)
             x =  (dble(i)+HALF)*dx(1) + prob_lo_x

             rho_Hext(i,j) = ey*(ONE + &
                  .00625_dp_t * sin(2*pi*x/L_x) &
                  + .01875_dp_t * sin((6*pi*x/L_x) + pi/3.d0) &
                  + .01250_dp_t * sin((8*pi*x/L_x) + pi/5.d0))*2.5d16

             Hmax = max(Hmax,rho_Hext(i,j))

             rho_Hext(i,j) = rho_Hext(i,j) * s(i,j,rho_comp)
          end do
       end do

       ! if (parallel_IOProcessor()) print *,'MAX VALUE OF H ',Hmax

    end if

  end subroutine get_rho_Hext_2d

  subroutine get_rho_Hext_3d(rho_Hext,s,lo,hi,ng,dx,time)

    use bl_constants_module
    use variables, only: rho_comp
    use probin_module, only: prob_lo_x, prob_lo_y, prob_lo_z

    integer, intent(in) :: lo(:), hi(:), ng
    real(kind=dp_t), intent(inout) :: rho_Hext(lo(1):,lo(2):,lo(3):)
    real(kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),time

    integer :: i,j,k
    real(kind=dp_t) :: x, y
    real(kind=dp_t) :: z,z_layer
    real(kind=dp_t) :: ez,Hmax
    real(kind=dp_t) :: pi,L_x,L_y


    rho_Hext = 0.0_dp_t
    Hmax = 0.0_dp_t

    ! HACK -- these are the domain sizes
    L_x = 2.5d8
    L_y = 2.5d8

    pi = 3.1415926535897932384626433d0

    if (time <= 200.0) then

       z_layer = 1.25d8

       do k = lo(3),hi(3)
          z = (dble(k)+HALF)*dx(3) + prob_lo_z
          ez = exp(-(z-z_layer)*(z-z_layer)/1.e14)

          do j = lo(2),hi(2)
             y = (dble(j)+HALF)*dx(2) + prob_lo_y

             do i = lo(1),hi(1)
                x =  (dble(i)+HALF)*dx(1) + prob_lo_x

                rho_Hext(i,j,k) = ez*(ONE + &
                     .00625_dp_t * sin(2*pi*x/L_x) * sin(2*pi*y/L_y) &
                     + .01875_dp_t * sin((6*pi*x/L_x) + pi/3.d0) * sin((6*pi*y/L_y) + pi/3.d0) &
                     + .01250_dp_t * sin((8*pi*x/L_x) + pi/5.d0) * sin((8*pi*y/L_y) + pi/5.d0))*2.5d16

                Hmax = max(Hmax,rho_Hext(i,j,k))

                rho_Hext(i,j,k) = rho_Hext(i,j,k) * s(i,j,k,rho_comp)
                
             end do
          end do
       enddo

       ! if (parallel_IOProcessor()) print *,'MAX VALUE OF H ',Hmax

    end if

  end subroutine get_rho_Hext_3d

end module heating_module
