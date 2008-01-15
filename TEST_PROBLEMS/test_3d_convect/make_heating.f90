! This routine returns the externally imposed (i.e. not reactions)
! heating source term to the enthalpy equation (actually rho * H
! is returned, where H has units of erg/g/s).

module heating_module

  use bl_types
  use bl_constants_module
  use multifab_module

  implicit none

contains

  subroutine get_H_2d (H,lo,hi,dx,time)

    use probin_module, only: prob_lo_x, prob_lo_y

    integer, intent(in) :: lo(:), hi(:)
    real(kind=dp_t), intent(inout) :: H(lo(1):,lo(2):)
    real(kind=dp_t), intent(in   ) :: dx(:),time

    integer :: i,j
    real(kind=dp_t) :: x,y,y_layer
    real(kind=dp_t) :: ey,Hmax
    real(kind=dp_t) :: pi,L_x

    H = 0.0_dp_t
    Hmax = 0.0_dp_t

    ! HACK -- this is the domain size
    L_x = 2.5d8
    pi = 3.1415926535897932384626433d0

    if (time <= 200.0) then

       y_layer = 1.25d7

       do j = lo(2),hi(2)
          y = prob_lo_y + (dble(j)+HALF)*dx(2)
          ey = exp(-(y-y_layer)*(y-y_layer)/1.e14)

          do i = lo(1),hi(1)
             x = prob_lo_x + (dble(i)+HALF)*dx(1)

             ! best so far
             H(i,j) = ey*(ONE + &
                  .00625_dp_t * sin(2*pi*x/L_x) &
                  + .01875_dp_t * sin((6*pi*x/L_x) + pi/3.d0) &
                  + .01250_dp_t * sin((8*pi*x/L_x) + pi/5.d0))*2.5d16


             Hmax = max(Hmax,H(i,j))
          end do
       end do

       !       if (parallel_IOProcessor()) print *,'MAX VALUE OF H ',Hmax

    end if

  end subroutine get_H_2d

  subroutine get_H_3d (H,lo,hi,dx,time)

    use probin_module, only: prob_lo_x, prob_lo_y, prob_lo_z

    integer, intent(in) :: lo(:), hi(:)
    real(kind=dp_t), intent(inout) :: H(lo(1):,lo(2):,lo(3):)
    real(kind=dp_t), intent(in   ) :: dx(:),time

    integer :: i,j,k

    real(kind=dp_t) :: x,y,z,z_layer
    real(kind=dp_t) :: ez,Hmax
    real(kind=dp_t) :: L_x, L_y, pi

    H = 0.0_dp_t
    Hmax = 0.0_dp_t

    ! HACK -- these are the domain sizes
    L_x = 2.5d8
    L_y = 2.5d8

    pi = 3.1415926535897932384626433d0

    if (time <= 200.0) then

       z_layer = 1.25d7

       do k = lo(3),hi(3)
          z = prob_lo_z + (dble(k)+HALF)*dx(3)
          ez = exp(-(z-z_layer)*(z-z_layer)/1.e14)

          do j = lo(2),hi(2)
             y = prob_lo_y + (dble(j)+HALF)*dx(2)

             do i = lo(1),hi(1)
                x = prob_lo_x + (dble(i)+HALF)*dx(1)

                H(i,j,k) = ez*(ONE + &
                     .00625_dp_t * sin(2*pi*x/L_x) * sin(2*pi*y/L_y) &
                     + .01875_dp_t * sin((6*pi*x/L_x) + pi/3.d0) * sin((6*pi*y/L_y) + pi/3.d0) &
                     + .01250_dp_t * sin((8*pi*x/L_x) + pi/5.d0) * sin((8*pi*y/L_y) + pi/5.d0))*2.5d16


                Hmax = max(Hmax,H(i,j,k))
             end do
          end do
       end do

       !       if (parallel_IOProcessor()) print *,'MAX VALUE OF H ',Hmax

    end if

  end subroutine get_H_3d
end module heating_module
