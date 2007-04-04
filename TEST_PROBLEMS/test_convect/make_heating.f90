module heating_module

  use bl_types
  use bl_constants_module
  use multifab_module

  implicit none

contains

   subroutine get_H_2d (H,lo,hi,dx,time)

      integer, intent(in) :: lo(:), hi(:)
      real(kind=dp_t), intent(inout) :: H(lo(1):,lo(2):)
      real(kind=dp_t), intent(in   ) :: dx(:),time

      integer :: i,j
      real(kind=dp_t) :: x,x0,x1,x2
      real(kind=dp_t) :: y,y0,y1,y2,y_layer
      real(kind=dp_t) :: r0,r1,r2
      real(kind=dp_t) :: ey,Hmax
      real(kind=dp_t) :: pi,L_x

      H = 0.0_dp_t
      Hmax = 0.0_dp_t

      L_x = 2.5d8
      pi = 3.1415926535897932384626433d0

      if (time <= 200.0) then

        ! First point at (0.5,.65)
        x0 = 5.0d7
        y0 = 6.5d7

        ! Second point at (1.2,..85)
        x1 = 1.2d8
        y1 = 8.5d7

        ! Third point at (2.0,.75)
        x2 = 2.d8
        y2 = 7.5d7

        y_layer = y2 + 5.0d7

        do j = lo(2),hi(2)
          y = (dble(j)+HALF)*dx(2)
          ey = exp(-(y-y_layer)*(y-y_layer)/1.e14)
          do i = lo(1),hi(1)
            x =  (dble(i)+HALF)*dx(1)

            r0 = sqrt( (x-x0)**2 +(y-y0)**2 ) / 2.5e6
            r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / 2.5e6
            r2 = sqrt( (x-x2)**2 +(y-y2)**2 ) / 2.5e6

!           H(i,j) = (ey &
!                     + .00625_dp_t * exp(-((x-x0)**2 +(y-y0)**2)/0.25e14) &
!                     + .01875_dp_t * exp(-((x-x1)**2 +(y-y1)**2)/0.25e14) &
!                     + .01250_dp_t * exp(-((x-x2)**2 +(y-y2)**2)/0.25e14) ) * 1.d17

!           H(i,j) = (  .00625_dp_t * exp(-((x-x0)**2 +(y-y0)**2)/0.25e14) &
!                     + .01875_dp_t * exp(-((x-x1)**2 +(y-y1)**2)/0.25e14) &
!                     + .01250_dp_t * exp(-((x-x2)**2 +(y-y2)**2)/0.25e14) ) * 1.d17

!            H(i,j) = ey*(.01875_dp_t * (1.d0 + sin(2*pi*x/L_x)) &
!                      +  .01250_dp_t * (1.d0 + sin((6*pi*x/L_x) + pi/3.d0)) &
!                      +  .00625_dp_t * (1.d0 + sin((8*pi*x/L_x) + pi/5.d0)) )*1.d17
 
! best so far
            H(i,j) = ey*(ONE + &
                        .00625_dp_t * sin(2*pi*x/L_x) &
                      + .01875_dp_t * sin((6*pi*x/L_x) + pi/3.d0) &
                      + .01250_dp_t * sin((8*pi*x/L_x) + pi/5.d0))*1.d17

!           H(i,j) = ey*(ONE + &
!                       .00625_dp_t * sin(2*pi*x/L_x) &
!                     + .0025_dp_t * sin(pi*x/L_x + .562) &
!                     + .01875_dp_t * sin((6*pi*x/L_x) + pi/3.d0) &
!                     + .01250_dp_t * sin((8*pi*x/L_x) + pi/5.d0))*1.d15

!           H(i,j) = ey*(ONE + &
!                       .00625_dp_t * sin(2*pi*x/L_x) &
!                     + .0025_dp_t * sin(pi*x/L_x + .562) &
!                     + .01875_dp_t * sin((6*pi*x/L_x) + pi/3.d0) &
!                     + .01250_dp_t * sin((8*pi*x/L_x) + pi/5.d0))*1.d14

!           H(i,j) = ey * 1.d17

!           ! HACK NO HEATING
!           H(i,j) = ZERO

            Hmax = max(Hmax,H(i,j))
          end do
        end do

!       if (parallel_IOProcessor()) print *,'MAX VALUE OF H ',Hmax

      end if

   end subroutine get_H_2d

   subroutine get_H_3d (H,lo,hi,dx,time)

      integer, intent(in) :: lo(:), hi(:)
      real(kind=dp_t), intent(inout) :: H(lo(1):,lo(2):,lo(3):)
      real(kind=dp_t), intent(in   ) :: dx(:),time

      integer :: i,j,k
      real(kind=dp_t) :: x,x0,x1,x2
      real(kind=dp_t) :: z,z0,z1,z2,z_layer
      real(kind=dp_t) :: r0,r1,r2
      real(kind=dp_t) :: ez,Hmax

      H = 0.0_dp_t
      Hmax = 0.0_dp_t

      if (time <= 2.0) then

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
          z = (dble(k)+HALF)*dx(3)
          ez = exp(-(z-z_layer)*(z-z_layer)/1.e14)
          do i = lo(1),hi(1)
!           y =  (dble(j)+HALF)*dx(2)
            x =  (dble(i)+HALF)*dx(1)

            r0 = sqrt( (x-x0)**2 +(z-z0)**2 ) / 2.5e6
            r1 = sqrt( (x-x1)**2 +(z-z1)**2 ) / 2.5e6
            r2 = sqrt( (x-x2)**2 +(z-z2)**2 ) / 2.5e6

!           H(i,j,k) = (ez &
!                       + .00625_dp_t * exp(-((x-x0)**2 +(z-z0)**2)/0.25e14) &
!                       + .01875_dp_t * exp(-((x-x1)**2 +(z-z1)**2)/0.25e14) &
!                       + .01250_dp_t * exp(-((x-x2)**2 +(z-z2)**2)/0.25e14) ) * 1.d17

!           H(i,j,k) = (  .00625_dp_t * exp(-((x-x0)**2 +(z-z0)**2)/0.25e14) &
!                       + .01875_dp_t * exp(-((x-x1)**2 +(z-z1)**2)/0.25e14) &
!                       + .01250_dp_t * exp(-((x-x2)**2 +(z-z2)**2)/0.25e14) ) * 1.d17

            H(i,j,k) = (  .00625_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r0))) &
                        + .01875_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r1))) &
                        + .01250_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r2))) ) * 1.d17
                      
            ! HACK NO HEATING
            H(i,j,k) = ZERO
                      
            Hmax = max(Hmax,H(i,j,k))
          end do
        end do
        end do

!       if (parallel_IOProcessor()) print *,'MAX VALUE OF H ',Hmax

      end if

   end subroutine get_H_3d
end module heating_module
