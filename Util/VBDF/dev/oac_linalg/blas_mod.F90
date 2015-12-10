        module blas_mod
        implicit none
        contains

        subroutine daxpy( n, a, x, incx, y, incy)
!$acc   routine seq
        implicit none
        integer n,incx,incy
        real*8 a, x(*),y(*)
        integer i,ix,iy

!!$acc   data pcopyin(x) pcopy(y)
        if ((incx.eq.1).and.(incy.eq.1)) then
!!$acc           kernels
!!$acc           loop independent vector  private(i)
                do i=1,n
                   y(i) = a*x(i) + y(i)
                enddo
!!$acc           end kernels
        else
                ix = 1
                iy = 1
!!$acc           kernels
!!$acc           loop independent vector private(i,ix,iy)
                do i=1,n
                   ix = 1+(i-1)*incx
                   iy = 1+(i-1)*incy
                   y(iy) = a*x(ix) + y(iy)
                enddo
!!$acc           end kernels
        endif
!!$acc   end data
        return
        end subroutine daxpy

        end module blas_mod
