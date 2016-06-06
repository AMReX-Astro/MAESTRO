!This code was made available from OLCF (specifically, from Ed D'Azevedo) and
!has not been cleared for sharing.  For now, just tinkering.

        subroutine daxpy( n, a, x, incx, y, incy)
!$acc   routine seq
        implicit none
        integer n,incx,incy
        real*8 a, x(*),y(*)
        integer i,ix,iy

        if ((incx.eq.1).and.(incy.eq.1)) then
                do i=1,n
                   y(i) = a*x(i) + y(i)
                enddo
        else
                ix = 1
                iy = 1
                do i=1,n
                   ix = 1+(i-1)*incx
                   iy = 1+(i-1)*incy
                   y(iy) = a*x(ix) + y(iy)
                enddo
        endif
        return
        end subroutine daxpy
