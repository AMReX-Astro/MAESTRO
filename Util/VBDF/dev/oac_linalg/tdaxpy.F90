        program tdaxpy
        use blas_mod
        implicit none

        integer n,m
        parameter(n=1024*16,m=64)
        integer t1,t2,count_rate


        integer it,ntimes
        integer i,j,ij,incx,incy
        real*8 x(n,m),y(n,m), a

        a = dble(1)
        do j=1,m
        do i=1,n
           ij = 1 + (j-1)*n
           x(i,j) = dble(ij)/dble(n)
           y(i,j) = dble(0)
        enddo
        enddo

        incx = 1
        incy = 1

        ntimes = 10
        call system_clock(t1,count_rate)

!$acc   data pcopyin(a,incx,incy,x) pcopy(y)
        do it=1,ntimes
!$acc   kernels
!$acc   loop independent gang
        do j=1,m
           call daxpy( n, a, x(1,j), incx, y(1,j), incy )
        enddo
!$acc   end kernels
        enddo
!$acc   end data

        call system_clock(t2,count_rate)

        write(*,*) 'n,m ', n,m
        write(*,*) 'time for daxpy ', real(t2-t1)/real(count_rate)
        write(*,*) 'sum(y) ', sum(y)


        stop
        end

