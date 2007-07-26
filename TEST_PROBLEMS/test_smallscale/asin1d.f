c-------------------------------------------------
c     
c     
c     -USAGE-            
c     
c     OVERVIEW: call asin1d (fname,xlo,xhi,y,Niny)   
c     computes y(xlo,xhi) by averaging data over  
c     the range [xlo,xhi]. Here x is a scalar, y  
c     a vector, with Niny components used. piece- 
c     constant extrapolation is applied to points 
c     laying outside the range of the orig data   
c     reads data from file fname
c     
c     CALLING PROGRAM NEEDS: none                 
c     PROCEDURE NEEDS:  none                      
c     COMPILER NEEDS: f90
c     SYSTEM NEEDS: none                          
c     
      subroutine asin1d(fname, xlo, xhi, y_vector, Niny, reverse)
      implicit none
c     
      character(len=*) fname
      integer Niny            
      double precision xlo, xhi
      double precision y_vector(Niny)
c     
      integer i, j, k
      integer lo_loside, lo_hiside       
      integer hi_loside, hi_hiside
      double precision ssum    
      double precision ylo, yhi, x1, y1, x2, y2, dydx   
      double precision, allocatable :: x_data(:)
      double precision, allocatable :: y_data(:,:)
      integer uni, ios
      integer N, M
      logical first
      logical reverse
      save N, M, first, y_data, x_data
      data first /.TRUE./

      if ( first ) then
         uni = 88
         first = .FALSE.
         open(unit = uni, file = fname, status = 'old', err = 7001)
         go to 7002
 7001    print *,'ASIN1D: Open failed on file ', fname
!         call bl_error("STOP")
 7002    read(unit = uni, fmt = *, iostat = ios) N, M
         if ( ios .ne. 0) then
!            call bl_error('ASIN1D:0: failed to read N, M')
         end if
         allocate(x_data(N))
         allocate(y_data(N,M))
         read(unit = uni, fmt = *, iostat = ios) (x_data(i),i=1,N)
         if ( ios .ne. 0 ) then
!            call bl_error('ASIN1D:1: failed to read x_data')
         end if
         if (reverse) then
            do i = N,1,-1
               read(unit = uni, fmt = *, iostat = ios) 
     &              (y_data(i,j),j=1,M)
               if ( ios .ne. 0 ) then
!                  call bl_error('ASIN1D:2: failed to read y_data')
               end if
            end do
         else
            do i = 1, N
               read(unit = uni, fmt = *, iostat = ios) 
     &              (y_data(i,j),j=1,M)
               if ( ios .ne. 0 ) then
!                  call bl_error('ASIN1D:2: failed to read y_data')
               end if
            end do
         endif

         close(unit = uni)
      end if
      if ( Niny .lt. M ) then
         print *, 'ASIN1D: Niny(', Niny, ') to small for N(', M, ')'
!         call bl_error("STOP")
      end if

      lo_loside = 0
      lo_hiside = 0
      hi_loside = 0
      hi_hiside = 0

      if (xlo .le. x_data(1)) then
         lo_loside = 1
         lo_hiside = 1
      end if
      if (xhi .le. x_data(1)) then
         hi_loside = 1
         hi_hiside = 1
      end if
      if (xlo .ge. x_data(N)) then
         lo_loside = N
         lo_hiside = N
      end if
      if (xhi .ge. x_data(N)) then
         hi_loside = N
         hi_hiside = N
      end if

      if (lo_loside.eq.0) then
         do i = 1, N-1                           
            if ( (xlo .ge. x_data(i)) .and.
     &           (xlo .lt. x_data(i+1)) ) then
               lo_loside  = i
               lo_hiside  = i+1
               exit
            end if
         end do
      end if

      if (hi_loside.eq.0) then            
         do i = 1, N-1                           
            if ( (xhi .ge. x_data(i)) .and.
     &           (xhi .lt. x_data(i+1)) ) then
               hi_loside = i
               hi_hiside = i + 1
               exit
            end if
         end do
      end if
      
      do j = 1, M
         x1 = x_data(lo_loside)
         y1 = y_data(lo_loside,j)
         x2 = x_data(lo_hiside)
         y2 = y_data(lo_hiside,j)
         if (lo_loside.eq.lo_hiside) then
            dydx = 0.d0
         else
            dydx = (y2-y1)/(x2-x1)
         end if
         ylo = y1 + dydx*(xlo - x1)
         if ( lo_loside .eq. hi_loside ) then
            yhi = y1 + dydx*(xhi - x1)
            y_vector(j) = 0.5d0*(ylo + yhi)
         else
            ssum = (x2 - xlo) * 0.5d0 * (ylo + y2)
            x1 = x_data(hi_loside)
            y1 = y_data(hi_loside,j)
            x2 = x_data(hi_hiside)
            y2 = y_data(hi_hiside,j)
            if ( hi_loside .eq. hi_hiside ) then
               dydx = 0.d0
            else
               dydx = (y2-y1)/(x2-x1)
            end if
            yhi = y1 + dydx*(xhi - x1)
            ssum = ssum + (xhi - x1)*0.5d0*(yhi+y1)
            do k = lo_hiside, hi_loside - 1
               ssum = ssum + (x_data(k+1) - x_data(k))
     &              * 0.5d0
     &              * (y_data(k,j) + y_data(k+1,j))
            end do
            y_vector(j) = ssum / (xhi - xlo)
         end if
      end do

      end                                         
