program grid_2d

  integer :: i, j

  integer, parameter :: nx = 16
  integer, parameter :: ny = 8
  integer, parameter :: nzonesx = 512
  integer, parameter :: nzonesy = 512
  integer, parameter :: iyoff = 272

  integer :: ix, iy
  integer :: nlevs, ngrids

 99 format(i1)
100 format('         ((',i4,',',i4,') (',i4,',',i4,') ('i4,',',i4,'))')
101 format('   ((',i4,',',i4,') (',i4,',',i4,') ('i4,',',i4,'))', i5)

  ix = 0
  nlevs = 1
  ngrids = nx*ny

  write (*, 99) nlevs
  write (*,101) 0,0,nzonesx-1,nzonesy-1,0,0,ngrids
  do i = 1, nx

     iy = 0
     do j = 1, ny

        write (*,100) ix,iy+iyoff,ix+nzonesx/nx-1,iy+nzonesy/ny-1+iyoff,0,0

        iy = iy + nzonesy/ny
     enddo

     ix = ix + nzonesx/nx
  enddo
        
end program grid_2d

