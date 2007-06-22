program grid

  integer :: i, j, k

  integer, parameter :: nx = 4
  integer, parameter :: ny = 4
  integer, parameter :: nz = 4
  integer, parameter :: nzonesx = 320
  integer, parameter :: nzonesy = 320
  integer, parameter :: nzonesz = 512

  integer :: ix, iy, iz
  integer :: nlevs, ngrids

 99 format(i1)
100 format('         ((',i3,',',i3,',',i3,') (',i3,',',i3,',',i3,') ('i3,',',i3,',',i3,'))')
101 format('   ((',i3,',',i3,',',i3,') (',i3,',',i3,',',i3,') ('i3,',',i3,',',i3,'))', i3)

  ix = 0
  nlevs = 1
  ngrids = nx*ny*nz

  write (*, 99) nlevs
  write (*,101) 0,0,0,nzonesx-1,nzonesy-1,nzonesz-1,0,0,0,ngrids
  do i = 1, nx

     iy = 0
     do j = 1, ny

        iz = 0
        do k = 1, nz

           write (*,100) ix,iy,iz, ix+nzonesx/nx-1,iy+nzonesy/ny-1,iz+nzonesz/nz-1,0,0,0

           iz = iz + nzonesz/nz
        enddo

        iy = iy + nzonesy/ny
     enddo

     ix = ix + nzonesx/nx
  enddo
        
end program grid

