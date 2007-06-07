program grid

  integer :: i, j, k

! integer, parameter :: nx = 4
! integer, parameter :: ny = 4
! integer, parameter :: nz = 4
! integer, parameter :: nzones = 256

  integer, parameter :: nx = 2
  integer, parameter :: ny = 2
  integer, parameter :: nz = 2
  integer, parameter :: nzones = 128

  integer :: ix, iy, iz
  integer :: nlevs, ngrids

 99 format(i1)
100 format('         ((',i3,',',i3,',',i3,') (',i3,',',i3,',',i3,') ('i3,',',i3,',',i3,'))')
101 format('   ((',i3,',',i3,',',i3,') (',i3,',',i3,',',i3,') ('i3,',',i3,',',i3,'))', i3)

  ix = 0
  nlevs = 1
  ngrids = nx*ny*nz

  write (*, 99) nlevs
  write (*,101) 0,0,0,nzones-1,nzones-1,nzones-1,0,0,0,ngrids
  do i = 1, nx

     iy = 0
     do j = 1, ny

        iz = 0
        do k = 1, nz

           write (*,100) ix,iy,iz, ix+nzones/nx-1,iy+nzones/ny-1,iz+nzones/nz-1,0,0,0

           iz = iz + nzones/nz
        enddo

        iy = iy + nzones/ny
     enddo

     ix = ix + nzones/nx
  enddo
        
end program grid

