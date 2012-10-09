  program wave

  implicit none

  integer i,iter,nx,nt
  parameter (nx = 512)
  real*8 :: x(nx), p0_0(nx), p0_T_yes(nx), phse_T_no(nx), phse_T_yes(nx)
  real*8 :: dummy

  open (11,file='NO_ETA/p00000')
  do i = 1,nx
     read(11,*) x(i), p0_0(i)
  end do
  close(11)

  open (91,file='NO_ETA/phse0222')
  do i = 1,nx
     read(91,*) dummy, phse_T_no(i)
  end do
  do i = 1,nx
     phse_T_no(i) = phse_T_no(i) / p0_0(i)
  end do
  close(91)

  open (13,file='YES_ETA/phse0221')
  do i = 1,nx
     read(13,*) dummy, phse_T_yes(i)
  end do
  do i = 1,nx
     phse_T_yes(i) = phse_T_yes(i) / p0_0(i)
  end do
  close(13)

  open (14,file='YES_ETA/p00221')
  do i = 1,nx
     read(14,*) dummy, p0_T_yes(i)
  end do
  do i = 1,nx
     print *,'P0 ',p0_T_yes(i)
     p0_T_yes(i) = p0_T_yes(i) / p0_0(i)
  end do
  close(14)

  open (12,file='RATIO')
  do i = 1,nx
     write(12,1000) x(i),  p0_T_yes(i), phse_T_yes(i), phse_T_no(i)
  end do
  close(12)

  open (12,file='ONE')
  do i = 1,nx
     write(12,*) x(i),  1.d0
  end do
  close(12)

1000 format(f12.2,2x,f10.6,2x,f10.6,2x,f10.6)

  end

