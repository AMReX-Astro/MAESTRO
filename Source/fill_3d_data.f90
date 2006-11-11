module fill_3d_module

  use bl_constants_module
  use bl_types
  use multifab_module
  use variables
  use geometry

  implicit none
  
contains

  subroutine fill_3d_data (data,s0,dx,ng)

    integer        , intent(in   ) :: ng
    real(kind=dp_t), intent(  out) :: data(-ng:,-ng:,-ng:)
    real(kind=dp_t), intent(in   ) ::   s0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer                  :: i,j,k,n,index
    integer                  :: nr,nx,ny,nz
    real(kind=dp_t)          :: x,y,z
    real(kind=dp_t)          :: radius

    nr = size(s0,dim=1)

    nx = size(data,dim=1)-2*ng
    ny = size(data,dim=2)-2*ng
    nz = size(data,dim=3)-2*ng

    do k = 0,nz-1
      z = (dble(k)+HALF)*dx(3) - center(3)
      do j = 0,nz-1
        y = (dble(j)+HALF)*dx(2) - center(2)
        do i = 0,nz-1
          x = (dble(i)+HALF)*dx(1) - center(1)
          radius = sqrt(x**2 + y**2 + z**2)
          index = radius / dr
          if (index .lt. 0 .or. index .gt. nr-1) then
            print *,'RADIUS ',radius
            print *,'BOGUS INDEX IN FILL_3D: ',index
            print *,'NOT IN RANGE 0 TO ',nr-1
            print *,'I J K ',i,j,k
            print *,'X Y Z ',x,y,z
            x = 1.0 / 0.0
            stop
          end if
          data(i,j,k) = s0(index)
        end do
      end do
    end do

  end subroutine fill_3d_data

  subroutine make_3d_normal (dx,normal,ng)

    integer        , intent(in   ) :: ng
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(  out) :: normal(-ng:,-ng:,-ng:,:)

    integer                  :: i,j,k,nx,ny,nz
    real(kind=dp_t)          :: x,y,z,radius

    nx = size(normal,dim=1)-2*ng
    ny = size(normal,dim=2)-2*ng
    nz = size(normal,dim=3)-2*ng

    do k = -ng,nz-1+ng
      z = (dble(k)-HALF)*dx(3) - center(3)
      do j = -ng,nz-1+ng
        y = (dble(j)-HALF)*dx(2) - center(2)
        do i = -ng,nz-1+ng
          x = (dble(i)-HALF)*dx(1) - center(1)

          radius = sqrt(x**2 + y**2 + z**2)

          normal(i,j,k,1) = x / radius
          normal(i,j,k,2) = y / radius
          normal(i,j,k,3) = z / radius

        end do
      end do
    end do

  end subroutine make_3d_normal

end module fill_3d_module
