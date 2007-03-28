module fill_3d_module

  use bl_constants_module
  use bl_types
  use multifab_module
  use variables
  use geometry

  implicit none
  
contains

  subroutine fill_3d_data (data,s0,lo,hi,dx,ng)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(  out) :: data(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(kind=dp_t), intent(in   ) ::   s0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer                  :: i,j,k,n,nr,index
    real(kind=dp_t)          :: x,y,z
    real(kind=dp_t)          :: radius

    nr = size(s0,dim=1)

    do k = lo(3),hi(3)
      z = (dble(k)+HALF)*dx(3) - center(3)
      do j = lo(2),hi(2)
        y = (dble(j)+HALF)*dx(2) - center(2)
        do i = lo(1),hi(1)
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

  subroutine make_3d_normal (normal,lo,hi,dx,ng)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(  out) :: normal(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)

    integer                  :: i,j,k
    real(kind=dp_t)          :: x,y,z,radius

    if (spherical .eq. 1) then
      do k = lo(3)-ng,hi(3)+ng
        z = (dble(k)-HALF)*dx(3) - center(3)
        do j = -lo(2)-ng,hi(2)+ng
          y = (dble(j)-HALF)*dx(2) - center(2)
          do i = -lo(1)-ng,hi(1)+ng
            x = (dble(i)-HALF)*dx(1) - center(1)
  
            radius = sqrt(x**2 + y**2 + z**2)
  
            normal(i,j,k,1) = x / radius
            normal(i,j,k,2) = y / radius
            normal(i,j,k,3) = z / radius
  
          end do
        end do
      end do
    else 
      print *,'SHOULDNT CALL MAKE_3D_NORMAL WITH SPHERICAL = 0'
      stop
    end if

  end subroutine make_3d_normal

  subroutine make_w0_cart(w0,w0_cart,normal,dx)

      real(kind=dp_t), intent(in   ) :: w0(:)
      type(multifab) , intent(inout) :: w0_cart
      type(multifab) , intent(in   ) :: normal
      real(kind=dp_t), intent(in   ) :: dx(:)

      ! Local variables
      integer :: i,lo(w0_cart%dim),hi(w0_cart%dim),dm,ng
      real(kind=dp_t), pointer :: wp(:,:,:,:)
      real(kind=dp_t), pointer :: np(:,:,:,:)

      dm = w0_cart%dim
      ng = w0_cart%ng

      call setval(w0_cart,ZERO,all=.true.)

      do i = 1, w0_cart%nboxes
         if ( multifab_remote(w0_cart, i) ) cycle
         wp => dataptr(w0_cart, i)
         lo =  lwb(get_box(w0_cart, i))
         hi =  upb(get_box(w0_cart, i))
         if (spherical .eq. 1) then
           np => dataptr(normal, i)
         else
           np => Null()
         end if
         call put_w0_on_3d_cells(w0,wp(:,:,:,:),np(:,:,:,:),lo,hi,dx,ng)
      end do

      call multifab_fill_boundary(w0_cart)

  end subroutine make_w0_cart

  subroutine put_w0_on_3d_cells (w0,w0_cell,normal,lo,hi,dx,ng)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(  out) :: w0_cell(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) ::  normal(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer                  :: i,j,k,n,nr,index
    real(kind=dp_t)          :: x,y,z
    real(kind=dp_t)          :: radius,rfac,w0_cell_val

    nr = size(w0,dim=1)

    if (spherical .eq. 1) then
     do k = lo(3),hi(3)
      z = (dble(k)+HALF)*dx(3) - center(3)
      do j = lo(2),hi(2)
        y = (dble(j)+HALF)*dx(2) - center(2)
        do i = lo(1),hi(1)
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
          rfac = (radius - dble(index)*dr) / dr
          if (rfac .lt. 0.0 .or. rfac .gt. 1.0) then
            print *,'BAD RFAC ',rfac
            print *,'RADIUS, INDEX*DR ',radius, dble(index)*dr
            x = 1.0 / 0.0
            stop
          end if
          w0_cell_val = rfac * w0(index) + (ONE-rfac) * w0(index+1)
          w0_cell(i,j,k,1) = w0_cell_val * normal(i,j,k,1)
          w0_cell(i,j,k,2) = w0_cell_val * normal(i,j,k,2)
          w0_cell(i,j,k,3) = w0_cell_val * normal(i,j,k,3)
        end do
      end do
     end do

    else 

      w0_cell = ZERO
      do k = lo(3),hi(3)
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
        w0_cell(i,j,k,3) = HALF * (w0(k) + w0(k+1))
      end do
      end do
      end do

    end if

  end subroutine put_w0_on_3d_cells

end module fill_3d_module
