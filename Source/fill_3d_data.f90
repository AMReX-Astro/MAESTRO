module fill_3d_module

  use bl_types
  use multifab_module

  implicit none

  private

  public :: fill_3d_data_wrapper, fill_3d_data
  public :: make_3d_normal, make_w0_cart, put_w0_on_3d_cells_sphr
  
contains

  subroutine fill_3d_data_wrapper(nlevs,s0_cart,s0,dx,in_comp)

    use bl_constants_module
    use geometry

    ! for spherical problems, this copies the base state onto a multifab
    ! sames as the function fill_3d_data_wrap, excpet we assume
    ! start_comp = 1, num_comp = 1, and the base state only has one component

    integer        , intent(in   )        :: nlevs
    type(multifab) , intent(inout)        :: s0_cart(:)
    real(kind=dp_t), intent(in   )        :: s0(:,0:)
    real(kind=dp_t), intent(in   )        :: dx(:,:)
    integer        , intent(in), optional :: in_comp

    ! local
    integer :: i,comp,ng,n
    integer :: lo(s0_cart(1)%dim),hi(s0_cart(1)%dim)
    real(kind=dp_t), pointer :: s0p(:,:,:,:)

    ng = s0_cart(1)%ng

    if ( present(in_comp) ) then
       comp = in_comp
    else
       comp = 1
    endif

    do n=1,nlevs
       do i=1,s0_cart(n)%nboxes
          if ( multifab_remote(s0_cart(n),i) ) cycle
          s0p => dataptr(s0_cart(n),i)
          lo = lwb(get_box(s0_cart(n),i))
          hi = upb(get_box(s0_cart(n),i))
          call fill_3d_data(n,s0p(:,:,:,comp),s0(n,:),lo,hi,dx(n,:),ng)
       end do
       
       call multifab_fill_boundary_c(s0_cart(n),comp,1)
    enddo

  end subroutine fill_3d_data_wrapper

  subroutine fill_3d_data(n,data,s0,lo,hi,dx,ng)

    use bl_constants_module
    use geometry, z_geometry => z ! Because there's a local z variable.
    
    integer        , intent(in   ) :: n,lo(:),hi(:),ng
    real(kind=dp_t), intent(  out) :: data(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(kind=dp_t), intent(in   ) ::   s0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    integer                  :: i,j,k,index
    real(kind=dp_t)          :: x,y,z
    real(kind=dp_t)          :: radius
    
    do k = lo(3),hi(3)
       z = (dble(k)+HALF)*dx(3) - center(3)
       do j = lo(2),hi(2)
          y = (dble(j)+HALF)*dx(2) - center(2)
          do i = lo(1),hi(1)
             x = (dble(i)+HALF)*dx(1) - center(1)
             radius = sqrt(x**2 + y**2 + z**2)
             index = int(radius / dr(n))
             if (index .lt. 0 .or. index .gt. nr(n)-1) then
                print *,'RADIUS ',radius
                print *,'BOGUS INDEX IN FILL_3D: ',index
                print *,'NOT IN RANGE 0 TO ',nr(n)-1
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

    use bl_constants_module
    use geometry, z_geometry => z ! Because there's a local z variable.
    
    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(  out) :: normal(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)

    integer                  :: i,j,k
    real(kind=dp_t)          :: x,y,z,radius

    if (spherical .eq. 1) then
      do k = lo(3)-ng,hi(3)+ng
        z = (dble(k)+HALF)*dx(3) - center(3)
        do j = lo(2)-ng,hi(2)+ng
          y = (dble(j)+HALF)*dx(2) - center(2)
          do i = lo(1)-ng,hi(1)+ng
            x = (dble(i)+HALF)*dx(1) - center(1)
  
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

  subroutine make_w0_cart(nlevs,w0,w0_cart,normal,dx)

    use bl_constants_module
    use geometry
    
    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(inout) :: w0_cart(:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    
    ! Local variables
    integer :: lo(w0_cart(1)%dim),hi(w0_cart(1)%dim)
    integer :: i,n,dm,ng
    real(kind=dp_t), pointer :: wp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)
    
    dm = w0_cart(1)%dim
    ng = w0_cart(1)%ng
    
    do n = 1, nlevs
    
       call setval(w0_cart(n),ZERO,all=.true.)
       
       do i = 1, w0_cart(n)%nboxes
          if ( multifab_remote(w0_cart(n), i) ) cycle
          wp => dataptr(w0_cart(n), i)
          lo = lwb(get_box(w0_cart(n), i))
          hi = upb(get_box(w0_cart(n), i))
          if (spherical .eq. 1) then
             np => dataptr(normal(n), i)
             call put_w0_on_3d_cells_sphr(n,w0(n,:),wp(:,:,:,:),np(:,:,:,:),lo,hi,dx(n,:),ng)
          else
             call put_w0_on_3d_cells_cart(n,w0(n,:),wp(:,:,:,:),lo,hi,dx(n,dm),ng)
          end if
       end do
       
       call multifab_fill_boundary(w0_cart(n))

    enddo
    
  end subroutine make_w0_cart
  
  subroutine put_w0_on_3d_cells_cart(n,w0,w0_cell,lo,hi,dz,ng)

    use bl_constants_module
    use geometry

    integer        , intent(in   ) :: n,lo(:),hi(:),ng
    real(kind=dp_t), intent(  out) :: w0_cell(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: dz

    integer         :: i,j,k
    integer         :: rr,klo,khi

    rr = int( dz / dr(n) + 1.d-12)

    w0_cell = ZERO
    do k = lo(3),hi(3)
    do j = lo(2),hi(2)
    do i = lo(1),hi(1)
      klo = rr*k
      khi = rr*(k+1)
      if(khi .gt. hi(3)) then
         khi = klo
      endif
      w0_cell(i,j,k,3) =  HALF * (w0(klo) + w0(khi))
    end do
    end do
    end do

  end subroutine put_w0_on_3d_cells_cart

  subroutine put_w0_on_3d_cells_sphr(n,w0,w0_cell,normal,lo,hi,dx,ng)

    use bl_constants_module
    use geometry, nr_geometry => nr, z_geometry => z ! Because nr & z are local variables.

    integer        , intent(in   ) :: n,lo(:),hi(:),ng
    real(kind=dp_t), intent(  out) :: w0_cell(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) ::  normal(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer                  :: i,j,k,nr,index
    real(kind=dp_t)          :: x,y,z
    real(kind=dp_t)          :: radius,rfac,w0_cell_val

    nr = size(w0,dim=1)

    do k = lo(3),hi(3)
      z = (dble(k)+HALF)*dx(3) - center(3)
      do j = lo(2),hi(2)
        y = (dble(j)+HALF)*dx(2) - center(2)
        do i = lo(1),hi(1)
          x = (dble(i)+HALF)*dx(1) - center(1)
          radius = sqrt(x**2 + y**2 + z**2)
          index = int(radius / dr(n))
          if (index .lt. 0 .or. index .gt. nr-1) then
            print *,'RADIUS ',radius
            print *,'BOGUS INDEX IN PUT_ON_CELLS: ',index
            print *,'NOT IN RANGE 0 TO ',nr-1
            print *,'I J K ',i,j,k
            print *,'X Y Z ',x,y,z
            x = 1.0 / 0.0
            stop
          end if
          rfac = (radius - dble(index)*dr(n)) / dr(n)
          if (rfac .lt. 0.0 .or. rfac .gt. 1.0) then
            print *,'BAD RFAC ',rfac
            print *,'RADIUS, INDEX*DR ',radius, dble(index)*dr(n)
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

  end subroutine put_w0_on_3d_cells_sphr

end module fill_3d_module
