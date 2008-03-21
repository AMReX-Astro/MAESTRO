module fill_3d_module

  use bl_types
  use multifab_module

  implicit none

  private

  public :: fill_3d_data_c, fill_3d_data
  public :: make_3d_normal, make_w0_cart, put_w0_on_3d_cells_sphr
  
contains

  subroutine fill_3d_data_c(nlevs,dx,the_bc_level,mla,s0_cart,s0,in_comp,bc_comp)

    use bl_prof_module
    use bl_constants_module
    use define_bc_module
    use ml_layout_module
    use multifab_physbc_module
    use ml_restriction_module, only: ml_cc_restriction_c
    use multifab_fill_ghost_module
    !
    ! for spherical problems, this copies the base state onto a multifab
    ! sames as the function fill_3d_data_wrap, except we assume
    ! start_comp = 1, num_comp = 1, and the base state only has one component
    !
    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(inout) :: s0_cart(:)
    real(kind=dp_t), intent(in   ) :: s0(:,0:)
    integer        , intent(in   ) :: in_comp,bc_comp

    integer                  :: i,ng,n
    integer                  :: lo(s0_cart(1)%dim),hi(s0_cart(1)%dim)
    real(kind=dp_t), pointer :: s0p(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "fill_3d_data_c")

    ng = s0_cart(1)%ng

    do n=1,nlevs
       do i=1,s0_cart(n)%nboxes
          if ( multifab_remote(s0_cart(n),i) ) cycle
          s0p => dataptr(s0_cart(n),i)
          lo = lwb(get_box(s0_cart(n),i))
          hi = upb(get_box(s0_cart(n),i))
          call fill_3d_data(n,s0p(:,:,:,in_comp),s0(n,:),lo,hi,dx(n,:),ng)
       end do
    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(s0_cart(nlevs),in_comp,1)

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s0_cart(nlevs),in_comp,bc_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(s0_cart(n-1),in_comp,s0_cart(n), &
                                   in_comp,mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s0_cart(n),s0_cart(n-1),ng,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n),in_comp,bc_comp,1)
    end do

 end if

    call destroy(bpt)

  end subroutine fill_3d_data_c

  subroutine fill_3d_data(n,data,s0,lo,hi,dx,ng)

    use bl_constants_module
    use geometry, only: center, dr, nr, base_cc_loc
    use bl_error_module
    
    integer        , intent(in   ) :: n,lo(:),hi(:),ng
    real(kind=dp_t), intent(  out) :: data(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(kind=dp_t), intent(in   ) ::   s0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    integer         :: i,j,k,index
    real(kind=dp_t) :: x,y,z
    real(kind=dp_t) :: radius
    logical         :: use_linear_interp

    use_linear_interp = .false.
    
    do k = lo(3),hi(3)
       z = (dble(k)+HALF)*dx(3) - center(3)
       do j = lo(2),hi(2)
          y = (dble(j)+HALF)*dx(2) - center(2)
          do i = lo(1),hi(1)
             x = (dble(i)+HALF)*dx(1) - center(1)
             radius = sqrt(x**2 + y**2 + z**2)
             index = int(radius / dr(n))

             if ( .false. ) then
                if (index .lt. 0 .or. index .gt. nr(n)-1) then
                   print *,'RADIUS ',radius
                   print *,'BOGUS INDEX IN FILL_3D: ',index
                   print *,'NOT IN RANGE 0 TO ',nr(n)-1
                   print *,'I J K ',i,j,k
                   print *,'X Y Z ',x,y,z
                   call bl_error(' ')
                end if
             end if

             if (use_linear_interp) then
                if(radius .ge. base_cc_loc(n,index)) then
                   if (index .eq. nr(n)-1) then
                      data(i,j,k) = s0(index)
                   else
                      data(i,j,k) = s0(index+1)*(radius-base_cc_loc(n,index))/dr(n) &
                           + s0(index)*(base_cc_loc(n,index+1)-radius)/dr(n)
                   end if
                else
                   if (index .eq. 0) then
                      data(i,j,k) = s0(index)
                   else
                      data(i,j,k) = s0(index)*(radius-base_cc_loc(n,index-1))/dr(n) &
                           + s0(index-1)*(base_cc_loc(n,index)-radius)/dr(n)
                   end if
                end if
             else
                data(i,j,k) = s0(index)
             end if

          end do
       end do
    end do
    
  end subroutine fill_3d_data
  
  subroutine make_3d_normal (normal,lo,hi,dx,ng)

    use bl_constants_module
    use geometry, only: spherical, center
    
    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(  out) :: normal(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)

    integer         :: i,j,k
    real(kind=dp_t) :: x,y,z,radius

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
      call bl_error('SHOULDNT CALL MAKE_3D_NORMAL WITH SPHERICAL = 0')
    end if

  end subroutine make_3d_normal

  subroutine make_w0_cart(nlevs,w0,w0_cart,normal,dx,the_bc_level,mla)

    use bl_prof_module
    use bl_constants_module
    use define_bc_module
    use geometry, only: spherical
    use ml_layout_module
    use multifab_physbc_module
    use ml_restriction_module, only: ml_cc_restriction_c
    use multifab_fill_ghost_module
    
    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(inout) :: w0_cart(:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla
    
    integer :: lo(w0_cart(1)%dim),hi(w0_cart(1)%dim)
    integer :: i,n,dm,ng
    real(kind=dp_t), pointer :: wp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_w0_cart")
    
    dm = w0_cart(1)%dim
    ng = w0_cart(1)%ng
    
    do n = 1, nlevs
    
       call setval(w0_cart(n),ZERO,all=.true.)
       
       do i = 1, w0_cart(n)%nboxes
          if ( multifab_remote(w0_cart(n), i) ) cycle
          wp => dataptr(w0_cart(n), i)
          lo =  lwb(get_box(w0_cart(n), i))
          hi =  upb(get_box(w0_cart(n), i))
          if (spherical .eq. 1) then
             np => dataptr(normal(n), i)
             call put_w0_on_3d_cells_sphr(n,w0(n,:),wp(:,:,:,:),np(:,:,:,:),lo,hi,dx(n,:),ng)
          else
             call put_w0_on_3d_cells_cart(n,w0(n,:),wp(:,:,:,:),lo,hi,dx(n,dm),ng)
          end if
       end do

    enddo

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(w0_cart(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(w0_cart(nlevs),1,dm,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(w0_cart(n-1),1,w0_cart(n),1,mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(w0_cart(n),w0_cart(n-1),ng,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n),1,dm,1)

       end do

    end if

    call destroy(bpt)
    
  end subroutine make_w0_cart
  
  subroutine put_w0_on_3d_cells_cart(n,w0,w0_cell,lo,hi,dz,ng)

    use bl_constants_module
    use geometry, only: dr

    integer        , intent(in   ) :: n,lo(:),hi(:),ng
    real(kind=dp_t), intent(  out) :: w0_cell(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: dz

    integer :: i,j,k
    integer :: rr,klo,khi

    rr = int( dz / dr(n) + 1.d-12)

    w0_cell = ZERO
    do k = lo(3),hi(3)
       klo = rr*k
       khi = rr*(k+1)
       if (khi .gt. hi(3)) khi = klo
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             w0_cell(i,j,k,3) =  HALF * (w0(klo) + w0(khi))
          end do
       end do
    end do

  end subroutine put_w0_on_3d_cells_cart

  subroutine put_w0_on_3d_cells_sphr(n,w0,w0_cell,normal,lo,hi,dx,ng)

    use bl_constants_module
    use geometry, only: center, dr

    integer        , intent(in   ) :: n,lo(:),hi(:),ng
    real(kind=dp_t), intent(  out) :: w0_cell(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) ::  normal(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer         :: i,j,k,nr,index
    real(kind=dp_t) :: x,y,z
    real(kind=dp_t) :: radius,rfac,w0_cell_val

    nr = size(w0,dim=1)

    do k = lo(3),hi(3)
      z = (dble(k)+HALF)*dx(3) - center(3)
      do j = lo(2),hi(2)
        y = (dble(j)+HALF)*dx(2) - center(2)
        do i = lo(1),hi(1)
          x = (dble(i)+HALF)*dx(1) - center(1)
          radius = sqrt(x**2 + y**2 + z**2)
          index  = int(radius / dr(n))

          if ( .false. ) then
             if (index .lt. 0 .or. index .gt. nr-1) then
                print *,'RADIUS ',radius
                print *,'BOGUS INDEX IN PUT_ON_CELLS: ',index
                print *,'NOT IN RANGE 0 TO ',nr-1
                print *,'I J K ',i,j,k
                print *,'X Y Z ',x,y,z
                call bl_error(' ')            
             end if
          end if

          rfac = (radius - dble(index)*dr(n)) / dr(n)

          if ( .false. ) then
             if (rfac .lt. 0.0 .or. rfac .gt. 1.0) then
                print *,'BAD RFAC ',rfac
                print *,'RADIUS, INDEX*DR ',radius, dble(index)*dr(n)
                call bl_error(' ')
             end if
          end if
          
          w0_cell_val      = rfac * w0(index) + (ONE-rfac) * w0(index+1)
          w0_cell(i,j,k,1) = w0_cell_val * normal(i,j,k,1)
          w0_cell(i,j,k,2) = w0_cell_val * normal(i,j,k,2)
          w0_cell(i,j,k,3) = w0_cell_val * normal(i,j,k,3)
        end do
      end do
    end do

  end subroutine put_w0_on_3d_cells_sphr

end module fill_3d_module
