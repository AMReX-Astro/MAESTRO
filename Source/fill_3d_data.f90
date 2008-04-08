module fill_3d_module

  use bl_types
  use multifab_module

  implicit none

  private
  
  public :: make_3d_normal, put_1d_array_on_cart
  public :: fill_3d_data_c, fill_3d_data
  public :: put_1d_vector_on_3d_cells, put_1d_vector_on_3d_cells_sphr
  
contains

  subroutine make_3d_normal(normal,lo,hi,dx,ng)

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


  subroutine fill_3d_data_c(nlevs,dx,the_bc_level,mla,s0_cart,s0,bc_comp)

    use bl_prof_module
    use bl_constants_module
    use define_bc_module
    use ml_layout_module
    use multifab_physbc_module
    use ml_restriction_module, only: ml_cc_restriction_c
    use multifab_fill_ghost_module
    
    ! for spherical problems, this copies the base state onto a multifab
    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla
    type(multifab) , intent(inout) :: s0_cart(:)
    real(kind=dp_t), intent(in   ) :: s0(:,0:)
    integer        , intent(in   ) :: bc_comp

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
          call fill_3d_data(n,s0p(:,:,:,1),s0(n,:),lo,hi,dx(n,:),ng)
       end do
    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(s0_cart(nlevs),1,1)

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s0_cart(nlevs),1,bc_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(s0_cart(n-1),1,s0_cart(n),1,mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s0_cart(n),s0_cart(n-1),ng,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n),1,bc_comp,1)
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
  


  subroutine put_1d_vector_on_3d_cells(nlevs,w0,w0_cart,normal,dx,bc_comp,is_edge_centered, &
                                       the_bc_level,mla)

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
    integer        , intent(in   ) :: bc_comp
    logical        , intent(in   ) :: is_edge_centered
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla
    
    integer :: lo(w0_cart(1)%dim),hi(w0_cart(1)%dim)
    integer :: i,n,dm,ng
    real(kind=dp_t), pointer :: wp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_1d_vector_on_3d_cells")
    
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
             call put_1d_vector_on_3d_cells_sphr(n,is_edge_centered,w0(n,:),wp(:,:,:,:), &
                                                 np(:,:,:,:),lo,hi,dx(n,:),ng)
          else
             call put_1d_vector_on_3d_cells_cart(n,is_edge_centered,w0(n,:),wp(:,:,:,:), &
                                                 lo,hi,dx(n,dm),ng)
          end if
       end do

    enddo

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(w0_cart(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(w0_cart(nlevs),1,bc_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(w0_cart(n-1),1,w0_cart(n),1,mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(w0_cart(n),w0_cart(n-1),ng,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n),1,bc_comp,1)

       end do

    end if

    call destroy(bpt)
    
  end subroutine put_1d_vector_on_3d_cells
  
  subroutine put_1d_vector_on_3d_cells_cart(n,is_edge_centered,w0,w0_cell,lo,hi,dz,ng)

    use bl_constants_module
    use geometry, only: dr

    integer        , intent(in   ) :: n,lo(:),hi(:),ng
    logical        , intent(in   ) :: is_edge_centered
    real(kind=dp_t), intent(  out) :: w0_cell(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: dz

    integer :: i,j,k
    integer :: rr,klo,khi

    rr = int( dz / dr(n) + 1.d-12)

    w0_cell = ZERO

    if (is_edge_centered) then

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

    else

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                w0_cell(i,j,k,3) = w0(k)
             end do
          end do
       end do

    end if


  end subroutine put_1d_vector_on_3d_cells_cart

  subroutine put_1d_vector_on_3d_cells_sphr(n,is_edge_centered,w0,w0_cell,normal,lo,hi,dx,ng)

    use bl_constants_module
    use geometry, only: center, dr, nr, base_cc_loc

    integer        , intent(in   ) :: n,lo(:),hi(:),ng
    logical        , intent(in   ) :: is_edge_centered
    real(kind=dp_t), intent(  out) :: w0_cell(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) ::  normal(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer         :: i,j,k,index
    real(kind=dp_t) :: x,y,z
    real(kind=dp_t) :: radius,rfac,w0_cell_val

    if (is_edge_centered) then

       ! use linear interpolation between two nearest edge-centered values

       do k = lo(3),hi(3)
          z = (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2),hi(2)
             y = (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1),hi(1)
                x = (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(n))
                
                rfac = (radius - dble(index)*dr(n)) / dr(n)
                 
                w0_cell_val      = rfac * w0(index) + (ONE-rfac) * w0(index+1)
                w0_cell(i,j,k,1) = w0_cell_val * normal(i,j,k,1)
                w0_cell(i,j,k,2) = w0_cell_val * normal(i,j,k,2)
                w0_cell(i,j,k,3) = w0_cell_val * normal(i,j,k,3)
             end do
          end do
       end do

    else

       ! use linear interpolation between two nearest cell-centered values

       do k = lo(3),hi(3)
          z = (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2),hi(2)
             y = (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1),hi(1)
                x = (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(n))
                
                if (radius .ge. base_cc_loc(n,index)) then
                   if (index .eq. nr(n)-1) then
                      w0_cell_val = w0(index)
                   else
                      w0_cell_val = w0(index+1)*(radius-base_cc_loc(n,index))/dr(n) &
                           + w0(index)*(base_cc_loc(n,index+1)-radius)/dr(n)
                   endif
                else
                   if (index .eq. 0) then
                      w0_cell_val = w0(index)
                   else
                      w0_cell_val = w0(index)*(radius-base_cc_loc(n,index-1))/dr(n) &
                           + w0(index-1)*(base_cc_loc(n,index)-radius)/dr(n)
                   end if
                end if

                w0_cell(i,j,k,1) = w0_cell_val * normal(i,j,k,1)
                w0_cell(i,j,k,2) = w0_cell_val * normal(i,j,k,2)
                w0_cell(i,j,k,3) = w0_cell_val * normal(i,j,k,3)
             end do
          end do
       end do

    end if

  end subroutine put_1d_vector_on_3d_cells_sphr


  subroutine put_1d_array_on_cart(nlevs,s0,s0_cart,bc_comp,is_edge_centered,is_vector, &
                                  dx,the_bc_level,mla,interp_type,normal)

    use bl_prof_module
    use bl_constants_module
    use define_bc_module
    use geometry, only: spherical
    use ml_layout_module
    use multifab_physbc_module
    use ml_restriction_module, only: ml_cc_restriction_c
    use multifab_fill_ghost_module
    
    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(in   ) :: s0(:,0:)
    type(multifab) , intent(inout) :: s0_cart(:)
    integer        , intent(in   ) :: bc_comp
    logical        , intent(in   ) :: is_edge_centered,is_vector
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(inout) :: mla
    integer        , intent(in   ), optional :: interp_type
    type(multifab) , intent(in   ), optional :: normal(:)
    
    integer :: lo(s0_cart(1)%dim)
    integer :: hi(s0_cart(1)%dim)
    integer :: i,n,dm,ng
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_1d_array_on_cart")

    if (spherical .eq. 1 .and. is_vector .and. (.not. present(normal)) ) then
       call bl_error('Error: Calling put_1d_array_on_cart for spherical with an input vector and without normal')
    end if

    if (spherical .eq. 1 .and. (.not. present(interp_type)) ) then
       call bl_error('Error: Calling put_1d_array_on_cart for spherical without an interp_type')
    end if
    
    dm = s0_cart(1)%dim
    ng = s0_cart(1)%ng
    
    do n=1,nlevs
       
       do i = 1, s0_cart(n)%nboxes
          if ( multifab_remote(s0_cart(n), i) ) cycle
          sp => dataptr(s0_cart(n), i)
          lo =  lwb(get_box(s0_cart(n), i))
          hi =  upb(get_box(s0_cart(n), i))
          select case (dm)
          case (2)
             call put_1d_array_on_cart_2d(n,is_edge_centered,is_vector, &
                                          s0(n,:),sp(:,:,1,:),lo,hi,ng)
          case (3)
             if (spherical .eq. 0) then
                call put_1d_array_on_cart_3d(n,is_edge_centered,is_vector, &
                                             s0(n,:),sp(:,:,:,:),lo,hi,ng)
             else
                if (is_vector) then
                   np => dataptr(normal(n), i)
                end if
!                call put_1d_array_on_cart_3d_sphr(n,is_edge_centered, &
!                                                  is_vector,interp_type, &
!                                                  s0(n,:),sp(:,:,:,:),np(:,:,:,:), &
!                                                  lo,hi,dx(n,:),ng)
             end if
          end select
       end do

    enddo

    
    ! Warning: these boundary conditions need to be reworked to handle the is_vector cases

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(s0_cart(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s0_cart(nlevs),1,bc_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(s0_cart(n-1),1,s0_cart(n),1,mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s0_cart(n),s0_cart(n-1),ng,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n),1,bc_comp,1)

       end do

    end if

    call destroy(bpt)
    
  end subroutine put_1d_array_on_cart

  subroutine put_1d_array_on_cart_2d(n,is_edge_centered,is_vector,s0,s0_cart,lo,hi,ng)

    use bl_constants_module
    use geometry, only: dr

    integer        , intent(in   ) :: n
    integer        , intent(in   ) :: lo(:),hi(:),ng
    logical        , intent(in   ) :: is_edge_centered,is_vector
    real(kind=dp_t), intent(in   ) :: s0(0:)
    real(kind=dp_t), intent(inout) :: s0_cart(lo(1)-ng:,lo(2)-ng:,:)

    integer :: i,j

    s0_cart = ZERO

    if (is_edge_centered) then

       if (is_vector) then

          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                s0_cart(i,j,2) = HALF * (s0(j) + s0(j+1))
             end do
          end do

       else

          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                s0_cart(i,j,1) = HALF * (s0(j) + s0(j+1))
             end do
          end do

       end if

    else

       if (is_vector) then

          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                s0_cart(i,j,2) = s0(j)
             end do
          end do

       else

          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                s0_cart(i,j,1) = s0(j)
             end do
          end do

       end if

    end if

  end subroutine put_1d_array_on_cart_2d

  subroutine put_1d_array_on_cart_3d(n,is_edge_centered,is_vector,s0,s0_cart,lo,hi,ng)

    use bl_constants_module
    use geometry, only: dr

    integer        , intent(in   ) :: n
    integer        , intent(in   ) :: lo(:),hi(:),ng
    logical        , intent(in   ) :: is_edge_centered,is_vector
    real(kind=dp_t), intent(in   ) :: s0(0:)
    real(kind=dp_t), intent(inout) :: s0_cart(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)

    integer :: i,j,k

    s0_cart = ZERO

    if (is_edge_centered) then

       if (is_vector) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   s0_cart(i,j,k,3) = HALF * (s0(k) + s0(k+1))
                end do
             end do
          end do
          
       else
          
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   s0_cart(i,j,k,1) = HALF * (s0(k) + s0(k+1))
                end do
             end do
          end do
          
       end if

    else

       if (is_vector) then

          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   s0_cart(i,j,k,3) = s0(k)
                end do
             end do
          end do
          
       else
          
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                do i=lo(1),hi(1)
                   s0_cart(i,j,k,1) = s0(k)
                end do
             end do
          end do
          
       end if

    end if

  end subroutine put_1d_array_on_cart_3d

end module fill_3d_module
