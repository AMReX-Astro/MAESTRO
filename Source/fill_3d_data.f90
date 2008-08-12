module fill_3d_module

  use bl_types
  use multifab_module

  implicit none

  private
  
  public :: put_1d_array_on_cart,  put_1d_array_on_cart_3d_sphr
  public :: put_w0_on_edges, put_w0_on_edges_3d_sphr
  public :: make_3d_normal
  
contains  

  subroutine put_1d_array_on_cart(nlevs,s0,s0_cart,bc_comp,is_edge_centered,is_vector, &
                                  dx,the_bc_level,mla,normal)

    use bl_prof_module
    use bl_constants_module
    use define_bc_module
    use geometry, only: spherical
    use ml_layout_module
    use multifab_physbc_module
    use ml_restriction_module, only: ml_cc_restriction_c
    use multifab_fill_ghost_module
    use variables, only: foextrap_comp
    
    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(in   ) :: s0(:,0:)
    type(multifab) , intent(inout) :: s0_cart(:)
    integer        , intent(in   ) :: bc_comp
    logical        , intent(in   ) :: is_edge_centered,is_vector
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ), optional :: normal(:)
    
    integer :: lo(s0_cart(1)%dim)
    integer :: hi(s0_cart(1)%dim)
    integer :: i,n,dm,ng_s,ng_n,comp
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_1d_array_on_cart")

    if (spherical .eq. 1 .and. is_vector .and. (.not. present(normal)) ) then
       call bl_error('Error: Calling put_1d_array_on_cart for spherical with is_vector=T and without normal')
    end if

    dm = s0_cart(1)%dim
    ng_s = s0_cart(1)%ng
    
    do n=1,nlevs
       
       do i = 1, s0_cart(n)%nboxes
          if ( multifab_remote(s0_cart(n), i) ) cycle
          sp => dataptr(s0_cart(n), i)
          lo =  lwb(get_box(s0_cart(n), i))
          hi =  upb(get_box(s0_cart(n), i))
          select case (dm)

          case (2)
             call put_1d_array_on_cart_2d(is_edge_centered,is_vector, &
                                          s0(n,:),sp(:,:,1,:),lo,hi,ng_s)

          case (3)
             if (spherical .eq. 0) then
                call put_1d_array_on_cart_3d(is_edge_centered,is_vector, &
                                             s0(n,:),sp(:,:,:,:),lo,hi,ng_s)
             else
                if (is_vector) then
                   np => dataptr(normal(n), i)
                   ng_n = normal(n)%ng
                
                   call put_1d_array_on_cart_3d_sphr(n,is_edge_centered,is_vector, &
                                                     s0(n,:),sp(:,:,:,:), &
                                                     lo,hi,dx(n,:),ng_s,ng_n,np(:,:,:,:))
                else
                   call put_1d_array_on_cart_3d_sphr(n,is_edge_centered,is_vector, &
                                                     s0(n,:),sp(:,:,:,:), &
                                                     lo,hi,dx(n,:),ng_s,ng_n)
                end if
             endif

          end select
       end do

    enddo

    
    if (is_vector) then

       if (bc_comp .eq. foextrap_comp) then

          ! Here we fill each of the dm components using foextrap
          do comp=1,dm
             if (nlevs .eq. 1) then
                call multifab_fill_boundary_c(s0_cart(nlevs),comp,1)
                call multifab_physbc(s0_cart(nlevs),comp,bc_comp,1,the_bc_level(nlevs))
             else
                do n=nlevs,2,-1
                   call ml_cc_restriction_c(s0_cart(n-1),comp,s0_cart(n),comp, &
                                            mla%mba%rr(n-1,:),1)
                   call multifab_fill_ghost_cells(s0_cart(n),s0_cart(n-1),ng_s, &
                                                  mla%mba%rr(n-1,:),the_bc_level(n-1), &
                                                  the_bc_level(n),comp,bc_comp,1)
                end do
             end if
          end do

       else

          ! Here we fill each of the dm components using bc_comp+comp
          if (nlevs .eq. 1) then
             call multifab_fill_boundary_c(s0_cart(nlevs),1,dm)
             call multifab_physbc(s0_cart(nlevs),1,bc_comp,dm,the_bc_level(nlevs))
          else
             do n=nlevs,2,-1
                call ml_cc_restriction_c(s0_cart(n-1),1,s0_cart(n),1,mla%mba%rr(n-1,:),dm)
                call multifab_fill_ghost_cells(s0_cart(n),s0_cart(n-1),ng_s, &
                                               mla%mba%rr(n-1,:),the_bc_level(n-1), &
                                               the_bc_level(n),1,bc_comp,dm)
             end do
          end if

       end if

    else

       ! Here will fill the one component using bc_comp
       if (nlevs .eq. 1) then
          call multifab_fill_boundary_c(s0_cart(nlevs),1,1)
          call multifab_physbc(s0_cart(nlevs),1,bc_comp,1,the_bc_level(nlevs))
       else
          do n=nlevs,2,-1
             call ml_cc_restriction_c(s0_cart(n-1),1,s0_cart(n),1,mla%mba%rr(n-1,:),1)
             call multifab_fill_ghost_cells(s0_cart(n),s0_cart(n-1),ng_s,mla%mba%rr(n-1,:), &
                                            the_bc_level(n-1),the_bc_level(n),1,bc_comp,1)
          end do
       end if

    end if

    call destroy(bpt)
    
  end subroutine put_1d_array_on_cart

  subroutine put_1d_array_on_cart_2d(is_edge_centered,is_vector,s0,s0_cart,lo,hi,ng_s)

    use bl_constants_module
    use geometry, only: dr

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    logical        , intent(in   ) :: is_edge_centered,is_vector
    real(kind=dp_t), intent(in   ) :: s0(0:)
    real(kind=dp_t), intent(inout) :: s0_cart(lo(1)-ng_s:,lo(2)-ng_s:,:)

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

  subroutine put_1d_array_on_cart_3d(is_edge_centered,is_vector,s0,s0_cart,lo,hi,ng_s)

    use bl_constants_module
    use geometry, only: dr

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    logical        , intent(in   ) :: is_edge_centered,is_vector
    real(kind=dp_t), intent(in   ) :: s0(0:)
    real(kind=dp_t), intent(inout) :: s0_cart(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)

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

  subroutine put_1d_array_on_cart_3d_sphr(n,is_edge_centered,is_vector,s0,s0_cart,lo,hi, &
                                          dx,ng_s,ng_n,normal)

    ! note: ng_n is required only to dimension normal.  Since normal is 
    ! optional, if you do not pass normal in, then you can use any dummy 
    ! value for ng_n

    use bl_constants_module
    use geometry, only: dr, center, r_cc_loc, nr_fine
    use probin_module, only: interp_type_radial_bin_to_cart

    integer        , intent(in   ) :: n
    integer        , intent(in   ) :: lo(:),hi(:),ng_s, ng_n
    logical        , intent(in   ) :: is_edge_centered,is_vector
    real(kind=dp_t), intent(in   ) :: s0(0:)
    real(kind=dp_t), intent(inout) :: s0_cart(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ), optional :: normal(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)

    integer         :: i,j,k,index
    real(kind=dp_t) :: x,y,z
    real(kind=dp_t) :: radius,rfac,s0_cart_val

    if (is_vector .and. (.not. present(normal)) ) then
       call bl_error('Error: Calling put_1d_array_on_cart_3d_sphr with is_vector=T and without normal')
    end if

    if (is_edge_centered) then

       ! interpolate from radial bin edge values

       do k = lo(3),hi(3)
          z = (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2),hi(2)
             y = (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1),hi(1)
                x = (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(n))
                
                rfac = (radius - dble(index)*dr(n)) / dr(n)
                s0_cart_val      = rfac * s0(index) + (ONE-rfac) * s0(index+1)

                if (is_vector) then
                   s0_cart(i,j,k,1) = s0_cart_val * normal(i,j,k,1)
                   s0_cart(i,j,k,2) = s0_cart_val * normal(i,j,k,2)
                   s0_cart(i,j,k,3) = s0_cart_val * normal(i,j,k,3)
                else
                   s0_cart(i,j,k,1) = s0_cart_val
                end if
             end do
          end do
       end do

    else

       ! interpolate from radial bin centered values

       do k = lo(3),hi(3)
          z = (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2),hi(2)
             y = (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1),hi(1)
                x = (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(n))
                
                if (interp_type_radial_bin_to_cart .eq. 1) then

                   s0_cart_val = s0(index)

                else if (interp_type_radial_bin_to_cart .eq. 2) then

                   if (radius .ge. r_cc_loc(n,index)) then
                      if (index .eq. nr_fine-1) then
                         s0_cart_val = s0(index)
                      else
                         s0_cart_val = s0(index+1)*(radius-r_cc_loc(n,index))/dr(n) &
                              + s0(index)*(r_cc_loc(n,index+1)-radius)/dr(n)
                      endif
                   else
                      if (index .eq. 0) then
                         s0_cart_val = s0(index)
                      else
                         s0_cart_val = s0(index)*(radius-r_cc_loc(n,index-1))/dr(n) &
                              + s0(index-1)*(r_cc_loc(n,index)-radius)/dr(n)
                      end if
                   end if

                else
                   call bl_error('Error: interp_type_radial_bin_to_cart not defined')
                end if

                if (is_vector) then
                   s0_cart(i,j,k,1) = s0_cart_val * normal(i,j,k,1)
                   s0_cart(i,j,k,2) = s0_cart_val * normal(i,j,k,2)
                   s0_cart(i,j,k,3) = s0_cart_val * normal(i,j,k,3)
                else
                   s0_cart(i,j,k,1) = s0_cart_val
                end if
             end do
          end do
       end do

    end if

  end subroutine put_1d_array_on_cart_3d_sphr

  subroutine put_w0_on_edges(nlevs,w0,w0mac,dx,normal)

    use geometry, only: spherical

    integer        , intent(in   )           :: nlevs
    real(kind=dp_t), intent(in   )           :: w0(:,0:)
    type(multifab) , intent(inout)           :: w0mac(:)
    real(kind=dp_t), intent(in   )           :: dx(:,:)
    type(multifab) , intent(in   ), optional :: normal(:)

    integer :: lo(w0mac(1)%dim)
    integer :: hi(w0mac(1)%dim)
    integer :: i,n,dm,ng_w0,ng_n
    real(kind=dp_t), pointer :: w0mp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)
    
    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_w0_on_edges")

    if (spherical .eq. 1 .and. (.not. present(normal)) ) then
       call bl_error('Error: Calling put_w0_on_edges for spherical without normal')
    end if

    dm = w0mac(1)%dim
    ng_w0 = w0mac(1)%ng
    ng_n = normal(1)%ng

    do n=1,nlevs
       do i=1,w0mac(n)%nboxes
          if ( multifab_remote(w0mac(n), i) ) cycle
          w0mp => dataptr(w0mac(n), i)
          np   => dataptr(normal(n), i)
          lo = lwb(get_box(w0mac(n), i))
          hi = upb(get_box(w0mac(n), i))
          select case (dm)
          case (2)
             call bl_error('Error: Should not have to call put_w0_on_edges in 2d')
          case (3)
             if (spherical .eq. 0) then
                call put_w0_on_edges_3d_sphr(n,w0(n,:),w0mp(:,:,:,:),ng_w0,np(:,:,:,:), &
                                             ng_n,lo,hi,dx(n,:))
             else
                call bl_error('Error: Should not have to call put_w0_on_edges in 3d cart')
             end if
          end select
       end do
    end do

  end subroutine put_w0_on_edges
  
  subroutine put_w0_on_edges_3d_sphr(n,w0,w0mac,ng_w0,normal,ng_n,lo,hi,dx)

    use bl_constants_module
    use geometry, only: dr, center

    integer        , intent(in   ) :: n,lo(:),hi(:),ng_w0,ng_n
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(inout) ::  w0mac(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:,:)
    real(kind=dp_t), intent(inout) :: normal(lo(1)-ng_n :,lo(2)-ng_n :,lo(3)-ng_n :,:)    
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer         :: i,j,k,index,w0mac_interp_type
    real(kind=dp_t) :: x,y,z
    real(kind=dp_t) :: radius,w0_cart_val,rfac
    real(kind=dp_t), allocatable :: w0_cc(:,:,:,:)
    real(kind=dp_t), allocatable :: w0_nodal(:,:,:,:)

    ! we currently have three different ideas for computing w0mac
    ! 1.  Interpolate w0 to cell centers, then average to edges
    ! 2.  Interpolate w0 to edges directly
    ! 3.  Interpolate w0 to nodes, then average to edges

    w0mac_interp_type = 1

    if (w0mac_interp_type .eq. 1) then

       allocate(w0_cc(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2,3))

       do k = lo(3)-2,hi(3)+2
          z = (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-2,hi(2)+2
             y = (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-2,hi(1)+2
                x = (dble(i)+HALF)*dx(1) - center(1)

                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(n))
                
                rfac = (radius - dble(index)*dr(n)) / dr(n)
                w0_cart_val = rfac * w0(index) + (ONE-rfac) * w0(index+1)

                w0_cc(i,j,k,1) = w0_cart_val * normal(i,j,k,1)
                w0_cc(i,j,k,2) = w0_cart_val * normal(i,j,k,2)
                w0_cc(i,j,k,3) = w0_cart_val * normal(i,j,k,3)

             end do
          end do
       end do

       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+2
                w0mac(i,j,k,1) = HALF* (w0_cc(i-1,j,k,1) + w0_cc(i,j,k,1))
             end do
          end do
       end do

       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+2
             do i=lo(1)-1,hi(1)+1
                w0mac(i,j,k,2) = HALF* (w0_cc(i,j-1,k,2) + w0_cc(i,j,k,2))
             end do
          end do
       end do

       do k=lo(3)-1,hi(3)+2
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                w0mac(i,j,k,3) = HALF* (w0_cc(i,j,k-1,3) + w0_cc(i,j,k,3))
             end do
          end do
       end do

    else if (w0mac_interp_type .eq. 2) then

       do k = lo(3)-1,hi(3)+1
          z = (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+2
                x = (dble(i)     )*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(n))

                rfac = (radius - dble(index)*dr(n)) / dr(n)
                w0_cart_val = rfac * w0(index) + (ONE-rfac) * w0(index+1)

                w0mac(i,j,k,1) = w0_cart_val * x / radius

             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+1
          z = (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+2
             y = (dble(j)     )*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(n))

                rfac = (radius - dble(index)*dr(n)) / dr(n)
                w0_cart_val = rfac * w0(index) + (ONE-rfac) * w0(index+1)

                w0mac(i,j,k,2) = w0_cart_val * y / radius

             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+2
          z = (dble(k)     )*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(n))

                rfac = (radius - dble(index)*dr(n)) / dr(n)
                w0_cart_val = rfac * w0(index) + (ONE-rfac) * w0(index+1)

                w0mac(i,j,k,3) = w0_cart_val * z / radius

             end do
          end do
       end do

    else if (w0mac_interp_type .eq. 3) then

       allocate(w0_nodal(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,lo(3)-1:hi(3)+2,3))

       do k = lo(3)-1,hi(3)+2
          z = (dble(k))*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+2
             y = (dble(j))*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+2
                x = (dble(i))*dx(1) - center(1)

                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(n))
                
                rfac = (radius - dble(index)*dr(n)) / dr(n)
                w0_cart_val = rfac * w0(index) + (ONE-rfac) * w0(index+1)

                w0_nodal(i,j,k,1) = w0_cart_val * x / radius
                w0_nodal(i,j,k,2) = w0_cart_val * y / radius
                w0_nodal(i,j,k,3) = w0_cart_val * z / radius

             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+1
          do j = lo(2)-1,hi(2)+1
             do i = lo(1)-1,hi(1)+2
                w0mac(i,j,k,1) = FOURTH*( w0_nodal(i,j,k,1) + w0_nodal(i,j+1,k,1) &
                                         +w0_nodal(i,j,k+1,1) + w0_nodal(i,j+1,k+1,1))
             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+1
          do j = lo(2)-1,hi(2)+2
             do i = lo(1)-1,hi(1)+1
                w0mac(i,j,k,2) = FOURTH*( w0_nodal(i,j,k,2) + w0_nodal(i+1,j,k,2) &
                                         +w0_nodal(i,j,k+1,2) + w0_nodal(i+1,j,k+1,2))
             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+2
          do j = lo(2)-1,hi(2)+1
             do i = lo(1)-1,hi(1)+1
                w0mac(i,j,k,2) = FOURTH*( w0_nodal(i,j,k,3) + w0_nodal(i+1,j,k,3) &
                                         +w0_nodal(i,j+1,k,3) + w0_nodal(i+1,j+1,k,3))
             end do
          end do
       end do

    else
       call bl_error('Error: fill_3d_data:w0mac_interp_type > 3')
    end if

  end subroutine put_w0_on_edges_3d_sphr

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

end module fill_3d_module
