module fill_3d_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use bl_constants_module
  use bl_prof_module

  implicit none

  private
  
  public :: put_1d_array_on_cart,  put_1d_array_on_cart_3d_sphr
  public :: put_w0_on_edges, put_s0_on_edges
  public :: make_normal, make_normal_3d_sphr
  public :: put_data_on_faces
  
contains  

  subroutine put_1d_array_on_cart(s0,s0_cart,bc_comp,is_input_edge_centered, &
                                  is_output_a_vector,dx,the_bc_level,mla,normal)

    use bl_constants_module
    use define_bc_module
    use geometry, only: spherical, dm, nlevs
    use ml_layout_module
    use multifab_physbc_module
    use ml_restriction_module, only: ml_cc_restriction_c
    use multifab_fill_ghost_module
    use variables, only: foextrap_comp
    
    real(kind=dp_t), intent(in   ) :: s0(:,0:)
    type(multifab) , intent(inout) :: s0_cart(:)
    integer        , intent(in   ) :: bc_comp
    logical        , intent(in   ) :: is_input_edge_centered,is_output_a_vector
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ), optional :: normal(:)
    
    integer :: lo(dm)
    integer :: hi(dm)
    integer :: i,n,ng_s,ng_n,comp
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_1d_array_on_cart")

    if (spherical .eq. 1 .and. is_output_a_vector .and. (.not. present(normal)) ) then
       call bl_error('Error: Calling put_1d_array_on_cart for spherical with is_output_a_vector=T and without normal')
    end if

    ng_s = s0_cart(1)%ng
    
    do n=1,nlevs
       
       do i = 1, s0_cart(n)%nboxes
          if ( multifab_remote(s0_cart(n), i) ) cycle
          sp => dataptr(s0_cart(n), i)
          lo =  lwb(get_box(s0_cart(n), i))
          hi =  upb(get_box(s0_cart(n), i))
          select case (dm)
          case (2)
             call put_1d_array_on_cart_2d(is_input_edge_centered,is_output_a_vector, &
                                          s0(n,:),sp(:,:,1,:),lo,hi,ng_s)
          case (3)
             if (spherical .eq. 0) then
                call put_1d_array_on_cart_3d(is_input_edge_centered,is_output_a_vector, &
                                             s0(n,:),sp(:,:,:,:),lo,hi,ng_s)
             else
                if (is_output_a_vector) then
                   np => dataptr(normal(n), i)
                   ng_n = normal(n)%ng
               
                   call put_1d_array_on_cart_3d_sphr(is_input_edge_centered, &
                                                     is_output_a_vector,s0(1,:), &
                                                     sp(:,:,:,:),lo,hi,dx(n,:),ng_s,ng_n, &
                                                     np(:,:,:,:))
                else
                   call put_1d_array_on_cart_3d_sphr(is_input_edge_centered, &
                                                     is_output_a_vector,s0(1,:), &
                                                     sp(:,:,:,:),lo,hi,dx(n,:),ng_s,ng_n)
                end if
             endif
          end select
       end do

    enddo
    
    if (is_output_a_vector) then

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
                                                  the_bc_level(n),comp,bc_comp,1, &
                                                  fill_crse_input=.false.)
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
                                               the_bc_level(n),1,bc_comp,dm, &
                                               fill_crse_input=.false.)
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
                                            the_bc_level(n-1),the_bc_level(n),1,bc_comp,1, &
                                            fill_crse_input=.false.)
          end do
       end if

    end if

    call destroy(bpt)
    
  end subroutine put_1d_array_on_cart

  subroutine put_1d_array_on_cart_2d(is_input_edge_centered,is_output_a_vector,s0,s0_cart, &
                                     lo,hi,ng_s)

    use bl_constants_module
    use geometry, only: dr

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    logical        , intent(in   ) :: is_input_edge_centered,is_output_a_vector
    real(kind=dp_t), intent(in   ) :: s0(0:)
    real(kind=dp_t), intent(inout) :: s0_cart(lo(1)-ng_s:,lo(2)-ng_s:,:)

    integer :: i,j

    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_1d_array_on_cart_2d")

    s0_cart = ZERO

    if (is_input_edge_centered) then

       if (is_output_a_vector) then

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

       if (is_output_a_vector) then

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

    call destroy(bpt)

  end subroutine put_1d_array_on_cart_2d

  subroutine put_1d_array_on_cart_3d(is_input_edge_centered,is_output_a_vector,s0,s0_cart, &
                                     lo,hi,ng_s)

    use bl_constants_module
    use geometry, only: dr

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    logical        , intent(in   ) :: is_input_edge_centered,is_output_a_vector
    real(kind=dp_t), intent(in   ) :: s0(0:)
    real(kind=dp_t), intent(inout) :: s0_cart(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)

    integer :: i,j,k

    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_1d_array_on_cart_3d")

    s0_cart = ZERO

    if (is_input_edge_centered) then

       if (is_output_a_vector) then

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

       if (is_output_a_vector) then

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

    call destroy(bpt)

  end subroutine put_1d_array_on_cart_3d

  subroutine put_1d_array_on_cart_3d_sphr(is_input_edge_centered,is_output_a_vector, &
                                          s0,s0_cart,lo,hi,dx,ng_s,ng_n,normal)

    ! note: ng_n is required only to dimension normal.  Since normal is 
    ! optional, if you do not pass normal in, then you can use any dummy 
    ! value for ng_n

    use bl_constants_module
    use geometry, only: dr, center, r_cc_loc, nr_fine
    use probin_module, only: s0_interp_type

    integer        , intent(in   ) :: lo(:),hi(:),ng_s, ng_n
    logical        , intent(in   ) :: is_input_edge_centered,is_output_a_vector
    real(kind=dp_t), intent(in   ) :: s0(0:)
    real(kind=dp_t), intent(inout) :: s0_cart(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ), optional :: normal(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)

    integer         :: i,j,k,index
    real(kind=dp_t) :: x,y,z
    real(kind=dp_t) :: radius,rfac,s0_cart_val

    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_1d_array_on_cart_3d_sphr")

    if (is_output_a_vector .and. (.not. present(normal)) ) then
       call bl_error('Error: Calling put_1d_array_on_cart_3d_sphr with is_output_a_vector=T and without normal')
    end if

    if (is_input_edge_centered) then

       ! interpolate from radial bin edge values

       do k = lo(3),hi(3)
          z = (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2),hi(2)
             y = (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1),hi(1)
                x = (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(1))
                
                rfac = (radius - dble(index)*dr(1)) / dr(1)
                s0_cart_val      = rfac * s0(index) + (ONE-rfac) * s0(index+1)

                if (is_output_a_vector) then
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
                index  = int(radius / dr(1))
                
                if (s0_interp_type .eq. 1) then

                   s0_cart_val = s0(index)

                else if (s0_interp_type .eq. 2) then

                   if (radius .ge. r_cc_loc(1,index)) then
                      if (index .eq. nr_fine-1) then
                         s0_cart_val = s0(index)
                      else
                         s0_cart_val = s0(index+1)*(radius-r_cc_loc(1,index))/dr(1) &
                              + s0(index)*(r_cc_loc(1,index+1)-radius)/dr(1)
                      endif
                   else
                      if (index .eq. 0) then
                         s0_cart_val = s0(index)
                      else
                         s0_cart_val = s0(index)*(radius-r_cc_loc(1,index-1))/dr(1) &
                              + s0(index-1)*(r_cc_loc(1,index)-radius)/dr(1)
                      end if
                   end if

                else
                   call bl_error('Error: s0_interp_type not defined')
                end if

                if (is_output_a_vector) then
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

    call destroy(bpt)

  end subroutine put_1d_array_on_cart_3d_sphr

  subroutine put_w0_on_edges(mla,w0,w0mac,dx)

    use bl_constants_module
    use geometry, only: spherical, nr_fine, dm, nlevs
    use probin_module, only: w0mac_interp_type
    use variables, only: foextrap_comp,press_comp
    use define_bc_module
    use macproject_module
    use fabio_module

    type(ml_layout), intent(in   ) :: mla
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(inout) :: w0mac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    ! Local variables
    integer         :: lo(dm),hi(dm)
    integer         :: i,n,ng_w0
    real(kind=dp_t) :: w0rhs(mla%nlevel,0:nr_fine-1)

    ! Local pointers
    real(kind=dp_t), pointer :: w0xp(:,:,:,:)
    real(kind=dp_t), pointer :: w0yp(:,:,:,:)
    real(kind=dp_t), pointer :: w0zp(:,:,:,:)
    
    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_w0_on_edges")

    if (dm .eq. 2 .or. spherical .eq. 0) then
       call bl_error('Error: only call put_w0_on_edges for spherical')
    end if

    ng_w0 = w0mac(1,1)%ng
    
    do n=1,nlevs
       do i=1,w0mac(n,1)%nboxes
          if ( multifab_remote(w0mac(n,1), i) ) cycle
          w0xp => dataptr(w0mac(n,1), i)
          w0yp => dataptr(w0mac(n,2), i)
          w0zp => dataptr(w0mac(n,3), i)
          lo = lwb(get_box(w0mac(n,1), i))
          hi = upb(get_box(w0mac(n,1), i))
          call put_w0_on_edges_3d_sphr(w0(1,:),w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1), &
                                       ng_w0,lo,hi,dx(n,:))
       end do
    end do

    call destroy(bpt)

  end subroutine put_w0_on_edges

  subroutine put_w0_on_edges_3d_sphr(w0,w0macx,w0macy,w0macz,ng_w0,lo,hi,dx)

    use bl_constants_module
    use geometry, only: dr, center, nr_fine
    use probin_module, only: w0mac_interp_type

    integer        , intent(in   ) :: lo(:),hi(:),ng_w0
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(inout) :: w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(inout) :: w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(inout) :: w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer         :: i,j,k,index
    real(kind=dp_t) :: x,y,z
    real(kind=dp_t) :: radius,w0_cart_val,rfac
    real(kind=dp_t), allocatable :: w0_cc(:,:,:,:)
    real(kind=dp_t), allocatable :: w0_nodal(:,:,:,:)

    if (ng_w0 .ne. 1) then
       call bl_error('Error: put_w0_on_edges_3d_sphr assumes one ghost cell')
    end if

    ! we currently have three different ideas for computing w0mac
    ! 1.  Interpolate w0 to cell centers, then average to edges
    ! 2.  Interpolate w0 to edges directly
    ! 3.  Interpolate w0 to nodes, then average to edges

    if (w0mac_interp_type .eq. 1) then

       allocate(w0_cc(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2,lo(3)-2:hi(3)+2,3))

       do k = lo(3)-2,hi(3)+2
          z = (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-2,hi(2)+2
             y = (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-2,hi(1)+2
                x = (dble(i)+HALF)*dx(1) - center(1)

                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(1))
                
                rfac = (radius - dble(index)*dr(1)) / dr(1)

                if (index .lt. nr_fine) then
                   w0_cart_val = rfac * w0(index) + (ONE-rfac) * w0(index+1)
                else
                   w0_cart_val = w0(nr_fine)
                end if

                w0_cc(i,j,k,1) = w0_cart_val * x / radius
                w0_cc(i,j,k,2) = w0_cart_val * y / radius
                w0_cc(i,j,k,3) = w0_cart_val * z / radius

             end do
          end do
       end do

       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+2
                w0macx(i,j,k) = HALF* (w0_cc(i-1,j,k,1) + w0_cc(i,j,k,1))
             end do
          end do
       end do

       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+2
             do i=lo(1)-1,hi(1)+1
                w0macy(i,j,k) = HALF* (w0_cc(i,j-1,k,2) + w0_cc(i,j,k,2))
             end do
          end do
       end do

       do k=lo(3)-1,hi(3)+2
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                w0macz(i,j,k) = HALF* (w0_cc(i,j,k-1,3) + w0_cc(i,j,k,3))
             end do
          end do
       end do

       deallocate(w0_cc)

    else if (w0mac_interp_type .eq. 2) then

       do k = lo(3)-1,hi(3)+1
          z = (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+2
                x = (dble(i)     )*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(1))

                rfac = (radius - dble(index)*dr(1)) / dr(1)

                if (index .lt. nr_fine) then
                   w0_cart_val = rfac * w0(index) + (ONE-rfac) * w0(index+1)
                else
                   w0_cart_val = w0(nr_fine)
                end if

                w0macx(i,j,k) = w0_cart_val * x / radius

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
                index  = int(radius / dr(1))

                rfac = (radius - dble(index)*dr(1)) / dr(1)

                if (index .lt. nr_fine) then
                   w0_cart_val = rfac * w0(index) + (ONE-rfac) * w0(index+1)
                else
                   w0_cart_val = w0(nr_fine)
                end if

                w0macy(i,j,k) = w0_cart_val * y / radius

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
                index  = int(radius / dr(1))

                rfac = (radius - dble(index)*dr(1)) / dr(1)

                if (index .lt. nr_fine) then
                   w0_cart_val = rfac * w0(index) + (ONE-rfac) * w0(index+1)
                else
                   w0_cart_val = w0(nr_fine)
                end if

                w0macz(i,j,k) = w0_cart_val * z / radius

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
                index  = int(radius / dr(1))
                
                rfac = (radius - dble(index)*dr(1)) / dr(1)

                if (index .lt. nr_fine) then
                   w0_cart_val = rfac * w0(index) + (ONE-rfac) * w0(index+1)
                else
                   w0_cart_val = w0(nr_fine)
                end if

                w0_nodal(i,j,k,1) = w0_cart_val * x / radius
                w0_nodal(i,j,k,2) = w0_cart_val * y / radius
                w0_nodal(i,j,k,3) = w0_cart_val * z / radius

             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+1
          do j = lo(2)-1,hi(2)+1
             do i = lo(1)-1,hi(1)+2
                w0macx(i,j,k) = FOURTH*( w0_nodal(i,j,k,1) + w0_nodal(i,j+1,k,1) &
                                        +w0_nodal(i,j,k+1,1) + w0_nodal(i,j+1,k+1,1))
             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+1
          do j = lo(2)-1,hi(2)+2
             do i = lo(1)-1,hi(1)+1
                w0macy(i,j,k) = FOURTH*( w0_nodal(i,j,k,2) + w0_nodal(i+1,j,k,2) &
                                        +w0_nodal(i,j,k+1,2) + w0_nodal(i+1,j,k+1,2))
             end do
          end do
       end do

       do k = lo(3)-1,hi(3)+2
          do j = lo(2)-1,hi(2)+1
             do i = lo(1)-1,hi(1)+1
                w0macz(i,j,k) = FOURTH*( w0_nodal(i,j,k,3) + w0_nodal(i+1,j,k,3) &
                                        +w0_nodal(i,j+1,k,3) + w0_nodal(i+1,j+1,k,3))
             end do
          end do
       end do

       deallocate(w0_nodal)

    else
       call bl_error('Error: fill_3d_data:w0mac_interp_type can only be 1,2 or 3')
    end if

  end subroutine put_w0_on_edges_3d_sphr

  subroutine put_s0_on_edges(mla,s0,s0mac,dx)

    use bl_constants_module
    use geometry, only: spherical, nr_fine, dm, nlevs
    use variables, only: foextrap_comp,press_comp
    use define_bc_module
    use macproject_module
    use fabio_module

    type(ml_layout), intent(in   ) :: mla
    real(kind=dp_t), intent(in   ) :: s0(:,0:)
    type(multifab) , intent(inout) :: s0mac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)

    ! Local variables
    integer         :: lo(dm),hi(dm)
    integer         :: i,n,ng_s0

    ! Local pointers
    real(kind=dp_t), pointer :: s0xp(:,:,:,:)
    real(kind=dp_t), pointer :: s0yp(:,:,:,:)
    real(kind=dp_t), pointer :: s0zp(:,:,:,:)
    
    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_s0_on_edges")

    if (dm .eq. 2 .or. spherical .eq. 0) then
       call bl_error('Error: only call put_s0_on_edges for spherical')
    end if

    ng_s0 = s0mac(1,1)%ng
    
    do n=1,nlevs
       do i=1,s0mac(n,1)%nboxes
          if ( multifab_remote(s0mac(n,1), i) ) cycle
          s0xp => dataptr(s0mac(n,1), i)
          s0yp => dataptr(s0mac(n,2), i)
          s0zp => dataptr(s0mac(n,3), i)
          lo = lwb(get_box(s0mac(n,1), i))
          hi = upb(get_box(s0mac(n,1), i))
          call put_s0_on_edges_3d_sphr(s0(1,:),s0xp(:,:,:,1),s0yp(:,:,:,1), &
                                       s0zp(:,:,:,1),ng_s0,lo,hi,dx(n,:))
       end do
    end do

    call destroy(bpt)

  end subroutine put_s0_on_edges

  subroutine put_s0_on_edges_3d_sphr(s0,s0macx,s0macy,s0macz,ng_s0,lo,hi,dx)

    use bl_constants_module
    use geometry, only: dr, center, nr_fine
    use probin_module, only: s0mac_interp_type

    integer        , intent(in   ) :: lo(:),hi(:),ng_s0
    real(kind=dp_t), intent(in   ) :: s0(0:)
    real(kind=dp_t), intent(inout) :: s0macx(lo(1)-ng_s0:,lo(2)-ng_s0:,lo(3)-ng_s0:)
    real(kind=dp_t), intent(inout) :: s0macy(lo(1)-ng_s0:,lo(2)-ng_s0:,lo(3)-ng_s0:)
    real(kind=dp_t), intent(inout) :: s0macz(lo(1)-ng_s0:,lo(2)-ng_s0:,lo(3)-ng_s0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer         :: i,j,k,index
    real(kind=dp_t) :: x,y,z
    real(kind=dp_t) :: radius,s0_cart_val,rfac
    real(kind=dp_t), allocatable :: s0_nodal(:,:,:,:)

    if (ng_s0 .ne. 1) then
       call bl_error('Error: put_s0_on_edges_3d_sphr assumes one ghost cell')
    end if

    ! we currently have three different ideas for computing s0mac
    ! 1.  Interpolate s0 to cell centers, then average to edges
    ! 2.  Interpolate s0 to edges directly
    ! 3.  Interpolate s0 to nodes, then average to edges

    if (s0mac_interp_type .eq. 1) then

       ! writing this is slightly more compicated since we need to create a 
       ! full multifab with a ghost cell in order to average to faces.
       ! we do this already in mkflux.f90, for example, without using this function.
       call bl_error('Error: put_s0_on_edges_3d_sphr with s0mac_interp_type=1 not written yet')
      
    else if (s0mac_interp_type .eq. 2) then

       do k = lo(3)-1,hi(3)+1
          z = (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+2
                x = (dble(i)     )*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int((radius-HALF*dr(1)) / dr(1))

                rfac = ((radius-HALF*dr(1)) - dble(index)*dr(1)) / dr(1)

                if (index .lt. nr_fine) then
                   s0macx(i,j,k) = rfac * s0(index) + (ONE-rfac) * s0(index+1)
                else
                   s0macx(i,j,k) = s0(nr_fine)
                end if

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
                index  = int((radius-HALF*dr(1)) / dr(1))

                rfac = ((radius-HALF*dr(1)) - dble(index)*dr(1)) / dr(1)

                if (index .lt. nr_fine) then
                   s0macy(i,j,k) = rfac * s0(index) + (ONE-rfac) * s0(index+1)
                else
                   s0macy(i,j,k) = s0(nr_fine)
                end if

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
                index  = int((radius-HALF*dr(1)) / dr(1))

                rfac = ((radius-HALF*dr(1)) - dble(index)*dr(1)) / dr(1)

                if (index .lt. nr_fine) then
                   s0macz(i,j,k) = rfac * s0(index) + (ONE-rfac) * s0(index+1)
                else
                   s0macz(i,j,k) = s0(nr_fine)
                end if

             end do
          end do
       end do

    else if (s0mac_interp_type .eq. 3) then

       call bl_error('Error: put_s0_on_edges_3d_sphr with s0mac_interp_type=3 not written yet')

    else

       call bl_error('Error: fill_3d_data:s0mac_interp_type can only be 1,2 or 3')

    end if

  end subroutine put_s0_on_edges_3d_sphr

  subroutine make_normal(normal,dx)

    use geometry, only: spherical, dm, nlevs

    type(multifab) , intent(inout) :: normal(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
        
    integer             :: lo(dm),hi(dm)
    integer             :: n,i,ng_n
    real(dp_t), pointer :: nop(:,:,:,:)
        
    ng_n = normal(1)%ng

    if (spherical .eq. 1) then
       do n = 1,nlevs
          do i = 1, normal(n)%nboxes
             if ( multifab_remote(normal(n), i) ) cycle
             nop => dataptr(normal(n), i)
             lo =  lwb(get_box(normal(n), i))
             hi =  upb(get_box(normal(n), i))
             call make_normal_3d_sphr(nop(:,:,:,:),lo,hi,dx(n,:),ng_n)
          end do
       end do
    end if

  end subroutine make_normal

  subroutine make_normal_3d_sphr(normal,lo,hi,dx,ng)

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

  end subroutine make_normal_3d_sphr


  subroutine put_data_on_faces(ccfab,comp,beta,harmonic_avg)

    use geometry, only: nlevs, dm

    type(multifab) , intent(in   ) :: ccfab(:)
    type(multifab) , intent(inout) :: beta(:)
    integer        , intent(in   ) :: comp
    logical        , intent(in   ) :: harmonic_avg

    ! local
    integer :: n,i,ng_cc,ng_fc
    integer :: lo(dm),hi(dm)

    real(kind=dp_t), pointer :: ccfabp(:,:,:,:),betap(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_data_on_faces")

    ng_cc = ccfab(1)%ng
    ng_fc = beta(1)%ng

    ! setup beta = ccfab on faces
    do n=1,nlevs
       do i=1,beta(n)%nboxes
          if (multifab_remote(beta(n),i)) cycle
          ccfabp => dataptr(ccfab(n),i)
          betap   => dataptr(beta(n),i)
          lo = lwb(get_box(beta(n),i))
          hi = upb(get_box(beta(n),i))
          select case (dm)
          case (2)
             call put_data_on_faces_2d(lo,ccfabp(:,:,1,comp),ng_cc,betap(:,:,1,:),ng_fc, &
                                       harmonic_avg)
          case (3)
             call put_data_on_faces_3d(lo,ccfabp(:,:,:,comp),ng_cc,betap(:,:,:,:),ng_fc, &
                                       harmonic_avg)
          end select
       end do
    enddo

    call destroy(bpt)

  end subroutine put_data_on_faces


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! put beta on faces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine put_data_on_faces_2d(lo,ccbeta,ng_cc,beta,ng_fc,harmonic_avg)

    integer        , intent(in   ) :: lo(:), ng_cc, ng_fc
    real(kind=dp_t), intent(in   ) :: ccbeta(lo(1)-ng_cc:,lo(2)-ng_cc:)
    real(kind=dp_t), intent(inout) ::   beta(lo(1)-ng_fc:,lo(2)-ng_fc:,:)
    logical        , intent(in   ) :: harmonic_avg

    ! Local
    integer :: i,j
    integer :: nx,ny

    nx = size(beta,dim=1) - 2
    ny = size(beta,dim=2) - 2

    if (harmonic_avg) then

       do j = lo(2),lo(2)+ny-1
          do i = lo(1),lo(1)+nx
             beta(i,j,1) = TWO*(ccbeta(i,j)*ccbeta(i-1,j))/(ccbeta(i,j) + ccbeta(i-1,j))
          end do
       end do

       do j = lo(2),lo(2)+ny
          do i = lo(1),lo(1)+nx-1
             beta(i,j,2) = TWO*(ccbeta(i,j)*ccbeta(i,j-1))/(ccbeta(i,j) + ccbeta(i,j-1))
          end do
       end do

    else

       call bl_error('put_data_on_faces_2d non-harmonic not written yet')

    end if

  end subroutine put_data_on_faces_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! put beta on faces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine put_data_on_faces_3d(lo,ccbeta,ng_cc,beta,ng_fc,harmonic_avg)

    integer        , intent(in   ) :: lo(:), ng_cc, ng_fc
    real(kind=dp_t), intent(in   ) :: ccbeta(lo(1)-ng_cc:,lo(2)-ng_cc:,lo(3)-ng_cc:)
    real(kind=dp_t), intent(inout) ::   beta(lo(1)-ng_fc:,lo(2)-ng_fc:,lo(3)-ng_fc:,:)
    logical        , intent(in   ) :: harmonic_avg

    ! Local
    integer :: i,j,k
    integer :: nx,ny,nz

    nx = size(beta,dim=1) - 2
    ny = size(beta,dim=2) - 2
    nz = size(beta,dim=3) - 2

    if (harmonic_avg) then

       do k = lo(3),lo(3)+nz-1
          do j = lo(2),lo(2)+ny-1
             do i = lo(1),lo(1)+nx
                beta(i,j,k,1) = TWO*(ccbeta(i,j,k)*ccbeta(i-1,j,k))/(ccbeta(i,j,k) &
                     + ccbeta(i-1,j,k))
             end do
          end do
       end do

       do k = lo(3),lo(3)+nz-1
          do j = lo(2),lo(2)+ny
             do i = lo(1),lo(1)+nx-1
                beta(i,j,k,2) = TWO*(ccbeta(i,j,k)*ccbeta(i,j-1,k))/(ccbeta(i,j,k) &
                     + ccbeta(i,j-1,k))
             end do
          end do
       end do

       do k = lo(3),lo(3)+nz
          do j = lo(2),lo(2)+ny-1
             do i = lo(1),lo(1)+nx-1
                beta(i,j,k,3) = TWO*(ccbeta(i,j,k)*ccbeta(i,j,k-1))/(ccbeta(i,j,k) &
                     + ccbeta(i,j,k-1))
             end do
          end do
       end do

    else

       call bl_error('put_data_on_faces_2d non-harmonic not written yet')

    end if

  end subroutine put_data_on_faces_3d

end module fill_3d_module
