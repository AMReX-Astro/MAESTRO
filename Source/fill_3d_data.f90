module fill_3d_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use bl_constants_module
  use bl_prof_module

  implicit none

  private
  
  public :: put_1d_array_on_cart,  put_1d_array_on_cart_3d_sphr
  public :: make_w0mac, make_s0mac
  public :: make_normal, make_normal_3d_sphr
  public :: put_data_on_faces
  public :: put_1d_array_on_cart_irreg
  
contains  

  subroutine put_1d_array_on_cart(s0,s0_cart,bc_comp,is_input_edge_centered, &
                                  is_output_a_vector,dx,the_bc_level,mla)

    use bl_constants_module
    use define_bc_module
    use geometry, only: spherical
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
    
    integer :: lo(mla%dim),hi(mla%dim),dm,nlevs
    integer :: i,n,ng_s,comp
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_1d_array_on_cart")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_s = nghost(s0_cart(1))
    
    do n=1,nlevs
       
       do i = 1, nfabs(s0_cart(n))
          sp => dataptr(s0_cart(n), i)
          lo =  lwb(get_box(s0_cart(n), i))
          hi =  upb(get_box(s0_cart(n), i))
          select case (dm)
          case (1)
             call put_1d_array_on_cart_1d(is_input_edge_centered, &
                                          s0(n,:),sp(:,1,1,:),lo,hi,ng_s)
          case (2)
             call put_1d_array_on_cart_2d(is_input_edge_centered,is_output_a_vector, &
                                          s0(n,:),sp(:,:,1,:),lo,hi,ng_s)
          case (3)
             if (spherical .eq. 0) then
                call put_1d_array_on_cart_3d(is_input_edge_centered,is_output_a_vector, &
                                             s0(n,:),sp(:,:,:,:),lo,hi,ng_s)
             else
                call put_1d_array_on_cart_3d_sphr(is_input_edge_centered, &
                                                  is_output_a_vector,s0(1,:), &
                                                  sp(:,:,:,:),lo,hi,dx(n,:),ng_s)
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

  subroutine put_1d_array_on_cart_1d(is_input_edge_centered,s0,s0_cart, &
                                     lo,hi,ng_s)

    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    logical        , intent(in   ) :: is_input_edge_centered
    real(kind=dp_t), intent(in   ) :: s0(0:)
    real(kind=dp_t), intent(inout) :: s0_cart(lo(1)-ng_s:,:)

    integer :: i

    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_1d_array_on_cart_1d")

    s0_cart = ZERO

    if (is_input_edge_centered) then

       ! we don't do anything different in 1-d if it is a vector
       do i=lo(1),hi(1)
          s0_cart(i,1) = HALF * (s0(i) + s0(i+1))
       end do

    else
       
       do i=lo(1),hi(1)
          s0_cart(i,1) = s0(i)
       end do

    end if

    call destroy(bpt)

  end subroutine put_1d_array_on_cart_1d

  subroutine put_1d_array_on_cart_2d(is_input_edge_centered,is_output_a_vector,s0,s0_cart, &
                                     lo,hi,ng_s)

    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    logical        , intent(in   ) :: is_input_edge_centered,is_output_a_vector
    real(kind=dp_t), intent(in   ) :: s0(0:)
    real(kind=dp_t), intent(inout) :: s0_cart(lo(1)-ng_s:,lo(2)-ng_s:,:)

    integer :: j

    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_1d_array_on_cart_2d")

    s0_cart = ZERO

    if (is_input_edge_centered) then

       if (is_output_a_vector) then

          do j=lo(2),hi(2)
             s0_cart(:,j,2) = HALF * (s0(j) + s0(j+1))
          end do

       else

          do j=lo(2),hi(2)
             s0_cart(:,j,1) = HALF * (s0(j) + s0(j+1))
          end do

       end if

    else

       if (is_output_a_vector) then

          do j=lo(2),hi(2)
             s0_cart(:,j,2) = s0(j)
          end do

       else

          do j=lo(2),hi(2)
             s0_cart(:,j,1) = s0(j)
          end do

       end if

    end if

    call destroy(bpt)

  end subroutine put_1d_array_on_cart_2d

  subroutine put_1d_array_on_cart_3d(is_input_edge_centered,is_output_a_vector,s0,s0_cart, &
                                     lo,hi,ng_s)

    use bl_constants_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    logical        , intent(in   ) :: is_input_edge_centered,is_output_a_vector
    real(kind=dp_t), intent(in   ) :: s0(0:)
    real(kind=dp_t), intent(inout) :: s0_cart(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)

    integer :: k

    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_1d_array_on_cart_3d")

    s0_cart = ZERO

    if (is_input_edge_centered) then

       if (is_output_a_vector) then

          do k=lo(3),hi(3)
             s0_cart(:,:,k,3) = HALF * (s0(k) + s0(k+1))
          end do
          
       else
          
          do k=lo(3),hi(3)
             s0_cart(:,:,k,1) = HALF * (s0(k) + s0(k+1))
          end do
          
       end if

    else

       if (is_output_a_vector) then

          do k=lo(3),hi(3)
             s0_cart(:,:,k,3) = s0(k)
          end do
          
       else
          
          do k=lo(3),hi(3)
             s0_cart(:,:,k,1) = s0(k)
          end do
          
       end if

    end if

    call destroy(bpt)

  end subroutine put_1d_array_on_cart_3d

  subroutine put_1d_array_on_cart_3d_sphr(is_input_edge_centered,is_output_a_vector, &
                                          s0,s0_cart,lo,hi,dx,ng_s)

    use bl_constants_module
    use geometry, only: dr, center, r_cc_loc, nr_fine, r_edge_loc
    use probin_module, only: s0_interp_type, w0_interp_type, prob_lo

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    logical        , intent(in   ) :: is_input_edge_centered,is_output_a_vector
    real(kind=dp_t), intent(in   ) :: s0(0:)
    real(kind=dp_t), intent(inout) :: s0_cart(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer         :: i,j,k,index
    real(kind=dp_t) :: x,y,z
    real(kind=dp_t) :: radius,rfac,s0_cart_val

    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_1d_array_on_cart_3d_sphr")

    if (is_input_edge_centered) then

       ! we currently have three different ideas for computing s0_cart, 
       ! where s0 is edge-centered.
       ! 1.  Piecewise constant
       ! 2.  Piecewise linear
       ! 3.  Quadratic

       if (w0_interp_type .eq. 1) then

          !$OMP PARALLEL DO PRIVATE(i,j,k,x,y,z,radius,index,rfac,s0_cart_val)
          do k = lo(3),hi(3)
             z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
             do j = lo(2),hi(2)
                y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
                do i = lo(1),hi(1)
                   x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                   radius = sqrt(x**2 + y**2 + z**2)
                   index  = int(radius / dr(1))

                   rfac = (radius - dble(index)*dr(1)) / dr(1)

                   if (rfac .gt. 0.5d0) then
                      s0_cart_val = s0(index+1)
                   else
                      s0_cart_val = s0(index)
                   end if

                   if (is_output_a_vector) then
                      s0_cart(i,j,k,1) = s0_cart_val * x / radius
                      s0_cart(i,j,k,2) = s0_cart_val * y / radius
                      s0_cart(i,j,k,3) = s0_cart_val * z / radius
                   else
                      s0_cart(i,j,k,1) = s0_cart_val
                   end if

                end do
             end do
          end do
          !$OMP END PARALLEL DO

       else if (w0_interp_type .eq. 2) then

          !$OMP PARALLEL DO PRIVATE(i,j,k,x,y,z,radius,index,rfac,s0_cart_val)
          do k = lo(3),hi(3)
             z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
             do j = lo(2),hi(2)
                y = prob_lo(2) +(dble(j)+HALF)*dx(2) - center(2)
                do i = lo(1),hi(1)
                   x = prob_lo(1) +(dble(i)+HALF)*dx(1) - center(1)
                   radius = sqrt(x**2 + y**2 + z**2)
                   index  = int(radius / dr(1))

                   rfac = (radius - dble(index)*dr(1)) / dr(1)

                   if (index .lt. nr_fine) then
                      s0_cart_val = rfac * s0(index+1) + (ONE-rfac) * s0(index)
                   else
                      s0_cart_val = s0(nr_fine)
                   end if

                   if (is_output_a_vector) then
                      s0_cart(i,j,k,1) = s0_cart_val * x / radius
                      s0_cart(i,j,k,2) = s0_cart_val * y / radius
                      s0_cart(i,j,k,3) = s0_cart_val * z / radius
                   else
                      s0_cart(i,j,k,1) = s0_cart_val
                   end if

                end do
             end do
          end do
          !$OMP END PARALLEL DO

       else if (w0_interp_type .eq. 3) then

          !$OMP PARALLEL DO PRIVATE(i,j,k,x,y,z,radius,index,s0_cart_val)
          do k = lo(3),hi(3)
             z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
             do j = lo(2),hi(2)
                y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
                do i = lo(1),hi(1)
                   x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                   radius = sqrt(x**2 + y**2 + z**2)
                   index  = int(radius / dr(1))

                   ! index refers to the lo point in the quadratic stencil
                   if (index .le. 0) then
                      index = 0
                   else if (index .ge. nr_fine-1) then
                      index = nr_fine-2
                   else if (radius-r_edge_loc(1,index) .lt. r_edge_loc(1,index+1)) then
                      index = index-1
                   end if

                   call quad_interp(radius, &
                                    r_edge_loc(1,index),r_edge_loc(1,index+1), &
                                    r_edge_loc(1,index+2), &
                                    s0_cart_val, &
                                    s0(index),s0(index+1),s0(index+2))

                   if (is_output_a_vector) then
                      s0_cart(i,j,k,1) = s0_cart_val * x / radius
                      s0_cart(i,j,k,2) = s0_cart_val * y / radius
                      s0_cart(i,j,k,3) = s0_cart_val * z / radius
                   else
                      s0_cart(i,j,k,1) = s0_cart_val
                   end if

                end do
             end do
          end do
          !$OMP END PARALLEL DO

       else
          call bl_error('Error: w0_interp_type not defined')
       end if

    else

       ! we currently have three different ideas for computing s0_cart, 
       ! where s0 is bin-centered.
       ! 1.  Piecewise constant
       ! 2.  Piecewise linear
       ! 3.  Quadratic
       
       if (s0_interp_type .eq. 1) then

          !$OMP PARALLEL DO PRIVATE(i,j,k,x,y,z,radius,index,s0_cart_val)
          do k = lo(3),hi(3)
             z = prob_lo(3) +(dble(k)+HALF)*dx(3) - center(3)
             do j = lo(2),hi(2)
                y = prob_lo(2) +(dble(j)+HALF)*dx(2) - center(2)
                do i = lo(1),hi(1)
                   x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                   radius = sqrt(x**2 + y**2 + z**2)
                   index  = int(radius / dr(1))

                   s0_cart_val = s0(index)

                   if (is_output_a_vector) then
                      s0_cart(i,j,k,1) = s0_cart_val * x / radius
                      s0_cart(i,j,k,2) = s0_cart_val * y / radius
                      s0_cart(i,j,k,3) = s0_cart_val * z / radius
                   else
                      s0_cart(i,j,k,1) = s0_cart_val
                   end if

                end do
             end do
          end do
          !$OMP END PARALLEL DO

       else if (s0_interp_type .eq. 2) then

          !$OMP PARALLEL DO PRIVATE(i,j,k,x,y,z,radius,index,s0_cart_val)
          do k = lo(3),hi(3)
             z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
             do j = lo(2),hi(2)
                y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
                do i = lo(1),hi(1)
                   x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                   radius = sqrt(x**2 + y**2 + z**2)
                   index  = int(radius / dr(1))

                   if (radius .ge. r_cc_loc(1,index)) then
                      if (index .ge. nr_fine-1) then
                         s0_cart_val = s0(nr_fine-1)
                      else
                         s0_cart_val = s0(index+1)*(radius-r_cc_loc(1,index))/dr(1) &
                              + s0(index)*(r_cc_loc(1,index+1)-radius)/dr(1)
                      endif
                   else
                      if (index .eq. 0) then
                         s0_cart_val = s0(index)
                      else if (index .gt. nr_fine-1) then
                         s0_cart_val = s0(nr_fine-1)
                      else
                         s0_cart_val = s0(index)*(radius-r_cc_loc(1,index-1))/dr(1) &
                              + s0(index-1)*(r_cc_loc(1,index)-radius)/dr(1)
                      end if
                   end if

                   if (is_output_a_vector) then
                      s0_cart(i,j,k,1) = s0_cart_val * x / radius
                      s0_cart(i,j,k,2) = s0_cart_val * y / radius
                      s0_cart(i,j,k,3) = s0_cart_val * z / radius
                   else
                      s0_cart(i,j,k,1) = s0_cart_val
                   end if

                end do
             end do
          end do
          !$OMP END PARALLEL DO

       else if (s0_interp_type .eq. 3) then

          !$OMP PARALLEL DO PRIVATE(i,j,k,x,y,z,radius,index,s0_cart_val)
          do k = lo(3),hi(3)
             z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
             do j = lo(2),hi(2)
                y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
                do i = lo(1),hi(1)
                   x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                   radius = sqrt(x**2 + y**2 + z**2)
                   index  = int(radius / dr(1))

                   ! index refers to the center point in the quadratic stencil.
                   ! we need to modify this if we're too close to the edge
                   if (index .eq. 0) then
                      index = 1
                   else if (index .ge. nr_fine-1) then
                      index = nr_fine-2
                   end if

                   call quad_interp(radius, &
                                    r_cc_loc(1,index-1),r_cc_loc(1,index), &
                                    r_cc_loc(1,index+1), &
                                    s0_cart_val, &
                                    s0(index-1),s0(index),s0(index+1))

                   if (is_output_a_vector) then
                      s0_cart(i,j,k,1) = s0_cart_val * x / radius
                      s0_cart(i,j,k,2) = s0_cart_val * y / radius
                      s0_cart(i,j,k,3) = s0_cart_val * z / radius
                   else
                      s0_cart(i,j,k,1) = s0_cart_val
                   end if

                end do
             end do
          end do
          !$OMP END PARALLEL DO

       else
          call bl_error('Error: s0_interp_type not defined')
       end if

    end if

    call destroy(bpt)

  end subroutine put_1d_array_on_cart_3d_sphr

  subroutine quad_interp(x,x0,x1,x2,y,y0,y1,y2)

    real(kind=dp_t), intent(in   ) :: x,x0,x1,x2,y0,y1,y2
    real(kind=dp_t), intent(  out) :: y
    
    y = y0 + (y1-y0)/(x1-x0)*(x-x0) &
           + ((y2-y1)/(x2-x1)-(y1-y0)/(x1-x0))/(x2-x0)*(x-x0)*(x-x1)

    if (y .gt. max(y0,y1,y2)) y = max(y0,y1,y2)
    if (y .lt. min(y0,y1,y2)) y = min(y0,y1,y2)

  end subroutine quad_interp

  subroutine make_w0mac(mla,w0,w0mac,dx,the_bc_level)

    use bl_constants_module
    use geometry, only: spherical
    use probin_module, only: w0mac_interp_type
    use define_bc_module

    type(ml_layout), intent(in   ) :: mla
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(inout) :: w0mac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local variables
    integer         :: lo(mla%dim),hi(mla%dim)
    integer         :: i,n,ng_w0,ng_wc,dm,nlevs

    ! Local pointers
    real(kind=dp_t), pointer :: w0xp(:,:,:,:)
    real(kind=dp_t), pointer :: w0yp(:,:,:,:)
    real(kind=dp_t), pointer :: w0zp(:,:,:,:)
    real(kind=dp_t), pointer :: w0p(:,:,:,:)

    type(multifab) :: w0_cart(mla%nlevel)
    
    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_w0mac")

    if (spherical .eq. 0) then
       call bl_error('Error: only call make_w0mac for spherical')
    end if

    dm = mla%dim
    nlevs = mla%nlevel

    ! we first need to construct a cart version of w0
    do n=1,nlevs
       call build(w0_cart(n),mla%la(n),dm,2)
    end do

    if (w0mac_interp_type .eq. 1) then
       call put_1d_array_on_cart(w0,w0_cart,1,.true.,.true.,dx,the_bc_level,mla)
    end if

    ng_w0 = nghost(w0mac(1,1))
    ng_wc = nghost(w0_cart(1))

    if (ng_w0 .ne. 1) then
       call bl_error('Error: make_w0mac_3d_sphr assumes one ghost cell')
    end if
    
    do n=1,nlevs
       do i=1, nfabs(w0mac(n,1))
          w0xp => dataptr(w0mac(n,1), i)
          w0yp => dataptr(w0mac(n,2), i)
          w0zp => dataptr(w0mac(n,3), i)
          w0p  => dataptr(w0_cart(n), i)
          lo = lwb(get_box(w0mac(n,1), i))
          hi = upb(get_box(w0mac(n,1), i))
          call make_w0mac_3d_sphr(w0(1,:),w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1), &
                                  ng_w0,w0p(:,:,:,:),ng_wc,lo,hi,dx(n,:))
       end do
    end do

    do n=1,nlevs
       call destroy(w0_cart(n))
    end do

    call destroy(bpt)

  end subroutine make_w0mac

  subroutine make_w0mac_3d_sphr(w0,w0macx,w0macy,w0macz,ng_w0,w0_cart,ng_wc,lo,hi,dx)

    use bl_constants_module
    use geometry, only: dr, center, nr_fine, r_edge_loc
    use probin_module, only: w0mac_interp_type, prob_lo

    integer        , intent(in   ) :: lo(:),hi(:),ng_w0,ng_wc
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(inout) ::  w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(inout) ::  w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(inout) ::  w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(inout) :: w0_cart(lo(1)-ng_wc:,lo(2)-ng_wc:,lo(3)-ng_wc:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer         :: i,j,k,index
    real(kind=dp_t) :: x,y,z
    real(kind=dp_t) :: radius,w0_cart_val,rfac
    real(kind=dp_t), allocatable :: w0_nodal(:,:,:,:)

    ! we currently have three different ideas for computing w0mac
    ! 1.  Interpolate w0 to cell centers, then average to edges
    ! 2.  Interpolate w0 to edges directly using linear interpolation
    ! 3.  Interpolate w0 to edges directly using quadratic interpolation
    ! 4.  Interpolate w0 to nodes, then average to edges

    if (w0mac_interp_type .eq. 1) then

       !$OMP PARALLEL PRIVATE(i,j,k)

       !$OMP DO
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+2
                w0macx(i,j,k) = HALF* (w0_cart(i-1,j,k,1) + w0_cart(i,j,k,1))
             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k=lo(3)-1,hi(3)+1
          do j=lo(2)-1,hi(2)+2
             do i=lo(1)-1,hi(1)+1
                w0macy(i,j,k) = HALF* (w0_cart(i,j-1,k,2) + w0_cart(i,j,k,2))
             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k=lo(3)-1,hi(3)+2
          do j=lo(2)-1,hi(2)+1
             do i=lo(1)-1,hi(1)+1
                w0macz(i,j,k) = HALF* (w0_cart(i,j,k-1,3) + w0_cart(i,j,k,3))
             end do
          end do
       end do
       !$OMP END DO

       !$OMP END PARALLEL

    else if (w0mac_interp_type .eq. 2) then

       !$OMP PARALLEL PRIVATE(i,j,k,x,y,z,radius,index,rfac,w0_cart_val)

       !$OMP DO
       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+2
                x = prob_lo(1) + (dble(i)     )*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(1))

                rfac = (radius - dble(index)*dr(1)) / dr(1)

                if (index .lt. nr_fine) then
                   w0_cart_val = rfac * w0(index+1) + (ONE-rfac) * w0(index)
                else
                   w0_cart_val = w0(nr_fine)
                end if

                w0macx(i,j,k) = w0_cart_val * x / radius

             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+2
             y = prob_lo(2) + (dble(j)     )*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(1))

                rfac = (radius - dble(index)*dr(1)) / dr(1)

                if (index .lt. nr_fine) then
                   w0_cart_val = rfac * w0(index+1) + (ONE-rfac) * w0(index)
                else
                   w0_cart_val = w0(nr_fine)
                end if

                w0macy(i,j,k) = w0_cart_val * y / radius

             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k = lo(3)-1,hi(3)+2
          z = prob_lo(3) + (dble(k)     )*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(1))

                rfac = (radius - dble(index)*dr(1)) / dr(1)

                if (index .lt. nr_fine) then
                   w0_cart_val = rfac * w0(index+1) + (ONE-rfac) * w0(index)
                else
                   w0_cart_val = w0(nr_fine)
                end if

                w0macz(i,j,k) = w0_cart_val * z / radius

             end do
          end do
       end do
       !$OMP END DO

       !$OMP END PARALLEL

    else if (w0mac_interp_type .eq. 3) then

       !$OMP PARALLEL PRIVATE(i,j,k,x,y,z,radius,index,w0_cart_val)

       !$OMP DO
       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+2
                x = prob_lo(1) + (dble(i)     )*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(1))

                ! index refers to the lo point in the quadratic stencil
                if (index .le. 0) then
                   index = 0
                else if (index .ge. nr_fine-1) then
                   index = nr_fine-2
                else if (radius-r_edge_loc(1,index) .lt. r_edge_loc(1,index+1)) then
                   index = index-1
                end if

                call quad_interp(radius, &
                                 r_edge_loc(1,index),r_edge_loc(1,index+1), &
                                 r_edge_loc(1,index+2), &
                                 w0_cart_val, &
                                 w0(index),w0(index+1),w0(index+2))

                w0macx(i,j,k) = w0_cart_val * x / radius

             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+2
             y = prob_lo(2) + (dble(j)     )*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(1))

                ! index refers to the lo point in the quadratic stencil
                if (index .le. 0) then
                   index = 0
                else if (index .ge. nr_fine-1) then
                   index = nr_fine-2
                else if (radius-r_edge_loc(1,index) .lt. r_edge_loc(1,index+1)) then
                   index = index-1
                end if

                call quad_interp(radius, &
                                 r_edge_loc(1,index),r_edge_loc(1,index+1), &
                                 r_edge_loc(1,index+2), &
                                 w0_cart_val, &
                                 w0(index),w0(index+1),w0(index+2))

                w0macy(i,j,k) = w0_cart_val * y / radius

             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k = lo(3)-1,hi(3)+2
          z = prob_lo(3) + (dble(k)     )*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(1))

                ! index refers to the lo point in the quadratic stencil
                if (index .le. 0) then
                   index = 0
                else if (index .ge. nr_fine-1) then
                   index = nr_fine-2
                else if (radius-r_edge_loc(1,index) .lt. r_edge_loc(1,index+1)) then
                   index = index-1
                end if

                call quad_interp(radius, &
                                 r_edge_loc(1,index),r_edge_loc(1,index+1), &
                                 r_edge_loc(1,index+2), &
                                 w0_cart_val, &
                                 w0(index),w0(index+1),w0(index+2))

                w0macz(i,j,k) = w0_cart_val * z / radius

             end do
          end do
       end do
       !$OMP END DO

       !$OMP END PARALLEL

    else if (w0mac_interp_type .eq. 4) then

       allocate(w0_nodal(lo(1)-1:hi(1)+2,lo(2)-1:hi(2)+2,lo(3)-1:hi(3)+2,3))

       !$OMP PARALLEL DO PRIVATE(i,j,k,x,y,z,radius,index,rfac,w0_cart_val)
       do k = lo(3)-1,hi(3)+2
          z = prob_lo(3) + (dble(k))*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+2
             y = prob_lo(2) + (dble(j))*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+2
                x = prob_lo(1) + (dble(i))*dx(1) - center(1)

                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(1))
                
                rfac = (radius - dble(index)*dr(1)) / dr(1)

                if (index .lt. nr_fine) then
                   w0_cart_val = rfac * w0(index+1) + (ONE-rfac) * w0(index)
                else
                   w0_cart_val = w0(nr_fine)
                end if

                w0_nodal(i,j,k,1) = w0_cart_val * x / radius
                w0_nodal(i,j,k,2) = w0_cart_val * y / radius
                w0_nodal(i,j,k,3) = w0_cart_val * z / radius

             end do
          end do
       end do
       !$OMP END PARALLEL DO

       !$OMP PARALLEL PRIVATE(i,j,k)

       !$OMP DO
       do k = lo(3)-1,hi(3)+1
          do j = lo(2)-1,hi(2)+1
             do i = lo(1)-1,hi(1)+2
                w0macx(i,j,k) = FOURTH*( w0_nodal(i,j,k  ,1) + w0_nodal(i,j+1,k  ,1) &
                                        +w0_nodal(i,j,k+1,1) + w0_nodal(i,j+1,k+1,1))
             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k = lo(3)-1,hi(3)+1
          do j = lo(2)-1,hi(2)+2
             do i = lo(1)-1,hi(1)+1
                w0macy(i,j,k) = FOURTH*( w0_nodal(i,j,k  ,2) + w0_nodal(i+1,j,k  ,2) &
                                        +w0_nodal(i,j,k+1,2) + w0_nodal(i+1,j,k+1,2))
             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k = lo(3)-1,hi(3)+2
          do j = lo(2)-1,hi(2)+1
             do i = lo(1)-1,hi(1)+1
                w0macz(i,j,k) = FOURTH*( w0_nodal(i,j  ,k,3) + w0_nodal(i+1,j  ,k,3) &
                                        +w0_nodal(i,j+1,k,3) + w0_nodal(i+1,j+1,k,3))
             end do
          end do
       end do
       !$OMP END DO

       !$OMP END PARALLEL

       deallocate(w0_nodal)

    else
       call bl_error('Error: w0mac_interp_type not defined')
    end if

  end subroutine make_w0mac_3d_sphr

  subroutine make_s0mac(mla,s0,s0mac,dx,bccomp,the_bc_level)

    use bl_constants_module
    use geometry, only: spherical
    use define_bc_module
    use probin_module, only: s0mac_interp_type

    type(ml_layout), intent(in   ) :: mla
    real(kind=dp_t), intent(in   ) :: s0(:,0:)
    type(multifab) , intent(inout) :: s0mac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: bccomp
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! Local variables
    integer         :: lo(mla%dim),hi(mla%dim)
    integer         :: i,n,ng_sm,ng_s0,dm,nlevs

    ! Local pointers
    real(kind=dp_t), pointer :: s0xp(:,:,:,:)
    real(kind=dp_t), pointer :: s0yp(:,:,:,:)
    real(kind=dp_t), pointer :: s0zp(:,:,:,:)
    real(kind=dp_t), pointer :: s0p(:,:,:,:)

    type(multifab) :: s0_cart(mla%nlevel)
    
    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_s0mac")

    dm = mla%dim
    nlevs = mla%nlevel

    if (spherical .eq. 0) then
       call bl_error('Error: only call make_s0mac for spherical')
    end if

    ! we first need to construct a cart version of s0
    do n=1,nlevs
       call build(s0_cart(n),mla%la(n),1,2)
    end do

    if (s0mac_interp_type .eq. 1) then
       call put_1d_array_on_cart(s0,s0_cart,bccomp,.false.,.false.,dx,the_bc_level,mla)
    end if

    ng_sm = nghost(s0mac(1,1))
    ng_s0 = nghost(s0_cart(1))

    if (ng_sm .ne. 1) then
       call bl_error('Error: make_s0mac assumes one ghost cell in s0mac')
    end if

    do n=1,nlevs
       do i=1, nfabs(s0mac(n,1))
          s0xp => dataptr(s0mac(n,1), i)
          s0yp => dataptr(s0mac(n,2), i)
          s0zp => dataptr(s0mac(n,3), i)
          s0p  => dataptr(s0_cart(n), i)
          lo = lwb(get_box(s0mac(n,1), i))
          hi = upb(get_box(s0mac(n,1), i))
          call make_s0mac_3d_sphr(s0(1,:),s0xp(:,:,:,1),s0yp(:,:,:,1), &
                                  s0zp(:,:,:,1),ng_sm,s0p(:,:,:,1),ng_s0, &
                                  lo,hi,dx(n,:))
       end do
    end do

    do n=1,nlevs
       call destroy(s0_cart(n))
    end do

    call destroy(bpt)

  end subroutine make_s0mac

  subroutine make_s0mac_3d_sphr(s0,s0macx,s0macy,s0macz,ng_sm,s0_cart,ng_s0,lo,hi,dx)

    use bl_constants_module
    use geometry, only: dr, center, nr_fine, r_cc_loc
    use probin_module, only: s0mac_interp_type, prob_lo

    integer        , intent(in   ) :: lo(:),hi(:),ng_sm,ng_s0
    real(kind=dp_t), intent(in   ) :: s0(0:)
    real(kind=dp_t), intent(inout) ::  s0macx(lo(1)-ng_sm:,lo(2)-ng_sm:,lo(3)-ng_sm:)
    real(kind=dp_t), intent(inout) ::  s0macy(lo(1)-ng_sm:,lo(2)-ng_sm:,lo(3)-ng_sm:)
    real(kind=dp_t), intent(inout) ::  s0macz(lo(1)-ng_sm:,lo(2)-ng_sm:,lo(3)-ng_sm:)
    real(kind=dp_t), intent(inout) :: s0_cart(lo(1)-ng_s0:,lo(2)-ng_s0:,lo(3)-ng_s0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer         :: i,j,k,index
    real(kind=dp_t) :: x,y,z
    real(kind=dp_t) :: radius

    ! we currently have three different ideas for computing s0mac
    ! 1.  Interpolate s0 to cell centers, then average to edges
    ! 2.  Interpolate s0 to edges directly using linear interpolation
    ! 3.  Interpolate s0 to edges directly using quadratic interpolation
    ! 4.  Interpolate s0 to nodes, then average to edges

    if (s0mac_interp_type .eq. 1) then

       !$OMP PARALLEL PRIVATE(i,j,k)

       !$OMP DO
       do k = lo(3)-1,hi(3)+1
          do j = lo(2)-1,hi(2)+1
             do i = lo(1)-1,hi(1)+2
                s0macx(i,j,k) = HALF*(s0_cart(i,j,k)+s0_cart(i-1,j,k))
             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k = lo(3)-1,hi(3)+1
          do j = lo(2)-1,hi(2)+2
             do i = lo(1)-1,hi(1)+1
                s0macy(i,j,k) = HALF*(s0_cart(i,j,k)+s0_cart(i,j-1,k))
             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k = lo(3)-1,hi(3)+2
          do j = lo(2)-1,hi(2)+1
             do i = lo(1)-1,hi(1)+1
                s0macz(i,j,k) = HALF*(s0_cart(i,j,k)+s0_cart(i,j,k-1))
             end do
          end do
       end do
       !$OMP END DO

       !$OMP END PARALLEL
      
    else if (s0mac_interp_type .eq. 2) then

       !$OMP PARALLEL PRIVATE(i,j,k,x,y,z,radius,index)

       !$OMP DO
       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+2
                x = prob_lo(1) + (dble(i)     )*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(1))

                if (radius .ge. r_cc_loc(1,index)) then
                   if (index .ge. nr_fine-1) then
                      s0macx(i,j,k) = s0(nr_fine-1)
                   else
                      s0macx(i,j,k) = s0(index+1)*(radius-r_cc_loc(1,index))/dr(1) &
                           + s0(index)*(r_cc_loc(1,index+1)-radius)/dr(1)
                   endif
                else
                   if (index .eq. 0) then
                      s0macx(i,j,k) = s0(index)
                   else if (index .gt. nr_fine-1) then
                      s0macx(i,j,k) = s0(nr_fine-1)
                   else
                      s0macx(i,j,k) = s0(index)*(radius-r_cc_loc(1,index-1))/dr(1) &
                           + s0(index-1)*(r_cc_loc(1,index)-radius)/dr(1)
                   end if
                end if

             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+2
             y = prob_lo(2) + (dble(j)     )*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(1))

                if (radius .ge. r_cc_loc(1,index)) then
                   if (index .ge. nr_fine-1) then
                      s0macy(i,j,k) = s0(nr_fine-1)
                   else
                      s0macy(i,j,k) = s0(index+1)*(radius-r_cc_loc(1,index))/dr(1) &
                           + s0(index)*(r_cc_loc(1,index+1)-radius)/dr(1)
                   endif
                else
                   if (index .eq. 0) then
                      s0macy(i,j,k) = s0(index)
                   else if (index .gt. nr_fine-1) then
                      s0macy(i,j,k) = s0(nr_fine-1)
                   else
                      s0macy(i,j,k) = s0(index)*(radius-r_cc_loc(1,index-1))/dr(1) &
                           + s0(index-1)*(r_cc_loc(1,index)-radius)/dr(1)
                   end if
                end if

             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k = lo(3)-1,hi(3)+2
          z = prob_lo(3) + (dble(k)     )*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(1))

                if (radius .ge. r_cc_loc(1,index)) then
                   if (index .ge. nr_fine-1) then
                      s0macz(i,j,k) = s0(nr_fine-1)
                   else
                      s0macz(i,j,k) = s0(index+1)*(radius-r_cc_loc(1,index))/dr(1) &
                           + s0(index)*(r_cc_loc(1,index+1)-radius)/dr(1)
                   endif
                else
                   if (index .eq. 0) then
                      s0macz(i,j,k) = s0(index)
                   else if (index .gt. nr_fine-1) then
                      s0macz(i,j,k) = s0(nr_fine-1)
                   else
                      s0macz(i,j,k) = s0(index)*(radius-r_cc_loc(1,index-1))/dr(1) &
                           + s0(index-1)*(r_cc_loc(1,index)-radius)/dr(1)
                   end if
                end if

             end do
          end do
       end do
       !$OMP END DO

       !$OMP END PARALLEL

    else if (s0mac_interp_type .eq. 3) then

       !$OMP PARALLEL PRIVATE(i,j,k,x,y,z,radius,index)

       !$OMP DO
       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+2
                x = prob_lo(1) + (dble(i)     )*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(1))

                ! index refers to the center point in the quadratic stencil.
                ! we need to modify this if we're too close to the edge
                if (index .eq. 0) then
                   index = 1
                else if (index .ge. nr_fine-1) then
                   index = nr_fine-2
                end if

                call quad_interp(radius, &
                                 r_cc_loc(1,index-1),r_cc_loc(1,index), &
                                 r_cc_loc(1,index+1), &
                                 s0macx(i,j,k), &
                                 s0(index-1),s0(index),s0(index+1))
             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k = lo(3)-1,hi(3)+1
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+2
             y = prob_lo(2) + (dble(j)     )*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(1))

                ! index refers to the center point in the quadratic stencil.
                ! we need to modify this if we're too close to the edge
                if (index .eq. 0) then
                   index = 1
                else if (index .ge. nr_fine-1) then
                   index = nr_fine-2
                end if

                call quad_interp(radius, &
                                 r_cc_loc(1,index-1),r_cc_loc(1,index), &
                                 r_cc_loc(1,index+1), &
                                 s0macy(i,j,k), &
                                 s0(index-1),s0(index),s0(index+1))
             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k = lo(3)-1,hi(3)+2
          z = prob_lo(3) + (dble(k)     )*dx(3) - center(3)
          do j = lo(2)-1,hi(2)+1
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-1,hi(1)+1
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
                radius = sqrt(x**2 + y**2 + z**2)
                index  = int(radius / dr(1))

                ! index refers to the center point in the quadratic stencil.
                ! we need to modify this if we're too close to the edge
                if (index .eq. 0) then
                   index = 1
                else if (index .ge. nr_fine-1) then
                   index = nr_fine-2
                end if

                call quad_interp(radius, &
                                 r_cc_loc(1,index-1),r_cc_loc(1,index), &
                                 r_cc_loc(1,index+1), &
                                 s0macz(i,j,k), &
                                 s0(index-1),s0(index),s0(index+1))
             end do
          end do
       end do
       !$OMP END DO

       !$OMP END PARALLEL

    else

       call bl_error('Error: s0mac_interp_type not defined')

    end if

  end subroutine make_s0mac_3d_sphr

  subroutine make_normal(normal,dx)

    use geometry, only: spherical

    type(multifab) , intent(inout) :: normal(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
        
    integer             :: lo(get_dim(normal(1))),hi(get_dim(normal(1)))
    integer             :: n,i,ng_n,dm,nlevs
    real(dp_t), pointer :: nop(:,:,:,:)
        
    ng_n = nghost(normal(1))
    dm = get_dim(normal(1))
    nlevs = size(normal)

    if (spherical .eq. 1) then
       do n = 1,nlevs
          do i = 1, nfabs(normal(n))
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
    use probin_module, only: prob_lo
    
    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(  out) :: normal(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)

    integer         :: i,j,k
    real(kind=dp_t) :: x,y,z,radius

    ! normal is the unit vector in the radial direction (e_r) in spherical
    ! coordinates.
    !
    ! in terms of Cartesian coordinates, with unit vectors e_x, e_y, e_z,
    !    e_r = sin(theta)cos(phi) e_x + sin(theta)sin(phi) e_y + cos(theta) e_z
    ! or
    !    e_r = (x/R) e_x + (y/R) e_y + (z/R) e_z

    if (spherical .eq. 1) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,x,y,z,radius)
       do k = lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
          do j = lo(2)-ng,hi(2)+ng
             y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
             do i = lo(1)-ng,hi(1)+ng
                x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)

                radius = sqrt(x**2 + y**2 + z**2)

                normal(i,j,k,1) = x / radius
                normal(i,j,k,2) = y / radius
                normal(i,j,k,3) = z / radius

             end do
          end do
       end do
       !$OMP END PARALLEL DO

    else 
       call bl_error('SHOULDNT CALL MAKE_3D_NORMAL WITH SPHERICAL = 0')
    end if

  end subroutine make_normal_3d_sphr

  subroutine put_data_on_faces(mla,ccfab,comp,beta,harmonic_avg)

    use ml_restriction_module, only: ml_edge_restriction

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: ccfab(:)
    type(multifab) , intent(inout) :: beta(:,:)
    integer        , intent(in   ) :: comp
    logical        , intent(in   ) :: harmonic_avg

    ! local
    integer :: n,i,ng_c,ng_b,dm,nlevs
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: ccfabp(:,:,:,:)
    real(kind=dp_t), pointer :: bxp(:,:,:,:)
    real(kind=dp_t), pointer :: byp(:,:,:,:)
    real(kind=dp_t), pointer :: bzp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_data_on_faces")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_c = nghost(ccfab(1))
    ng_b = nghost(beta(1,1))

    ! setup beta = ccfab on faces
    do n=1,nlevs
       do i=1, nfabs(ccfab(n))
          ccfabp => dataptr(ccfab(n),i)
          bxp   => dataptr(beta(n,1),i)
          lo = lwb(get_box(ccfab(n),i))
          hi = upb(get_box(ccfab(n),i))
          select case (dm)
          case (1)
             call put_data_on_faces_1d(lo,hi,ccfabp(:,1,1,comp),ng_c, &
                                       bxp(:,1,1,1),ng_b, &
                                       harmonic_avg)
          case (2)
             byp   => dataptr(beta(n,2),i)
             call put_data_on_faces_2d(lo,hi,ccfabp(:,:,1,comp),ng_c, &
                                       bxp(:,:,1,1),byp(:,:,1,1),ng_b, &
                                       harmonic_avg)
          case (3)
             byp   => dataptr(beta(n,2),i)
             bzp   => dataptr(beta(n,3),i)
             call put_data_on_faces_3d(lo,hi,ccfabp(:,:,:,comp),ng_c, &
                                       bxp(:,:,:,1),byp(:,:,:,1),bzp(:,:,:,1),ng_b, &
                                       harmonic_avg)
          end select
       end do
    enddo

    ! Make sure that the fine edges average down onto the coarse edges.
    do n = nlevs,2,-1
       do i = 1,dm
          call ml_edge_restriction(beta(n-1,i),beta(n,i),mla%mba%rr(n-1,:),i)
       end do
    end do

    call destroy(bpt)

  end subroutine put_data_on_faces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! put beta on faces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine put_data_on_faces_1d(lo,hi,ccbeta,ng_c,betax,ng_b,harmonic_avg)

    integer        , intent(in   ) :: lo(:), hi(:), ng_c, ng_b
    real(kind=dp_t), intent(in   ) :: ccbeta(lo(1)-ng_c:)
    real(kind=dp_t), intent(inout) ::  betax(lo(1)-ng_b:)
    logical        , intent(in   ) :: harmonic_avg

    ! Local
    integer         :: i
    real(kind=dp_t) :: denom

    if (harmonic_avg) then

       do i = lo(1),hi(1)+1
          denom = (ccbeta(i) + ccbeta(i-1))
          if (denom .ne. 0.d0) then
             betax(i) = TWO*(ccbeta(i)*ccbeta(i-1)) / denom
          else
             betax(i) = HALF*(ccbeta(i)+ccbeta(i-1))
          end if
       end do

    else

       do i = lo(1),hi(1)+1
          betax(i) = HALF*(ccbeta(i)+ccbeta(i-1))
       end do

    end if

  end subroutine put_data_on_faces_1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! put beta on faces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine put_data_on_faces_2d(lo,hi,ccbeta,ng_c,betax,betay,ng_b,harmonic_avg)

    integer        , intent(in   ) :: lo(:), hi(:), ng_c, ng_b
    real(kind=dp_t), intent(in   ) :: ccbeta(lo(1)-ng_c:,lo(2)-ng_c:)
    real(kind=dp_t), intent(inout) ::  betax(lo(1)-ng_b:,lo(2)-ng_b:)
    real(kind=dp_t), intent(inout) ::  betay(lo(1)-ng_b:,lo(2)-ng_b:)
    logical        , intent(in   ) :: harmonic_avg

    ! Local
    integer         :: i,j
    real(kind=dp_t) :: denom

    if (harmonic_avg) then

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             denom = (ccbeta(i,j) + ccbeta(i-1,j))
             if (denom .ne. 0.d0) then
                betax(i,j) = TWO*(ccbeta(i,j)*ccbeta(i-1,j)) / denom
             else
                betax(i,j) = HALF*(ccbeta(i,j)+ccbeta(i-1,j))
             end if
          end do
       end do

       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
             denom = (ccbeta(i,j) + ccbeta(i,j-1))
             if (denom .ne. 0.d0) then
                betay(i,j) = TWO*(ccbeta(i,j)*ccbeta(i,j-1)) / denom
             else
                betay(i,j) = HALF*(ccbeta(i,j)+ccbeta(i,j-1))
             end if
          end do
       end do

    else

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1
             betax(i,j) = HALF*(ccbeta(i,j)+ccbeta(i-1,j))
          end do
       end do

       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)
             betay(i,j) = HALF*(ccbeta(i,j)+ccbeta(i,j-1))
          end do
       end do

    end if

  end subroutine put_data_on_faces_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! put beta on faces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine put_data_on_faces_3d(lo,hi,ccbeta,ng_c,betax,betay,betaz,ng_b,harmonic_avg)

    integer        , intent(in   ) :: lo(:), hi(:), ng_c, ng_b
    real(kind=dp_t), intent(in   ) :: ccbeta(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
    real(kind=dp_t), intent(inout) ::  betax(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real(kind=dp_t), intent(inout) ::  betay(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real(kind=dp_t), intent(inout) ::  betaz(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    logical        , intent(in   ) :: harmonic_avg

    ! Local
    integer         :: i,j,k
    real(kind=dp_t) :: denom

    if (harmonic_avg) then

       !$OMP PARALLEL PRIVATE(i,j,k,denom)

       !$OMP DO
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)+1
                denom = (ccbeta(i,j,k) + ccbeta(i-1,j,k))
                if (denom .ne. 0.d0) then
                   betax(i,j,k) = TWO*(ccbeta(i,j,k) * ccbeta(i-1,j,k)) / denom
                else
                   betax(i,j,k) = HALF*denom
                end if
             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)
                denom = (ccbeta(i,j,k) + ccbeta(i,j-1,k))
                if (denom .ne. 0.d0) then
                   betay(i,j,k) = TWO*(ccbeta(i,j,k) * ccbeta(i,j-1,k)) / denom
                else
                   betay(i,j,k) = HALF*denom
                end if
             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k = lo(3),hi(3)+1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                denom = (ccbeta(i,j,k) + ccbeta(i,j,k-1))
                if (denom .ne. 0.d0) then
                   betaz(i,j,k) = TWO*(ccbeta(i,j,k) * ccbeta(i,j,k-1)) / denom
                else
                   betaz(i,j,k) = HALF*denom
                end if
             end do
          end do
       end do
       !$OMP END DO

       !$OMP END PARALLEL

    else

       !$OMP PARALLEL PRIVATE(i,j,k)

       !$OMP DO
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)+1
                betax(i,j,k) = HALF*(ccbeta(i,j,k)+ccbeta(i-1,j,k))
             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)+1
             do i = lo(1),hi(1)
                betay(i,j,k) = HALF*(ccbeta(i,j,k)+ccbeta(i,j-1,k))
             end do
          end do
       end do
       !$OMP END DO NOWAIT

       !$OMP DO
       do k = lo(3),hi(3)+1
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                betaz(i,j,k) = HALF*(ccbeta(i,j,k)+ccbeta(i,j,k-1))
             end do
          end do
       end do
       !$OMP END DO

       !$OMP END PARALLEL

    end if

  end subroutine put_data_on_faces_3d

  subroutine put_1d_array_on_cart_irreg(s0,s0_cart,bc_comp,dx,the_bc_level,mla)

    use bl_constants_module
    use define_bc_module
    use geometry, only: spherical, nr_irreg
    use ml_layout_module
    use multifab_physbc_module
    use ml_restriction_module, only: ml_cc_restriction_c
    use multifab_fill_ghost_module
    
    real(kind=dp_t), intent(in   ) :: s0(:,0:)
    type(multifab) , intent(inout) :: s0_cart(:)
    integer        , intent(in   ) :: bc_comp
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(in   ) :: mla
    
    integer :: lo(mla%dim),hi(mla%dim)
    integer :: i,n,r,ng_s,dm,nlevs
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    real(kind=dp_t), allocatable :: radii(:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_1d_array_on_cart")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_s = nghost(s0_cart(1))
    
    ! radii contains every possible distance that a cell-center at the finest
    ! level can map into
    allocate(radii(0:nr_irreg))

    !$OMP PARALLEL DO PRIVATE(r)
    do r=0,nr_irreg
       radii(r) = sqrt(0.75d0+2.d0*r)*dx(nlevs,1)
    end do
    !$OMP END PARALLEL DO

    do n=1,nlevs
       
       do i = 1, nfabs(s0_cart(n))
          sp => dataptr(s0_cart(n), i)
          lo =  lwb(get_box(s0_cart(n), i))
          hi =  upb(get_box(s0_cart(n), i))
          select case (dm)
          case (2)
             call bl_error("Only call put_1d_array_on_cart_irreg for 3D spherical!")
          case (3)
             if (spherical .eq. 0) then
                call bl_error("Only call put_1d_array_on_cart_irreg for 3D spherical!")
             else
                call put_1d_array_on_cart_irreg_sphr(s0(1,:),radii,sp(:,:,:,:),lo,hi, &
                                                     dx(n,:),ng_s)
             endif
          end select
       end do

    enddo

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

    call destroy(bpt)
    
  end subroutine put_1d_array_on_cart_irreg

  subroutine put_1d_array_on_cart_irreg_sphr(s0,radii,s0_cart,lo,hi,dx,ng_s)

    use bl_constants_module
    use geometry, only: center, nr_irreg
    use probin_module, only: prob_lo

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(in   ) :: s0(0:),radii(0:)
    real(kind=dp_t), intent(inout) :: s0_cart(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer         :: i,j,k,index
    real(kind=dp_t) :: x,y,z
    real(kind=dp_t) :: radius

    type(bl_prof_timer), save :: bpt

    call build(bpt, "put_1d_array_on_cart_irreg_3d_sphr")

    !$OMP PARALLEL DO PRIVATE(i,j,k,x,y,z,radius,index)
    do k = lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+HALF)*dx(3) - center(3)
       do j = lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+HALF)*dx(2) - center(2)
          do i = lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+HALF)*dx(1) - center(1)
             radius = sqrt(x**2 + y**2 + z**2)

             ! figure out which radii index this point maps into
             index = ((radius / dx(1))**2 - 0.75d0) / 2.d0
                
             ! due to roundoff error, need to ensure that we are in the proper radial bin
             if (index .lt. nr_irreg) then
                if (abs(radius-radii(index)) .gt. abs(radius-radii(index+1))) then
                   index = index+1
                end if
             end if

             s0_cart(i,j,k,1) = s0(index)

          end do
       end do
    end do
    !$OMP END PARALLEL DO

    call destroy(bpt)

  end subroutine put_1d_array_on_cart_irreg_sphr

end module fill_3d_module
