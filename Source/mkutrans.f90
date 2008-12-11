module mkutrans_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use define_bc_module

  implicit none

  private

  public :: mkutrans

contains

  subroutine mkutrans(u,utrans,w0,w0mac,dx,dt,the_bc_level)

    use bl_prof_module
    use create_umac_grown_module
    use geometry, only: dm, nlevs, spherical

    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(inout) :: utrans(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: utp(:,:,:,:)
    real(kind=dp_t), pointer :: vtp(:,:,:,:)
    real(kind=dp_t), pointer :: wtp(:,:,:,:)
    real(kind=dp_t), pointer :: w0xp(:,:,:,:)
    real(kind=dp_t), pointer :: w0yp(:,:,:,:)
    real(kind=dp_t), pointer :: w0zp(:,:,:,:)
    integer                  :: lo(dm),hi(dm)
    integer                  :: i,n,ng_u,ng_ut,ng_w0

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mkutrans")

    ng_u  = u(1)%ng
    ng_ut = utrans(1,1)%ng
    ng_w0 = w0mac(1,1)%ng

    do n=1,nlevs

       do i=1,u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          up => dataptr(u(n),i)
          utp => dataptr(utrans(n,1),i)
          vtp => dataptr(utrans(n,2),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call mkutrans_2d(n,up(:,:,1,:), ng_u, &
                              utp(:,:,1,1), vtp(:,:,1,1), ng_ut, w0(n,:), &
                              lo,hi,dx(n,:),dt,&
                              the_bc_level(n)%adv_bc_level_array(i,:,:,:), &
                              the_bc_level(n)%phys_bc_level_array(i,:,:))
          case (3)
             wtp => dataptr(utrans(n,3), i)
             w0xp => dataptr(w0mac(n,1), i)
             w0yp => dataptr(w0mac(n,2), i)
             w0zp => dataptr(w0mac(n,3), i)
             if (spherical .eq. 1) then
                call mkutrans_3d(n, up(:,:,:,:), ng_u, &
                                 utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), ng_ut, &
                                 w0(1,:), w0xp(:,:,:,1), w0yp(:,:,:,1), w0zp(:,:,:,1),&
                                 ng_w0, lo, hi, dx(n,:), dt, &
                                 the_bc_level(n)%adv_bc_level_array(i,:,:,:), &
                                 the_bc_level(n)%phys_bc_level_array(i,:,:))
             else
                call mkutrans_3d(n, up(:,:,:,:), ng_u, &
                                 utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), ng_ut, &
                                 w0(n,:), w0xp(:,:,:,1), w0yp(:,:,:,1), w0zp(:,:,:,1),&
                                 ng_w0, lo, hi, dx(n,:), dt, &
                                 the_bc_level(n)%adv_bc_level_array(i,:,:,:), &
                                 the_bc_level(n)%phys_bc_level_array(i,:,:))
             end if
          end select
       end do

    end do

    if (nlevs .gt. 1) then
       do n=2,nlevs
          call create_umac_grown(n,utrans(n,:),utrans(n-1,:))
       end do
    else
       do n=1,nlevs
          do i=1,dm
             call multifab_fill_boundary(utrans(n,i))
          enddo
       end do
    end if
    
    ! we don't need calls to multifab_physbc or multifab_fill_ghost cells since the boundary 
    ! conditions are handled within mkutrans_2d and _3d.
    ! I don't think a call to ml_edge_restriction makes sense here

    call destroy(bpt)

  end subroutine mkutrans

  subroutine mkutrans_2d(n,u,ng_u,utrans,vtrans,ng_ut,w0,lo,hi,dx,dt,adv_bc,phys_bc)

    use bc_module
    use slope_module
    use geometry, only: nr
    use variables, only: rel_eps
    use probin_module, only: use_ppm
    use ppm_module

    integer,         intent(in   ) :: n,lo(:),hi(:),ng_u,ng_ut
    real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u :,lo(2)-ng_u :,:)
    real(kind=dp_t), intent(inout) :: utrans(lo(1)-ng_ut:,lo(2)-ng_ut:)
    real(kind=dp_t), intent(inout) :: vtrans(lo(1)-ng_ut:,lo(2)-ng_ut:)
    real(kind=dp_t), intent(in   ) :: w0(0:)    
    real(kind=dp_t), intent(in   ) :: dt,dx(:)
    integer        , intent(in   ) :: adv_bc(:,:,:)
    integer        , intent(in   ) :: phys_bc(:,:)
    
    real(kind=dp_t) :: slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1)
    real(kind=dp_t) :: slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1)

    real(kind=dp_t), allocatable :: Ipx(:,:,:)
    real(kind=dp_t), allocatable :: Imx(:,:,:)
    real(kind=dp_t), allocatable :: Ipy(:,:,:)
    real(kind=dp_t), allocatable :: Imy(:,:,:)
    
    real(kind=dp_t), allocatable :: ulx(:,:),urx(:,:)
    real(kind=dp_t), allocatable :: vly(:,:),vry(:,:)

    real(kind=dp_t) hx,hy,dt2,uavg,vlo,vhi

    integer :: i,j,is,js,ie,je

    logical :: test
    
    allocate(ulx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(urx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1))

    allocate(vly(lo(1)-1:hi(1)+1,lo(2):hi(2)+1))
    allocate(vry(lo(1)-1:hi(1)+1,lo(2):hi(2)+1))

    allocate(Ipx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(Imx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(Ipy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(Imy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))

    is = lo(1)
    js = lo(2)
    ie = hi(1)
    je = hi(2)
    
    dt2 = HALF*dt
    
    hx = dx(1)
    hy = dx(2)
    
    if (use_ppm) then
       call ppm_2d(u(:,:,1),ng_u,u,ng_u,Ipx,Imx,lo,hi,adv_bc(:,:,1:),dx,dt)
       call ppm_2d(u(:,:,2),ng_u,u,ng_u,Ipy,Imy,lo,hi,adv_bc(:,:,2:),dx,dt)
    else
       call slopex_2d(u(:,:,1:),slopex,lo,hi,ng_u,1,adv_bc(:,:,1:))
       call slopey_2d(u(:,:,2:),slopey,lo,hi,ng_u,1,adv_bc(:,:,2:))
    end if

    !******************************************************************
    ! create utrans
    !******************************************************************

    if (use_ppm) then
       do j=js,je
          do i=is,ie+1
             ! extrapolate to edges
             if (u(i-1,j,1) .gt. ZERO) then
                ulx(i,j) = u(i-1,j,1) + Ipx(i-1,j,1)
             else
                ulx(i,j) = u(i-1,j,1)
             end if
             if (u(i,j,1) .lt. ZERO) then
                urx(i,j) = u(i  ,j,1) + Imx(i  ,j,1)
             else
                urx(i,j) = u(i  ,j,1)
             end if
          end do
       end do
    else
       do j=js,je
          do i=is,ie+1
             ! extrapolate to edges
             ulx(i,j) = u(i-1,j,1) + (HALF-(dt2/hx)*max(ZERO,u(i-1,j,1)))*slopex(i-1,j,1)
             urx(i,j) = u(i  ,j,1) - (HALF+(dt2/hx)*min(ZERO,u(i  ,j,1)))*slopex(i  ,j,1)
          end do
       end do
    end if

    ! impose lo side bc's
    if (phys_bc(1,1) .eq. INLET) then
       ulx(is,js:je) = u(is-1,js:je,1)
       urx(is,js:je) = u(is-1,js:je,1)
    else if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
       ulx(is,js:je) = ZERO
       urx(is,js:je) = ZERO
    else if (phys_bc(1,1) .eq. OUTLET) then
       ulx(is,js:je) = min(urx(is,js:je),ZERO)
       urx(is,js:je) = min(urx(is,js:je),ZERO)
    end if

    ! impose hi side bc's    
    if (phys_bc(1,2) .eq. INLET) then
       ulx(ie+1,js:je) = u(ie+1,js:je,1)
       urx(ie+1,js:je) = u(ie+1,js:je,1)
    else if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
       ulx(ie+1,js:je) = ZERO
       urx(ie+1,js:je) = ZERO
    else if (phys_bc(1,2) .eq. OUTLET) then
       ulx(ie+1,js:je) = max(ulx(ie+1,js:je),ZERO)
       urx(ie+1,js:je) = max(ulx(ie+1,js:je),ZERO)
    end if

    do j=js,je
       do i=is,ie+1
          ! upwind
          uavg = HALF*(ulx(i,j)+urx(i,j))
          test = ((ulx(i,j) .le. ZERO .and. urx(i,j) .ge. ZERO) .or. &
               (abs(ulx(i,j)+urx(i,j)) .lt. rel_eps))
          utrans(i,j) = merge(ulx(i,j),urx(i,j),uavg .gt. ZERO)
          utrans(i,j) = merge(ZERO,utrans(i,j),test)
       end do
    end do

    !******************************************************************
    ! create vtrans
    !******************************************************************

    if (use_ppm) then
       do j=js,je+1
          ! compute effect of w0
          if (j .eq. 0) then
             vlo = w0(j)
             vhi = HALF*(w0(j)+w0(j+1))
          else if (j .eq. nr(n)) then
             vlo = HALF*(w0(j-1)+w0(j))
             vhi = w0(j)
          else
             vlo = HALF*(w0(j-1)+w0(j))
             vhi = HALF*(w0(j)+w0(j+1))
          end if
          do i=is,ie
             ! extrapolate to edges
             if (u(i,j-1,2)+vlo .gt. ZERO) then
                vly(i,j) = u(i,j-1,2) + Ipy(i,j-1,2)
             else
                vly(i,j) = u(i,j-1,2)
             end if
             if (u(i,j,2)+vhi .lt. ZERO) then
                vry(i,j) = u(i,j  ,2) + Imy(i,j,2)
             else
                vry(i,j) = u(i,j  ,2)
             end if
          end do
       end do
    else
       do j=js,je+1
          ! compute effect of w0
          if (j .eq. 0) then
             vlo = w0(j)
             vhi = HALF*(w0(j)+w0(j+1))
          else if (j .eq. nr(n)) then
             vlo = HALF*(w0(j-1)+w0(j))
             vhi = w0(j)
          else
             vlo = HALF*(w0(j-1)+w0(j))
             vhi = HALF*(w0(j)+w0(j+1))
          end if
          do i=is,ie
             ! extrapolate to edges
             vly(i,j) = u(i,j-1,2) + (HALF-(dt2/hy)*max(ZERO,u(i,j-1,2)+vlo))*slopey(i,j-1,1)
             vry(i,j) = u(i,j  ,2) - (HALF+(dt2/hy)*min(ZERO,u(i,j  ,2)+vhi))*slopey(i,j  ,1)
          end do
       end do
    end if

    ! impose lo side bc's
    if (phys_bc(2,1) .eq. INLET) then
       vly(is:ie,js) = u(is:ie,js-1,2)
       vry(is:ie,js) = u(is:ie,js-1,2)
    else if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
       vly(is:ie,js) = ZERO
       vry(is:ie,js) = ZERO
    else if (phys_bc(2,1) .eq. OUTLET) then
       vly(is:ie,js) = min(vry(is:ie,js),ZERO)
       vry(is:ie,js) = min(vry(is:ie,js),ZERO)
    end if

    ! impose hi side bc's
    if (phys_bc(2,2) .eq. INLET) then
       vly(is:ie,je+1) = u(is:ie,je+1,2)
       vry(is:ie,je+1) = u(is:ie,je+1,2)
    else if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
       vly(is:ie,je+1) = ZERO
       vry(is:ie,je+1) = ZERO
    else if (phys_bc(2,2) .eq. OUTLET) then
       vly(is:ie,je+1) = max(vly(is:ie,je+1),ZERO)
       vry(is:ie,je+1) = max(vly(is:ie,je+1),ZERO)
    end if

    do j=js,je+1
       do i=is,ie
          ! upwind using full v
          uavg = HALF*(vly(i,j)+vry(i,j))
          test = ((vly(i,j)+w0(j) .le. ZERO .and. vry(i,j)+w0(j) .ge. ZERO) .or. &
               (abs(vly(i,j)+vry(i,j)+TWO*w0(j)) .lt. rel_eps))
          vtrans(i,j) = merge(vly(i,j),vry(i,j),uavg+w0(j) .gt. ZERO)
          vtrans(i,j) = merge(ZERO,vtrans(i,j),test)
       enddo
    enddo

    deallocate(ulx,urx,vly,vry)
    deallocate(Ipx,Imx,Ipy,Imy)

  end subroutine mkutrans_2d
  
  subroutine mkutrans_3d(n,u,ng_u,utrans,vtrans,wtrans,ng_ut,w0,w0macx,w0macy,w0macz, &
                         ng_w0,lo,hi,dx,dt,adv_bc,phys_bc)

    use bc_module
    use slope_module
    use geometry, only: nr, spherical
    use variables, only: rel_eps
    
    integer,         intent(in)    :: n,lo(:),hi(:),ng_u,ng_ut,ng_w0    
    real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :,:)
    real(kind=dp_t), intent(inout) :: utrans(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(inout) :: vtrans(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(inout) :: wtrans(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: dt,dx(:)
    integer        , intent(in   ) :: adv_bc(:,:,:)
    integer        , intent(in   ) :: phys_bc(:,:)
    
    real(kind=dp_t) :: slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1)
    real(kind=dp_t) :: slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1)
    real(kind=dp_t) :: slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1)
    
    real(kind=dp_t) hx,hy,hz,dt2,uavg,uhi,ulo,vhi,vlo,whi,wlo
    
    logical :: test

    integer :: i,j,k,is,js,ks,ie,je,ke
    
    real(kind=dp_t), allocatable:: ulx(:,:,:),urx(:,:,:)
    real(kind=dp_t), allocatable:: vly(:,:,:),vry(:,:,:)
    real(kind=dp_t), allocatable:: wlz(:,:,:),wrz(:,:,:)

    allocate(ulx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(urx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    allocate(vly(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(vry(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))

    allocate(wlz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(wrz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))

    is = lo(1)
    js = lo(2)
    ks = lo(3)
    ie = hi(1)
    je = hi(2)
    ke = hi(3)
    
    dt2 = HALF*dt
    
    hx = dx(1)
    hy = dx(2)
    hz = dx(3)
    
    do k = lo(3)-1,hi(3)+1
       call slopex_2d(u(:,:,k,1:),slopex(:,:,k,:),lo,hi,ng_u,1,adv_bc(:,:,1:))
       call slopey_2d(u(:,:,k,2:),slopey(:,:,k,:),lo,hi,ng_u,1,adv_bc(:,:,2:))
    end do
    call slopez_3d(u(:,:,:,3:),slopez,lo,hi,ng_u,1,adv_bc(:,:,3:))
    
    !******************************************************************
    ! create utrans
    !******************************************************************

    do k=ks,ke
       do j=js,je
          do i=is,ie+1
             ! compute effect of w0
             if (spherical .eq. 1) then
                ulo = u(i-1,j,k,1) + HALF * (w0macx(i-1,j,k)+w0macx(i  ,j,k))
                uhi = u(i  ,j,k,1) + HALF * (w0macx(i  ,j,k)+w0macx(i+1,j,k))
             else
                ulo = u(i-1,j,k,1)
                uhi = u(i  ,j,k,1)
             end if
             
             ! extrapolate to edges
             ulx(i,j,k) = u(i-1,j,k,1) + (HALF - dt2*max(ZERO,ulo)/hx)*slopex(i-1,j,k,1)
             urx(i,j,k) = u(i  ,j,k,1) - (HALF + dt2*min(ZERO,uhi)/hx)*slopex(i  ,j,k,1)
          end do
       end do
    end do

    ! impose lo side bc's
    if (phys_bc(1,1) .eq. INLET) then
       ulx(is,js:je,ks:ke) = u(is-1,js:je,ks:ke,1)
       urx(is,js:je,ks:ke) = u(is-1,js:je,ks:ke,1)
    else if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
       ulx(is,js:je,ks:ke) = ZERO
       urx(is,js:je,ks:ke) = ZERO
    else if (phys_bc(1,1) .eq. OUTLET) then
       ulx(is,js:je,ks:ke) = min(urx(is,js:je,ks:ke),ZERO)
       urx(is,js:je,ks:ke) = min(urx(is,js:je,ks:ke),ZERO)
    end if

    ! impose hi side bc's
    if (phys_bc(1,2) .eq. INLET) then
       ulx(ie+1,js:je,ks:ke) = u(ie+1,js:je,ks:ke,1)
       urx(ie+1,js:je,ks:ke) = u(ie+1,js:je,ks:ke,1)
    else if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
       ulx(ie+1,js:je,ks:ke) = ZERO
       urx(ie+1,js:je,ks:ke) = ZERO
    else if (phys_bc(1,2) .eq. OUTLET) then
       ulx(ie+1,js:je,ks:ke) = max(ulx(ie+1,js:je,ks:ke),ZERO)
       urx(ie+1,js:je,ks:ke) = max(ulx(ie+1,js:je,ks:ke),ZERO)
    end if

    do k=ks,ke
       do j=js,je
          do i=is,ie+1
             if (spherical .eq. 1) then
                ! upwind using full u
                uavg = HALF*(ulx(i,j,k)+urx(i,j,k))
                test = ((ulx(i,j,k)+w0macx(i,j,k) .le. ZERO .and. &
                     urx(i,j,k)+w0macx(i,j,k) .ge. ZERO) .or. &
                     (abs(ulx(i,j,k)+urx(i,j,k)+TWO*w0macx(i,j,k)) .lt. rel_eps))
                utrans(i,j,k) = merge(ulx(i,j,k),urx(i,j,k),uavg+w0macx(i,j,k) .gt. ZERO)
                utrans(i,j,k) = merge(ZERO,utrans(i,j,k),test)
             else
                ! upwind
                uavg = HALF*(ulx(i,j,k)+urx(i,j,k))
                test = ((ulx(i,j,k) .le. ZERO .and. urx(i,j,k) .ge. ZERO) .or. &
                     (abs(ulx(i,j,k)+urx(i,j,k)) .lt. rel_eps))
                utrans(i,j,k) = merge(ulx(i,j,k),urx(i,j,k),uavg .gt. ZERO)
                utrans(i,j,k) = merge(ZERO,utrans(i,j,k),test)
             end if
          enddo
       enddo
    enddo

    !******************************************************************
    ! create vtrans
    !******************************************************************

    do k=ks,ke
       do j=js,je+1
          do i=is,ie
             ! compute effect of w0
             if (spherical .eq. 1) then
                vlo = u(i,j-1,k,2) + HALF * (w0macy(i,j-1,k)+w0macy(i,j  ,k))
                vhi = u(i,j  ,k,2) + HALF * (w0macy(i,j  ,k)+w0macy(i,j+1,k))
             else
                vlo = u(i,j-1,k,2)
                vhi = u(i,j  ,k,2)
             end if

             ! extrapolate to edges
             vly(i,j,k) = u(i,j-1,k,2) + (HALF - dt2*max(ZERO,vlo)/hy)*slopey(i,j-1,k,1)
             vry(i,j,k) = u(i,j  ,k,2) - (HALF + dt2*min(ZERO,vhi)/hy)*slopey(i,j  ,k,1)
         enddo
       enddo
    enddo

    ! impose lo side bc's
    if (phys_bc(2,1) .eq. INLET) then
       vly(is:ie,js,ks:ke) = u(is:ie,js-1,ks:ke,2)
       vry(is:ie,js,ks:ke) = u(is:ie,js-1,ks:ke,2)
    else if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
       vly(is:ie,js,ks:ke) = ZERO
       vry(is:ie,js,ks:ke) = ZERO
    else if (phys_bc(2,1) .eq. OUTLET) then
       vly(is:ie,js,ks:ke) = min(vry(is:ie,js,ks:ke),ZERO)
       vry(is:ie,js,ks:ke) = min(vry(is:ie,js,ks:ke),ZERO)
    end if

    ! impose hi side bc's
    if (phys_bc(2,2) .eq. INLET) then
       vly(is:ie,je+1,ks:ke) = u(is:ie,je+1,ks:ke,2)
       vry(is:ie,je+1,ks:ke) = u(is:ie,je+1,ks:ke,2)
    else if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
       vly(is:ie,je+1,ks:ke) = ZERO
       vry(is:ie,je+1,ks:ke) = ZERO
    else if (phys_bc(2,2) .eq. OUTLET) then
       vly(is:ie,je+1,ks:ke) = max(vly(is:ie,je+1,ks:ke),ZERO)
       vry(is:ie,je+1,ks:ke) = max(vly(is:ie,je+1,ks:ke),ZERO)
    end if

    do k=ks,ke
       do j=js,je+1
          do i=is,ie
             if (spherical .eq. 1) then
                ! upwind using full v
                uavg = HALF*(vly(i,j,k)+vry(i,j,k))
                test = ((vly(i,j,k)+w0macy(i,j,k) .le. ZERO .and. &
                     vry(i,j,k)+w0macy(i,j,k) .ge. ZERO) .or. &
                     (abs(vly(i,j,k)+vry(i,j,k)+TWO*w0macy(i,j,k)) .lt. rel_eps))
                vtrans(i,j,k) = merge(vly(i,j,k),vry(i,j,k),uavg+w0macy(i,j,k) .gt. ZERO)
                vtrans(i,j,k) = merge(ZERO,vtrans(i,j,k),test)
             else
                ! upwind
                uavg = HALF*(vly(i,j,k)+vry(i,j,k))
                test = ((vly(i,j,k) .le. ZERO .and. vry(i,j,k) .ge. ZERO) .or. &
                     (abs(vly(i,j,k)+vry(i,j,k)) .lt. rel_eps))
                vtrans(i,j,k) = merge(vly(i,j,k),vry(i,j,k),uavg .gt. ZERO)
                vtrans(i,j,k) = merge(ZERO,vtrans(i,j,k),test)
             end if
          enddo
       enddo
    enddo

    !******************************************************************
    ! create wtrans
    !******************************************************************

    do k=ks,ke+1
       do j=js,je
          do i=is,ie
             ! compute effect of w0
             if (spherical .eq. 1) then
                wlo = u(i,j,k-1,3) + HALF * (w0macz(i,j,k-1)+w0macz(i,j,k  ))
                whi = u(i,j,k  ,3) + HALF * (w0macz(i,j,k  )+w0macz(i,j,k+1))
             else
                if (k .eq. 0) then
                   wlo = u(i,j,k-1,3) + w0(k)
                   whi = u(i,j,k  ,3) + HALF*(w0(k)+w0(k+1))
                else if (k .eq. nr(n)) then
                   wlo = u(i,j,k-1,3) + HALF*(w0(k-1)+w0(k))
                   whi = u(i,j,k  ,3) + w0(k)
                else
                   wlo = u(i,j,k-1,3) + HALF*(w0(k-1)+w0(k))
                   whi = u(i,j,k  ,3) + HALF*(w0(k)+w0(k+1))
                end if
             end if

             ! extrapolate to edges
             wlz(i,j,k) = u(i,j,k-1,3) + (HALF - dt2*max(ZERO,wlo)/hz)*slopez(i,j,k-1,1)
             wrz(i,j,k) = u(i,j,k  ,3) - (HALF + dt2*min(ZERO,whi)/hz)*slopez(i,j,k  ,1)
          end do
       end do
    end do
    
    ! impose lo side bc's
    if (phys_bc(3,1) .eq. INLET) then
       wlz(is:ie,js:je,ks) = u(is:is,js:je,ks-1,3)
       wrz(is:ie,js:je,ks) = u(is:is,js:je,ks-1,3)
    else if (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
       wlz(is:ie,js:je,ks) = ZERO
       wrz(is:ie,js:je,ks) = ZERO
    else if (phys_bc(3,1) .eq. OUTLET) then
       wlz(is:ie,js:je,ks) = min(wrz(is:ie,js:je,ks),ZERO)
       wrz(is:ie,js:je,ks) = min(wrz(is:ie,js:je,ks),ZERO)
    end if

    ! impose hi side bc's
    if (phys_bc(3,2) .eq. INLET) then
       wlz(is:ie,js:je,ke+1) = u(is:ie,js:je,ke+1,3)
       wrz(is:ie,js:je,ke+1) = u(is:ie,js:je,ke+1,3)
    else if (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
       wlz(is:ie,js:je,ke+1) = ZERO
       wrz(is:ie,js:je,ke+1) = ZERO
    else if (phys_bc(3,2) .eq. OUTLET) then
       wlz(is:ie,js:je,ke+1) = max(wlz(is:ie,js:je,ke+1),ZERO)
       wrz(is:ie,js:je,ke+1) = max(wlz(is:ie,js:je,ke+1),ZERO)
    end if

    do k=ks,ke+1
       do j=js,je
          do i=is,ie
             if (spherical .eq. 1) then
                ! upwind using full w
                uavg = HALF*(wlz(i,j,k)+wrz(i,j,k))
                test = ((wlz(i,j,k)+w0macz(i,j,k) .le. ZERO .and. &
                     wrz(i,j,k)+w0macz(i,j,k) .ge. ZERO) .or. &
                     (abs(wlz(i,j,k)+wrz(i,j,k)+TWO*w0macz(i,j,k)) .lt. rel_eps))
                wtrans(i,j,k) = merge(wlz(i,j,k),wrz(i,j,k),uavg+w0macz(i,j,k) .gt. ZERO)
                wtrans(i,j,k) = merge(ZERO,wtrans(i,j,k),test)
             else
                ! upwind using full w
                uavg = HALF*(wlz(i,j,k)+wrz(i,j,k))
                test = ((wlz(i,j,k)+w0(k).le.ZERO .and. wrz(i,j,k)+w0(k).ge.ZERO) .or. &
                     (abs(wlz(i,j,k)+wrz(i,j,k)+TWO*w0(k)) .lt. rel_eps))
                wtrans(i,j,k) = merge(wlz(i,j,k),wrz(i,j,k),uavg+w0(k) .gt. ZERO)
                wtrans(i,j,k) = merge(ZERO,wtrans(i,j,k),test)
             end if
          enddo
       enddo
    enddo

    deallocate(ulx,urx,vly,vry,wlz,wrz)

  end subroutine mkutrans_3d
  
end module mkutrans_module
