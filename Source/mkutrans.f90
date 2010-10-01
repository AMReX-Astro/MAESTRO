module mkutrans_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: mkutrans

contains

  subroutine mkutrans(u,ufull,utrans,w0,w0mac,dx,dt,the_bc_level,mla)

    use bl_prof_module
    use create_umac_grown_module
    use geometry, only: dm, nlevs, spherical
    use ml_restriction_module, only: ml_edge_restriction_c
    use multifab_physbc_edgevel_module, only: multifab_physbc_edgevel

    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(in   ) :: ufull(:)
    type(multifab) , intent(inout) :: utrans(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(in   ) :: mla

    ! local
    real(kind=dp_t), pointer :: up(:,:,:,:)
    real(kind=dp_t), pointer :: ufp(:,:,:,:)
    real(kind=dp_t), pointer :: utp(:,:,:,:)
    real(kind=dp_t), pointer :: vtp(:,:,:,:)
    real(kind=dp_t), pointer :: wtp(:,:,:,:)
    real(kind=dp_t), pointer :: w0xp(:,:,:,:)
    real(kind=dp_t), pointer :: w0yp(:,:,:,:)
    real(kind=dp_t), pointer :: w0zp(:,:,:,:)
    integer                  :: lo(dm),hi(dm)
    integer                  :: i,n,n_1d,ng_u,ng_uf,ng_ut,ng_w0

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mkutrans")

    ng_u  = nghost(u(1))
    ng_uf = nghost(ufull(1))
    ng_ut = nghost(utrans(1,1))
    ng_w0 = nghost(w0mac(1,1))

    do n=1,nlevs

       do i=1, nboxes(u(n))
          if ( multifab_remote(u(n),i) ) cycle
          up => dataptr(u(n),i)
          ufp => dataptr(ufull(n),i)
          utp => dataptr(utrans(n,1),i)
          vtp => dataptr(utrans(n,2),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call mkutrans_2d(up(:,:,1,:), ng_u, &
                              ufp(:,:,1,:), ng_uf, &
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
                n_1d = 1
             else
                n_1d = n
             end if
             call mkutrans_3d(up(:,:,:,:), ng_u, &
                              ufp(:,:,:,:), ng_uf, &
                              utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), ng_ut, &
                              w0(n_1d,:), w0xp(:,:,:,1), w0yp(:,:,:,1), w0zp(:,:,:,1),&
                              ng_w0, lo, hi, dx(n,:), dt, &
                              the_bc_level(n)%adv_bc_level_array(i,:,:,:), &
                              the_bc_level(n)%phys_bc_level_array(i,:,:))
          end select
       end do

    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
      do i=1,dm
          call multifab_fill_boundary(utrans(nlevs,i))
       enddo

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc_edgevel(utrans(nlevs,:),the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          do i=1,dm
             call ml_edge_restriction_c(utrans(n-1,i),1,utrans(n,i),1,mla%mba%rr(n-1,:),i,1)
          enddo

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc_edgevel are called for
          ! level n and level 1 (if n=2)
          call create_umac_grown(n,utrans(n,:),utrans(n-1,:),the_bc_level(n-1),the_bc_level(n))

       end do

    end if

    call destroy(bpt)

  end subroutine mkutrans

  subroutine mkutrans_2d(u,ng_u,ufull,ng_uf,utrans,vtrans,ng_ut,w0, &
                         lo,hi,dx,dt,adv_bc,phys_bc)

    use bc_module
    use slope_module
    use variables, only: rel_eps
    use probin_module, only: ppm_type
    use ppm_module

    integer,         intent(in   ) :: lo(:),hi(:),ng_u,ng_uf,ng_ut
    real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u :,lo(2)-ng_u :,:)
    real(kind=dp_t), intent(in   ) ::  ufull(lo(1)-ng_u :,lo(2)-ng_u :,:)
    real(kind=dp_t), intent(inout) :: utrans(lo(1)-ng_ut:,lo(2)-ng_ut:)
    real(kind=dp_t), intent(inout) :: vtrans(lo(1)-ng_ut:,lo(2)-ng_ut:)
    real(kind=dp_t), intent(in   ) :: w0(0:)    
    real(kind=dp_t), intent(in   ) :: dt,dx(:)
    integer        , intent(in   ) :: adv_bc(:,:,:)
    integer        , intent(in   ) :: phys_bc(:,:)
    
    real(kind=dp_t) :: slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1)
    real(kind=dp_t) :: slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1)

    real(kind=dp_t), allocatable :: Ip(:,:,:)
    real(kind=dp_t), allocatable :: Im(:,:,:)
    
    real(kind=dp_t), allocatable :: ulx(:,:),urx(:,:)
    real(kind=dp_t), allocatable :: vly(:,:),vry(:,:)

    real(kind=dp_t) hx,hy,dt2,uavg

    integer :: i,j,is,js,ie,je

    logical :: test
    
    allocate(ulx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(urx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1))

    allocate(vly(lo(1)-1:hi(1)+1,lo(2):hi(2)+1))
    allocate(vry(lo(1)-1:hi(1)+1,lo(2):hi(2)+1))

    allocate(Ip(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(Im(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))

    is = lo(1)
    js = lo(2)
    ie = hi(1)
    je = hi(2)
    
    dt2 = HALF*dt
    
    hx = dx(1)
    hy = dx(2)
    
    if (ppm_type .eq. 0) then
       call slopex_2d(u(:,:,1:),slopex,lo,hi,ng_u,1,adv_bc(:,:,1:))
       call slopey_2d(u(:,:,2:),slopey,lo,hi,ng_u,1,adv_bc(:,:,2:))
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_2d(u(:,:,1),ng_u,ufull,ng_uf,Ip,Im,lo,hi,adv_bc(:,:,1),dx,dt)
    end if

    !******************************************************************
    ! create utrans
    !******************************************************************

    if (ppm_type .eq. 0) then
       do j=js,je
          do i=is,ie+1
             ! extrapolate to edges
             ulx(i,j) = u(i-1,j,1) + (HALF-(dt2/hx)*max(ZERO,ufull(i-1,j,1)))*slopex(i-1,j,1)
             urx(i,j) = u(i  ,j,1) - (HALF+(dt2/hx)*min(ZERO,ufull(i  ,j,1)))*slopex(i  ,j,1)
          end do
       end do
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do j=js,je
          do i=is,ie+1
             ! extrapolate to edges
             ulx(i,j) = Ip(i-1,j,1)
             urx(i,j) = Im(i  ,j,1)
          end do
       end do
    end if

    ! impose lo i side bc's
    select case(phys_bc(1,1))
    case (INLET)
       ulx(is,js:je) = u(is-1,js:je,1)
       urx(is,js:je) = u(is-1,js:je,1)
    case (SLIP_WALL, NO_SLIP_WALL, SYMMETRY)
       ulx(is,js:je) = ZERO
       urx(is,js:je) = ZERO
    case (OUTLET)
       ulx(is,js:je) = min(urx(is,js:je),ZERO)
       urx(is,js:je) = ulx(is,js:je)
    case (INTERIOR, PERIODIC) 
    case  default
       call bl_error("mkutrans_2d: invalid boundary type phys_bc(1,1)")
    end select

    ! impose hi i side bc's    
    select case(phys_bc(1,2))
    case (INLET)
       ulx(ie+1,js:je) = u(ie+1,js:je,1)
       urx(ie+1,js:je) = u(ie+1,js:je,1)
    case (SLIP_WALL, NO_SLIP_WALL, SYMMETRY)
       ulx(ie+1,js:je) = ZERO
       urx(ie+1,js:je) = ZERO
    case (OUTLET)
       ulx(ie+1,js:je) = max(ulx(ie+1,js:je),ZERO)
       urx(ie+1,js:je) = ulx(ie+1,js:je)
    case (INTERIOR, PERIODIC) 
    case  default
       call bl_error("mkutrans_2d: invalid boundary type phys_bc(1,2)")
    end select

    do j=js,je
       do i=is,ie+1
          ! solve Riemann problem using full velocity
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

    if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_2d(u(:,:,2),ng_u,ufull,ng_uf,Ip,Im,lo,hi,adv_bc(:,:,2),dx,dt)
    end if
       
    if (ppm_type .eq. 0) then
       do j=js,je+1
          do i=is,ie
             ! extrapolate to edges
             vly(i,j) = u(i,j-1,2) + (HALF-(dt2/hy)*max(ZERO,ufull(i,j-1,2)))*slopey(i,j-1,1)
             vry(i,j) = u(i,j  ,2) - (HALF+(dt2/hy)*min(ZERO,ufull(i,j  ,2)))*slopey(i,j  ,1)
          end do
       end do

    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do j=js,je+1
          do i=is,ie
             ! extrapolate to edges
             vly(i,j) = Ip(i,j-1,2)
             vry(i,j) = Im(i,j  ,2)
          end do
       end do
    end if

    ! impose lo side bc's
    select case(phys_bc(2,1))
    case (INLET)
       vly(is:ie,js) = u(is:ie,js-1,2)
       vry(is:ie,js) = u(is:ie,js-1,2)
    case (SLIP_WALL, NO_SLIP_WALL, SYMMETRY)
       vly(is:ie,js) = ZERO
       vry(is:ie,js) = ZERO
    case (OUTLET)
       vly(is:ie,js) = min(vry(is:ie,js),ZERO)
       vry(is:ie,js) = vly(is:ie,js)
    case (INTERIOR, PERIODIC) 
    case  default
       call bl_error("mkutrans_2d: invalid boundary type phys_bc(2,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(2,2))
    case (INLET)
       vly(is:ie,je+1) = u(is:ie,je+1,2)
       vry(is:ie,je+1) = u(is:ie,je+1,2)
    case (SLIP_WALL, NO_SLIP_WALL, SYMMETRY)
       vly(is:ie,je+1) = ZERO
       vry(is:ie,je+1) = ZERO
    case (OUTLET)
       vly(is:ie,je+1) = max(vly(is:ie,je+1),ZERO)
       vry(is:ie,je+1) = vly(is:ie,je+1)
    case (INTERIOR, PERIODIC) 
    case  default
       call bl_error("mkutrans_2d: invalid boundary type phys_bc(2,2)")
    end select

    do j=js,je+1
       do i=is,ie
          ! solve Riemann problem using full velocity
          uavg = HALF*(vly(i,j)+vry(i,j))
          test = ((vly(i,j)+w0(j) .le. ZERO .and. vry(i,j)+w0(j) .ge. ZERO) .or. &
               (abs(vly(i,j)+vry(i,j)+TWO*w0(j)) .lt. rel_eps))
          vtrans(i,j) = merge(vly(i,j),vry(i,j),uavg+w0(j) .gt. ZERO)
          vtrans(i,j) = merge(ZERO,vtrans(i,j),test)
       enddo
    enddo

    deallocate(ulx,urx,vly,vry)
    deallocate(Ip,Im)

  end subroutine mkutrans_2d
  
  subroutine mkutrans_3d(u,ng_u,ufull,ng_uf,utrans,vtrans,wtrans,ng_ut, &
                         w0,w0macx,w0macy,w0macz,ng_w0,lo,hi,dx,dt,adv_bc,phys_bc)

    use bc_module
    use slope_module
    use geometry, only: spherical
    use variables, only: rel_eps
    use probin_module, only: ppm_type
    use ppm_module
    
    integer,         intent(in)    :: lo(:),hi(:),ng_u,ng_uf,ng_ut,ng_w0    
    real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :,:)
    real(kind=dp_t), intent(in   ) ::  ufull(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :,:)
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
    
    real(kind=dp_t), allocatable :: slopex(:,:,:,:)
    real(kind=dp_t), allocatable :: slopey(:,:,:,:)
    real(kind=dp_t), allocatable :: slopez(:,:,:,:)
    
    real(kind=dp_t), allocatable :: Ip(:,:,:,:)
    real(kind=dp_t), allocatable :: Im(:,:,:,:)

    real(kind=dp_t) hx,hy,hz,dt2,uavg
    
    logical :: test

    integer :: i,j,k,is,js,ks,ie,je,ke
    
    real(kind=dp_t), allocatable:: ulx(:,:,:),urx(:,:,:)
    real(kind=dp_t), allocatable:: vly(:,:,:),vry(:,:,:)
    real(kind=dp_t), allocatable:: wlz(:,:,:),wrz(:,:,:)

    allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))
    allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))
    allocate(slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))

    allocate(Ip(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(Im(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

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
    
    if (ppm_type .eq. 0) then
       do k = lo(3)-1,hi(3)+1
          call slopex_2d(u(:,:,k,1:),slopex(:,:,k,:),lo,hi,ng_u,1,adv_bc(:,:,1:))
          call slopey_2d(u(:,:,k,2:),slopey(:,:,k,:),lo,hi,ng_u,1,adv_bc(:,:,2:))
       end do
       call slopez_3d(u(:,:,:,3:),slopez,lo,hi,ng_u,1,adv_bc(:,:,3:))
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_3d(u(:,:,:,1),ng_u,ufull,ng_uf,Ip,Im,lo,hi,adv_bc(:,:,1),dx,dt)
    end if
    
    !******************************************************************
    ! create utrans
    !******************************************************************

    allocate(ulx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(urx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    if (ppm_type .eq. 0) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks,ke
          do j=js,je
             do i=is,ie+1
                ! extrapolate to edges
                ulx(i,j,k) = u(i-1,j,k,1) &
                     + (HALF-(dt2/hx)*max(ZERO,ufull(i-1,j,k,1)))*slopex(i-1,j,k,1)
                urx(i,j,k) = u(i  ,j,k,1) &
                     - (HALF+(dt2/hx)*min(ZERO,ufull(i  ,j,k,1)))*slopex(i  ,j,k,1)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do k=ks,ke
          do j=js,je
             do i=is,ie+1
                ! extrapolate to edges
                ulx(i,j,k) = Ip(i-1,j,k,1)
                urx(i,j,k) = Im(i  ,j,k,1)
             end do
          end do
       end do
    end if

    deallocate(slopex)

    ! impose lo side bc's
    select case(phys_bc(1,1))
    case (INLET)
       ulx(is,js:je,ks:ke) = u(is-1,js:je,ks:ke,1)
       urx(is,js:je,ks:ke) = u(is-1,js:je,ks:ke,1)
    case (SLIP_WALL, NO_SLIP_WALL, SYMMETRY)
       ulx(is,js:je,ks:ke) = ZERO
       urx(is,js:je,ks:ke) = ZERO
    case (OUTLET)
       ulx(is,js:je,ks:ke) = min(urx(is,js:je,ks:ke),ZERO)
       urx(is,js:je,ks:ke) = ulx(is,js:je,ks:ke)
    case (INTERIOR, PERIODIC) 
    case  default
       call bl_error("mkutrans_3d: invalid boundary type phys_bc(1,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(1,2))
    case (INLET)
       ulx(ie+1,js:je,ks:ke) = u(ie+1,js:je,ks:ke,1)
       urx(ie+1,js:je,ks:ke) = u(ie+1,js:je,ks:ke,1)
    case (SLIP_WALL, NO_SLIP_WALL, SYMMETRY)
       ulx(ie+1,js:je,ks:ke) = ZERO
       urx(ie+1,js:je,ks:ke) = ZERO
    case (OUTLET)
       ulx(ie+1,js:je,ks:ke) = max(ulx(ie+1,js:je,ks:ke),ZERO)
       urx(ie+1,js:je,ks:ke) = ulx(ie+1,js:je,ks:ke)
    case (INTERIOR, PERIODIC) 
    case  default
       call bl_error("mkutrans_3d: invalid boundary type phys_bc(1,2)")
    end select

    if (spherical .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke
          do j=js,je
             do i=is,ie+1
                ! solve Riemann problem using full velocity
                uavg = HALF*(ulx(i,j,k)+urx(i,j,k))
                test = ((ulx(i,j,k)+w0macx(i,j,k) .le. ZERO .and. &
                     urx(i,j,k)+w0macx(i,j,k) .ge. ZERO) .or. &
                     (abs(ulx(i,j,k)+urx(i,j,k)+TWO*w0macx(i,j,k)) .lt. rel_eps))
                utrans(i,j,k) = merge(ulx(i,j,k),urx(i,j,k),uavg+w0macx(i,j,k) .gt. ZERO)
                utrans(i,j,k) = merge(ZERO,utrans(i,j,k),test)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke
          do j=js,je
             do i=is,ie+1
                ! solve Riemann problem using full velocity
                uavg = HALF*(ulx(i,j,k)+urx(i,j,k))
                test = ((ulx(i,j,k) .le. ZERO .and. urx(i,j,k) .ge. ZERO) .or. &
                     (abs(ulx(i,j,k)+urx(i,j,k)) .lt. rel_eps))
                utrans(i,j,k) = merge(ulx(i,j,k),urx(i,j,k),uavg .gt. ZERO)
                utrans(i,j,k) = merge(ZERO,utrans(i,j,k),test)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    end if

    deallocate(ulx,urx)

    !******************************************************************
    ! create vtrans
    !******************************************************************

    if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_3d(u(:,:,:,2),ng_u,ufull,ng_uf,Ip,Im,lo,hi,adv_bc(:,:,2),dx,dt)
    end if

    allocate(vly(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(vry(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))

    if (ppm_type .eq. 0) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks,ke
          do j=js,je+1
             do i=is,ie
                ! extrapolate to edges
                vly(i,j,k) = u(i,j-1,k,2) &
                     + (HALF-(dt2/hy)*max(ZERO,ufull(i,j-1,k,2)))*slopey(i,j-1,k,1)
                vry(i,j,k) = u(i,j  ,k,2) &
                     - (HALF+(dt2/hy)*min(ZERO,ufull(i,j  ,k,2)))*slopey(i,j  ,k,1)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do k=ks,ke
          do j=js,je+1
             do i=is,ie
                ! extrapolate to edges
                vly(i,j,k) = Ip(i,j-1,k,2)
                vry(i,j,k) = Im(i,j  ,k,2)
             enddo
          enddo
       enddo
    end if

    deallocate(slopey)

    ! impose lo side bc's
    select case(phys_bc(2,1))
    case (INLET)
       vly(is:ie,js,ks:ke) = u(is:ie,js-1,ks:ke,2)
       vry(is:ie,js,ks:ke) = u(is:ie,js-1,ks:ke,2)
    case (SLIP_WALL, NO_SLIP_WALL, SYMMETRY)
       vly(is:ie,js,ks:ke) = ZERO
       vry(is:ie,js,ks:ke) = ZERO
    case (OUTLET)
       vly(is:ie,js,ks:ke) = min(vry(is:ie,js,ks:ke),ZERO)
       vry(is:ie,js,ks:ke) = vly(is:ie,js,ks:ke)
    case (INTERIOR, PERIODIC) 
    case  default
       call bl_error("mkutrans_3d: invalid boundary type phys_bc(2,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(2,2))
    case (INLET)
       vly(is:ie,je+1,ks:ke) = u(is:ie,je+1,ks:ke,2)
       vry(is:ie,je+1,ks:ke) = u(is:ie,je+1,ks:ke,2)
    case (SLIP_WALL, NO_SLIP_WALL, SYMMETRY)
       vly(is:ie,je+1,ks:ke) = ZERO
       vry(is:ie,je+1,ks:ke) = ZERO
    case (OUTLET)
       vly(is:ie,je+1,ks:ke) = max(vly(is:ie,je+1,ks:ke),ZERO)
       vry(is:ie,je+1,ks:ke) = vly(is:ie,je+1,ks:ke)
    case (INTERIOR, PERIODIC) 
    case  default
       call bl_error("mkutrans_3d: invalid boundary type phys_bc(2,2)")
    end select
    
    if (spherical .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke
          do j=js,je+1
             do i=is,ie
                ! solve Riemann problem using full velocity
                uavg = HALF*(vly(i,j,k)+vry(i,j,k))
                test = ((vly(i,j,k)+w0macy(i,j,k) .le. ZERO .and. &
                     vry(i,j,k)+w0macy(i,j,k) .ge. ZERO) .or. &
                     (abs(vly(i,j,k)+vry(i,j,k)+TWO*w0macy(i,j,k)) .lt. rel_eps))
                vtrans(i,j,k) = merge(vly(i,j,k),vry(i,j,k),uavg+w0macy(i,j,k) .gt. ZERO)
                vtrans(i,j,k) = merge(ZERO,vtrans(i,j,k),test)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke
          do j=js,je+1
             do i=is,ie
                ! solve Riemann problem using full velocity
                uavg = HALF*(vly(i,j,k)+vry(i,j,k))
                test = ((vly(i,j,k) .le. ZERO .and. vry(i,j,k) .ge. ZERO) .or. &
                     (abs(vly(i,j,k)+vry(i,j,k)) .lt. rel_eps))
                vtrans(i,j,k) = merge(vly(i,j,k),vry(i,j,k),uavg .gt. ZERO)
                vtrans(i,j,k) = merge(ZERO,vtrans(i,j,k),test)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    end if

    deallocate(vly,vry)

    !******************************************************************
    ! create wtrans
    !******************************************************************

    if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_3d(u(:,:,:,3),ng_u,ufull,ng_uf,Ip,Im,lo,hi,adv_bc(:,:,3),dx,dt)
    end if

    allocate(wlz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(wrz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))

    if (ppm_type .eq. 0) then
       !$OMP PARALLEL DO PRIVATE(i,j,k)
       do k=ks,ke+1
          do j=js,je
             do i=is,ie
                ! extrapolate to edges
                wlz(i,j,k) = u(i,j,k-1,3) &
                     + (HALF-(dt2/hz)*max(ZERO,ufull(i,j,k-1,3)))*slopez(i,j,k-1,1)
                wrz(i,j,k) = u(i,j,k  ,3) &
                     - (HALF+(dt2/hz)*min(ZERO,ufull(i,j,k  ,3)))*slopez(i,j,k  ,1)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do k=ks,ke+1
          do j=js,je
             do i=is,ie
                ! extrapolate to edges
                wlz(i,j,k) = Ip(i,j,k-1,3)
                wrz(i,j,k) = Im(i,j,k  ,3)
             end do
          end do
       end do
    end if

    deallocate(slopez,Ip,Im)
    
    ! impose lo side bc's
    select case(phys_bc(3,1))
    case (INLET)
       wlz(is:ie,js:je,ks) = u(is:is,js:je,ks-1,3)
       wrz(is:ie,js:je,ks) = u(is:is,js:je,ks-1,3)
    case (SLIP_WALL, NO_SLIP_WALL, SYMMETRY)
       wlz(is:ie,js:je,ks) = ZERO
       wrz(is:ie,js:je,ks) = ZERO
    case (OUTLET)
       wlz(is:ie,js:je,ks) = min(wrz(is:ie,js:je,ks),ZERO)
       wrz(is:ie,js:je,ks) = wlz(is:ie,js:je,ks)
    case (INTERIOR, PERIODIC) 
    case  default
       call bl_error("mkutrans_3d: invalid boundary type phys_bc(3,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(3,2))
    case (INLET)
       wlz(is:ie,js:je,ke+1) = u(is:ie,js:je,ke+1,3)
       wrz(is:ie,js:je,ke+1) = u(is:ie,js:je,ke+1,3)
    case (SLIP_WALL, NO_SLIP_WALL, SYMMETRY)
       wlz(is:ie,js:je,ke+1) = ZERO
       wrz(is:ie,js:je,ke+1) = ZERO
    case (OUTLET)
       wlz(is:ie,js:je,ke+1) = max(wlz(is:ie,js:je,ke+1),ZERO)
       wrz(is:ie,js:je,ke+1) = wlz(is:ie,js:je,ke+1)
    case (INTERIOR, PERIODIC) 
    case  default
       call bl_error("mkutrans_3d: invalid boundary type phys_bc(3,2)")
    end select
    
    if (spherical .eq. 1) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke+1
          do j=js,je
             do i=is,ie
                ! solve Riemann problem using full velocity
                uavg = HALF*(wlz(i,j,k)+wrz(i,j,k))
                test = ((wlz(i,j,k)+w0macz(i,j,k) .le. ZERO .and. &
                     wrz(i,j,k)+w0macz(i,j,k) .ge. ZERO) .or. &
                     (abs(wlz(i,j,k)+wrz(i,j,k)+TWO*w0macz(i,j,k)) .lt. rel_eps))
                wtrans(i,j,k) = merge(wlz(i,j,k),wrz(i,j,k),uavg+w0macz(i,j,k) .gt. ZERO)
                wtrans(i,j,k) = merge(ZERO,wtrans(i,j,k),test)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    else
       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke+1
          do j=js,je
             do i=is,ie
                ! solve Riemann problem using full velocity
                uavg = HALF*(wlz(i,j,k)+wrz(i,j,k))
                test = ((wlz(i,j,k)+w0(k).le.ZERO .and. wrz(i,j,k)+w0(k).ge.ZERO) .or. &
                     (abs(wlz(i,j,k)+wrz(i,j,k)+TWO*w0(k)) .lt. rel_eps))
                wtrans(i,j,k) = merge(wlz(i,j,k),wrz(i,j,k),uavg+w0(k) .gt. ZERO)
                wtrans(i,j,k) = merge(ZERO,wtrans(i,j,k),test)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    end if

    deallocate(wlz,wrz)

  end subroutine mkutrans_3d
  
end module mkutrans_module
