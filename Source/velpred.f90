module velpred_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: velpred

contains

  subroutine velpred(u,umac,utrans,force,normal,w0,w0mac,dx,dt,the_bc_level,mla)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: nr_fine, dr, spherical, dm, nlevs
    use variables, only: foextrap_comp
    use fill_3d_module
    use multifab_physbc_module
    use ml_restriction_module, only : ml_edge_restriction_c

    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: utrans(:,:),force(:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(in   ) :: mla

    integer                  :: i,r,n
    integer                  :: ng_u,ng_um,ng_ut,ng_f,ng_w0,ng_gw,ng_n
    integer                  :: lo(dm), hi(dm)
    real(kind=dp_t), pointer :: uop(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: utp(:,:,:,:)
    real(kind=dp_t), pointer :: vtp(:,:,:,:)
    real(kind=dp_t), pointer :: wtp(:,:,:,:)
    real(kind=dp_t), pointer :: w0xp(:,:,:,:)
    real(kind=dp_t), pointer :: w0yp(:,:,:,:)
    real(kind=dp_t), pointer :: w0zp(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)
    real(kind=dp_t), pointer :: gw0p(:,:,:,:)

    real(kind=dp_t), allocatable :: gradw0_rad(:)
    type(multifab) :: gradw0_cart

    type(bl_prof_timer), save :: bpt

    call build(bpt, "velpred")

    ng_u = u(1)%ng
    ng_um = umac(1,1)%ng
    ng_ut = utrans(1,1)%ng
    ng_f = force(1)%ng
    ng_w0 = w0mac(1,1)%ng
    ng_n = normal(1)%ng

    if (spherical .eq. 1) then
       allocate(gradw0_rad(0:nr_fine-1))
       ! NOTE: here we are doing the computation at the finest level
       do r=0,nr_fine-1
          gradw0_rad(r) = (w0(nlevs,r+1) - w0(nlevs,r)) / dr(nlevs)
       enddo
    endif

    do n=1,nlevs

       call multifab_build(gradw0_cart, u(n)%la,1,1)

       if (spherical .eq. 1) then

          do i = 1, gradw0_cart%nboxes
             if ( multifab_remote(u(n),i) ) cycle
             gw0p => dataptr(gradw0_cart, i)
             lo = lwb(get_box(gradw0_cart,i))
             hi = upb(get_box(gradw0_cart,i))

             call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,gradw0_rad,gw0p, &
                                               lo,hi,dx(n,:),gradw0_cart%ng,0)

          enddo


          ! fill ghost cells for two adjacent grids at the same level
          ! this includes periodic domain boundary ghost cells 
          call multifab_fill_boundary(gradw0_cart)

          ! fill non-periodic domain boundary ghost cells.
          ! NOTE: not sure what the BC should be for gradw0_cart.  Right
          ! now I am just using foextrap_comp.
          call multifab_physbc(gradw0_cart,1,foextrap_comp,1,the_bc_level(n))

       else
          call setval(gradw0_cart, ZERO, all=.true.)
       endif

       ng_gw = gradw0_cart%ng

       do i = 1, u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          uop  => dataptr(u(n),i)
          ump  => dataptr(umac(n,1),i)
          vmp  => dataptr(umac(n,2),i)
          utp  => dataptr(utrans(n,1),i)
          vtp  => dataptr(utrans(n,2),i)
          fp   => dataptr(force(n),i)
          lo   =  lwb(get_box(u(n),i))
          hi   =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call velpred_2d(n, uop(:,:,1,:), ng_u, &
                             utp(:,:,1,1), vtp(:,:,1,1), ng_ut, &
                             ump(:,:,1,1), vmp(:,:,1,1), ng_um, &
                             fp(:,:,1,:), ng_f, w0(n,:), lo, hi, dx(n,:), dt, &
                             the_bc_level(n)%phys_bc_level_array(i,:,:), &
                             the_bc_level(n)%adv_bc_level_array(i,:,:,:))

          case (3)
             wmp  => dataptr(  umac(n,3),i)
             wtp  => dataptr(utrans(n,3),i)
             w0xp  => dataptr(w0mac(n,1),i)
             w0yp  => dataptr(w0mac(n,2),i)
             w0zp  => dataptr(w0mac(n,3),i)
             gw0p => dataptr(gradw0_cart,i)
             np => dataptr(normal(n),i)
             call velpred_3d(n, uop(:,:,:,:), ng_u, &
                             ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                             utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), ng_ut, &
                             fp(:,:,:,:), ng_f, np(:,:,:,:), ng_n, &
                             w0(n,:),w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1), &
                             ng_w0, gw0p(:,:,:,1), ng_gw, lo, hi, dx(n,:), dt, &
                             the_bc_level(n)%phys_bc_level_array(i,:,:), &
                             the_bc_level(n)%adv_bc_level_array(i,:,:,:))
          end select
       end do

       call destroy(gradw0_cart)

    end do

    do n = nlevs,2,-1
       do i = 1, dm
          call ml_edge_restriction_c(umac(n-1,i),1,umac(n,i),1,mla%mba%rr(n-1,:),i,1)
       enddo
    enddo

    if (spherical .eq. 1) then
       deallocate(gradw0_rad)
    endif

    call destroy(bpt)

  end subroutine velpred

  subroutine velpred_2d(n,u,ng_u,utrans,vtrans,ng_ut,umac,vmac,ng_um,force,ng_f, &
                        w0,lo,hi,dx,dt,phys_bc,adv_bc)

    use geometry, only: nr
    use bc_module
    use slope_module
    use bl_constants_module
    use variables, only: rel_eps

    integer        , intent(in   ) :: n,lo(:),hi(:),ng_u,ng_um,ng_ut,ng_f
    real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u :,lo(2)-ng_u :,:)
    real(kind=dp_t), intent(in   ) :: utrans(lo(1)-ng_ut:,lo(2)-ng_ut:)
    real(kind=dp_t), intent(in   ) :: vtrans(lo(1)-ng_ut:,lo(2)-ng_ut:)
    real(kind=dp_t), intent(inout) ::   umac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(inout) ::   vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,:)
    real(kind=dp_t), intent(in   ) ::     w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    integer        , intent(in   ) :: phys_bc(:,:)
    integer        , intent(in   ) :: adv_bc(:,:,:)

    ! Local variables
    real(kind=dp_t) :: slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2)
    real(kind=dp_t) :: slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2)

    ! these correspond to u_L^x, etc.
    real(kind=dp_t), allocatable :: ulx(:,:,:),urx(:,:,:),uimhx(:,:,:)
    real(kind=dp_t), allocatable :: uly(:,:,:),ury(:,:,:),uimhy(:,:,:)

    ! these correspond to umac_L, etc.
    real(kind=dp_t), allocatable :: umacl(:,:),umacr(:,:)
    real(kind=dp_t), allocatable :: vmacl(:,:),vmacr(:,:)

    real(kind=dp_t) :: hx, hy, dt2, dt4, uavg, abs_eps, umax, vlo, vhi

    integer :: i,j,is,js,ie,je

    logical :: test

    allocate(  ulx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(  urx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(uimhx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,2))

    allocate(  uly(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))
    allocate(  ury(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))
    allocate(uimhy(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,2))

    allocate(umacl(lo(1):hi(1)+1,lo(2):hi(2)))
    allocate(umacr(lo(1):hi(1)+1,lo(2):hi(2)))

    allocate(vmacl(lo(1):hi(1),lo(2):hi(2)+1))
    allocate(vmacr(lo(1):hi(1),lo(2):hi(2)+1))

    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)

    abs_eps = 1.d-8
    
    ! Compute rel_eps, which is relative to the max velocity
    umax = abs(u(is,js,1))
    do j = js,je; do i = is,ie
       umax = max(umax,abs(u(i,j,1)))
    end do; end do
    do j = js,je; do i = is,ie
       umax = max(umax,abs(u(i,j,2)+HALF*(w0(j)+w0(j+1))))
    end do; end do
    if (umax .eq. 0.d0) then
       rel_eps = abs_eps
    else
       rel_eps = abs_eps * umax
    endif

    dt2 = HALF*dt
    dt4 = dt/4.0d0

    hx = dx(1)
    hy = dx(2)

    ! compute velocity slopes
    call slopex_2d(u,slopex,lo,hi,ng_u,2,adv_bc)
    call slopey_2d(u,slopey,lo,hi,ng_u,2,adv_bc)

    !******************************************************************
    ! Create u_{\i-\half\e_x}^x, etc.
    !******************************************************************

    do j=js-1,je+1
       do i=is,ie+1
          ! extrapolate both components of velocity to left face
          ulx(i,j,1) = u(i-1,j,1) + (HALF - dt2*max(ZERO,u(i-1,j,1)/hx))*slopex(i-1,j,1)
          ulx(i,j,2) = u(i-1,j,2) + (HALF - dt2*max(ZERO,u(i-1,j,1)/hx))*slopex(i-1,j,2)

          ! extrapolate both components of velocity to right face
          urx(i,j,1) = u(i  ,j,1) - (HALF + dt2*min(ZERO,u(i  ,j,1)/hx))*slopex(i  ,j,1)
          urx(i,j,2) = u(i  ,j,2) - (HALF + dt2*min(ZERO,u(i  ,j,1)/hx))*slopex(i  ,j,2)

          ! impose lo side bc's
          if(i .eq. is) then
             ulx(i,j,1) = merge(u(is-1,j,1),ulx(i,j,1),phys_bc(1,1) .eq. INLET)
             urx(i,j,1) = merge(u(is-1,j,1),urx(i,j,1),phys_bc(1,1) .eq. INLET)
             ulx(i,j,2) = merge(u(is-1,j,2),ulx(i,j,2),phys_bc(1,1) .eq. INLET)
             urx(i,j,2) = merge(u(is-1,j,2),urx(i,j,2),phys_bc(1,1) .eq. INLET)
             if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                ulx(i,j,1) = ZERO
                urx(i,j,1) = ZERO
                ulx(i,j,2) = merge(ZERO,urx(i,j,2),phys_bc(1,1) .eq. NO_SLIP_WALL)
                urx(i,j,2) = merge(ZERO,urx(i,j,2),phys_bc(1,1) .eq. NO_SLIP_WALL)
             endif
          endif

          ! impose hi side bc's
          if(i .eq. ie+1) then
             ulx(i,j,1) = merge(u(ie+1,j,1),ulx(i,j,1),phys_bc(1,2) .eq. INLET)
             urx(i,j,1) = merge(u(ie+1,j,1),urx(i,j,1),phys_bc(1,2) .eq. INLET)
             ulx(i,j,2) = merge(u(ie+1,j,2),ulx(i,j,2),phys_bc(1,2) .eq. INLET)
             urx(i,j,2) = merge(u(ie+1,j,2),urx(i,j,2),phys_bc(1,2) .eq. INLET)
             if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                ulx(i,j,1) = ZERO
                urx(i,j,1) = ZERO
                ulx(i,j,2) = merge(ZERO,ulx(i,j,2),phys_bc(1,2) .eq. NO_SLIP_WALL)
                urx(i,j,2) = merge(ZERO,ulx(i,j,2),phys_bc(1,2) .eq. NO_SLIP_WALL)
             endif
          endif

          ! make normal component of uimhx by first solving a normal Riemann problem
          ! uimhx(1) = { ulx(1) if uavg > 0,
          !              urx(1) if uavg < 0,
          !              0   if (ulx(1) < 0 and urx(1) > 0) or |ulx(1)+urx(1)|<rel_eps }
          uavg = HALF*(ulx(i,j,1)+urx(i,j,1))
          test = ((ulx(i,j,1) .le. ZERO .and. urx(i,j,1) .ge. ZERO) .or. &
               (abs(ulx(i,j,1)+urx(i,j,1)) .lt. rel_eps))
          uimhx(i,j,1) = merge(ulx(i,j,1),urx(i,j,1),uavg .gt. ZERO)
          uimhx(i,j,1) = merge(ZERO,uimhx(i,j,1),test)

          ! now upwind to get transverse component of uimhx
          ! uimhx(2) = { ulx(2) if uimhx(1) > rel_eps
          !              urx(2) if uimhx(1) < rel_eps
          !              uavg   if |uimhx(1)| < rel_eps }
          uimhx(i,j,2) = merge(ulx(i,j,2),urx(i,j,2),uimhx(i,j,1).gt.ZERO)
          uavg = HALF*(ulx(i,j,2)+urx(i,j,2))
          uimhx(i,j,2) = merge(uavg,uimhx(i,j,2),abs(uimhx(i,j,1)).lt.rel_eps)
       enddo
    enddo

    do j=js,je+1
       do i=is-1,ie+1
          if (j .eq. 0) then
             vlo = u(i,j-1,2) + w0(j)
             vhi = u(i,j  ,2) + HALF*(w0(j)+w0(j+1))
          else if (j .eq. nr(n)) then
             vlo = u(i,j-1,2) + HALF*(w0(j-1)+w0(j))
             vhi = u(i,j  ,2) + w0(j)
          else
             vlo = u(i,j-1,2) + HALF*(w0(j-1)+w0(j))
             vhi = u(i,j  ,2) + HALF*(w0(j)+w0(j+1))
          end if

          ! extrapolate both components of velocity to left face
          uly(i,j,1) = u(i,j-1,1) + (HALF - dt2*max(ZERO,vlo/hy))*slopey(i,j-1,1)
          uly(i,j,2) = u(i,j-1,2) + (HALF - dt2*max(ZERO,vlo/hy))*slopey(i,j-1,2)

          ! extrapolate both components of velocity to right face
          ury(i,j,1) = u(i,j  ,1) - (HALF + dt2*min(ZERO,vhi/hy))*slopey(i,j  ,1)
          ury(i,j,2) = u(i,j  ,2) - (HALF + dt2*min(ZERO,vhi/hy))*slopey(i,j  ,2)

          ! impose lo side bc's
          if(j .eq. js) then
             uly(i,j,1) = merge(u(i,js-1,1),uly(i,j,1),phys_bc(2,1) .eq. INLET)
             ury(i,j,1) = merge(u(i,js-1,1),ury(i,j,1),phys_bc(2,1) .eq. INLET)
             uly(i,j,2) = merge(u(i,js-1,2),uly(i,j,2),phys_bc(2,1) .eq. INLET)
             ury(i,j,2) = merge(u(i,js-1,2),ury(i,j,2),phys_bc(2,1) .eq. INLET)
             if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                uly(i,j,1) = merge(ZERO,ury(i,j,1),phys_bc(2,1) .eq. NO_SLIP_WALL)
                ury(i,j,1) = merge(ZERO,ury(i,j,1),phys_bc(2,1) .eq. NO_SLIP_WALL)
                uly(i,j,2) = ZERO
                ury(i,j,2) = ZERO
             endif
          endif

          ! impose hi side bc's
          if(j .eq. je+1) then
             uly(i,j,1) = merge(u(i,je+1,1),uly(i,j,1),phys_bc(2,2) .eq. INLET)
             ury(i,j,1) = merge(u(i,je+1,1),ury(i,j,1),phys_bc(2,2) .eq. INLET)
             uly(i,j,2) = merge(u(i,je+1,2),uly(i,j,2),phys_bc(2,2) .eq. INLET)
             ury(i,j,2) = merge(u(i,je+1,2),ury(i,j,2),phys_bc(2,2) .eq. INLET)
             if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                uly(i,j,1) = merge(ZERO,uly(i,j,1),phys_bc(2,2) .eq. NO_SLIP_WALL)
                ury(i,j,1) = merge(ZERO,uly(i,j,1),phys_bc(2,2) .eq. NO_SLIP_WALL)
                uly(i,j,2) = ZERO
                ury(i,j,2) = ZERO
             endif
          endif

          ! make normal component of uimhy by first solving a normal Riemann problem
          ! uimhy(2) = { uly(2) if uavg > 0,
          !              ury(2) if uavg < 0,
          !              0   if (uly(2) < 0 and ury(2) > 0) or |uly(2)+ury(2)|<rel_eps }
          uavg = HALF*(uly(i,j,2)+ury(i,j,2))
          test = ((uly(i,j,2)+w0(j) .le. ZERO .and. ury(i,j,2)+w0(j) .ge. ZERO) .or. &
               (abs(uly(i,j,2)+ury(i,j,2)+TWO*w0(j)) .lt. rel_eps))
          uimhy(i,j,2) = merge(uly(i,j,2),ury(i,j,2),uavg+w0(j) .gt. ZERO)
          uimhy(i,j,2) = merge(ZERO,uimhy(i,j,2),test)

          ! now upwind to get transverse component of uimhy
          ! uimhy(1) = { uly(1) if uimhy(2) > rel_eps
          !              ury(1) if uimhy(2) < rel_eps
          !              uavg   if |uimhx(2)| < rel_eps }
          uimhy(i,j,1) = merge(uly(i,j,1),ury(i,j,1),uimhy(i,j,2)+w0(j).gt.ZERO)
          uavg = HALF*(uly(i,j,1)+ury(i,j,1))
          uimhy(i,j,1) = merge(uavg,uimhy(i,j,1),abs(uimhy(i,j,2)+w0(j)).lt.rel_eps)
       enddo
    enddo

    !******************************************************************
    ! Create umac and vmac
    !******************************************************************

    do j=js,je
       do i=is,ie+1
          ! extrapolate to edges
          umacl(i,j) = ulx(i,j,1) &
               - (dt4/hy)*(uimhy(i-1,j+1,2)+w0(j+1)+uimhy(i-1,j,2)+w0(j)) &
               * (uimhy(i-1,j+1,1)-uimhy(i-1,j,1)) + dt2*force(i-1,j,1)
          umacr(i,j) = urx(i,j,1) &
               - (dt4/hy)*(uimhy(i  ,j+1,2)+w0(j+1)+uimhy(i  ,j,2)+w0(j)) &
               * (uimhy(i  ,j+1,1)-uimhy(i  ,j,1)) + dt2*force(i  ,j,1)

          ! solve Riemann problem
          uavg = HALF*(umacl(i,j)+umacr(i,j))
          test = ((umacl(i,j) .le. ZERO .and. umacr(i,j) .ge. ZERO) .or. &
               (abs(umacl(i,j)+umacr(i,j)) .lt. rel_eps))
          umac(i,j) = merge(umacl(i,j),umacr(i,j),uavg .gt. ZERO)
          umac(i,j) = merge(uavg,umac(i,j),test) ! varden uses ZERO instead of uavg
       enddo
    enddo

    ! Apply boundary conditions
    do j=js,je
       ! lo side
       if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
          umac(is,j) = ZERO
       elseif (phys_bc(1,1) .eq. INLET) then
          umac(is,j) = u(is-1,j,1)
       elseif (phys_bc(1,1) .eq. OUTLET) then
          umac(is,j) = min(umacr(is,j),ZERO)
       endif

       ! hi side
       if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
          umac(ie+1,j) = ZERO
       elseif (phys_bc(1,2) .eq. INLET) then
          umac(ie+1,j) = u(ie+1,j,1)
       elseif (phys_bc(1,2) .eq. OUTLET) then
          umac(ie+1,j) = max(umacl(ie+1,j),ZERO)
       endif
    enddo

    do j=js,je+1
       do i=is,ie
          ! extrapolate to edges
          vmacl(i,j) = uly(i,j,2) &
               - (dt4/hx)*(uimhx(i+1,j-1,1)+uimhx(i,j-1,1)) &
               * (uimhx(i+1,j-1,2)-uimhx(i,j-1,2)) + dt2*force(i,j-1,2)
          vmacr(i,j) = ury(i,j,2) &
               - (dt4/hx)*(uimhx(i+1,j  ,1)+uimhx(i,j  ,1)) &
               * (uimhx(i+1,j  ,2)-uimhx(i,j  ,2)) + dt2*force(i,j  ,2)

          ! add the (Utilde . e_r) d w_0 /dr e_r term here
          if (j .eq. 0) then
             ! vmacl unchanged since dw_0 / dr = 0
             vmacr(i,j) = vmacr(i,j)-(dt4/hy)*(uimhy(i,j+1,2)+uimhy(i,j  ,2))*(w0(j+1)-w0(j))
          else if (j .eq. nr(n)) then
             vmacl(i,j) = vmacl(i,j)-(dt4/hy)*(uimhy(i,j  ,2)+uimhy(i,j-1,2))*(w0(j)-w0(j-1))
             ! vmacr unchanged since dw_0 / dr = 0
          else
             vmacl(i,j) = vmacl(i,j)-(dt4/hy)*(uimhy(i,j  ,2)+uimhy(i,j-1,2))*(w0(j)-w0(j-1))
             vmacr(i,j) = vmacr(i,j)-(dt4/hy)*(uimhy(i,j+1,2)+uimhy(i,j  ,2))*(w0(j+1)-w0(j))
          end if

          ! solve Riemann problem
          uavg = HALF*(vmacl(i,j)+vmacr(i,j))
          test = ((vmacl(i,j)+w0(j) .le. ZERO .and. vmacr(i,j)+w0(j) .ge. ZERO) .or. &
               (abs(vmacl(i,j)+vmacr(i,j)+TWO*w0(j)) .lt. rel_eps))
          vmac(i,j) = merge(vmacl(i,j),vmacr(i,j),uavg+w0(j) .gt. ZERO)
          vmac(i,j) = merge(uavg,vmac(i,j),test) ! varden uses ZERO instead of uavg
       enddo
    enddo

    ! Apply boundary conditions
    do i=is,ie
       ! lo side
       if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
          vmac(i,js) = ZERO
       elseif (phys_bc(2,1) .eq. INLET) then
          vmac(i,js) = u(i,js-1,2)
       elseif (phys_bc(2,1) .eq. OUTLET) then
          vmac(i,js) = min(vmacr(i,js),ZERO)
       endif

       ! hi side
       if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
          vmac(i,je+1) = ZERO
       elseif (phys_bc(2,2) .eq. INLET) then
          vmac(i,je+1) = u(i,je+1,2)
       elseif (phys_bc(2,2) .eq. OUTLET) then
          vmac(i,je+1) = max(vmacl(i,je+1),ZERO)
       endif
    enddo

    deallocate(ulx,urx,uimhx,uly,ury,uimhy,umacl,umacr,vmacl,vmacr)

  end subroutine velpred_2d

  subroutine velpred_3d(n,u,ng_u, &
                        umac,vmac,wmac,ng_um,utrans,vtrans,wtrans,ng_ut, &
                        force,ng_f,normal,ng_n,w0,w0macx,w0macy,w0macz,ng_w0, &
                        gradw0_cart,ng_gw,lo,hi,dx,dt,phys_bc,adv_bc)

    use bc_module
    use slope_module
    use bl_constants_module
    use geometry, only: spherical, nr
    use variables, only: rel_eps

    integer        , intent(in   ) :: n,lo(:),hi(:)
    integer        , intent(in   ) :: ng_u,ng_um,ng_ut,ng_f,ng_n,ng_w0,ng_gw
    real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :,:)
    real(kind=dp_t), intent(inout) ::   umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(inout) ::   vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(inout) ::   wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) :: utrans(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(in   ) :: vtrans(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(in   ) :: wtrans(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    real(kind=dp_t), intent(in   ) :: normal(lo(1)-ng_n :,lo(2)-ng_n :,lo(3)-ng_n :,:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)    
    real(kind=dp_t), intent(in   ) :: gradw0_cart(lo(1)-ng_gw:,lo(2)-ng_gw:,lo(3)-ng_gw:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    integer        , intent(in   ) :: phys_bc(:,:)
    integer        , intent(in   ) :: adv_bc(:,:,:)

    ! local variables only needed to non-corner coupling code
    real(kind=dp_t) :: s_l(lo(1)-1:hi(1)+2)
    real(kind=dp_t) :: s_r(lo(1)-1:hi(1)+2)
    real(kind=dp_t) :: s_b(lo(2)-1:hi(2)+2)
    real(kind=dp_t) :: s_t(lo(2)-1:hi(2)+2)
    real(kind=dp_t) :: s_u(lo(3)-1:hi(3)+2)
    real(kind=dp_t) :: s_d(lo(3)-1:hi(3)+2)
    real(kind=dp_t) :: ubardt2, vbardt2, wbardt2
    real(kind=dp_t) :: splus, sminus
    real(kind=dp_t) :: savg,st
    real(kind=dp_t) :: sptop,spbot,smtop,smbot,splft,sprgt,smlft,smrgt

    ! local variables
    real(kind=dp_t) :: slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3)
    real(kind=dp_t) :: slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3)
    real(kind=dp_t) :: slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3)

    ! these correspond to u_L^x, etc.
    real(kind=dp_t), allocatable:: ulx(:,:,:,:),urx(:,:,:,:),uimhx(:,:,:,:)
    real(kind=dp_t), allocatable:: uly(:,:,:,:),ury(:,:,:,:),uimhy(:,:,:,:)
    real(kind=dp_t), allocatable:: ulz(:,:,:,:),urz(:,:,:,:),uimhz(:,:,:,:)

    ! these correspond to u_L^{y|z}, etc.
    real(kind=dp_t), allocatable:: ulyz(:,:,:)
    real(kind=dp_t), allocatable:: uryz(:,:,:)
    real(kind=dp_t), allocatable:: uimhyz(:,:,:)

    real(kind=dp_t), allocatable:: ulzy(:,:,:)
    real(kind=dp_t), allocatable:: urzy(:,:,:)
    real(kind=dp_t), allocatable:: uimhzy(:,:,:)

    real(kind=dp_t), allocatable:: vlxz(:,:,:)
    real(kind=dp_t), allocatable:: vrxz(:,:,:)
    real(kind=dp_t), allocatable:: vimhxz(:,:,:)

    real(kind=dp_t), allocatable:: vlzx(:,:,:)
    real(kind=dp_t), allocatable:: vrzx(:,:,:)
    real(kind=dp_t), allocatable:: vimhzx(:,:,:)

    real(kind=dp_t), allocatable:: wlxy(:,:,:)
    real(kind=dp_t), allocatable:: wrxy(:,:,:)
    real(kind=dp_t), allocatable:: wimhxy(:,:,:)

    real(kind=dp_t), allocatable:: wlyx(:,:,:)
    real(kind=dp_t), allocatable:: wryx(:,:,:)
    real(kind=dp_t), allocatable:: wimhyx(:,:,:)

    ! these correspond to umac_L, etc.
    real(kind=dp_t), allocatable:: umacl(:,:,:),umacr(:,:,:)
    real(kind=dp_t), allocatable:: vmacl(:,:,:),vmacr(:,:,:)
    real(kind=dp_t), allocatable:: wmacl(:,:,:),wmacr(:,:,:)

    real(kind=dp_t) :: hx, hy, hz, dt2, dt4, dt6, uavg, abs_eps, umax
    real(kind=dp_t) :: ulo, uhi, vlo, vhi, wlo, whi, Ut_dot_er

    integer :: i,j,k,is,js,ie,je,ks,ke

    logical :: test

    ! normal predictor states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    allocate(ulx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(urx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(uimhx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

    allocate(uly(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(ury(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(uimhy(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1,3))

    allocate(ulz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1,3))
    allocate(urz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1,3))
    allocate(uimhz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1,3))

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    allocate(ulyz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(uryz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(uimhyz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))

    allocate(ulzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))
    allocate(urzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))
    allocate(uimhzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))

    allocate(vlxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))
    allocate(vrxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))
    allocate(vimhxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))

    allocate(vlzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(vrzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(vimhzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))

    allocate(wlxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))
    allocate(wrxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))
    allocate(wimhxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))

    allocate(wlyx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(wryx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(wimhyx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))

    ! mac states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo:hi in the transverse direction
    allocate(umacl(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)))
    allocate(umacr(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)))
    allocate(vmacl(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(vmacr(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(wmacl(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1))
    allocate(wmacr(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1))

    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)
    ks = lo(3)
    ke = hi(3)

    abs_eps = 1.d-8

    ! Compute rel_eps, which is relative to the max velocity
    if (spherical .eq. 1) then
       umax = abs(u(is,js,ks,1))
       do k = ks,ke; do j = js,je; do i = is,ie
          umax = max(umax,abs(u(i,j,k,1)+HALF*(w0macx(i,j,k)+w0macx(i+1,j,k))))
       end do; end do; end do
       do k = ks,ke; do j = js,je; do i = is,ie
          umax = max(umax,abs(u(i,j,k,2)+HALF*(w0macy(i,j,k)+w0macy(i,j+1,k))))
       end do; end do; end do
       do k = ks,ke; do j = js,je; do i = is,ie
          umax = max(umax,abs(u(i,j,k,3)+HALF*(w0macz(i,j,k)+w0macz(i,j,k+1))))
       end do; end do; end do
    else
       umax = abs(u(is,js,ks,1))
       do k = ks,ke; do j = js,je; do i = is,ie
          umax = max(umax,abs(u(i,j,k,1)))
       end do; end do; end do
       do k = ks,ke; do j = js,je; do i = is,ie
          umax = max(umax,abs(u(i,j,k,2)))
       end do; end do; end do
       do k = ks,ke; do j = js,je; do i = is,ie
          umax = max(umax,abs(u(i,j,k,3)+HALF*(w0(k)+w0(k+1))))
       end do; end do; end do
    end if
    if (umax .eq. 0.d0) then
       rel_eps = abs_eps
    else
       rel_eps = abs_eps * umax
    endif

    dt2 = HALF*dt
    dt4 = dt/4.0d0
    dt6 = dt/6.0d0

    hx = dx(1)
    hy = dx(2)
    hz = dx(3)

    do k = lo(3)-1,hi(3)+1
       call slopex_2d(u(:,:,k,:),slopex(:,:,k,:),lo,hi,ng_u,3,adv_bc)
       call slopey_2d(u(:,:,k,:),slopey(:,:,k,:),lo,hi,ng_u,3,adv_bc)
    end do
    call slopez_3d(u,slopez,lo,hi,ng_u,3,adv_bc)

    ! CORNER-COUPLING CODE (NOT FINISHED YET)
    if (.false.) then

    !******************************************************************
    ! Create u_{\i-\half\e_x}^x, etc.
    !******************************************************************

    do k=ks-1,ke+1
       do j=js-1,je+1
          do i=is,ie+1
             if (spherical .eq. 1) then
                ulo = u(i-1,j,k,1) + HALF * (w0macx(i-1,j,k)+w0macx(i  ,j,k))
                uhi = u(i  ,j,k,1) + HALF * (w0macx(i  ,j,k)+w0macx(i+1,j,k))
             else
                ulo = u(i-1,j,k,1)
                uhi = u(i  ,j,k,1)
             end if
             
             ! extrapolate all components of velocity to left face
             ulx(i,j,k,1) = u(i-1,j,k,1) + (HALF - dt2*max(ZERO,ulo)/hx)*slopex(i-1,j,k,1)
             ulx(i,j,k,2) = u(i-1,j,k,2) + (HALF - dt2*max(ZERO,ulo)/hx)*slopex(i-1,j,k,2)
             ulx(i,j,k,3) = u(i-1,j,k,3) + (HALF - dt2*max(ZERO,ulo)/hx)*slopex(i-1,j,k,3)

             ! extrapolate all components of velocity to right face
             urx(i,j,k,1) = u(i,j,k,1) - (HALF + dt2*min(ZERO,uhi)/hx)*slopex(i,j,k,1)
             urx(i,j,k,2) = u(i,j,k,2) - (HALF + dt2*min(ZERO,uhi)/hx)*slopex(i,j,k,2)
             urx(i,j,k,3) = u(i,j,k,3) - (HALF + dt2*min(ZERO,uhi)/hx)*slopex(i,j,k,3)

             ! impose lo side bc's
             if(i .eq. is) then
                ulx(i,j,k,1) = merge(u(is-1,j,k,1),ulx(i,j,k,1),phys_bc(1,1) .eq. INLET)
                urx(i,j,k,1) = merge(u(is-1,j,k,1),urx(i,j,k,1),phys_bc(1,1) .eq. INLET)
                ulx(i,j,k,2) = merge(u(is-1,j,k,2),ulx(i,j,k,2),phys_bc(1,1) .eq. INLET)
                urx(i,j,k,2) = merge(u(is-1,j,k,2),urx(i,j,k,2),phys_bc(1,1) .eq. INLET)
                ulx(i,j,k,3) = merge(u(is-1,j,k,3),ulx(i,j,k,3),phys_bc(1,1) .eq. INLET)
                urx(i,j,k,3) = merge(u(is-1,j,k,3),urx(i,j,k,3),phys_bc(1,1) .eq. INLET)
                if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                   ulx(i,j,k,1) = ZERO
                   urx(i,j,k,1) = ZERO
                   ulx(i,j,k,2) = merge(ZERO,urx(i,j,k,2),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   urx(i,j,k,2) = merge(ZERO,urx(i,j,k,2),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   ulx(i,j,k,3) = merge(ZERO,urx(i,j,k,3),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   urx(i,j,k,3) = merge(ZERO,urx(i,j,k,3),phys_bc(1,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(i .eq. ie+1) then
                ulx(i,j,k,1) = merge(u(ie+1,j,k,1),ulx(i,j,k,1),phys_bc(1,2) .eq. INLET)
                urx(i,j,k,1) = merge(u(ie+1,j,k,1),urx(i,j,k,1),phys_bc(1,2) .eq. INLET)
                ulx(i,j,k,2) = merge(u(ie+1,j,k,2),ulx(i,j,k,2),phys_bc(1,2) .eq. INLET)
                urx(i,j,k,2) = merge(u(ie+1,j,k,2),urx(i,j,k,2),phys_bc(1,2) .eq. INLET)
                ulx(i,j,k,3) = merge(u(ie+1,j,k,3),ulx(i,j,k,3),phys_bc(1,2) .eq. INLET)
                urx(i,j,k,3) = merge(u(ie+1,j,k,3),urx(i,j,k,3),phys_bc(1,2) .eq. INLET)
                if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                   ulx(i,j,k,1) = ZERO
                   urx(i,j,k,1) = ZERO
                   ulx(i,j,k,2) = merge(ZERO,ulx(i,j,k,2),phys_bc(1,2) .eq. NO_SLIP_WALL)
                   urx(i,j,k,2) = merge(ZERO,ulx(i,j,k,2),phys_bc(1,2) .eq. NO_SLIP_WALL)
                   ulx(i,j,k,3) = merge(ZERO,ulx(i,j,k,3),phys_bc(1,2) .eq. NO_SLIP_WALL)
                   urx(i,j,k,3) = merge(ZERO,ulx(i,j,k,3),phys_bc(1,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! make normal component of uimhx by first solving a normal Riemann problem
             uavg = HALF*(ulx(i,j,k,1)+urx(i,j,k,1))
             test = ((ulx(i,j,k,1) .le. ZERO .and. urx(i,j,k,1) .ge. ZERO) .or. &
                  (abs(ulx(i,j,k,1)+urx(i,j,k,1)) .lt. rel_eps))
             uimhx(i,j,k,1) = merge(ulx(i,j,k,1),urx(i,j,k,1),uavg .gt. ZERO)
             uimhx(i,j,k,1) = merge(ZERO,uimhx(i,j,k,1),test)

             ! now upwind to get transverse components of uimhx
             uimhx(i,j,k,2) = merge(ulx(i,j,k,2),urx(i,j,k,2),uimhx(i,j,k,1).gt.ZERO)
             uavg = HALF*(ulx(i,j,k,2)+urx(i,j,k,2))
             uimhx(i,j,k,2) = merge(uavg,uimhx(i,j,k,2),abs(uimhx(i,j,k,1)).lt.rel_eps)

             uimhx(i,j,k,3) = merge(ulx(i,j,k,3),urx(i,j,k,3),uimhx(i,j,k,1).gt.ZERO)
             uavg = HALF*(ulx(i,j,k,3)+urx(i,j,k,3))
             uimhx(i,j,k,3) = merge(uavg,uimhx(i,j,k,3),abs(uimhx(i,j,k,1)).lt.rel_eps)
          enddo
       enddo
    enddo

    do k=ks-1,ke+1
       do j=js,je+1
          do i=is-1,ie+1
             if (spherical .eq. 1) then
                vlo = u(i,j-1,k,2) + HALF * (w0macy(i,j-1,k)+w0macy(i,j  ,k))
                vhi = u(i,j  ,k,2) + HALF * (w0macy(i,j  ,k)+w0macy(i,j+1,k))
             else
                vlo = u(i,j-1,k,2)
                vhi = u(i,j  ,k,2)
             end if

             ! extrapolate all components of velocity to left face
             uly(i,j,k,1) = u(i,j-1,k,1) + (HALF - dt2*max(ZERO,vlo)/hy)*slopey(i,j-1,k,1)
             uly(i,j,k,2) = u(i,j-1,k,2) + (HALF - dt2*max(ZERO,vlo)/hy)*slopey(i,j-1,k,2)
             uly(i,j,k,3) = u(i,j-1,k,3) + (HALF - dt2*max(ZERO,vlo)/hy)*slopey(i,j-1,k,3)

             ! extrapolate all components of velocity to right face
             ury(i,j,k,1) = u(i,j,k,1) - (HALF + dt2*min(ZERO,vhi)/hy)*slopey(i,j,k,1)
             ury(i,j,k,2) = u(i,j,k,2) - (HALF + dt2*min(ZERO,vhi)/hy)*slopey(i,j,k,2)
             ury(i,j,k,3) = u(i,j,k,3) - (HALF + dt2*min(ZERO,vhi)/hy)*slopey(i,j,k,3)

             ! impose lo side bc's
             if(j .eq. js) then
                uly(i,j,k,1) = merge(u(i,js-1,k,1),uly(i,j,k,1),phys_bc(2,1) .eq. INLET)
                ury(i,j,k,1) = merge(u(i,js-1,k,1),ury(i,j,k,1),phys_bc(2,1) .eq. INLET)
                uly(i,j,k,2) = merge(u(i,js-1,k,2),uly(i,j,k,2),phys_bc(2,1) .eq. INLET)
                ury(i,j,k,2) = merge(u(i,js-1,k,2),ury(i,j,k,2),phys_bc(2,1) .eq. INLET)
                uly(i,j,k,3) = merge(u(i,js-1,k,3),uly(i,j,k,3),phys_bc(2,1) .eq. INLET)
                ury(i,j,k,3) = merge(u(i,js-1,k,3),ury(i,j,k,3),phys_bc(2,1) .eq. INLET)
                if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                   uly(i,j,k,1) = merge(ZERO,ury(i,j,k,1),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   ury(i,j,k,1) = merge(ZERO,ury(i,j,k,1),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   uly(i,j,k,2) = ZERO
                   ury(i,j,k,2) = ZERO
                   uly(i,j,k,3) = merge(ZERO,ury(i,j,k,3),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   ury(i,j,k,3) = merge(ZERO,ury(i,j,k,3),phys_bc(2,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(j .eq. je+1) then
                uly(i,j,k,1) = merge(u(i,je+1,k,1),uly(i,j,k,1),phys_bc(2,2) .eq. INLET)
                ury(i,j,k,1) = merge(u(i,je+1,k,1),ury(i,j,k,1),phys_bc(2,2) .eq. INLET)
                uly(i,j,k,2) = merge(u(i,je+1,k,2),uly(i,j,k,2),phys_bc(2,2) .eq. INLET)
                ury(i,j,k,2) = merge(u(i,je+1,k,2),ury(i,j,k,2),phys_bc(2,2) .eq. INLET)
                uly(i,j,k,3) = merge(u(i,je+1,k,3),uly(i,j,k,3),phys_bc(2,2) .eq. INLET)
                ury(i,j,k,3) = merge(u(i,je+1,k,3),ury(i,j,k,3),phys_bc(2,2) .eq. INLET)
                if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                   uly(i,j,k,1) = merge(ZERO,uly(i,j,k,1),phys_bc(2,2) .eq. NO_SLIP_WALL)
                   ury(i,j,k,1) = merge(ZERO,uly(i,j,k,1),phys_bc(2,2) .eq. NO_SLIP_WALL)
                   uly(i,j,k,2) = ZERO
                   ury(i,j,k,2) = ZERO
                   uly(i,j,k,3) = merge(ZERO,uly(i,j,k,3),phys_bc(2,2) .eq. NO_SLIP_WALL)
                   ury(i,j,k,3) = merge(ZERO,uly(i,j,k,3),phys_bc(2,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! make normal component of uimhy by first solving a normal Riemann problem
             uavg = HALF*(uly(i,j,k,2)+ury(i,j,k,2))
             test = ((uly(i,j,k,2) .le. ZERO .and. ury(i,j,k,2) .ge. ZERO) .or. &
                  (abs(uly(i,j,k,2)+ury(i,j,k,2)) .lt. rel_eps))
             uimhy(i,j,k,2) = merge(uly(i,j,k,2),ury(i,j,k,2),uavg .gt. ZERO)
             uimhy(i,j,k,2) = merge(ZERO,uimhy(i,j,k,2),test)

             ! now upwind to get transverse components of uimhy
             uimhy(i,j,k,1) = merge(uly(i,j,k,1),ury(i,j,k,1),uimhy(i,j,k,2).gt.ZERO)
             uavg = HALF*(uly(i,j,k,1)+ury(i,j,k,1))
             uimhy(i,j,k,1) = merge(uavg,uimhy(i,j,k,1),abs(uimhy(i,j,k,2)).lt.rel_eps)

             uimhy(i,j,k,3) = merge(uly(i,j,k,3),ury(i,j,k,3),uimhy(i,j,k,2).gt.ZERO)
             uavg = HALF*(uly(i,j,k,3)+ury(i,j,k,3))
             uimhy(i,j,k,3) = merge(uavg,uimhy(i,j,k,3),abs(uimhy(i,j,k,2)).lt.rel_eps)
          enddo
       enddo
    enddo

    do k=ks,ke+1
       do j=js-1,je+1
          do i=is-1,ie+1
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

             ! extrapolate all components of velocity to left face
             ulz(i,j,k,1) = u(i,j,k-1,1) + (HALF - dt2*max(ZERO,wlo)/hz)*slopez(i,j,k-1,1)
             ulz(i,j,k,2) = u(i,j,k-1,2) + (HALF - dt2*max(ZERO,wlo)/hz)*slopez(i,j,k-1,2)
             ulz(i,j,k,3) = u(i,j,k-1,3) + (HALF - dt2*max(ZERO,wlo)/hz)*slopez(i,j,k-1,3)

             ! extrapolate all components of velocity to right face
             urz(i,j,k,1) = u(i,j,k,1) - (HALF + dt2*min(ZERO,whi)/hz)*slopez(i,j,k,1)
             urz(i,j,k,2) = u(i,j,k,2) - (HALF + dt2*min(ZERO,whi)/hz)*slopez(i,j,k,2)
             urz(i,j,k,3) = u(i,j,k,3) - (HALF + dt2*min(ZERO,whi)/hz)*slopez(i,j,k,3)

             ! impose lo side bc's
             if(k .eq. ks) then
                ulz(i,j,k,1) = merge(u(i,j,ks-1,1),ulz(i,j,k,1),phys_bc(3,1) .eq. INLET)
                urz(i,j,k,1) = merge(u(i,j,ks-1,1),urz(i,j,k,1),phys_bc(3,1) .eq. INLET)
                ulz(i,j,k,2) = merge(u(i,j,ks-1,2),ulz(i,j,k,2),phys_bc(3,1) .eq. INLET)
                urz(i,j,k,2) = merge(u(i,j,ks-1,2),urz(i,j,k,2),phys_bc(3,1) .eq. INLET)
                ulz(i,j,k,3) = merge(u(i,j,ks-1,3),ulz(i,j,k,3),phys_bc(3,1) .eq. INLET)
                urz(i,j,k,3) = merge(u(i,j,ks-1,3),urz(i,j,k,3),phys_bc(3,1) .eq. INLET)
                if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                   ulz(i,j,k,1) = merge(ZERO,urz(i,j,k,1),phys_bc(3,1) .eq. NO_SLIP_WALL)
                   urz(i,j,k,1) = merge(ZERO,urz(i,j,k,1),phys_bc(3,1) .eq. NO_SLIP_WALL)
                   ulz(i,j,k,2) = merge(ZERO,urz(i,j,k,2),phys_bc(3,1) .eq. NO_SLIP_WALL)
                   urz(i,j,k,2) = merge(ZERO,urz(i,j,k,2),phys_bc(3,1) .eq. NO_SLIP_WALL)
                   ulz(i,j,k,3) = ZERO
                   urz(i,j,k,3) = ZERO
                endif
             endif

             ! impose hi side bc's
             if(k .eq. ke+1) then
                ulz(i,j,k,1) = merge(u(i,j,ke+1,1),ulz(i,j,k,1),phys_bc(3,2) .eq. INLET)
                urz(i,j,k,1) = merge(u(i,j,ke+1,1),urz(i,j,k,1),phys_bc(3,2) .eq. INLET)
                ulz(i,j,k,2) = merge(u(i,j,ke+1,2),ulz(i,j,k,2),phys_bc(3,2) .eq. INLET)
                urz(i,j,k,2) = merge(u(i,j,ke+1,2),urz(i,j,k,2),phys_bc(3,2) .eq. INLET)
                ulz(i,j,k,3) = merge(u(i,j,ke+1,3),ulz(i,j,k,3),phys_bc(3,2) .eq. INLET)
                urz(i,j,k,3) = merge(u(i,j,ke+1,3),urz(i,j,k,3),phys_bc(3,2) .eq. INLET)
                if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                   ulz(i,j,k,1) = merge(ZERO,ulz(i,j,k,1),phys_bc(3,2) .eq. NO_SLIP_WALL)
                   urz(i,j,k,1) = merge(ZERO,ulz(i,j,k,1),phys_bc(3,2) .eq. NO_SLIP_WALL)
                   ulz(i,j,k,2) = merge(ZERO,ulz(i,j,k,2),phys_bc(3,2) .eq. NO_SLIP_WALL)
                   urz(i,j,k,2) = merge(ZERO,ulz(i,j,k,2),phys_bc(3,2) .eq. NO_SLIP_WALL)
                   ulz(i,j,k,3) = ZERO
                   urz(i,j,k,3) = ZERO
                endif
             endif

             ! make normal component of uimhz by first solving a normal Riemann problem
             uavg = HALF*(ulz(i,j,k,3)+urz(i,j,k,3))
             test = ((ulz(i,j,k,3) .le. ZERO .and. urz(i,j,k,3) .ge. ZERO) .or. &
                  (abs(ulz(i,j,k,3)+urz(i,j,k,3)) .lt. rel_eps))
             uimhz(i,j,k,3) = merge(ulz(i,j,k,3),urz(i,j,k,3),uavg .gt. ZERO)
             uimhz(i,j,k,3) = merge(ZERO,uimhz(i,j,k,3),test)

             ! now upwind to get transverse components of uimhz
             uimhz(i,j,k,1) = merge(ulz(i,j,k,1),urz(i,j,k,1),uimhz(i,j,k,3).gt.ZERO)
             uavg = HALF*(ulz(i,j,k,1)+urz(i,j,k,1))
             uimhz(i,j,k,1) = merge(uavg,uimhz(i,j,k,1),abs(uimhz(i,j,k,3)).lt.rel_eps)

             uimhz(i,j,k,2) = merge(ulz(i,j,k,2),urz(i,j,k,2),uimhz(i,j,k,3).gt.ZERO)
             uavg = HALF*(ulz(i,j,k,2)+urz(i,j,k,2))
             uimhz(i,j,k,2) = merge(uavg,uimhz(i,j,k,2),abs(uimhz(i,j,k,3)).lt.rel_eps)
          enddo
       enddo
    enddo

    !******************************************************************
    ! Create u_{\i-\half\e_y}^{y|z}, etc.
    !******************************************************************

    ! uimhyz loop
    do k=ks,ke
       do j=js,je+1
          do i=is-1,ie+1
             ! extrapolate to faces
             ulyz(i,j,k) = uly(i,j,k,1) &
                  - (dt6/hz)*(uimhz(i,j-1,k+1,3)+uimhz(i,j-1,k,3))*(uimhz(i,j-1,k+1,1)-uimhz(i,j-1,k,1))
             uryz(i,j,k) = ury(i,j,k,1) &
                  - (dt6/hz)*(uimhz(i,j  ,k+1,3)+uimhz(i,j  ,k,3))*(uimhz(i,j  ,k+1,1)-uimhz(i,j  ,k,1))

             ! impose lo side bc's
             if(j .eq. js) then
                ulyz(i,j,k) = merge(u(i,js-1,k,1),ulyz(i,j,k),phys_bc(2,1) .eq. INLET)
                uryz(i,j,k) = merge(u(i,js-1,k,1),uryz(i,j,k),phys_bc(2,1) .eq. INLET)
                if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                   ulyz(i,j,k) = merge(ZERO,uryz(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   uryz(i,j,k) = merge(ZERO,uryz(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(j .eq. je+1) then
                ulyz(i,j,k) = merge(u(i,je+1,k,1),ulyz(i,j,k),phys_bc(2,2) .eq. INLET)
                uryz(i,j,k) = merge(u(i,je+1,k,1),uryz(i,j,k),phys_bc(2,2) .eq. INLET)
                if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                   ulyz(i,j,k) = merge(ZERO,ulyz(i,j,k),phys_bc(2,2) .eq. NO_SLIP_WALL)
                   uryz(i,j,k) = merge(ZERO,ulyz(i,j,k),phys_bc(2,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! upwind
             uimhyz(i,j,k) = merge(ulyz(i,j,k),uryz(i,j,k),uimhy(i,j,k,2).gt.ZERO)
             uavg = HALF*(ulyz(i,j,k)+uryz(i,j,k))
             uimhyz(i,j,k) = merge(uavg,uimhyz(i,j,k),abs(uimhy(i,j,k,2)).lt.rel_eps)
          enddo
       enddo
    enddo

    ! uimhzy loop
    do k=ks,ke+1
       do j=js,je
          do i=is-1,ie+1
             ! extrapolate to faces
             ulzy(i,j,k) = ulz(i,j,k,1) &
                  - (dt6/hy)*(uimhy(i,j+1,k-1,2)+uimhy(i,j,k-1,2))*(uimhy(i,j+1,k-1,1)-uimhy(i,j,k-1,1))
             urzy(i,j,k) = urz(i,j,k,1) &
                  - (dt6/hy)*(uimhy(i,j+1,k  ,2)+uimhy(i,j,k  ,2))*(uimhy(i,j+1,k  ,1)-uimhy(i,j,k  ,1))

             ! impose lo side bc's
             if(k .eq. ks) then
                ulzy(i,j,k) = merge(u(i,j,ks-1,1),ulzy(i,j,k),phys_bc(3,1) .eq. INLET)
                urzy(i,j,k) = merge(u(i,j,ks-1,1),urzy(i,j,k),phys_bc(3,1) .eq. INLET)
                if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                   ulzy(i,j,k) = merge(ZERO,urzy(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                   urzy(i,j,k) = merge(ZERO,urzy(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(k .eq. ke+1) then
                ulzy(i,j,k) = merge(u(i,j,ke+1,1),ulzy(i,j,k),phys_bc(3,2) .eq. INLET)
                urzy(i,j,k) = merge(u(i,j,ke+1,1),urzy(i,j,k),phys_bc(3,2) .eq. INLET)
                if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                   ulzy(i,j,k) = merge(ZERO,ulzy(i,j,k),phys_bc(3,2) .eq. NO_SLIP_WALL)
                   urzy(i,j,k) = merge(ZERO,ulzy(i,j,k),phys_bc(3,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! upwind
             uimhzy(i,j,k) = merge(ulzy(i,j,k),urzy(i,j,k),uimhz(i,j,k,3).gt.ZERO)
             uavg = HALF*(ulzy(i,j,k)+urzy(i,j,k))
             uimhzy(i,j,k) = merge(uavg,uimhzy(i,j,k),abs(uimhz(i,j,k,3)).lt.rel_eps)
          enddo
       enddo
    enddo

    ! vimhxz loop
    do k=ks,ke
       do j=js-1,je+1
          do i=is,ie+1
             ! extrapolate to faces
             vlxz(i,j,k) = ulx(i,j,k,2) &
                  - (dt6/hz)*(uimhz(i-1,j,k+1,3)+uimhz(i-1,j,k,3))*(uimhz(i-1,j,k+1,2)-uimhz(i-1,j,k,2))
             vrxz(i,j,k) = urx(i,j,k,2) &
                  - (dt6/hz)*(uimhz(i  ,j,k+1,3)+uimhz(i  ,j,k,3))*(uimhz(i  ,j,k+1,2)-uimhz(i  ,j,k,2))

             ! impose lo side bc's
             if(i .eq. is) then
                vlxz(i,j,k) = merge(u(is-1,j,k,2),vlxz(i,j,k),phys_bc(1,1) .eq. INLET)
                vrxz(i,j,k) = merge(u(is-1,j,k,2),vrxz(i,j,k),phys_bc(1,1) .eq. INLET)
                if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                   vlxz(i,j,k) = merge(ZERO,vrxz(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   vrxz(i,j,k) = merge(ZERO,vrxz(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(i .eq. ie+1) then
                vlxz(i,j,k) = merge(u(ie+1,j,k,2),vlxz(i,j,k),phys_bc(1,2) .eq. INLET)
                vrxz(i,j,k) = merge(u(ie+1,j,k,2),vrxz(i,j,k),phys_bc(1,2) .eq. INLET)
                if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                   vlxz(i,j,k) = merge(ZERO,vlxz(i,j,k),phys_bc(1,2) .eq. NO_SLIP_WALL)
                   vrxz(i,j,k) = merge(ZERO,vlxz(i,j,k),phys_bc(1,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! upwind
             vimhxz(i,j,k) = merge(vlxz(i,j,k),vrxz(i,j,k),uimhx(i,j,k,1).gt.ZERO)
             uavg = HALF*(vlxz(i,j,k)+vrxz(i,j,k))
             vimhxz(i,j,k) = merge(uavg,vimhxz(i,j,k),abs(uimhx(i,j,k,1)).lt.rel_eps)
          enddo
       enddo
    enddo

    ! vimhzx loop
    do k=ks,ke+1
       do j=js-1,je+1
          do i=is,ie
             ! extrapolate to faces
             vlzx(i,j,k) = ulz(i,j,k,2) &
                  - (dt6/hx)*(uimhx(i+1,j,k-1,1)+uimhx(i,j,k-1,1))*(uimhx(i+1,j,k-1,2)-uimhx(i,j,k-1,2))
             vrzx(i,j,k) = urz(i,j,k,2) &
                  - (dt6/hx)*(uimhx(i+1,j,k  ,1)+uimhx(i,j,k  ,1))*(uimhx(i+1,j,k  ,2)-uimhx(i,j,k  ,2))

             ! impose lo side bc's
             if(k .eq. ks) then
                vlzx(i,j,k) = merge(u(i,j,ks-1,1),vlzx(i,j,k),phys_bc(3,1) .eq. INLET)
                vrzx(i,j,k) = merge(u(i,j,ks-1,1),vrzx(i,j,k),phys_bc(3,1) .eq. INLET)
                if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                   vlzx(i,j,k) = merge(ZERO,vrzx(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                   vrzx(i,j,k) = merge(ZERO,vrzx(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(k .eq. ke+1) then
                vlzx(i,j,k) = merge(u(i,j,ke+1,1),vlzx(i,j,k),phys_bc(3,2) .eq. INLET)
                vrzx(i,j,k) = merge(u(i,j,ke+1,1),vrzx(i,j,k),phys_bc(3,2) .eq. INLET)
                if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                   vlzx(i,j,k) = merge(ZERO,vlzx(i,j,k),phys_bc(3,2) .eq. NO_SLIP_WALL)
                   vrzx(i,j,k) = merge(ZERO,vlzx(i,j,k),phys_bc(3,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! upwind
             vimhzx(i,j,k) = merge(vlzx(i,j,k),vrzx(i,j,k),uimhz(i,j,k,3).gt.ZERO)
             uavg = HALF*(vlzx(i,j,k)+vrzx(i,j,k))
             vimhzx(i,j,k) = merge(uavg,vimhzx(i,j,k),abs(uimhz(i,j,k,3)).lt.rel_eps)
          enddo
       enddo
    enddo

    ! wimhxy loop
    do k=ks-1,ke+1
       do j=js,je
          do i=is,ie+1
             ! extrapolate to faces
             wlxy(i,j,k) = ulx(i,j,k,3) &
                  - (dt6/hy)*(uimhy(i-1,j+1,k,2)+uimhy(i-1,j,k,2))*(uimhy(i-1,j+1,k,3)-uimhy(i-1,j,k,3))
             wrxy(i,j,k) = urx(i,j,k,3) &
                  - (dt6/hy)*(uimhy(i  ,j+1,k,2)+uimhy(i  ,j,k,2))*(uimhy(i  ,j+1,k,3)-uimhy(i  ,j,k,3))

             ! impose lo side bc's
             if(i .eq. is) then
                wlxy(i,j,k) = merge(u(is-1,j,k,3),wlxy(i,j,k),phys_bc(1,1) .eq. INLET)
                wrxy(i,j,k) = merge(u(is-1,j,k,3),wrxy(i,j,k),phys_bc(1,1) .eq. INLET)
                if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                   wlxy(i,j,k) = merge(ZERO,wrxy(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   wrxy(i,j,k) = merge(ZERO,wrxy(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(i .eq. ie+1) then
                wlxy(i,j,k) = merge(u(ie+1,j,k,3),wlxy(i,j,k),phys_bc(1,2) .eq. INLET)
                wrxy(i,j,k) = merge(u(ie+1,j,k,3),wrxy(i,j,k),phys_bc(1,2) .eq. INLET)
                if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                   wlxy(i,j,k) = merge(ZERO,wlxy(i,j,k),phys_bc(1,2) .eq. NO_SLIP_WALL)
                   wrxy(i,j,k) = merge(ZERO,wlxy(i,j,k),phys_bc(1,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! upwind
             wimhxy(i,j,k) = merge(wlxy(i,j,k),wrxy(i,j,k),uimhx(i,j,k,1).gt.ZERO)
             uavg = HALF*(wlxy(i,j,k)+wrxy(i,j,k))
             wimhxy(i,j,k) = merge(uavg,wimhxy(i,j,k),abs(uimhx(i,j,k,1)).lt.rel_eps)
          enddo
       enddo
    enddo

    ! wimhyx loop
    do k=ks-1,ke+1
       do j=js,je+1
          do i=is,ie
             ! extrapolate to faces
             wlyx(i,j,k) = uly(i,j,k,3) &
                  - (dt6/hx)*(uimhx(i+1,j-1,k,1)+uimhx(i,j-1,k,1))*(uimhx(i+1,j-1,k,3)-uimhx(i,j-1,k,3))
             wryx(i,j,k) = ury(i,j,k,3) &
                  - (dt6/hx)*(uimhx(i+1,j  ,k,1)+uimhx(i,j  ,k,1))*(uimhx(i+1,j  ,k,3)-uimhx(i,j  ,k,3))

             ! impose lo side bc's
             if(j .eq. js) then
                wlyx(i,j,k) = merge(u(i,js-1,k,3),wlyx(i,j,k),phys_bc(2,1) .eq. INLET)
                wryx(i,j,k) = merge(u(i,js-1,k,3),wryx(i,j,k),phys_bc(2,1) .eq. INLET)
                if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                   wlyx(i,j,k) = merge(ZERO,wryx(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   wryx(i,j,k) = merge(ZERO,wryx(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                endif
             endif

             ! impose hi side bc's
             if(j .eq. je+1) then
                wlyx(i,j,k) = merge(u(i,je+1,k,3),wlyx(i,j,k),phys_bc(2,2) .eq. INLET)
                wryx(i,j,k) = merge(u(i,je+1,k,3),wryx(i,j,k),phys_bc(2,2) .eq. INLET)
                if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                   wlyx(i,j,k) = merge(ZERO,wlyx(i,j,k),phys_bc(2,2) .eq. NO_SLIP_WALL)
                   wryx(i,j,k) = merge(ZERO,wlyx(i,j,k),phys_bc(2,2) .eq. NO_SLIP_WALL)
                endif
             endif

             ! upwind
             wimhyx(i,j,k) = merge(wlyx(i,j,k),wryx(i,j,k),uimhy(i,j,k,2).gt.ZERO)
             uavg = HALF*(wlyx(i,j,k)+wryx(i,j,k))
             wimhyx(i,j,k) = merge(uavg,wimhyx(i,j,k),abs(uimhy(i,j,k,2)).lt.rel_eps)
          enddo
       enddo
    enddo


    !******************************************************************
    ! Create umac, etc.
    !******************************************************************

    do k=ks,ke
       do j=js,je
          do i=is,ie+1
             ! extrapolate to edges
             umacl(i,j,k) = ulx(i,j,k,1) &
                  - (dt4/hy)*(uimhy(i-1,j+1,k  ,2)+uimhy(i-1,j,k,2))*(uimhyz(i-1,j+1,k  )-uimhyz(i-1,j,k)) &
                  - (dt4/hz)*(uimhz(i-1,j  ,k+1,3)+uimhz(i-1,j,k,3))*(uimhzy(i-1,j  ,k+1)-uimhzy(i-1,j,k)) &
                  + dt2*force(i-1,j,k,1)
             umacr(i,j,k) = urx(i,j,k,1) &
                  - (dt4/hy)*(uimhy(i  ,j+1,k  ,2)+uimhy(i  ,j,k,2))*(uimhyz(i  ,j+1,k  )-uimhyz(i  ,j,k)) &
                  - (dt4/hz)*(uimhz(i  ,j  ,k+1,3)+uimhz(i  ,j,k,3))*(uimhzy(i  ,j  ,k+1)-uimhzy(i  ,j,k)) &
                  + dt2*force(i  ,j,k,1)

             ! solve Riemann problem
             uavg = HALF*(umacl(i,j,k)+umacr(i,j,k))
             test = ((umacl(i,j,k) .le. ZERO .and. umacr(i,j,k) .ge. ZERO) .or. &
                  (abs(umacl(i,j,k)+umacr(i,j,k)) .lt. rel_eps))
             umac(i,j,k) = merge(umacl(i,j,k),umacr(i,j,k),uavg .gt. ZERO)
             umac(i,j,k) = merge(uavg,umac(i,j,k),test) ! varden uses ZERO instead of uavg
          enddo
       enddo
    enddo

    ! Apply boundary conditions
    do k=ks,ke
       do j=js,je
          ! lo side
          if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
             umac(is,j,k) = ZERO
          elseif (phys_bc(1,1) .eq. INLET) then
             umac(is,j,k) = u(is-1,j,k,1)
          elseif (phys_bc(1,1) .eq. OUTLET) then
             umac(is,j,k) = min(umacr(is,j,k),ZERO)
          endif

          ! hi side
          if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
             umac(ie+1,j,k) = ZERO
          elseif (phys_bc(1,2) .eq. INLET) then
             umac(ie+1,j,k) = u(ie+1,j,k,1)
          elseif (phys_bc(1,2) .eq. OUTLET) then
             umac(ie+1,j,k) = max(umacl(ie+1,j,k),ZERO)
          endif
       enddo
    enddo

    do k=ks,ke
       do j=js,je+1
          do i=is,ie
             ! extrapolate to edges
             vmacl(i,j,k) = uly(i,j,k,2) &
                  - (dt4/hx)*(uimhx(i+1,j-1,k  ,1)+uimhx(i,j-1,k,1))*(vimhxz(i+1,j-1,k  )-vimhxz(i,j-1,k)) &
                  - (dt4/hz)*(uimhz(i  ,j-1,k+1,3)+uimhz(i,j-1,k,3))*(vimhzx(i  ,j-1,k+1)-vimhzx(i,j-1,k)) &
                  + dt2*force(i,j-1,k,2)
             vmacr(i,j,k) = ury(i,j,k,2) &
                  - (dt4/hx)*(uimhx(i+1,j  ,k  ,1)+uimhx(i,j  ,k,1))*(vimhxz(i+1,j  ,k  )-vimhxz(i,j  ,k)) &
                  - (dt4/hz)*(uimhz(i  ,j  ,k+1,3)+uimhz(i,j  ,k,3))*(vimhzx(i  ,j  ,k+1)-vimhzx(i,j  ,k)) &
                  + dt2*force(i,j  ,k,2)

             ! solve Riemann problem
             uavg = HALF*(vmacl(i,j,k)+vmacr(i,j,k))
             test = ((vmacl(i,j,k) .le. ZERO .and. vmacr(i,j,k) .ge. ZERO) .or. &
                  (abs(vmacl(i,j,k)+vmacr(i,j,k)) .lt. rel_eps))
             vmac(i,j,k) = merge(vmacl(i,j,k),vmacr(i,j,k),uavg .gt. ZERO)
             vmac(i,j,k) = merge(uavg,vmac(i,j,k),test) ! varden uses ZERO instead of uavg
          enddo
       enddo
    enddo

    ! Apply boundary conditions
    do k=ks,ke
       do i=is,ie
          ! lo side
          if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
             vmac(i,js,k) = ZERO
          elseif (phys_bc(2,1) .eq. INLET) then
             vmac(i,js,k) = u(i,js-1,k,2)
          elseif (phys_bc(2,1) .eq. OUTLET) then
             vmac(i,js,k) = min(vmacr(i,js,k),ZERO)
          endif

          ! hi side
          if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
             vmac(i,je+1,k) = ZERO
          elseif (phys_bc(2,2) .eq. INLET) then
             vmac(i,je+1,k) = u(i,je+1,k,2)
          elseif (phys_bc(2,2) .eq. OUTLET) then
             vmac(i,je+1,k) = max(vmacl(i,je+1,k),ZERO)
          endif
       enddo
    enddo

    do k=ks,ke+1
       do j=js,je
          do i=is,ie
             ! extrapolate to edges
             wmacl(i,j,k) = ulz(i,j,k,3) &
                  - (dt4/hx)*(uimhx(i+1,j  ,k-1,1)+uimhx(i,j,k-1,1))*(wimhxy(i+1,j  ,k-1)-wimhxy(i,j,k-1)) &
                  - (dt4/hy)*(uimhy(i  ,j+1,k-1,2)+uimhy(i,j,k-1,2))*(wimhyx(i  ,j+1,k-1)-wimhyx(i,j,k-1)) &
                  + dt2*force(i,j,k-1,3)
             wmacr(i,j,k) = urz(i,j,k,3) &
                  - (dt4/hx)*(uimhx(i+1,j  ,k  ,1)+uimhx(i,j,k  ,1))*(wimhxy(i+1,j  ,k  )-wimhxy(i,j,k  )) &
                  - (dt4/hy)*(uimhy(i  ,j+1,k  ,2)+uimhy(i,j,k  ,2))*(wimhyx(i  ,j+1,k  )-wimhyx(i,j,k  )) &
                  + dt2*force(i,j,k  ,3)

             ! solve Riemann problem
             uavg = HALF*(wmacl(i,j,k)+wmacr(i,j,k))
             test = ((wmacl(i,j,k) .le. ZERO .and. wmacr(i,j,k) .ge. ZERO) .or. &
                  (abs(wmacl(i,j,k)+wmacr(i,j,k)) .lt. rel_eps))
             wmac(i,j,k) = merge(wmacl(i,j,k),wmacr(i,j,k),uavg .gt. ZERO)
             wmac(i,j,k) = merge(uavg,wmac(i,j,k),test) ! varden uses ZERO instead of uavg
          enddo
       enddo
    enddo

    ! Apply boundary conditions
    do j=js,je
       do i=is,ie
          ! lo side
          if (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
             wmac(i,j,ks) = ZERO
          elseif (phys_bc(3,1) .eq. INLET) then
             wmac(i,j,ks) = u(i,j,ks-1,3)
          elseif (phys_bc(3,1) .eq. OUTLET) then
             wmac(i,j,ks) = min(wmacr(i,j,ks),ZERO)
          endif

          ! hi side
          if (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
             wmac(i,j,ke+1) = ZERO
          elseif (phys_bc(3,2) .eq. INLET) then
             wmac(i,j,ke+1) = u(i,j,ke+1,3)
          elseif (phys_bc(3,2) .eq. OUTLET) then
             wmac(i,j,ke+1) = max(wmacl(i,j,ke+1),ZERO)
          endif
       enddo
    enddo

    end if ! end CORNER-COUPLING CODE

    ! NON-CORNER-COUPLING CODE

    if (.true.) then

    !********************************
    ! Loop for edge states on x-edges.
    !********************************

    do k = ks,ke 
       do j = js,je 
          do i = is-1,ie+1 

             ! Do transverse in j direction
             if (spherical .eq. 1) then
                vlo = u(i,j  ,k,2) + HALF*(w0macy(i,j  ,k)+w0macy(i,j+1,k))
                vhi = u(i,j+1,k,2) + HALF*(w0macy(i,j+1,k)+w0macy(i,j+2,k))
             else
                vlo = u(i,j  ,k,2)
                vhi = u(i,j+1,k,2)
             end if

             spbot = u(i,j  ,k,1) + (HALF - dt2*max(ZERO,vlo)/hy)*slopey(i,j  ,k,1)
             sptop = u(i,j+1,k,1) - (HALF + dt2*min(ZERO,vhi)/hy)*slopey(i,j+1,k,1)

             sptop = merge(u(i,je+1,k,1),sptop,j.eq.je .and. phys_bc(2,2) .eq. INLET)
             spbot = merge(u(i,je+1,k,1),spbot,j.eq.je .and. phys_bc(2,2) .eq. INLET)

             if (j .eq. je .and. &
                  (phys_bc(2,2).eq.SLIP_WALL.or.phys_bc(2,2).eq.NO_SLIP_WALL)) then
                sptop = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                spbot = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
             endif

             ! upwind based on full vtrans
             if (spherical .eq. 1) then
                splus = merge(spbot,sptop,vtrans(i,j+1,k)+w0macy(i,j+1,k).gt.ZERO)
                savg  = HALF * (spbot + sptop)
                splus = merge(splus, savg, abs(vtrans(i,j+1,k)+w0macy(i,j+1,k)) .gt. rel_eps)
             else
                splus = merge(spbot,sptop,vtrans(i,j+1,k).gt.ZERO)
                savg  = HALF * (spbot + sptop)
                splus = merge(splus, savg, abs(vtrans(i,j+1,k)) .gt. rel_eps)
             end if

             if (spherical .eq. 1) then
                vlo = u(i,j-1,k,2) + HALF*(w0macy(i,j-1,k)+w0macy(i,j  ,k))
                vhi = u(i,j  ,k,2) + HALF*(w0macy(i,j  ,k)+w0macy(i,j+1,k))
             else
                vlo = u(i,j-1,k,2)
                vhi = u(i,j  ,k,2)
             end if

             smbot = u(i,j-1,k,1) + (HALF - dt2*max(ZERO,vlo)/hy)*slopey(i,j-1,k,1)
             smtop = u(i,j  ,k,1) - (HALF + dt2*min(ZERO,vhi)/hy)*slopey(i,j  ,k,1)

             smtop = merge(u(i,js-1,k,1),smtop,j.eq.js .and. phys_bc(2,1) .eq. INLET)
             smbot = merge(u(i,js-1,k,1),smbot,j.eq.js .and. phys_bc(2,1) .eq. INLET)

             if (j .eq. js .and. &
                  (phys_bc(2,1).eq.SLIP_WALL.or.phys_bc(2,1).eq.NO_SLIP_WALL)) then
                smbot = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                smtop = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
             endif

             ! upwind based on full vtrans
             if (spherical .eq. 1) then
                sminus = merge(smbot,smtop,vtrans(i,j,k)+w0macy(i,j,k).gt.ZERO)
                savg   = HALF * (smbot + smtop)
                sminus = merge(sminus, savg, abs(vtrans(i,j,k)+w0macy(i,j,k)) .gt. rel_eps)
             else
                sminus = merge(smbot,smtop,vtrans(i,j,k).gt.ZERO)
                savg   = HALF * (smbot + smtop)
                sminus = merge(sminus, savg, abs(vtrans(i,j,k)) .gt. rel_eps)
             end if

             if (spherical .eq. 1) then
                st = force(i,j,k,1) - HALF &
                     * (vtrans(i,j,k)+w0macy(i,j,k)+vtrans(i,j+1,k)+w0macy(i,j+1,k)) &
                     * (splus-sminus) / hy
             else
                st = force(i,j,k,1)-HALF*(vtrans(i,j,k)+vtrans(i,j+1,k))*(splus-sminus) / hy
             end if

             ! Do transverse in k direction
             if (spherical .eq. 1) then
                whi = u(i,j,k+1,3) + HALF*(w0macz(i,j,k+1)+w0macz(i,j,k+2))
                wlo = u(i,j,k  ,3) + HALF*(w0macz(i,j,k  )+w0macz(i,j,k+1))
             else
                wlo = u(i,j,k,3) + HALF * (w0(k)+w0(k+1))
                if (k .eq. nr(n)-1) then
                   whi = u(i,j,k+1,3) + w0(k+1)
                else
                   whi = u(i,j,k+1,3) + HALF* (w0(k+1)+w0(k+2))
                end if
             end if

             spbot = u(i,j,k  ,1) + (HALF - dt2*max(ZERO,wlo)/hz)*slopez(i,j,k  ,1)
             sptop = u(i,j,k+1,1) - (HALF + dt2*min(ZERO,whi)/hz)*slopez(i,j,k+1,1)

             sptop = merge(u(i,j,ke+1,1),sptop,k.eq.ke .and. phys_bc(3,2) .eq. INLET)
             spbot = merge(u(i,j,ke+1,1),spbot,k.eq.ke .and. phys_bc(3,2) .eq. INLET)

             if (k .eq. ke .and. &
                  (phys_bc(3,2).eq.SLIP_WALL.or.phys_bc(3,2).eq.NO_SLIP_WALL)) then
                sptop = merge(ZERO,spbot,phys_bc(3,2).eq.NO_SLIP_WALL)
                spbot = merge(ZERO,spbot,phys_bc(3,2).eq.NO_SLIP_WALL)
             endif

             ! upwind based on full wtrans
             if (spherical .eq. 1) then
                splus = merge(spbot,sptop,wtrans(i,j,k+1)+w0macz(i,j,k+1).gt.ZERO)
                savg  = HALF * (spbot + sptop)
                splus = merge(splus, savg, abs(wtrans(i,j,k+1)+w0macz(i,j,k+1)) .gt. rel_eps)
             else
                splus = merge(spbot,sptop,wtrans(i,j,k+1)+w0(k+1).gt.ZERO)
                savg  = HALF * (spbot + sptop)
                splus = merge(splus, savg, abs(wtrans(i,j,k+1)+w0(k+1)) .gt. rel_eps)
             end if

             if (spherical .eq. 1) then
                whi = u(i,j,k  ,3) + HALF*(w0macz(i,j,k  )+w0macz(i,j,k+1))
                wlo = u(i,j,k-1,3) + HALF*(w0macz(i,j,k-1)+w0macz(i,j,k  ))
             else
                if (k .eq. 0) then
                   wlo = u(i,j,k-1,3) + w0(k)
                else
                   wlo = u(i,j,k-1,3) + HALF* (w0(k-1)+w0(k))
                end if
                whi = u(i,j,k,3) + HALF* (w0(k)+w0(k+1))
             end if

             smtop = u(i,j,k  ,1) - (HALF + dt2*min(ZERO,whi)/hz)*slopez(i,j,k  ,1)
             smbot = u(i,j,k-1,1) + (HALF - dt2*max(ZERO,wlo)/hz)*slopez(i,j,k-1,1)

             smtop = merge(u(i,j,ks-1,1),smtop,k.eq.ks .and. phys_bc(3,1) .eq. INLET)
             smbot = merge(u(i,j,ks-1,1),smbot,k.eq.ks .and. phys_bc(3,1) .eq. INLET)

             if (k .eq. ks .and. &
                  (phys_bc(3,1).eq.SLIP_WALL.or.phys_bc(3,1).eq.NO_SLIP_WALL)) then
                smbot = merge(ZERO,smtop,phys_bc(3,1).eq.NO_SLIP_WALL)
                smtop = merge(ZERO,smtop,phys_bc(3,1).eq.NO_SLIP_WALL)
             endif

             ! upwind based on full wtrans
             if (spherical .eq. 1) then
                sminus = merge(smbot,smtop,wtrans(i,j,k)+w0macz(i,j,k).gt.ZERO)
                savg   = HALF * (smbot + smtop)
                sminus = merge(sminus, savg, abs(wtrans(i,j,k)+w0macz(i,j,k)) .gt. rel_eps)
             else
                sminus = merge(smbot,smtop,wtrans(i,j,k)+w0(k).gt.ZERO)
                savg   = HALF * (smbot + smtop)
                sminus = merge(sminus, savg, abs(wtrans(i,j,k)+w0(k)) .gt. rel_eps)
             end if

             if (spherical .eq. 1) then
                st = st - HALF &
                     * (wtrans(i,j,k)+w0macz(i,j,k)+wtrans(i,j,k+1)+w0macz(i,j,k+1))&
                     * (splus-sminus) / hz
             else
                st = st - HALF &
                     * (wtrans(i,j,k)+w0(k)+wtrans(i,j,k+1)+w0(k+1)) &
                     * (splus-sminus) / hz
             end if

             ! add the (Utilde . e_r) d w_0 /dr e_r term here
             if (spherical .eq. 1) then

                Ut_dot_er = HALF*(utrans(i,j,k) + utrans(i+1,j,k))*normal(i,j,k,1) + &
                            HALF*(vtrans(i,j,k) + vtrans(i,j+1,k))*normal(i,j,k,2) + &
                            HALF*(wtrans(i,j,k) + wtrans(i,j,k+1))*normal(i,j,k,3)

                st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,1)

             endif

             if (spherical .eq. 1) then
                ubardt2 = dt2/hx * ( u(i,j,k,1) + HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)) )
             else
                ubardt2 = dt2/hx * u(i,j,k,1)
             end if

             s_l(i+1)= u(i,j,k,1) + (HALF-max(ZERO,ubardt2))*slopex(i,j,k,1) + dt2*st
             s_r(i  )= u(i,j,k,1) - (HALF+min(ZERO,ubardt2))*slopex(i,j,k,1) + dt2*st

          enddo

          ! upwind based on full umac
          do i = is, ie+1 
             savg = HALF*(s_r(i) + s_l(i))
             test = ( (s_l(i).le.ZERO.and.s_r(i).ge.ZERO) .or. (abs(s_l(i)+s_r(i)).lt.rel_eps) )
             umac(i,j,k)=merge(s_l(i),s_r(i),savg.gt.ZERO)
             umac(i,j,k)=merge(savg,umac(i,j,k),test)
          enddo

          if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
             umac(is,j,k) = ZERO
          elseif (phys_bc(1,1) .eq. INLET) then
             umac(is,j,k) = u(is-1,j,k,1)
          elseif (phys_bc(1,1) .eq. OUTLET) then
             umac(is,j,k) = MIN(s_r(is),ZERO)
          endif
          if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
             umac(ie+1,j,k) = ZERO
          elseif (phys_bc(1,2) .eq. INLET) then
             umac(ie+1,j,k) = u(ie+1,j,k,1)
          elseif (phys_bc(1,2) .eq. OUTLET) then
             umac(ie+1,j,k) = MAX(s_l(ie+1),ZERO)
          endif

       enddo
    enddo

    !********************************
    ! Loop for edge states on y-edges.
    !********************************

    do k = ks, ke 
       do i = is, ie 
          do j = js-1, je+1 

             ! Do transverse in i direction
             if (spherical .eq. 1) then
                ulo = u(i  ,j,k,1) + HALF*(w0macx(i  ,j,k)+w0macx(i+1,j,k))
                uhi = u(i+1,j,k,1) + HALF*(w0macx(i+1,j,k)+w0macx(i+2,j,k))
             else
                ulo = u(i  ,j,k,1)
                uhi = u(i+1,j,k,1)
             end if

             splft = u(i  ,j,k,2) + (HALF - dt2*max(ZERO,ulo)/hx)*slopex(i  ,j,k,2)
             sprgt = u(i+1,j,k,2) - (HALF + dt2*min(ZERO,uhi)/hx)*slopex(i+1,j,k,2)

             sprgt = merge(u(ie+1,j,k,2),sprgt,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
             splft = merge(u(ie+1,j,k,2),splft,i.eq.ie .and. phys_bc(1,2) .eq. INLET)

             if (i .eq. ie .and. &
                  (phys_bc(1,2).eq.SLIP_WALL.or.phys_bc(1,2).eq.NO_SLIP_WALL)) then
                sprgt = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                splft = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
             endif

             ! upwind based on full utrans
             if (spherical .eq. 1) then
                splus = merge(splft,sprgt,utrans(i+1,j,k)+w0macx(i+1,j,k).gt.ZERO)
                savg  = HALF * (splft + sprgt)
                splus = merge(splus, savg, abs(utrans(i+1,j,k)+w0macx(i+1,j,k)) .gt. rel_eps)
             else
                splus = merge(splft,sprgt,utrans(i+1,j,k).gt.ZERO)
                savg  = HALF * (splft + sprgt)
                splus = merge(splus, savg, abs(utrans(i+1,j,k)) .gt. rel_eps)
             end if

             if (spherical .eq. 1) then
                ulo = u(i-1,j,k,1) + HALF*(w0macx(i-1,j,k)+w0macx(i  ,j,k))
                uhi = u(i  ,j,k,1) + HALF*(w0macx(i  ,j,k)+w0macx(i+1,j,k))
             else
                ulo = u(i-1,j,k,1)
                uhi = u(i  ,j,k,1)
             end if

             smlft = u(i-1,j,k,2) + (HALF - dt2*max(ZERO,ulo)/hx)*slopex(i-1,j,k,2)
             smrgt = u(i  ,j,k,2) - (HALF + dt2*min(ZERO,uhi)/hx)*slopex(i  ,j,k,2)

             smrgt = merge(u(is-1,j,k,2),smrgt,i.eq.is .and. phys_bc(1,1) .eq. INLET)
             smlft = merge(u(is-1,j,k,2),smlft,i.eq.is .and. phys_bc(1,1) .eq. INLET)

             if (i .eq. is .and. &
                  (phys_bc(1,1).eq.SLIP_WALL.or.phys_bc(1,1).eq.NO_SLIP_WALL)) then
                smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
             endif

             ! upwind based on full utrans
             if (spherical .eq. 1) then
                sminus = merge(smlft,smrgt,utrans(i,j,k)+w0macx(i,j,k).gt.ZERO)
                savg   = HALF * (smlft + smrgt)
                sminus = merge(sminus, savg, abs(utrans(i,j,k)+w0macx(i,j,k)) .gt. rel_eps)
             else
                sminus = merge(smlft,smrgt,utrans(i,j,k).gt.ZERO)
                savg   = HALF * (smlft + smrgt)
                sminus = merge(sminus, savg, abs(utrans(i,j,k)) .gt. rel_eps)
             end if

             if (spherical .eq. 1) then
                st = force(i,j,k,2) - HALF &
                     * (utrans(i,j,k)+w0macx(i,j,k)+utrans(i+1,j,k)+w0macx(i+1,j,k)) &
                     * (splus-sminus) / hx
             else
                st = force(i,j,k,2) - HALF*(utrans(i,j,k)+utrans(i+1,j,k))*(splus-sminus) / hx
             end if

             ! Do transverse in k direction
             if (spherical .eq. 1) then
                whi = u(i,j,k+1,3) + HALF*(w0macz(i,j,k+1)+w0macz(i,j,k+2))
                wlo = u(i,j,k  ,3) + HALF*(w0macz(i,j,k  )+w0macz(i,j,k+1))
             else

                wlo = u(i,j,k,3) + HALF* (w0(k)+w0(k+1))
                if (k .eq. nr(n)-1) then
                   whi = u(i,j,k+1,3) + w0(k+1)
                else
                   whi = u(i,j,k+1,3) + HALF* (w0(k+1)+w0(k+2))
                end if

             end if

             splft = u(i,j,k  ,2) + (HALF - dt2*max(ZERO,wlo)/hz)*slopez(i,j,k  ,2)
             sprgt = u(i,j,k+1,2) - (HALF + dt2*min(ZERO,whi)/hz)*slopez(i,j,k+1,2)

             sprgt = merge(u(i,j,ke+1,2),sprgt,k.eq.ke .and. phys_bc(3,2) .eq. INLET)
             splft = merge(u(i,j,ke+1,2),splft,k.eq.ke .and. phys_bc(3,2) .eq. INLET)

             if (k .eq. ke .and. &
                  (phys_bc(3,2).eq.SLIP_WALL.or.phys_bc(3,2).eq.NO_SLIP_WALL)) then
                sprgt = merge(ZERO,splft,phys_bc(3,2).eq.NO_SLIP_WALL)
                splft = merge(ZERO,splft,phys_bc(3,2).eq.NO_SLIP_WALL)
             endif

             ! upwind based on full wtrans
             if (spherical .eq. 1) then
                splus = merge(splft,sprgt,wtrans(i,j,k+1)+w0macz(i,j,k+1).gt.ZERO)
                savg  = HALF * (splft + sprgt)
                splus = merge(splus, savg, abs(wtrans(i,j,k+1)+w0macz(i,j,k+1)) .gt. rel_eps)
             else
                splus = merge(splft,sprgt,wtrans(i,j,k+1)+w0(k+1).gt.ZERO)
                savg  = HALF * (splft + sprgt)
                splus = merge(splus, savg, abs(wtrans(i,j,k+1)+w0(k+1)) .gt. rel_eps)
             end if

             if (spherical .eq. 1) then
                whi = u(i,j,k  ,3) + HALF*(w0macz(i,j,k  )+w0macz(i,j,k+1))
                wlo = u(i,j,k-1,3) + HALF*(w0macz(i,j,k-1)+w0macz(i,j,k  ))
             else

                if (k .eq. 0) then
                   wlo = u(i,j,k-1,3) + w0(k)
                else
                   wlo = u(i,j,k-1,3) + HALF* (w0(k-1)+w0(k))
                end if
                whi = u(i,j,k,3) + HALF* (w0(k)+w0(k+1))

             end if

             smrgt = u(i,j,k  ,2) - (HALF + dt2*min(ZERO,whi)/hz)*slopez(i,j,k  ,2)
             smlft = u(i,j,k-1,2) + (HALF - dt2*max(ZERO,wlo)/hz)*slopez(i,j,k-1,2)

             smrgt = merge(u(i,j,ks-1,2),smrgt,k.eq.ks .and. phys_bc(3,1) .eq. INLET)
             smlft = merge(u(i,j,ks-1,2),smlft,k.eq.ks .and. phys_bc(3,1) .eq. INLET)

             if (k .eq. ks .and. &
                  (phys_bc(3,1).eq.SLIP_WALL.or.phys_bc(3,1).eq.NO_SLIP_WALL)) then
                smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
             endif

             ! upwind based on full wtrans
             if (spherical .eq. 1) then
                sminus = merge(smlft,smrgt,wtrans(i,j,k)+w0macz(i,j,k).gt.ZERO)
                savg   = HALF * (smlft + smrgt)
                sminus = merge(sminus, savg, abs(wtrans(i,j,k)+w0macz(i,j,k)) .gt. rel_eps)
             else
                sminus = merge(smlft,smrgt,wtrans(i,j,k)+w0(k).gt.ZERO)
                savg   = HALF * (smlft + smrgt)
                sminus = merge(sminus, savg, abs(wtrans(i,j,k)+w0(k)) .gt. rel_eps)
             end if

             if (spherical .eq. 1) then
                st = st - HALF &
                     * (wtrans(i,j,k)+w0macz(i,j,k)+wtrans(i,j,k+1)+w0macz(i,j,k+1)) &
                     * (splus-sminus) / hz
             else
                st = st - HALF &
                     * (wtrans(i,j,k)+w0(k)+wtrans(i,j,k+1)+w0(k+1)) &
                     * (splus-sminus) / hz
             end if

             ! add the (Utilde . e_r) d w_0 /dr e_r term here
             if (spherical .eq. 1) then

                Ut_dot_er = HALF*(utrans(i,j,k) + utrans(i+1,j,k))*normal(i,j,k,1) + &
                            HALF*(vtrans(i,j,k) + vtrans(i,j+1,k))*normal(i,j,k,2) + &
                            HALF*(wtrans(i,j,k) + wtrans(i,j,k+1))*normal(i,j,k,3)

                st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,2)

             endif

             if (spherical .eq. 1) then
                vbardt2 = dt2/hy * ( u(i,j,k,2) + HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)) )
             else
                vbardt2 = dt2/hy * u(i,j,k,2)
             end if

             s_b(j+1)= u(i,j,k,2) + (HALF-max(ZERO,vbardt2))*slopey(i,j,k,2) + dt2*st
             s_t(j  )= u(i,j,k,2) - (HALF+min(ZERO,vbardt2))*slopey(i,j,k,2) + dt2*st

          enddo

          ! upwind based on full vmac
          do j = js, je+1 
             savg = HALF*(s_b(j) + s_t(j))
             test = ( (s_b(j).le.ZERO.and.s_t(j).ge.ZERO) .or. (abs(s_b(j)+s_t(j)).lt.rel_eps) )
             vmac(i,j,k)=merge(s_b(j),s_t(j),savg.gt.ZERO)
             vmac(i,j,k)=merge(savg,vmac(i,j,k),test)
          enddo

          if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
             vmac(i,js,k) = ZERO
          elseif (phys_bc(2,1) .eq. INLET) then
             vmac(i,js,k) = u(i,js-1,k,2)
          elseif (phys_bc(2,1) .eq. OUTLET) then
             vmac(i,js,k) = MIN(s_t(js),ZERO)
          endif

          if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
             vmac(i,je+1,k) = ZERO
          elseif (phys_bc(2,2) .eq. INLET) then
             vmac(i,je+1,k) = u(i,je+1,k,2)
          elseif (phys_bc(2,2) .eq. OUTLET) then
             vmac(i,je+1,k) = MAX(s_b(je+1),ZERO)
          endif

       enddo
    enddo

    !********************************
    ! Loop for edge states on z-edges.
    !********************************

    do j = js, je 
       do i = is, ie 
          do k = ks-1,ke+1

             ! Do transverse in i direction
             if (spherical .eq. 1) then
                ulo = u(i  ,j,k,1) + HALF*(w0macx(i  ,j,k)+w0macx(i+1,j,k))
                uhi = u(i+1,j,k,1) + HALF*(w0macx(i+1,j,k)+w0macx(i+2,j,k))
             else
                ulo = u(i  ,j,k,1)
                uhi = u(i+1,j,k,1)
             end if

             splft = u(i  ,j,k,3) + (HALF - dt2*max(ZERO,ulo)/hx)*slopex(i  ,j,k,3)
             sprgt = u(i+1,j,k,3) - (HALF + dt2*min(ZERO,uhi)/hx)*slopex(i+1,j,k,3)

             sprgt = merge(u(ie+1,j,k,3),sprgt,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
             splft = merge(u(ie+1,j,k,3),splft,i.eq.ie .and. phys_bc(1,2) .eq. INLET)

             if (i .eq. ie .and. &
                  (phys_bc(1,2).eq.SLIP_WALL.or.phys_bc(1,2).eq.NO_SLIP_WALL)) then
                sprgt = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                splft = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
             endif

             ! upwind based on full utrans
             if (spherical .eq. 1) then
                splus = merge(splft,sprgt,utrans(i+1,j,k)+w0macx(i+1,j,k).gt.ZERO)
                savg  = HALF * (splft + sprgt)
                splus = merge(splus, savg, abs(utrans(i+1,j,k)+w0macx(i+1,j,k)) .gt. rel_eps)
             else
                splus = merge(splft,sprgt,utrans(i+1,j,k).gt.ZERO)
                savg  = HALF * (splft + sprgt)
                splus = merge(splus, savg, abs(utrans(i+1,j,k)) .gt. rel_eps)
             end if

             if (spherical .eq. 1) then
                ulo = u(i-1,j,k,1) + HALF*(w0macx(i-1,j,k)+w0macx(i  ,j,k))
                uhi = u(i  ,j,k,1) + HALF*(w0macx(i  ,j,k)+w0macx(i+1,j,k))
             else
                ulo = u(i-1,j,k,1)
                uhi = u(i  ,j,k,1)
             end if

             smlft = u(i-1,j,k,3) + (HALF - dt2*max(ZERO,ulo)/hx)*slopex(i-1,j,k,3)
             smrgt = u(i  ,j,k,3) - (HALF + dt2*min(ZERO,uhi)/hx)*slopex(i  ,j,k,3)

             smrgt = merge(u(is-1,j,k,3),smrgt,i.eq.is .and. phys_bc(1,1) .eq. INLET)
             smlft = merge(u(is-1,j,k,3),smlft,i.eq.is .and. phys_bc(1,1) .eq. INLET)

             if (i .eq. is .and. &
                  (phys_bc(1,1).eq.SLIP_WALL.or.phys_bc(1,1).eq.NO_SLIP_WALL)) then
                smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
             endif

             ! upwind based on full utrans
             if (spherical .eq. 1) then
                sminus = merge(smlft,smrgt,utrans(i,j,k)+w0macx(i,j,k).gt.ZERO)
                savg   = HALF * (smlft + smrgt)
                sminus = merge(sminus, savg, abs(utrans(i,j,k)+w0macx(i,j,k)) .gt. rel_eps)
             else
                sminus = merge(smlft,smrgt,utrans(i,j,k).gt.ZERO)
                savg   = HALF * (smlft + smrgt)
                sminus = merge(sminus, savg, abs(utrans(i,j,k)) .gt. rel_eps)
             end if

             if (spherical .eq. 1) then
                st = force(i,j,k,3) - HALF &
                     * (utrans(i,j,k)+w0macx(i,j,k)+utrans(i+1,j,k)+w0macx(i+1,j,k)) &
                     * (splus-sminus) / hx
             else
                st = force(i,j,k,3) - HALF*(utrans(i,j,k)+utrans(i+1,j,k))*(splus-sminus) / hx
             end if

             ! Do transverse in j direction
             if (spherical .eq. 1) then
                vlo = u(i,j  ,k,2) + HALF*(w0macy(i,j  ,k)+w0macy(i,j+1,k))
                vhi = u(i,j+1,k,2) + HALF*(w0macy(i,j+1,k)+w0macy(i,j+2,k))
             else
                vlo = u(i,j  ,k,2)
                vhi = u(i,j+1,k,2)
             end if

             spbot = u(i,j  ,k,3) + (HALF - dt2*max(ZERO,vlo)/hy)*slopey(i,j  ,k,3)
             sptop = u(i,j+1,k,3) - (HALF + dt2*min(ZERO,vhi)/hy)*slopey(i,j+1,k,3)

             sptop = merge(u(i,je+1,k,3),sptop,j.eq.je .and. phys_bc(2,2) .eq. INLET)
             spbot = merge(u(i,je+1,k,3),spbot,j.eq.je .and. phys_bc(2,2) .eq. INLET)

             if (j .eq. je .and. &
                  (phys_bc(2,2).eq.SLIP_WALL.or.phys_bc(2,2).eq.NO_SLIP_WALL)) then
                sptop = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                spbot = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
             endif

             ! upwind based on full vtrans
             if (spherical .eq. 1) then
                splus = merge(spbot,sptop,vtrans(i,j+1,k)+w0macy(i,j+1,k).gt.ZERO)
                savg  = HALF * (spbot + sptop)
                splus = merge(splus, savg, abs(vtrans(i,j+1,k)+w0macy(i,j+1,k)) .gt. rel_eps)
             else
                splus = merge(spbot,sptop,vtrans(i,j+1,k).gt.ZERO)
                savg  = HALF * (spbot + sptop)
                splus = merge(splus, savg, abs(vtrans(i,j+1,k)) .gt. rel_eps)
             end if

             if (spherical .eq. 1) then
                vlo = u(i,j-1,k,2) + HALF*(w0macy(i,j-1,k)+w0macy(i,j  ,k))
                vhi = u(i,j  ,k,2) + HALF*(w0macy(i,j  ,k)+w0macy(i,j+1,k))
             else
                vlo = u(i,j-1,k,2)
                vhi = u(i,j  ,k,2)
             end if

             smbot = u(i,j-1,k,3) + (HALF - dt2*max(ZERO,vlo)/hy)*slopey(i,j-1,k,3)
             smtop = u(i,j  ,k,3) - (HALF + dt2*min(ZERO,vhi)/hy)*slopey(i,j  ,k,3)

             smtop = merge(u(i,js-1,k,3),smtop,j.eq.js .and. phys_bc(2,1) .eq. INLET)
             smbot = merge(u(i,js-1,k,3),smbot,j.eq.js .and. phys_bc(2,1) .eq. INLET)

             if (j .eq. js .and. &
                  (phys_bc(2,1).eq.SLIP_WALL.or.phys_bc(2,1).eq.NO_SLIP_WALL)) then
                smbot = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                smtop = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
             endif

             ! upwind based on full vtrans
             if (spherical .eq. 1) then
                sminus = merge(smbot,smtop,vtrans(i,j,k)+w0macy(i,j,k).gt.ZERO)
                savg   = HALF * (smbot + smtop)
                sminus = merge(sminus, savg, abs(vtrans(i,j,k)+w0macy(i,j,k)) .gt. rel_eps)
             else
                sminus = merge(smbot,smtop,vtrans(i,j,k).gt.ZERO)
                savg   = HALF * (smbot + smtop)
                sminus = merge(sminus, savg, abs(vtrans(i,j,k)) .gt. rel_eps)
             end if

             if (spherical .eq. 1) then
                st = st - HALF &
                     * (vtrans(i,j,k)+w0macy(i,j,k)+vtrans(i,j+1,k)+w0macy(i,j+1,k)) &
                     * (splus - sminus) / hy
             else
                st = st - HALF * (vtrans(i,j,k)+vtrans(i,j+1,k))*(splus - sminus) / hy
             end if

             ! add the (Utilde . e_r) d w_0 /dr e_r term here
             if (spherical .eq. 0) then

                if (k .ge. 0 .and. k .le. nr(n)-1) then
                   st = st - HALF * (wtrans(i,j,k)+wtrans(i,j,k+1))*(w0(k+1)-w0(k)) / hz
                else
                   ! dw0/dr=0 and therefore st is unchanged
                end if

             else if (spherical .eq. 1) then

                Ut_dot_er = HALF*(utrans(i,j,k) + utrans(i+1,j,k))*normal(i,j,k,1) + &
                            HALF*(vtrans(i,j,k) + vtrans(i,j+1,k))*normal(i,j,k,2) + &
                            HALF*(wtrans(i,j,k) + wtrans(i,j,k+1))*normal(i,j,k,3)

                st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,3)

             endif

             if (spherical .eq. 1) then
                wbardt2 = dt2/hz * ( u(i,j,k,3) + HALF*(w0macz(i,j,k)+w0macz(i,j,k+1)) )
             else
                if (k .lt. 0) then
                   wbardt2 = dt2/hz * ( u(i,j,k,3) + w0(k+1) )
                else if(k .gt. nr(n)-1) then
                   wbardt2 = dt2/hz * ( u(i,j,k,3) + w0(k) )
                else
                   wbardt2 = dt2/hz * ( u(i,j,k,3) + HALF*(w0(k)+w0(k+1)) )
                end if
             end if

             s_d(k+1)= u(i,j,k,3) + (HALF-max(ZERO,wbardt2))*slopez(i,j,k,3) + dt2*st
             s_u(k  )= u(i,j,k,3) - (HALF+min(ZERO,wbardt2))*slopez(i,j,k,3) + dt2*st

          enddo

          ! upwind based on full wmac
          if (spherical .eq. 1) then
             do k = ks, ke+1 
                savg = HALF*(s_d(k) + s_u(k))
                test = ( (s_d(k)+w0macz(i,j,k).le.ZERO.and.s_u(k)+w0macz(i,j,k).ge.ZERO) &
                     .or. (abs(s_d(k)+s_u(k)+TWO*w0macz(i,j,k)).lt.rel_eps) )
                wmac(i,j,k)=merge(s_d(k),s_u(k),savg+w0macz(i,j,k).gt.ZERO)
                wmac(i,j,k)=merge(savg,wmac(i,j,k),test)
             enddo
          else
             do k = ks, ke+1 
                savg = HALF*(s_d(k) + s_u(k))
                test = ( (s_d(k)+w0(k).le.ZERO.and.s_u(k)+w0(k).ge.ZERO) &
                     .or. (abs(s_d(k)+s_u(k)+TWO*w0(k)).lt.rel_eps) )
                wmac(i,j,k)=merge(s_d(k),s_u(k),savg+w0(k).gt.ZERO)
                wmac(i,j,k)=merge(savg,wmac(i,j,k),test)
             enddo
          end if

          if (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
             wmac(i,j,ks) = ZERO
          elseif (phys_bc(3,1) .eq. INLET) then
             wmac(i,j,ks) = u(i,j,ks-1,3)
          elseif (phys_bc(3,1) .eq. OUTLET) then
             wmac(i,j,ks) = MIN(s_u(ks),ZERO)
          endif

          if (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
             wmac(i,j,ke+1 ) = ZERO
          elseif (phys_bc(3,2) .eq. INLET) then
             wmac(i,j,ke+1) = u(i,j,ke+1,3)
          elseif (phys_bc(3,2) .eq. OUTLET) then
             wmac(i,j,ke+1) = MAX(s_d(ke+1),ZERO)
          endif

       enddo
    enddo

    end if ! end NON-CORNER-COUPLING CODE

    deallocate(ulx,urx,uimhx,uly,ury,uimhy,ulz,urz,uimhz)
    deallocate(ulyz,uryz,uimhyz,ulzy,urzy,uimhzy)
    deallocate(vlxz,vrxz,vimhxz,vlzx,vrzx,vimhzx)
    deallocate(wlxy,wrxy,wimhxy,wlyx,wryx,wimhyx)
    deallocate(umacl,umacr,vmacl,vmacr,wmacl,wmacr)

  end subroutine velpred_3d

end module velpred_module
