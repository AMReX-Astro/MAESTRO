! velpred is called by advance_premac -- it is used to predict the
! normal velocities to the interfaces.  We don't care about the
! transverse velocities here.  The prediction is done piecewise linear
! or with PPM depending on ppm_type.

module velpred_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: velpred

contains

  subroutine velpred(u,ufull,umac,utrans,force,w0,w0mac,dx,dt,the_bc_level,mla)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: spherical
    use fill_3d_module
    use multifab_physbc_module
    use ml_cc_restriction_module, only : ml_edge_restriction_c

    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(in   ) :: ufull(:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: utrans(:,:),force(:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(in   ) :: mla

    integer                  :: i,n,n_1d
    integer                  :: ng_u,ng_uf,ng_um,ng_ut,ng_f,ng_w0
    integer                  :: lo(mla%dim), hi(mla%dim),dm,nlevs
    real(kind=dp_t), pointer :: uop(:,:,:,:)
    real(kind=dp_t), pointer :: ufp(:,:,:,:)
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

    type(bl_prof_timer), save :: bpt

    call build(bpt, "velpred")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_u  = nghost(u(1))
    ng_uf = nghost(ufull(1))
    ng_um = nghost(umac(1,1))
    ng_ut = nghost(utrans(1,1))
    ng_f  = nghost(force(1))
    ng_w0 = nghost(w0mac(1,1))

    do n=1,nlevs
       do i = 1, nfabs(u(n))
          uop  => dataptr(u(n),i)
          ufp  => dataptr(ufull(n),i)
          ump  => dataptr(umac(n,1),i)
          utp  => dataptr(utrans(n,1),i)
          fp   => dataptr(force(n),i)
          lo   =  lwb(get_box(u(n),i))
          hi   =  upb(get_box(u(n),i))
          select case (dm)
          case (1)
             call velpred_1d(uop(:,1,1,:), ng_u, &
                             ufp(:,1,1,:), ng_uf, &
                             ump(:,1,1,1), ng_um, &
                             fp(:,1,1,1), ng_f, w0(n,:), lo, hi, dx(n,:), dt, &
                             the_bc_level(n)%phys_bc_level_array(i,:,:), &
                             the_bc_level(n)%adv_bc_level_array(i,:,:,:))
          case (2)
             vtp  => dataptr(utrans(n,2),i)
             vmp  => dataptr(  umac(n,2),i)
             call velpred_2d(uop(:,:,1,:), ng_u, &
                             ufp(:,:,1,:), ng_uf, &
                             utp(:,:,1,1), vtp(:,:,1,1), ng_ut, &
                             ump(:,:,1,1), vmp(:,:,1,1), ng_um, &
                             fp(:,:,1,:), ng_f, w0(n,:), lo, hi, dx(n,:), dt, &
                             the_bc_level(n)%phys_bc_level_array(i,:,:), &
                             the_bc_level(n)%adv_bc_level_array(i,:,:,:))

          case (3)
             vmp  => dataptr(  umac(n,2),i)
             wmp  => dataptr(  umac(n,3),i)
             vtp  => dataptr(utrans(n,2),i)
             wtp  => dataptr(utrans(n,3),i)
             w0xp  => dataptr(w0mac(n,1),i)
             w0yp  => dataptr(w0mac(n,2),i)
             w0zp  => dataptr(w0mac(n,3),i)
             if (spherical .eq. 1) then
                n_1d = 1
             else
                n_1d = n
             end if
             call velpred_3d(uop(:,:,:,:), ng_u, &
                             ufp(:,:,:,:), ng_uf, &
                             ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                             utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), ng_ut, &
                             fp(:,:,:,:), ng_f, &
                             w0(n_1d,:),w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1), &
                             ng_w0, lo, hi, dx(n,:), dt, &
                             the_bc_level(n)%phys_bc_level_array(i,:,:), &
                             the_bc_level(n)%adv_bc_level_array(i,:,:,:))
          end select
       end do
    end do

    do n = nlevs,2,-1
       do i = 1, dm
          call ml_edge_restriction_c(umac(n-1,i),1,umac(n,i),1,mla%mba%rr(n-1,:),i,1)
       enddo
    enddo

    call destroy(bpt)

  end subroutine velpred

  subroutine velpred_1d(u,ng_u,ufull,ng_uf,umac,ng_um,force,ng_f, &
                        w0,lo,hi,dx,dt,phys_bc,adv_bc)

    use bc_module
    use slope_module
    use bl_constants_module
    use variables, only: rel_eps
    use probin_module, only: ppm_type, ppm_trace_forces
    use ppm_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_uf,ng_um,ng_f
    real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u :,:)
    real(kind=dp_t), intent(in   ) ::  ufull(lo(1)-ng_uf:,:)
    real(kind=dp_t), intent(inout) ::   umac(lo(1)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)-ng_f :)
    real(kind=dp_t), intent(in   ) ::     w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    integer        , intent(in   ) :: phys_bc(:,:)
    integer        , intent(in   ) :: adv_bc(:,:,:)

    ! Local variables
    real(kind=dp_t) :: slopex(lo(1)-1:hi(1)+1,1)

    real(kind=dp_t), allocatable :: Ipu(:), Ipf(:)
    real(kind=dp_t), allocatable :: Imu(:), Imf(:)

    ! these correspond to umac_L, etc.
    real(kind=dp_t), allocatable :: umacl(:),umacr(:)

    real(kind=dp_t) :: hx, dt2, dt4, uavg

    integer :: i,is,ie

    logical :: test

    allocate(Ipu(lo(1)-1:hi(1)+1))
    allocate(Imu(lo(1)-1:hi(1)+1))

    allocate(Ipf(lo(1)-1:hi(1)+1))
    allocate(Imf(lo(1)-1:hi(1)+1))

    allocate(umacl(lo(1):hi(1)+1))
    allocate(umacr(lo(1):hi(1)+1))

    is = lo(1)
    ie = hi(1)

    dt2 = HALF*dt
    dt4 = dt/4.0d0

    hx = dx(1)

    if (ppm_type .eq. 0) then
       call slopex_1d(u,slopex,lo,hi,ng_u,1,adv_bc)

    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_1d(u(:,1),ng_u,ufull(:,1),ng_uf,Ipu,Imu,lo,hi,adv_bc(:,:,1),dx,dt,.false.)
       if (ppm_trace_forces .eq. 1) then
          call ppm_1d(force(:),ng_f,ufull(:,1),ng_uf,Ipf,Imf,lo,hi,adv_bc(:,:,1),dx,dt,.false.)
       endif
    end if

    !******************************************************************
    ! Create umac 
    !******************************************************************

    if (ppm_type .eq. 0) then
       do i=is,ie+1
          ! extrapolate velocity to left face
          umacl(i) = u(i-1,1) + (HALF-(dt2/hx)*max(ZERO,ufull(i-1,1)))*slopex(i-1,1) &
               + dt2*force(i-1)
          ! extrapolate velocity to right face
          umacr(i) = u(i  ,1) - (HALF+(dt2/hx)*min(ZERO,ufull(i  ,1)))*slopex(i  ,1) &
               + dt2*force(i  )
       end do
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       if (ppm_trace_forces .eq. 0) then
          do i=is,ie+1
             ! extrapolate velocity to left face
             umacl(i) = Ipu(i-1) + dt2*force(i-1)
             ! extrapolate velocity to right face
             umacr(i) = Imu(i  ) + dt2*force(i  )
          end do
       else 
          do i=is,ie+1
             ! extrapolate velocity to left face
             umacl(i) = Ipu(i-1) + dt2*Ipf(i-1)
             ! extrapolate velocity to right face
             umacr(i) = Imu(i  ) + dt2*Imf(i  )
          end do
       endif
    end if

    do i=is,ie+1
       ! solve Riemann problem using full velocity
       uavg = HALF*(umacl(i)+umacr(i))
       test = ((umacl(i)+w0(i) .le. ZERO .and. umacr(i)+w0(i) .ge. ZERO) .or. &
           (abs(umacl(i)+umacr(i)+TWO*w0(i)) .lt. rel_eps))
       umac(i) = merge(umacl(i),umacr(i),uavg+w0(i) .gt. ZERO)
       umac(i) = merge(ZERO,umac(i),test)
    enddo

    ! impose lo side bc's
    select case(phys_bc(1,1))
    case (INLET)
       umac(is) = u(is-1,1)
    case (SLIP_WALL, NO_SLIP_WALL, SYMMETRY)
       umac(is) = ZERO
    case (OUTLET)
       umac(is) = min(umacr(is),ZERO)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_1d: invalid boundary type phys_bc(1,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(1,2))
    case (INLET)
       umac(ie+1) = u(ie+1,1)
    case (SLIP_WALL, NO_SLIP_WALL, SYMMETRY)
       umac(ie+1) = ZERO
    case (OUTLET)
       umac(ie+1) = max(umacl(ie+1),ZERO)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_1d: invalid boundary type phys_bc(1,2)")
    end select

    deallocate(umacl,umacr)
    deallocate(Ipu,Imu)

  end subroutine velpred_1d

  subroutine velpred_2d(u,ng_u,ufull,ng_uf,utrans,vtrans,ng_ut,umac,vmac,ng_um,force,ng_f, &
                        w0,lo,hi,dx,dt,phys_bc,adv_bc)

    use bc_module
    use slope_module
    use bl_constants_module
    use variables, only: rel_eps
    use probin_module, only: ppm_type, ppm_trace_forces
    use ppm_module

    integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_uf,ng_um,ng_ut,ng_f
    real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u :,lo(2)-ng_u :,:)
    real(kind=dp_t), intent(in   ) ::  ufull(lo(1)-ng_uf:,lo(2)-ng_uf:,:)
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

    real(kind=dp_t), allocatable :: Ipu(:,:,:), Ipfx(:,:,:)
    real(kind=dp_t), allocatable :: Imu(:,:,:), Imfx(:,:,:)
    real(kind=dp_t), allocatable :: Ipv(:,:,:), Ipfy(:,:,:)
    real(kind=dp_t), allocatable :: Imv(:,:,:), Imfy(:,:,:)

    ! these correspond to u_L^x, etc.
    real(kind=dp_t), allocatable :: ulx(:,:,:),urx(:,:,:),uimhx(:,:,:)
    real(kind=dp_t), allocatable :: uly(:,:,:),ury(:,:,:),uimhy(:,:,:)

    ! these correspond to umac_L, etc.
    real(kind=dp_t), allocatable :: umacl(:,:),umacr(:,:)
    real(kind=dp_t), allocatable :: vmacl(:,:),vmacr(:,:)

    real(kind=dp_t) :: hx, hy, dt2, dt4, uavg, maxu, minu
    real(kind=dp_t) :: fl, fr

    integer :: i,j,is,js,ie,je

    logical :: test

    allocate(Ipu(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(Imu(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(Ipv(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(Imv(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))

    allocate(Ipfx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(Imfx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(Ipfy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(Imfy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))

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

    dt2 = HALF*dt
    dt4 = dt/4.0d0

    hx = dx(1)
    hy = dx(2)

    if (ppm_type .eq. 0) then
       call slopex_2d(u,slopex,lo,hi,ng_u,2,adv_bc)
       call slopey_2d(u,slopey,lo,hi,ng_u,2,adv_bc)
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_2d(u(:,:,1),ng_u, &
                   ufull(:,:,1),ufull(:,:,2),ng_uf, &
                   Ipu,Imu,lo,hi,adv_bc(:,:,1),dx,dt,.false.)
       call ppm_2d(u(:,:,2),ng_u, &
                   ufull(:,:,1),ufull(:,:,2),ng_uf, &
                   Ipv,Imv,lo,hi,adv_bc(:,:,2),dx,dt,.false.)

       ! trace forces, if necessary.  Note by default the ppm routines
       ! will trace each component to each interface in all coordinate
       ! directions, but we really only need the force traced along
       ! its respective dimension.  This should be simplified later.
       if (ppm_trace_forces == 1) then
          call ppm_2d(force(:,:,1),ng_f, &
                      ufull(:,:,1),ufull(:,:,2),ng_uf, &
                      Ipfx,Imfx,lo,hi,adv_bc(:,:,1),dx,dt,.false.)
          call ppm_2d(force(:,:,2),ng_f, &
                      ufull(:,:,1),ufull(:,:,2),ng_uf, &
                      Ipfy,Imfy,lo,hi,adv_bc(:,:,2),dx,dt,.false.)
       endif
    end if
       
    !******************************************************************
    ! Create u_{\i-\half\e_x}^x, etc.
    !******************************************************************

    if (ppm_type .eq. 0) then
       do j=js-1,je+1
          do i=is,ie+1
             maxu = max(ZERO,ufull(i-1,j,1))
             minu = min(ZERO,ufull(i  ,j,1))
             ! extrapolate both components of velocity to left face
             ulx(i,j,1) = u(i-1,j,1) + (HALF - (dt2/hx)*maxu)*slopex(i-1,j,1)
             ulx(i,j,2) = u(i-1,j,2) + (HALF - (dt2/hx)*maxu)*slopex(i-1,j,2)
             ! extrapolate both components of velocity to right face
             urx(i,j,1) = u(i  ,j,1) - (HALF + (dt2/hx)*minu)*slopex(i  ,j,1)
             urx(i,j,2) = u(i  ,j,2) - (HALF + (dt2/hx)*minu)*slopex(i  ,j,2)
          end do
       end do
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do j=js-1,je+1
          do i=is,ie+1
             ! extrapolate both components of velocity to left face
             ulx(i,j,1) = Ipu(i-1,j,1)
             ulx(i,j,2) = Ipv(i-1,j,1)
             ! extrapolate both components of velocity to right face
             urx(i,j,1) = Imu(i,j,1)
             urx(i,j,2) = Imv(i,j,1)
          end do
       end do
    end if
    
    ! impose lo side bc's
    select case(phys_bc(1,1))
    case (INLET)
       ulx(is,js-1:je+1,1:2) = u(is-1,js-1:je+1,1:2)
       urx(is,js-1:je+1,1:2) = u(is-1,js-1:je+1,1:2)
    case (SLIP_WALL, SYMMETRY)
       ulx(is,js-1:je+1,1) = ZERO
       urx(is,js-1:je+1,1) = ZERO
       ulx(is,js-1:je+1,2) = urx(is,js-1:je+1,2)
    case (NO_SLIP_WALL)
       ulx(is,js-1:je+1,1:2) = ZERO
       urx(is,js-1:je+1,1:2) = ZERO
    case (OUTLET)
       ulx(is,js-1:je+1,1) = min(urx(is,js-1:je+1,1),ZERO)
       urx(is,js-1:je+1,1) = ulx(is,js-1:je+1,1)
       ulx(is,js-1:je+1,2) = urx(is,js-1:je+1,2)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_2d: invalid boundary type phys_bc(1,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(1,2))
    case (INLET)
       ulx(ie+1,js-1:je+1,1:2) = u(ie+1,js-1:je+1,1:2)
       urx(ie+1,js-1:je+1,1:2) = u(ie+1,js-1:je+1,1:2)
    case (SLIP_WALL, SYMMETRY)
       ulx(ie+1,js-1:je+1,1) = ZERO
       urx(ie+1,js-1:je+1,1) = ZERO
       urx(ie+1,js-1:je+1,2) = ulx(ie+1,js-1:je+1,2)
    case (NO_SLIP_WALL)
       ulx(ie+1,js-1:je+1,1:2) = ZERO
       urx(ie+1,js-1:je+1,1:2) = ZERO
    case (OUTLET)
       ulx(ie+1,js-1:je+1,1) = max(ulx(ie+1,js-1:je+1,1),ZERO)
       urx(ie+1,js-1:je+1,1) = ulx(ie+1,js-1:je+1,2)
       urx(ie+1,js-1:je+1,2) = ulx(ie+1,js-1:je+1,2)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_2d: invalid boundary type phys_bc(1,2)")
    end select

    do j=js-1,je+1
       do i=is,ie+1
          ! No need to compute uimhx(:,:,1) since it's equal to utrans-w0
          ! upwind using full velocity to get transverse component of uimhx
          ! Note: utrans already contains w0
          uimhx(i,j,2) = merge(ulx(i,j,2),urx(i,j,2),utrans(i,j).gt.ZERO)
          uavg = HALF*(ulx(i,j,2)+urx(i,j,2))
          uimhx(i,j,2) = merge(uavg,uimhx(i,j,2),abs(utrans(i,j)).lt.rel_eps)
       enddo
    enddo

    if (ppm_type .eq. 0) then
       do j=js,je+1
          do i=is-1,ie+1
             maxu = max(ZERO,ufull(i,j-1,2))
             minu = min(ZERO,ufull(i,j  ,2))
             ! extrapolate both components of velocity to left face
             uly(i,j,1) = u(i,j-1,1) + (HALF-(dt2/hy)*maxu)*slopey(i,j-1,1)
             uly(i,j,2) = u(i,j-1,2) + (HALF-(dt2/hy)*maxu)*slopey(i,j-1,2)
             ! extrapolate both components of velocity to right face
             ury(i,j,1) = u(i,j  ,1) - (HALF+(dt2/hy)*minu)*slopey(i,j  ,1)
             ury(i,j,2) = u(i,j  ,2) - (HALF+(dt2/hy)*minu)*slopey(i,j  ,2)
          end do
       end do
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do j=js,je+1
          do i=is-1,ie+1
             ! extrapolate both components of velocity to left face
             uly(i,j,1) = Ipu(i,j-1,2)
             uly(i,j,2) = Ipv(i,j-1,2)
             ! extrapolate both components of velocity to right face
             ury(i,j,1) = Imu(i,j,2)
             ury(i,j,2) = Imv(i,j,2)
          end do
       end do
    end if

    ! impose lo side bc's
    select case(phys_bc(2,1))
    case (INLET)
       uly(is-1:ie+1,js,1:2) = u(is-1:ie+1,js-1,1:2)
       ury(is-1:ie+1,js,1:2) = u(is-1:ie+1,js-1,1:2)
    case (SLIP_WALL, SYMMETRY)
       uly(is-1:ie+1,js,1) = ury(is-1:ie+1,js,1)
       uly(is-1:ie+1,js,2) = ZERO
       ury(is-1:ie+1,js,2) = ZERO
    case (NO_SLIP_WALL)
       uly(is-1:ie+1,js,1:2) = ZERO
       ury(is-1:ie+1,js,1:2) = ZERO
    case (OUTLET)
       uly(is-1:ie+1,js,1) = ury(is-1:ie+1,js,1)
       uly(is-1:ie+1,js,2) = min(ury(is-1:ie+1,js,2),ZERO)
       ury(is-1:ie+1,js,2) = uly(is-1:ie+1,js,2)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_2d: invalid boundary type phys_bc(2,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(2,2))
    case (INLET)
       uly(is-1:ie+1,je+1,1:2) = u(is-1:ie+1,je+1,1:2)
       ury(is-1:ie+1,je+1,1:2) = u(is-1:ie+1,je+1,1:2)
    case (SLIP_WALL, SYMMETRY)
       ury(is-1:ie+1,je+1,1) = uly(is-1:ie+1,je+1,1)
       uly(is-1:ie+1,je+1,2) = ZERO
       ury(is-1:ie+1,je+1,2) = ZERO
    case (NO_SLIP_WALL)
       uly(is-1:ie+1,je+1,1:2) = ZERO
       ury(is-1:ie+1,je+1,1:2) = ZERO
    case (OUTLET)
       ury(is-1:ie+1,je+1,1) = uly(is-1:ie+1,je+1,1)
       uly(is-1:ie+1,je+1,2) = max(uly(is-1:ie+1,je+1,2),ZERO)
       ury(is-1:ie+1,je+1,2) = uly(is-1:ie+1,je+1,2)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_2d: invalid boundary type phys_bc(2,2)")
    end select

    do j=js,je+1
       do i=is-1,ie+1
          ! No need to compute uimhy(:,:,2) since it's equal to vtrans-w0
          ! upwind using full velocity to get transverse component of uimhy
          ! Note: vtrans already contains w0
          uimhy(i,j,1) = merge(uly(i,j,1),ury(i,j,1),vtrans(i,j).gt.ZERO)
          uavg = HALF*(uly(i,j,1)+ury(i,j,1))
          uimhy(i,j,1) = merge(uavg,uimhy(i,j,1),abs(vtrans(i,j)).lt.rel_eps)
       enddo
    enddo

    !******************************************************************
    ! Create umac and vmac
    !******************************************************************

    do j=js,je
       do i=is,ie+1
          ! use the traced force if ppm_trace_forces = 1
          fl = merge(force(i-1,j,1), Ipfx(i-1,j,1), ppm_trace_forces == 0)
          fr = merge(force(i,j  ,1), Imfx(i,  j,1), ppm_trace_forces == 0)

          ! extrapolate to edges
          umacl(i,j) = ulx(i,j,1) &
               - (dt4/hy)*(vtrans(i-1,j+1)+vtrans(i-1,j)) &
               * (uimhy(i-1,j+1,1)-uimhy(i-1,j,1)) + dt2*fl
          umacr(i,j) = urx(i,j,1) &
               - (dt4/hy)*(vtrans(i  ,j+1)+vtrans(i  ,j)) &
               * (uimhy(i  ,j+1,1)-uimhy(i  ,j,1)) + dt2*fr

          ! solve Riemann problem using full velocity
          uavg = HALF*(umacl(i,j)+umacr(i,j))
          test = ((umacl(i,j) .le. ZERO .and. umacr(i,j) .ge. ZERO) .or. &
              (abs(umacl(i,j)+umacr(i,j)) .lt. rel_eps))
          umac(i,j) = merge(umacl(i,j),umacr(i,j),uavg .gt. ZERO)
          umac(i,j) = merge(ZERO,umac(i,j),test)
       enddo
    enddo

    ! impose lo side bc's
    select case(phys_bc(1,1))
    case (INLET)
       umac(is,js:je) = u(is-1,js:je,1)
    case (SLIP_WALL, NO_SLIP_WALL, SYMMETRY)
       umac(is,js:je) = ZERO
    case (OUTLET)
       umac(is,js:je) = min(umacr(is,js:je),ZERO)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_2d: invalid boundary type phys_bc(1,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(1,2))
    case (INLET)
       umac(ie+1,js:je) = u(ie+1,js:je,1)
    case (SLIP_WALL, NO_SLIP_WALL, SYMMETRY)
       umac(ie+1,js:je) = ZERO
    case (OUTLET)
       umac(ie+1,js:je) = max(umacl(ie+1,js:je),ZERO)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_2d: invalid boundary type phys_bc(1,2)")
    end select


    do j=js,je+1
       do i=is,ie
          ! use the traced force if ppm_trace_forces = 1
          fl = merge(force(i,j-1,2), Ipfy(i,j-1,2), ppm_trace_forces == 0)
          fr = merge(force(i,j  ,2), Imfy(i,j  ,2), ppm_trace_forces == 0)
          
          ! extrapolate to edges
          vmacl(i,j) = uly(i,j,2) &
               - (dt4/hx)*(utrans(i+1,j-1)+utrans(i,j-1)) &
               * (uimhx(i+1,j-1,2)-uimhx(i,j-1,2)) + dt2*fl
          vmacr(i,j) = ury(i,j,2) &
               - (dt4/hx)*(utrans(i+1,j  )+utrans(i,j  )) &
               * (uimhx(i+1,j  ,2)-uimhx(i,j  ,2)) + dt2*fr
          
          ! solve Riemann problem using full velocity
          uavg = HALF*(vmacl(i,j)+vmacr(i,j))
          test = ((vmacl(i,j)+w0(j) .le. ZERO .and. vmacr(i,j)+w0(j) .ge. ZERO) .or. &
               (abs(vmacl(i,j)+vmacr(i,j)+TWO*w0(j)) .lt. rel_eps))
          vmac(i,j) = merge(vmacl(i,j),vmacr(i,j),uavg+w0(j) .gt. ZERO)
          vmac(i,j) = merge(ZERO,vmac(i,j),test)
       enddo
    enddo


    ! impose lo side bc's
    select case(phys_bc(2,1))
    case (INLET)
       vmac(is:ie,js) = u(is:ie,js-1,2)
    case (SLIP_WALL, NO_SLIP_WALL, SYMMETRY)
       vmac(is:ie,js) = ZERO
    case (OUTLET)
       vmac(is:ie,js) = min(vmacr(is:ie,js),ZERO)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_2d: invalid boundary type phys_bc(2,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(2,2))
    case (INLET)
       vmac(is:ie,je+1) = u(is:ie,je+1,2)
    case (SLIP_WALL, NO_SLIP_WALL, SYMMETRY)
       vmac(is:ie,je+1) = ZERO
    case (OUTLET)
       vmac(is:ie,je+1) = max(vmacl(is:ie,je+1),ZERO)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_2d: invalid boundary type phys_bc(2,2)")
    end select

    deallocate(ulx,urx,uimhx,uly,ury,uimhy,umacl,umacr,vmacl,vmacr)
    deallocate(Ipu,Imu,Ipv,Imv)

  end subroutine velpred_2d

  subroutine velpred_3d(u,ng_u,ufull,ng_uf, &
                        umac,vmac,wmac,ng_um,utrans,vtrans,wtrans,ng_ut, &
                        force,ng_f,w0,w0macx,w0macy,w0macz,ng_w0, &
                        lo,hi,dx,dt,phys_bc,adv_bc)

    use bc_module
    use slope_module
    use bl_constants_module
    use geometry, only: spherical
    use variables, only: rel_eps
    use probin_module, only: ppm_type, ppm_trace_forces
    use ppm_module

    integer        , intent(in   ) :: lo(:),hi(:)
    integer        , intent(in   ) :: ng_u,ng_uf,ng_um,ng_ut,ng_f,ng_w0
    real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :,:)
    real(kind=dp_t), intent(in   ) ::  ufull(lo(1)-ng_uf:,lo(2)-ng_uf:,lo(3)-ng_uf:,:)
    real(kind=dp_t), intent(inout) ::   umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(inout) ::   vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(inout) ::   wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) :: utrans(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(in   ) :: vtrans(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(in   ) :: wtrans(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)    
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    integer        , intent(in   ) :: phys_bc(:,:)
    integer        , intent(in   ) :: adv_bc(:,:,:)

    ! local variables
    real(kind=dp_t), allocatable :: slopex(:,:,:,:)
    real(kind=dp_t), allocatable :: slopey(:,:,:,:)
    real(kind=dp_t), allocatable :: slopez(:,:,:,:)

    real(kind=dp_t), allocatable :: Ipu(:,:,:,:), Ipfx(:,:,:,:)
    real(kind=dp_t), allocatable :: Imu(:,:,:,:), Imfx(:,:,:,:)
    real(kind=dp_t), allocatable :: Ipv(:,:,:,:), Ipfy(:,:,:,:)
    real(kind=dp_t), allocatable :: Imv(:,:,:,:), Imfy(:,:,:,:)
    real(kind=dp_t), allocatable :: Ipw(:,:,:,:), Ipfz(:,:,:,:)
    real(kind=dp_t), allocatable :: Imw(:,:,:,:), Imfz(:,:,:,:)

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

    real(kind=dp_t) :: hx, hy, hz, dt2, dt4, dt6, uavg, maxu, minu
    real(kind=dp_t) :: fl, fr

    integer :: i,j,k,is,js,ie,je,ks,ke,ung

    logical :: test

    allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

    allocate(Ipu(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(Imu(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(Ipv(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(Imv(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(Ipw(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(Imw(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

    allocate(Ipfx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(Imfx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(Ipfy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(Imfy(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(Ipfz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(Imfz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)
    ks = lo(3)
    ke = hi(3)

    dt2 = HALF*dt
    dt4 = dt/4.0d0
    dt6 = dt/6.0d0

    hx = dx(1)
    hy = dx(2)
    hz = dx(3)

    if (ppm_type .eq. 0) then
       ung = ng_u
       !$OMP PARALLEL DO PRIVATE(k) FIRSTPRIVATE(ung)
       do k = lo(3)-1,hi(3)+1
          call slopex_2d(u(:,:,k,:),slopex(:,:,k,:),lo,hi,ung,3,adv_bc)
          call slopey_2d(u(:,:,k,:),slopey(:,:,k,:),lo,hi,ung,3,adv_bc)
       end do
       !$OMP END PARALLEL DO
       call slopez_3d(u,slopez,lo,hi,ng_u,3,adv_bc)
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       call ppm_3d(u(:,:,:,1),ng_u, &
                   ufull(:,:,:,1),ufull(:,:,:,2),ufull(:,:,:,3),ng_uf, &
                   Ipu,Imu,lo,hi,adv_bc(:,:,1),dx,dt,.false.)
       call ppm_3d(u(:,:,:,2),ng_u, &
                   ufull(:,:,:,1),ufull(:,:,:,2),ufull(:,:,:,3),ng_uf, &
                   Ipv,Imv,lo,hi,adv_bc(:,:,2),dx,dt,.false.)
       call ppm_3d(u(:,:,:,3),ng_u, &
                   ufull(:,:,:,1),ufull(:,:,:,2),ufull(:,:,:,3),ng_uf, &
                   Ipw,Imw,lo,hi,adv_bc(:,:,3),dx,dt,.false.)

       ! trace forces, if necessary.  Note by default the ppm routines
       ! will trace each component to each interface in all coordinate
       ! directions, but we really only need the force traced along
       ! its respective dimension.  This should be simplified later.
       if (ppm_trace_forces == 1) then
          call ppm_3d(force(:,:,:,1),ng_u, &
                      ufull(:,:,:,1),ufull(:,:,:,2),ufull(:,:,:,3),ng_uf, &
                      Ipfx,Imfx,lo,hi,adv_bc(:,:,1),dx,dt,.false.)
          call ppm_3d(force(:,:,:,2),ng_u, &
                      ufull(:,:,:,1),ufull(:,:,:,2),ufull(:,:,:,3),ng_uf, &
                      Ipfy,Imfy,lo,hi,adv_bc(:,:,2),dx,dt,.false.)
          call ppm_3d(force(:,:,:,3),ng_u, &
                      ufull(:,:,:,1),ufull(:,:,:,2),ufull(:,:,:,3),ng_uf, &
                      Ipfz,Imfz,lo,hi,adv_bc(:,:,3),dx,dt,.false.)
       endif
    end if

    !******************************************************************
    ! Create u_{\i-\half\e_x}^x, etc.
    !******************************************************************

    ! normal predictor states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    allocate(ulx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(urx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

    if (ppm_type .eq. 0) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,maxu,minu)
       do k=ks-1,ke+1
          do j=js-1,je+1
             do i=is,ie+1
                maxu = (HALF - dt2*max(ZERO,ufull(i-1,j,k,1))/hx)
                minu = (HALF + dt2*min(ZERO,ufull(i  ,j,k,1))/hx)

                ! extrapolate all components of velocity to left face
                ulx(i,j,k,1) = u(i-1,j,k,1) + maxu * slopex(i-1,j,k,1)
                ulx(i,j,k,2) = u(i-1,j,k,2) + maxu * slopex(i-1,j,k,2)
                ulx(i,j,k,3) = u(i-1,j,k,3) + maxu * slopex(i-1,j,k,3)

                ! extrapolate all components of velocity to right face
                urx(i,j,k,1) = u(i,j,k,1) - minu * slopex(i,j,k,1)
                urx(i,j,k,2) = u(i,j,k,2) - minu * slopex(i,j,k,2)
                urx(i,j,k,3) = u(i,j,k,3) - minu * slopex(i,j,k,3)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do k=ks-1,ke+1
          do j=js-1,je+1
             do i=is,ie+1
                ! extrapolate all components of velocity to left face
                ulx(i,j,k,1) = Ipu(i-1,j,k,1)
                ulx(i,j,k,2) = Ipv(i-1,j,k,1)
                ulx(i,j,k,3) = Ipw(i-1,j,k,1)

                ! extrapolate all components of velocity to right face
                urx(i,j,k,1) = Imu(i,j,k,1)
                urx(i,j,k,2) = Imv(i,j,k,1)
                urx(i,j,k,3) = Imw(i,j,k,1)
             end do
          end do
       end do
    end if

    deallocate(slopex)

    ! impose lo side bc's
    select case(phys_bc(1,1))
    case (INLET)
       ulx(is,js-1:je+1,ks-1:ke+1,1:3) = u(is-1,js-1:je+1,ks-1:ke+1,1:3)
       urx(is,js-1:je+1,ks-1:ke+1,1:3) = u(is-1,js-1:je+1,ks-1:ke+1,1:3)
    case (SLIP_WALL, SYMMETRY)
       ulx(is,js-1:je+1,ks-1:ke+1,1) = ZERO
       urx(is,js-1:je+1,ks-1:ke+1,1) = ZERO
       ulx(is,js-1:je+1,ks-1:ke+1,2) = urx(is,js-1:je+1,ks-1:ke+1,2)
       ulx(is,js-1:je+1,ks-1:ke+1,3) = urx(is,js-1:je+1,ks-1:ke+1,3)
    case (NO_SLIP_WALL)
       ulx(is,js-1:je+1,ks-1:ke+1,1:3) = ZERO
       urx(is,js-1:je+1,ks-1:ke+1,1:3) = ZERO
    case (OUTLET)
       ulx(is,js-1:je+1,ks-1:ke+1,1) = min(urx(is,js-1:je+1,ks-1:ke+1,1),ZERO)
       urx(is,js-1:je+1,ks-1:ke+1,1) = ulx(is,js-1:je+1,ks-1:ke+1,1)
       ulx(is,js-1:je+1,ks-1:ke+1,2) = urx(is,js-1:je+1,ks-1:ke+1,2)
       ulx(is,js-1:je+1,ks-1:ke+1,3) = urx(is,js-1:je+1,ks-1:ke+1,3)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(1,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(1,2))
    case (INLET)
       ulx(ie+1,js-1:je+1,ks-1:ke+1,1:3) = u(ie+1,js-1:je+1,ks-1:ke+1,1:)
       urx(ie+1,js-1:je+1,ks-1:ke+1,1:3) = u(ie+1,js-1:je+1,ks-1:ke+1,1:3)
    case (SLIP_WALL, SYMMETRY)
       ulx(ie+1,js-1:je+1,ks-1:ke+1,1) = ZERO
       urx(ie+1,js-1:je+1,ks-1:ke+1,1) = ZERO
       urx(ie+1,js-1:je+1,ks-1:ke+1,2) = ulx(ie+1,js-1:je+1,ks-1:ke+1,2)
       urx(ie+1,js-1:je+1,ks-1:ke+1,3) = ulx(ie+1,js-1:je+1,ks-1:ke+1,3)
    case (NO_SLIP_WALL)
       ulx(ie+1,js-1:je+1,ks-1:ke+1,1:3) = ZERO
       urx(ie+1,js-1:je+1,ks-1:ke+1,1:3) = ZERO
    case (OUTLET)
       ulx(ie+1,js-1:je+1,ks-1:ke+1,1) = max(ulx(ie+1,js-1:je+1,ks-1:ke+1,1),ZERO)
       urx(ie+1,js-1:je+1,ks-1:ke+1,1) = ulx(ie+1,js-1:je+1,ks-1:ke+1,1)
       urx(ie+1,js-1:je+1,ks-1:ke+1,2) = ulx(ie+1,js-1:je+1,ks-1:ke+1,2)
       urx(ie+1,js-1:je+1,ks-1:ke+1,3) = ulx(ie+1,js-1:je+1,ks-1:ke+1,3)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(1,2)")
    end select

    allocate(uimhx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks-1,ke+1
       do j=js-1,je+1
          do i=is,ie+1
             ! No need to compute uimhx(:,:,:,1) since it's equal to utrans-w0
             ! upwind using full velocity to get transverse components of uimhx
             ! Note: utrans already contains w0
             uimhx(i,j,k,2) = merge(ulx(i,j,k,2),urx(i,j,k,2),utrans(i,j,k).gt.ZERO)
             uavg = HALF*(ulx(i,j,k,2)+urx(i,j,k,2))
             uimhx(i,j,k,2) = merge(uavg,uimhx(i,j,k,2),abs(utrans(i,j,k)).lt.rel_eps)
             
             uimhx(i,j,k,3) = merge(ulx(i,j,k,3),urx(i,j,k,3),utrans(i,j,k).gt.ZERO)
             uavg = HALF*(ulx(i,j,k,3)+urx(i,j,k,3))
             uimhx(i,j,k,3) = merge(uavg,uimhx(i,j,k,3),abs(utrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! normal predictor states

    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    allocate(uly(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(ury(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1,3))

    if (ppm_type .eq. 0) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,minu,maxu)
       do k=ks-1,ke+1
          do j=js,je+1
             do i=is-1,ie+1
                maxu = (HALF - dt2*max(ZERO,ufull(i,j-1,k,2))/hy)
                minu = (HALF + dt2*min(ZERO,ufull(i,j  ,k,2))/hy)

                ! extrapolate all components of velocity to left face
                uly(i,j,k,1) = u(i,j-1,k,1) + maxu * slopey(i,j-1,k,1)
                uly(i,j,k,2) = u(i,j-1,k,2) + maxu * slopey(i,j-1,k,2)
                uly(i,j,k,3) = u(i,j-1,k,3) + maxu * slopey(i,j-1,k,3)

                ! extrapolate all components of velocity to right face
                ury(i,j,k,1) = u(i,j,k,1) - minu * slopey(i,j,k,1)
                ury(i,j,k,2) = u(i,j,k,2) - minu * slopey(i,j,k,2)
                ury(i,j,k,3) = u(i,j,k,3) - minu * slopey(i,j,k,3)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do k=ks-1,ke+1
          do j=js,je+1
             do i=is-1,ie+1
                ! extrapolate all components of velocity to left face
                uly(i,j,k,1) = Ipu(i,j-1,k,2)
                uly(i,j,k,2) = Ipv(i,j-1,k,2)
                uly(i,j,k,3) = Ipw(i,j-1,k,2)

                ! extrapolate all components of velocity to right face
                ury(i,j,k,1) = Imu(i,j,k,2)
                ury(i,j,k,2) = Imv(i,j,k,2)
                ury(i,j,k,3) = Imw(i,j,k,2)
             enddo
          enddo
       enddo
    end if

    deallocate(slopey)

    ! impose lo side bc's
    select case(phys_bc(2,1))
    case (INLET)
       uly(is-1:ie+1,js,ks-1:ke+1,1:3) = u(is-1:ie+1,js-1,ks-1:ke+1,1:3)
       ury(is-1:ie+1,js,ks-1:ke+1,1:3) = u(is-1:ie+1,js-1,ks-1:ke+1,1:3)
    case (SLIP_WALL, SYMMETRY)
       uly(is-1:ie+1,js,ks-1:ke+1,1) = ury(is-1:ie+1,js,ks-1:ke+1,1)
       uly(is-1:ie+1,js,ks-1:ke+1,2) = ZERO
       ury(is-1:ie+1,js,ks-1:ke+1,2) = ZERO
       uly(is-1:ie+1,js,ks-1:ke+1,3) = ury(is-1:ie+1,js,ks-1:ke+1,3) 
    case (NO_SLIP_WALL)
       uly(is-1:ie+1,js,ks-1:ke+1,1:3) = ZERO
       ury(is-1:ie+1,js,ks-1:ke+1,1:3) = ZERO
    case (OUTLET)
       uly(is-1:ie+1,js,ks-1:ke+1,1) = ury(is-1:ie+1,js,ks-1:ke+1,1)
       uly(is-1:ie+1,js,ks-1:ke+1,2) = min(ury(is-1:ie+1,js,ks-1:ke+1,2),ZERO)
       ury(is-1:ie+1,js,ks-1:ke+1,2) = uly(is-1:ie+1,js,ks-1:ke+1,2)
       uly(is-1:ie+1,js,ks-1:ke+1,3) = ury(is-1:ie+1,js,ks-1:ke+1,3) 
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(2,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(2,2))
    case (INLET)
       uly(is-1:ie+1,je+1,ks-1:ke+1,1:3) = u(is-1:ie+1,je+1,ks-1:ke+1,1:3)
       ury(is-1:ie+1,je+1,ks-1:ke+1,1:3) = u(is-1:ie+1,je+1,ks-1:ke+1,1:3)
    case (SLIP_WALL, SYMMETRY)
       ury(is-1:ie+1,je+1,ks-1:ke+1,1) = uly(is-1:ie+1,je+1,ks-1:ke+1,1)
       uly(is-1:ie+1,je+1,ks-1:ke+1,2) = ZERO
       ury(is-1:ie+1,je+1,ks-1:ke+1,2) = ZERO
       ury(is-1:ie+1,je+1,ks-1:ke+1,3) = uly(is-1:ie+1,je+1,ks-1:ke+1,3)
    case (NO_SLIP_WALL)
       uly(is-1:ie+1,je+1,ks-1:ke+1,1:3) = ZERO
       ury(is-1:ie+1,je+1,ks-1:ke+1,1:3) = ZERO
    case (OUTLET)
       ury(is-1:ie+1,je+1,ks-1:ke+1,1) = uly(is-1:ie+1,je+1,ks-1:ke+1,1)
       uly(is-1:ie+1,je+1,ks-1:ke+1,2) = max(uly(is-1:ie+1,je+1,ks-1:ke+1,2),ZERO)
       ury(is-1:ie+1,je+1,ks-1:ke+1,2) = uly(is-1:ie+1,je+1,ks-1:ke+1,2)
       ury(is-1:ie+1,je+1,ks-1:ke+1,3) = uly(is-1:ie+1,je+1,ks-1:ke+1,3)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(2,2)")
    end select

    allocate(uimhy(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1,3))

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks-1,ke+1
       do j=js,je+1
          do i=is-1,ie+1
             ! No need to compute uimhy(:,:,:,2) since it's equal to vtrans-w0
             ! upwind using full velocity to get transverse components of uimhy
             ! Note: vtrans already contains w0
             uimhy(i,j,k,1) = merge(uly(i,j,k,1),ury(i,j,k,1),vtrans(i,j,k).gt.ZERO)
             uavg = HALF*(uly(i,j,k,1)+ury(i,j,k,1))
             uimhy(i,j,k,1) = merge(uavg,uimhy(i,j,k,1),abs(vtrans(i,j,k)).lt.rel_eps)
             
             uimhy(i,j,k,3) = merge(uly(i,j,k,3),ury(i,j,k,3),vtrans(i,j,k).gt.ZERO)
             uavg = HALF*(uly(i,j,k,3)+ury(i,j,k,3))
             uimhy(i,j,k,3) = merge(uavg,uimhy(i,j,k,3),abs(vtrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! normal predictor states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    allocate(ulz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1,3))
    allocate(urz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1,3))

    if (ppm_type .eq. 0) then
       !$OMP PARALLEL DO PRIVATE(i,j,k,minu,maxu)
       do k=ks,ke+1
          do j=js-1,je+1
             do i=is-1,ie+1
                maxu = (HALF - dt2*max(ZERO,ufull(i,j,k-1,3))/hz)
                minu = (HALF + dt2*min(ZERO,ufull(i,j,k  ,3))/hz)

                ! extrapolate all components of velocity to left face
                ulz(i,j,k,1) = u(i,j,k-1,1) + maxu * slopez(i,j,k-1,1)
                ulz(i,j,k,2) = u(i,j,k-1,2) + maxu * slopez(i,j,k-1,2)
                ulz(i,j,k,3) = u(i,j,k-1,3) + maxu * slopez(i,j,k-1,3)

                ! extrapolate all components of velocity to right face
                urz(i,j,k,1) = u(i,j,k,1) - minu * slopez(i,j,k,1)
                urz(i,j,k,2) = u(i,j,k,2) - minu * slopez(i,j,k,2)
                urz(i,j,k,3) = u(i,j,k,3) - minu * slopez(i,j,k,3)
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else if (ppm_type .eq. 1 .or. ppm_type .eq. 2) then
       do k=ks,ke+1
          do j=js-1,je+1
             do i=is-1,ie+1
                ! extrapolate all components of velocity to left face
                ulz(i,j,k,1) = Ipu(i,j,k-1,3)
                ulz(i,j,k,2) = Ipv(i,j,k-1,3)
                ulz(i,j,k,3) = Ipw(i,j,k-1,3)

                ! extrapolate all components of velocity to right face
                urz(i,j,k,1) = Imu(i,j,k,3)
                urz(i,j,k,2) = Imv(i,j,k,3)
                urz(i,j,k,3) = Imw(i,j,k,3)
             end do
          end do
       end do
    end if

    deallocate(slopez,Ipu,Imu,Ipv,Imv,Ipw,Imw)

    ! impose lo side bc's
    select case(phys_bc(3,1))
    case (INLET)
       ulz(is-1:ie+1,js-1:je+1,ks,1:3) = u(is-1:ie+1,js-1:je+1,ks-1,1:3)
       urz(is-1:ie+1,js-1:je+1,ks,1:3) = u(is-1:ie+1,js-1:je+1,ks-1,1:3)
    case (SLIP_WALL, SYMMETRY)
       ulz(is-1:ie+1,js-1:je+1,ks,1) = urz(is-1:ie+1,js-1:je+1,ks,1)
       ulz(is-1:ie+1,js-1:je+1,ks,2) = urz(is-1:ie+1,js-1:je+1,ks,2)
       ulz(is-1:ie+1,js-1:je+1,ks,3) = ZERO
       urz(is-1:ie+1,js-1:je+1,ks,3) = ZERO
    case (NO_SLIP_WALL)
       ulz(is-1:ie+1,js-1:je+1,ks,1:3) = ZERO
       urz(is-1:ie+1,js-1:je+1,ks,1:3) = ZERO
    case (OUTLET)
       ulz(is-1:ie+1,js-1:je+1,ks,1) = urz(is-1:ie+1,js-1:je+1,ks,1)
       ulz(is-1:ie+1,js-1:je+1,ks,2) = urz(is-1:ie+1,js-1:je+1,ks,2)
       ulz(is-1:ie+1,js-1:je+1,ks,3) = min(urz(is-1:ie+1,js-1:je+1,ks,3),ZERO)
       urz(is-1:ie+1,js-1:je+1,ks,3) = ulz(is-1:ie+1,js-1:je+1,ks,3)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(3,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(3,2))
    case (INLET)
       ulz(is-1:ie+1,js-1:je+1,ke+1,1:3) = u(is-1:ie+1,js-1:je+1,ke+1,1:3)
       urz(is-1:ie+1,js-1:je+1,ke+1,1:3) = u(is-1:ie+1,js-1:je+1,ke+1,1:3)
    case (SLIP_WALL, SYMMETRY)
       urz(is-1:ie+1,js-1:je+1,ke+1,1) = ulz(is-1:ie+1,js-1:je+1,ke+1,1)
       urz(is-1:ie+1,js-1:je+1,ke+1,2) = ulz(is-1:ie+1,js-1:je+1,ke+1,2)
       ulz(is-1:ie+1,js-1:je+1,ke+1,3) = ZERO
       urz(is-1:ie+1,js-1:je+1,ke+1,3) = ZERO
    case (NO_SLIP_WALL)
       ulz(is-1:ie+1,js-1:je+1,ke+1,1:3) = ZERO
       urz(is-1:ie+1,js-1:je+1,ke+1,1:3) = ZERO
    case (OUTLET)
       urz(is-1:ie+1,js-1:je+1,ke+1,1) = ulz(is-1:ie+1,js-1:je+1,ke+1,1)
       urz(is-1:ie+1,js-1:je+1,ke+1,2) = ulz(is-1:ie+1,js-1:je+1,ke+1,2)
       ulz(is-1:ie+1,js-1:je+1,ke+1,3) = max(ulz(is-1:ie+1,js-1:je+1,ke+1,3),ZERO)
       urz(is-1:ie+1,js-1:je+1,ke+1,3) = ulz(is-1:ie+1,js-1:je+1,ke+1,3)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(3,2)")
    end select

    allocate(uimhz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1,3))

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks,ke+1
       do j=js-1,je+1
          do i=is-1,ie+1
             ! No need to compute uimhz(:,:,:,3) since it's equal to wtrans-w0
             ! upwind using full velocity to get transverse components of uimhz
             ! Note: wtrans already contains w0
             uimhz(i,j,k,1) = merge(ulz(i,j,k,1),urz(i,j,k,1),wtrans(i,j,k).gt.ZERO)
             uavg = HALF*(ulz(i,j,k,1)+urz(i,j,k,1))
             uimhz(i,j,k,1) = merge(uavg,uimhz(i,j,k,1),abs(wtrans(i,j,k)).lt.rel_eps)
             
             uimhz(i,j,k,2) = merge(ulz(i,j,k,2),urz(i,j,k,2),wtrans(i,j,k).gt.ZERO)
             uavg = HALF*(ulz(i,j,k,2)+urz(i,j,k,2))
             uimhz(i,j,k,2) = merge(uavg,uimhz(i,j,k,2),abs(wtrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    !******************************************************************
    ! Create u_{\i-\half\e_y}^{y|z}, etc.
    !******************************************************************

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    allocate(ulyz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(uryz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(uimhyz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))

    ! uimhyz loop
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=ks,ke
       do j=js,je+1
          do i=is-1,ie+1
             ! extrapolate to faces
                ulyz(i,j,k) = uly(i,j,k,1) - (dt6/hz)*(wtrans(i,j-1,k+1)+wtrans(i,j-1,k)) &
                     * (uimhz(i,j-1,k+1,1)-uimhz(i,j-1,k,1))
                uryz(i,j,k) = ury(i,j,k,1) - (dt6/hz)*(wtrans(i,j  ,k+1)+wtrans(i,j  ,k)) &
                     * (uimhz(i,j  ,k+1,1)-uimhz(i,j  ,k,1))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! impose lo side bc's
    select case(phys_bc(2,1))
    case (INLET)
       ulyz(is-1:ie+1,js,ks:ke) = u(is-1:ie+1,js-1,ks:ke,1)
       uryz(is-1:ie+1,js,ks:ke) = u(is-1:ie+1,js-1,ks:ke,1)
    case (SLIP_WALL, SYMMETRY, OUTLET)
       ulyz(is-1:ie+1,js,ks:ke) = uryz(is-1:ie+1,js,ks:ke)
    case (NO_SLIP_WALL)
       ulyz(is-1:ie+1,js,ks:ke) = ZERO
       uryz(is-1:ie+1,js,ks:ke) = ZERO
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(2,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(2,2))
    case (INLET)
       ulyz(is-1:ie+1,je+1,ks:ke) = u(is-1:ie+1,je+1,ks:ke,1)
       uryz(is-1:ie+1,je+1,ks:ke) = u(is-1:ie+1,je+1,ks:ke,1)
    case (SLIP_WALL, SYMMETRY, OUTLET)
       uryz(is-1:ie+1,je+1,ks:ke) = ulyz(is-1:ie+1,je+1,ks:ke)
    case (NO_SLIP_WALL)
       ulyz(is-1:ie+1,je+1,ks:ke) = ZERO
       uryz(is-1:ie+1,je+1,ks:ke) = ZERO
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(2,2)")
    end select

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks,ke
       do j=js,je+1
          do i=is-1,ie+1
             ! upwind using full velocity
             uimhyz(i,j,k) = merge(ulyz(i,j,k),uryz(i,j,k),vtrans(i,j,k).gt.ZERO)
             uavg = HALF*(ulyz(i,j,k)+uryz(i,j,k))
             uimhyz(i,j,k) = merge(uavg,uimhyz(i,j,k),abs(vtrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(ulyz,uryz)

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    allocate(ulzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))
    allocate(urzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))
    allocate(uimhzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))

    ! uimhzy loop
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=ks,ke+1
       do j=js,je
          do i=is-1,ie+1
             ! extrapolate to faces
                ulzy(i,j,k) = ulz(i,j,k,1) - (dt6/hy)*(vtrans(i,j+1,k-1)+vtrans(i,j,k-1)) &
                     * (uimhy(i,j+1,k-1,1)-uimhy(i,j,k-1,1))
                urzy(i,j,k) = urz(i,j,k,1) - (dt6/hy)*(vtrans(i,j+1,k  )+vtrans(i,j,k  )) &
                     * (uimhy(i,j+1,k  ,1)-uimhy(i,j,k  ,1))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! impose lo side bc's
    select case(phys_bc(3,1))
    case (INLET)
       ulzy(is-1:ie+1,js:je,ks) = u(is-1:ie+1,js:je,ks-1,1)
       urzy(is-1:ie+1,js:je,ks) = u(is-1:ie+1,js:je,ks-1,1)
    case (SLIP_WALL, SYMMETRY, OUTLET)
       ulzy(is-1:ie+1,js:je,ks) = urzy(is-1:ie+1,js:je,ks)
    case (NO_SLIP_WALL)
       ulzy(is-1:ie+1,js:je,ks) = ZERO
       urzy(is-1:ie+1,js:je,ks) = ZERO
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(3,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(3,2))
    case (INLET)
       ulzy(is-1:ie+1,js:je,ke+1) = u(is-1:ie+1,js:je,ke+1,1)
       urzy(is-1:ie+1,js:je,ke+1) = u(is-1:ie+1,js:je,ke+1,1)
    case (SLIP_WALL, SYMMETRY, OUTLET)
       urzy(is-1:ie+1,js:je,ke+1) = ulzy(is-1:ie+1,js:je,ke+1)
    case (NO_SLIP_WALL)
       ulzy(is-1:ie+1,js:je,ke+1) = ZERO
       urzy(is-1:ie+1,js:je,ke+1) = ZERO
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(3,2)")
    end select

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks,ke+1
       do j=js,je
          do i=is-1,ie+1
             ! upwind using full velocity
             uimhzy(i,j,k) = merge(ulzy(i,j,k),urzy(i,j,k),wtrans(i,j,k).gt.ZERO)
             uavg = HALF*(ulzy(i,j,k)+urzy(i,j,k))
             uimhzy(i,j,k) = merge(uavg,uimhzy(i,j,k),abs(wtrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(ulzy,urzy)

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    allocate(vlxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))
    allocate(vrxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))

    ! vimhxz loop
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=ks,ke
       do j=js-1,je+1
          do i=is,ie+1
             ! extrapolate to faces
             vlxz(i,j,k) = ulx(i,j,k,2) - (dt6/hz)*(wtrans(i-1,j,k+1)+wtrans(i-1,j,k)) &
                  * (uimhz(i-1,j,k+1,2)-uimhz(i-1,j,k,2))
             vrxz(i,j,k) = urx(i,j,k,2) - (dt6/hz)*(wtrans(i  ,j,k+1)+wtrans(i  ,j,k)) &
                  * (uimhz(i  ,j,k+1,2)-uimhz(i  ,j,k,2))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(uimhz)

    ! impose lo side bc's
    select case(phys_bc(1,1))
    case (INLET)
       vlxz(is,js-1:je+1,ks:ke) = u(is-1,js-1:je+1,ks:ke,2)
       vrxz(is,js-1:je+1,ks:ke) = u(is-1,js-1:je+1,ks:ke,2)
    case (SLIP_WALL, SYMMETRY, OUTLET)
       vlxz(is,js-1:je+1,ks:ke) = vrxz(is,js-1:je+1,ks:ke)
    case (NO_SLIP_WALL)
       vlxz(is,js-1:je+1,ks:ke) = ZERO
       vrxz(is,js-1:je+1,ks:ke) = ZERO       
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(1,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(1,2))
    case (INLET)
       vlxz(ie+1,js-1:je+1,ks:ke) = u(ie+1,js-1:je+1,ks:ke,2)
       vrxz(ie+1,js-1:je+1,ks:ke) = u(ie+1,js-1:je+1,ks:ke,2)
    case (SLIP_WALL, SYMMETRY, OUTLET)
       vrxz(ie+1,js-1:je+1,ks:ke) = vlxz(ie+1,js-1:je+1,ks:ke)
    case (NO_SLIP_WALL)
       vlxz(ie+1,js-1:je+1,ks:ke) = ZERO
       vrxz(ie+1,js-1:je+1,ks:ke) = ZERO
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(1,2)")
    end select

    allocate(vimhxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks,ke
       do j=js-1,je+1
          do i=is,ie+1
             ! upwind using full velocity
             vimhxz(i,j,k) = merge(vlxz(i,j,k),vrxz(i,j,k),utrans(i,j,k).gt.ZERO)
             uavg = HALF*(vlxz(i,j,k)+vrxz(i,j,k))
             vimhxz(i,j,k) = merge(uavg,vimhxz(i,j,k),abs(utrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(vlxz,vrxz)

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    allocate(vlzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(vrzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(vimhzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))

    ! vimhzx loop
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=ks,ke+1
       do j=js-1,je+1
          do i=is,ie
             ! extrapolate to faces
             vlzx(i,j,k) = ulz(i,j,k,2) - (dt6/hx)*(utrans(i+1,j,k-1)+utrans(i,j,k-1)) &
                  * (uimhx(i+1,j,k-1,2)-uimhx(i,j,k-1,2))
             vrzx(i,j,k) = urz(i,j,k,2) - (dt6/hx)*(utrans(i+1,j,k  )+utrans(i,j,k  )) &
                  * (uimhx(i+1,j,k  ,2)-uimhx(i,j,k  ,2))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    ! impose lo side bc's
    select case(phys_bc(3,1))
    case (INLET)
       vlzx(is:ie,js-1:je+1,ks) = u(is:ie,js-1:je+1,ks-1,2)
       vrzx(is:ie,js-1:je+1,ks) = u(is:ie,js-1:je+1,ks-1,2)
    case (SLIP_WALL, SYMMETRY, OUTLET)
       vlzx(is:ie,js-1:je+1,ks) = vrzx(is:ie,js-1:je+1,ks)
    case (NO_SLIP_WALL)
       vlzx(is:ie,js-1:je+1,ks) = ZERO
       vrzx(is:ie,js-1:je+1,ks) = ZERO       
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(3,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(3,2))
    case (INLET)
       vlzx(is:ie,js-1:je+1,ke+1) = u(is:ie,js-1:je+1,ke+1,2)
       vrzx(is:ie,js-1:je+1,ke+1) = u(is:ie,js-1:je+1,ke+1,2)
    case (SLIP_WALL, SYMMETRY, OUTLET)
       vrzx(is:ie,js-1:je+1,ke+1) = vlzx(is:ie,js-1:je+1,ke+1)
    case (NO_SLIP_WALL)
       vlzx(is:ie,js-1:je+1,ke+1) = ZERO
       vrzx(is:ie,js-1:je+1,ke+1) = ZERO
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(3,2)")
    end select

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks,ke+1
       do j=js-1,je+1
          do i=is,ie
             ! upwind using full velocity
             vimhzx(i,j,k) = merge(vlzx(i,j,k),vrzx(i,j,k),wtrans(i,j,k).gt.ZERO)
             uavg = HALF*(vlzx(i,j,k)+vrzx(i,j,k))
             vimhzx(i,j,k) = merge(uavg,vimhzx(i,j,k),abs(wtrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(vlzx,vrzx)

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    allocate(wlxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))
    allocate(wrxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))

    ! wimhxy loop
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=ks-1,ke+1
       do j=js,je
          do i=is,ie+1
             ! extrapolate to faces
             wlxy(i,j,k) = ulx(i,j,k,3) - (dt6/hy)*(vtrans(i-1,j+1,k)+vtrans(i-1,j,k)) &
                  * (uimhy(i-1,j+1,k,3)-uimhy(i-1,j,k,3))
             wrxy(i,j,k) = urx(i,j,k,3) - (dt6/hy)*(vtrans(i  ,j+1,k)+vtrans(i  ,j,k)) &
                  * (uimhy(i  ,j+1,k,3)-uimhy(i  ,j,k,3))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(uimhy)

    ! impose lo side bc's
    select case(phys_bc(1,1))
    case (INLET)
       wlxy(is,js:je,ks-1:ke+1) = u(is-1,js:je,ks-1:ke+1,3)
       wrxy(is,js:je,ks-1:ke+1) = u(is-1,js:je,ks-1:ke+1,3)
    case (SLIP_WALL, SYMMETRY, OUTLET)
       wlxy(is,js:je,ks-1:ke+1) = wrxy(is,js:je,ks-1:ke+1)
    case (NO_SLIP_WALL)
       wlxy(is,js:je,ks-1:ke+1) = ZERO
       wrxy(is,js:je,ks-1:ke+1) = ZERO
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(1,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(1,2))
    case (INLET)
       wlxy(ie+1,js:je,ks-1:ke+1) = u(ie+1,js:je,ks-1:ke+1,3)
       wrxy(ie+1,js:je,ks-1:ke+1) = u(ie+1,js:je,ks-1:ke+1,3)
    case (SLIP_WALL, SYMMETRY, OUTLET)
       wrxy(ie+1,js:je,ks-1:ke+1) = wlxy(ie+1,js:je,ks-1:ke+1)
    case (NO_SLIP_WALL)
       wlxy(ie+1,js:je,ks-1:ke+1) = ZERO
       wrxy(ie+1,js:je,ks-1:ke+1) = ZERO
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(1,2)")
    end select

    allocate(wimhxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks-1,ke+1
       do j=js,je
          do i=is,ie+1
             ! upwind using full velocity
             wimhxy(i,j,k) = merge(wlxy(i,j,k),wrxy(i,j,k),utrans(i,j,k).gt.ZERO)
             uavg = HALF*(wlxy(i,j,k)+wrxy(i,j,k))
             wimhxy(i,j,k) = merge(uavg,wimhxy(i,j,k),abs(utrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(wlxy,wrxy)

    ! transverse states
    ! lo-1:hi+1 in base direction
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    allocate(wlyx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(wryx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))

    ! wimhyx loop
    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=ks-1,ke+1
       do j=js,je+1
          do i=is,ie
             ! extrapolate to faces
             wlyx(i,j,k) = uly(i,j,k,3) - (dt6/hx)*(utrans(i+1,j-1,k)+utrans(i,j-1,k)) &
                  * (uimhx(i+1,j-1,k,3)-uimhx(i,j-1,k,3))
             wryx(i,j,k) = ury(i,j,k,3) - (dt6/hx)*(utrans(i+1,j  ,k)+utrans(i,j  ,k)) &
                  * (uimhx(i+1,j  ,k,3)-uimhx(i,j  ,k,3))
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(uimhx)

    ! impose lo side bc's
    select case(phys_bc(2,1))
    case (INLET)
       wlyx(is:ie,js,ks-1:ke+1) = u(is:ie,js-1,ks-1:ke+1,3)
       wryx(is:ie,js,ks-1:ke+1) = u(is:ie,js-1,ks-1:ke+1,3)
    case (SLIP_WALL, SYMMETRY, OUTLET)
       wlyx(is:ie,js,ks-1:ke+1) = wryx(is:ie,js,ks-1:ke+1)
    case (NO_SLIP_WALL)
       wlyx(is:ie,js,ks-1:ke+1) = ZERO
       wryx(is:ie,js,ks-1:ke+1) = ZERO
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(2,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(2,2))
    case (INLET)
       wlyx(is:ie,je+1,ks-1:ke+1) = u(is:ie,je+1,ks-1:ke+1,3)
       wryx(is:ie,je+1,ks-1:ke+1) = u(is:ie,je+1,ks-1:ke+1,3)
    case (SLIP_WALL, SYMMETRY, OUTLET)
       wryx(is:ie,je+1,ks-1:ke+1) = wlyx(is:ie,je+1,ks-1:ke+1)
    case (NO_SLIP_WALL)
       wlyx(is:ie,je+1,ks-1:ke+1) = ZERO
       wryx(is:ie,je+1,ks-1:ke+1) = ZERO
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(2,2)")
    end select

    allocate(wimhyx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))

    !$OMP PARALLEL DO PRIVATE(i,j,k,uavg)
    do k=ks-1,ke+1
       do j=js,je+1
          do i=is,ie
             ! upwind using full velocity
             wimhyx(i,j,k) = merge(wlyx(i,j,k),wryx(i,j,k),vtrans(i,j,k).gt.ZERO)
             uavg = HALF*(wlyx(i,j,k)+wryx(i,j,k))
             wimhyx(i,j,k) = merge(uavg,wimhyx(i,j,k),abs(vtrans(i,j,k)).lt.rel_eps)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(wlyx,wryx)

    !******************************************************************
    ! Create umac, etc.
    !******************************************************************

    ! mac states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo:hi in the transverse direction
    allocate(umacl(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)))
    allocate(umacr(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)))

    !$OMP PARALLEL DO PRIVATE(i,j,k,fl,fr)
    do k=ks,ke
       do j=js,je
          do i=is,ie+1
             ! use the traced force if ppm_trace_forces = 1
             fl = merge(force(i-1,j,k,1), Ipfx(i-1,j,k,1), ppm_trace_forces == 0)
             fr = merge(force(i,j  ,k,1), Imfx(i,  j,k,1), ppm_trace_forces == 0)
             
             ! extrapolate to edges
             umacl(i,j,k) = ulx(i,j,k,1) &
                  - (dt4/hy)*(vtrans(i-1,j+1,k  )+vtrans(i-1,j,k)) &
                  * (uimhyz(i-1,j+1,k  )-uimhyz(i-1,j,k)) &
                  - (dt4/hz)*(wtrans(i-1,j  ,k+1)+wtrans(i-1,j,k)) &
                  * (uimhzy(i-1,j  ,k+1)-uimhzy(i-1,j,k)) &
                  + dt2*fl
             umacr(i,j,k) = urx(i,j,k,1) &
                  - (dt4/hy)*(vtrans(i  ,j+1,k  )+vtrans(i  ,j,k)) &
                  * (uimhyz(i  ,j+1,k  )-uimhyz(i  ,j,k)) &
                  - (dt4/hz)*(wtrans(i  ,j  ,k+1)+wtrans(i  ,j,k)) &
                  * (uimhzy(i  ,j  ,k+1)-uimhzy(i  ,j,k)) &
                  + dt2*fr
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    if (spherical .eq. 1) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke
          do j=js,je
             do i=is,ie+1
                ! solve Riemann problem using full velocity
                uavg = HALF*(umacl(i,j,k)+umacr(i,j,k))
                test = ((umacl(i,j,k)+w0macx(i,j,k) .le. ZERO .and. &
                         umacr(i,j,k)+w0macx(i,j,k) .ge. ZERO) .or. &
                    (abs(umacl(i,j,k)+umacr(i,j,k)+TWO*w0macx(i,j,k)) .lt. rel_eps))
                umac(i,j,k) = merge(umacl(i,j,k),umacr(i,j,k),uavg+w0macx(i,j,k) .gt. ZERO)
                umac(i,j,k) = merge(ZERO,umac(i,j,k),test)
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
                uavg = HALF*(umacl(i,j,k)+umacr(i,j,k))
                test = ((umacl(i,j,k) .le. ZERO .and. umacr(i,j,k) .ge. ZERO) .or. &
                    (abs(umacl(i,j,k)+umacr(i,j,k)) .lt. rel_eps))
                umac(i,j,k) = merge(umacl(i,j,k),umacr(i,j,k),uavg .gt. ZERO)
                umac(i,j,k) = merge(ZERO,umac(i,j,k),test)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

    end if

    deallocate(ulx,urx,uimhyz,uimhzy)

    ! impose lo side bc's
    select case(phys_bc(1,1))
    case (INLET)
       umac(is,js:je,ks:ke) = u(is-1,js:je,ks:ke,1)
    case (SLIP_WALL, NO_SLIP_WALL, SYMMETRY)
       umac(is,js:je,ks:ke) = ZERO
    case (OUTLET)
       umac(is,js:je,ks:ke) = min(umacr(is,js:je,ks:ke),ZERO)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(1,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(1,2))
    case (INLET)
       umac(ie+1,js:je,ks:ke) = u(ie+1,js:je,ks:ke,1)
    case (SLIP_WALL, SYMMETRY, NO_SLIP_WALL)
       umac(ie+1,js:je,ks:ke) = ZERO
    case (OUTLET)
       umac(ie+1,js:je,ks:ke) = max(umacl(ie+1,js:je,ks:ke),ZERO)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(1,2)")
    end select

    deallocate(umacl,umacr)

    ! mac states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo:hi in the transverse direction
    allocate(vmacl(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(vmacr(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3)))

    !$OMP PARALLEL DO PRIVATE(i,j,k,fl,fr)
    do k=ks,ke
       do j=js,je+1
          do i=is,ie
             ! use the traced force if ppm_trace_forces = 1
             fl = merge(force(i,j-1,k,2), Ipfy(i,j-1,k,2), ppm_trace_forces == 0)
             fr = merge(force(i,j  ,k,2), Imfy(i,j  ,k,2), ppm_trace_forces == 0)

             ! extrapolate to edges
             vmacl(i,j,k) = uly(i,j,k,2) &
                  - (dt4/hx)*(utrans(i+1,j-1,k  )+utrans(i,j-1,k)) &
                  * (vimhxz(i+1,j-1,k  )-vimhxz(i,j-1,k)) &
                  - (dt4/hz)*(wtrans(i  ,j-1,k+1)+wtrans(i,j-1,k)) &
                  * (vimhzx(i  ,j-1,k+1)-vimhzx(i,j-1,k)) &
                  + dt2*fl
             vmacr(i,j,k) = ury(i,j,k,2) &
                  - (dt4/hx)*(utrans(i+1,j  ,k  )+utrans(i,j  ,k)) &
                  * (vimhxz(i+1,j  ,k  )-vimhxz(i,j  ,k)) &
                  - (dt4/hz)*(wtrans(i  ,j  ,k+1)+wtrans(i,j  ,k)) &
                  * (vimhzx(i  ,j  ,k+1)-vimhzx(i,j  ,k)) &
                  + dt2*fr
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    if (spherical .eq. 1) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke
          do j=js,je+1
             do i=is,ie
                ! solve Riemann problem using full velocity
                uavg = HALF*(vmacl(i,j,k)+vmacr(i,j,k))
                test = ((vmacl(i,j,k)+w0macy(i,j,k) .le. ZERO .and. &
                         vmacr(i,j,k)+w0macy(i,j,k) .ge. ZERO) .or. &
                    (abs(vmacl(i,j,k)+vmacr(i,j,k)+TWO*w0macy(i,j,k)) .lt. rel_eps))
                vmac(i,j,k) = merge(vmacl(i,j,k),vmacr(i,j,k),uavg+w0macy(i,j,k) .gt. ZERO)
                vmac(i,j,k) = merge(ZERO,vmac(i,j,k),test)
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
                uavg = HALF*(vmacl(i,j,k)+vmacr(i,j,k))
                test = ((vmacl(i,j,k) .le. ZERO .and. vmacr(i,j,k) .ge. ZERO) .or. &
                    (abs(vmacl(i,j,k)+vmacr(i,j,k)) .lt. rel_eps))
                vmac(i,j,k) = merge(vmacl(i,j,k),vmacr(i,j,k),uavg .gt. ZERO)
                vmac(i,j,k) = merge(ZERO,vmac(i,j,k),test)
             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

    end if

    deallocate(uly,ury,vimhxz,vimhzx)

    ! impose lo side bc's
    select case(phys_bc(2,1))
    case (INLET)
       vmac(is:ie,js,ks:ke) = u(is:ie,js-1,ks:ke,2)
    case (SLIP_WALL, SYMMETRY, NO_SLIP_WALL)
       vmac(is:ie,js,ks:ke) = ZERO
    case (OUTLET)
       vmac(is:ie,js,ks:ke) = min(vmacr(is:ie,js,ks:ke),ZERO)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(2,1)")
    end select

    ! impose hi side bc's
    select case(phys_bc(2,2))
    case (INLET)
       vmac(is:ie,je+1,ks:ke) = u(is:ie,je+1,ks:ke,2)
    case (SLIP_WALL, SYMMETRY, NO_SLIP_WALL)
       vmac(is:ie,je+1,ks:ke) = ZERO
    case (OUTLET)
       vmac(is:ie,je+1,ks:ke) = max(vmacl(is:ie,je+1,ks:ke),ZERO)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(2,2)")
    end select

    deallocate(vmacl,vmacr)

    ! mac states
    ! Allocated from lo:hi+1 in the normal direction
    ! lo:hi in the transverse direction
    allocate(wmacl(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1))
    allocate(wmacr(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1))

    !$OMP PARALLEL DO PRIVATE(i,j,k,fl,fr)
    do k=ks,ke+1
       do j=js,je
          do i=is,ie
             ! use the traced force if ppm_trace_forces = 1
             fl = merge(force(i,j,k-1,3), Ipfz(i,j,k-1,3), ppm_trace_forces == 0)
             fr = merge(force(i,j,k  ,3), Imfz(i,j,k  ,3), ppm_trace_forces == 0)

             ! extrapolate to edges
             wmacl(i,j,k) = ulz(i,j,k,3) &
                  - (dt4/hx)*(utrans(i+1,j  ,k-1)+utrans(i,j,k-1)) &
                  * (wimhxy(i+1,j  ,k-1)-wimhxy(i,j,k-1)) &
                  - (dt4/hy)*(vtrans(i  ,j+1,k-1)+vtrans(i,j,k-1)) &
                  * (wimhyx(i  ,j+1,k-1)-wimhyx(i,j,k-1)) &
                  + dt2*fl
             wmacr(i,j,k) = urz(i,j,k,3) &
                  - (dt4/hx)*(utrans(i+1,j  ,k  )+utrans(i,j,k  )) &
                  * (wimhxy(i+1,j  ,k  )-wimhxy(i,j,k  )) &
                  - (dt4/hy)*(vtrans(i  ,j+1,k  )+vtrans(i,j,k  )) &
                  * (wimhyx(i  ,j+1,k  )-wimhyx(i,j,k  )) &
                  + dt2*fr
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    if (spherical .eq. 1) then

       !$OMP PARALLEL DO PRIVATE(i,j,k,uavg,test)
       do k=ks,ke+1
          do j=js,je
             do i=is,ie
                ! solve Riemann problem using full velocity
                uavg = HALF*(wmacl(i,j,k)+wmacr(i,j,k))
                test = ((wmacl(i,j,k)+w0macz(i,j,k) .le. ZERO .and. &
                         wmacr(i,j,k)+w0macz(i,j,k) .ge. ZERO) .or. &
                    (abs(wmacl(i,j,k)+wmacr(i,j,k)+TWO*w0macz(i,j,k)) .lt. rel_eps))
                wmac(i,j,k) = merge(wmacl(i,j,k),wmacr(i,j,k),uavg+w0macz(i,j,k) .gt. ZERO)
                wmac(i,j,k) = merge(ZERO,wmac(i,j,k),test)
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
                uavg = HALF*(wmacl(i,j,k)+wmacr(i,j,k))
                test = ((wmacl(i,j,k)+w0(k) .le. ZERO .and. &
                         wmacr(i,j,k)+w0(k) .ge. ZERO) .or. &
                    (abs(wmacl(i,j,k)+wmacr(i,j,k)+TWO*w0(k)) .lt. rel_eps))
                wmac(i,j,k) = merge(wmacl(i,j,k),wmacr(i,j,k),uavg+w0(k) .gt. ZERO)
                wmac(i,j,k) = merge(ZERO,wmac(i,j,k),test)

             enddo
          enddo
       enddo
       !$OMP END PARALLEL DO

    end if

    deallocate(ulz,urz,wimhxy,wimhyx)

    ! impose hi side bc's
    select case(phys_bc(3,1))
    case (INLET)
       wmac(is:ie,js:je,ks) = u(is:ie,js:je,ks-1,3)
    case (SLIP_WALL, SYMMETRY, NO_SLIP_WALL)
       wmac(is:ie,js:je,ks) = ZERO
    case (OUTLET)
       wmac(is:ie,js:je,ks) = min(wmacr(is:ie,js:je,ks),ZERO)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(3,1)")
    end select

    ! impose lo side bc's
    select case(phys_bc(3,2))
    case (INLET)
       wmac(is:ie,js:je,ke+1) = u(is:ie,js:je,ke+1,3)
    case (SLIP_WALL, SYMMETRY, NO_SLIP_WALL)
       wmac(is:ie,js:je,ke+1) = ZERO
    case (OUTLET)
       wmac(is:ie,js:je,ke+1) = max(wmacl(is:ie,js:je,ke+1),ZERO)
    case (INTERIOR, PERIODIC)
    case  default
       call bl_error("velpred_3d: invalid boundary type phys_bc(3,2)")
    end select

    deallocate(wmacl,wmacr)

  end subroutine velpred_3d

end module velpred_module
