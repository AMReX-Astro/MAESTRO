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
    use geometry, only: dm, nlevs

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
             call mkutrans_3d(n, up(:,:,:,:), ng_u, &
                              utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), ng_ut, &
                              w0(n,:), w0xp(:,:,:,1), w0yp(:,:,:,1), w0zp(:,:,:,1),&
                              ng_w0, lo, hi, dx(n,:), dt, &
                              the_bc_level(n)%adv_bc_level_array(i,:,:,:), &
                              the_bc_level(n)%phys_bc_level_array(i,:,:))
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

    real(kind=dp_t) hx,hy,dt2,uavg,vlo,vhi

    integer :: i,j,is,js,ie,je

    logical :: test
    
    real(kind=dp_t), allocatable :: ulx(:,:),urx(:,:)
    real(kind=dp_t), allocatable :: vly(:,:),vry(:,:)
    
    allocate(ulx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(urx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1))

    allocate(vly(lo(1)-1:hi(1)+1,lo(2):hi(2)+1))
    allocate(vry(lo(1)-1:hi(1)+1,lo(2):hi(2)+1))

    is = lo(1)
    js = lo(2)
    ie = hi(1)
    je = hi(2)
    
    dt2 = HALF*dt
    
    hx = dx(1)
    hy = dx(2)
    
    call slopex_2d(u(:,:,1:),slopex,lo,hi,ng_u,1,adv_bc(:,:,1:))
    call slopey_2d(u(:,:,2:),slopey,lo,hi,ng_u,1,adv_bc(:,:,2:))

    !******************************************************************
    ! create utrans
    !******************************************************************

    do j = js,je
       do i = is,ie+1
          ! extrapolate to edges
          ulx(i,j) = u(i-1,j,1) + (HALF - (dt2/hx)*max(ZERO,u(i-1,j,1)))*slopex(i-1,j,1)
          urx(i,j) = u(i  ,j,1) - (HALF + (dt2/hx)*min(ZERO,u(i  ,j,1)))*slopex(i  ,j,1)
       end do
    end do

    do j=js,je
       ! impose lo side bc's
       ulx(is,j) = merge(u(is-1,j,1),ulx(is,j),phys_bc(1,1) .eq. INLET)
       urx(is,j) = merge(u(is-1,j,1),urx(is,j),phys_bc(1,1) .eq. INLET)
       if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
          ulx(is,j) = ZERO
          urx(is,j) = ZERO
       endif

       ! impose hi side bc's
       ulx(ie+1,j) = merge(u(ie+1,j,1),ulx(ie+1,j),phys_bc(1,2) .eq. INLET)
       urx(ie+1,j) = merge(u(ie+1,j,1),urx(ie+1,j),phys_bc(1,2) .eq. INLET)
       if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
          ulx(ie+1,j) = ZERO
          urx(ie+1,j) = ZERO
       endif
    end do

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
          vly(i,j) = u(i,j-1,2) + (HALF - (dt2/hy)*max(ZERO,u(i,j-1,2)+vlo))*slopey(i,j-1,1)
          vry(i,j) = u(i,j  ,2) - (HALF + (dt2/hy)*min(ZERO,u(i,j  ,2)+vhi))*slopey(i,j  ,1)
       end do
    end do

    do i=is,ie
       ! impose lo side bc's
       vly(i,js) = merge(u(i,js-1,2),vly(i,js),phys_bc(2,1) .eq. INLET)
       vry(i,js) = merge(u(i,js-1,2),vry(i,js),phys_bc(2,1) .eq. INLET)
       if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
          vly(i,js) = ZERO
          vry(i,js) = ZERO
       endif

       ! impose hi side bc's
       vly(i,je+1) = merge(u(i,je+1,2),vly(i,je+1),phys_bc(2,2) .eq. INLET)
       vry(i,je+1) = merge(u(i,je+1,2),vry(i,je+1),phys_bc(2,2) .eq. INLET)
       if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
          vly(i,je+1) = ZERO
          vry(i,je+1) = ZERO
       endif
    end do

    do j=js,je+1
       do i=is,ie
          ! upwind
          uavg = HALF*(vly(i,j)+vry(i,j))
          test = ((vly(i,j)+w0(j) .le. ZERO .and. vry(i,j)+w0(j) .ge. ZERO) .or. &
               (abs(vly(i,j)+vry(i,j)+TWO*w0(j)) .lt. rel_eps))
          vtrans(i,j) = merge(vly(i,j),vry(i,j),uavg+w0(j) .gt. ZERO)
          vtrans(i,j) = merge(ZERO,vtrans(i,j),test)
       enddo
    enddo

    deallocate(ulx,urx,vly,vry)

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
    
    real(kind=dp_t) :: slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3)
    real(kind=dp_t) :: slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3)
    real(kind=dp_t) :: slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3)
    
    real(kind=dp_t) ulft,urgt,vbot,vtop,wbot,wtop
    real(kind=dp_t) uhi,ulo,vhi,vlo,whi,wlo
    real(kind=dp_t) hx, hy, hz, dt2
    
    logical :: test
    integer :: i,j,k,is,js,ks,ie,je,ke
    
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
       call slopex_2d(u(:,:,k,:),slopex(:,:,k,:),lo,hi,ng_u,3,adv_bc)
       call slopey_2d(u(:,:,k,:),slopey(:,:,k,:),lo,hi,ng_u,3,adv_bc)
    end do
    call slopez_3d(u(:,:,:,:),slopez,lo,hi,ng_u,3,adv_bc)
    
    ! Create the x-velocity to be used for transverse derivatives.
    do k = ks,ke
       do j = js,je
          do i = is,ie+1

             if (spherical .eq. 1) then
                uhi = u(i  ,j,k,1) + HALF* (w0macx(i  ,j,k)+w0macx(i+1,j,k))
                ulo = u(i-1,j,k,1) + HALF* (w0macx(i-1,j,k)+w0macx(i  ,j,k))
             else
                uhi = u(i  ,j,k,1)
                ulo = u(i-1,j,k,1)
             end if
             
             urgt = u(i,j,k  ,1) - (HALF + dt2*min(ZERO,uhi)/hx) * slopex(i  ,j,k,1)
             ulft = u(i-1,j,k,1) + (HALF - dt2*max(ZERO,ulo)/hx) * slopex(i-1,j,k,1)

             urgt = merge(u(is-1,j,k,1),urgt,i.eq.is   .and. phys_bc(1,1) .eq. INLET)
             urgt = merge(u(ie+1,j,k,1),urgt,i.eq.ie+1 .and. phys_bc(1,2) .eq. INLET)
             urgt = merge(ZERO           ,urgt,i.eq.is   .and. &
                  (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL))
             urgt = merge(ZERO           ,urgt,i.eq.ie+1 .and. &
                  (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL))
             
             ulft = merge(u(is-1,j,k,1),ulft,i.eq.is   .and. phys_bc(1,1) .eq. INLET)
             ulft = merge(u(ie+1,j,k,1),ulft,i.eq.ie+1 .and. phys_bc(1,2) .eq. INLET)
             ulft = merge(ZERO           ,ulft,i.eq.is   .and. &
                  (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL))
             ulft = merge(ZERO           ,ulft,i.eq.ie+1 .and. &
                  (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL))
             
             if (spherical .eq. 1) then
                ! upwind based on u, not utilde
                utrans(i,j,k) = merge(ulft,urgt,(ulft+urgt+TWO*w0macx(i,j,k)).gt.ZERO)
                test=( (ulft+w0macx(i,j,k).le. ZERO.and.urgt+w0macx(i,j,k).ge.ZERO) .or. &
                     (abs(ulft+urgt+TWO*w0macx(i,j,k)) .lt. rel_eps) )
                utrans(i,j,k) = merge(ZERO,utrans(i,j,k),test)
             else
                utrans(i,j,k) = merge(ulft,urgt,(ulft+urgt).gt.ZERO)
                test=( (ulft .le. ZERO  .and.  urgt .ge. ZERO)  .or. &
                     (abs(ulft+urgt) .lt. rel_eps) )
                utrans(i,j,k) = merge(ZERO,utrans(i,j,k),test)
             end if

          enddo
       enddo
    enddo
    
    ! Create the y-velocity to be used for transverse derivatives.
    do j = js,je+1
       do k = ks,ke
          do i = is,ie

             if (spherical .eq. 1) then
                vhi = u(i,j  ,k,2) + HALF* (w0macy(i,j  ,k)+w0macy(i,j+1,k))
                vlo = u(i,j-1,k,2) + HALF* (w0macy(i,j-1,k)+w0macy(i,j  ,k))
             else
                vhi = u(i,j  ,k,2)
                vlo = u(i,j-1,k,2)
             end if
             
             vtop = u(i,j  ,k,2) - (HALF + dt2*min(ZERO,vhi)/hy) * slopey(i,j  ,k,2)
             vbot = u(i,j-1,k,2) + (HALF - dt2*max(ZERO,vlo)/hy) * slopey(i,j-1,k,2)

             vtop = merge(u(i,js-1,k,2),vtop,j.eq.js   .and. phys_bc(2,1) .eq. INLET)
             vtop = merge(u(i,je+1,k,2),vtop,j.eq.je+1 .and. phys_bc(2,2) .eq. INLET)
             vtop = merge(ZERO           ,vtop,j.eq.js   .and. &
                  (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL))
             vtop = merge(ZERO           ,vtop,j.eq.je+1 .and. &
                  (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL))
             
             vbot = merge(u(i,js-1,k,2),vbot,j.eq.js   .and. phys_bc(2,1) .eq. INLET)
             vbot = merge(u(i,je+1,k,2),vbot,j.eq.je+1 .and. phys_bc(2,2) .eq. INLET)
             vbot = merge(ZERO           ,vbot,j.eq.js   .and. &
                  (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL))
             vbot = merge(ZERO           ,vbot,j.eq.je+1 .and. &
                  (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL))
             
             if (spherical .eq. 1) then
                ! upwind based on v, not vtilde
                vtrans(i,j,k)=merge(vbot,vtop,(vbot+vtop+TWO*w0macy(i,j,k)).gt.ZERO)
                test = ( (vbot+w0macy(i,j,k).le.ZERO.and.vtop+w0macy(i,j,k).ge.ZERO) .or. &
                     (abs(vbot+vtop+TWO*w0macy(i,j,k)) .lt. rel_eps) )
                vtrans(i,j,k) = merge(ZERO,vtrans(i,j,k),test)
             else
                vtrans(i,j,k)=merge(vbot,vtop,(vbot+vtop).gt.ZERO)
                test = ( (vbot .le. ZERO  .and.  vtop .ge. ZERO)  .or. &
                     (abs(vbot+vtop) .lt. rel_eps))
                vtrans(i,j,k) = merge(ZERO,vtrans(i,j,k),test)
             end if
             
          enddo
       enddo
    enddo
    
    ! Create the z-velocity to be used for transverse derivatives.
    do k = ks,ke+1
       do j = js,je
          do i = is,ie

             if (spherical .eq. 1) then
                whi = u(i,j,k  ,3) + HALF* (w0macz(i,j,k  )+w0macz(i,j,k+1))
                wlo = u(i,j,k-1,3) + HALF* (w0macz(i,j,k-1)+w0macz(i,j,k  ))
             else
                if (k .eq. nr(n)) then
                   whi = u(i,j,k,3) + w0(k)
                else
                   whi = u(i,j,k,3) + HALF* (w0(k)+w0(k+1))
                end if

                if (k .eq. ZERO) then
                   wlo = u(i,j,k-1,3) + w0(k)
                else
                   wlo = u(i,j,k-1,3) + HALF* (w0(k-1)+w0(k))
                end if
             end if

             wtop = u(i,j,k  ,3) - (HALF + dt2*min(ZERO,whi)/hz) * slopez(i,j,k  ,3)
             wbot = u(i,j,k-1,3) + (HALF - dt2*max(ZERO,wlo)/hz) * slopez(i,j,k-1,3)
             
             wtop = merge(u(i,j,ks-1,3),wtop,k.eq.ks   .and. phys_bc(3,1) .eq. INLET)
             wtop = merge(u(i,j,ke+1,3),wtop,k.eq.ke+1 .and. phys_bc(3,2) .eq. INLET)
             wtop = merge(ZERO           ,wtop,k.eq.ks   .and. &
                  (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL))
             wtop = merge(ZERO           ,wtop,k.eq.ke+1 .and. &
                  (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL))
             
             wbot = merge(u(i,j,ks-1,3),wbot,k.eq.ks   .and. phys_bc(3,1) .eq. INLET)
             wbot = merge(u(i,j,ke+1,3),wbot,k.eq.ke+1 .and. phys_bc(3,2) .eq. INLET)
             wbot = merge(ZERO           ,wbot,k.eq.ks   .and. &
                  (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL))
             wbot = merge(ZERO           ,wbot,k.eq.ke+1 .and. &
                  (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL))
             
             ! upwind based on w, not wtilde
             if (spherical .eq. 1) then
                wtrans(i,j,k)=merge(wbot,wtop,(wbot+wtop+TWO*w0macz(i,j,k)).gt.ZERO)
                test = ( (wbot+w0macz(i,j,k).le.ZERO.and.wtop+w0macz(i,j,k).ge.ZERO) .or. &
                     (abs(wbot+wtop+TWO*w0macz(i,j,k)) .lt. rel_eps) )
                wtrans(i,j,k) = merge(ZERO,wtrans(i,j,k),test)
             else
                wtrans(i,j,k)=merge(wbot,wtop,(wbot+wtop+TWO*w0(k)).gt.ZERO)
                test = ( (wbot+w0(k) .le. ZERO  .and.  wtop+w0(k) .ge. ZERO)  .or. &
                     (abs(wbot+wtop+TWO*w0(k)) .lt. rel_eps))
                wtrans(i,j,k) = merge(ZERO,wtrans(i,j,k),test)
             end if
             
          enddo
       enddo
    enddo

  end subroutine mkutrans_3d
  
end module mkutrans_module
