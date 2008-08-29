module mkutrans_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use define_bc_module

  implicit none

  private

  public :: mkutrans

contains

  subroutine mkutrans(nlevs,u,utrans,w0,w0mac,dx,dt,the_bc_level)

    use bl_prof_module
    use create_umac_grown_module
    use geometry, only: dm

    integer        , intent(in   ) :: nlevs
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
    integer                  :: lo(dm)
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
          select case (dm)
          case (2)
             call mkutrans_2d(n,up(:,:,1,:), ng_u, &
                              utp(:,:,1,1), vtp(:,:,1,1), ng_ut, w0(n,:), &
                              lo,dx(n,:),dt,&
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
                              ng_w0, lo, dx(n,:), dt, &
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

  subroutine mkutrans_2d(n,vel,ng_u,utrans,vtrans,ng_ut,w0,lo,dx,dt,adv_bc,phys_bc)

    use bc_module
    use slope_module
    use geometry, only: nr
    use probin_module, only: use_new_godunov

    integer,         intent(in   ) :: n,lo(2),ng_u,ng_ut
    real(kind=dp_t), intent(in   ) ::    vel(lo(1)-ng_u :,lo(2)-ng_u :,:)
    real(kind=dp_t), intent(inout) :: utrans(lo(1)-ng_ut:,lo(2)-ng_ut:)
    real(kind=dp_t), intent(inout) :: vtrans(lo(1)-ng_ut:,lo(2)-ng_ut:)
    real(kind=dp_t), intent(in   ) :: w0(0:)    
    real(kind=dp_t), intent(in   ) :: dt,dx(:)
    integer        , intent(in   ) :: adv_bc(:,:,:)
    integer        , intent(in   ) :: phys_bc(:,:)
    
    real(kind=dp_t), allocatable ::  velx(:,:,:), vely(:,:,:)
    
    real(kind=dp_t) hx, hy, dth, umax, dw0drhi, dw0drlo
    real(kind=dp_t) ulft,urgt,vbot,vtop,vlo,vhi,eps,abs_eps

    integer :: hi(2), i,j,is,js,ie,je
    
    logical :: test
    
    hi(1) = lo(1) + size(vel,dim=1) - (2*ng_u+1)
    hi(2) = lo(2) + size(vel,dim=2) - (2*ng_u+1)
    
    is = lo(1)
    js = lo(2)
    ie = hi(1)
    je = hi(2)
    
    abs_eps = 1.d-8
    
    ! Compute eps, which is relative to the max velocity
    umax = abs(vel(is,js,1))
    do j = js,je; do i = is,ie
       umax = max(umax,abs(vel(i,j,1)))
    end do; end do
    do j = js,je; do i = is,ie
       umax = max(umax,abs(vel(i,j,2)))
    end do; end do
    
    if (umax .eq. 0.d0) then
       eps = abs_eps
    else
       eps = abs_eps * umax
    endif
    
    dth = HALF * dt
    
    hx = dx(1)
    hy = dx(2)
    
    allocate(velx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1))
    allocate(vely(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1))
    
    call slopex_2d(vel(:,:,1:),velx,lo,ng_u,1,adv_bc)
    call slopey_2d(vel(:,:,2:),vely,lo,ng_u,1,adv_bc)
    
    ! Create the x-velocity to be used for transverse derivatives.
    do j = js,je
       do i = is,ie+1 
          
          if (use_new_godunov) then
             urgt = vel(i  ,j,1) - (HALF + dth*min(ZERO,vel(i  ,j,1))/hx) * velx(i  ,j,1)
             ulft = vel(i-1,j,1) + (HALF - dth*max(ZERO,vel(i-1,j,1))/hx) * velx(i-1,j,1)
          else
             urgt = vel(i  ,j,1) - (HALF + dth*vel(i  ,j,1)/hx) * velx(i  ,j,1)
             ulft = vel(i-1,j,1) + (HALF - dth*vel(i-1,j,1)/hx) * velx(i-1,j,1)
          end if
          
          urgt = merge(vel(is-1,j,1),urgt,i.eq.is   .and. phys_bc(1,1) .eq. INLET)
          urgt = merge(vel(ie+1,j,1),urgt,i.eq.ie+1 .and. phys_bc(1,2) .eq. INLET)
          urgt = merge(ZERO     ,urgt,i.eq.is   .and. &
               (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL))
          urgt = merge(ZERO     ,urgt,i.eq.ie+1 .and. &
               (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL))
          
          ulft = merge(vel(is-1,j,1),ulft,i.eq.is   .and. phys_bc(1,1) .eq. INLET)
          ulft = merge(vel(ie+1,j,1),ulft,i.eq.ie+1 .and. phys_bc(1,2) .eq. INLET)
          ulft = merge(ZERO     ,ulft,i.eq.is   .and. &
               (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL))
          ulft = merge(ZERO     ,ulft,i.eq.ie+1 .and. &
               (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL))
          
          utrans(i,j) = merge(ulft,urgt,(ulft+urgt).gt.ZERO)
          test = ( (ulft .le. ZERO  .and.  urgt .ge. ZERO)  .or.  &
               (abs(ulft+urgt) .lt. eps) )
          utrans(i,j) = merge(ZERO,utrans(i,j),test)
          
       enddo
    enddo
    
    ! Create the y-velocity to be used for transverse derivatives.
    do j = js,je+1 
       do i = is,ie

          if (j .eq. nr(n)) then
             vhi = vel(i,j,2) + w0(j)
             dw0drhi = (w0(j)-w0(j-1))/dx(2)
          else
             vhi = vel(i,j,2) + HALF*(w0(j+1) + w0(j))
             dw0drhi = (w0(j+1)-w0(j))/dx(2)
          end if

          if (j .eq. ZERO) then
             vlo = vel(i,j-1,2) + w0(j)
             dw0drlo = (w0(j+1)-w0(j))/dx(2)
          else
             vlo = vel(i,j-1,2) + HALF*(w0(j) + w0(j-1))
             dw0drlo = (w0(j)-w0(j-1))/dx(2)
          end if
          
          if (use_new_godunov) then
             vtop = vel(i,j  ,2) - (HALF + dth*min(ZERO,vhi)/hy) * vely(i,j  ,1) &
                  - dth*vel(i,j  ,2)*dw0drhi
             vbot = vel(i,j-1,2) + (HALF - dth*max(ZERO,vlo)/hy) * vely(i,j-1,1) &
                  - dth*vel(i,j-1,2)*dw0drlo
          else
             vtop = vel(i,j  ,2) - (HALF + dth*vhi/hy) * vely(i,j  ,1)
             vbot = vel(i,j-1,2) + (HALF - dth*vlo/hy) * vely(i,j-1,1)
          end if

          vtop = merge(vel(i,js-1,2),vtop,j.eq.js   .and. phys_bc(2,1) .eq. INLET)
          vtop = merge(vel(i,je+1,2),vtop,j.eq.je+1 .and. phys_bc(2,2) .eq. INLET)
          vtop = merge(ZERO     ,vtop,j.eq.js   .and. &
               (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL))
          vtop = merge(ZERO     ,vtop,j.eq.je+1 .and. &
               (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL))
          
          vbot = merge(vel(i,js-1,2),vbot,j.eq.js   .and. phys_bc(2,1) .eq. INLET)
          vbot = merge(vel(i,je+1,2),vbot,j.eq.je+1 .and. phys_bc(2,2) .eq. INLET)
          vbot = merge(ZERO     ,vbot,j.eq.js   .and. &
                       (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL))
          vbot = merge(ZERO     ,vbot,j.eq.je+1 .and. &
                       (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL))
          
          vtrans(i,j)=merge(vbot,vtop,(vbot+vtop).gt.ZERO)
          test = ( (vbot .le. ZERO  .and.  vtop .ge. ZERO)  .or.  &
               (abs(vbot+vtop) .lt. eps))
          vtrans(i,j) = merge(ZERO,vtrans(i,j),test)
       enddo
    enddo

  end subroutine mkutrans_2d
  
  subroutine mkutrans_3d(n,vel,ng_u,utrans,vtrans,wtrans,ng_ut,w0,w0macx,w0macy,w0macz, &
                         ng_w0,lo,dx,dt,adv_bc,phys_bc)

    use bc_module
    use slope_module
    use geometry, only: nr, spherical
    
    integer, intent(in)            :: n,lo(3),ng_u,ng_ut,ng_w0    
    real(kind=dp_t), intent(in   ) ::         vel(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :,:)
    real(kind=dp_t), intent(inout) ::      utrans(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(inout) ::      vtrans(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(inout) ::      wtrans(lo(1)-ng_ut:,lo(2)-ng_ut:,lo(3)-ng_ut:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: dt,dx(:)
    integer        , intent(in   ) :: adv_bc(:,:,:)
    integer        , intent(in   ) :: phys_bc(:,:)
    
    real(kind=dp_t), allocatable::  velx(:,:,:,:),vely(:,:,:,:),velz(:,:,:,:)
    
    real(kind=dp_t) ulft,urgt,vbot,vtop,wbot,wtop
    real(kind=dp_t) uhi,ulo,vhi,vlo,whi,wlo
    real(kind=dp_t) hx, hy, hz, dth, umax, eps, abs_eps
    
    logical :: test
    
    integer :: hi(3)
    integer :: i,j,k,is,js,ks,ie,je,ke
    
    hi(1) = lo(1) + size(vel,dim=1) - (2*ng_u+1)
    hi(2) = lo(2) + size(vel,dim=2) - (2*ng_u+1)
    hi(3) = lo(3) + size(vel,dim=3) - (2*ng_u+1)
    
    is = lo(1)
    js = lo(2)
    ks = lo(3)
    ie = hi(1)
    je = hi(2)
    ke = hi(3)
    
    abs_eps = 1.d-8
    
    ! Compute eps, which is relative to the max velocity
    umax = abs(vel(is,js,ks,1))
    do k = ks,ke; do j = js,je; do i = is,ie
       umax = max(umax,abs(vel(i,j,k,1)))
    end do; end do; end do
    do k = ks,ke; do j = js,je; do i = is,ie
       umax = max(umax,abs(vel(i,j,k,2)))
    end do; end do; end do
    do k = ks,ke; do j = js,je; do i = is,ie
       umax = max(umax,abs(vel(i,j,k,3)))
    end do; end do; end do
    
    if (umax .eq. 0.d0) then
       eps = abs_eps
    else
       eps = abs_eps * umax
    endif
    
    dth = HALF * dt
    
    hx = dx(1)
    hy = dx(2)
    hz = dx(3)
    
    allocate(velx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))
    allocate(vely(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))
    allocate(velz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))
    
    do k = lo(3)-1,hi(3)+1
       call slopex_2d(vel(:,:,k,1:),velx(:,:,k,:),lo,ng_u,1,adv_bc)
       call slopey_2d(vel(:,:,k,2:),vely(:,:,k,:),lo,ng_u,1,adv_bc)
    end do
    call slopez_3d(vel(:,:,:,3:),velz,lo,ng_u,1,adv_bc)
    
    ! Create the x-velocity to be used for transverse derivatives.
    do k = ks,ke
       do j = js,je
          do i = is,ie+1

             if (spherical .eq. 1) then
                uhi = vel(i  ,j,k,1) + HALF* (w0macx(i  ,j,k)+w0macx(i+1,j,k))
                ulo = vel(i-1,j,k,1) + HALF* (w0macx(i-1,j,k)+w0macx(i  ,j,k))
             else
                uhi = vel(i  ,j,k,1)
                ulo = vel(i-1,j,k,1)
             end if
             
             urgt = vel(i,j,k  ,1) - (HALF + dth*uhi/hx) * velx(i  ,j,k,1)
             ulft = vel(i-1,j,k,1) + (HALF - dth*ulo/hx) * velx(i-1,j,k,1)

             urgt = merge(vel(is-1,j,k,1),urgt,i.eq.is   .and. phys_bc(1,1) .eq. INLET)
             urgt = merge(vel(ie+1,j,k,1),urgt,i.eq.ie+1 .and. phys_bc(1,2) .eq. INLET)
             urgt = merge(ZERO           ,urgt,i.eq.is   .and. &
                  (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL))
             urgt = merge(ZERO           ,urgt,i.eq.ie+1 .and. &
                  (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL))
             
             ulft = merge(vel(is-1,j,k,1),ulft,i.eq.is   .and. phys_bc(1,1) .eq. INLET)
             ulft = merge(vel(ie+1,j,k,1),ulft,i.eq.ie+1 .and. phys_bc(1,2) .eq. INLET)
             ulft = merge(ZERO           ,ulft,i.eq.is   .and. &
                  (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL))
             ulft = merge(ZERO           ,ulft,i.eq.ie+1 .and. &
                  (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL))
             
             utrans(i,j,k) = merge(ulft,urgt,(ulft+urgt).gt.ZERO)
             test=( (ulft .le. ZERO  .and.  urgt .ge. ZERO)  .or. &
                  (abs(ulft+urgt) .lt. eps) )
             utrans(i,j,k) = merge(ZERO,utrans(i,j,k),test)
          enddo
       enddo
    enddo
    
    ! Create the y-velocity to be used for transverse derivatives.
    do j = js,je+1
       do k = ks,ke
          do i = is,ie

             if (spherical .eq. 1) then
                vhi = vel(i,j  ,k,2) + HALF* (w0macy(i,j  ,k)+w0macy(i,j+1,k))
                vlo = vel(i,j-1,k,2) + HALF* (w0macy(i,j-1,k)+w0macy(i,j  ,k))
             else
                vhi = vel(i,j  ,k,2)
                vlo = vel(i,j-1,k,2)
             end if
             
             vtop = vel(i,j  ,k,2) - (HALF + dth*vhi/hy) * vely(i,j  ,k,1)
             vbot = vel(i,j-1,k,2) + (HALF - dth*vlo/hy) * vely(i,j-1,k,1)

             vtop = merge(vel(i,js-1,k,2),vtop,j.eq.js   .and. phys_bc(2,1) .eq. INLET)
             vtop = merge(vel(i,je+1,k,2),vtop,j.eq.je+1 .and. phys_bc(2,2) .eq. INLET)
             vtop = merge(ZERO           ,vtop,j.eq.js   .and. &
                  (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL))
             vtop = merge(ZERO           ,vtop,j.eq.je+1 .and. &
                  (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL))
             
             vbot = merge(vel(i,js-1,k,2),vbot,j.eq.js   .and. phys_bc(2,1) .eq. INLET)
             vbot = merge(vel(i,je+1,k,2),vbot,j.eq.je+1 .and. phys_bc(2,2) .eq. INLET)
             vbot = merge(ZERO           ,vbot,j.eq.js   .and. &
                  (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL))
             vbot = merge(ZERO           ,vbot,j.eq.je+1 .and. &
                  (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL))
             
             vtrans(i,j,k)=merge(vbot,vtop,(vbot+vtop).gt.ZERO)
             test = ( (vbot .le. ZERO  .and.  vtop .ge. ZERO)  .or. &
                  (abs(vbot+vtop) .lt. eps))
             vtrans(i,j,k) = merge(ZERO,vtrans(i,j,k),test)
             
          enddo
       enddo
    enddo
    
    ! Create the z-velocity to be used for transverse derivatives.
    do k = ks,ke+1
       do j = js,je
          do i = is,ie

             if (spherical .eq. 1) then
                whi = vel(i,j,k  ,3) + HALF* (w0macz(i,j,k  )+w0macz(i,j,k+1))
                wlo = vel(i,j,k-1,3) + HALF* (w0macz(i,j,k-1)+w0macz(i,j,k  ))
             else

                if (k .eq. nr(n)) then
                   whi = vel(i,j,k,3) + w0(k)
                else
                   whi = vel(i,j,k,3) + HALF* (w0(k)+w0(k+1))
                end if

                if (k .eq. ZERO) then
                   wlo = vel(i,j,k-1,3) + w0(k)
                else
                   wlo = vel(i,j,k-1,3) + HALF* (w0(k-1)+w0(k))
                end if

             end if

             wtop = vel(i,j,k  ,3) - (HALF + dth*whi/hz) * velz(i,j,k  ,1)
             wbot = vel(i,j,k-1,3) + (HALF - dth*wlo/hz) * velz(i,j,k-1,1)
             
             wtop = merge(vel(i,j,ks-1,3),wtop,k.eq.ks   .and. phys_bc(3,1) .eq. INLET)
             wtop = merge(vel(i,j,ke+1,3),wtop,k.eq.ke+1 .and. phys_bc(3,2) .eq. INLET)
             wtop = merge(ZERO           ,wtop,k.eq.ks   .and. &
                  (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL))
             wtop = merge(ZERO           ,wtop,k.eq.ke+1 .and. &
                  (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL))
             
             wbot = merge(vel(i,j,ks-1,3),wbot,k.eq.ks   .and. phys_bc(3,1) .eq. INLET)
             wbot = merge(vel(i,j,ke+1,3),wbot,k.eq.ke+1 .and. phys_bc(3,2) .eq. INLET)
             wbot = merge(ZERO           ,wbot,k.eq.ks   .and. &
                  (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL))
             wbot = merge(ZERO           ,wbot,k.eq.ke+1 .and. &
                  (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL))
             
             wtrans(i,j,k)=merge(wbot,wtop,(wbot+wtop).gt.ZERO)
             test = ( (wbot .le. ZERO  .and.  wtop .ge. ZERO)  .or. &
                  (abs(wbot+wtop) .lt. eps))
             wtrans(i,j,k) = merge(ZERO,wtrans(i,j,k),test)
             
          enddo
       enddo
    enddo

  end subroutine mkutrans_3d
  
end module mkutrans_module
