module velpred_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: velpred
  
contains

  subroutine velpred(s,u,umac,utrans,force,normal,w0,w0mac,dx,dt,the_bc_level,mla)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: nr_fine, dr, spherical, dm, nlevs
    use variables, only: foextrap_comp
    use fill_3d_module
    use multifab_physbc_module
    use ml_restriction_module, only : ml_edge_restriction_c

    type(multifab) , intent(in   ) :: s(:),u(:)
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(in   ) :: utrans(:,:),force(:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(ml_layout), intent(in   ) :: mla

    integer                  :: i,r,comp,n
    integer                  :: ng_s,ng_u,ng_se,ng_um,ng_ut,ng_f,ng_w0,ng_gw,ng_n
    integer                  :: lo(dm), hi(dm)
    real(kind=dp_t), pointer :: sop(:,:,:,:)
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

    ng_s = s(1)%ng
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

       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n),i) ) cycle
          sop  => dataptr(s(n),i)
          uop  => dataptr(u(n),i)
          ump  => dataptr(umac(n,1),i)
          vmp  => dataptr(umac(n,2),i)
          utp  => dataptr(utrans(n,1),i)
          vtp  => dataptr(utrans(n,2),i)
          fp   => dataptr(force(n),i)
          lo   =  lwb(get_box(s(n),i))
          hi   =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             do comp = 1,dm
                call velpred_2d(n,sop(:,:,1,:), ng_s, uop(:,:,1,:), ng_u, &
                                ump(:,:,1,1), vmp(:,:,1,1), ng_um, &
                                utp(:,:,1,1), vtp(:,:,1,1), ng_ut, &
                                fp(:,:,1,:), ng_f, w0(n,:), &
                                lo, hi, dx(n,:), dt, &
                                the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                the_bc_level(n)%adv_bc_level_array(i,:,:,comp:), &
                                comp)
             end do

          case (3)
             wmp  => dataptr(  umac(n,3),i)
             wtp  => dataptr(utrans(n,3),i)
             w0xp  => dataptr(w0mac(n,1),i)
             w0yp  => dataptr(w0mac(n,2),i)
             w0zp  => dataptr(w0mac(n,3),i)
             gw0p => dataptr(gradw0_cart,i)
             np => dataptr(normal(n),i)
             do comp = 1, dm
                call velpred_3d(n,sop(:,:,:,:), ng_s, uop(:,:,:,:), ng_u, &
                                ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), ng_ut, &
                                fp(:,:,:,:), ng_f, np(:,:,:,:), ng_n, &
                                w0(n,:),w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1), &
                                ng_w0, gw0p(:,:,:,1), ng_gw, &
                                lo, hi, dx(n,:), dt, &
                                the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                the_bc_level(n)%adv_bc_level_array(i,:,:,comp:), &
                                comp)
             end do
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

    subroutine velpred_2d(n,s,ng_s,u,ng_u,umac,vmac,ng_um,utrans, &
                          vtrans,ng_ut,force,ng_f,w0,lo,hi,dx,dt,phys_bc,adv_bc,comp)

    use geometry, only: nr
    use bc_module
    use slope_module
    use bl_constants_module
    use probin_module, only: use_new_godunov

    integer        , intent(in   ) :: n,lo(:),hi(:),ng_s,ng_u,ng_um,ng_ut,ng_f
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s :,lo(2)-ng_s :,:)
    real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u :,lo(2)-ng_u :,:)
    real(kind=dp_t), intent(inout) ::   umac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(inout) ::   vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) :: utrans(lo(1)-ng_ut:,lo(2)-ng_ut:)
    real(kind=dp_t), intent(in   ) :: vtrans(lo(1)-ng_ut:,lo(2)-ng_ut:)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,:)
    real(kind=dp_t), intent(in   ) ::     w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    integer        , intent(in   ) :: phys_bc(:,:)
    integer        , intent(in   ) :: adv_bc(:,:,:)
    integer        , intent(in   ) :: comp
    
    ! Local variables
    real(kind=dp_t) :: slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1)
    real(kind=dp_t) :: slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1)
    real(kind=dp_t) :: s_l(lo(1)-1:hi(1)+2)
    real(kind=dp_t) :: s_r(lo(1)-1:hi(1)+2)
    real(kind=dp_t) :: s_b(lo(2)-1:hi(2)+2)
    real(kind=dp_t) :: s_t(lo(2)-1:hi(2)+2)
    
    real(kind=dp_t) :: ubardth, vbardth
    real(kind=dp_t) :: hx, hy, dth
    real(kind=dp_t) :: splus,sminus
    real(kind=dp_t) :: savg,st
    real(kind=dp_t) :: vlo,vhi
    real(kind=dp_t) :: sptop,spbot,smtop,smbot,splft,sprgt,smlft,smrgt
    real(kind=dp_t) :: abs_eps, eps, umax
    
    integer :: i,j,is,js,ie,je
    logical :: test
    
    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)
    
    call slopex_2d(s(:,:,comp:),slopex,lo,hi,ng_s,1,adv_bc)
    call slopey_2d(s(:,:,comp:),slopey,lo,hi,ng_s,1,adv_bc)
    
    abs_eps = 1.0d-8
    
    dth = HALF*dt
    
    hx = dx(1)
    hy = dx(2)
    
    umax = abs(utrans(is,js))
    
    do j = js,je; do i = is,ie+1
       umax = max(umax,abs(utrans(i,j)))
    end do; end do
    do j = js,je+1; do i = is,ie
       umax = max(umax,abs(vtrans(i,j)))
    end do; end do
    
    if (umax .eq. 0.d0) then
       eps = abs_eps
    else
       eps = abs_eps * umax
    endif
    
    !********************************
    ! Loop for edge states on x-edges.
    !********************************

    if (comp .eq. 1) then
       do j = js,je 
          do i = is-1,ie+1 
             
             vlo = u(i,j,2) + HALF * (w0(j)+w0(j+1))
             if (j .ne. nr(n)-1) then
                vhi = u(i,j+1,2) + HALF * (w0(j+1)+w0(j+2))
             else
                vhi = u(i,j+1,2) + w0(j+1)
             end if
             
             if (use_new_godunov) then
                spbot = s(i,j  ,comp) + (HALF - dth*max(ZERO,vlo)/hy) * slopey(i,j  ,1)
                sptop = s(i,j+1,comp) - (HALF + dth*min(ZERO,vhi)/hy) * slopey(i,j+1,1)
             else
                spbot = s(i,j  ,comp) + (HALF - dth*vlo/hy) * slopey(i,j  ,1)
                sptop = s(i,j+1,comp) - (HALF + dth*vhi/hy) * slopey(i,j+1,1)
             end if
             
             sptop = merge(s(i,je+1,comp),sptop,j.eq.je .and. phys_bc(2,2) .eq. INLET)
             spbot = merge(s(i,je+1,comp),spbot,j.eq.je .and. phys_bc(2,2) .eq. INLET)
             
             if (j .eq. je .and. &
                  (phys_bc(2,2).eq.SLIP_WALL.or.phys_bc(2,2).eq.NO_SLIP_WALL)) then
                if (comp .eq. 2) then
                   sptop = ZERO
                   spbot = ZERO
                else
                   sptop = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                   spbot = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                endif
             endif
             
             splus = merge(spbot,sptop,vtrans(i,j+1).gt.ZERO)
             savg  = HALF * (spbot + sptop)
             splus = merge(splus, savg, abs(vtrans(i,j+1)) .gt. eps)
             
             if (j .ne. 0) then
                vlo = u(i,j-1,2) + HALF * (w0(j-1)+w0(j  ))
             else
                vlo = u(i,j-1,2) + w0(j)
             end if
             vhi = u(i,j  ,2) + HALF * (w0(j  )+w0(j+1))
             
             if (use_new_godunov) then
                smtop = s(i,j  ,comp) - (HALF + dth*min(ZERO,vhi)/hy) * slopey(i,j  ,1)
                smbot = s(i,j-1,comp) + (HALF - dth*max(ZERO,vlo)/hy) * slopey(i,j-1,1)
             else
                smtop = s(i,j  ,comp) - (HALF + dth*vhi/hy) * slopey(i,j  ,1)
                smbot = s(i,j-1,comp) + (HALF - dth*vlo/hy) * slopey(i,j-1,1)
             end if
             
             smtop = merge(s(i,js-1,comp),smtop,j.eq.js .and. phys_bc(2,1) .eq. INLET)
             smbot = merge(s(i,js-1,comp),smbot,j.eq.js .and. phys_bc(2,1) .eq. INLET)
             
             if (j .eq. js .and. &
                  (phys_bc(2,1).eq.SLIP_WALL.or.phys_bc(2,1).eq.NO_SLIP_WALL)) then
                if (comp .eq. 2) then
                   smtop = ZERO
                   smbot = ZERO
                else
                   smbot = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                   smtop = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                endif
             endif
             
             sminus = merge(smbot,smtop,vtrans(i,j).gt.ZERO)
             savg   = HALF * (smbot + smtop)
             sminus = merge(sminus, savg, abs(vtrans(i,j)) .gt. eps)
             
             st = force(i,j,comp) - HALF * (vtrans(i,j)+vtrans(i,j+1))*(splus - sminus) / hy
             
             if (comp.eq.2) then
                ! vtrans contains w0 so we need to subtract it off
                st = st - HALF * (vtrans(i,j)+vtrans(i,j+1)-w0(j+1)-w0(j))*(w0(j+1)-w0(j))/hy
             end if
             
             ubardth = dth*u(i,j,1)/hx
             
             if (use_new_godunov) then
                s_l(i+1)= s(i,j,comp) + (HALF-max(ZERO,ubardth))*slopex(i,j,1) + dth*st
                s_r(i  )= s(i,j,comp) - (HALF+min(ZERO,ubardth))*slopex(i,j,1) + dth*st
             else
                s_l(i+1)= s(i,j,comp) + (HALF-ubardth)*slopex(i,j,1) + dth*st
                s_r(i  )= s(i,j,comp) - (HALF+ubardth)*slopex(i,j,1) + dth*st
             end if
             
          enddo
          
          do i = is, ie+1 
             savg = HALF*(s_r(i) + s_l(i))
             test = ( (s_l(i) .le. ZERO  .and. s_r(i) .ge. ZERO) .or. &
                  (abs(s_l(i) + s_r(i)) .lt. eps) )
             umac(i,j)=merge(s_l(i),s_r(i),savg.gt.ZERO)
             umac(i,j)=merge(savg,umac(i,j),test)
          enddo
          
          if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
             if (comp .eq. 1) then
                umac(is,j) = ZERO
             else
                umac(is,j) = merge(ZERO,s_r(is),phys_bc(1,1).eq.NO_SLIP_WALL)
             endif
          elseif (phys_bc(1,1) .eq. INLET) then
             umac(is,j) = s(is-1,j,comp)
          elseif (phys_bc(1,1) .eq. OUTLET) then
             if (comp.eq.1) then
                umac(is,j) = MIN(s_r(is),ZERO)
             else
                umac(is,j) = s_r(is)
             end if
          endif
          if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
             if (comp .eq. 1) then
                umac(ie+1,j) = ZERO
             else
                umac(ie+1,j) = merge(ZERO,s_l(ie+1),phys_bc(1,2).eq.NO_SLIP_WALL)
             endif
          elseif (phys_bc(1,2) .eq. INLET) then
             umac(ie+1,j) = s(ie+1,j,comp)
          elseif (phys_bc(1,2) .eq. OUTLET) then
             if (comp.eq.1) then
                umac(ie+1,j) = MAX(s_l(ie+1),ZERO)
             else
                umac(ie+1,j) = s_l(ie+1)
             end if
          endif
          
       enddo
    endif
    
    !********************************
    ! Loop for edge states on y-edges.
    !********************************

    if (comp .eq. 2) then
       do i = is, ie 
          do j = js-1, je+1 
             
             if (use_new_godunov) then
                splft = s(i,j  ,comp) + (HALF - dth*max(ZERO,u(i  ,j,1))/hx)*slopex(i  ,j,1)
                sprgt = s(i+1,j,comp) - (HALF + dth*min(ZERO,u(i+1,j,1))/hx)*slopex(i+1,j,1)
             else
                splft = s(i,j  ,comp) + (HALF - dth*u(i  ,j,1)/hx) * slopex(i  ,j,1)
                sprgt = s(i+1,j,comp) - (HALF + dth*u(i+1,j,1)/hx) * slopex(i+1,j,1)
             end if
             
             sprgt = merge(s(ie+1,j,comp),sprgt,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
             splft = merge(s(ie+1,j,comp),splft,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
             
             if (i .eq. ie .and. &
                  (phys_bc(1,2).eq.SLIP_WALL.or.phys_bc(1,2).eq.NO_SLIP_WALL)) then
                if (comp .eq. 1) then
                   splft = ZERO
                   sprgt = ZERO
                else
                   sprgt = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                   splft = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                endif
             endif
             
             splus = merge(splft,sprgt,utrans(i+1,j).gt.ZERO)
             savg  = HALF * (splft + sprgt)
             splus = merge(splus, savg, abs(utrans(i+1,j)) .gt. eps)
             
             if (use_new_godunov) then
                smrgt = s(i  ,j,comp) - (HALF + dth*min(ZERO,u(i  ,j,1))/hx)*slopex(i  ,j,1)
                smlft = s(i-1,j,comp) + (HALF - dth*max(ZERO,u(i-1,j,1))/hx)*slopex(i-1,j,1)
             else
                smrgt = s(i  ,j,comp) - (HALF + dth*u(i  ,j,1)/hx) * slopex(i  ,j,1)
                smlft = s(i-1,j,comp) + (HALF - dth*u(i-1,j,1)/hx) * slopex(i-1,j,1)
             end if
             
             smrgt = merge(s(is-1,j,comp),smrgt,i.eq.is .and. phys_bc(1,1) .eq. INLET)
             smlft = merge(s(is-1,j,comp),smlft,i.eq.is .and. phys_bc(1,1) .eq. INLET)
             
             if (i .eq. is .and. &
                  (phys_bc(1,1).eq.SLIP_WALL.or.phys_bc(1,1).eq.NO_SLIP_WALL)) then
                if (comp .eq. 1) then
                   smlft = ZERO
                   smrgt = ZERO
                else
                   smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                   smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                endif
             endif
             
             sminus = merge(smlft,smrgt,utrans(i,j).gt.ZERO)
             savg   = HALF * (smlft + smrgt)
             sminus = merge(sminus, savg, abs(utrans(i,j)) .gt. eps)
             
             st = force(i,j,comp) - HALF * (utrans(i,j)+utrans(i+1,j))*(splus - sminus) / hx
             
             if (comp .eq. 2) then
                ! vtrans contains w0 so we need to subtract it off
                if (j .ge. 0 .and. j .le. nr(n)-1) then
                   st = st - HALF * &
                        (vtrans(i,j)+vtrans(i,j+1)-w0(j+1)-w0(j))*(w0(j+1)-w0(j))/hy
                end if
             end if

             if (j .ge. 0 .and. j .le. nr(n)-1) then
                vbardth = dth / hy * ( u(i,j,2) + HALF * (w0(j)+w0(j+1)) )
             else
                vbardth = dth / hy * u(i,j,2) 
             end if
             
             if (use_new_godunov) then
                s_b(j+1)= s(i,j,comp) + (HALF-max(ZERO,vbardth))*slopey(i,j,1) + dth*st
                s_t(j  )= s(i,j,comp) - (HALF+min(ZERO,vbardth))*slopey(i,j,1) + dth*st
             else
                s_b(j+1)= s(i,j,comp) + (HALF-vbardth)*slopey(i,j,1) + dth*st
                s_t(j  )= s(i,j,comp) - (HALF+vbardth)*slopey(i,j,1) + dth*st
             end if
             
          enddo
          
          do j = js, je+1 
             savg = HALF*(s_b(j) + s_t(j))
             test = ( (s_b(j) .le. ZERO .and. s_t(j) .ge. ZERO) .or. &
                  (abs(s_b(j) + s_t(j)) .lt. eps) )
             vmac(i,j)=merge(s_b(j),s_t(j),savg.gt.ZERO)
             vmac(i,j)=merge(savg,vmac(i,j),test)
          enddo
             
          if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
             if (comp .eq. 2) then
                vmac(i,js) = ZERO
             else
                vmac(i,js) = merge(ZERO,s_t(js),phys_bc(2,1).eq.NO_SLIP_WALL)
             endif
          elseif (phys_bc(2,1) .eq. INLET) then
             vmac(i,js) = s(i,js-1,comp)
          elseif (phys_bc(2,1) .eq. OUTLET) then
             if (comp.eq.2) then
                vmac(i,js) = MIN(s_t(js),ZERO)
             else
                vmac(i,js) = s_t(js)
             end if
          endif
          
          if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
             if (comp .eq. 2) then
                vmac(i,je+1) = ZERO
             else
                vmac(i,je+1) = merge(ZERO,s_b(je+1),phys_bc(2,2).eq.NO_SLIP_WALL)
             endif
          elseif (phys_bc(2,2) .eq. INLET) then
             vmac(i,je+1) = s(i,je+1,comp)
          elseif (phys_bc(2,2) .eq. OUTLET) then
             if (comp.eq.2) then
                vmac(i,je+1) = MAX(s_b(je+1),ZERO)
             else
                vmac(i,je+1) = s_b(je+1)
             end if
          endif
          
       enddo
    end if

  end subroutine velpred_2d
    
  subroutine velpred_3d(n,s,ng_s,u,ng_u, &
                        umac,vmac,wmac,ng_um,utrans,vtrans,wtrans,ng_ut, &
                        force,ng_f,normal,ng_n,w0,w0macx,w0macy,w0macz,ng_w0, &
                        gradw0_cart,ng_gw,lo,hi,dx,dt,phys_bc,adv_bc,comp)

    use bc_module
    use slope_module
    use bl_constants_module
    use geometry, only: spherical, nr
    use probin_module, only: use_new_godunov

    integer        , intent(in   ) :: n,lo(:),hi(:)
    integer        , intent(in   ) :: ng_s,ng_u,ng_um,ng_ut,ng_f,ng_n,ng_w0,ng_gw
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)
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
    integer        , intent(in   ) :: comp

    real(kind=dp_t) :: slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1)
    real(kind=dp_t) :: slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1)
    real(kind=dp_t) :: slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1)
    real(kind=dp_t) :: s_l(lo(1)-1:hi(1)+2)
    real(kind=dp_t) :: s_r(lo(1)-1:hi(1)+2)
    real(kind=dp_t) :: s_b(lo(2)-1:hi(2)+2)
    real(kind=dp_t) :: s_t(lo(2)-1:hi(2)+2)
    real(kind=dp_t) :: s_u(lo(3)-1:hi(3)+2)
    real(kind=dp_t) :: s_d(lo(3)-1:hi(3)+2)
    
    real(kind=dp_t) :: ubardth, vbardth, wbardth
    real(kind=dp_t) :: hx, hy, hz, dth, splus, sminus
    real(kind=dp_t) :: savg,st,ulo,uhi,vlo,vhi,wlo,whi
    real(kind=dp_t) :: sptop,spbot,smtop,smbot,splft,sprgt,smlft,smrgt
    real(kind=dp_t) :: abs_eps, eps, umax

    real(kind=dp_t) :: Ut_dot_er

    integer :: i,j,k,is,js,ie,je,ks,ke
    logical :: test
    
    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)
    ks = lo(3)
    ke = hi(3)
    
    do k = lo(3)-1,hi(3)+1
       call slopex_2d(s(:,:,k,comp:),slopex(:,:,k,:),lo,hi,ng_s,1,adv_bc)
       call slopey_2d(s(:,:,k,comp:),slopey(:,:,k,:),lo,hi,ng_s,1,adv_bc)
    end do
    call slopez_3d(s(:,:,:,comp:),slopez,lo,hi,ng_s,1,adv_bc)
    
    abs_eps = 1.0d-8
    
    dth = HALF*dt
    
    hx = dx(1)
    hy = dx(2)
    hz = dx(3)
    
    umax = abs(utrans(is,js,ks))
    
    do k = ks,ke; do j = js,je; do i = is,ie+1
       umax = max(umax,abs(utrans(i,j,k)))
    end do; end do; end do
    do k = ks,ke; do j = js,je+1; do i = is,ie
       umax = max(umax,abs(vtrans(i,j,k)))
    end do; end do; end do
    do k = ks,ke+1; do j = js,je; do i = is,ie
       umax = max(umax,abs(wtrans(i,j,k)))
    end do; end do; end do
    
    if (umax .eq. 0.d0) then
       eps = abs_eps
    else
       eps = abs_eps * umax
    endif
    
    !********************************
    ! Loop for edge states on x-edges.
    !********************************

     if (comp .eq. 1) then
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

                 if (use_new_godunov) then
                    spbot = s(i,j  ,k,comp) + (HALF - dth*max(ZERO,vlo)/hy)*slopey(i,j  ,k,1)
                    sptop = s(i,j+1,k,comp) - (HALF + dth*min(ZERO,vhi)/hy)*slopey(i,j+1,k,1)
                 else
                    spbot = s(i,j  ,k,comp) + (HALF - dth*vlo/hy) * slopey(i,j  ,k,1)
                    sptop = s(i,j+1,k,comp) - (HALF + dth*vhi/hy) * slopey(i,j+1,k,1)
                 end if
                 
                 sptop = merge(s(i,je+1,k,comp),sptop,j.eq.je .and. phys_bc(2,2) .eq. INLET)
                 spbot = merge(s(i,je+1,k,comp),spbot,j.eq.je .and. phys_bc(2,2) .eq. INLET)
                 
                 if (j .eq. je .and. &
                      (phys_bc(2,2).eq.SLIP_WALL.or.phys_bc(2,2).eq.NO_SLIP_WALL)) then
                    if (comp.eq.2) then
                       sptop = ZERO
                       spbot = ZERO
                    else
                       sptop = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                       spbot = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                    endif
                 endif
                 
                 splus = merge(spbot,sptop,vtrans(i,j+1,k).gt.ZERO)
                 savg  = HALF * (spbot + sptop)
                 splus = merge(splus, savg, abs(vtrans(i,j+1,k)) .gt. eps)
                 
                 if (spherical .eq. 1) then
                    vlo = u(i,j-1,k,2) + HALF*(w0macy(i,j-1,k)+w0macy(i,j  ,k))
                    vhi = u(i,j  ,k,2) + HALF*(w0macy(i,j  ,k)+w0macy(i,j+1,k))
                 else
                    vlo = u(i,j-1,k,2)
                    vhi = u(i,j  ,k,2)
                 end if

                 if (use_new_godunov) then
                    smbot = s(i,j-1,k,comp) + (HALF - dth*max(ZERO,vlo)/hy)*slopey(i,j-1,k,1)
                    smtop = s(i,j  ,k,comp) - (HALF + dth*min(ZERO,vhi)/hy)*slopey(i,j  ,k,1)
                 else
                    smbot = s(i,j-1,k,comp) + (HALF - dth*vlo/hy) * slopey(i,j-1,k,1)
                    smtop = s(i,j  ,k,comp) - (HALF + dth*vhi/hy) * slopey(i,j  ,k,1)
                 end if
                 
                 smtop = merge(s(i,js-1,k,comp),smtop,j.eq.js .and. phys_bc(2,1) .eq. INLET)
                 smbot = merge(s(i,js-1,k,comp),smbot,j.eq.js .and. phys_bc(2,1) .eq. INLET)
                 
                 if (j .eq. js .and. &
                      (phys_bc(2,1).eq.SLIP_WALL.or.phys_bc(2,1).eq.NO_SLIP_WALL)) then
                    if (comp .eq. 2) then
                       smtop = ZERO
                       smbot = ZERO
                    else
                       smbot = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                       smtop = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                    endif
                 endif
                 
                 sminus = merge(smbot,smtop,vtrans(i,j,k).gt.ZERO)
                 savg   = HALF * (smbot + smtop)
                 sminus = merge(sminus, savg, abs(vtrans(i,j,k)) .gt. eps)
                 
                 st = force(i,j,k,comp) - &
                      HALF * (vtrans(i,j,k)+vtrans(i,j+1,k))*(splus - sminus) / hy
                 
                 ! Do transverse in k direction
                 if (spherical .eq. 1) then
                    whi = u(i,j,k+1,3) + HALF*(w0macz(i,j,k+1)+w0macz(i,j,k+2))
                    wlo = u(i,j,k  ,3) + HALF*(w0macz(i,j,k  )+w0macz(i,j,k+1))
                 else

                    wlo = u(i,j,k,3) + HALF * (w0(k)+w0(k+1))
                    if (k .ne. nr(n)-1) then
                       whi = u(i,j,k+1,3) + HALF* (w0(k+1)+w0(k+2))
                    else
                       whi = u(i,j,k+1,3) + w0(k+1)
                    end if
                    
                 end if

                 if (use_new_godunov) then
                    spbot = s(i,j,k  ,comp) + (HALF - dth*max(ZERO,wlo)/hz)*slopez(i,j,k  ,1)
                    sptop = s(i,j,k+1,comp) - (HALF + dth*min(ZERO,whi)/hz)*slopez(i,j,k+1,1)
                 else
                    spbot = s(i,j,k  ,comp) + (HALF - dth*wlo/hz) * slopez(i,j,k  ,1)
                    sptop = s(i,j,k+1,comp) - (HALF + dth*whi/hz) * slopez(i,j,k+1,1)
                 end if
                 
                 sptop = merge(s(i,j,ke+1,comp),sptop,k.eq.ke .and. phys_bc(3,2) .eq. INLET)
                 spbot = merge(s(i,j,ke+1,comp),spbot,k.eq.ke .and. phys_bc(3,2) .eq. INLET)
                 
                 if (k .eq. ke .and. &
                      (phys_bc(3,2).eq.SLIP_WALL.or.phys_bc(3,2).eq.NO_SLIP_WALL)) then
                    if (comp .eq. 3) then
                       sptop = ZERO
                       spbot = ZERO
                    else
                       sptop = merge(ZERO,spbot,phys_bc(3,2).eq.NO_SLIP_WALL)
                       spbot = merge(ZERO,spbot,phys_bc(3,2).eq.NO_SLIP_WALL)
                    endif
                 endif
                 
                 splus = merge(spbot,sptop,wtrans(i,j,k+1).gt.ZERO)
                 savg  = HALF * (spbot + sptop)
                 splus = merge(splus, savg, abs(wtrans(i,j,k+1)) .gt. eps)
                 
                 if (spherical .eq. 1) then
                    whi = u(i,j,k  ,3) + HALF*(w0macz(i,j,k  )+w0macz(i,j,k+1))
                    wlo = u(i,j,k-1,3) + HALF*(w0macz(i,j,k-1)+w0macz(i,j,k  ))
                 else

                    if (k .ne. 0) then
                       wlo = u(i,j,k-1,3) + HALF* (w0(k-1)+w0(k))
                    else
                       wlo = u(i,j,k-1,3) + w0(k)
                    end if
                    whi = u(i,j,k,3) + HALF* (w0(k)+w0(k+1))

                 end if
                 
                 if (use_new_godunov) then
                    smtop = s(i,j,k  ,comp) - (HALF + dth*min(ZERO,whi)/hz)*slopez(i,j,k  ,1)
                    smbot = s(i,j,k-1,comp) + (HALF - dth*max(ZERO,wlo)/hz)*slopez(i,j,k-1,1)
                 else
                    smtop = s(i,j,k  ,comp) - (HALF + dth*whi/hz) * slopez(i,j,k  ,1)
                    smbot = s(i,j,k-1,comp) + (HALF - dth*wlo/hz) * slopez(i,j,k-1,1)
                 end if
                 
                 smtop = merge(s(i,j,ks-1,comp),smtop,k.eq.ks .and. phys_bc(3,1) .eq. INLET)
                 smbot = merge(s(i,j,ks-1,comp),smbot,k.eq.ks .and. phys_bc(3,1) .eq. INLET)
                 
                 if (k .eq. ks .and. &
                      (phys_bc(3,1).eq.SLIP_WALL.or.phys_bc(3,1).eq.NO_SLIP_WALL)) then
                    if (comp.eq.3) then
                       smtop = ZERO
                       smbot = ZERO
                    else
                       smbot = merge(ZERO,smtop,phys_bc(3,1).eq.NO_SLIP_WALL)
                       smtop = merge(ZERO,smtop,phys_bc(3,1).eq.NO_SLIP_WALL)
                    endif
                 endif
                 
                 sminus = merge(smbot,smtop,wtrans(i,j,k).gt.ZERO)
                 savg   = HALF * (smbot + smtop)
                 sminus = merge(sminus, savg, abs(wtrans(i,j,k)) .gt. eps)
                 
                 st = st - HALF * (wtrans(i,j,k)+wtrans(i,j,k+1))*(splus - sminus) / hz
                 
                 ! add the (Utilde . e_r) d w_0 /dr e_r term here
                 if (spherical .eq. 0 .and. comp.eq.3) then
                    
                    ! wtrans contains w0 so we need to subtract it off
                    st = st - HALF*(wtrans(i,j,k)+wtrans(i,j,k+1)-w0(k+1)-w0(k)) &
                         * (w0(k+1)-w0(k)) / hz
                    
                 else if (spherical .eq. 1) then
                    
                    ! u/v/wtrans contain w0, so we need to subtract it off.  
                    Ut_dot_er = (HALF*(utrans(i,j,k) + utrans(i+1,j,k)) - &
                                 HALF*(w0macx(i,j,k)+w0macx(i+1,j,k))*normal(i,j,k,1)) + &
                                 (HALF*(vtrans(i,j,k) + vtrans(i,j+1,k)) - &
                                 HALF*(w0macy(i,j,k)+w0macy(i,j+1,k))*normal(i,j,k,2)) + &
                                 (HALF*(wtrans(i,j,k) + wtrans(i,j,k+1)) - &
                                 HALF*(w0macz(i,j,k)+w0macz(i,j,k+1))*normal(i,j,k,3))
                    
                    if (comp .eq. 1) then
                       st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,1)
                    else if (comp .eq. 2) then
                       st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,2)
                    else if (comp .eq. 3) then
                       st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,3)
                    endif
                    
                 endif
                 
                 if (spherical .eq. 1) then
                    ubardth = dth/hx * ( u(i,j,k,1) + HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)) )
                 else
                    ubardth = dth/hx * u(i,j,k,1)
                 end if
                 
                 if (use_new_godunov) then
                    s_l(i+1)= s(i,j,k,comp) + (HALF-max(ZERO,ubardth))*slopex(i,j,k,1) &
                         + dth*st
                    s_r(i  )= s(i,j,k,comp) - (HALF+min(ZERO,ubardth))*slopex(i,j,k,1) &
                         + dth*st
                 else
                    s_l(i+1)= s(i,j,k,comp) + (HALF-ubardth)*slopex(i,j,k,1) + dth*st
                    s_r(i  )= s(i,j,k,comp) - (HALF+ubardth)*slopex(i,j,k,1) + dth*st
                 end if
                 
              enddo
              
              do i = is, ie+1 
                 savg = HALF*(s_r(i) + s_l(i))
                 test = ( (s_l(i) .le. ZERO .and. s_r(i) .ge. ZERO) .or. &
                      (abs(s_l(i) + s_r(i)) .lt. eps) )
                 umac(i,j,k)=merge(s_l(i),s_r(i),savg.gt.ZERO)
                 umac(i,j,k)=merge(savg,umac(i,j,k),test)
              enddo
              
              if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                 if (comp .eq. 1) then
                    umac(is,j,k) = ZERO
                 else
                    umac(is,j,k) = merge(ZERO,s_r(is),phys_bc(1,1).eq.NO_SLIP_WALL)
                 endif
              elseif (phys_bc(1,1) .eq. INLET) then
                 umac(is,j,k) = s(is-1,j,k,comp)
              elseif (phys_bc(1,1) .eq. OUTLET) then
                 if (comp.eq.1) then
                    umac(is,j,k) = MIN(s_r(is),ZERO)
                 else
                    umac(is,j,k) = s_r(is)
                 end if
              endif
              if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                 if (comp .eq. 1) then
                    umac(ie+1,j,k) = ZERO
                 else
                    umac(ie+1,j,k) = merge(ZERO,s_l(ie+1),phys_bc(1,2).eq.NO_SLIP_WALL)
                 endif
              elseif (phys_bc(1,2) .eq. INLET) then
                 umac(ie+1,j,k) = s(ie+1,j,k,comp)
              elseif (phys_bc(1,2) .eq. OUTLET) then
                 if (comp.eq.1) then
                    umac(ie+1,j,k) = MAX(s_l(ie+1),ZERO)
                 else
                    umac(ie+1,j,k) = s_l(ie+1)
                 end if
              endif
              
           enddo
        enddo
     endif

     !********************************
     ! Loop for edge states on y-edges.
     !********************************

     if (comp .eq. 2) then
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

                 if (use_new_godunov) then
                    splft = s(i  ,j,k,comp) + (HALF - dth*max(ZERO,ulo)/hx)*slopex(i  ,j,k,1)
                    sprgt = s(i+1,j,k,comp) - (HALF + dth*min(ZERO,uhi)/hx)*slopex(i+1,j,k,1)
                 else
                    splft = s(i  ,j,k,comp) + (HALF - dth*ulo/hx) * slopex(i  ,j,k,1)
                    sprgt = s(i+1,j,k,comp) - (HALF + dth*uhi/hx) * slopex(i+1,j,k,1)
                 end if
                 
                 sprgt = merge(s(ie+1,j,k,comp),sprgt,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
                 splft = merge(s(ie+1,j,k,comp),splft,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
                 
                 if (i .eq. ie .and. &
                      (phys_bc(1,2).eq.SLIP_WALL.or.phys_bc(1,2).eq.NO_SLIP_WALL)) then
                    if (comp .eq. 1) then
                       splft = ZERO
                       sprgt = ZERO
                    else
                       sprgt = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                       splft = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                    endif
                 endif
                 
                 splus = merge(splft,sprgt,utrans(i+1,j,k).gt.ZERO)
                 savg  = HALF * (splft + sprgt)
                 splus = merge(splus, savg, abs(utrans(i+1,j,k)) .gt. eps)
                 
                 if (spherical .eq. 1) then
                    ulo = u(i-1,j,k,1) + HALF*(w0macx(i-1,j,k)+w0macx(i  ,j,k))
                    uhi = u(i  ,j,k,1) + HALF*(w0macx(i  ,j,k)+w0macx(i+1,j,k))
                 else
                    ulo = u(i-1,j,k,1)
                    uhi = u(i  ,j,k,1)
                 end if

                 if (use_new_godunov) then
                    smlft = s(i-1,j,k,comp) + (HALF - dth*max(ZERO,ulo)/hx)*slopex(i-1,j,k,1)
                    smrgt = s(i  ,j,k,comp) - (HALF + dth*min(ZERO,uhi)/hx)*slopex(i  ,j,k,1)
                 else
                    smlft = s(i-1,j,k,comp) + (HALF - dth*ulo/hx) * slopex(i-1,j,k,1)
                    smrgt = s(i  ,j,k,comp) - (HALF + dth*uhi/hx) * slopex(i  ,j,k,1)
                 end if
                 
                 smrgt = merge(s(is-1,j,k,comp),smrgt,i.eq.is .and. phys_bc(1,1) .eq. INLET)
                 smlft = merge(s(is-1,j,k,comp),smlft,i.eq.is .and. phys_bc(1,1) .eq. INLET)
                 
                 if (i .eq. is .and. &
                      (phys_bc(1,1).eq.SLIP_WALL.or.phys_bc(1,1).eq.NO_SLIP_WALL)) then
                    if (comp .eq. 1) then
                       smlft = ZERO
                       smrgt = ZERO
                    else
                       smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                       smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                    endif
                 endif
                 
                 sminus = merge(smlft,smrgt,utrans(i,j,k).gt.ZERO)
                 savg   = HALF * (smlft + smrgt)
                 sminus = merge(sminus, savg, abs(utrans(i,j,k)) .gt. eps)
                 
                 st = force(i,j,k,comp) - &
                      HALF * (utrans(i,j,k)+utrans(i+1,j,k))*(splus - sminus) / hx
                 
                 ! Do transverse in k direction
                 if (spherical .eq. 1) then
                    whi = u(i,j,k+1,3) + HALF*(w0macz(i,j,k+1)+w0macz(i,j,k+2))
                    wlo = u(i,j,k  ,3) + HALF*(w0macz(i,j,k  )+w0macz(i,j,k+1))
                 else

                    wlo = u(i,j,k,3) + HALF* (w0(k)+w0(k+1))
                    if (k .ne. nr(n)-1) then
                       whi = u(i,j,k+1,3) + HALF* (w0(k+1)+w0(k+2))
                    else
                       whi = u(i,j,k+1,3) + w0(k+1)
                    end if
                    
                 end if

                 if (use_new_godunov) then
                    splft = s(i,j,k  ,comp) + (HALF - dth*max(ZERO,wlo)/hz)*slopez(i,j,k  ,1)
                    sprgt = s(i,j,k+1,comp) - (HALF + dth*min(ZERO,whi)/hz)*slopez(i,j,k+1,1)
                 else
                    splft = s(i,j,k  ,comp) + (HALF - dth*wlo/hz) * slopez(i,j,k  ,1)
                    sprgt = s(i,j,k+1,comp) - (HALF + dth*whi/hz) * slopez(i,j,k+1,1)
                 end if
                 
                 sprgt = merge(s(i,j,ke+1,comp),sprgt,k.eq.ke .and. phys_bc(3,2) .eq. INLET)
                 splft = merge(s(i,j,ke+1,comp),splft,k.eq.ke .and. phys_bc(3,2) .eq. INLET)
                 
                 if (k .eq. ke .and. &
                      (phys_bc(3,2).eq.SLIP_WALL.or.phys_bc(3,2).eq.NO_SLIP_WALL)) then
                    if (comp .eq. 3) then
                       splft = ZERO
                       sprgt = ZERO
                    else
                       sprgt = merge(ZERO,splft,phys_bc(3,2).eq.NO_SLIP_WALL)
                       splft = merge(ZERO,splft,phys_bc(3,2).eq.NO_SLIP_WALL)
                    endif
                 endif
                 
                 splus = merge(splft,sprgt,wtrans(i,j,k+1).gt.ZERO)
                 savg  = HALF * (splft + sprgt)
                 splus = merge(splus, savg, abs(wtrans(i,j,k+1)) .gt. eps)
                 
                 if (spherical .eq. 1) then
                    whi = u(i,j,k  ,3) + HALF*(w0macz(i,j,k  )+w0macz(i,j,k+1))
                    wlo = u(i,j,k-1,3) + HALF*(w0macz(i,j,k-1)+w0macz(i,j,k  ))
                 else

                    if (k .ne. 0) then
                       wlo = u(i,j,k-1,3) + HALF* (w0(k-1)+w0(k))
                    else
                       wlo = u(i,j,k-1,3) + w0(k)
                    end if
                    whi = u(i,j,k,3) + HALF* (w0(k)+w0(k+1))

                 end if
                 
                 if (use_new_godunov) then
                    smrgt = s(i,j,k  ,comp) - (HALF + dth*min(ZERO,whi)/hz)*slopez(i,j,k  ,1)
                    smlft = s(i,j,k-1,comp) + (HALF - dth*max(ZERO,wlo)/hz)*slopez(i,j,k-1,1)
                 else
                    smrgt = s(i,j,k  ,comp) - (HALF + dth*whi/hz) * slopez(i,j,k  ,1)
                    smlft = s(i,j,k-1,comp) + (HALF - dth*wlo/hz) * slopez(i,j,k-1,1)
                 end if
                 
                 smrgt = merge(s(i,j,ks-1,comp),smrgt,k.eq.ks .and. phys_bc(3,1) .eq. INLET)
                 smlft = merge(s(i,j,ks-1,comp),smlft,k.eq.ks .and. phys_bc(3,1) .eq. INLET)
                 
                 if (k .eq. ks .and. &
                      (phys_bc(3,1).eq.SLIP_WALL.or.phys_bc(3,1).eq.NO_SLIP_WALL)) then
                    if (comp .eq. 3) then
                       smlft = ZERO
                       smrgt = ZERO
                    else
                       smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                       smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                    endif
                 endif
                 
                 sminus = merge(smlft,smrgt,wtrans(i,j,k).gt.ZERO)
                 savg   = HALF * (smlft + smrgt)
                 sminus = merge(sminus, savg, abs(wtrans(i,j,k)) .gt. eps)
                 
                 st = st - HALF * (wtrans(i,j,k)+wtrans(i,j,k+1))*(splus - sminus) / hz

                 ! add the (Utilde . e_r) d w_0 /dr e_r term here
                 if (spherical .eq. 0 .and. comp.eq.3) then
                    
                    ! wtrans contains w0 so we need to subtract it off
                    st = st - HALF*(wtrans(i,j,k)+wtrans(i,j,k+1)-w0(k+1)-w0(k)) &
                         * (w0(k+1)-w0(k)) / hz
                    
                 else if (spherical .eq. 1) then
                    
                    ! u/v/wtrans contain w0, so we need to subtract it off.  
                    Ut_dot_er = (HALF*(utrans(i,j,k) + utrans(i+1,j,k)) - &
                                 HALF*(w0macx(i,j,k)+w0macx(i+1,j,k))*normal(i,j,k,1)) + &
                                 (HALF*(vtrans(i,j,k) + vtrans(i,j+1,k)) - &
                                 HALF*(w0macy(i,j,k)+w0macy(i,j+1,k))*normal(i,j,k,2)) + &
                                 (HALF*(wtrans(i,j,k) + wtrans(i,j,k+1)) - &
                                 HALF*(w0macz(i,j,k)+w0macz(i,j,k+1))*normal(i,j,k,3))
                    
                    if (comp .eq. 1) then
                       st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,1)
                    else if (comp .eq. 2) then
                       st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,2)
                    else if (comp .eq. 3) then
                       st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,3)
                    endif
                    
                 endif

                 if (spherical .eq. 1) then
                    vbardth = dth/hy * ( u(i,j,k,2) + HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)) )
                 else
                    vbardth = dth/hy * u(i,j,k,2)
                 end if
                                  
                 if (use_new_godunov) then
                    s_b(j+1)= s(i,j,k,comp) + (HALF-max(ZERO,vbardth))*slopey(i,j,k,1) &
                         + dth*st
                    s_t(j  )= s(i,j,k,comp) - (HALF+min(ZERO,vbardth))*slopey(i,j,k,1) &
                         + dth*st
                 else
                    s_b(j+1)= s(i,j,k,comp) + (HALF-vbardth)*slopey(i,j,k,1) + dth*st
                    s_t(j  )= s(i,j,k,comp) - (HALF+vbardth)*slopey(i,j,k,1) + dth*st
                 end if
                 
              enddo
              
              do j = js, je+1 
                 savg = HALF*(s_b(j) + s_t(j))
                 test = ( (s_b(j) .le. ZERO .and. s_t(j) .ge. ZERO) .or. &
                      (abs(s_b(j) + s_t(j)) .lt. eps) )
                 vmac(i,j,k)=merge(s_b(j),s_t(j),savg.gt.ZERO)
                 vmac(i,j,k)=merge(savg,vmac(i,j,k),test)
              enddo
              
              if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                 if (comp .eq. 2) then
                    vmac(i,js,k) = ZERO
                 else
                    vmac(i,js,k) = merge(ZERO,s_t(js),phys_bc(2,1).eq.NO_SLIP_WALL)
                 endif
              elseif (phys_bc(2,1) .eq. INLET) then
                 vmac(i,js,k) = s(i,js-1,k,comp)
              elseif (phys_bc(2,1) .eq. OUTLET) then
                 if (comp.eq.2) then
                    vmac(i,js,k) = MIN(s_t(js),ZERO)
                 else
                    vmac(i,js,k) = s_t(js)
                 end if
              endif
              
              if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                 if (comp .eq. 2) then
                    vmac(i,je+1,k) = ZERO
                 else
                    vmac(i,je+1,k) = merge(ZERO,s_b(je+1),phys_bc(2,2).eq.NO_SLIP_WALL)
                 endif
              elseif (phys_bc(2,2) .eq. INLET) then
                 vmac(i,je+1,k) = s(i,je+1,k,comp)
              elseif (phys_bc(2,2) .eq. OUTLET) then
                 if (comp.eq.2) then
                    vmac(i,je+1,k) = MAX(s_b(je+1),ZERO)
                 else
                    vmac(i,je+1,k) = s_b(je+1)
                 end if
              endif
              
           enddo
        enddo
     endif

     !********************************
     ! Loop for edge states on z-edges.
     !********************************

     if (comp .eq. 3) then
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
                 
                 if (use_new_godunov) then
                    splft = s(i  ,j,k,comp) + (HALF - dth*max(ZERO,ulo)/hx)*slopex(i  ,j,k,1)
                    sprgt = s(i+1,j,k,comp) - (HALF + dth*min(ZERO,uhi)/hx)*slopex(i+1,j,k,1)
                 else
                    splft = s(i  ,j,k,comp) + (HALF - dth*ulo/hx) * slopex(i  ,j,k,1)
                    sprgt = s(i+1,j,k,comp) - (HALF + dth*uhi/hx) * slopex(i+1,j,k,1)
                 end if
                 
                 sprgt = merge(s(ie+1,j,k,comp),sprgt,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
                 splft = merge(s(ie+1,j,k,comp),splft,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
                 
                 if (i .eq. ie .and. &
                      (phys_bc(1,2).eq.SLIP_WALL.or.phys_bc(1,2).eq.NO_SLIP_WALL)) then
                    if (comp .eq. 1) then
                       splft = ZERO
                       sprgt = ZERO
                    else
                       sprgt = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                       splft = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                    endif
                 endif
                 
                 splus = merge(splft,sprgt,utrans(i+1,j,k).gt.ZERO)
                 savg  = HALF * (splft + sprgt)
                 splus = merge(splus, savg, abs(utrans(i+1,j,k)) .gt. eps)
                 
                 if (spherical .eq. 1) then
                    ulo = u(i-1,j,k,1) + HALF*(w0macx(i-1,j,k)+w0macx(i  ,j,k))
                    uhi = u(i  ,j,k,1) + HALF*(w0macx(i  ,j,k)+w0macx(i+1,j,k))
                 else
                    ulo = u(i-1,j,k,1)
                    uhi = u(i  ,j,k,1)
                 end if
                 
                 if (use_new_godunov) then
                    smlft = s(i-1,j,k,comp) + (HALF - dth*max(ZERO,ulo)/hx)*slopex(i-1,j,k,1)
                    smrgt = s(i  ,j,k,comp) - (HALF + dth*min(ZERO,uhi)/hx)*slopex(i  ,j,k,1)
                 else
                    smlft = s(i-1,j,k,comp) + (HALF - dth*ulo/hx) * slopex(i-1,j,k,1)
                    smrgt = s(i  ,j,k,comp) - (HALF + dth*uhi/hx) * slopex(i  ,j,k,1)
                 end if
                 
                 smrgt = merge(s(is-1,j,k,comp),smrgt,i.eq.is .and. phys_bc(1,1) .eq. INLET)
                 smlft = merge(s(is-1,j,k,comp),smlft,i.eq.is .and. phys_bc(1,1) .eq. INLET)
                 
                 if (i .eq. is .and. &
                      (phys_bc(1,1).eq.SLIP_WALL.or.phys_bc(1,1).eq.NO_SLIP_WALL)) then
                    if (comp .eq. 1) then
                       smlft = ZERO
                       smrgt = ZERO
                    else
                       smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                       smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                    endif
                 endif
                 
                 sminus = merge(smlft,smrgt,utrans(i,j,k).gt.ZERO)
                 savg   = HALF * (smlft + smrgt)
                 sminus = merge(sminus, savg, abs(utrans(i,j,k)) .gt. eps)
                 
                 st = force(i,j,k,comp) - &
                      HALF * (utrans(i,j,k)+utrans(i+1,j,k))*(splus - sminus) / hx
                 
                 ! Do transverse in j direction
                 if (spherical .eq. 1) then
                    vlo = u(i,j  ,k,2) + HALF*(w0macy(i,j  ,k)+w0macy(i,j+1,k))
                    vhi = u(i,j+1,k,2) + HALF*(w0macy(i,j+1,k)+w0macy(i,j+2,k))
                 else
                    vlo = u(i,j  ,k,2)
                    vhi = u(i,j+1,k,2)
                 end if
                 
                 if (use_new_godunov) then
                    spbot = s(i,j  ,k,comp) + (HALF - dth*max(ZERO,vlo)/hy)*slopey(i,j  ,k,1)
                    sptop = s(i,j+1,k,comp) - (HALF + dth*min(ZERO,vhi)/hy)*slopey(i,j+1,k,1)
                 else
                    spbot = s(i,j  ,k,comp) + (HALF - dth*vlo/hy) * slopey(i,j  ,k,1)
                    sptop = s(i,j+1,k,comp) - (HALF + dth*vhi/hy) * slopey(i,j+1,k,1)
                 end if
                 
                 sptop = merge(s(i,je+1,k,comp),sptop,j.eq.je .and. phys_bc(2,2) .eq. INLET)
                 spbot = merge(s(i,je+1,k,comp),spbot,j.eq.je .and. phys_bc(2,2) .eq. INLET)
                 
                 if (j .eq. je .and. &
                      (phys_bc(2,2).eq.SLIP_WALL.or.phys_bc(2,2).eq.NO_SLIP_WALL)) then
                    if (comp.eq.2) then
                       sptop = ZERO
                       spbot = ZERO
                    else
                       sptop = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                       spbot = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                    endif
                 endif
                 
                 splus = merge(spbot,sptop,vtrans(i,j+1,k).gt.ZERO)
                 savg  = HALF * (spbot + sptop)
                 splus = merge(splus, savg, abs(vtrans(i,j+1,k)) .gt. eps)
                 
                 if (spherical .eq. 1) then
                    vlo = u(i,j-1,k,2) + HALF*(w0macy(i,j-1,k)+w0macy(i,j  ,k))
                    vhi = u(i,j  ,k,2) + HALF*(w0macy(i,j  ,k)+w0macy(i,j+1,k))
                 else
                    vlo = u(i,j-1,k,2)
                    vhi = u(i,j  ,k,2)
                 end if
                 
                 if (use_new_godunov) then
                    smbot = s(i,j-1,k,comp) + (HALF - dth*max(ZERO,vlo)/hy)*slopey(i,j-1,k,1)
                    smtop = s(i,j  ,k,comp) - (HALF + dth*min(ZERO,vhi)/hy)*slopey(i,j  ,k,1)
                 else
                    smbot = s(i,j-1,k,comp) + (HALF - dth*vlo/hy) * slopey(i,j-1,k,1)
                    smtop = s(i,j  ,k,comp) - (HALF + dth*vhi/hy) * slopey(i,j  ,k,1)
                 end if
                 
                 smtop = merge(s(i,js-1,k,comp),smtop,j.eq.js .and. phys_bc(2,1) .eq. INLET)
                 smbot = merge(s(i,js-1,k,comp),smbot,j.eq.js .and. phys_bc(2,1) .eq. INLET)
                 
                 if (j .eq. js .and. &
                      (phys_bc(2,1).eq.SLIP_WALL.or.phys_bc(2,1).eq.NO_SLIP_WALL)) then
                    if (comp .eq. 2) then
                       smtop = ZERO
                       smbot = ZERO
                    else
                       smbot = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                       smtop = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                    endif
                 endif
                 
                 sminus = merge(smbot,smtop,vtrans(i,j,k).gt.ZERO)
                 savg   = HALF * (smbot + smtop)
                 sminus = merge(sminus, savg, abs(vtrans(i,j,k)) .gt. eps)
                 
                 st = st - HALF * (vtrans(i,j,k)+vtrans(i,j+1,k))*(splus - sminus) / hy
                 
                 ! add the (Utilde . e_r) d w_0 /dr e_r term here
                 if (spherical .eq. 0 .and. comp .eq. 3) then
                    
                    if (k .ge. 0 .and. k .le. nr(n)-1) then
                       st = st - HALF*(wtrans(i,j,k)+wtrans(i,j,k+1)-w0(k+1)-w0(k)) &
                            * (w0(k+1)-w0(k)) / hz
                    end if
                    
                 else if (spherical .eq. 1) then
                    
                    ! u/v/wtrans contain w0, so we need to subtract it off.  
                    Ut_dot_er = (HALF*(utrans(i,j,k) + utrans(i+1,j,k)) - &
                                 HALF*(w0macx(i,j,k)+w0macx(i+1,j,k))*normal(i,j,k,1)) + &
                                 (HALF*(vtrans(i,j,k) + vtrans(i,j+1,k)) - &
                                 HALF*(w0macy(i,j,k)+w0macy(i,j+1,k))*normal(i,j,k,2)) + &
                                 (HALF*(wtrans(i,j,k) + wtrans(i,j,k+1)) - &
                                 HALF*(w0macz(i,j,k)+w0macz(i,j,k+1))*normal(i,j,k,3))
                    
                    if (comp .eq. 1) then
                       st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,1)
                    else if (comp .eq. 2) then
                       st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,2)
                    else if (comp .eq. 3) then
                       st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,3)
                    endif
                    
                 endif

                 if (spherical .eq. 1) then
                    wbardth = dth/hz * ( u(i,j,k,3) + HALF*(w0macz(i,j,k)+w0macz(i,j,k+1)) )
                 else
                    if (k .ge. 0 .and. k .le. nr(n)-1) then
                       wbardth = dth/hz * ( u(i,j,k,3) + HALF*(w0(k)+w0(k+1)) )
                    else
                       wbardth = dth/hz * u(i,j,k,3)
                    end if
                 end if

                 if (use_new_godunov) then
                    s_d(k+1)= s(i,j,k,comp) + (HALF-max(ZERO,wbardth))*slopez(i,j,k,1) &
                         + dth*st
                    s_u(k  )= s(i,j,k,comp) - (HALF+min(ZERO,wbardth))*slopez(i,j,k,1) &
                         + dth*st
                 else
                    s_d(k+1)= s(i,j,k,comp) + (HALF-wbardth)*slopez(i,j,k,1) + dth*st
                    s_u(k  )= s(i,j,k,comp) - (HALF+wbardth)*slopez(i,j,k,1) + dth*st
                 end if
                 
              enddo
              
              do k = ks, ke+1 
                 savg = HALF*(s_d(k) + s_u(k))
                 test = ( (s_d(k) .le. ZERO .and. s_u(k) .ge. ZERO) .or. &
                      (abs(s_d(k) + s_u(k)) .lt. eps) )
                 wmac(i,j,k)=merge(s_d(k),s_u(k),savg.gt.ZERO)
                 wmac(i,j,k)=merge(savg,wmac(i,j,k),test)
              enddo
              
              if (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                 if (comp .eq. 2) then
                    wmac(i,j,ks) = ZERO
                 else
                    wmac(i,j,ks) = merge(ZERO,s_u(ks),phys_bc(3,1).eq.NO_SLIP_WALL)
                 endif
              elseif (phys_bc(3,1) .eq. INLET) then
                 wmac(i,j,ks) = s(i,j,ks-1,comp)
              elseif (phys_bc(3,1) .eq. OUTLET) then
                 if (comp.eq.3) then
                    wmac(i,j,ks) = MIN(s_u(ks),ZERO)
                 else
                    wmac(i,j,ks) = s_u(ks)
                 end if
              endif
              
              if (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                 if (comp .eq. 2) then
                    wmac(i,j,ke+1 ) = ZERO
                 else
                    wmac(i,j,ke+1) = merge(ZERO,s_d(ke+1),phys_bc(3,2).eq.NO_SLIP_WALL)
                 endif
              elseif (phys_bc(3,2) .eq. INLET) then
                 wmac(i,j,ke+1) = s(i,j,ke+1,comp)
              elseif (phys_bc(3,2) .eq. OUTLET) then
                 if (comp.eq.3) then
                    wmac(i,j,ke+1) = MAX(s_d(ke+1),ZERO)
                 else
                    wmac(i,j,ke+1) = s_d(ke+1)
                 end if
              endif
              
           enddo
        enddo
     endif

   end subroutine velpred_3d
  
 end module velpred_module
