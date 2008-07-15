! make_edge_state constructs the edge state of a variable, using a 
! second-order Taylor expansion in space (through dx/2) and time 
! (though dt/2). 
!
! If is_vel = .true., then we are computing the edge states for the
! velocity.
!
! If velpred = 1 (then is_vel should also be 1), and we are computing 
! the edge states for the MAC.  In this case, we only need the normal 
! edge states (i.e. the x-edge state for u, the y-edge state for v, ...)
!
! If velpred = 0, then we are computing all edge states for each 
! variable.  This is what is done for the final updates of the state 
! variables and velocity.

module make_edge_state_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: make_edge_state, make_edge_state_1d
  
contains

  subroutine make_edge_state(nlevs,s,u,sedge,umac,utrans,force, &
                             normal,w0,w0_cart_vec, &
                             dx,dt,is_vel,the_bc_level,velpred, &
                             start_scomp,start_bccomp,num_comp,mla)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: nr_fine, r_start_coord, r_end_coord, dr, spherical
    use variables, only: foextrap_comp
    use fill_3d_module
    use multifab_physbc_module
    use ml_restriction_module, only : ml_edge_restriction_c

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(in   ) :: s(:),u(:)
    type(multifab) , intent(inout) :: sedge(:,:),umac(:,:)
    type(multifab) , intent(in   ) :: utrans(:,:),force(:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0_cart_vec(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    logical        , intent(in   ) :: is_vel
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    integer        , intent(in   ) :: velpred,start_scomp,start_bccomp,num_comp
    type(ml_layout), intent(in   ) :: mla

    integer                  :: i,r,scomp,bccomp,ng,dm,n
    integer                  :: lo(u(1)%dim), hi(u(1)%dim)
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: uop(:,:,:,:)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: utp(:,:,:,:)
    real(kind=dp_t), pointer :: vtp(:,:,:,:)
    real(kind=dp_t), pointer :: wtp(:,:,:,:)
    real(kind=dp_t), pointer :: w0p(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)
    real(kind=dp_t), pointer :: gw0p(:,:,:,:)

    real(kind=dp_t), allocatable :: gradw0_rad(:)
    type(multifab) :: gradw0_cart

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_edge_state")

    dm = u(1)%dim
    ng = s(1)%ng

    if (spherical .eq. 1) then
       allocate (gradw0_rad(0:nr_fine-1))

       ! NOTE: here we are doing the computation at the finest level
       do r=r_start_coord(nlevs),r_end_coord(nlevs)
          gradw0_rad(r) = (w0(nlevs,r+1) - w0(nlevs,r)) / dr(nlevs)
       enddo
    endif

    do n=1,nlevs

       call multifab_build(gradw0_cart, u(n)%la,1,1)

       if (spherical .eq. 1 .and. is_vel) then

          do i = 1, gradw0_cart%nboxes
             if ( multifab_remote(u(n),i) ) cycle
             gw0p => dataptr(gradw0_cart, i)
             lo = lwb(get_box(gradw0_cart,i))
             hi = upb(get_box(gradw0_cart,i))
             
             call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,gradw0_rad,gw0p, &
                                               lo,hi,dx(n,:),gradw0_cart%ng)
             
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
       

       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n),i) ) cycle
          sop  => dataptr(s(n),i)
          uop  => dataptr(u(n),i)
          sepx => dataptr(sedge(n,1),i)
          sepy => dataptr(sedge(n,2),i)
          ump  => dataptr(umac(n,1),i)
          vmp  => dataptr(umac(n,2),i)
          utp  => dataptr(utrans(n,1),i)
          vtp  => dataptr(utrans(n,2),i)
          fp   => dataptr(force(n),i)
          lo   =  lwb(get_box(s(n),i))
          select case (dm)
          case (2)
             do scomp = start_scomp, start_scomp + num_comp - 1
                bccomp = start_bccomp + scomp - start_scomp
                call make_edge_state_2d(n,sop(:,:,1,:), uop(:,:,1,:), &
                                        sepx(:,:,1,:), sepy(:,:,1,:), &
                                        ump(:,:,1,1), vmp(:,:,1,1), &
                                        utp(:,:,1,1), vtp(:,:,1,1), &
                                        fp(:,:,1,:), w0(n,:), &
                                        lo, dx(n,:), dt, is_vel, &
                                        the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                        the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp:), &
                                        velpred, ng, scomp)
             end do

          case (3)
             wmp  => dataptr(  umac(n,3),i)
             wtp  => dataptr(utrans(n,3),i)
             sepz => dataptr( sedge(n,3),i)
             w0p  => dataptr(w0_cart_vec(n),i)
             gw0p => dataptr(gradw0_cart,i)
             if (spherical .eq. 1) then
                np => dataptr(normal(n),i)
             end if
             do scomp = start_scomp, start_scomp + num_comp - 1
                bccomp = start_bccomp + scomp - start_scomp
                call make_edge_state_3d(n,sop(:,:,:,:), uop(:,:,:,:), &
                                        sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                        ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                        utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), &
                                        fp(:,:,:,:), np(:,:,:,:), &
                                        w0p(:,:,:,:), gw0p(:,:,:,1), &
                                        lo, dx(n,:), dt, is_vel, &
                                        the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                        the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp:), &
                                        velpred, ng, scomp)
             end do
          end select
       end do

       call destroy(gradw0_cart)

    end do
    !
    ! We call ml_edge_restriction for the output velocity if is_vel .eq. .true.
    ! we do not call ml_edge_restriction for scalars because instead we will call
    ! ml_edge_restriction on the fluxes in mkflux.
    !
    if (is_vel .and. velpred .eq. 1) then
       do n = nlevs,2,-1
          do i = 1, dm
             call ml_edge_restriction_c(umac(n-1,i),1,umac(n,i),1,mla%mba%rr(n-1,:),i,1)
          enddo
       enddo
    end if

    if (is_vel .and. velpred .eq. 0) then
       do n = nlevs,2,-1
          do i = 1, dm
             call ml_edge_restriction_c(sedge(n-1,i),1,sedge(n,i),1,mla%mba%rr(n-1,:),i,dm)
          enddo
       enddo
    end if

    if (spherical .eq. 1 .and. is_vel) then
       deallocate(gradw0_rad)
    endif

    call destroy(bpt)
    
  end subroutine make_edge_state

  
  subroutine make_edge_state_2d(n,s,u,sedgex,sedgey,umac,vmac,utrans, &
                                vtrans,force,w0,lo,dx,dt,is_vel,phys_bc,adv_bc,velpred, &
                                ng,comp)

    use geometry, only: nr, r_start_coord, r_end_coord
    use bc_module
    use slope_module
    use bl_constants_module
    use probin_module, only: use_new_godunov

    integer        , intent(in   ) :: n,lo(:)
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(inout) :: sedgex(lo(1)   :,lo(2)   :,:)
    real(kind=dp_t), intent(inout) :: sedgey(lo(1)   :,lo(2)   :,:)
    real(kind=dp_t), intent(inout) ::   umac(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(inout) ::   vmac(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(in   ) :: utrans(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(in   ) :: vtrans(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)- 1:,lo(2)- 1:,:)
    real(kind=dp_t), intent(in   ) ::     w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    logical        , intent(in   ) :: is_vel
    integer        , intent(in   ) :: phys_bc(:,:)
    integer        , intent(in   ) :: adv_bc(:,:,:)
    integer        , intent(in   ) :: velpred
    integer        , intent(in   ) :: ng,comp
    
    ! Local variables
    real(kind=dp_t), allocatable :: slopex(:,:,:),slopey(:,:,:)
    real(kind=dp_t), allocatable :: s_l(:),s_r(:),s_b(:),s_t(:)
    
    real(kind=dp_t) :: ubardth, vbardth
    real(kind=dp_t) :: hx, hy, dth
    real(kind=dp_t) :: splus,sminus
    real(kind=dp_t) :: savg,st
    real(kind=dp_t) :: vlo,vhi
    real(kind=dp_t) :: sptop,spbot,smtop,smbot,splft,sprgt,smlft,smrgt
    real(kind=dp_t) :: abs_eps, eps, umax
    
    integer :: hi(2)
    integer :: i,j,is,js,ie,je
    logical :: test
    
    hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)
    
    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)
    
    allocate(s_l(lo(1)-1:hi(1)+2))
    allocate(s_r(lo(1)-1:hi(1)+2))
    allocate(s_b(lo(2)-1:hi(2)+2))
    allocate(s_t(lo(2)-1:hi(2)+2))
    
    allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1))
    allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1))
    
    call slopex_2d(s(:,:,comp:),slopex,lo,ng,1,adv_bc)
    call slopey_2d(s(:,:,comp:),slopey,lo,ng,1,adv_bc)
    
    abs_eps = 1.0d-8
    
    dth = HALF*dt
    
    hx = dx(1)
    hy = dx(2)
    
    if (velpred .eq. 1) then
       
       umax = abs(utrans(is,js))

       do j = js,je; do i = is,ie+1
          umax = max(umax,abs(utrans(i,j)))
       end do; end do
       do j = js,je+1; do i = is,ie
          umax = max(umax,abs(vtrans(i,j)))
       end do; end do
       
    else 
       
       umax = abs(umac(is,js))

       do j = js,je; do i = is,ie+1
          umax = max(umax,abs(umac(i,j)))
       end do; end do
       do j = js,je+1; do i = is,ie
          umax = max(umax,abs(vmac(i,j)))
       end do; end do

    end if
    
    if (umax .eq. 0.d0) then
       eps = abs_eps
    else
       eps = abs_eps * umax
    endif
    
    !********************************
    ! Loop for edge states on x-edges.
    !********************************

    if (velpred .eq. 0 .or. comp .eq. 1) then
       do j = js,je 
          do i = is-1,ie+1 
             
             vlo = u(i,j  ,2) + HALF * (w0(j  )+w0(j+1))
             if ((j+2) .le. nr(n)) then
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
                if (is_vel .and. comp .eq. 2) then
                   sptop = ZERO
                   spbot = ZERO
                elseif (is_vel .and. comp .eq. 1) then
                   sptop = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                   spbot = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                else
                   sptop = spbot
                endif
             endif
             
             splus = merge(spbot,sptop,vtrans(i,j+1).gt.ZERO)
             savg  = HALF * (spbot + sptop)
             splus = merge(splus, savg, abs(vtrans(i,j+1)) .gt. eps)
             
             if (j.ge.1) then
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
                if (is_vel .and. (comp .eq. 2)) then
                   smtop = ZERO
                   smbot = ZERO
                elseif (is_vel .and. (comp .ne. 2)) then
                   smbot = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                   smtop = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                else
                   smbot = smtop
                endif
             endif
             
             sminus = merge(smbot,smtop,vtrans(i,j).gt.ZERO)
             savg   = HALF * (smbot + smtop)
             sminus = merge(sminus, savg, abs(vtrans(i,j)) .gt. eps)
             
             st = force(i,j,comp) - HALF * (vtrans(i,j)+vtrans(i,j+1))*(splus - sminus) / hy
             
             if (is_vel .and. comp.eq.2) then
                ! vtrans contains w0 so we need to subtract it off
                st = st - HALF * (vtrans(i,j)+vtrans(i,j+1)-w0(j+1)-w0(j))*(w0(j+1)-w0(j))/hy
             end if
             
             ubardth = dth*u(i,j,1)/hx
             
             if(velpred .eq. 1) then
                if (use_new_godunov) then
                   s_l(i+1)= s(i,j,comp) + (HALF-max(ZERO,ubardth))*slopex(i,j,1) + dth*st
                   s_r(i  )= s(i,j,comp) - (HALF+min(ZERO,ubardth))*slopex(i,j,1) + dth*st
                else
                   s_l(i+1)= s(i,j,comp) + (HALF-ubardth)*slopex(i,j,1) + dth*st
                   s_r(i  )= s(i,j,comp) - (HALF+ubardth)*slopex(i,j,1) + dth*st
                end if
             else
                s_l(i+1)= s(i,j,comp) + (HALF-dth*umac(i+1,j)/hx)*slopex(i,j,1) + dth*st
                s_r(i  )= s(i,j,comp) - (HALF+dth*umac(i  ,j)/hx)*slopex(i,j,1) + dth*st
             endif
             
          enddo
          
          if (velpred .eq. 1) then
             do i = is, ie+1 
                savg = HALF*(s_r(i) + s_l(i))
                test = ( (s_l(i) .le. ZERO  .and. s_r(i) .ge. ZERO) .or. &
                     (abs(s_l(i) + s_r(i)) .lt. eps) )
                sedgex(i,j,comp)=merge(s_l(i),s_r(i),savg.gt.ZERO)
                sedgex(i,j,comp)=merge(savg,sedgex(i,j,comp),test)
             enddo
          else
             do i = is, ie+1 
                sedgex(i,j,comp)=merge(s_l(i),s_r(i),umac(i,j).gt.ZERO)
                savg = HALF*(s_r(i) + s_l(i))
                sedgex(i,j,comp)=merge(savg,sedgex(i,j,comp),abs(umac(i,j)) .lt. eps)
             enddo
          endif
          
          if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
             if (is_vel .and. comp .eq. 1) then
                sedgex(is,j,comp) = ZERO
             elseif (is_vel .and. comp .ne. 1) then
                sedgex(is,j,comp) = merge(ZERO,s_r(is),phys_bc(1,1).eq.NO_SLIP_WALL)
             else 
                sedgex(is,j,comp) = s_r(is)
             endif
          elseif (phys_bc(1,1) .eq. INLET) then
             sedgex(is,j,comp) = s(is-1,j,comp)
          elseif (phys_bc(1,1) .eq. OUTLET) then
             if (is_vel .and. comp.eq.1) then
                sedgex(is,j,comp) = MIN(s_r(is),ZERO)
             else
                sedgex(is,j,comp) = s_r(is)
             end if
          endif
          if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
             if (is_vel .and. comp .eq. 1) then
                sedgex(ie+1,j,comp) = ZERO
             else if (is_vel .and. comp .ne. 1) then
                sedgex(ie+1,j,comp) = merge(ZERO,s_l(ie+1),phys_bc(1,2).eq.NO_SLIP_WALL)
             else 
                sedgex(ie+1,j,comp) = s_l(ie+1)
             endif
          elseif (phys_bc(1,2) .eq. INLET) then
             sedgex(ie+1,j,comp) = s(ie+1,j,comp)
          elseif (phys_bc(1,2) .eq. OUTLET) then
             if (is_vel .and. comp.eq.1) then
                sedgex(ie+1,j,comp) = MAX(s_l(ie+1),ZERO)
             else
                sedgex(ie+1,j,comp) = s_l(ie+1)
             end if
          endif
          
          if (velpred .eq. 1) then
             do i = is, ie+1 
                umac(i,j) = sedgex(i,j,1)
             enddo
          endif
       enddo
    endif
    
    !********************************
    ! Loop for edge states on y-edges.
    !********************************

    if (velpred .eq. 0 .or. comp .eq. 2) then
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
                if (is_vel .and. comp .eq. 1) then
                   splft = ZERO
                   sprgt = ZERO
                elseif (is_vel .and. comp .ne. 1) then
                   sprgt = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                   splft = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                else
                   sprgt = splft
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
                if (is_vel .and. comp .eq. 1) then
                   smlft = ZERO
                   smrgt = ZERO
                elseif (is_vel .and. comp .ne. 1) then
                   smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                   smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                else
                   smlft = smrgt
                endif
             endif
             
             sminus = merge(smlft,smrgt,utrans(i,j).gt.ZERO)
             savg   = HALF * (smlft + smrgt)
             sminus = merge(sminus, savg, abs(utrans(i,j)) .gt. eps)
             
             st = force(i,j,comp) - HALF * (utrans(i,j)+utrans(i+1,j))*(splus - sminus) / hx
             
             if (is_vel .and. comp .eq. 2) then
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
             
             if(velpred .eq. 1) then
                if (use_new_godunov) then
                   s_b(j+1)= s(i,j,comp) + (HALF-max(ZERO,vbardth))*slopey(i,j,1) + dth*st
                   s_t(j  )= s(i,j,comp) - (HALF+min(ZERO,vbardth))*slopey(i,j,1) + dth*st
                else
                   s_b(j+1)= s(i,j,comp) + (HALF-vbardth)*slopey(i,j,1) + dth*st
                   s_t(j  )= s(i,j,comp) - (HALF+vbardth)*slopey(i,j,1) + dth*st
                end if
             else
                s_b(j+1)= s(i,j,comp) + (HALF-dth*vmac(i,j+1)/hy)*slopey(i,j,1) + dth*st
                s_t(j  )= s(i,j,comp) - (HALF+dth*vmac(i,j  )/hy)*slopey(i,j,1) + dth*st
             endif
             
          enddo
          
          if (velpred .eq. 1) then
             do j = js, je+1 
                savg = HALF*(s_b(j) + s_t(j))
                test = ( (s_b(j) .le. ZERO .and. s_t(j) .ge. ZERO) .or. &
                     (abs(s_b(j) + s_t(j)) .lt. eps) )
                sedgey(i,j,comp)=merge(s_b(j),s_t(j),savg.gt.ZERO)
                sedgey(i,j,comp)=merge(savg,sedgey(i,j,comp),test)
             enddo
             
          else
             
             do j = js, je+1 
                sedgey(i,j,comp)=merge(s_b(j),s_t(j),vmac(i,j).gt.ZERO)
                savg = HALF*(s_b(j) + s_t(j))
                sedgey(i,j,comp)=merge(savg,sedgey(i,j,comp),abs(vmac(i,j)) .lt. eps)
             enddo

          endif
          
          if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
             if (is_vel .and. comp .eq. 2) then
                sedgey(i,js,comp) = ZERO
             elseif (is_vel .and. comp .ne. 2) then
                sedgey(i,js,comp) = merge(ZERO,s_t(js),phys_bc(2,1).eq.NO_SLIP_WALL)
             else 
                sedgey(i,js,comp) = s_t(js)
             endif
          elseif (phys_bc(2,1) .eq. INLET) then
             sedgey(i,js,comp) = s(i,js-1,comp)
          elseif (phys_bc(2,1) .eq. OUTLET) then
             if (is_vel .and. comp.eq.2) then
                sedgey(i,js,comp) = MIN(s_t(js),ZERO)
             else
                sedgey(i,js,comp) = s_t(js)
             end if
          endif
          
          if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
             if (is_vel .and. comp .eq. 2) then
                sedgey(i,je+1,comp) = ZERO
             elseif (is_vel .and. comp .ne. 2) then
                sedgey(i,je+1,comp) = merge(ZERO,s_b(je+1),phys_bc(2,2).eq.NO_SLIP_WALL)
             else 
                sedgey(i,je+1,comp) = s_b(je+1)
             endif
          elseif (phys_bc(2,2) .eq. INLET) then
             sedgey(i,je+1,comp) = s(i,je+1,comp)
          elseif (phys_bc(2,2) .eq. OUTLET) then
             if (is_vel .and. comp.eq.2) then
                sedgey(i,je+1,comp) = MAX(s_b(je+1),ZERO)
             else
                sedgey(i,je+1,comp) = s_b(je+1)
             end if
          endif
          
          if (velpred .eq. 1) then
             do j = js, je+1 
                vmac(i,j) = sedgey(i,j,2)
             enddo
          end if
          
       enddo
    end if
    
  end subroutine make_edge_state_2d
  
  
  subroutine make_edge_state_3d(n, s,u, &
                                sedgex,sedgey,sedgez, &
                                umac,vmac,wmac, &
                                utrans,vtrans,wtrans, &
                                force, normal, &
                                w0_cart_vec,gradw0_cart, &
                                lo,dx,dt,is_vel, &
                                phys_bc,adv_bc, &
                                velpred,ng,comp)

    use bc_module
    use slope_module
    use bl_constants_module
    use geometry, only: spherical, r_start_coord, r_end_coord, nr

    integer        , intent(in   ) :: n, lo(:)
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(inout) :: sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t), intent(inout) :: sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t), intent(inout) :: sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t), intent(inout) ::   umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real(kind=dp_t), intent(inout) ::   vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real(kind=dp_t), intent(inout) ::   wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real(kind=dp_t), intent(in   ) :: utrans(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real(kind=dp_t), intent(in   ) :: vtrans(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real(kind=dp_t), intent(in   ) :: wtrans(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real(kind=dp_t), intent(in   ) :: normal(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real(kind=dp_t), intent(in   ) :: w0_cart_vec(lo(1)- 2:,lo(2)- 2:,lo(3)- 2:,:)
    real(kind=dp_t), intent(in   ) :: gradw0_cart(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    logical        , intent(in   ) :: is_vel
    integer        , intent(in   ) :: phys_bc(:,:)
    integer        , intent(in   ) :: adv_bc(:,:,:)
    integer        , intent(in   ) :: velpred
    integer        , intent(in   ) :: ng,comp

    real(kind=dp_t), allocatable :: slopex(:,:,:,:),slopey(:,:,:,:),slopez(:,:,:,:)
    real(kind=dp_t), allocatable :: s_l(:),s_r(:),s_b(:),s_t(:),s_u(:),s_d(:)
    
    real(kind=dp_t) :: ubardth, vbardth, wbardth
    real(kind=dp_t) :: hx, hy, hz, dth, splus, sminus
    real(kind=dp_t) :: savg,st,ulo,uhi,vlo,vhi,wlo,whi
    real(kind=dp_t) :: sptop,spbot,smtop,smbot,splft,sprgt,smlft,smrgt
    real(kind=dp_t) :: abs_eps, eps, umax

    real(kind=dp_t) :: Ut_dot_er

    integer :: hi(3),i,j,k,is,js,ie,je,ks,ke
    logical :: test
    
    hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
    hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)
    hi(3) = lo(3) + size(s,dim=3) - (2*ng+1)
    
    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)
    ks = lo(3)
    ke = hi(3)
    
    allocate(s_l(lo(1)-1:hi(1)+2))
    allocate(s_r(lo(1)-1:hi(1)+2))
    allocate(s_b(lo(2)-1:hi(2)+2))
    allocate(s_t(lo(2)-1:hi(2)+2))
    allocate(s_d(lo(3)-1:hi(3)+2))
    allocate(s_u(lo(3)-1:hi(3)+2))

    allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))
    allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))
    allocate(slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))
    
    do k = lo(3)-1,hi(3)+1
       call slopex_2d(s(:,:,k,comp:),slopex(:,:,k,:),lo,ng,1,adv_bc)
       call slopey_2d(s(:,:,k,comp:),slopey(:,:,k,:),lo,ng,1,adv_bc)
    end do
    call slopez_3d(s(:,:,:,comp:),slopez,lo,ng,1,adv_bc)
    
    abs_eps = 1.0d-8
    
    dth = HALF*dt
    
    hx = dx(1)
    hy = dx(2)
    hz = dx(3)
    
    if (velpred .eq. 1) then
       
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
        
     else 
        
        umax = abs(umac(is,js,ks))

        do k = ks,ke; do j = js,je; do i = is,ie+1
           umax = max(umax,abs(umac(i,j,k)))
        end do; end do; end do
        do k = ks,ke; do j = js,je+1; do i = is,ie
           umax = max(umax,abs(vmac(i,j,k)))
        end do; end do; end do
        do k = ks,ke+1; do j = js,je; do i = is,ie
           umax = max(umax,abs(wmac(i,j,k)))
        end do; end do; end do

     end if
     
     if (umax .eq. 0.d0) then
        eps = abs_eps
     else
        eps = abs_eps * umax
     endif
     
     !********************************
     ! Loop for edge states on x-edges.
     !********************************

     if (velpred .eq. 0 .or. comp .eq. 1) then
        do k = ks,ke 
           do j = js,je 
              do i = is-1,ie+1 
                 
                 ! Do transverse in j direction
                 vlo = u(i,j  ,k,2) + w0_cart_vec(i,j  ,k,2)
                 vhi = u(i,j+1,k,2) + w0_cart_vec(i,j+1,k,2)
                 
                 spbot = s(i,j  ,k,comp) + (HALF - dth*vlo/hy) * slopey(i,j  ,k,1)
                 sptop = s(i,j+1,k,comp) - (HALF + dth*vhi/hy) * slopey(i,j+1,k,1)
                 
                 sptop = merge(s(i,je+1,k,comp),sptop,j.eq.je .and. phys_bc(2,2) .eq. INLET)
                 spbot = merge(s(i,je+1,k,comp),spbot,j.eq.je .and. phys_bc(2,2) .eq. INLET)
                 
                 if (j .eq. je .and. &
                      (phys_bc(2,2).eq.SLIP_WALL.or.phys_bc(2,2).eq.NO_SLIP_WALL)) then
                    if (is_vel .and. comp.eq.2) then
                       sptop = ZERO
                       spbot = ZERO
                    elseif (is_vel .and. comp.ne.2) then
                       sptop = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                       spbot = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                    else
                       sptop = spbot
                    endif
                 endif
                 
                 splus = merge(spbot,sptop,vtrans(i,j+1,k).gt.ZERO)
                 savg  = HALF * (spbot + sptop)
                 splus = merge(splus, savg, abs(vtrans(i,j+1,k)) .gt. eps)
                 
                 vlo = u(i,j-1,k,2) + w0_cart_vec(i,j-1,k,2)
                 vhi = u(i,j  ,k,2) + w0_cart_vec(i,j  ,k,2)
                 
                 smbot = s(i,j-1,k,comp) + (HALF - dth*vlo/hy) * slopey(i,j-1,k,1)
                 smtop = s(i,j  ,k,comp) - (HALF + dth*vhi/hy) * slopey(i,j  ,k,1)
                 
                 smtop = merge(s(i,js-1,k,comp),smtop,j.eq.js .and. phys_bc(2,1) .eq. INLET)
                 smbot = merge(s(i,js-1,k,comp),smbot,j.eq.js .and. phys_bc(2,1) .eq. INLET)
                 
                 if (j .eq. js .and. &
                      (phys_bc(2,1).eq.SLIP_WALL.or.phys_bc(2,1).eq.NO_SLIP_WALL)) then
                    if (is_vel .and. (comp .eq. 2)) then
                       smtop = ZERO
                       smbot = ZERO
                    elseif (is_vel .and. (comp .ne. 2)) then
                       smbot = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                       smtop = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                    else
                       smbot = smtop
                    endif
                 endif
                 
                 sminus = merge(smbot,smtop,vtrans(i,j,k).gt.ZERO)
                 savg   = HALF * (smbot + smtop)
                 sminus = merge(sminus, savg, abs(vtrans(i,j,k)) .gt. eps)
                 
                 st = force(i,j,k,comp) - &
                      HALF * (vtrans(i,j,k)+vtrans(i,j+1,k))*(splus - sminus) / hy
                 
                 ! Do transverse in k direction
                 wlo = u(i,j,k  ,3) + w0_cart_vec(i,j,k  ,3)
                 whi = u(i,j,k+1,3) + w0_cart_vec(i,j,k+1,3)
                 
                 spbot = s(i,j,k  ,comp) + (HALF - dth*wlo/hz) * slopez(i,j,k  ,1)
                 sptop = s(i,j,k+1,comp) - (HALF + dth*whi/hz) * slopez(i,j,k+1,1)
                 
                 sptop = merge(s(i,j,ke+1,comp),sptop,k.eq.ke .and. phys_bc(3,2) .eq. INLET)
                 spbot = merge(s(i,j,ke+1,comp),spbot,k.eq.ke .and. phys_bc(3,2) .eq. INLET)
                 
                 if (k .eq. ke .and. &
                      (phys_bc(3,2).eq.SLIP_WALL.or.phys_bc(3,2).eq.NO_SLIP_WALL)) then
                    if (is_vel .and. comp .eq. 3) then
                       sptop = ZERO
                       spbot = ZERO
                    elseif (is_vel .and. (comp.ne.3)) then
                       sptop = merge(ZERO,spbot,phys_bc(3,2).eq.NO_SLIP_WALL)
                       spbot = merge(ZERO,spbot,phys_bc(3,2).eq.NO_SLIP_WALL)
                    else
                       sptop = spbot
                    endif
                 endif
                 
                 splus = merge(spbot,sptop,wtrans(i,j,k+1).gt.ZERO)
                 savg  = HALF * (spbot + sptop)
                 splus = merge(splus, savg, abs(wtrans(i,j,k+1)) .gt. eps)
                 
                 wlo = u(i,j,k-1,3) + w0_cart_vec(i,j,k-1,3)
                 whi = u(i,j,k  ,3) + w0_cart_vec(i,j,k  ,3)
                 
                 smtop = s(i,j,k  ,comp) - (HALF + dth*whi/hz) * slopez(i,j,k  ,1)
                 smbot = s(i,j,k-1,comp) + (HALF - dth*wlo/hz) * slopez(i,j,k-1,1)
                 
                 smtop = merge(s(i,j,ks-1,comp),smtop,k.eq.ks .and. phys_bc(3,1) .eq. INLET)
                 smbot = merge(s(i,j,ks-1,comp),smbot,k.eq.ks .and. phys_bc(3,1) .eq. INLET)
                 
                 if (k .eq. ks .and. &
                      (phys_bc(3,1).eq.SLIP_WALL.or.phys_bc(3,1).eq.NO_SLIP_WALL)) then
                    if (is_vel .and. (comp.eq.3)) then
                       smtop = ZERO
                       smbot = ZERO
                    elseif (is_vel .and. (comp.ne.3)) then
                       smbot = merge(ZERO,smtop,phys_bc(3,1).eq.NO_SLIP_WALL)
                       smtop = merge(ZERO,smtop,phys_bc(3,1).eq.NO_SLIP_WALL)
                    else
                       smbot = smtop
                    endif
                 endif
                 
                 sminus = merge(smbot,smtop,wtrans(i,j,k).gt.ZERO)
                 savg   = HALF * (smbot + smtop)
                 sminus = merge(sminus, savg, abs(wtrans(i,j,k)) .gt. eps)
                 
                 st = st - HALF * (wtrans(i,j,k)+wtrans(i,j,k+1))*(splus - sminus) / hz
                 
                 if (is_vel) then
                    ! add the (Utilde . e_r) d w_0 /dr e_r term here

                    if (spherical .eq. 0 .and. comp.eq.3) then

                       ! wtrans contains w0 so we need to subtract it off
                       st = st - HALF * (wtrans(i,j,k)+wtrans(i,j,k+1)- &
                            w0_cart_vec(i,j,k+1,3)-w0_cart_vec(i,j,k,3))* &
                            (w0_cart_vec(i,j,k+1,3)-w0_cart_vec(i,j,k,3))/hz

                    else if (spherical .eq. 1) then

                       ! u/v/wtrans contain w0, so we need to subtract it off.  
                       ! Note w0_cart_vec is cell-centered
                       Ut_dot_er = (HALF*(utrans(i,j,k) + utrans(i+1,j,k)) - &
                                          w0_cart_vec(i,j,k,1))*normal(i,j,k,1) + &
                                   (HALF*(vtrans(i,j,k) + vtrans(i,j+1,k)) - &
                                          w0_cart_vec(i,j,k,2))*normal(i,j,k,2) + &
                                   (HALF*(wtrans(i,j,k) + wtrans(i,j,k+1)) - &
                                          w0_cart_vec(i,j,k,3))*normal(i,j,k,3)
                       
                       if (comp .eq. 1) then
                          st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,1)
                       else if (comp .eq. 2) then
                          st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,2)
                       else if (comp .eq. 3) then
                          st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,3)
                       endif
                                        
                    endif

                 end if
                 
                 ubardth = dth/hx * ( u(i,j,k,1) + w0_cart_vec(i,j,k,1))
                 
                 if (velpred .eq. 1) then
                    s_l(i+1)= s(i,j,k,comp) + (HALF-ubardth)*slopex(i,j,k,1) + dth*st
                    s_r(i  )= s(i,j,k,comp) - (HALF+ubardth)*slopex(i,j,k,1) + dth*st
                 else
                    s_l(i+1)= s(i,j,k,comp) + (HALF-dth*umac(i+1,j,k)/hx)*slopex(i,j,k,1) &
                         + dth*st
                    s_r(i  )= s(i,j,k,comp) - (HALF+dth*umac(i  ,j,k)/hx)*slopex(i,j,k,1) &
                         + dth*st
                 endif
                 
              enddo
              
              if (velpred .eq. 1) then
                 do i = is, ie+1 
                    savg = HALF*(s_r(i) + s_l(i))
                    test = ( (s_l(i) .le. ZERO .and. s_r(i) .ge. ZERO) .or. &
                         (abs(s_l(i) + s_r(i)) .lt. eps) )
                    sedgex(i,j,k,comp)=merge(s_l(i),s_r(i),savg.gt.ZERO)
                    sedgex(i,j,k,comp)=merge(savg,sedgex(i,j,k,comp),test)
                 enddo
              else
                 do i = is, ie+1 
                    sedgex(i,j,k,comp)=merge(s_l(i),s_r(i),umac(i,j,k).gt.ZERO)
                    savg = HALF*(s_r(i) + s_l(i))
                    sedgex(i,j,k,comp)=merge(savg,sedgex(i,j,k,comp),abs(umac(i,j,k)) .lt. eps)
                 enddo
              endif
              
              if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                 if (is_vel .and. comp .eq. 1) then
                    sedgex(is,j,k,comp) = ZERO
                 elseif (is_vel .and. comp .ne. 1) then
                    sedgex(is,j,k,comp) = merge(ZERO,s_r(is),phys_bc(1,1).eq.NO_SLIP_WALL)
                 else 
                    sedgex(is,j,k,comp) = s_r(is)
                 endif
              elseif (phys_bc(1,1) .eq. INLET) then
                 sedgex(is,j,k,comp) = s(is-1,j,k,comp)
              elseif (phys_bc(1,1) .eq. OUTLET) then
                 if (is_vel .and. comp.eq.1) then
                    sedgex(is,j,k,comp) = MIN(s_r(is),ZERO)
                 else
                    sedgex(is,j,k,comp) = s_r(is)
                 end if
              endif
              if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                 if (is_vel .and. comp .eq. 1) then
                    sedgex(ie+1,j,k,comp) = ZERO
                 else if (is_vel .and. comp .ne. 1) then
                    sedgex(ie+1,j,k,comp) = merge(ZERO,s_l(ie+1),phys_bc(1,2).eq.NO_SLIP_WALL)
                 else 
                    sedgex(ie+1,j,k,comp) = s_l(ie+1)
                 endif
              elseif (phys_bc(1,2) .eq. INLET) then
                 sedgex(ie+1,j,k,comp) = s(ie+1,j,k,comp)
              elseif (phys_bc(1,2) .eq. OUTLET) then
                 if (is_vel .and. comp.eq.1) then
                    sedgex(ie+1,j,k,comp) = MAX(s_l(ie+1),ZERO)
                 else
                    sedgex(ie+1,j,k,comp) = s_l(ie+1)
                 end if
              endif
              
              if (velpred .eq. 1) then
                 do i = is, ie+1 
                    umac(i,j,k) = sedgex(i,j,k,1)
                 enddo
              endif
              
           enddo
        enddo
     endif

     !********************************
     ! Loop for edge states on y-edges.
     !********************************

     if (velpred .eq. 0 .or. comp .eq. 2) then
        do k = ks, ke 
           do i = is, ie 
              do j = js-1, je+1 
                 
                 ! Do transverse in i direction
                 ulo = u(i  ,j,k,1) + w0_cart_vec(i  ,j,k,1)
                 uhi = u(i+1,j,k,1) + w0_cart_vec(i+1,j,k,1)
                 
                 splft = s(i  ,j,k,comp) + (HALF - dth*ulo/hx) * slopex(i  ,j,k,1)
                 sprgt = s(i+1,j,k,comp) - (HALF + dth*uhi/hx) * slopex(i+1,j,k,1)
                 
                 sprgt = merge(s(ie+1,j,k,comp),sprgt,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
                 splft = merge(s(ie+1,j,k,comp),splft,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
                 
                 if (i .eq. ie .and. &
                      (phys_bc(1,2).eq.SLIP_WALL.or.phys_bc(1,2).eq.NO_SLIP_WALL)) then
                    if (is_vel .and. comp .eq. 1) then
                       splft = ZERO
                       sprgt = ZERO
                    elseif (is_vel .and. comp .ne. 1) then
                       sprgt = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                       splft = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                    else
                       sprgt = splft
                    endif
                 endif
                 
                 splus = merge(splft,sprgt,utrans(i+1,j,k).gt.ZERO)
                 savg  = HALF * (splft + sprgt)
                 splus = merge(splus, savg, abs(utrans(i+1,j,k)) .gt. eps)
                 
                 ulo = u(i-1,j,k,1) + w0_cart_vec(i-1,j,k,1)
                 uhi = u(i  ,j,k,1) + w0_cart_vec(i  ,j,k,1)
                 
                 smlft = s(i-1,j,k,comp) + (HALF - dth*ulo/hx) * slopex(i-1,j,k,1)
                 smrgt = s(i  ,j,k,comp) - (HALF + dth*uhi/hx) * slopex(i  ,j,k,1)
                 
                 smrgt = merge(s(is-1,j,k,comp),smrgt,i.eq.is .and. phys_bc(1,1) .eq. INLET)
                 smlft = merge(s(is-1,j,k,comp),smlft,i.eq.is .and. phys_bc(1,1) .eq. INLET)
                 
                 if (i .eq. is .and. &
                      (phys_bc(1,1).eq.SLIP_WALL.or.phys_bc(1,1).eq.NO_SLIP_WALL)) then
                    if (is_vel .and. comp .eq. 1) then
                       smlft = ZERO
                       smrgt = ZERO
                    elseif (is_vel .and. comp .ne. 1) then
                       smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                       smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                    else
                       smlft = smrgt
                    endif
                 endif
                 
                 sminus = merge(smlft,smrgt,utrans(i,j,k).gt.ZERO)
                 savg   = HALF * (smlft + smrgt)
                 sminus = merge(sminus, savg, abs(utrans(i,j,k)) .gt. eps)
                 
                 st = force(i,j,k,comp) - &
                      HALF * (utrans(i,j,k)+utrans(i+1,j,k))*(splus - sminus) / hx
                 
                 ! Do transverse in k direction
                 wlo = u(i,j,k  ,3) + w0_cart_vec(i,j,k  ,3)
                 whi = u(i,j,k+1,3) + w0_cart_vec(i,j,k+1,3)
                 
                 splft = s(i,j,k  ,comp) + (HALF - dth*wlo/hz) * slopez(i,j,k  ,1)
                 sprgt = s(i,j,k+1,comp) - (HALF + dth*whi/hz) * slopez(i,j,k+1,1)
                 
                 sprgt = merge(s(i,j,ke+1,comp),sprgt,k.eq.ke .and. phys_bc(3,2) .eq. INLET)
                 splft = merge(s(i,j,ke+1,comp),splft,k.eq.ke .and. phys_bc(3,2) .eq. INLET)
                 
                 if (k .eq. ke .and. &
                      (phys_bc(3,2).eq.SLIP_WALL.or.phys_bc(3,2).eq.NO_SLIP_WALL)) then
                    if (is_vel .and. comp .eq. 3) then
                       splft = ZERO
                       sprgt = ZERO
                    elseif (is_vel .and. comp .ne. 3) then
                       sprgt = merge(ZERO,splft,phys_bc(3,2).eq.NO_SLIP_WALL)
                       splft = merge(ZERO,splft,phys_bc(3,2).eq.NO_SLIP_WALL)
                    else
                       sprgt = splft
                    endif
                 endif
                 
                 splus = merge(splft,sprgt,wtrans(i,j,k+1).gt.ZERO)
                 savg  = HALF * (splft + sprgt)
                 splus = merge(splus, savg, abs(wtrans(i,j,k+1)) .gt. eps)
                 
                 wlo = u(i,j,k-1,3) + w0_cart_vec(i,j,k-1,3)
                 whi = u(i,j,k  ,3) + w0_cart_vec(i,j,k  ,3)
                 
                 smrgt = s(i,j,k  ,comp) - (HALF + dth*whi/hz) * slopez(i,j,k  ,1)
                 smlft = s(i,j,k-1,comp) + (HALF - dth*wlo/hz) * slopez(i,j,k-1,1)
                 
                 smrgt = merge(s(i,j,ks-1,comp),smrgt,k.eq.ks .and. phys_bc(3,1) .eq. INLET)
                 smlft = merge(s(i,j,ks-1,comp),smlft,k.eq.ks .and. phys_bc(3,1) .eq. INLET)
                 
                 if (k .eq. ks .and. &
                      (phys_bc(3,1).eq.SLIP_WALL.or.phys_bc(3,1).eq.NO_SLIP_WALL)) then
                    if (is_vel .and. comp .eq. 3) then
                       smlft = ZERO
                       smrgt = ZERO
                    elseif (is_vel .and. comp .ne. 3) then
                       smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                       smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                    else
                       smlft = smrgt
                    endif
                 endif
                 
                 sminus = merge(smlft,smrgt,wtrans(i,j,k).gt.ZERO)
                 savg   = HALF * (smlft + smrgt)
                 sminus = merge(sminus, savg, abs(wtrans(i,j,k)) .gt. eps)
                 
                 st = st - HALF * (wtrans(i,j,k)+wtrans(i,j,k+1))*(splus - sminus) / hz
                 

                 if (is_vel) then
                    ! add the (Utilde . e_r) d w_0 /dr e_r term here

                    if (spherical .eq. 0 .and. comp.eq.3) then

                       ! wtrans contains w0 so we need to subtract it off
                       st = st - HALF * (wtrans(i,j,k)+wtrans(i,j,k+1)- &
                            w0_cart_vec(i,j,k+1,3)-w0_cart_vec(i,j,k,3))* &
                            (w0_cart_vec(i,j,k+1,3)-w0_cart_vec(i,j,k,3))/hz

                    else if (spherical .eq. 1) then

                       ! u/v/wtrans contain w0, so we need to subtract it off.  
                       ! Note w0_cart_vec is cell-centered
                       Ut_dot_er = (HALF*(utrans(i,j,k) + utrans(i+1,j,k)) - &
                                          w0_cart_vec(i,j,k,1))*normal(i,j,k,1) + &
                                   (HALF*(vtrans(i,j,k) + vtrans(i,j+1,k)) - &
                                          w0_cart_vec(i,j,k,2))*normal(i,j,k,2) + &
                                   (HALF*(wtrans(i,j,k) + wtrans(i,j,k+1)) - &
                                          w0_cart_vec(i,j,k,3))*normal(i,j,k,3)

                       if (comp .eq. 1) then
                          st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,1)
                       else if (comp .eq. 2) then
                          st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,2)
                       else if (comp .eq. 3) then
                          st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,3)
                       endif
                                        
                    endif

                 end if

                 vbardth = dth/hy * ( u(i,j,k,2) + w0_cart_vec(i,j,k,2))
                 
                 if(velpred .eq. 1) then
                    s_b(j+1)= s(i,j,k,comp) + (HALF-vbardth)*slopey(i,j,k,1) + dth*st
                    s_t(j  )= s(i,j,k,comp) - (HALF+vbardth)*slopey(i,j,k,1) + dth*st
                 else
                    s_b(j+1)= s(i,j,k,comp) + (HALF-dth*vmac(i,j+1,k)/hy)*slopey(i,j,k,1) &
                         + dth*st
                    s_t(j  )= s(i,j,k,comp) - (HALF+dth*vmac(i,j,  k)/hy)*slopey(i,j,k,1) &
                         + dth*st
                 endif
                 
              enddo
              
              if (velpred .eq. 1) then
                 do j = js, je+1 
                    savg = HALF*(s_b(j) + s_t(j))
                    test = ( (s_b(j) .le. ZERO .and. s_t(j) .ge. ZERO) .or. &
                         (abs(s_b(j) + s_t(j)) .lt. eps) )
                    sedgey(i,j,k,comp)=merge(s_b(j),s_t(j),savg.gt.ZERO)
                    sedgey(i,j,k,comp)=merge(savg,sedgey(i,j,k,comp),test)
                 enddo
              else
                 do j = js, je+1 
                    sedgey(i,j,k,comp)=merge(s_b(j),s_t(j),vmac(i,j,k).gt.ZERO)
                    savg = HALF*(s_b(j) + s_t(j))
                    sedgey(i,j,k,comp)=merge(savg,sedgey(i,j,k,comp),abs(vmac(i,j,k)) .lt. eps)
                 enddo
              endif
              
              if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                 if (is_vel .and. comp .eq. 2) then
                    sedgey(i,js,k,comp) = ZERO
                 elseif (is_vel .and. comp .ne. 2) then
                    sedgey(i,js,k,comp) = merge(ZERO,s_t(js),phys_bc(2,1).eq.NO_SLIP_WALL)
                 else 
                    sedgey(i,js,k,comp) = s_t(js)
                 endif
              elseif (phys_bc(2,1) .eq. INLET) then
                 sedgey(i,js,k,comp) = s(i,js-1,k,comp)
              elseif (phys_bc(2,1) .eq. OUTLET) then
                 if (is_vel .and. comp.eq.2) then
                    sedgey(i,js,k,comp) = MIN(s_t(js),ZERO)
                 else
                    sedgey(i,js,k,comp) = s_t(js)
                 end if
              endif
              
              if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                 if (is_vel .and. comp .eq. 2) then
                    sedgey(i,je+1,k,comp) = ZERO
                 elseif (is_vel .and. comp .ne. 2) then
                    sedgey(i,je+1,k,comp) = merge(ZERO,s_b(je+1),phys_bc(2,2).eq.NO_SLIP_WALL)
                 else 
                    sedgey(i,je+1,k,comp) = s_b(je+1)
                 endif
              elseif (phys_bc(2,2) .eq. INLET) then
                 sedgey(i,je+1,k,comp) = s(i,je+1,k,comp)
              elseif (phys_bc(2,2) .eq. OUTLET) then
                 if (is_vel .and. comp.eq.2) then
                    sedgey(i,je+1,k,comp) = MAX(s_b(je+1),ZERO)
                 else
                    sedgey(i,je+1,k,comp) = s_b(je+1)
                 end if
              endif
              
              if (velpred .eq. 1) then
                 do j = js, je+1 
                    vmac(i,j,k) = sedgey(i,j,k,2)
                 enddo
              endif
              
           enddo
        enddo
     endif

     !********************************
     ! Loop for edge states on z-edges.
     !********************************

     if (velpred .eq. 0 .or. comp .eq. 3) then
        do j = js, je 
           do i = is, ie 
              do k = ks-1,ke+1
                 
                 ! Do transverse in i direction
                 ulo = u(i  ,j,k,1) + w0_cart_vec(i  ,j,k,1)
                 uhi = u(i+1,j,k,1) + w0_cart_vec(i+1,j,k,1)
                 
                 splft = s(i  ,j,k,comp) + (HALF - dth*ulo/hx) * slopex(i  ,j,k,1)
                 sprgt = s(i+1,j,k,comp) - (HALF + dth*uhi/hx) * slopex(i+1,j,k,1)
                 
                 sprgt = merge(s(ie+1,j,k,comp),sprgt,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
                 splft = merge(s(ie+1,j,k,comp),splft,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
                 
                 if (i .eq. ie .and. &
                      (phys_bc(1,2).eq.SLIP_WALL.or.phys_bc(1,2).eq.NO_SLIP_WALL)) then
                    if (is_vel .and. comp .eq. 1) then
                       splft = ZERO
                       sprgt = ZERO
                    elseif (is_vel .and. comp .ne. 1) then
                       sprgt = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                       splft = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
                    else
                       sprgt = splft
                    endif
                 endif
                 
                 splus = merge(splft,sprgt,utrans(i+1,j,k).gt.ZERO)
                 savg  = HALF * (splft + sprgt)
                 splus = merge(splus, savg, abs(utrans(i+1,j,k)) .gt. eps)
                 
                 ulo = u(i-1,j,k,1) + w0_cart_vec(i-1,j,k,1)
                 uhi = u(i  ,j,k,1) + w0_cart_vec(i  ,j,k,1)
                 
                 smlft = s(i-1,j,k,comp) + (HALF - dth*ulo/hx) * slopex(i-1,j,k,1)
                 smrgt = s(i  ,j,k,comp) - (HALF + dth*uhi/hx) * slopex(i  ,j,k,1)
                 
                 smrgt = merge(s(is-1,j,k,comp),smrgt,i.eq.is .and. phys_bc(1,1) .eq. INLET)
                 smlft = merge(s(is-1,j,k,comp),smlft,i.eq.is .and. phys_bc(1,1) .eq. INLET)
                 
                 if (i .eq. is .and. &
                      (phys_bc(1,1).eq.SLIP_WALL.or.phys_bc(1,1).eq.NO_SLIP_WALL)) then
                    if (is_vel .and. comp .eq. 1) then
                       smlft = ZERO
                       smrgt = ZERO
                    elseif (is_vel .and. comp .ne. 1) then
                       smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                       smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
                    else
                       smlft = smrgt
                    endif
                 endif
                 
                 sminus = merge(smlft,smrgt,utrans(i,j,k).gt.ZERO)
                 savg   = HALF * (smlft + smrgt)
                 sminus = merge(sminus, savg, abs(utrans(i,j,k)) .gt. eps)
                 
                 st = force(i,j,k,comp) - &
                      HALF * (utrans(i,j,k)+utrans(i+1,j,k))*(splus - sminus) / hx
                 
                 ! Do transverse in j direction
                 vlo = u(i,j  ,k,2) + w0_cart_vec(i,j  ,k,2)
                 vhi = u(i,j+1,k,2) + w0_cart_vec(i,j+1,k,2)
                 
                 spbot = s(i,j  ,k,comp) + (HALF - dth*vlo/hy) * slopey(i,j  ,k,1)
                 sptop = s(i,j+1,k,comp) - (HALF + dth*vhi/hy) * slopey(i,j+1,k,1)
                 
                 sptop = merge(s(i,je+1,k,comp),sptop,j.eq.je .and. phys_bc(2,2) .eq. INLET)
                 spbot = merge(s(i,je+1,k,comp),spbot,j.eq.je .and. phys_bc(2,2) .eq. INLET)
                 
                 if (j .eq. je .and. &
                      (phys_bc(2,2).eq.SLIP_WALL.or.phys_bc(2,2).eq.NO_SLIP_WALL)) then
                    if (is_vel .and. comp.eq.2) then
                       sptop = ZERO
                       spbot = ZERO
                    elseif (is_vel .and. comp.ne.2) then
                       sptop = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                       spbot = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
                    else
                       sptop = spbot
                    endif
                 endif
                 
                 splus = merge(spbot,sptop,vtrans(i,j+1,k).gt.ZERO)
                 savg  = HALF * (spbot + sptop)
                 splus = merge(splus, savg, abs(vtrans(i,j+1,k)) .gt. eps)
                 
                 vlo = u(i,j-1,k,2) + w0_cart_vec(i,j-1,k,2)
                 vhi = u(i,j  ,k,2) + w0_cart_vec(i,j  ,k,2)
                 
                 smbot = s(i,j-1,k,comp) + (HALF - dth*vlo/hy) * slopey(i,j-1,k,1)
                 smtop = s(i,j  ,k,comp) - (HALF + dth*vhi/hy) * slopey(i,j  ,k,1)
                 
                 smtop = merge(s(i,js-1,k,comp),smtop,j.eq.js .and. phys_bc(2,1) .eq. INLET)
                 smbot = merge(s(i,js-1,k,comp),smbot,j.eq.js .and. phys_bc(2,1) .eq. INLET)
                 
                 if (j .eq. js .and. &
                      (phys_bc(2,1).eq.SLIP_WALL.or.phys_bc(2,1).eq.NO_SLIP_WALL)) then
                    if (is_vel .and. (comp .eq. 2)) then
                       smtop = ZERO
                       smbot = ZERO
                    elseif (is_vel .and. (comp .ne. 2)) then
                       smbot = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                       smtop = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
                    else
                       smbot = smtop
                    endif
                 endif
                 
                 sminus = merge(smbot,smtop,vtrans(i,j,k).gt.ZERO)
                 savg   = HALF * (smbot + smtop)
                 sminus = merge(sminus, savg, abs(vtrans(i,j,k)) .gt. eps)
                 
                 st = st - HALF * (vtrans(i,j,k)+vtrans(i,j+1,k))*(splus - sminus) / hy
                 
                 if (is_vel) then
                    ! add the (Utilde . e_r) d w_0 /dr e_r term here

                    if (spherical .eq. 0 .and. comp .eq. 3) then

                       if (k .ge. 0 .and. k .le. nr(n)) then
                          ! wtrans contains w0 so we need to subtract it off
                          st = st - HALF * (wtrans(i,j,k)+wtrans(i,j,k+1)- &
                               w0_cart_vec(i,j,k+1,3)-w0_cart_vec(i,j,k,3))* &
                               (w0_cart_vec(i,j,k+1,3)-w0_cart_vec(i,j,k,3))/hz
                       end if

                    else if (spherical .eq. 1) then

                       ! u/v/wtrans contain w0, so we need to subtract it off.  
                       ! Note w0_cart_vec is cell-centered
                       Ut_dot_er = (HALF*(utrans(i,j,k) + utrans(i+1,j,k)) - &
                                          w0_cart_vec(i,j,k,1))*normal(i,j,k,1) + &
                                   (HALF*(vtrans(i,j,k) + vtrans(i,j+1,k)) - &
                                          w0_cart_vec(i,j,k,2))*normal(i,j,k,2) + &
                                   (HALF*(wtrans(i,j,k) + wtrans(i,j,k+1)) - &
                                          w0_cart_vec(i,j,k,3))*normal(i,j,k,3)
                       
                       if (comp .eq. 1) then
                          st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,1)
                       else if (comp .eq. 2) then
                          st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,2)
                       else if (comp .eq. 3) then
                          st = st - Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,3)
                       endif
                                        
                    endif

                 end if

                 wbardth = dth/hz * ( u(i,j,k,3) + w0_cart_vec(i,j,k,3))
                 
                 if(velpred .eq. 1) then
                    s_d(k+1)= s(i,j,k,comp) + (HALF-wbardth)*slopez(i,j,k,1) + dth*st
                    s_u(k  )= s(i,j,k,comp) - (HALF+wbardth)*slopez(i,j,k,1) + dth*st
                 else
                    s_d(k+1)= s(i,j,k,comp) + (HALF-dth*wmac(i,j,k+1)/hz)*slopez(i,j,k,1)&
                         + dth*st
                    s_u(k  )= s(i,j,k,comp) - (HALF+dth*wmac(i,j,k  )/hz)*slopez(i,j,k,1) &
                         + dth*st
                 endif
                 
              enddo
              
              if (velpred .eq. 1) then
                 do k = ks, ke+1 
                    savg = HALF*(s_d(k) + s_u(k))
                    test = ( (s_d(k) .le. ZERO .and. s_u(k) .ge. ZERO) .or. &
                         (abs(s_d(k) + s_u(k)) .lt. eps) )
                    sedgez(i,j,k,comp)=merge(s_d(k),s_u(k),savg.gt.ZERO)
                    sedgez(i,j,k,comp)=merge(savg,sedgez(i,j,k,comp),test)
                 enddo
              else
                 do k = ks, ke+1 
                    sedgez(i,j,k,comp)=merge(s_d(k),s_u(k),wmac(i,j,k).gt.ZERO)
                    savg = HALF*(s_d(k) + s_u(k))
                    sedgez(i,j,k,comp)=merge(savg,sedgez(i,j,k,comp),abs(wmac(i,j,k)) .lt. eps)
                 enddo
              endif
              
              if (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                 if (is_vel .and. comp .eq. 2) then
                    sedgez(i,j,ks,comp) = ZERO
                 elseif (is_vel .and. comp .ne. 2) then
                    sedgez(i,j,ks,comp) = merge(ZERO,s_u(ks),phys_bc(3,1).eq.NO_SLIP_WALL)
                 else 
                    sedgez(i,j,ks,comp) = s_u(ks)
                 endif
              elseif (phys_bc(3,1) .eq. INLET) then
                 sedgez(i,j,ks,comp) = s(i,j,ks-1,comp)
              elseif (phys_bc(3,1) .eq. OUTLET) then
                 if (is_vel .and. comp.eq.3) then
                    sedgez(i,j,ks,comp) = MIN(s_u(ks),ZERO)
                 else
                    sedgez(i,j,ks,comp) = s_u(ks)
                 end if
              endif
              
              if (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                 if (is_vel .and. comp .eq. 2) then
                    sedgez(i,j,ke+1,comp) = ZERO
                 elseif (is_vel .and. comp .ne. 2) then
                    sedgez(i,j,ke+1,comp) = merge(ZERO,s_d(ke+1),phys_bc(3,2).eq.NO_SLIP_WALL)
                 else 
                    sedgez(i,j,ke+1,comp) = s_d(ke+1)
                 endif
              elseif (phys_bc(3,2) .eq. INLET) then
                 sedgez(i,j,ke+1,comp) = s(i,j,ke+1,comp)
              elseif (phys_bc(3,2) .eq. OUTLET) then
                 if (is_vel .and. comp.eq.3) then
                    sedgez(i,j,ke+1,comp) = MAX(s_d(ke+1),ZERO)
                 else
                    sedgez(i,j,ke+1,comp) = s_d(ke+1)
                 end if
              endif
              
              if (velpred .eq. 1) then
                 do k = ks, ke+1 
                    wmac(i,j,k) = sedgez(i,j,k,3)
                 enddo
              endif
              
           enddo
        enddo
     endif
     
   end subroutine make_edge_state_3d
   
   subroutine make_edge_state_1d(nlevs,s,sedgex,umac,force,dx,dt)

     use geometry, only: r_start_coord, r_end_coord, nr_fine, nr
     use probin_module, only: slope_order
     use bl_constants_module
     
     integer        , intent(in   ) :: nlevs
     real(kind=dp_t), intent(in   ) ::      s(:,0:)
     real(kind=dp_t), intent(inout) :: sedgex(:,0:)
     real(kind=dp_t), intent(in   ) ::   umac(:,0:)
     real(kind=dp_t), intent(in   ) ::  force(:,0:)
     real(kind=dp_t), intent(in   ) :: dx(:),dt
     
     real(kind=dp_t), allocatable::  slopex(:,:)
     real(kind=dp_t), allocatable::  s_l(:,:),s_r(:,:)
     real(kind=dp_t), allocatable:: dxscr(:,:,:)
     real(kind=dp_t) :: dmin,dpls,ds,del,slim,sflag
     real(kind=dp_t) :: ubardth, dth, savg
     real(kind=dp_t) :: abs_eps, eps, umax, u
     
     integer :: r,lo,hi,n
     integer        , parameter :: cen = 1, lim = 2, flag = 3, fromm = 4
     real(kind=dp_t), parameter :: fourthirds = 4.0_dp_t / 3.0_dp_t
        
     allocate(s_l   (nlevs,-1:nr_fine+1))
     allocate(s_r   (nlevs,-1:nr_fine+1))
     allocate(slopex(nlevs, 0:nr_fine-1))
     allocate(dxscr (nlevs, 0:nr_fine-1,4))

     abs_eps = 1.0d-8
     dth = HALF*dt

     ! compute eps based on umax
     umax = ZERO
     do n=1,nlevs
        lo = r_start_coord(n)
        hi = r_end_coord(n)
        do r = lo,hi+1
           umax = max(umax,abs(umac(n,r)))
        end do
     end do
     eps = abs_eps * umax
     
     ! compute slopes
     do n=1,nlevs

        lo = r_start_coord(n)
        hi = r_end_coord(n)
        
        if (slope_order .eq. 0) then
           
           slopex(n,:) = ZERO
           
        else if (slope_order .eq. 2) then
           
           if (n .eq. 1) then ! second order slopes for coarsest level
              
              ! do standard limiting on interior cells
              do r = lo+1,hi-1
                 del = half*(s(n,r+1) - s(n,r-1))
                 dpls = two*(s(n,r+1) - s(n,r  ))
                 dmin = two*(s(n,r  ) - s(n,r-1))
                 slim = min(abs(dpls), abs(dmin))
                 slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                 sflag = sign(one,del)
                 slopex(n,r)= sflag*min(slim,abs(del))
              enddo

              ! set slopes next to domain boundaries to zero
              slopex(n,lo) = ZERO
              slopex(n,hi) = ZERO

           else ! second order slopes for non-coarsest levels
                            
              do r = lo,hi
                 if (r .eq. 0 .or. r .eq. nr(n)-1) then
                    ! set slopes next to domain boundaries to zero
                    slopex(n,r) = ZERO
                 else
                    ! do standard limiting on interior cells
                    del = half*(s(n,r+1) - s(n,r-1))
                    dpls = two*(s(n,r+1) - s(n,r  ))
                    dmin = two*(s(n,r  ) - s(n,r-1))
                    slim = min(abs(dpls), abs(dmin))
                    slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                    sflag = sign(one,del)
                    slopex(n,r)= sflag*min(slim,abs(del))
                 end if
              enddo

           end if
           
        else if (slope_order .eq. 4) then
           
           if (n .eq. 1) then ! fourth order slopes for coarsest level
              
              do r = lo+1,hi-1
                 dxscr(n,r,cen) = half*(s(n,r+1)-s(n,r-1))
                 dpls = two*(s(n,r+1)-s(n,r  ))
                 dmin = two*(s(n,r  )-s(n,r-1))
                 dxscr(n,r,lim)= min(abs(dmin),abs(dpls))
                 dxscr(n,r,lim) = merge(dxscr(n,r,lim),zero,dpls*dmin.gt.ZERO)
                 dxscr(n,r,flag) = sign(one,dxscr(n,r,cen))
                 dxscr(n,r,fromm)= dxscr(n,r,flag)*min(dxscr(n,r,lim),abs(dxscr(n,r,cen)))
              enddo
              
              dxscr(n,lo,fromm) = ZERO
              dxscr(n,hi,fromm) = ZERO
              
              ! fourth order limited slopes on interior
              do r = lo+1,hi-1
                 ds = fourthirds * dxscr(n,r,cen) - &
                      sixth * (dxscr(n,r+1,fromm) + dxscr(n,r-1,fromm))
                 slopex(n,r) = dxscr(n,r,flag)*min(abs(ds),dxscr(n,r,lim))
              enddo

              ! set slopes adjacent to domain boundaries to zero
              slopex(n,lo) = ZERO
              slopex(n,hi) = ZERO

           else ! fourth order slopes for non-coarsest levels

              do r=lo,hi

                 if (r .eq. 0 .or. r .eq. nr(n)-1) then
                    dxscr(n,r,fromm) = ZERO
                 else
                    dxscr(n,r,cen) = half*(s(n,r+1)-s(n,r-1))
                    dpls = two*(s(n,r+1)-s(n,r  ))
                    dmin = two*(s(n,r  )-s(n,r-1))
                    dxscr(n,r,lim)= min(abs(dmin),abs(dpls))
                    dxscr(n,r,lim) = merge(dxscr(n,r,lim),zero,dpls*dmin.gt.ZERO)
                    dxscr(n,r,flag) = sign(one,dxscr(n,r,cen))
                    dxscr(n,r,fromm)= dxscr(n,r,flag)*min(dxscr(n,r,lim),abs(dxscr(n,r,cen)))
                 end if

              end do
              
              do r=lo,hi

                 if (r .eq. 0 .or. r .eq. nr(n)-1) then
                    ! set slopes adjacent to domain boundaries to zero
                    slopex(n,r) = ZERO
                 else if (r .eq. r_start_coord(n) .or. r .eq. r_end_coord(n)) then
                    ! drop order to second-order limited differences
                    del = half*(s(n,r+1) - s(n,r-1))
                    dpls = two*(s(n,r+1) - s(n,r  ))
                    dmin = two*(s(n,r  ) - s(n,r-1))
                    slim = min(abs(dpls), abs(dmin))
                    slim = merge(slim, zero, dpls*dmin.gt.ZERO)
                    sflag = sign(one,del)
                    slopex(n,r)= sflag*min(slim,abs(del))
                 else
                    ! fourth-order limited slopes on interior
                    ds = fourthirds * dxscr(n,r,cen) - &
                         sixth * (dxscr(n,r+1,fromm) + dxscr(n,r-1,fromm))
                    slopex(n,r) = dxscr(n,r,flag)*min(abs(ds),dxscr(n,r,lim))
                 end if
                    
              end do

           end if ! which level
           
        end if ! slope order

     end do ! end compute slopes

     ! compute s_l and s_r
     do n=1,nlevs

        lo = r_start_coord(n)
        hi = r_end_coord(n)

        do r = lo,hi

           u = HALF * (umac(n,r) + umac(n,r+1))
           ubardth = dth*u/dx(n)
        
           s_l(n,r+1)= s(n,r) + (HALF-ubardth)*slopex(n,r) + dth * force(n,r)
           s_r(n,r  )= s(n,r) - (HALF+ubardth)*slopex(n,r) + dth * force(n,r)
           
        end do
        
     end do ! end compute s_l and s_r

     ! compute edge states from s_l and s_r
     do n=1,nlevs

        lo = r_start_coord(n)
        hi = r_end_coord(n)
        
        ! if we are not at the finest level
        ! copy in the s_r and s_l states from the next finer level at the c-f interface
        if (n .ne. nlevs) then
           s_r(n,r_start_coord(n+1)/2) = s_r(n+1,r_start_coord(n+1))
           s_l(n,(r_end_coord(n+1)+1)/2) = s_l(n+1,r_end_coord(n+1)+1)
        end if

        ! if we are not at the coarsest level
        ! copy in the s_l and s_r states from the next coarser level at the c-f interface
        if (n .ne. 1) then
           s_l(n,lo) = s_l(n-1,lo/2)
           s_r(n,hi+1) = s_r(n-1,(hi+1)/2)
        end if

        do r=lo,hi+1
           if (r .eq. 0) then
              sedgex(n,r) = s_r(n,r)
           else if (r .eq. nr(n)) then
              sedgex(n,r) = s_l(n,r)
           else
              sedgex(n,r)=merge(s_l(n,r),s_r(n,r),umac(n,r).gt.ZERO)
              savg = HALF*(s_r(n,r) + s_l(n,r))
              sedgex(n,r)=merge(savg,sedgex(n,r),abs(umac(n,r)) .lt. eps)
           end if
        end do
        
     end do ! end compute edge state from s_l and s_r
     
   end subroutine make_edge_state_1d
   
 end module make_edge_state_module
