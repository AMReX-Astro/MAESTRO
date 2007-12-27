module make_edge_state_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: make_edge_state, make_edge_state_1d
  
contains

  subroutine make_edge_state(nlevs,s,u,sedge,umac,utrans,force,w0,w0_cart_vec,dx, &
                             dt,is_vel,the_bc_level,velpred,start_scomp,start_bccomp, &
                             num_comp,mla)

    use bl_prof_module
    use bl_constants_module
    use ml_restriction_module, only : ml_edge_restriction_c

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(in   ) :: s(:),u(:)
    type(multifab) , intent(inout) :: sedge(:,:),umac(:,:)
    type(multifab) , intent(in   ) :: utrans(:,:),force(:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0_cart_vec(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    logical        , intent(in   ) :: is_vel
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    integer        , intent(in   ) :: velpred,start_scomp,start_bccomp,num_comp
    type(ml_layout), intent(inout) :: mla

    ! local
    integer                  :: i,scomp,bccomp,ng,dm,n
    integer                  :: lo(u(1)%dim)
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

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_edge_state")

    dm = u(1)%dim
    ng = s(1)%ng
    
    do n=1,nlevs

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
          lo =  lwb(get_box(s(n),i))
          select case (dm)
          case (2)
             do scomp = start_scomp, start_scomp + num_comp - 1
                bccomp = start_bccomp + scomp - start_scomp
                call make_edge_state_2d(n,sop(:,:,1,:), uop(:,:,1,:), &
                                        sepx(:,:,1,:), sepy(:,:,1,:), &
                                        ump(:,:,1,1), vmp(:,:,1,1), &
                                        utp(:,:,1,1), vtp(:,:,1,1), fp(:,:,1,:), w0(n,:), &
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
             do scomp = start_scomp, start_scomp + num_comp - 1
                bccomp = start_bccomp + scomp - start_scomp
                call make_edge_state_3d(n,sop(:,:,:,:), uop(:,:,:,:), &
                                        sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                        ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                        utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), &
                                        fp(:,:,:,:), &
                                        w0(n,:), w0p(:,:,:,:), &
                                        lo, dx(n,:), dt, is_vel, &
                                        the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                        the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp:), &
                                        velpred, ng, scomp)
             end do
          end select
       end do

    end do ! end loop over levels

    ! we call ml_edge_restriction for the output velocity if is_vel .eq. .true.
    ! we do not call ml_edge_restriction for scalars because instead we will call 
    ! ml_edge_restriction on the fluxes in mkflux
    if(is_vel .and. velpred .eq. 1) then
       do n = nlevs,2,-1
          do i = 1, dm
             call ml_edge_restriction_c(umac(n-1,i),1,umac(n,i),1,mla%mba%rr(n-1,:),i,1)
          enddo
       enddo
    end if
    if(is_vel .and. velpred .eq. 0) then
       do n = nlevs,2,-1
          do i = 1, dm
             call ml_edge_restriction_c(sedge(n-1,i),1,sedge(n,i),1,mla%mba%rr(n-1,:),i,dm)
          enddo
       enddo
    end if

    call destroy(bpt)
    
  end subroutine make_edge_state

  
  subroutine make_edge_state_2d(n,s,u,sedgex,sedgey,umac,vmac,utrans, &
                                vtrans,force,w0,lo,dx,dt,is_vel,phys_bc,adv_bc,velpred, &
                                ng,comp)

    use geometry, only: nr
    use bc_module
    use slope_module
    use bl_constants_module

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
    
    integer :: hi(2),slope_order
    integer :: i,j,is,js,ie,je
    logical :: test
    
    slope_order = 4
    
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
    
    call slopex_2d(s(:,:,comp:),slopex,lo,ng,1,adv_bc,slope_order)
    call slopey_2d(s(:,:,comp:),slopey,lo,ng,1,adv_bc,slope_order)
    
    abs_eps = 1.0d-8
    
    dth = HALF*dt
    
    hx = dx(1)
    hy = dx(2)
    
    if (velpred .eq. 1) then
       
       umax = abs(utrans(is,js))
       do j = js,je
          do i = is,ie+1
             umax = max(umax,abs(utrans(i,j)))
          end do
       end do
       do j = js,je+1
          do i = is,ie
             umax = max(umax,abs(vtrans(i,j)))
          end do
       end do
       
    else 
       
       umax = abs(umac(is,js))
       do j = js,je
          do i = is,ie+1
             umax = max(umax,abs(umac(i,j)))
          end do
       end do
       do j = js,je+1
          do i = is,ie
             umax = max(umax,abs(vmac(i,j)))
          end do
       end do

    end if
    
    if(umax .eq. 0.d0) then
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
             if ((j+2).le.nr(n)) then
                vhi = u(i,j+1,2) + HALF * (w0(j+1)+w0(j+2))
             else
                vhi = u(i,j+1,2) + w0(j+1)
             end if
             
             spbot = s(i,j  ,comp) + (HALF - dth*vlo/hy) * slopey(i,j  ,1)
             sptop = s(i,j+1,comp) - (HALF + dth*vhi/hy) * slopey(i,j+1,1)
             
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
             
             smtop = s(i,j  ,comp) - (HALF + dth*vhi/hy) * slopey(i,j  ,1)
             smbot = s(i,j-1,comp) + (HALF - dth*vlo/hy) * slopey(i,j-1,1)
             
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
                st = st - HALF * (vtrans(i,j)+vtrans(i,j+1))*(w0(j+1)-w0(j))/hy
             end if
             
             ubardth = dth*u(i,j,1)/hx
             
             if(velpred .eq. 1) then
                s_l(i+1)= s(i,j,comp) + (HALF-ubardth)*slopex(i,j,1) + dth*st
                s_r(i  )= s(i,j,comp) - (HALF+ubardth)*slopex(i,j,1) + dth*st
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
             
             splft = s(i,j  ,comp) + (HALF - dth*u(i  ,j,1)/hx) * slopex(i  ,j,1)
             sprgt = s(i+1,j,comp) - (HALF + dth*u(i+1,j,1)/hx) * slopex(i+1,j,1)
             
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
             
             smrgt = s(i  ,j,comp) - (HALF + dth*u(i  ,j,1)/hx) * slopex(i  ,j,1)
             smlft = s(i-1,j,comp) + (HALF - dth*u(i-1,j,1)/hx) * slopex(i-1,j,1)
             
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
             
             if (is_vel .and. comp.eq.2 .and. j.ge.0 .and. j.lt.nr(n)) then
                st = st - HALF * (vtrans(i,j)+vtrans(i,j+1))*(w0(j+1)-w0(j))/hy
             end if
             
             if (j .ge. 0 .and. j.lt.nr(n)) then
                vbardth = dth / hy * ( u(i,j,2) + HALF * (w0(j)+w0(j+1)) )
             else
                vbardth = dth / hy * u(i,j,2) 
             end if
             
             if(velpred .eq. 1) then
                s_b(j+1)= s(i,j,comp) + (HALF-vbardth)*slopey(i,j,1) + dth*st
                s_t(j  )= s(i,j,comp) - (HALF+vbardth)*slopey(i,j,1) + dth*st
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
    
    deallocate(s_l)
    deallocate(s_r)
    deallocate(s_b)
    deallocate(s_t)
    
    deallocate(slopex)
    deallocate(slopey)
    
  end subroutine make_edge_state_2d
  
  
  subroutine make_edge_state_3d(n,s,u,sedgex,sedgey,sedgez,umac, &
                                vmac,wmac,utrans,vtrans,wtrans,force,w0,w0_cart_vec,lo, &
                                dx,dt,is_vel,phys_bc,adv_bc,velpred,ng,comp)

    use bc_module
    use slope_module
    use bl_constants_module
    use geometry, only: spherical, nr

    integer        , intent(in   ) :: n,lo(:)
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
    real(kind=dp_t), intent(in   ) ::     w0(0:)
    real(kind=dp_t), intent(in   ) :: w0_cart_vec(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    logical        , intent(in   ) :: is_vel
    integer        , intent(in   ) :: phys_bc(:,:)
    integer        , intent(in   ) :: adv_bc(:,:,:)
    integer        , intent(in   ) :: velpred
    integer        , intent(in   ) :: ng,comp

    ! Local variables
    real(kind=dp_t), allocatable :: slopex(:,:,:,:),slopey(:,:,:,:),slopez(:,:,:,:)
    real(kind=dp_t), allocatable :: s_l(:),s_r(:),s_b(:),s_t(:),s_u(:),s_d(:)
    
    real(kind=dp_t) :: ubardth, vbardth, wbardth
    real(kind=dp_t) :: hx, hy, hz, dth
    real(kind=dp_t) :: splus,sminus
    real(kind=dp_t) :: savg,st
    real(kind=dp_t) :: ulo,uhi,vlo,vhi,wlo,whi
    real(kind=dp_t) :: sptop,spbot,smtop,smbot,splft,sprgt,smlft,smrgt
    real(kind=dp_t) :: abs_eps, eps, umax
    
    integer :: hi(3),slope_order
    integer :: i,j,k,is,js,ie,je,ks,ke
    logical :: test
    
    slope_order = 4
    
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
       call slopex_2d(s(:,:,k,comp:),slopex(:,:,k,:),lo,ng,1,adv_bc,slope_order)
       call slopey_2d(s(:,:,k,comp:),slopey(:,:,k,:),lo,ng,1,adv_bc,slope_order)
    end do
    call slopez_3d(s(:,:,:,comp:),slopez,lo,ng,1,adv_bc,slope_order)
    
    abs_eps = 1.0d-8
    
    dth = HALF*dt
    
    hx = dx(1)
    hy = dx(2)
    hz = dx(3)
    
    if (velpred .eq. 1) then
       
       umax = abs(utrans(is,js,ks))
       do k = ks,ke
          do j = js,je
             do i = is,ie+1
                umax = max(umax,abs(utrans(i,j,k)))
             end do
          end do
       end do
       do k = ks,ke
          do j = js,je+1
             do i = is,ie
                umax = max(umax,abs(vtrans(i,j,k)))
             end do
          end do
        end do
        do k = ks,ke+1
           do j = js,je
              do i = is,ie
                 umax = max(umax,abs(wtrans(i,j,k)))
              end do
           end do
        end do
        
     else 
        
        umax = abs(umac(is,js,ks))
        do k = ks,ke
           do j = js,je
              do i = is,ie+1
                 umax = max(umax,abs(umac(i,j,k)))
              end do
           end do
        end do
        do k = ks,ke
           do j = js,je+1
              do i = is,ie
                 umax = max(umax,abs(vmac(i,j,k)))
              end do
           end do
        end do
        do k = ks,ke+1
           do j = js,je
              do i = is,ie
                 umax = max(umax,abs(wmac(i,j,k)))
              end do
           end do
        end do

     end if
     
     if(umax .eq. 0.d0) then
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
                 
                 ! NOTE NOTE : THIS IS WRONG FOR SPHERICAL !!
                 if (spherical .eq. 0 .and. is_vel .and. comp.eq.3) then
                    st = st - HALF * (wtrans(i,j,k)+wtrans(i,j,k+1))*(w0(k+1)-w0(k))/hz
                 end if
                 
                 ubardth = dth/hx * ( u(i,j,k,1) + w0_cart_vec(i,j,k,1))
                 
                 if(velpred .eq. 1) then
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
                 
                 ! NOTE NOTE : THIS IS WRONG FOR SPHERICAL !!
                 if (spherical .eq. 0 .and. is_vel .and. comp.eq.3) then
                    st = st - HALF * (wtrans(i,j,k)+wtrans(i,j,k+1))*(w0(k+1)-w0(k))/hz
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
                 
                 ! NOTE NOTE : THIS IS WRONG FOR SPHERICAL !!
                 if (spherical.eq.0.and.is_vel.and.comp.eq.3.and.k.ge.0.and.k.lt.nr(n)) then
                    st = st - HALF * (wtrans(i,j,k)+wtrans(i,j,k+1))*(w0(k+1)-w0(k))/hz
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
     
     deallocate(s_l)
     deallocate(s_r)
     deallocate(s_b)
     deallocate(s_t)
     deallocate(s_d)
     deallocate(s_u)
     
     deallocate(slopex)
     deallocate(slopey)
     deallocate(slopez)
     
   end subroutine make_edge_state_3d
   
   subroutine make_edge_state_1d(n,s,sedgex,umac,force,lo,dx,dt)

     use geometry, only: nr
     use bl_constants_module
     
     integer        , intent(in   ) :: n, lo
     real(kind=dp_t), intent(in   ) ::      s(lo:)
     real(kind=dp_t), intent(inout) :: sedgex(lo:)
     real(kind=dp_t), intent(in   ) ::   umac(lo:)
     real(kind=dp_t), intent(in   ) ::  force(lo:)
     real(kind=dp_t), intent(in   ) :: dx,dt
     
     ! Local variables
     real(kind=dp_t), allocatable::  slopex(:)
     real(kind=dp_t), allocatable::  s_l(:),s_r(:)
     real(kind=dp_t), allocatable:: dxscr(:,:)
     real(kind=dp_t) :: dmin,dpls,ds
     real(kind=dp_t) :: ubardth, dth, savg
     real(kind=dp_t) :: abs_eps, eps, umax, u
     real(kind=dp_t) :: fourthirds
     
     integer :: i,is,ie,hi
     integer, parameter :: cen = 1, lim = 2, flag = 3, fromm = 4
     
     hi = lo + nr(n) - 1
     
     allocate(s_l(lo-1:hi+2),s_r(lo-1:hi+2))
     allocate(slopex(lo:hi))
     
     allocate(dxscr(lo:hi,4))
     
     abs_eps = 1.0d-8
     
     dth = HALF*dt
     
     is = lo
     ie = lo + nr(n) - 1
     
     umax = ZERO
     do i = is,ie+1
        umax = max(umax,abs(umac(i)))
     end do
     
     eps = abs_eps * umax
     
     ! Compute fourth-order slopes
     do i = is+1,ie-1
        dxscr(i,cen) = half*(s(i+1)-s(i-1))
        dmin = two*(s(i  )-s(i-1))
        dpls = two*(s(i+1)-s(i  ))
        dxscr(i,lim)= min(abs(dmin),abs(dpls))
        dxscr(i,lim) = merge(dxscr(i,lim),zero,dpls*dmin.gt.ZERO)
        dxscr(i,flag) = sign(one,dxscr(i,cen))
        dxscr(i,fromm)= dxscr(i,flag)*min(dxscr(i,lim),abs(dxscr(i,cen)))
     enddo
     
     dxscr(is,fromm) = ZERO
     dxscr(ie,fromm) = ZERO
     
     fourthirds = 4.0_dp_t / 3.0_dp_t
     
     do i = is+1,ie-1
        ds = fourthirds * dxscr(i,cen) - sixth * (dxscr(i+1,fromm) + dxscr(i-1,fromm))
        slopex(i) = dxscr(i,flag)*min(abs(ds),dxscr(i,lim))
     enddo
     
     slopex(is) = ZERO
     slopex(ie) = ZERO
     
     ! Use fourth-order slopes to compute edge values
     do i = is,ie
        
        u = HALF * (umac(i) + umac(i+1))
        ubardth = dth*u/dx
        
        s_l(i+1)= s(i) + (HALF-ubardth)*slopex(i) + dth * force(i)
        s_r(i  )= s(i) - (HALF+ubardth)*slopex(i) + dth * force(i)
        
     enddo
     
     sedgex(is  ) = s_r(is  )
     sedgex(ie+1) = s_l(ie+1)
     
     do i = is+1, ie 
        sedgex(i)=merge(s_l(i),s_r(i),umac(i).gt.ZERO)
        savg = HALF*(s_r(i) + s_l(i))
        sedgex(i)=merge(savg,sedgex(i),abs(umac(i)) .lt. eps)
     enddo
     
     deallocate(s_l)
     deallocate(s_r)
     deallocate(slopex)
     deallocate(dxscr)
     
   end subroutine make_edge_state_1d
   
 end module make_edge_state_module
