! make_edge_scal constructs the edge state of a scalar, using a 
! second-order Taylor expansion in space (through dx/2) and time 
! (though dt/2).   We use only MAC-projected edge velocities in this
! prediction.
!
! We are computing all edge states for each 
! variable.  This is what is done for the final updates of the state 
! variables and velocity.

module make_edge_scal_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_layout_module

  implicit none

  private

  public :: make_edge_scal
  
contains

  subroutine make_edge_scal(nlevs,s,sedge,umac,force,w0,w0_cart_vec,dx, &
                             dt,is_vel,the_bc_level,start_scomp,start_bccomp, &
                             num_comp,mla)

    use bl_prof_module
    use bl_constants_module
    use ml_restriction_module, only : ml_edge_restriction_c

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: sedge(:,:),umac(:,:)
    type(multifab) , intent(in   ) :: force(:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0_cart_vec(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    logical        , intent(in   ) :: is_vel
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    integer        , intent(in   ) :: start_scomp,start_bccomp,num_comp
    type(ml_layout), intent(inout) :: mla

    integer                  :: i,scomp,bccomp,ng,dm,n,lo(s(1)%dim)
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: wtp(:,:,:,:)
    real(kind=dp_t), pointer :: w0p(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_edge_scal")

    dm = s(1)%dim
    ng = s(1)%ng
    
    do n=1,nlevs

       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n),i) ) cycle
          sop  => dataptr(s(n),i)
          sepx => dataptr(sedge(n,1),i)
          sepy => dataptr(sedge(n,2),i)
          ump  => dataptr(umac(n,1),i)
          vmp  => dataptr(umac(n,2),i)
          fp   => dataptr(force(n),i)
          lo   =  lwb(get_box(s(n),i))
          select case (dm)
          case (2)
             do scomp = start_scomp, start_scomp + num_comp - 1
                bccomp = start_bccomp + scomp - start_scomp
                call make_edge_scal_2d(n,sop(:,:,1,:), &
                                        sepx(:,:,1,:), sepy(:,:,1,:), &
                                        ump(:,:,1,1), vmp(:,:,1,1), &
                                        fp(:,:,1,:), w0(n,:), &
                                        lo, dx(n,:), dt, is_vel, &
                                        the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                        the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp:), &
                                        ng, scomp)
             end do
          case (3)
!            wmp  => dataptr(  umac(n,3),i)
!            sepz => dataptr( sedge(n,3),i)
!            w0p  => dataptr(w0_cart_vec(n),i)
!            do scomp = start_scomp, start_scomp + num_comp - 1
!               bccomp = start_bccomp + scomp - start_scomp
!               call make_edge_scal_3d(n,sop(:,:,:,:), &
!                                       sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
!                                       ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
!                                       utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), &
!                                       fp(:,:,:,:), &
!                                       w0(n,:), w0p(:,:,:,:), &
!                                       lo, dx(n,:), dt, &
!                                       the_bc_level(n)%phys_bc_level_array(i,:,:), &
!                                       the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp:), &
!                                       ng, scomp)
!            end do
          end select
       end do

    end do
    !
    ! We call ml_edge_restriction for the output velocity if is_vel .eq. .true.
    ! we do not call ml_edge_restriction for scalars because instead we will call
    ! ml_edge_restriction on the fluxes in mkflux.
    !
    if (is_vel) then
       do n = nlevs,2,-1
          do i = 1, dm
             call ml_edge_restriction_c(sedge(n-1,i),1,sedge(n,i),1,mla%mba%rr(n-1,:),i,dm)
          enddo
       enddo
    end if

    call destroy(bpt)
    
  end subroutine make_edge_scal

  
  subroutine make_edge_scal_2d(n,s,sedgex,sedgey,umac,vmac,force,w0,lo,dx,dt,is_vel, &
                               phys_bc,adv_bc,ng,comp)

    use geometry, only: nr
    use bc_module
    use slope_module
    use bl_constants_module

    integer        , intent(in   ) :: n,lo(:)
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(inout) :: sedgex(lo(1)   :,lo(2)   :,:)
    real(kind=dp_t), intent(inout) :: sedgey(lo(1)   :,lo(2)   :,:)
    real(kind=dp_t), intent(inout) ::   umac(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(inout) ::   vmac(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)- 1:,lo(2)- 1:,:)
    real(kind=dp_t), intent(in   ) ::     w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    logical        , intent(in   ) :: is_vel
    integer        , intent(in   ) :: phys_bc(:,:)
    integer        , intent(in   ) :: adv_bc(:,:,:)
    integer        , intent(in   ) :: ng,comp
    
    ! Local variables
    real(kind=dp_t), allocatable :: slopex(:,:,:),slopey(:,:,:)
    real(kind=dp_t), allocatable :: s_l(:),s_r(:),s_b(:),s_t(:)
    
    real(kind=dp_t) :: hx, hy, dth
    real(kind=dp_t) :: splus,sminus
    real(kind=dp_t) :: savg,st
    real(kind=dp_t) :: uedge,vedge
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
    
       
    umax = abs(umac(is,js))

    do j = js,je; do i = is,ie+1
       umax = max(umax,abs(umac(i,j)))
    end do; end do
    do j = js,je+1; do i = is,ie
       umax = max(umax,abs(vmac(i,j)))
    end do; end do
    
    if (umax .eq. 0.d0) then
       eps = abs_eps
    else
       eps = abs_eps * umax
    endif
    
    !********************************
    ! Loop for edge states on x-edges.
    !********************************

    do j = js,je 
       do i = is-1,ie+1 
             
             spbot = s(i,j  ,comp) + (HALF - dth*vmac(i,j+1)/hy) * slopey(i,j  ,1)
             sptop = s(i,j+1,comp) - (HALF + dth*vmac(i,j+1)/hy) * slopey(i,j+1,1)
             
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
             
             splus = merge(spbot,sptop,vmac(i,j+1).gt.ZERO)
             savg  = HALF * (spbot + sptop)
             splus = merge(splus, savg, abs(vmac(i,j+1)) .gt. eps)
             
             smtop = s(i,j  ,comp) - (HALF + dth*vmac(i,j)/hy) * slopey(i,j  ,1)
             smbot = s(i,j-1,comp) + (HALF - dth*vmac(i,j)/hy) * slopey(i,j-1,1)
             
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

             sminus = merge(smbot,smtop,vmac(i,j).gt.ZERO)
             savg   = HALF * (smbot + smtop)
             sminus = merge(sminus, savg, abs(vmac(i,j)) .gt. eps)
             
             st = force(i,j,comp) - HALF * (vmac(i,j)+vmac(i,j+1))*(splus - sminus) / hy
             
             if (is_vel .and. comp .eq. 2) then
                st = st - HALF * (vmac(i,j)+vmac(i,j+1)-w0(j+1)-w0(j))*(w0(j+1)-w0(j))/hy
             end if

             s_l(i+1)= s(i,j,comp) + (HALF-dth*umac(i+1,j)/hx)*slopex(i,j,1) + dth*st
             s_r(i  )= s(i,j,comp) - (HALF+dth*umac(i  ,j)/hx)*slopex(i,j,1) + dth*st
             
       enddo
          
       do i = is, ie+1 
          sedgex(i,j,comp)=merge(s_l(i),s_r(i),umac(i,j).gt.ZERO)
          savg = HALF*(s_r(i) + s_l(i))
          sedgex(i,j,comp)=merge(savg,sedgex(i,j,comp),abs(umac(i,j)) .lt. eps)
       enddo

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
       
    enddo
    
    !********************************
    ! Loop for edge states on y-edges.
    !********************************

    do i = is, ie 
       do j = js-1, je+1 
             
          splft = s(i,j  ,comp) + (HALF - dth*umac(i+1,j)/hx) * slopex(i  ,j,1)
          sprgt = s(i+1,j,comp) - (HALF + dth*umac(i+1,j)/hx) * slopex(i+1,j,1)
          
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
          
          splus = merge(splft,sprgt,umac(i+1,j).gt.ZERO)
          savg  = HALF * (splft + sprgt)
          splus = merge(splus, savg, abs(umac(i+1,j)) .gt. eps)
          
          smrgt = s(i  ,j,comp) - (HALF + dth*umac(i,j)/hx) * slopex(i  ,j,1)
          smlft = s(i-1,j,comp) + (HALF - dth*umac(i,j)/hx) * slopex(i-1,j,1)
          
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
          
          sminus = merge(smlft,smrgt,umac(i,j).gt.ZERO)
          savg   = HALF * (smlft + smrgt)
          sminus = merge(sminus, savg, abs(umac(i,j)) .gt. eps)
          
          st = force(i,j,comp) - HALF * (umac(i,j)+umac(i+1,j))*(splus - sminus) / hx
          
          if (is_vel .and. comp .eq. 2 .and. j .ge. 0 .and. j .lt. nr(n)) then
             st = st - HALF * s(i,j,comp)*(w0(j+1)-w0(j))/hy
          end if
          
          s_b(j+1)= s(i,j,comp) + (HALF-dth*vmac(i,j+1)/hy)*slopey(i,j,1) + dth*st
          s_t(j  )= s(i,j,comp) - (HALF+dth*vmac(i,j  )/hy)*slopey(i,j,1) + dth*st
          
       enddo
       
       do j = js, je+1 
          sedgey(i,j,comp)=merge(s_b(j),s_t(j),vmac(i,j).gt.ZERO)
          savg = HALF*(s_b(j) + s_t(j))
          sedgey(i,j,comp)=merge(savg,sedgey(i,j,comp),abs(vmac(i,j)) .lt. eps)
       enddo
          
       if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
          sedgey(i,js,comp) = s_t(js)
       elseif (phys_bc(2,1) .eq. INLET) then
          sedgey(i,js,comp) = s(i,js-1,comp)
       elseif (phys_bc(2,1) .eq. OUTLET) then
          sedgey(i,js,comp) = s_t(js)
       endif
          
       if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
          sedgey(i,je+1,comp) = s_b(je+1)
       elseif (phys_bc(2,2) .eq. INLET) then
          sedgey(i,je+1,comp) = s(i,je+1,comp)
       elseif (phys_bc(2,2) .eq. OUTLET) then
          sedgey(i,je+1,comp) = s_b(je+1)
       endif
       
    enddo
    
  end subroutine make_edge_scal_2d
   
 end module make_edge_scal_module
