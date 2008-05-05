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

  subroutine make_edge_scal(nlevs,s,sedge,umac,force, &
                            normal,w0,w0_cart_vec, &
                            dx,dt,is_vel,the_bc_level, &
                            start_scomp,start_bccomp,num_comp,is_conservative,mla)

    use bl_prof_module
    use bl_constants_module
    use geometry
    use variables, only: foextrap_comp
    use fill_3d_module
    use multifab_physbc_module
    use ml_restriction_module, only : ml_edge_restriction_c

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: sedge(:,:)
    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: force(:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0_cart_vec(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    logical        , intent(in   ) :: is_vel
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    integer        , intent(in   ) :: start_scomp,start_bccomp,num_comp
    logical        , intent(in   ) :: is_conservative
    type(ml_layout), intent(inout) :: mla

    integer                  :: i,r,scomp,bccomp,ng,dm,n
    integer                  :: lo(s(1)%dim), hi(s(1)%dim)
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: w0p(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: gw0p(:,:,:,:)

    real(kind=dp_t), allocatable :: gradw0_rad(:)
    type(multifab) :: gradw0_cart


    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_edge_scal")

    dm = s(1)%dim
    ng = s(1)%ng

    if (spherical .eq. 1) then
       allocate (gradw0_rad(0:nr(nlevs)-1))

       ! NOTE: here we are doing the computation at the finest level
       do r = 0, nr(nlevs)-1
          gradw0_rad(r) = (w0(nlevs,r+1) - w0(nlevs,r)) / dr(nlevs)
       enddo
    endif

    do n=1,nlevs

       if (spherical .eq. 1 .and. is_vel) then
          call multifab_build(gradw0_cart, s(n)%la,1,1)

          do i = 1, gradw0_cart%nboxes
             if ( multifab_remote(s(n),i) ) cycle
             gw0p => dataptr(gradw0_cart, i)
             lo = lwb(get_box(gradw0_cart,i))
             hi = upb(get_box(gradw0_cart,i))
             
             call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,1,gradw0_rad,gw0p, &
                                               lo,hi,dx(n,:),gradw0_cart%ng)
             
          enddo
          
          
          ! fill ghost cells for two adjacent grids at the same level
          ! this includes periodic domain boundary ghost cells 
          call multifab_fill_boundary(gradw0_cart)

          ! fill non-periodic domain boundary ghost cells.
          ! NOTE: not sure what the BC should be for gradw0_cart.  Right
          ! now I am just using foextrap_comp.
          call multifab_physbc(gradw0_cart,1,foextrap_comp,1,the_bc_level(n))
       endif

 
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
                call make_edge_scal_2d(n, sop(:,:,1,:), &
                                       sepx(:,:,1,:), sepy(:,:,1,:), &
                                       ump(:,:,1,1), vmp(:,:,1,1), &
                                       fp(:,:,1,:), w0(n,:), &
                                       lo, dx(n,:), dt, is_vel, &
                                       the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                       the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp:), &
                                       ng, scomp, is_conservative)
             end do

          case (3)
            wmp  => dataptr(  umac(n,3),i)
            sepz => dataptr( sedge(n,3),i)
            w0p  => dataptr(w0_cart_vec(n),i)
            do scomp = start_scomp, start_scomp + num_comp - 1
               bccomp = start_bccomp + scomp - start_scomp
               call make_edge_scal_3d(n, sop(:,:,:,:), &
                                      sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                      ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                      fp(:,:,:,:), w0p(:,:,:,:), &
                                      lo, dx(n,:), dt, is_vel, &
                                      the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                      the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp:), &
                                      ng, scomp, is_conservative)
            end do
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

    if (spherical .eq. 1 .and. is_vel) then
       deallocate(gradw0_rad)
       call destroy(gradw0_cart)
    endif

    call destroy(bpt)
    
  end subroutine make_edge_scal

  
  subroutine make_edge_scal_2d(n,s,sedgex,sedgey,umac,vmac, &
                               force,w0,lo,dx,dt,is_vel,phys_bc,adv_bc, &
                               ng,comp,is_conservative)

    use geometry, only: nr
    use bc_module
    use slope_module
    use bl_constants_module
    use probin_module, only: use_new_godunov

    integer        , intent(in   ) :: n,lo(:)
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(inout) :: sedgex(lo(1)   :,lo(2)   :,:)
    real(kind=dp_t), intent(inout) :: sedgey(lo(1)   :,lo(2)   :,:)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(in   ) ::   vmac(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)- 1:,lo(2)- 1:,:)
    real(kind=dp_t), intent(in   ) ::     w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    logical        , intent(in   ) :: is_vel
    integer        , intent(in   ) :: phys_bc(:,:)
    integer        , intent(in   ) :: adv_bc(:,:,:)
    integer        , intent(in   ) :: ng,comp
    logical        , intent(in   ) :: is_conservative

    ! Local variables
    real(kind=dp_t), allocatable :: slopex(:,:,:),slopey(:,:,:)
    real(kind=dp_t), allocatable :: s_l(:),s_r(:),s_b(:),s_t(:)

    real(kind=dp_t) :: hx,hy,dth,splus,sminus
    real(kind=dp_t) :: savg,st
    real(kind=dp_t) :: sptop,spbot,smtop,smbot,splft,sprgt,smlft,smrgt
    real(kind=dp_t) :: abs_eps,eps,umax,dw0drhi,dw0drlo,vtilde

    integer :: hi(2)
    integer :: i,j,is,js,ie,je

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

          if (use_new_godunov .and. is_vel .and. comp .eq. 2) then

             if ((j+2).le.nr(n)) then
                dw0drhi = (w0(j+2)-w0(j+1))/dx(2)
             else
                dw0drhi = (w0(j+1)-w0(j))/dx(2)
             end if

             if(j .ge. 0) then
                dw0drlo = (w0(j+1)-w0(j))/dx(2)
             else
                dw0drlo = (w0(j+2)-w0(j+1))/dx(2)
             end if

             vtilde = vmac(i,j+1) - w0(j+1)

             spbot = spbot - dth*vtilde*dw0drlo
             sptop = sptop - dth*vtilde*dw0drhi

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

          splus = merge(spbot,sptop,vmac(i,j+1).gt.ZERO)
          savg  = HALF * (spbot + sptop)
          splus = merge(splus, savg, abs(vmac(i,j+1)) .gt. eps)

          smtop = s(i,j  ,comp) - (HALF + dth*vmac(i,j)/hy) * slopey(i,j  ,1)
          smbot = s(i,j-1,comp) + (HALF - dth*vmac(i,j)/hy) * slopey(i,j-1,1)

          if (use_new_godunov .and. is_vel .and. comp .eq. 2) then

             if ((j+1).le.nr(n)) then
                dw0drhi = (w0(j+1)-w0(j))/dx(2)
             else
                dw0drhi = (w0(j+2)-w0(j+1))/dx(2)
             end if

             if(j-1 .ge. 0) then
                dw0drlo = (w0(j)-w0(j-1))/dx(2)
             else
                dw0drlo = (w0(j+1)-w0(j))/dx(2)
             end if

             vtilde = vmac(i,j) - w0(j)

             smbot = smbot - dth*vtilde*dw0drlo
             smtop = smtop - dth*vtilde*dw0drhi

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

          sminus = merge(smbot,smtop,vmac(i,j).gt.ZERO)
          savg   = HALF * (smbot + smtop)
          sminus = merge(sminus, savg, abs(vmac(i,j)) .gt. eps)

          if (is_conservative) then
             st = force(i,j,comp) - ( vmac(i,j+1)*splus-vmac(i,j)*sminus ) / hy &
                  - s(i,j,comp)*(umac(i+1,j)-umac(i,j)) / hx
          else
             st = force(i,j,comp) - HALF * (vmac(i,j)+vmac(i,j+1))*(splus - sminus) / hy
          end if

          if (is_vel .and. comp .eq. 2) then
             ! vmac contains w0 so we need to subtract it off
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
          
          if (is_conservative) then
             st = force(i,j,comp) - ( umac(i+1,j)*splus-umac(i,j)*sminus ) /hx &
                  - s(i,j,comp)*(vmac(i,j+1)-vmac(i,j)) / hy
          else
             st = force(i,j,comp) - HALF * (umac(i,j)+umac(i+1,j))*(splus - sminus) / hx
          end if
          
          if (is_vel .and. comp .eq. 2 .and. j .ge. 0 .and. j .lt. nr(n)) then
             ! vmac contains w0 so we need to subtract it off
             st = st - HALF * (vmac(i,j)+vmac(i,j+1)-w0(j+1)-w0(j))*(w0(j+1)-w0(j))/hy
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

  subroutine make_edge_scal_3d(n,s,sedgex,sedgey,sedgez,umac,vmac,wmac,force,w0_cart_vec, &
                               lo,dx,dt,is_vel,phys_bc,adv_bc,ng,comp,is_conservative)

    use geometry, only: nr, spherical
    use bc_module
    use slope_module
    use bl_constants_module

    integer        , intent(in   ) :: n,lo(:)
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(inout) :: sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t), intent(inout) :: sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t), intent(inout) :: sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real(kind=dp_t), intent(in   ) ::   vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real(kind=dp_t), intent(in   ) ::   wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
    real(kind=dp_t), intent(in   ) :: w0_cart_vec(lo(1)-2:,lo(2)-2:,lo(3)-2:,:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    logical        , intent(in   ) :: is_vel
    integer        , intent(in   ) :: phys_bc(:,:)
    integer        , intent(in   ) :: adv_bc(:,:,:)
    integer        , intent(in   ) :: ng,comp
    logical        , intent(in   ) :: is_conservative

    ! Local variables
    real(kind=dp_t), allocatable :: slopex(:,:,:,:),slopey(:,:,:,:),slopez(:,:,:,:)
    real(kind=dp_t), allocatable :: s_l(:),s_r(:),s_b(:),s_t(:),s_u(:),s_d(:)

    real(kind=dp_t) :: hx,hy,hz,dth,splus,sminus
    real(kind=dp_t) :: savg,st
    real(kind=dp_t) :: sptop,spbot,smtop,smbot,splft,sprgt,smlft,smrgt
    real(kind=dp_t) :: abs_eps,eps,umax

    integer :: hi(3)
    integer :: i,j,k,is,js,ks,ie,je,ke

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

    if (umax .eq. 0.d0) then
       eps = abs_eps
    else
       eps = abs_eps * umax
    endif

    !********************************
    ! Loop for edge states on x-edges.
    !********************************

    do k = ks, ke
       do j= js, je
          do i = is-1, ie+1

             spbot = s(i,j  ,k,comp) + (HALF - dth*vmac(i,j+1,k)/hy) * slopey(i,j  ,k,1)
             sptop = s(i,j+1,k,comp) - (HALF + dth*vmac(i,j+1,k)/hy) * slopey(i,j+1,k,1)

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

             splus = merge(spbot,sptop,vmac(i,j+1,k).gt.ZERO)
             savg  = HALF * (spbot + sptop)
             splus = merge(splus, savg, abs(vmac(i,j+1,k)) .gt. eps)

             smtop = s(i,j  ,k,comp) - (HALF + dth*vmac(i,j,k)/hy) * slopey(i,j  ,k,1)
             smbot = s(i,j-1,k,comp) + (HALF - dth*vmac(i,j,k)/hy) * slopey(i,j-1,k,1)

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

             sminus = merge(smbot,smtop,vmac(i,j,k).gt.ZERO)
             savg   = HALF * (smbot + smtop)
             sminus = merge(sminus, savg, abs(vmac(i,j,k)) .gt. eps)

             if (is_conservative) then
                st = force(i,j,k,comp) - ( vmac(i,j+1,k)*splus-vmac(i,j,k)*sminus )  / hy
             else
                st = force(i,j,k,comp) &
                     - HALF * (vmac(i,j,k)+vmac(i,j+1,k))*(splus - sminus) / hy
             end if

             spbot = s(i,j,k  ,comp) + (HALF - dth*wmac(i,j,k+1)/hz) * slopez(i,j,k  ,1)
             sptop = s(i,j,k+1,comp) - (HALF + dth*wmac(i,j,k+1)/hz) * slopez(i,j,k+1,1)
             
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
                 
             splus = merge(spbot,sptop,wmac(i,j,k+1).gt.ZERO)
             savg  = HALF * (spbot + sptop)
             splus = merge(splus, savg, abs(wmac(i,j,k+1)) .gt. eps)
                 
             smtop = s(i,j,k  ,comp) - (HALF + dth*wmac(i,j,k)/hz) * slopez(i,j,k  ,1)
             smbot = s(i,j,k-1,comp) + (HALF - dth*wmac(i,j,k)/hz) * slopez(i,j,k-1,1)
                 
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
             
             sminus = merge(smbot,smtop,wmac(i,j,k).gt.ZERO)
             savg   = HALF * (smbot + smtop)
             sminus = merge(sminus, savg, abs(wmac(i,j,k)) .gt. eps)

             if (is_conservative) then
                st = st - ( wmac(i,j,k+1)*splus-wmac(i,j,k)*sminus ) / hz &
                     - s(i,j,k,comp)*(umac(i+1,j,k)-umac(i,j,k)) / hx
             else
                st = st - HALF * (wmac(i,j,k)+wmac(i,j,k+1))*(splus - sminus) / hz
             end if

             ! NOTE NOTE : THIS IS WRONG FOR SPHERICAL !!
             if (spherical .eq. 0 .and. is_vel .and. comp .eq. 3) then
                ! wmac contains w0 so we need to subtract it off
                st = st - HALF*(wmac(i,j,k)+wmac(i,j,k+1)- &
                     w0_cart_vec(i,j,k+1,3)-w0_cart_vec(i,j,k,3)) &
                     * (w0_cart_vec(i,j,k+1,3)-w0_cart_vec(i,j,k,3))/hz
             end if

             s_l(i+1) = s(i,j,k,comp) + (HALF-dth*umac(i+1,j,k)/hx)*slopex(i,j,k,1) + dth*st
             s_r(i  ) = s(i,j,k,comp) - (HALF+dth*umac(i  ,j,k)/hx)*slopex(i,j,k,1) + dth*st

          end do

          do i = is, ie+1 
             sedgex(i,j,k,comp)=merge(s_l(i),s_r(i),umac(i,j,k).gt.ZERO)
             savg = HALF*(s_r(i) + s_l(i))
             sedgex(i,j,k,comp)=merge(savg,sedgex(i,j,k,comp),abs(umac(i,j,k)) .lt. eps)
          enddo
              
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
          
       end do
    end do

    !********************************
    ! Loop for edge states on y-edges.
    !********************************

    do k = ks, ke 
       do i = is, ie 
          do j = js-1, je+1 

             splft = s(i,j  ,k,comp) + (HALF - dth*umac(i+1,j,k)/hx) * slopex(i  ,j,k,1)
             sprgt = s(i+1,j,k,comp) - (HALF + dth*umac(i+1,j,k)/hx) * slopex(i+1,j,k,1)
             
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
             
             splus = merge(splft,sprgt,umac(i+1,j,k).gt.ZERO)
             savg  = HALF * (splft + sprgt)
             splus = merge(splus, savg, abs(umac(i+1,j,k)) .gt. eps)
             
             smrgt = s(i  ,j,k,comp) - (HALF + dth*umac(i,j,k)/hx) * slopex(i  ,j,k,1)
             smlft = s(i-1,j,k,comp) + (HALF - dth*umac(i,j,k)/hx) * slopex(i-1,j,k,1)
             
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
             
             sminus = merge(smlft,smrgt,umac(i,j,k).gt.ZERO)
             savg   = HALF * (smlft + smrgt)
             sminus = merge(sminus, savg, abs(umac(i,j,k)) .gt. eps)

             if (is_conservative) then
                st = force(i,j,k,comp) - ( umac(i+1,j,k)*splus-umac(i,j,k)*sminus ) / hx
             else
                st = force(i,j,k,comp) &
                     - HALF*(umac(i,j,k)+umac(i+1,j,k)) * (splus - sminus) / hx
             end if

             splft = s(i,j,k  ,comp) + (HALF - dth*wmac(i,j,k+1)/hz) * slopez(i,j,k  ,1)
             sprgt = s(i,j,k+1,comp) - (HALF + dth*wmac(i,j,k+1)/hz) * slopez(i,j,k+1,1)
                 
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
             
             splus = merge(splft,sprgt,wmac(i,j,k+1).gt.ZERO)
             savg  = HALF * (splft + sprgt)
             splus = merge(splus, savg, abs(wmac(i,j,k+1)) .gt. eps)
             
             smrgt = s(i,j,k  ,comp) - (HALF + dth*wmac(i,j,k)/hz) * slopez(i,j,k  ,1)
             smlft = s(i,j,k-1,comp) + (HALF - dth*wmac(i,j,k)/hz) * slopez(i,j,k-1,1)
                 
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
             
             sminus = merge(smlft,smrgt,wmac(i,j,k).gt.ZERO)
             savg   = HALF * (smlft + smrgt)
             sminus = merge(sminus, savg, abs(wmac(i,j,k)) .gt. eps)
                 
             if (is_conservative) then
                st = st - ( wmac(i,j,k+1)*splus-wmac(i,j,k)*sminus ) / hz &
                  - s(i,j,k,comp)*(vmac(i,j+1,k)-vmac(i,j,k)) / hy
             else
                st = st - HALF * (wmac(i,j,k)+wmac(i,j,k+1))*(splus - sminus) / hz
             end if
                 
             ! NOTE NOTE : THIS IS WRONG FOR SPHERICAL !!
             if (spherical .eq. 0 .and. is_vel .and. comp.eq.3) then
                ! wmac contains w0 so we need to subtract it off
                st = st - HALF * (wmac(i,j,k)+wmac(i,j,k+1)- &
                     w0_cart_vec(i,j,k+1,3)-w0_cart_vec(i,j,k,3))* &
                     (w0_cart_vec(i,j,k+1,3)-w0_cart_vec(i,j,k,3))/hz
             end if
                 
             s_b(j+1)= s(i,j,k,comp) + (HALF-dth*vmac(i,j+1,k)/hy)*slopey(i,j,k,1) + dth*st
             s_t(j  )= s(i,j,k,comp) - (HALF+dth*vmac(i,j,  k)/hy)*slopey(i,j,k,1) + dth*st

          end do
         
          do j = js, je+1 
             sedgey(i,j,k,comp)=merge(s_b(j),s_t(j),vmac(i,j,k).gt.ZERO)
             savg = HALF*(s_b(j) + s_t(j))
             sedgey(i,j,k,comp)=merge(savg,sedgey(i,j,k,comp),abs(vmac(i,j,k)) .lt. eps)
          enddo
          
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

       end do
    end do

    !********************************
    ! Loop for edge states on z-edges.
    !********************************

    do j = js, je 
       do i = is, ie 
          do k = ks-1, ke+1

             splft = s(i  ,j,k,comp) + (HALF - dth*umac(i+1,j,k)/hx) * slopex(i  ,j,k,1)
             sprgt = s(i+1,j,k,comp) - (HALF + dth*umac(i+1,j,k)/hx) * slopex(i+1,j,k,1)
                 
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
                 
             splus = merge(splft,sprgt,umac(i+1,j,k).gt.ZERO)
             savg  = HALF * (splft + sprgt)
             splus = merge(splus, savg, abs(umac(i+1,j,k)) .gt. eps)
             
             smlft = s(i-1,j,k,comp) + (HALF - dth*umac(i,j,k)/hx) * slopex(i-1,j,k,1)
             smrgt = s(i  ,j,k,comp) - (HALF + dth*umac(i,j,k)/hx) * slopex(i  ,j,k,1)
             
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
             
             sminus = merge(smlft,smrgt,umac(i,j,k).gt.ZERO)
             savg   = HALF * (smlft + smrgt)
             sminus = merge(sminus, savg, abs(umac(i,j,k)) .gt. eps)
                 
             if (is_conservative) then
                st = force(i,j,k,comp) - ( umac(i+1,j,k)*splus - umac(i,j,k)*sminus ) / hx
             else
                st = force(i,j,k,comp) &
                     - HALF * (umac(i,j,k)+umac(i+1,j,k))*(splus - sminus) / hx
             end if
                 
             spbot = s(i,j  ,k,comp) + (HALF - dth*vmac(i,j+1,k)/hy) * slopey(i,j  ,k,1)
             sptop = s(i,j+1,k,comp) - (HALF + dth*vmac(i,j+1,k)/hy) * slopey(i,j+1,k,1)
                 
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
             
             splus = merge(spbot,sptop,vmac(i,j+1,k).gt.ZERO)
             savg  = HALF * (spbot + sptop)
             splus = merge(splus, savg, abs(vmac(i,j+1,k)) .gt. eps)
                 
             smbot = s(i,j-1,k,comp) + (HALF - dth*vmac(i,j,k)/hy) * slopey(i,j-1,k,1)
             smtop = s(i,j  ,k,comp) - (HALF + dth*vmac(i,j,k)/hy) * slopey(i,j  ,k,1)
                 
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
                 
             sminus = merge(smbot,smtop,vmac(i,j,k).gt.ZERO)
             savg   = HALF * (smbot + smtop)
             sminus = merge(sminus, savg, abs(vmac(i,j,k)) .gt. eps)
                 
             if (is_conservative) then
                st = st - ( vmac(i,j+1,k)*splus-vmac(i,j,k)*sminus ) / hy &
                  - s(i,j,k,comp)*(wmac(i,j,k+1)-wmac(i,j,k)) / hz
             else
                st = st - HALF * (vmac(i,j,k)+vmac(i,j+1,k))*(splus - sminus) / hy
             end if
             
             ! NOTE NOTE : THIS IS WRONG FOR SPHERICAL !!
             if (spherical .eq. 0 .and. is_vel .and. comp .eq. 3) then
                ! wmac contains w0 so we need to subtract it off
                st = st - HALF*(wmac(i,j,k)+wmac(i,j,k+1)- &
                     w0_cart_vec(i,j,k+1,3)-w0_cart_vec(i,j,k,3)) * &
                     (w0_cart_vec(i,j,k+1,3)-w0_cart_vec(i,j,k,3))/hz
             end if
                 
             s_d(k+1)= s(i,j,k,comp) + (HALF-dth*wmac(i,j,k+1)/hz)*slopez(i,j,k,1) + dth*st
             s_u(k  )= s(i,j,k,comp) - (HALF+dth*wmac(i,j,k  )/hz)*slopez(i,j,k,1) + dth*st

          end do

          do k = ks, ke+1 
             sedgez(i,j,k,comp)=merge(s_d(k),s_u(k),wmac(i,j,k).gt.ZERO)
             savg = HALF*(s_d(k) + s_u(k))
             sedgez(i,j,k,comp)=merge(savg,sedgez(i,j,k,comp),abs(wmac(i,j,k)) .lt. eps)
          enddo
          
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
          
       end do
    end do

  end subroutine make_edge_scal_3d

end module make_edge_scal_module
 
