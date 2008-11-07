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

  subroutine make_edge_scal(s,sedge,umac,force,normal,w0,w0mac,dx,dt,is_vel,the_bc_level, &
                            start_scomp,start_bccomp,num_comp,is_conservative,mla)

    use bl_prof_module
    use bl_constants_module
    use geometry, only: spherical, nr_fine, dr, dm, nlevs
    use variables, only: foextrap_comp
    use fill_3d_module
    use multifab_physbc_module
    use ml_restriction_module, only : ml_edge_restriction_c

    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: sedge(:,:)
    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: force(:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    logical        , intent(in   ) :: is_vel
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    integer        , intent(in   ) :: start_scomp,start_bccomp,num_comp
    logical        , intent(in   ) :: is_conservative
    type(ml_layout), intent(in   ) :: mla

    integer                  :: i,r,scomp,bccomp,n
    integer                  :: lo(dm), hi(dm)
    integer                  :: ng_s,ng_se,ng_um,ng_f,ng_w0,ng_n,ng_gw
    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: w0xp(:,:,:,:)
    real(kind=dp_t), pointer :: w0yp(:,:,:,:)
    real(kind=dp_t), pointer :: w0zp(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: gw0p(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)

    real(kind=dp_t), allocatable :: gradw0_rad(:)
    type(multifab) :: gradw0_cart

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_edge_scal")

    ng_s = s(1)%ng
    ng_se = sedge(1,1)%ng
    ng_um = umac(1,1)%ng
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

       call multifab_build(gradw0_cart, s(n)%la,1,1)

       if (spherical .eq. 1 .and. is_vel) then

          do i = 1, gradw0_cart%nboxes
             if ( multifab_remote(s(n),i) ) cycle
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
          hi   =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             do scomp = start_scomp, start_scomp + num_comp - 1
                bccomp = start_bccomp + scomp - start_scomp
                call make_edge_scal_2d(n, sop(:,:,1,:), ng_s, &
                                       sepx(:,:,1,:), sepy(:,:,1,:), ng_se, &
                                       ump(:,:,1,1), vmp(:,:,1,1), ng_um, &
                                       fp(:,:,1,:), ng_f, w0(n,:), &
                                       lo, hi, dx(n,:), dt, is_vel, &
                                       the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                       the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp:), &
                                       scomp, is_conservative)
             end do

          case (3)
            wmp  => dataptr(  umac(n,3),i)
            sepz => dataptr( sedge(n,3),i)
            w0xp => dataptr(w0mac(n,1),i)
            w0yp => dataptr(w0mac(n,2),i)
            w0zp => dataptr(w0mac(n,3),i)
            np   => dataptr(normal(n), i)
            gw0p => dataptr(gradw0_cart, i)
            ng_gw = gradw0_cart%ng
            do scomp = start_scomp, start_scomp + num_comp - 1
               bccomp = start_bccomp + scomp - start_scomp
               call make_edge_scal_3d(n, sop(:,:,:,:), ng_s, &
                                      sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), ng_se, &
                                      ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_um, &
                                      fp(:,:,:,:), ng_f, np(:,:,:,:), ng_n, w0(n,:), &
                                      w0xp(:,:,:,1), w0yp(:,:,:,1), w0zp(:,:,:,1), ng_w0, &
                                      gw0p(:,:,:,1), ng_gw, lo, hi, dx(n,:), dt, is_vel, &
                                      the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                      the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp:), &
                                      scomp, is_conservative)
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
    if (is_vel) then
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
    
  end subroutine make_edge_scal

  
  subroutine make_edge_scal_2d(n,s,ng_s,sedgex,sedgey,ng_se,umac,vmac,ng_um, &
                               force,ng_f,w0,lo,hi,dx,dt,is_vel,phys_bc,adv_bc, &
                               comp,is_conservative)

    use geometry, only: nr
    use bc_module
    use slope_module
    use bl_constants_module
    use variables, only: rel_eps

    integer        , intent(in   ) :: lo(:),hi(:),n,ng_s,ng_se,ng_um,ng_f
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s :,lo(2)-ng_s :,:)
    real(kind=dp_t), intent(inout) :: sedgex(lo(1)-ng_se:,lo(2)-ng_se:,:)
    real(kind=dp_t), intent(inout) :: sedgey(lo(1)-ng_se:,lo(2)-ng_se:,:)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::   vmac(lo(1)-ng_um:,lo(2)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,:)
    real(kind=dp_t), intent(in   ) ::     w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    logical        , intent(in   ) :: is_vel
    integer        , intent(in   ) :: phys_bc(:,:)
    integer        , intent(in   ) :: adv_bc(:,:,:)
    integer        , intent(in   ) :: comp
    logical        , intent(in   ) :: is_conservative

    ! Local variables
    real(kind=dp_t) :: slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1)
    real(kind=dp_t) :: slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1)
    real(kind=dp_t) :: s_l(lo(1)-1:hi(1)+2)
    real(kind=dp_t) :: s_r(lo(1)-1:hi(1)+2)
    real(kind=dp_t) :: s_b(lo(2)-1:hi(2)+2)
    real(kind=dp_t) :: s_t(lo(2)-1:hi(2)+2)

    real(kind=dp_t) :: hx,hy,dt2,dt4,savg

    integer :: i,j,is,js,ie,je

    ! these correspond to s_L^x, etc.
    real(kind=dp_t), allocatable:: slx(:,:),srx(:,:)
    real(kind=dp_t), allocatable:: sly(:,:),sry(:,:)

    ! these correspond to s_{\i-\half\e_x}^x, etc.
    real(kind=dp_t), allocatable:: simhx(:,:),simhy(:,:)

    ! these correspond to \mathrm{sedge}_L^x, etc.
    real(kind=dp_t), allocatable:: sedgelx(:,:),sedgerx(:,:)
    real(kind=dp_t), allocatable:: sedgely(:,:),sedgery(:,:)

    ! Normal predictor states.
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse direction
    allocate(slx  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(srx  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(simhx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1))

    allocate(sly  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1))
    allocate(sry  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1))
    allocate(simhy(lo(1)-1:hi(1)+1,lo(2):hi(2)+1))

    ! Final edge states.
    ! lo:hi+1 in the normal direction
    ! lo:hi in the transverse direction
    allocate(sedgelx(lo(1):hi(1)+1,lo(2):hi(2)))
    allocate(sedgerx(lo(1):hi(1)+1,lo(2):hi(2)))
    allocate(sedgely(lo(1):hi(1),lo(2):hi(2)+1))
    allocate(sedgery(lo(1):hi(1),lo(2):hi(2)+1))

    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)

    call slopex_2d(s(:,:,comp:),slopex,lo,hi,ng_s,1,adv_bc)
    call slopey_2d(s(:,:,comp:),slopey,lo,hi,ng_s,1,adv_bc)

    dt2 = HALF*dt
    dt4 = dt/4.0d0

    hx = dx(1)
    hy = dx(2)
    
    !******************************************************************
    ! Create s_{\i-\half\e_x}^x, etc.
    !******************************************************************
    
    ! loop over appropriate x-faces
    do j=js-1,je+1
       do i=is,ie+1
          ! make slx, srx with 1D extrapolation
          slx(i,j) = s(i-1,j,comp) + (HALF - dt2*umac(i,j)/hx)*slopex(i-1,j,1)
          srx(i,j) = s(i  ,j,comp) - (HALF + dt2*umac(i,j)/hx)*slopex(i  ,j,1)

          ! impose lo side bc's
          if(i .eq. is) then
             slx(i,j) = merge(s(is-1,j,comp),slx(i,j),phys_bc(1,1) .eq. INLET)
             srx(i,j) = merge(s(is-1,j,comp),srx(i,j),phys_bc(1,1) .eq. INLET)
             if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                if(is_vel .and. comp .eq. 1) then
                   slx(i,j) = ZERO
                   srx(i,j) = ZERO
                else if(is_vel .and. comp .ne. 1) then
                   slx(i,j) = merge(ZERO,srx(i,j),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   srx(i,j) = merge(ZERO,srx(i,j),phys_bc(1,1) .eq. NO_SLIP_WALL)
                else
                   slx(i,j) = srx(i,j)
                endif
             endif
          endif

          ! impose hi side bc's
          if(i .eq. ie+1) then
             slx(i,j) = merge(s(ie+1,j,comp),slx(i,j),phys_bc(1,2) .eq. INLET)
             srx(i,j) = merge(s(ie+1,j,comp),srx(i,j),phys_bc(1,2) .eq. INLET)
             if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                if (is_vel .and. comp .eq. 1) then
                   slx(i,j) = ZERO
                   srx(i,j) = ZERO
                else if (is_vel .and. comp .ne. 1) then
                   slx(i,j) = merge(ZERO,slx(i,j),phys_bc(1,2).eq.NO_SLIP_WALL)
                   srx(i,j) = merge(ZERO,slx(i,j),phys_bc(1,2).eq.NO_SLIP_WALL)
                else
                   srx(i,j) = slx(i,j)
                endif
             endif
          endif

          ! make simhx by solving Riemann problem
          simhx(i,j) = merge(slx(i,j),srx(i,j),umac(i,j) .gt. ZERO)
          savg = HALF*(slx(i,j)+srx(i,j))
          simhx(i,j) = merge(simhx(i,j),savg,abs(umac(i,j)) .gt. rel_eps)
       enddo
    enddo

    ! loop over appropriate y-faces
    do j=js,je+1
       do i=is-1,ie+1
          ! make sly, sry with 1D extrapolation
          sly(i,j) = s(i,j-1,comp) + (HALF - dt2*vmac(i,j)/hy)*slopey(i,j-1,1)
          sry(i,j) = s(i,j  ,comp) - (HALF + dt2*vmac(i,j)/hy)*slopey(i,j  ,1)

          ! impose lo side bc's
          if(j .eq. js) then
             sly(i,j) = merge(s(is,j-1,comp),sly(i,j),phys_bc(2,1) .eq. INLET)
             sry(i,j) = merge(s(is,j-1,comp),sry(i,j),phys_bc(2,1) .eq. INLET)
             if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                if(is_vel .and. comp .eq. 2) then
                   sly(i,j) = ZERO
                   sry(i,j) = ZERO
                else if(is_vel .and. comp .ne. 2) then
                   sly(i,j) = merge(ZERO,sry(i,j),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   sry(i,j) = merge(ZERO,sry(i,j),phys_bc(2,1) .eq. NO_SLIP_WALL)
                else
                   sly(i,j) = sry(i,j)
                endif
             endif
          endif

          ! impose hi side bc's
          if(j .eq. je+1) then
             sly(i,j) = merge(s(i,je+1,comp),sly(i,j),phys_bc(2,2) .eq. INLET)
             sry(i,j) = merge(s(i,je+1,comp),sry(i,j),phys_bc(2,2) .eq. INLET)
             if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                if (is_vel .and. comp .eq. 2) then
                   sly(i,j) = ZERO
                   sry(i,j) = ZERO
                else if (is_vel .and. comp .ne. 2) then
                   sly(i,j) = merge(ZERO,sly(i,j),phys_bc(2,2).eq.NO_SLIP_WALL)
                   sry(i,j) = merge(ZERO,sly(i,j),phys_bc(2,2).eq.NO_SLIP_WALL)
                else
                   sry(i,j) = sly(i,j)
                endif
             endif
          endif

          ! make simhy by solving Riemann problem
          simhy(i,j) = merge(sly(i,j),sry(i,j),vmac(i,j) .gt. ZERO)
          savg = HALF*(sly(i,j)+sry(i,j))
          simhy(i,j) = merge(simhy(i,j),savg,abs(vmac(i,j)) .gt. rel_eps)
       enddo
    enddo

    !******************************************************************
    ! Create sedgelx, etc.
    !******************************************************************

    ! loop over appropriate x-faces
    do j=js,je
       do i=is,ie+1
          ! make sedgelx, sedgerx
          if(is_conservative) then
             sedgelx(i,j) = slx(i,j) &
                  - (dt2/hy)*(simhy(i-1,j+1)*vmac(i-1,j+1) - simhy(i-1,j)*vmac(i-1,j)) &
                  - (dt2/hx)*s(i-1,j,comp)*(umac(i  ,j)-umac(i-1,j)) &
                  + dt2*force(i-1,j,comp)
             sedgerx(i,j) = srx(i,j) &
                  - (dt2/hy)*(simhy(i  ,j+1)*vmac(i  ,j+1) - simhy(i  ,j)*vmac(i  ,j)) &
                  - (dt2/hx)*s(i  ,j,comp)*(umac(i+1,j)-umac(i  ,j)) &
                  + dt2*force(i  ,j,comp)
          else
             sedgelx(i,j) = slx(i,j) &
                  - (dt4/hy)*(vmac(i-1,j+1)+vmac(i-1,j))*(simhy(i-1,j+1)-simhy(i-1,j)) &
                  + dt2*force(i-1,j,comp)
             sedgerx(i,j) = srx(i,j) &
                  - (dt4/hy)*(vmac(i  ,j+1)+vmac(i  ,j))*(simhy(i  ,j+1)-simhy(i  ,j)) &
                  + dt2*force(i  ,j,comp)
          endif

          ! add the (Utilde . e_r) d w_0 /dr e_r term here
          ! vmac contains w0 so we need to subtract it off
          if (is_vel .and. comp .eq. 2) then
             sedgelx(i,j) = sedgelx(i,j) &
                  -(dt4/hy)*(vmac(i-1,j)-w0(j)+vmac(i-1,j+1)-w0(j+1))*(w0(j+1)-w0(j))
             sedgerx(i,j) = sedgerx(i,j) &
                  -(dt4/hy)*(vmac(i  ,j)-w0(j)+vmac(i  ,j+1)-w0(j+1))*(w0(j+1)-w0(j))
          end if

          ! make sedgex by solving Riemann problem
          ! boundary conditions enforced outside of i,j loop
          sedgex(i,j,comp) = merge(sedgelx(i,j),sedgerx(i,j),umac(i,j) .gt. ZERO)
          savg = HALF*(sedgelx(i,j)+sedgerx(i,j))
          sedgex(i,j,comp) = merge(sedgex(i,j,comp),savg,abs(umac(i,j)) .gt. rel_eps)
       enddo
    enddo

    ! sedgex boundary conditions
    do j=js,je
       ! lo side
       if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
          if (is_vel .and. comp .eq. 1) then
             sedgex(is,j,comp) = ZERO
          elseif (is_vel .and. comp .ne. 1) then
             sedgex(is,j,comp) = merge(ZERO,sedgerx(is,j),phys_bc(1,1).eq.NO_SLIP_WALL)
          else
             sedgex(is,j,comp) = sedgerx(is,j)
          endif
       elseif (phys_bc(1,1) .eq. INLET) then
          sedgex(is,j,comp) = s(is-1,j,comp)
       elseif (phys_bc(1,1) .eq. OUTLET) then
          if (is_vel .and. comp.eq.1) then
             sedgex(is,j,comp) = MIN(sedgerx(is,j),ZERO)
          else
             sedgex(is,j,comp) = sedgerx(is,j)
          end if
       endif

       ! hi side
       if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
          if (is_vel .and. comp .eq. 1) then
             sedgex(ie+1,j,comp) = ZERO
          else if (is_vel .and. comp .ne. 1) then
             sedgex(ie+1,j,comp) = merge(ZERO,sedgelx(ie+1,j),phys_bc(1,2).eq.NO_SLIP_WALL)
          else 
             sedgex(ie+1,j,comp) = sedgelx(ie+1,j)
          endif
       elseif (phys_bc(1,2) .eq. INLET) then
          sedgex(ie+1,j,comp) = s(ie+1,j,comp)
       elseif (phys_bc(1,2) .eq. OUTLET) then
          if (is_vel .and. comp.eq.1) then
             sedgex(ie+1,j,comp) = MAX(sedgelx(ie+1,j),ZERO)
          else
             sedgex(ie+1,j,comp) = sedgelx(ie+1,j)
          end if
       endif
    enddo

    ! loop over appropriate y-faces
    do j=js,je+1
       do i=is,ie
          ! make sedgely, sedgery
          if(is_conservative) then
             sedgely(i,j) = sly(i,j) &
                  - (dt2/hx)*(simhx(i+1,j-1)*umac(i+1,j-1) - simhx(i,j-1)*umac(i,j-1)) &
                  - (dt2/hy)*s(i,j-1,comp)*(vmac(i,j  )-vmac(i,j-1)) &
                  + dt2*force(i,j-1,comp)
             sedgery(i,j) = sry(i,j) &
                  - (dt2/hx)*(simhx(i+1,j  )*umac(i+1,j  ) - simhx(i,j  )*umac(i,j  )) &
                  - (dt2/hy)*s(i,j  ,comp)*(vmac(i,j+1)-vmac(i,j  )) &
                  + dt2*force(i,j  ,comp)
          else
             sedgely(i,j) = sly(i,j) &
                  - (dt4/hx)*(umac(i+1,j-1)+umac(i,j-1))*(simhx(i+1,j-1)-simhx(i,j-1)) &
                  + dt2*force(i,j-1,comp)
             sedgery(i,j) = sry(i,j) &
                  - (dt4/hx)*(umac(i+1,j  )+umac(i,j  ))*(simhx(i+1,j  )-simhx(i,j  )) &
                  + dt2*force(i,j  ,comp)
          endif

          ! add the (Utilde . e_r) d w_0 /dr e_r term here
          ! vmac contains w0 so we need to subtract it off
          if (is_vel .and. comp .eq. 2) then
             if (j .eq. 0) then
                ! sedgely unchanged since dw_0 / dr = 0
                sedgery(i,j) = sedgery(i,j) &
                     -(dt4/hy)*(vmac(i,j+1)-w0(j+1)+vmac(i,j  )-w0(j  ))*(w0(j+1)-w0(j  ))
             else if (j .eq. nr(n)) then
                sedgely(i,j) = sedgely(i,j) &
                     -(dt4/hy)*(vmac(i,j  )-w0(j  )+vmac(i,j-1)-w0(j-1))*(w0(j  )-w0(j-1))
                ! sedgery unchanged since dw_0 / dr = 0
             else
                sedgely(i,j) = sedgely(i,j) &
                     -(dt4/hy)*(vmac(i,j  )-w0(j  )+vmac(i,j-1)-w0(j-1))*(w0(j  )-w0(j-1))
                sedgery(i,j) = sedgery(i,j) &
                     -(dt4/hy)*(vmac(i,j+1)-w0(j+1)+vmac(i,j  )-w0(j  ))*(w0(j+1)-w0(j  ))
             end if
          end if

          ! make sedgey by solving Riemann problem
          ! boundary conditions enforced outside of i,j loop
          sedgey(i,j,comp) = merge(sedgely(i,j),sedgery(i,j),vmac(i,j) .gt. ZERO)
          savg = HALF*(sedgely(i,j)+sedgery(i,j))
          sedgey(i,j,comp) = merge(sedgey(i,j,comp),savg,abs(vmac(i,j)) .gt. rel_eps)
       enddo
    enddo

    ! sedgey boundary conditions
    do i=is,ie
       ! lo side
       if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
          if (is_vel .and. comp .eq. 2) then
             sedgey(i,js,comp) = ZERO
          elseif (is_vel .and. comp .ne. 2) then
             sedgey(i,js,comp) = merge(ZERO,sedgery(i,js),phys_bc(2,1).eq.NO_SLIP_WALL)
          else 
             sedgey(i,js,comp) = sedgery(i,js)
          endif
       elseif (phys_bc(2,1) .eq. INLET) then
          sedgey(i,js,comp) = s(i,js-1,comp)
       elseif (phys_bc(2,1) .eq. OUTLET) then
          if (is_vel .and. comp.eq.2) then
             sedgey(i,js,comp) = MIN(sedgery(i,js),ZERO)
          else
             sedgey(i,js,comp) = sedgery(i,js)
          end if
       endif

       ! hi side
       if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
          if (is_vel .and. comp .eq. 2) then
             sedgey(i,je+1,comp) = ZERO
          elseif (is_vel .and. comp .ne. 2) then
             sedgey(i,je+1,comp) = merge(ZERO,sedgely(i,je+1),phys_bc(2,2).eq.NO_SLIP_WALL)
          else 
             sedgey(i,je+1,comp) = sedgely(i,je+1)
          endif
       elseif (phys_bc(2,2) .eq. INLET) then
          sedgey(i,je+1,comp) = s(i,je+1,comp)
       elseif (phys_bc(2,2) .eq. OUTLET) then
          if (is_vel .and. comp.eq.2) then
             sedgey(i,je+1,comp) = MAX(sedgely(i,je+1),ZERO)
          else
             sedgey(i,je+1,comp) = sedgely(i,je+1)
          end if
       endif
    enddo

    deallocate(slx,srx,sly,sry,simhx,simhy,sedgelx,sedgerx,sedgely,sedgery)

  end subroutine make_edge_scal_2d

  subroutine make_edge_scal_3d(n,s,ng_s,sedgex,sedgey,sedgez,ng_se,umac,vmac,wmac,ng_um, &
                               force,ng_f,normal,ng_n,w0,w0macx,w0macy,w0macz,ng_w0, &
                               gradw0_cart,ng_gw,lo,hi,dx,dt,is_vel,phys_bc,adv_bc,comp, &
                               is_conservative)

    use geometry, only: spherical, nr
    use bc_module
    use slope_module
    use bl_constants_module
    use variables, only: rel_eps

    integer        , intent(in   ) :: n,lo(:),hi(:),ng_s,ng_se,ng_um,ng_f,ng_w0,ng_n,ng_gw
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)
    real(kind=dp_t), intent(inout) :: sedgex(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(inout) :: sedgey(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(inout) :: sedgez(lo(1)-ng_se:,lo(2)-ng_se:,lo(3)-ng_se:,:)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::   vmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::   wmac(lo(1)-ng_um:,lo(2)-ng_um:,lo(3)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    real(kind=dp_t), intent(in   ) :: normal(lo(1)-ng_n :,lo(2)-ng_n :,lo(3)-ng_n :,:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real(kind=dp_t), intent(in   ) :: gradw0_cart(lo(1)-ng_gw:,lo(2)-ng_gw:,lo(3)-ng_gw:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    logical        , intent(in   ) :: is_vel
    integer        , intent(in   ) :: phys_bc(:,:)
    integer        , intent(in   ) :: adv_bc(:,:,:)
    integer        , intent(in   ) :: comp
    logical        , intent(in   ) :: is_conservative

    ! Local variables
    real(kind=dp_t) :: slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1)
    real(kind=dp_t) :: slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1)
    real(kind=dp_t) :: slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1)
    real(kind=dp_t) :: s_l(lo(1)-1:hi(1)+2)
    real(kind=dp_t) :: s_r(lo(1)-1:hi(1)+2)
    real(kind=dp_t) :: s_b(lo(2)-1:hi(2)+2)
    real(kind=dp_t) :: s_t(lo(2)-1:hi(2)+2)
    real(kind=dp_t) :: s_u(lo(3)-1:hi(3)+2)
    real(kind=dp_t) :: s_d(lo(3)-1:hi(3)+2)

    real(kind=dp_t) :: hx,hy,hz,dt2,dt3,dt4,dt6,splus,sminus
    real(kind=dp_t) :: savg,st
    real(kind=dp_t) :: sptop,spbot,smtop,smbot,splft,sprgt,smlft,smrgt
    real(kind=dp_t) :: Ut_dot_er

    integer :: i,j,k,is,js,ks,ie,je,ke

    ! these correspond to s_L^x, etc.
    real(kind=dp_t), allocatable:: slx(:,:,:),srx(:,:,:)
    real(kind=dp_t), allocatable:: sly(:,:,:),sry(:,:,:)
    real(kind=dp_t), allocatable:: slz(:,:,:),srz(:,:,:)

    ! these correspond to s_{\i-\half\e_x}^x, etc.
    real(kind=dp_t), allocatable:: simhx(:,:,:),simhy(:,:,:),simhz(:,:,:)

    ! these correspond to s_L^{x|y}, etc.
    real(kind=dp_t), allocatable:: slxy(:,:,:),srxy(:,:,:),slxz(:,:,:),srxz(:,:,:)
    real(kind=dp_t), allocatable:: slyx(:,:,:),sryx(:,:,:),slyz(:,:,:),sryz(:,:,:)
    real(kind=dp_t), allocatable:: slzx(:,:,:),srzx(:,:,:),slzy(:,:,:),srzy(:,:,:)

    ! these correspond to s_{\i-\half\e_x}^{x|y}, etc.
    real(kind=dp_t), allocatable:: simhxy(:,:,:),simhxz(:,:,:)
    real(kind=dp_t), allocatable:: simhyx(:,:,:),simhyz(:,:,:)
    real(kind=dp_t), allocatable:: simhzx(:,:,:),simhzy(:,:,:)

    ! these correspond to \mathrm{sedge}_L^x, etc.
    real(kind=dp_t), allocatable:: sedgelx(:,:,:),sedgerx(:,:,:)
    real(kind=dp_t), allocatable:: sedgely(:,:,:),sedgery(:,:,:)
    real(kind=dp_t), allocatable:: sedgelz(:,:,:),sedgerz(:,:,:)

    ! Normal predictor states.
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    allocate(slx  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(srx  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(simhx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

    allocate(sly  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(sry  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(simhy(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))

    allocate(slz  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(srz  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(simhz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))

    ! These are transverse terms.  The size allocation is tricky.
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    ! lo-1:hi+1 in unused direction
    allocate(slxy  (lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))
    allocate(srxy  (lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))
    allocate(simhxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))

    allocate(slxz  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))
    allocate(srxz  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))
    allocate(simhxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))

    allocate(slyx  (lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(sryx  (lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(simhyx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))

    allocate(slyz  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(sryz  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(simhyz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))

    allocate(slzx  (lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(srzx  (lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(simhzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))

    allocate(slzy  (lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))
    allocate(srzy  (lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))
    allocate(simhzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))

    ! Final edge states.
    ! lo:hi+1 in the normal direction
    ! lo:hi in the transverse directions
    allocate(sedgelx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)))
    allocate(sedgerx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)))
    allocate(sedgely(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(sedgery(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(sedgelz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1))
    allocate(sedgerz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1))

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

    dt2 = HALF*dt
    dt3 = dt/3.0d0
    dt4 = dt/4.0d0
    dt6 = dt/6.0d0
    
    hx = dx(1)
    hy = dx(2)
    hz = dx(3)
    
    !******************************************************************
    ! Create s_{\i-\half\e_x}^x, etc.
    !******************************************************************
    
    ! loop over appropriate x-faces
    do k=ks-1,ke+1
       do j=js-1,je+1
          do i=is,ie+1
             ! make slx, srx with 1D extrapolation
             slx(i,j,k) = s(i-1,j,k,comp) + (HALF - dt2*umac(i,j,k)/hx)*slopex(i-1,j,k,1)
             srx(i,j,k) = s(i  ,j,k,comp) - (HALF + dt2*umac(i,j,k)/hx)*slopex(i  ,j,k,1)

             ! impose lo side bc's
             if(i .eq. is) then
                slx(i,j,k) = merge(s(is-1,j,k,comp),slx(i,j,k),phys_bc(1,1) .eq. INLET)
                srx(i,j,k) = merge(s(is-1,j,k,comp),srx(i,j,k),phys_bc(1,1) .eq. INLET)
                if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                   if(is_vel .and. comp .eq. 1) then
                      slx(i,j,k) = ZERO
                      srx(i,j,k) = ZERO
                   else if(is_vel .and. comp .ne. 1) then
                      slx(i,j,k) = merge(ZERO,srx(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                      srx(i,j,k) = merge(ZERO,srx(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   else
                      slx(i,j,k) = srx(i,j,k)
                   endif
                endif
             endif

             ! impose hi side bc's
             if(i .eq. ie+1) then
                slx(i,j,k) = merge(s(ie+1,j,k,comp),slx(i,j,k),phys_bc(1,2) .eq. INLET)
                srx(i,j,k) = merge(s(ie+1,j,k,comp),srx(i,j,k),phys_bc(1,2) .eq. INLET)
                if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                   if (is_vel .and. comp .eq. 1) then
                      slx(i,j,k) = ZERO
                      srx(i,j,k) = ZERO
                   else if (is_vel .and. comp .ne. 1) then
                      slx(i,j,k) = merge(ZERO,slx(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                      srx(i,j,k) = merge(ZERO,slx(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                   else
                      srx(i,j,k) = slx(i,j,k)
                   endif
                endif
             endif

             ! make simhx by solving Riemann problem
             simhx(i,j,k) = merge(slx(i,j,k),srx(i,j,k),umac(i,j,k) .gt. ZERO)
             savg = HALF*(slx(i,j,k)+srx(i,j,k))
             simhx(i,j,k) = merge(simhx(i,j,k),savg,abs(umac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo

    ! loop over appropriate y-faces
    do k=ks-1,ke+1
       do j=js,je+1
          do i=is-1,ie+1
             ! make sly, sry with 1D extrapolation
             sly(i,j,k) = s(i,j-1,k,comp) + (HALF - dt2*vmac(i,j,k)/hy)*slopey(i,j-1,k,1)
             sry(i,j,k) = s(i,j  ,k,comp) - (HALF + dt2*vmac(i,j,k)/hy)*slopey(i,j  ,k,1)

             ! impose lo side bc's
             if(j .eq. js) then
                sly(i,j,k) = merge(s(is,j-1,k,comp),sly(i,j,k),phys_bc(2,1) .eq. INLET)
                sry(i,j,k) = merge(s(is,j-1,k,comp),sry(i,j,k),phys_bc(2,1) .eq. INLET)
                if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                   if(is_vel .and. comp .eq. 2) then
                      sly(i,j,k) = ZERO
                      sry(i,j,k) = ZERO
                   else if(is_vel .and. comp .ne. 2) then
                      sly(i,j,k) = merge(ZERO,sry(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                      sry(i,j,k) = merge(ZERO,sry(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   else
                      sly(i,j,k) = sry(i,j,k)
                   endif
                endif
             endif

             ! impose hi side bc's
             if(j .eq. je+1) then
                sly(i,j,k) = merge(s(i,je+1,k,comp),sly(i,j,k),phys_bc(2,2) .eq. INLET)
                sry(i,j,k) = merge(s(i,je+1,k,comp),sry(i,j,k),phys_bc(2,2) .eq. INLET)
                if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                   if (is_vel .and. comp .eq. 2) then
                      sly(i,j,k) = ZERO
                      sry(i,j,k) = ZERO
                   else if (is_vel .and. comp .ne. 2) then
                      sly(i,j,k) = merge(ZERO,sly(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                      sry(i,j,k) = merge(ZERO,sly(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                   else
                      sry(i,j,k) = sly(i,j,k)
                   endif
                endif
             endif

             ! make simhy by solving Riemann problem
             simhy(i,j,k) = merge(sly(i,j,k),sry(i,j,k),vmac(i,j,k) .gt. ZERO)
             savg = HALF*(sly(i,j,k)+sry(i,j,k))
             simhy(i,j,k) = merge(simhy(i,j,k),savg,abs(vmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo

    ! loop over appropriate z-faces
    do k=ks,ke+1
       do j=js-1,je+1
          do i=is-1,ie+1
             ! make slz, srz with 1D extrapolation
             slz(i,j,k) = s(i,j,k-1,comp) + (HALF - dt2*wmac(i,j,k)/hz)*slopez(i,j,k-1,1)
             srz(i,j,k) = s(i,j,k  ,comp) - (HALF + dt2*wmac(i,j,k)/hz)*slopez(i,j,k  ,1)

             ! impose lo side bc's
             if(k .eq. ks) then
                slz(i,j,k) = merge(s(is,j,k-1,comp),slz(i,j,k),phys_bc(3,1) .eq. INLET)
                srz(i,j,k) = merge(s(is,j,k-1,comp),srz(i,j,k),phys_bc(3,1) .eq. INLET)
                if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                   if(is_vel .and. comp .eq. 3) then
                      slz(i,j,k) = ZERO
                      srz(i,j,k) = ZERO
                   else if(is_vel .and. comp .ne. 3) then
                      slz(i,j,k) = merge(ZERO,srz(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                      srz(i,j,k) = merge(ZERO,srz(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                   else
                      slz(i,j,k) = srz(i,j,k)
                   endif
                endif
             endif

             ! impose hi side bc's
             if(k .eq. ke+1) then
                slz(i,j,k) = merge(s(i,j,ke+1,comp),slz(i,j,k),phys_bc(3,2) .eq. INLET)
                srz(i,j,k) = merge(s(i,j,ke+1,comp),srz(i,j,k),phys_bc(3,2) .eq. INLET)
                if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                   if (is_vel .and. comp .eq. 3) then
                      slz(i,j,k) = ZERO
                      srz(i,j,k) = ZERO
                   else if (is_vel .and. comp .ne. 3) then
                      slz(i,j,k) = merge(ZERO,slz(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                      srz(i,j,k) = merge(ZERO,slz(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                   else
                      srz(i,j,k) = slz(i,j,k)
                   endif
                endif
             endif

             ! make simhz by solving Riemann problem
             simhz(i,j,k) = merge(slz(i,j,k),srz(i,j,k),wmac(i,j,k) .gt. ZERO)
             savg = HALF*(slz(i,j,k)+srz(i,j,k))
             simhz(i,j,k) = merge(simhz(i,j,k),savg,abs(wmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo

    !******************************************************************
    ! Create s_{\i-\half\e_x}^{x|y}, etc.
    !******************************************************************

    ! loop over appropriate xy faces
    do k=ks-1,ke+1
       do j=js,je
          do i=is,ie+1
             ! make slxy, srxy by updating 1D extrapolation
             if(is_conservative) then
                slxy(i,j,k) = slx(i,j,k) &
                     - (dt3/hy)*(simhy(i-1,j+1,k)*vmac(i-1,j+1,k) - simhy(i-1,j,k)*vmac(i-1,j,k))
                srxy(i,j,k) = srx(i,j,k) &
                     - (dt3/hy)*(simhy(i  ,j+1,k)*vmac(i  ,j+1,k) - simhy(i  ,j,k)*vmac(i  ,j,k))
             else
                slxy(i,j,k) = slx(i,j,k) &
                     - (dt6/hy)*(vmac(i-1,j+1,k)+vmac(i-1,j,k))*(simhy(i-1,j+1,k)-simhy(i-1,j,k))
                srxy(i,j,k) = srx(i,j,k) &
                     - (dt6/hy)*(vmac(i  ,j+1,k)+vmac(i  ,j,k))*(simhy(i  ,j+1,k)-simhy(i  ,j,k))
             endif

             ! impose lo side bc's
             if(i .eq. is) then
                slxy(i,j,k) = merge(s(is-1,j,k,comp),slxy(i,j,k),phys_bc(1,1) .eq. INLET)
                srxy(i,j,k) = merge(s(is-1,j,k,comp),srxy(i,j,k),phys_bc(1,1) .eq. INLET)
                if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                   if(is_vel .and. comp .eq. 1) then
                      slxy(i,j,k) = ZERO
                      srxy(i,j,k) = ZERO
                   else if(is_vel .and. comp .ne. 1) then
                      slxy(i,j,k) = merge(ZERO,srxy(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                      srxy(i,j,k) = merge(ZERO,srxy(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   else
                      slxy(i,j,k) = srxy(i,j,k)
                   endif
                endif
             endif

             ! impose hi side bc's
             if(i .eq. ie+1) then
                slxy(i,j,k) = merge(s(ie+1,j,k,comp),slxy(i,j,k),phys_bc(1,2) .eq. INLET)
                srxy(i,j,k) = merge(s(ie+1,j,k,comp),srxy(i,j,k),phys_bc(1,2) .eq. INLET)
                if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                   if (is_vel .and. comp .eq. 1) then
                      slxy(i,j,k) = ZERO
                      srxy(i,j,k) = ZERO
                   else if (is_vel .and. comp .ne. 1) then
                      slxy(i,j,k) = merge(ZERO,slxy(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                      srxy(i,j,k) = merge(ZERO,slxy(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                   else
                      srxy(i,j,k) = slxy(i,j,k)
                   endif
                endif
             endif

             ! make simhxy by solving Riemann problem
             simhxy(i,j,k) = merge(slxy(i,j,k),srxy(i,j,k),umac(i,j,k) .gt. ZERO)
             savg = HALF*(slxy(i,j,k)+srxy(i,j,k))
             simhxy(i,j,k) = merge(simhxy(i,j,k),savg,abs(umac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo

    ! loop over appropriate xz faces
    do k=ks,ke
       do j=js-1,je+1
          do i=is,ie+1
             ! make slxz, srxz by updating 1D extrapolation
             if(is_conservative) then
                slxz(i,j,k) = slx(i,j,k) &
                     - (dt3/hz)*(simhz(i-1,j,k+1)*wmac(i-1,j,k+1) - simhz(i-1,j,k)*wmac(i-1,j,k))
                srxz(i,j,k) = srx(i,j,k) &
                     - (dt3/hz)*(simhz(i  ,j,k+1)*wmac(i  ,j,k+1) - simhz(i  ,j,k)*wmac(i  ,j,k))
             else
                slxz(i,j,k) = slx(i,j,k) &
                     - (dt6/hz)*(wmac(i-1,j,k+1)+wmac(i-1,j,k))*(simhz(i-1,j,k+1)-simhz(i-1,j,k))
                srxz(i,j,k) = srx(i,j,k) &
                     - (dt6/hz)*(wmac(i  ,j,k+1)+wmac(i  ,j,k))*(simhz(i  ,j,k+1)-simhz(i  ,j,k))
             endif

             ! impose lo side bc's
             if(i .eq. is) then
                slxz(i,j,k) = merge(s(is-1,j,k,comp),slxz(i,j,k),phys_bc(1,1) .eq. INLET)
                srxz(i,j,k) = merge(s(is-1,j,k,comp),srxz(i,j,k),phys_bc(1,1) .eq. INLET)
                if(phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
                   if(is_vel .and. comp .eq. 1) then
                      slxz(i,j,k) = ZERO
                      srxz(i,j,k) = ZERO
                   else if(is_vel .and. comp .ne. 1) then
                      slxz(i,j,k) = merge(ZERO,srxz(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                      srxz(i,j,k) = merge(ZERO,srxz(i,j,k),phys_bc(1,1) .eq. NO_SLIP_WALL)
                   else
                      slxz(i,j,k) = srxz(i,j,k)
                   endif
                endif
             endif

             ! impose hi side bc's
             if(i .eq. ie+1) then
                slxz(i,j,k) = merge(s(ie+1,j,k,comp),slxz(i,j,k),phys_bc(1,2) .eq. INLET)
                srxz(i,j,k) = merge(s(ie+1,j,k,comp),srxz(i,j,k),phys_bc(1,2) .eq. INLET)
                if(phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
                   if (is_vel .and. comp .eq. 1) then
                      slxz(i,j,k) = ZERO
                      srxz(i,j,k) = ZERO
                   else if (is_vel .and. comp .ne. 1) then
                      slxz(i,j,k) = merge(ZERO,slxz(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                      srxz(i,j,k) = merge(ZERO,slxz(i,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
                   else
                      srxz(i,j,k) = slxz(i,j,k)
                   endif
                endif
             endif

             ! make simhxz by solving Riemann problem
             simhxz(i,j,k) = merge(slxz(i,j,k),srxz(i,j,k),umac(i,j,k) .gt. ZERO)
             savg = HALF*(slxz(i,j,k)+srxz(i,j,k))
             simhxz(i,j,k) = merge(simhxz(i,j,k),savg,abs(umac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo

    ! loop over appropriate yx faces
    do k=ks-1,ke+1
       do j=js,je+1
          do i=is,ie
             ! make slyx, sryx by updating 1D extrapolation
             if(is_conservative) then
                slyx(i,j,k) = sly(i,j,k) &
                     - (dt3/hx)*(simhx(i+1,j-1,k)*umac(i+1,j-1,k) - simhx(i,j-1,k)*umac(i,j-1,k))
                sryx(i,j,k) = sry(i,j,k) &
                     - (dt3/hx)*(simhx(i+1,j  ,k)*umac(i+1,j  ,k) - simhx(i,j  ,k)*umac(i,j  ,k))
             else
                slyx(i,j,k) = sly(i,j,k) &
                     - (dt6/hx)*(umac(i+1,j-1,k)+umac(i,j-1,k))*(simhx(i+1,j-1,k)-simhx(i,j-1,k))
                sryx(i,j,k) = sry(i,j,k) &
                     - (dt6/hx)*(umac(i+1,j  ,k)+umac(i,j  ,k))*(simhx(i+1,j  ,k)-simhx(i,j  ,k))
             endif

             ! impose lo side bc's
             if(j .eq. js) then
                slyx(i,j,k) = merge(s(is,j-1,k,comp),slyx(i,j,k),phys_bc(2,1) .eq. INLET)
                sryx(i,j,k) = merge(s(is,j-1,k,comp),sryx(i,j,k),phys_bc(2,1) .eq. INLET)
                if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                   if(is_vel .and. comp .eq. 2) then
                      slyx(i,j,k) = ZERO
                      sryx(i,j,k) = ZERO
                   else if(is_vel .and. comp .ne. 2) then
                      slyx(i,j,k) = merge(ZERO,sryx(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                      sryx(i,j,k) = merge(ZERO,sryx(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   else
                      slyx(i,j,k) = sryx(i,j,k)
                   endif
                endif
             endif

             ! impose hi side bc's
             if(j .eq. je+1) then
                slyx(i,j,k) = merge(s(i,je+1,k,comp),slyx(i,j,k),phys_bc(2,2) .eq. INLET)
                sryx(i,j,k) = merge(s(i,je+1,k,comp),sryx(i,j,k),phys_bc(2,2) .eq. INLET)
                if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                   if (is_vel .and. comp .eq. 2) then
                      slyx(i,j,k) = ZERO
                      sryx(i,j,k) = ZERO
                   else if (is_vel .and. comp .ne. 2) then
                      slyx(i,j,k) = merge(ZERO,slyx(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                      sryx(i,j,k) = merge(ZERO,slyx(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                   else
                      sryx(i,j,k) = slyx(i,j,k)
                   endif
                endif
             endif

             ! make simhyx by solving Riemann problem
             simhyx(i,j,k) = merge(slyx(i,j,k),sryx(i,j,k),vmac(i,j,k) .gt. ZERO)
             savg = HALF*(slyx(i,j,k)+sryx(i,j,k))
             simhyx(i,j,k) = merge(simhyx(i,j,k),savg,abs(vmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo

    ! loop over appropriate yz faces
    do k=ks,ke
       do j=js,je+1
          do i=is-1,ie+1
             ! make slyz, sryz by updating 1D extrapolation
             if(is_conservative) then
                slyz(i,j,k) = sly(i,j,k) &
                     - (dt3/hz)*(simhz(i,j-1,k+1)*wmac(i,j-1,k+1) - simhz(i,j-1,k)*wmac(i,j-1,k))
                sryz(i,j,k) = sry(i,j,k) &
                     - (dt3/hz)*(simhz(i,j  ,k+1)*wmac(i,j  ,k+1) - simhz(i,j  ,k)*wmac(i,j  ,k))
             else
                slyz(i,j,k) = sly(i,j,k) &
                     - (dt6/hz)*(wmac(i,j-1,k+1)+wmac(i,j-1,k))*(simhz(i,j-1,k+1)-simhz(i,j-1,k))
                sryz(i,j,k) = sry(i,j,k) &
                     - (dt6/hz)*(wmac(i,j  ,k+1)+wmac(i,j  ,k))*(simhz(i,j  ,k+1)-simhz(i,j  ,k))
             endif

             ! impose lo side bc's
             if(j .eq. js) then
                slyz(i,j,k) = merge(s(is,j-1,k,comp),slyz(i,j,k),phys_bc(2,1) .eq. INLET)
                sryz(i,j,k) = merge(s(is,j-1,k,comp),sryz(i,j,k),phys_bc(2,1) .eq. INLET)
                if(phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
                   if(is_vel .and. comp .eq. 2) then
                      slyz(i,j,k) = ZERO
                      sryz(i,j,k) = ZERO
                   else if(is_vel .and. comp .ne. 2) then
                      slyz(i,j,k) = merge(ZERO,sryz(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                      sryz(i,j,k) = merge(ZERO,sryz(i,j,k),phys_bc(2,1) .eq. NO_SLIP_WALL)
                   else
                      slyz(i,j,k) = sryz(i,j,k)
                   endif
                endif
             endif

             ! impose hi side bc's
             if(j .eq. je+1) then
                slyz(i,j,k) = merge(s(i,je+1,k,comp),slyz(i,j,k),phys_bc(2,2) .eq. INLET)
                sryz(i,j,k) = merge(s(i,je+1,k,comp),sryz(i,j,k),phys_bc(2,2) .eq. INLET)
                if(phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
                   if (is_vel .and. comp .eq. 2) then
                      slyz(i,j,k) = ZERO
                      sryz(i,j,k) = ZERO
                   else if (is_vel .and. comp .ne. 2) then
                      slyz(i,j,k) = merge(ZERO,slyz(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                      sryz(i,j,k) = merge(ZERO,slyz(i,j,k),phys_bc(2,2).eq.NO_SLIP_WALL)
                   else
                      sryz(i,j,k) = slyz(i,j,k)
                   endif
                endif
             endif

             ! make simhyz by solving Riemann problem
             simhyz(i,j,k) = merge(slyz(i,j,k),sryz(i,j,k),vmac(i,j,k) .gt. ZERO)
             savg = HALF*(slyz(i,j,k)+sryz(i,j,k))
             simhyz(i,j,k) = merge(simhyz(i,j,k),savg,abs(vmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo

    ! loop over appropriate zx faces
    do k=ks,ke+1
       do j=js-1,je+1
          do i=is,ie
             ! make slzx, srzx by updating 1D extrapolation
             if(is_conservative) then
                slzx(i,j,k) = slz(i,j,k) &
                     - (dt3/hx)*(simhx(i+1,j,k-1)*umac(i+1,j,k-1) - simhx(i,j,k-1)*umac(i,j,k-1))
                srzx(i,j,k) = srz(i,j,k) &
                     - (dt3/hx)*(simhx(i+1,j,k  )*umac(i+1,j,k  ) - simhx(i,j,k  )*umac(i,j,k  ))
             else
                slzx(i,j,k) = slz(i,j,k) &
                     - (dt6/hx)*(umac(i+1,j,k-1)+umac(i,j,k-1))*(simhx(i+1,j,k-1)-simhx(i,j,k-1))
                srzx(i,j,k) = srz(i,j,k) &
                     - (dt6/hx)*(umac(i+1,j,k  )+umac(i,j,k  ))*(simhx(i+1,j,k  )-simhx(i,j,k  ))
             endif

             ! impose lo side bc's
             if(k .eq. ks) then
                slzx(i,j,k) = merge(s(is,j,k-1,comp),slzx(i,j,k),phys_bc(3,1) .eq. INLET)
                srzx(i,j,k) = merge(s(is,j,k-1,comp),srzx(i,j,k),phys_bc(3,1) .eq. INLET)
                if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                   if(is_vel .and. comp .eq. 3) then
                      slzx(i,j,k) = ZERO
                      srzx(i,j,k) = ZERO
                   else if(is_vel .and. comp .ne. 3) then
                      slzx(i,j,k) = merge(ZERO,srzx(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                      srzx(i,j,k) = merge(ZERO,srzx(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                   else
                      slzx(i,j,k) = srzx(i,j,k)
                   endif
                endif
             endif

             ! impose hi side bc's
             if(k .eq. ke+1) then
                slzx(i,j,k) = merge(s(i,j,ke+1,comp),slzx(i,j,k),phys_bc(3,2) .eq. INLET)
                srzx(i,j,k) = merge(s(i,j,ke+1,comp),srzx(i,j,k),phys_bc(3,2) .eq. INLET)
                if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                   if (is_vel .and. comp .eq. 3) then
                      slzx(i,j,k) = ZERO
                      srzx(i,j,k) = ZERO
                   else if (is_vel .and. comp .ne. 3) then
                      slzx(i,j,k) = merge(ZERO,slzx(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                      srzx(i,j,k) = merge(ZERO,slzx(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                   else
                      srzx(i,j,k) = slzx(i,j,k)
                   endif
                endif
             endif

             ! make simhzx by solving Riemann problem
             simhzx(i,j,k) = merge(slzx(i,j,k),srzx(i,j,k),wmac(i,j,k) .gt. ZERO)
             savg = HALF*(slzx(i,j,k)+srzx(i,j,k))
             simhzx(i,j,k) = merge(simhzx(i,j,k),savg,abs(wmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo

    ! loop over appropriate zy faces
    do k=ks,ke+1
       do j=js,je
          do i=is-1,ie+1
             ! make slzy, srzy by updating 1D extrapolation
             if(is_conservative) then
                slzy(i,j,k) = slz(i,j,k) &
                     - (dt3/hy)*(simhy(i,j+1,k-1)*vmac(i,j+1,k-1) - simhy(i,j,k-1)*vmac(i,j,k-1))
                srzy(i,j,k) = srz(i,j,k) &
                     - (dt3/hy)*(simhy(i,j+1,k  )*vmac(i,j+1,k  ) - simhy(i,j,k  )*vmac(i,j,k  ))
             else
                slzy(i,j,k) = slz(i,j,k) &
                     - (dt6/hy)*(vmac(i,j+1,k-1)+vmac(i,j,k-1))*(simhy(i,j+1,k-1)-simhy(i,j,k-1))
                srzy(i,j,k) = srz(i,j,k) &
                     - (dt6/hy)*(vmac(i,j+1,k  )+vmac(i,j,k  ))*(simhy(i,j+1,k  )-simhy(i,j,k  ))
             endif

             ! impose lo side bc's
             if(k .eq. ks) then
                slzy(i,j,k) = merge(s(is,j,k-1,comp),slzy(i,j,k),phys_bc(3,1) .eq. INLET)
                srzy(i,j,k) = merge(s(is,j,k-1,comp),srzy(i,j,k),phys_bc(3,1) .eq. INLET)
                if(phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
                   if(is_vel .and. comp .eq. 3) then
                      slzy(i,j,k) = ZERO
                      srzy(i,j,k) = ZERO
                   else if(is_vel .and. comp .ne. 3) then
                      slzy(i,j,k) = merge(ZERO,srzy(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                      srzy(i,j,k) = merge(ZERO,srzy(i,j,k),phys_bc(3,1) .eq. NO_SLIP_WALL)
                   else
                      slzy(i,j,k) = srzy(i,j,k)
                   endif
                endif
             endif

             ! impose hi side bc's
             if(k .eq. ke+1) then
                slzy(i,j,k) = merge(s(i,j,ke+1,comp),slzy(i,j,k),phys_bc(3,2) .eq. INLET)
                srzy(i,j,k) = merge(s(i,j,ke+1,comp),srzy(i,j,k),phys_bc(3,2) .eq. INLET)
                if(phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
                   if (is_vel .and. comp .eq. 3) then
                      slzy(i,j,k) = ZERO
                      srzy(i,j,k) = ZERO
                   else if (is_vel .and. comp .ne. 3) then
                      slzy(i,j,k) = merge(ZERO,slzy(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                      srzy(i,j,k) = merge(ZERO,slzy(i,j,k),phys_bc(3,2).eq.NO_SLIP_WALL)
                   else
                      srzy(i,j,k) = slzy(i,j,k)
                   endif
                endif
             endif

             ! make simhzy by solving Riemann problem
             simhzy(i,j,k) = merge(slzy(i,j,k),srzy(i,j,k),wmac(i,j,k) .gt. ZERO)
             savg = HALF*(slzy(i,j,k)+srzy(i,j,k))
             simhzy(i,j,k) = merge(simhzy(i,j,k),savg,abs(wmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo

    !******************************************************************
    ! Create sedgelx, etc.
    !******************************************************************

    ! loop over appropriate x-faces
    do k=ks,ke
       do j=js,je
          do i=is,ie+1
             ! make sedgelx, sedgerx
             if(is_conservative) then
                sedgelx(i,j,k) = slx(i,j,k) &
                     - (dt2/hy)*(simhyz(i-1,j+1,k  )*vmac(i-1,j+1,k  ) - simhyz(i-1,j,k)*vmac(i-1,j,k)) &
                     - (dt2/hz)*(simhzy(i-1,j  ,k+1)*wmac(i-1,j  ,k+1) - simhzy(i-1,j,k)*wmac(i-1,j,k)) &
                     - (dt2/hx)*s(i-1,j,k,comp)*(umac(i  ,j,k)-umac(i-1,j,k)) &
                     + dt2*force(i-1,j,k,comp)
                sedgerx(i,j,k) = srx(i,j,k) &
                     - (dt2/hy)*(simhyz(i  ,j+1,k  )*vmac(i  ,j+1,  k) - simhyz(i  ,j,k)*vmac(i  ,j,k)) &
                     - (dt2/hz)*(simhzy(i  ,j  ,k+1)*wmac(i  ,j  ,k+1) - simhzy(i  ,j,k)*wmac(i  ,j,k)) &
                     - (dt2/hx)*s(i  ,j,k,comp)*(umac(i+1,j,k)-umac(i  ,j,k)) &
                     + dt2*force(i  ,j,k,comp)
             else
                sedgelx(i,j,k) = slx(i,j,k) &
                     - (dt4/hy)*(vmac(i-1,j+1,k  )+vmac(i-1,j,k))*(simhyz(i-1,j+1,k  )-simhyz(i-1,j,k)) &
                     - (dt4/hz)*(wmac(i-1,j  ,k+1)+wmac(i-1,j,k))*(simhzy(i-1,j  ,k+1)-simhzy(i-1,j,k)) &
                     + dt2*force(i-1,j,k,comp)
                sedgerx(i,j,k) = srx(i,j,k) &
                     - (dt4/hy)*(vmac(i  ,j+1,k  )+vmac(i  ,j,k))*(simhyz(i  ,j+1,k  )-simhyz(i  ,j,k)) &
                     - (dt4/hz)*(wmac(i  ,j  ,k+1)+wmac(i  ,j,k))*(simhzy(i  ,j  ,k+1)-simhzy(i  ,j,k)) &
                     + dt2*force(i  ,j,k,comp)
             endif

             ! make sedgex by solving Riemann problem
             ! boundary conditions enforced outside of i,j,k loop
             sedgex(i,j,k,comp) = merge(sedgelx(i,j,k),sedgerx(i,j,k),umac(i,j,k) .gt. ZERO)
             savg = HALF*(sedgelx(i,j,k)+sedgerx(i,j,k))
             sedgex(i,j,k,comp) = merge(sedgex(i,j,k,comp),savg,abs(umac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo

    ! sedgex boundary conditions
    do k=ks,ke
       do j=js,je
          ! lo side
          if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
             if (is_vel .and. comp .eq. 1) then
                sedgex(is,j,k,comp) = ZERO
             elseif (is_vel .and. comp .ne. 1) then
                sedgex(is,j,k,comp) = merge(ZERO,sedgerx(is,j,k),phys_bc(1,1).eq.NO_SLIP_WALL)
             else
                sedgex(is,j,k,comp) = sedgerx(is,j,k)
             endif
          elseif (phys_bc(1,1) .eq. INLET) then
             sedgex(is,j,k,comp) = s(is-1,j,k,comp)
          elseif (phys_bc(1,1) .eq. OUTLET) then
             if (is_vel .and. comp.eq.1) then
                sedgex(is,j,k,comp) = MIN(sedgerx(is,j,k),ZERO)
             else
                sedgex(is,j,k,comp) = sedgerx(is,j,k)
             end if
          endif

          ! hi side
          if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
             if (is_vel .and. comp .eq. 1) then
                sedgex(ie+1,j,k,comp) = ZERO
             else if (is_vel .and. comp .ne. 1) then
                sedgex(ie+1,j,k,comp) = merge(ZERO,sedgelx(ie+1,j,k),phys_bc(1,2).eq.NO_SLIP_WALL)
             else 
                sedgex(ie+1,j,k,comp) = sedgelx(ie+1,j,k)
             endif
          elseif (phys_bc(1,2) .eq. INLET) then
             sedgex(ie+1,j,k,comp) = s(ie+1,j,k,comp)
          elseif (phys_bc(1,2) .eq. OUTLET) then
             if (is_vel .and. comp.eq.1) then
                sedgex(ie+1,j,k,comp) = MAX(sedgelx(ie+1,j,k),ZERO)
             else
                sedgex(ie+1,j,k,comp) = sedgelx(ie+1,j,k)
             end if
          endif
       enddo
    enddo

    ! loop over appropriate y-faces
    do k=ks,ke
       do j=js,je+1
          do i=is,ie
             ! make sedgely, sedgery
             if(is_conservative) then
                sedgely(i,j,k) = sly(i,j,k) &
                     - (dt2/hx)*(simhxz(i+1,j-1,k  )*umac(i+1,j-1,k  ) - simhxz(i,j-1,k)*umac(i,j-1,k)) &
                     - (dt2/hz)*(simhzx(i  ,j-1,k+1)*wmac(i  ,j-1,k+1) - simhzx(i,j-1,k)*wmac(i,j-1,k)) &
                     - (dt2/hy)*s(i,j-1,k,comp)*(vmac(i,j  ,k)-vmac(i,j-1,k)) &
                     + dt2*force(i,j-1,k,comp)
                sedgery(i,j,k) = sry(i,j,k) &
                     - (dt2/hx)*(simhxz(i+1,j  ,k  )*umac(i+1,j  ,k  ) - simhxz(i,j  ,k)*umac(i,j  ,k)) &
                     - (dt2/hz)*(simhzx(i  ,j  ,k+1)*wmac(i  ,j  ,k+1) - simhzx(i,j  ,k)*wmac(i,j  ,k)) &
                     - (dt2/hy)*s(i,j  ,k,comp)*(vmac(i,j+1,k)-vmac(i,j  ,k)) &
                     + dt2*force(i,j  ,k,comp)
             else
                sedgely(i,j,k) = sly(i,j,k) &
                     - (dt4/hx)*(umac(i+1,j-1,k  )+umac(i,j-1,k))*(simhxz(i+1,j-1,k  )-simhxz(i,j-1,k)) &
                     - (dt4/hz)*(wmac(i  ,j-1,k+1)+wmac(i,j-1,k))*(simhzx(i  ,j-1,k+1)-simhzx(i,j-1,k)) &
                     + dt2*force(i,j-1,k,comp)
                sedgery(i,j,k) = sry(i,j,k) &
                     - (dt4/hx)*(umac(i+1,j  ,k  )+umac(i,j  ,k))*(simhxz(i+1,j  ,k  )-simhxz(i,j  ,k)) &
                     - (dt4/hz)*(wmac(i  ,j  ,k+1)+wmac(i,j  ,k))*(simhzx(i  ,j  ,k+1)-simhzx(i,j  ,k)) &
                     + dt2*force(i,j  ,k,comp)
             endif

             ! make sedgey by solving Riemann problem
             ! boundary conditions enforced outside of i,j,k loop
             sedgey(i,j,k,comp) = merge(sedgely(i,j,k),sedgery(i,j,k),vmac(i,j,k) .gt. ZERO)
             savg = HALF*(sedgely(i,j,k)+sedgery(i,j,k))
             sedgey(i,j,k,comp) = merge(sedgey(i,j,k,comp),savg,abs(vmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo

    ! sedgey boundary conditions
    do k=ks,ke
       do i=is,ie
          ! lo side
          if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
             if (is_vel .and. comp .eq. 2) then
                sedgey(i,js,k,comp) = ZERO
             elseif (is_vel .and. comp .ne. 2) then
                sedgey(i,js,k,comp) = merge(ZERO,sedgery(i,js,k),phys_bc(2,1).eq.NO_SLIP_WALL)
             else 
                sedgey(i,js,k,comp) = sedgery(i,js,k)
             endif
          elseif (phys_bc(2,1) .eq. INLET) then
             sedgey(i,js,k,comp) = s(i,js-1,k,comp)
          elseif (phys_bc(2,1) .eq. OUTLET) then
             if (is_vel .and. comp.eq.2) then
                sedgey(i,js,k,comp) = MIN(sedgery(i,js,k),ZERO)
             else
                sedgey(i,js,k,comp) = sedgery(i,js,k)
             end if
          endif

          ! hi side
          if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
             if (is_vel .and. comp .eq. 2) then
                sedgey(i,je+1,k,comp) = ZERO
             elseif (is_vel .and. comp .ne. 2) then
                sedgey(i,je+1,k,comp) = merge(ZERO,sedgely(i,je+1,k),phys_bc(2,2).eq.NO_SLIP_WALL)
             else 
                sedgey(i,je+1,k,comp) = sedgely(i,je+1,k)
             endif
          elseif (phys_bc(2,2) .eq. INLET) then
             sedgey(i,je+1,k,comp) = s(i,je+1,k,comp)
          elseif (phys_bc(2,2) .eq. OUTLET) then
             if (is_vel .and. comp.eq.2) then
                sedgey(i,je+1,k,comp) = MAX(sedgely(i,je+1,k),ZERO)
             else
                sedgey(i,je+1,k,comp) = sedgely(i,je+1,k)
             end if
          endif
       enddo
    enddo

    ! loop over appropriate z-faces
    do k=ks,ke+1
       do j=js,je
          do i=is,ie
             ! make sedgelz, sedgerz
             if(is_conservative) then
                sedgelz(i,j,k) = slz(i,j,k) &
                     - (dt2/hx)*(simhxy(i+1,j  ,k-1)*umac(i+1,j  ,k-1) - simhxy(i,j,k-1)*umac(i,j,k-1)) &
                     - (dt2/hy)*(simhyx(i  ,j+1,k-1)*vmac(i  ,j+1,k-1) - simhyx(i,j,k-1)*vmac(i,j,k-1)) &
                     - (dt2/hz)*s(i,j,k-1,comp)*(wmac(i,j,k  )-wmac(i,j,k-1)) &
                     + dt2*force(i,j,k-1,comp)
                sedgerz(i,j,k) = srz(i,j,k) &
                     - (dt2/hx)*(simhxy(i+1,j  ,k  )*umac(i+1,j  ,k  ) - simhxy(i,j,k  )*umac(i,j,k  )) &
                     - (dt2/hy)*(simhyx(i  ,j+1,k  )*vmac(i  ,j+1,k  ) - simhyx(i,j,k  )*vmac(i,j,k  )) &
                     - (dt2/hz)*s(i,j,k  ,comp)*(wmac(i,j,k+1)-wmac(i,j,k  )) &
                     + dt2*force(i,j,k  ,comp)
             else
                sedgelz(i,j,k) = slz(i,j,k) &
                     - (dt4/hx)*(umac(i+1,j  ,k-1)+umac(i,j,k-1))*(simhxy(i+1,j  ,k-1)-simhxy(i,j,k-1)) &
                     - (dt4/hy)*(vmac(i  ,j+1,k-1)+vmac(i,j,k-1))*(simhyx(i  ,j+1,k-1)-simhyx(i,j,k-1)) &
                     + dt2*force(i,j,k-1,comp)
                sedgerz(i,j,k) = srz(i,j,k) &
                     - (dt4/hx)*(umac(i+1,j  ,k  )+umac(i,j,k  ))*(simhxy(i+1,j  ,k  )-simhxy(i,j,k  )) &
                     - (dt4/hy)*(vmac(i  ,j+1,k  )+vmac(i,j,k  ))*(simhyx(i  ,j+1,k  )-simhyx(i,j,k  )) &
                     + dt2*force(i,j,k  ,comp)
             endif

             ! make sedgez by solving Riemann problem
             ! boundary conditions enforced outside of i,j,k loop
             sedgez(i,j,k,comp) = merge(sedgelz(i,j,k),sedgerz(i,j,k),wmac(i,j,k) .gt. ZERO)
             savg = HALF*(sedgelz(i,j,k)+sedgerz(i,j,k))
             sedgez(i,j,k,comp) = merge(sedgez(i,j,k,comp),savg,abs(wmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo

    ! sedgez boundary conditions
    do j=js,je
       do i=is,ie
          ! lo side
          if (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
             if (is_vel .and. comp .eq. 2) then
                sedgez(i,j,ks,comp) = ZERO
             elseif (is_vel .and. comp .ne. 2) then
                sedgez(i,j,ks,comp) = merge(ZERO,sedgerz(i,j,ks),phys_bc(3,1).eq.NO_SLIP_WALL)
             else 
                sedgez(i,j,ks,comp) = sedgerz(i,j,ks)
             endif
          elseif (phys_bc(3,1) .eq. INLET) then
             sedgez(i,j,ks,comp) = s(i,j,ks-1,comp)
          elseif (phys_bc(3,1) .eq. OUTLET) then
             if (is_vel .and. comp.eq.3) then
                sedgez(i,j,ks,comp) = MIN(sedgerz(i,j,ks),ZERO)
             else
                sedgez(i,j,ks,comp) = sedgerz(i,j,ks)
             end if
          endif

          ! hi side
          if (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
             if (is_vel .and. comp .eq. 2) then
                sedgez(i,j,ke+1,comp) = ZERO
             elseif (is_vel .and. comp .ne. 2) then
                sedgez(i,j,ke+1,comp) = merge(ZERO,sedgelz(i,j,ke+1),phys_bc(3,2).eq.NO_SLIP_WALL)
             else 
                sedgez(i,j,ke+1,comp) = sedgelz(i,j,ke+1)
             endif
          elseif (phys_bc(3,2) .eq. INLET) then
             sedgez(i,j,ke+1,comp) = s(i,j,ke+1,comp)
          elseif (phys_bc(3,2) .eq. OUTLET) then
             if (is_vel .and. comp.eq.3) then
                sedgez(i,j,ke+1,comp) = MAX(sedgelz(i,j,ke+1),ZERO)
             else
                sedgez(i,j,ke+1,comp) = sedgelz(i,j,ke+1)
             end if
          endif
       enddo
    enddo

    deallocate(slx,srx,simhx,sly,sry,simhy,slz,srz,simhz)
    deallocate(slxy,srxy,simhxy,slxz,srxz,simhxz,slyx,sryx,simhyx)
    deallocate(slyz,sryz,simhyz,slzx,srzx,simhzx,slzy,srzy,simhzy)
    deallocate(sedgelx,sedgerx,sedgely,sedgery,sedgelz,sedgerz)    

  end subroutine make_edge_scal_3d

end module make_edge_scal_module
 
