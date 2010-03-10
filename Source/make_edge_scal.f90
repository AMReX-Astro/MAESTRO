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

    real(kind=dp_t) :: gradw0_rad(1,0:nr_fine-1)
    type(multifab)  :: gradw0_cart(nlevs)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_edge_scal")

    do n=1,nlevs
       call multifab_build(gradw0_cart(n),s(n)%la,1,1)
    end do

    ng_s = s(1)%ng
    ng_se = sedge(1,1)%ng
    ng_um = umac(1,1)%ng
    ng_f = force(1)%ng
    ng_w0 = w0mac(1,1)%ng
    ng_n = normal(1)%ng
    ng_gw = gradw0_cart(1)%ng

    ! make a Cartesian version of dw0/dr
    if (spherical .eq. 1 .and. is_vel) then
       do r=0,nr_fine-1
          gradw0_rad(1,r) = (w0(1,r+1) - w0(1,r)) / dr(1)
       enddo
       call put_1d_array_on_cart(gradw0_rad,gradw0_cart,foextrap_comp,.false.,.false.,dx, &
                                 the_bc_level,mla)
    else
       do n=1,nlevs
          call setval(gradw0_cart(n),ZERO,all=.true.)
       end do
    endif

    do n=1,nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n),i) ) cycle
          sop  => dataptr(s(n),i)
          sepx => dataptr(sedge(n,1),i)
          ump  => dataptr(umac(n,1),i)
          fp   => dataptr(force(n),i)
          lo   =  lwb(get_box(s(n),i))
          hi   =  upb(get_box(s(n),i))
          select case (dm)
          case (1)
             do scomp = start_scomp, start_scomp + num_comp - 1
                bccomp = start_bccomp + scomp - start_scomp
                call make_edge_scal_1d(n, sop(:,1,1,:), ng_s, &
                                       sepx(:,1,1,:), ng_se, &
                                        ump(:,1,1,1), ng_um, &
                                       fp(:,1,1,:), ng_f, w0(n,:), &
                                       lo, hi, dx(n,:), dt, is_vel, &
                                       the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                       the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp:), &
                                       scomp, is_conservative)
             end do

          case (2)
             vmp  => dataptr(umac(n,2),i)
             sepy => dataptr(sedge(n,2),i)
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
             vmp  => dataptr(umac(n,2),i)
             wmp  => dataptr(  umac(n,3),i)
             sepy => dataptr(sedge(n,2),i)
             sepz => dataptr( sedge(n,3),i)
             w0xp => dataptr(w0mac(n,1),i)
             w0yp => dataptr(w0mac(n,2),i)
             w0zp => dataptr(w0mac(n,3),i)
             np   => dataptr(normal(n), i)
             gw0p => dataptr(gradw0_cart(n), i)
             do scomp = start_scomp, start_scomp + num_comp - 1
               bccomp = start_bccomp + scomp - start_scomp
               if (spherical .eq. 1) then
                  call make_edge_scal_3d(n, sop(:,:,:,:), ng_s, &
                                         sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                         ng_se, ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                         ng_um, fp(:,:,:,:), ng_f, np(:,:,:,:), ng_n, &
                                         w0(1,:), w0xp(:,:,:,1), w0yp(:,:,:,1), &
                                         w0zp(:,:,:,1), ng_w0, gw0p(:,:,:,1), ng_gw, lo, &
                                         hi, dx(n,:), dt, is_vel, &
                                         the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                         the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp:), &
                                         scomp, is_conservative)
               else
                  call make_edge_scal_3d(n, sop(:,:,:,:), ng_s, &
                                         sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                         ng_se, ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                         ng_um, fp(:,:,:,:), ng_f, np(:,:,:,:), ng_n, &
                                         w0(n,:), w0xp(:,:,:,1), w0yp(:,:,:,1), &
                                         w0zp(:,:,:,1), ng_w0, gw0p(:,:,:,1), ng_gw, lo, &
                                         hi, dx(n,:), dt, is_vel, &
                                         the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                         the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp:), &
                                         scomp, is_conservative)
               end if
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

    do n=1,nlevs
       call destroy(gradw0_cart(n))
    end do

    call destroy(bpt)
    
  end subroutine make_edge_scal
  
  subroutine make_edge_scal_1d(n,s,ng_s,sedgex,ng_se,umac,ng_um, &
                               force,ng_f,w0,lo,hi,dx,dt,is_vel,phys_bc,adv_bc, &
                               comp,is_conservative)

    use geometry, only: nr
    use bc_module
    use slope_module
    use bl_constants_module
    use variables, only: rel_eps
    use ppm_module
    use probin_module, only: ppm_type

    integer        , intent(in   ) :: lo(:),hi(:),n,ng_s,ng_se,ng_um,ng_f
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s :,:)
    real(kind=dp_t), intent(inout) :: sedgex(lo(1)-ng_se:,:)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)-ng_um:)
    real(kind=dp_t), intent(in   ) ::  force(lo(1)-ng_f :,:)
    real(kind=dp_t), intent(in   ) ::     w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:),dt
    logical        , intent(in   ) :: is_vel
    integer        , intent(in   ) :: phys_bc(:,:)
    integer        , intent(in   ) :: adv_bc(:,:,:)
    integer        , intent(in   ) :: comp
    logical        , intent(in   ) :: is_conservative

    ! Local variables
    real(kind=dp_t) :: slopex(lo(1)-1:hi(1)+1,1)

    real(kind=dp_t) :: hx,dt2,dt4,savg

    integer :: i,is,ie

    real(kind=dp_t), allocatable :: Ip(:)
    real(kind=dp_t), allocatable :: Im(:)

    ! these correspond to s_L^x, etc.
    real(kind=dp_t), allocatable:: slx(:),srx(:)

    ! these correspond to s_{\i-\half\e_x}^x, etc.
    real(kind=dp_t), allocatable:: simhx(:)

    ! these correspond to \mathrm{sedge}_L^x, etc.
    real(kind=dp_t), allocatable:: sedgelx(:),sedgerx(:)

    allocate(Ip(lo(1)-1:hi(1)+1))
    allocate(Im(lo(1)-1:hi(1)+1))

    ! Normal predictor states.
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse direction
    allocate(slx  (lo(1):hi(1)+1))
    allocate(srx  (lo(1):hi(1)+1))
    allocate(simhx(lo(1):hi(1)+1))

    ! Final edge states.
    ! lo:hi+1 in the normal direction
    ! lo:hi in the transverse direction
    allocate(sedgelx(lo(1):hi(1)+1))
    allocate(sedgerx(lo(1):hi(1)+1))

    is = lo(1)
    ie = hi(1)

    if (ppm_type .gt. 0) then
       call ppm_fpu_1d(n,s(:,comp),ng_s,umac,ng_um,Ip,Im,w0,lo,hi,adv_bc(:,:,1),dx,dt)
    else
       call slopex_1d(s(:,comp:),slopex,lo,hi,ng_s,1,adv_bc)
    end if

    dt2 = HALF*dt
    dt4 = dt/4.0d0

    hx = dx(1)
    
    !******************************************************************
    ! Create s_{\i-\half\e_x}^x, etc.
    !******************************************************************

    ! loop over appropriate x-faces    
    if (ppm_type .gt. 0) then
       do i=is,ie+1
          ! make slx, srx with 1D extrapolation
          slx(i) = Ip(i-1)
          srx(i) = Im(i  )
       end do
    else
       do i=is,ie+1
          ! make slx, srx with 1D extrapolation
          slx(i) = s(i-1,comp) + (HALF - dt2*umac(i)/hx)*slopex(i-1,1)
          srx(i) = s(i  ,comp) - (HALF + dt2*umac(i)/hx)*slopex(i  ,1)
       enddo
    end if

    ! impose lo side bc's
    if (phys_bc(1,1) .eq. INLET) then
       slx(is) = s(is-1,comp)
       srx(is) = s(is-1,comp)
    else if (phys_bc(1,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 1) then
          slx(is) = ZERO
          srx(is) = ZERO
       else
          slx(is) = srx(is)
       end if
    else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slx(is) = ZERO
          srx(is) = ZERO
       else
          slx(is) = srx(is)
       end if
    else if (phys_bc(1,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 1) then
          slx(is) = min(srx(is),ZERO)
          srx(is) = min(srx(is),ZERO)
       else
          slx(is) = srx(is)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(1,2) .eq. INLET) then
       slx(ie+1) = s(ie+1,comp)
       srx(ie+1) = s(ie+1,comp)
    else if (phys_bc(1,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 1) then
          slx(ie+1) = ZERO
          srx(ie+1) = ZERO
       else
          srx(ie+1) = slx(ie+1)
       end if
    else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slx(ie+1) = ZERO
          srx(ie+1) = ZERO
       else
          srx(ie+1) = slx(ie+1)
       end if
    else if (phys_bc(1,2) .eq. OUTLET) then       
       if (is_vel .and. comp .eq. 1) then
          slx(ie+1) = max(slx(ie+1),ZERO)
          srx(ie+1) = max(slx(ie+1),ZERO)
       else
          srx(ie+1) = slx(ie+1)
       end if
    end if

    do i=is,ie+1
       ! make simhx by solving Riemann problem
       simhx(i) = merge(slx(i),srx(i),umac(i) .gt. ZERO)
       savg = HALF*(slx(i)+srx(i))
       simhx(i) = merge(simhx(i),savg,abs(umac(i)) .gt. rel_eps)
    enddo

    !******************************************************************
    ! Create sedgelx, etc.
    !******************************************************************

    ! loop over appropriate x-faces
    do i=is,ie+1
       ! make sedgelx, sedgerx
       if(is_conservative) then
          sedgelx(i) = slx(i) &
               - (dt2/hx)*s(i-1,comp)*(umac(i  )-umac(i-1)) &
               + dt2*force(i-1,comp)
          sedgerx(i) = srx(i) &
               - (dt2/hx)*s(i  ,comp)*(umac(i+1)-umac(i  )) &
               + dt2*force(i  ,comp)
       else
          sedgelx(i) = slx(i) + dt2*force(i-1,comp)
          sedgerx(i) = srx(i) + dt2*force(i  ,comp)
       end if

       ! add the (Utilde . e_r) d w_0 /dr e_r term here
       ! umac contains w0 so we need to subtract it off

       if (is_vel .and. comp .eq. 1) then
          if (i .eq. 0) then
             ! sedgelx unchanged since dw_0 / dr = 0
             sedgerx(i) = sedgerx(i) &
                  -(dt4/hx)*(umac(i+1)-w0(i+1)+umac(i)-w0(i))*(w0(i+1)-w0(i))
          else if (i .eq. nr(n)) then
             sedgelx(i) = sedgelx(i) &
                  -(dt4/hx)*(umac(i)-w0(i)+umac(i-1)-w0(i-1))*(w0(i)-w0(i-1))
             ! sedgerx unchanged since dw_0 / dr = 0
          else
             sedgelx(i) = sedgelx(i) &
                  -(dt4/hx)*(umac(i)-w0(i)+umac(i-1)-w0(i-1))*(w0(i)-w0(i-1))
             sedgerx(i) = sedgerx(i) &
                  -(dt4/hx)*(umac(i+1)-w0(i+1)+umac(i)-w0(i))*(w0(i+1)-w0(i))
          end if
       end if

       ! make sedgex by solving Riemann problem
       ! boundary conditions enforced outside of i loop
       sedgex(i,comp) = merge(sedgelx(i),sedgerx(i),umac(i) .gt. ZERO)
       savg = HALF*(sedgelx(i)+sedgerx(i))
       sedgex(i,comp) = merge(sedgex(i,comp),savg,abs(umac(i)) .gt. rel_eps)
    enddo
 
    ! impose lo side bc's
    if (phys_bc(1,1) .eq. INLET) then
       sedgex(is,comp) = s(is-1,comp)
    else if (phys_bc(1,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 1) then
          sedgex(is,comp) = ZERO
       else
          sedgex(is,comp) = sedgerx(is)
       end if
    else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          sedgex(is,comp) = ZERO
       else
          sedgex(is,comp) = sedgerx(is)
       end if
    else if (phys_bc(1,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 1) then
          sedgex(is,comp) = min(sedgerx(is),ZERO)
       else
          sedgex(is,comp) = sedgerx(is)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(1,2) .eq. INLET) then
       sedgex(ie+1,comp) = s(ie+1,comp)
    else if (phys_bc(1,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 1) then
          sedgex(ie+1,comp) = ZERO
       else
          sedgex(ie+1,comp) = sedgelx(ie+1)
       end if
    else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          sedgex(ie+1,comp) = ZERO
       else
          sedgex(ie+1,comp) = sedgelx(ie+1)
       end if
    else if (phys_bc(1,2) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 1) then
          sedgex(ie+1,comp) = max(sedgelx(ie+1),ZERO)
       else
          sedgex(ie+1,comp) = sedgelx(ie+1)
       end if
    end if

    deallocate(slx,srx,simhx,sedgelx,sedgerx)
    deallocate(Ip,Im)

  end subroutine make_edge_scal_1d
  
  subroutine make_edge_scal_2d(n,s,ng_s,sedgex,sedgey,ng_se,umac,vmac,ng_um, &
                               force,ng_f,w0,lo,hi,dx,dt,is_vel,phys_bc,adv_bc, &
                               comp,is_conservative)

    use geometry, only: nr
    use bc_module
    use slope_module
    use bl_constants_module
    use variables, only: rel_eps
    use ppm_module
    use probin_module, only: ppm_type

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

    real(kind=dp_t) :: hx,hy,dt2,dt4,savg

    integer :: i,j,is,js,ie,je

    real(kind=dp_t), allocatable :: Ip(:,:,:)
    real(kind=dp_t), allocatable :: Im(:,:,:)

    ! these correspond to s_L^x, etc.
    real(kind=dp_t), allocatable:: slx(:,:),srx(:,:)
    real(kind=dp_t), allocatable:: sly(:,:),sry(:,:)

    ! these correspond to s_{\i-\half\e_x}^x, etc.
    real(kind=dp_t), allocatable:: simhx(:,:),simhy(:,:)

    ! these correspond to \mathrm{sedge}_L^x, etc.
    real(kind=dp_t), allocatable:: sedgelx(:,:),sedgerx(:,:)
    real(kind=dp_t), allocatable:: sedgely(:,:),sedgery(:,:)

    allocate(Ip(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
    allocate(Im(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))

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

    if (ppm_type .gt. 0) then
       call ppm_fpu_2d(n,s(:,:,comp),ng_s,umac,vmac,ng_um,Ip,Im,w0,lo,hi,adv_bc(:,:,1),dx,dt)
    else
       call slopex_2d(s(:,:,comp:),slopex,lo,hi,ng_s,1,adv_bc)
       call slopey_2d(s(:,:,comp:),slopey,lo,hi,ng_s,1,adv_bc)
    end if

    dt2 = HALF*dt
    dt4 = dt/4.0d0

    hx = dx(1)
    hy = dx(2)
    
    !******************************************************************
    ! Create s_{\i-\half\e_x}^x, etc.
    !******************************************************************

    ! loop over appropriate x-faces    
    if (ppm_type .gt. 0) then
       do j=js-1,je+1
          do i=is,ie+1
             ! make slx, srx with 1D extrapolation
             slx(i,j) = Ip(i-1,j,1)
             srx(i,j) = Im(i  ,j,1)
          end do
       end do
    else
       do j=js-1,je+1
          do i=is,ie+1
             ! make slx, srx with 1D extrapolation
             slx(i,j) = s(i-1,j,comp) + (HALF - dt2*umac(i,j)/hx)*slopex(i-1,j,1)
             srx(i,j) = s(i  ,j,comp) - (HALF + dt2*umac(i,j)/hx)*slopex(i  ,j,1)
          enddo
       enddo
    end if

    ! impose lo side bc's
    if (phys_bc(1,1) .eq. INLET) then
       slx(is,js-1:je+1) = s(is-1,js-1:je+1,comp)
       srx(is,js-1:je+1) = s(is-1,js-1:je+1,comp)
    else if (phys_bc(1,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 1) then
          slx(is,js-1:je+1) = ZERO
          srx(is,js-1:je+1) = ZERO
       else
          slx(is,js-1:je+1) = srx(is,js-1:je+1)
       end if
    else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slx(is,js-1:je+1) = ZERO
          srx(is,js-1:je+1) = ZERO
       else
          slx(is,js-1:je+1) = srx(is,js-1:je+1)
       end if
    else if (phys_bc(1,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 1) then
          slx(is,js-1:je+1) = min(srx(is,js-1:je+1),ZERO)
          srx(is,js-1:je+1) = min(srx(is,js-1:je+1),ZERO)
       else
          slx(is,js-1:je+1) = srx(is,js-1:je+1)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(1,2) .eq. INLET) then
       slx(ie+1,js-1:je+1) = s(ie+1,js-1:je+1,comp)
       srx(ie+1,js-1:je+1) = s(ie+1,js-1:je+1,comp)
    else if (phys_bc(1,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 1) then
          slx(ie+1,js-1:je+1) = ZERO
          srx(ie+1,js-1:je+1) = ZERO
       else
          srx(ie+1,js-1:je+1) = slx(ie+1,js-1:je+1)
       end if
    else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slx(ie+1,js-1:je+1) = ZERO
          srx(ie+1,js-1:je+1) = ZERO
       else
          srx(ie+1,js-1:je+1) = slx(ie+1,js-1:je+1)
       end if
    else if (phys_bc(1,2) .eq. OUTLET) then       
       if (is_vel .and. comp .eq. 1) then
          slx(ie+1,js-1:je+1) = max(slx(ie+1,js-1:je+1),ZERO)
          srx(ie+1,js-1:je+1) = max(slx(ie+1,js-1:je+1),ZERO)
       else
          srx(ie+1,js-1:je+1) = slx(ie+1,js-1:je+1)
       end if
    end if

    do j=js-1,je+1
       do i=is,ie+1
          ! make simhx by solving Riemann problem
          simhx(i,j) = merge(slx(i,j),srx(i,j),umac(i,j) .gt. ZERO)
          savg = HALF*(slx(i,j)+srx(i,j))
          simhx(i,j) = merge(simhx(i,j),savg,abs(umac(i,j)) .gt. rel_eps)
       enddo
    enddo

    ! loop over appropriate y-faces
    if (ppm_type .gt. 0) then
       do j=js,je+1
          do i=is-1,ie+1
             ! make sly, sry with 1D extrapolation
             sly(i,j) = Ip(i,j-1,2)
             sry(i,j) = Im(i,j  ,2)
          enddo
       enddo
    else
       do j=js,je+1
          do i=is-1,ie+1
             ! make sly, sry with 1D extrapolation
             sly(i,j) = s(i,j-1,comp) + (HALF - dt2*vmac(i,j)/hy)*slopey(i,j-1,1)
             sry(i,j) = s(i,j  ,comp) - (HALF + dt2*vmac(i,j)/hy)*slopey(i,j  ,1)
          enddo
       enddo
    end if

    ! impose lo side bc's
    if (phys_bc(2,1) .eq. INLET) then
       sly(is-1:ie+1,js) = s(is-1:ie+1,js-1,comp)
       sry(is-1:ie+1,js) = s(is-1:ie+1,js-1,comp)
    else if (phys_bc(2,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 2) then
          sly(is-1:ie+1,js) = ZERO
          sry(is-1:ie+1,js) = ZERO
       else
          sly(is-1:ie+1,js) = sry(is-1:ie+1,js)
       end if
    else if (phys_bc(2,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          sly(is-1:ie+1,js) = ZERO
          sry(is-1:ie+1,js) = ZERO
       else
          sly(is-1:ie+1,js) = sry(is-1:ie+1,js)
       end if
    else if (phys_bc(2,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 2) then
          sly(is-1:ie+1,js) = min(sry(is-1:ie+1,js),ZERO)
          sry(is-1:ie+1,js) = min(sry(is-1:ie+1,js),ZERO)
       else
          sly(is-1:ie+1,js) = sry(is-1:ie+1,js)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(2,2) .eq. INLET) then
       sly(is-1:ie+1,je+1) = s(is-1:ie+1,je+1,comp)
       sry(is-1:ie+1,je+1) = s(is-1:ie+1,je+1,comp)
    else if (phys_bc(2,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 2) then
          sly(is-1:ie+1,je+1) = ZERO
          sry(is-1:ie+1,je+1) = ZERO
       else
          sry(is-1:ie+1,je+1) = sly(is-1:ie+1,je+1)
       end if
    else if (phys_bc(2,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          sly(is-1:ie+1,je+1) = ZERO
          sry(is-1:ie+1,je+1) = ZERO
       else
          sry(is-1:ie+1,je+1) = sly(is-1:ie+1,je+1)
       end if
    else if (phys_bc(2,2) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 2) then
          sly(is-1:ie+1,je+1) = max(sly(is-1:ie+1,je+1),ZERO)
          sry(is-1:ie+1,je+1) = max(sly(is-1:ie+1,je+1),ZERO)
       else
          sry(is-1:ie+1,je+1) = sly(is-1:ie+1,je+1)
       end if
    end if

    do j=js,je+1
       do i=is-1,ie+1
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
          end if

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
 
    ! impose lo side bc's
    if (phys_bc(1,1) .eq. INLET) then
       sedgex(is,js:je,comp) = s(is-1,js:je,comp)
    else if (phys_bc(1,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 1) then
          sedgex(is,js:je,comp) = ZERO
       else
          sedgex(is,js:je,comp) = sedgerx(is,js:je)
       end if
    else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          sedgex(is,js:je,comp) = ZERO
       else
          sedgex(is,js:je,comp) = sedgerx(is,js:je)
       end if
    else if (phys_bc(1,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 1) then
          sedgex(is,js:je,comp) = min(sedgerx(is,js:je),ZERO)
       else
          sedgex(is,js:je,comp) = sedgerx(is,js:je)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(1,2) .eq. INLET) then
       sedgex(ie+1,js:je,comp) = s(ie+1,js:je,comp)
    else if (phys_bc(1,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 1) then
          sedgex(ie+1,js:je,comp) = ZERO
       else
          sedgex(ie+1,js:je,comp) = sedgelx(ie+1,js:je)
       end if
    else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          sedgex(ie+1,js:je,comp) = ZERO
       else
          sedgex(ie+1,js:je,comp) = sedgelx(ie+1,js:je)
       end if
    else if (phys_bc(1,2) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 1) then
          sedgex(ie+1,js:je,comp) = max(sedgelx(ie+1,js:je),ZERO)
       else
          sedgex(ie+1,js:je,comp) = sedgelx(ie+1,js:je)
       end if
    end if

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
          end if

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

    ! impose lo side bc's
    if (phys_bc(2,1) .eq. INLET) then
       sedgey(is:ie,js,comp) = s(is:ie,js-1,comp)
    else if (phys_bc(2,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 2) then
          sedgey(is:ie,js,comp) = ZERO
       else
          sedgey(is:ie,js,comp) = sedgery(is:ie,js)
       end if
    else if (phys_bc(2,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          sedgey(is:ie,js,comp) = ZERO
       else
          sedgey(is:ie,js,comp) = sedgery(is:ie,js)
       end if
    else if (phys_bc(2,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 2) then
          sedgey(is:ie,js,comp) = min(sedgery(is:ie,js),ZERO)
       else
          sedgey(is:ie,js,comp) = sedgery(is:ie,js)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(2,2) .eq. INLET) then
       sedgey(is:ie,je+1,comp) = s(is:ie,je+1,comp)
    else if (phys_bc(2,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 2) then
          sedgey(is:ie,je+1,comp) = ZERO
       else
          sedgey(is:ie,je+1,comp) = sedgely(is:ie,je+1)
       end if
    else if (phys_bc(2,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          sedgey(is:ie,je+1,comp) = ZERO
       else
          sedgey(is:ie,je+1,comp) = sedgely(is:ie,je+1)
       end if
    else if (phys_bc(2,2) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 2) then
          sedgey(is:ie,je+1,comp) = max(sedgely(is:ie,je+1),ZERO)
       else
          sedgey(is:ie,je+1,comp) = sedgely(is:ie,je+1)
       end if
    end if

    deallocate(slx,srx,sly,sry,simhx,simhy,sedgelx,sedgerx,sedgely,sedgery)
    deallocate(Ip,Im)

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
    use ppm_module
    use probin_module, only: ppm_type

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
    real(kind=dp_t), allocatable :: slopex(:,:,:,:)
    real(kind=dp_t), allocatable :: slopey(:,:,:,:)
    real(kind=dp_t), allocatable :: slopez(:,:,:,:)

    real(kind=dp_t) :: hx,hy,hz,dt2,dt3,dt4,dt6
    real(kind=dp_t) :: savg
    real(kind=dp_t) :: Ut_dot_er

    integer :: i,j,k,is,js,ks,ie,je,ke

    real(kind=dp_t), allocatable :: Ip(:,:,:,:)
    real(kind=dp_t), allocatable :: Im(:,:,:,:)

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

    allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))
    allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))
    allocate(slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))

    allocate(Ip(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
    allocate(Im(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)
    ks = lo(3)
    ke = hi(3)

    if (ppm_type .gt. 0) then
       call ppm_fpu_3d(n,s(:,:,:,comp),ng_s,umac,vmac,wmac,ng_um,Ip,Im,w0, &
                       w0macx,w0macy,w0macz,ng_w0,lo,hi,adv_bc(:,:,1),dx,dt)
    else
       do k = lo(3)-1,hi(3)+1
          call slopex_2d(s(:,:,k,comp:),slopex(:,:,k,:),lo,hi,ng_s,1,adv_bc)
          call slopey_2d(s(:,:,k,comp:),slopey(:,:,k,:),lo,hi,ng_s,1,adv_bc)
       end do
       call slopez_3d(s(:,:,:,comp:),slopez,lo,hi,ng_s,1,adv_bc)
    end if

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

    ! Normal predictor states.
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    allocate(slx  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(srx  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(simhx(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
    
    ! loop over appropriate x-faces
    if (ppm_type .gt. 0) then
!$omp parallel do private(i,j,k)
       do k=ks-1,ke+1
          do j=js-1,je+1
             do i=is,ie+1
                ! make slx, srx with 1D extrapolation
                slx(i,j,k) = Ip(i-1,j,k,1)
                srx(i,j,k) = Im(i  ,j,k,1)
             end do
          end do
       end do
!$omp end parallel do
    else
!$omp parallel do private(i,j,k)
       do k=ks-1,ke+1
          do j=js-1,je+1
             do i=is,ie+1
                ! make slx, srx with 1D extrapolation
                slx(i,j,k) = s(i-1,j,k,comp) + (HALF - dt2*umac(i,j,k)/hx)*slopex(i-1,j,k,1)
                srx(i,j,k) = s(i  ,j,k,comp) - (HALF + dt2*umac(i,j,k)/hx)*slopex(i  ,j,k,1)
             enddo
          enddo
       enddo
!$omp end parallel do
    end if

    ! impose lo side bc's
    if (phys_bc(1,1) .eq. INLET) then
       slx(is,js-1:je+1,ks-1:ke+1) = s(is-1,js-1:je+1,ks-1:ke+1,comp)
       srx(is,js-1:je+1,ks-1:ke+1) = s(is-1,js-1:je+1,ks-1:ke+1,comp)
    else if (phys_bc(1,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 1) then
          slx(is,js-1:je+1,ks-1:ke+1) = ZERO
          srx(is,js-1:je+1,ks-1:ke+1) = ZERO
       else
          slx(is,js-1:je+1,ks-1:ke+1) = srx(is,js-1:je+1,ks-1:ke+1)
       end if
    else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slx(is,js-1:je+1,ks-1:ke+1) = ZERO
          srx(is,js-1:je+1,ks-1:ke+1) = ZERO
       else
          slx(is,js-1:je+1,ks-1:ke+1) = srx(is,js-1:je+1,ks-1:ke+1)
       end if
    else if (phys_bc(1,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 1) then
          slx(is,js-1:je+1,ks-1:ke+1) = min(srx(is,js-1:je+1,ks-1:ke+1),ZERO)
          srx(is,js-1:je+1,ks-1:ke+1) = min(srx(is,js-1:je+1,ks-1:ke+1),ZERO)
       else
          slx(is,js-1:je+1,ks-1:ke+1) = srx(is,js-1:je+1,ks-1:ke+1)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(1,2) .eq. INLET) then
       slx(ie+1,js-1:je+1,ks-1:ke+1) = s(ie+1,js-1:je+1,ks-1:ke+1,comp)
       srx(ie+1,js-1:je+1,ks-1:ke+1) = s(ie+1,js-1:je+1,ks-1:ke+1,comp)
    else if (phys_bc(1,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 1) then
          slx(ie+1,js-1:je+1,ks-1:ke+1) = ZERO
          srx(ie+1,js-1:je+1,ks-1:ke+1) = ZERO
       else
          srx(ie+1,js-1:je+1,ks-1:ke+1) = slx(ie+1,js-1:je+1,ks-1:ke+1)
       end if
    else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slx(ie+1,js-1:je+1,ks-1:ke+1) = ZERO
          srx(ie+1,js-1:je+1,ks-1:ke+1) = ZERO
       else
          srx(ie+1,js-1:je+1,ks-1:ke+1) = slx(ie+1,js-1:je+1,ks-1:ke+1)
       end if
    else if (phys_bc(1,2) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 1) then
          slx(ie+1,js-1:je+1,ks-1:ke+1) = max(slx(ie+1,js-1:je+1,ks-1:ke+1),ZERO)
          srx(ie+1,js-1:je+1,ks-1:ke+1) = max(slx(ie+1,js-1:je+1,ks-1:ke+1),ZERO)
       else
          srx(ie+1,js-1:je+1,ks-1:ke+1) = slx(ie+1,js-1:je+1,ks-1:ke+1)
       end if
    end if

!$omp parallel do private(i,j,k,savg)
    do k=ks-1,ke+1
       do j=js-1,je+1
          do i=is,ie+1
             ! make simhx by solving Riemann problem
             simhx(i,j,k) = merge(slx(i,j,k),srx(i,j,k),umac(i,j,k) .gt. ZERO)
             savg = HALF*(slx(i,j,k)+srx(i,j,k))
             simhx(i,j,k) = merge(simhx(i,j,k),savg,abs(umac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
!$omp end parallel do

    ! Normal predictor states.
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    allocate(sly  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(sry  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(simhy(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3)-1:hi(3)+1))

    ! loop over appropriate y-faces
    if (ppm_type .gt. 0) then
!$omp parallel do private(i,j,k)
       do k=ks-1,ke+1
          do j=js,je+1
             do i=is-1,ie+1
                ! make sly, sry with 1D extrapolation
                sly(i,j,k) = Ip(i,j-1,k,2)
                sry(i,j,k) = Im(i,j  ,k,2)
             enddo
          enddo
       enddo
!$omp end parallel do
    else
!$omp parallel do private(i,j,k)
       do k=ks-1,ke+1
          do j=js,je+1
             do i=is-1,ie+1
                ! make sly, sry with 1D extrapolation
                sly(i,j,k) = s(i,j-1,k,comp) + (HALF - dt2*vmac(i,j,k)/hy)*slopey(i,j-1,k,1)
                sry(i,j,k) = s(i,j  ,k,comp) - (HALF + dt2*vmac(i,j,k)/hy)*slopey(i,j  ,k,1)
             enddo
          enddo
       enddo
!$omp end parallel do
    end if

    ! impose lo side bc's
    if (phys_bc(2,1) .eq. INLET) then
       sly(is-1:ie+1,js,ks-1:ke+1) = s(is-1:ie+1,js-1,ks-1:ke+1,comp)
       sry(is-1:ie+1,js,ks-1:ke+1) = s(is-1:ie+1,js-1,ks-1:ke+1,comp)
    else if (phys_bc(2,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 2) then
          sly(is-1:ie+1,js,ks-1:ke+1) = ZERO
          sry(is-1:ie+1,js,ks-1:ke+1) = ZERO
       else
          sly(is-1:ie+1,js,ks-1:ke+1) = sry(is-1:ie+1,js,ks-1:ke+1)
       end if
    else if (phys_bc(2,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          sly(is-1:ie+1,js,ks-1:ke+1) = ZERO
          sry(is-1:ie+1,js,ks-1:ke+1) = ZERO
       else
          sly(is-1:ie+1,js,ks-1:ke+1) = sry(is-1:ie+1,js,ks-1:ke+1)
       end if
    else if (phys_bc(2,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 2) then
          sly(is-1:ie+1,js,ks-1:ke+1) = min(sry(is-1:ie+1,js,ks-1:ke+1),ZERO)
          sry(is-1:ie+1,js,ks-1:ke+1) = min(sry(is-1:ie+1,js,ks-1:ke+1),ZERO)
       else
          sly(is-1:ie+1,js,ks-1:ke+1) = sry(is-1:ie+1,js,ks-1:ke+1)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(2,2) .eq. INLET) then
       sly(is-1:ie+1,je+1,ks-1:ke+1) = s(is-1:ie+1,je+1,ks-1:ke+1,comp)
       sry(is-1:ie+1,je+1,ks-1:ke+1) = s(is-1:ie+1,je+1,ks-1:ke+1,comp)
    else if (phys_bc(2,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 2) then
          sly(is-1:ie+1,je+1,ks-1:ke+1) = ZERO
          sry(is-1:ie+1,je+1,ks-1:ke+1) = ZERO
       else
          sry(is-1:ie+1,je+1,ks-1:ke+1) = sly(is-1:ie+1,je+1,ks-1:ke+1)
       end if
    else if (phys_bc(2,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          sly(is-1:ie+1,je+1,ks-1:ke+1) = ZERO
          sry(is-1:ie+1,je+1,ks-1:ke+1) = ZERO
       else
          sry(is-1:ie+1,je+1,ks-1:ke+1) = sly(is-1:ie+1,je+1,ks-1:ke+1)
       end if
    else if (phys_bc(2,2) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 2) then
          sly(is-1:ie+1,je+1,ks-1:ke+1) = max(sly(is-1:ie+1,je+1,ks-1:ke+1),ZERO)
          sry(is-1:ie+1,je+1,ks-1:ke+1) = max(sly(is-1:ie+1,je+1,ks-1:ke+1),ZERO)
       else
          sry(is-1:ie+1,je+1,ks-1:ke+1) = sly(is-1:ie+1,je+1,ks-1:ke+1)
       end if
    end if

!$omp parallel do private(i,j,k,savg)
    do k=ks-1,ke+1
       do j=js,je+1
          do i=is-1,ie+1
             ! make simhy by solving Riemann problem
             simhy(i,j,k) = merge(sly(i,j,k),sry(i,j,k),vmac(i,j,k) .gt. ZERO)
             savg = HALF*(sly(i,j,k)+sry(i,j,k))
             simhy(i,j,k) = merge(simhy(i,j,k),savg,abs(vmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
!$omp end parallel do

    ! Normal predictor states.
    ! Allocated from lo:hi+1 in the normal direction
    ! lo-1:hi+1 in the transverse directions
    allocate(slz  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(srz  (lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(simhz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)+1))

    ! loop over appropriate z-faces
    if (ppm_type .gt. 0) then
!$omp parallel do private(i,j,k)
       do k=ks,ke+1
          do j=js-1,je+1
             do i=is-1,ie+1
                ! make slz, srz with 1D extrapolation
                slz(i,j,k) = Ip(i,j,k-1,3)
                srz(i,j,k) = Im(i,j,k  ,3)
             enddo
          enddo
       enddo
!$omp end parallel do
    else
!$omp parallel do private(i,j,k)
       do k=ks,ke+1
          do j=js-1,je+1
             do i=is-1,ie+1
                ! make slz, srz with 1D extrapolation
                slz(i,j,k) = s(i,j,k-1,comp) + (HALF - dt2*wmac(i,j,k)/hz)*slopez(i,j,k-1,1)
                srz(i,j,k) = s(i,j,k  ,comp) - (HALF + dt2*wmac(i,j,k)/hz)*slopez(i,j,k  ,1)
             enddo
          enddo
       enddo
!$omp end parallel do
    end if

    deallocate(slopex,slopey,slopez,Ip,Im)

    ! impose lo side bc's
    if (phys_bc(3,1) .eq. INLET) then
       slz(is-1:ie+1,js-1:je+1,ks) = s(is-1:ie+1,js-1:je+1,ks,comp)
       srz(is-1:ie+1,js-1:je+1,ks) = s(is-1:ie+1,js-1:je+1,ks,comp)
    else if (phys_bc(3,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 3) then
          slz(is-1:ie+1,js-1:je+1,ks) = ZERO
          srz(is-1:ie+1,js-1:je+1,ks) = ZERO
       else
          slz(is-1:ie+1,js-1:je+1,ks) = srz(is-1:ie+1,js-1:je+1,ks)
       end if
    else if (phys_bc(3,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slz(is-1:ie+1,js-1:je+1,ks) = ZERO
          srz(is-1:ie+1,js-1:je+1,ks) = ZERO
       else
          slz(is-1:ie+1,js-1:je+1,ks) = srz(is-1:ie+1,js-1:je+1,ks)
       end if
    else if (phys_bc(3,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 3) then
          slz(is-1:ie+1,js-1:je+1,ks) = min(srz(is-1:ie+1,js-1:je+1,ks),ZERO)
          srz(is-1:ie+1,js-1:je+1,ks) = min(srz(is-1:ie+1,js-1:je+1,ks),ZERO)
       else
          slz(is-1:ie+1,js-1:je+1,ks) = srz(is-1:ie+1,js-1:je+1,ks)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(3,2) .eq. INLET) then
       slz(is-1:ie+1,js-1:je+1,ke+1) = s(is-1:ie+1,js-1:je+1,ke+1,comp)
       srz(is-1:ie+1,js-1:je+1,ke+1) = s(is-1:ie+1,js-1:je+1,ke+1,comp)
    else if (phys_bc(3,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 3) then
          slz(is-1:ie+1,js-1:je+1,ke+1) = ZERO
          srz(is-1:ie+1,js-1:je+1,ke+1) = ZERO
       else
          srz(is-1:ie+1,js-1:je+1,ke+1) = slz(is-1:ie+1,js-1:je+1,ke+1)
       end if
    else if (phys_bc(3,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slz(is-1:ie+1,js-1:je+1,ke+1) = ZERO
          srz(is-1:ie+1,js-1:je+1,ke+1) = ZERO
       else
          srz(is-1:ie+1,js-1:je+1,ke+1) = slz(is-1:ie+1,js-1:je+1,ke+1)
       end if
    else if (phys_bc(3,2) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 3) then
          slz(is-1:ie+1,js-1:je+1,ke+1) = max(slz(is-1:ie+1,js-1:je+1,ke+1),ZERO)
          srz(is-1:ie+1,js-1:je+1,ke+1) = max(slz(is-1:ie+1,js-1:je+1,ke+1),ZERO)
       else
          srz(is-1:ie+1,js-1:je+1,ke+1) = slz(is-1:ie+1,js-1:je+1,ke+1)
       end if
    end if

!$omp parallel do private(i,j,k,savg)
    do k=ks,ke+1
       do j=js-1,je+1
          do i=is-1,ie+1
             ! make simhz by solving Riemann problem
             simhz(i,j,k) = merge(slz(i,j,k),srz(i,j,k),wmac(i,j,k) .gt. ZERO)
             savg = HALF*(slz(i,j,k)+srz(i,j,k))
             simhz(i,j,k) = merge(simhz(i,j,k),savg,abs(wmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
!$omp end parallel do

    !******************************************************************
    ! Create s_{\i-\half\e_x}^{x|y}, etc.
    !******************************************************************

    allocate(slxy  (lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))
    allocate(srxy  (lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))
    allocate(simhxy(lo(1):hi(1)+1,lo(2):hi(2),lo(3)-1:hi(3)+1))

    ! loop over appropriate xy faces
    if (is_conservative) then
!$omp parallel do private(i,j,k)
       do k=ks-1,ke+1
          do j=js,je
             do i=is,ie+1
                ! make slxy, srxy by updating 1D extrapolation
                slxy(i,j,k) = slx(i,j,k) &
                     - (dt3/hy)*(simhy(i-1,j+1,k)*vmac(i-1,j+1,k) &
                     - simhy(i-1,j,k)*vmac(i-1,j,k))
                srxy(i,j,k) = srx(i,j,k) &
                     - (dt3/hy)*(simhy(i  ,j+1,k)*vmac(i  ,j+1,k) &
                     - simhy(i  ,j,k)*vmac(i  ,j,k))
             enddo
          enddo
       enddo
!$omp end parallel do
    else
!$omp parallel do private(i,j,k)
       do k=ks-1,ke+1
          do j=js,je
             do i=is,ie+1
                ! make slxy, srxy by updating 1D extrapolation
                slxy(i,j,k) = slx(i,j,k) &
                     - (dt6/hy)*(vmac(i-1,j+1,k)+vmac(i-1,j,k)) &
                     *(simhy(i-1,j+1,k)-simhy(i-1,j,k))
                srxy(i,j,k) = srx(i,j,k) &
                     - (dt6/hy)*(vmac(i  ,j+1,k)+vmac(i  ,j,k)) &
                     *(simhy(i  ,j+1,k)-simhy(i  ,j,k))
             enddo
          enddo
       enddo
!$omp end parallel do
    end if

    ! impose lo side bc's
    if (phys_bc(1,1) .eq. INLET) then
       slxy(is,js:je,ks-1:ke+1) = s(is-1,js:je,ks-1:ke+1,comp)
       srxy(is,js:je,ks-1:ke+1) = s(is-1,js:je,ks-1:ke+1,comp)
    else if (phys_bc(1,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 1) then
          slxy(is,js:je,ks-1:ke+1) = ZERO
          srxy(is,js:je,ks-1:ke+1) = ZERO
       else
          slxy(is,js:je,ks-1:ke+1) = srxy(is,js:je,ks-1:ke+1)
       end if
    else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slxy(is,js:je,ks-1:ke+1) = ZERO
          srxy(is,js:je,ks-1:ke+1) = ZERO
       else
          slxy(is,js:je,ks-1:ke+1) = srxy(is,js:je,ks-1:ke+1)
       end if
    else if (phys_bc(1,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 1) then
          slxy(is,js:je,ks-1:ke+1) = min(srxy(is,js:je,ks-1:ke+1),ZERO)
          srxy(is,js:je,ks-1:ke+1) = min(srxy(is,js:je,ks-1:ke+1),ZERO)
       else
          slxy(is,js:je,ks-1:ke+1) = srxy(is,js:je,ks-1:ke+1)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(1,2) .eq. INLET) then
       slxy(ie+1,js:je,ks-1:ke+1) = s(ie+1,js:je,ks-1:ke+1,comp)
       srxy(ie+1,js:je,ks-1:ke+1) = s(ie+1,js:je,ks-1:ke+1,comp)
    else if (phys_bc(1,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 1) then
          slxy(ie+1,js:je,ks-1:ke+1) = ZERO
          srxy(ie+1,js:je,ks-1:ke+1) = ZERO
       else
          srxy(ie+1,js:je,ks-1:ke+1) = slxy(ie+1,js:je,ks-1:ke+1)
       end if
    else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slxy(ie+1,js:je,ks-1:ke+1) = ZERO
          srxy(ie+1,js:je,ks-1:ke+1) = ZERO
       else
          srxy(ie+1,js:je,ks-1:ke+1) = slxy(ie+1,js:je,ks-1:ke+1)
       end if
    else if (phys_bc(1,2) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 1) then
          slxy(ie+1,js:je,ks-1:ke+1) = max(slxy(ie+1,js:je,ks-1:ke+1),ZERO)
          srxy(ie+1,js:je,ks-1:ke+1) = max(slxy(ie+1,js:je,ks-1:ke+1),ZERO)
       else
          srxy(ie+1,js:je,ks-1:ke+1) = slxy(ie+1,js:je,ks-1:ke+1)
       end if
    end if

!$omp parallel do private(i,j,k,savg)
    do k=ks-1,ke+1
       do j=js,je
          do i=is,ie+1
             ! make simhxy by solving Riemann problem
             simhxy(i,j,k) = merge(slxy(i,j,k),srxy(i,j,k),umac(i,j,k) .gt. ZERO)
             savg = HALF*(slxy(i,j,k)+srxy(i,j,k))
             simhxy(i,j,k) = merge(simhxy(i,j,k),savg,abs(umac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
!$omp end parallel do

    deallocate(slxy,srxy)

    ! loop over appropriate xz faces

    ! These are transverse terms.  The size allocation is tricky.
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    ! lo-1:hi+1 in unused direction
    allocate(slxz  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))
    allocate(srxz  (lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))
    allocate(simhxz(lo(1):hi(1)+1,lo(2)-1:hi(2)+1,lo(3):hi(3)))

    if (is_conservative) then
!$omp parallel do private(i,j,k)
       do k=ks,ke
          do j=js-1,je+1
             do i=is,ie+1
                ! make slxz, srxz by updating 1D extrapolation
                slxz(i,j,k) = slx(i,j,k) &
                     - (dt3/hz)*(simhz(i-1,j,k+1)*wmac(i-1,j,k+1) &
                     - simhz(i-1,j,k)*wmac(i-1,j,k))
                srxz(i,j,k) = srx(i,j,k) &
                     - (dt3/hz)*(simhz(i  ,j,k+1)*wmac(i  ,j,k+1) &
                     - simhz(i  ,j,k)*wmac(i  ,j,k))
             enddo
          enddo
       enddo
!$omp end parallel do
    else
!$omp parallel do private(i,j,k)
       do k=ks,ke
          do j=js-1,je+1
             do i=is,ie+1
                ! make slxz, srxz by updating 1D extrapolation
                slxz(i,j,k) = slx(i,j,k) &
                     - (dt6/hz)*(wmac(i-1,j,k+1)+wmac(i-1,j,k)) &
                     *(simhz(i-1,j,k+1)-simhz(i-1,j,k))
                srxz(i,j,k) = srx(i,j,k) &
                     - (dt6/hz)*(wmac(i  ,j,k+1)+wmac(i  ,j,k)) &
                     *(simhz(i  ,j,k+1)-simhz(i  ,j,k))
             enddo
          enddo
       enddo
!$omp end parallel do
    end if

    ! impose lo side bc's
    if (phys_bc(1,1) .eq. INLET) then
       slxz(is,js-1:je+1,ks:ke) = s(is-1,js-1:je+1,ks:ke,comp)
       srxz(is,js-1:je+1,ks:ke) = s(is-1,js-1:je+1,ks:ke,comp)
    else if (phys_bc(1,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 1) then
          slxz(is,js-1:je+1,ks:ke) = ZERO
          srxz(is,js-1:je+1,ks:ke) = ZERO
       else
          slxz(is,js-1:je+1,ks:ke) = srxz(is,js-1:je+1,ks:ke)
       end if
    else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slxz(is,js-1:je+1,ks:ke) = ZERO
          srxz(is,js-1:je+1,ks:ke) = ZERO
       else
          slxz(is,js-1:je+1,ks:ke) = srxz(is,js-1:je+1,ks:ke)
       end if
    else if (phys_bc(1,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 1) then
          slxz(is,js-1:je+1,ks:ke) = min(srxz(is,js-1:je+1,ks:ke),ZERO)
          srxz(is,js-1:je+1,ks:ke) = min(srxz(is,js-1:je+1,ks:ke),ZERO)
       else
          slxz(is,js-1:je+1,ks:ke) = srxz(is,js-1:je+1,ks:ke)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(1,2) .eq. INLET) then
       slxz(ie+1,js-1:je+1,ks:ke) = s(ie+1,js-1:je+1,ks:ke,comp)
       srxz(ie+1,js-1:je+1,ks:ke) = s(ie+1,js-1:je+1,ks:ke,comp)
    else if (phys_bc(1,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 1) then
          slxz(ie+1,js-1:je+1,ks:ke) = ZERO
          srxz(ie+1,js-1:je+1,ks:ke) = ZERO
       else
          srxz(ie+1,js-1:je+1,ks:ke) = slxz(ie+1,js-1:je+1,ks:ke)
       end if
    else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slxz(ie+1,js-1:je+1,ks:ke) = ZERO
          srxz(ie+1,js-1:je+1,ks:ke) = ZERO
       else
          srxz(ie+1,js-1:je+1,ks:ke) = slxz(ie+1,js-1:je+1,ks:ke)
       end if
    else if (phys_bc(1,2) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 1) then
          slxz(ie+1,js-1:je+1,ks:ke) = max(slxz(ie+1,js-1:je+1,ks:ke),ZERO)
          srxz(ie+1,js-1:je+1,ks:ke) = max(slxz(ie+1,js-1:je+1,ks:ke),ZERO)
       else
          srxz(ie+1,js-1:je+1,ks:ke) = slxz(ie+1,js-1:je+1,ks:ke)
       end if
    end if

!$omp parallel do private(i,j,k,savg)
    do k=ks,ke
       do j=js-1,je+1
          do i=is,ie+1
             ! make simhxz by solving Riemann problem
             simhxz(i,j,k) = merge(slxz(i,j,k),srxz(i,j,k),umac(i,j,k) .gt. ZERO)
             savg = HALF*(slxz(i,j,k)+srxz(i,j,k))
             simhxz(i,j,k) = merge(simhxz(i,j,k),savg,abs(umac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
!$omp end parallel do

    deallocate(slxz,srxz)

    ! loop over appropriate yx faces

    ! These are transverse terms.  The size allocation is tricky.
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    ! lo-1:hi+1 in unused direction
    allocate(slyx  (lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(sryx  (lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))
    allocate(simhyx(lo(1):hi(1),lo(2):hi(2)+1,lo(3)-1:hi(3)+1))

    if (is_conservative) then
!$omp parallel do private(i,j,k)
       do k=ks-1,ke+1
          do j=js,je+1
             do i=is,ie
                ! make slyx, sryx by updating 1D extrapolation
                slyx(i,j,k) = sly(i,j,k) &
                     - (dt3/hx)*(simhx(i+1,j-1,k)*umac(i+1,j-1,k) &
                     - simhx(i,j-1,k)*umac(i,j-1,k))
                sryx(i,j,k) = sry(i,j,k) &
                     - (dt3/hx)*(simhx(i+1,j  ,k)*umac(i+1,j  ,k) &
                     - simhx(i,j  ,k)*umac(i,j  ,k))
             enddo
          enddo
       enddo
!$omp end parallel do
    else
!$omp parallel do private(i,j,k)
       do k=ks-1,ke+1
          do j=js,je+1
             do i=is,ie
                ! make slyx, sryx by updating 1D extrapolation
                slyx(i,j,k) = sly(i,j,k) &
                     - (dt6/hx)*(umac(i+1,j-1,k)+umac(i,j-1,k)) &
                     *(simhx(i+1,j-1,k)-simhx(i,j-1,k))
                sryx(i,j,k) = sry(i,j,k) &
                     - (dt6/hx)*(umac(i+1,j  ,k)+umac(i,j  ,k)) &
                     *(simhx(i+1,j  ,k)-simhx(i,j  ,k))
             enddo
          enddo
       enddo
!$omp end parallel do
    end if

    ! impose lo side bc's
    if (phys_bc(2,1) .eq. INLET) then
       slyx(is:ie,js,ks-1:ke+1) = s(is:ie,js-1,ks-1:ke+1,comp)
       sryx(is:ie,js,ks-1:ke+1) = s(is:ie,js-1,ks-1:ke+1,comp)
    else if (phys_bc(2,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 2) then
          slyx(is:ie,js,ks-1:ke+1) = ZERO
          sryx(is:ie,js,ks-1:ke+1) = ZERO
       else
          slyx(is:ie,js,ks-1:ke+1) = sryx(is:ie,js,ks-1:ke+1)
       end if
    else if (phys_bc(2,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slyx(is:ie,js,ks-1:ke+1) = ZERO
          sryx(is:ie,js,ks-1:ke+1) = ZERO
       else
          slyx(is:ie,js,ks-1:ke+1) = sryx(is:ie,js,ks-1:ke+1)
       end if
    else if (phys_bc(2,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 2) then
          slyx(is:ie,js,ks-1:ke+1) = min(sryx(is:ie,js,ks-1:ke+1),ZERO)
          sryx(is:ie,js,ks-1:ke+1) = min(sryx(is:ie,js,ks-1:ke+1),ZERO)
       else
          slyx(is:ie,js,ks-1:ke+1) = sryx(is:ie,js,ks-1:ke+1)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(2,2) .eq. INLET) then
       slyx(is:ie,je+1,ks-1:ke+1) = s(is:ie,je+1,ks-1:ke+1,comp)
       sryx(is:ie,je+1,ks-1:ke+1) = s(is:ie,je+1,ks-1:ke+1,comp)
    else if (phys_bc(2,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 2) then
          slyx(is:ie,je+1,ks-1:ke+1) = ZERO
          sryx(is:ie,je+1,ks-1:ke+1) = ZERO
       else
          sryx(is:ie,je+1,ks-1:ke+1) = slyx(is:ie,je+1,ks-1:ke+1)
       end if
    else if (phys_bc(2,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slyx(is:ie,je+1,ks-1:ke+1) = ZERO
          sryx(is:ie,je+1,ks-1:ke+1) = ZERO
       else
          sryx(is:ie,je+1,ks-1:ke+1) = slyx(is:ie,je+1,ks-1:ke+1)
       end if
    else if (phys_bc(2,2) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 2) then
          slyx(is:ie,je+1,ks-1:ke+1) = max(slyx(is:ie,je+1,ks-1:ke+1),ZERO)
          sryx(is:ie,je+1,ks-1:ke+1) = max(slyx(is:ie,je+1,ks-1:ke+1),ZERO)
       else
          sryx(is:ie,je+1,ks-1:ke+1) = slyx(is:ie,je+1,ks-1:ke+1)
       end if
    end if

!$omp parallel do private(i,j,k,savg)
    do k=ks-1,ke+1
       do j=js,je+1
          do i=is,ie
             ! make simhyx by solving Riemann problem
             simhyx(i,j,k) = merge(slyx(i,j,k),sryx(i,j,k),vmac(i,j,k) .gt. ZERO)
             savg = HALF*(slyx(i,j,k)+sryx(i,j,k))
             simhyx(i,j,k) = merge(simhyx(i,j,k),savg,abs(vmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
!$omp end parallel do

    deallocate(slyx,sryx)

    ! loop over appropriate yz faces

    ! These are transverse terms.  The size allocation is tricky.
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    ! lo-1:hi+1 in unused direction
    allocate(slyz  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(sryz  (lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(simhyz(lo(1)-1:hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)))

    if (is_conservative) then
!$omp parallel do private(i,j,k)
       do k=ks,ke
          do j=js,je+1
             do i=is-1,ie+1
                ! make slyz, sryz by updating 1D extrapolation
                slyz(i,j,k) = sly(i,j,k) &
                     - (dt3/hz)*(simhz(i,j-1,k+1)*wmac(i,j-1,k+1) &
                     - simhz(i,j-1,k)*wmac(i,j-1,k))
                sryz(i,j,k) = sry(i,j,k) &
                     - (dt3/hz)*(simhz(i,j  ,k+1)*wmac(i,j  ,k+1) &
                     - simhz(i,j  ,k)*wmac(i,j  ,k))
             enddo
          enddo
       enddo
!$omp end parallel do
    else
!$omp parallel do private(i,j,k)
       do k=ks,ke
          do j=js,je+1
             do i=is-1,ie+1
                ! make slyz, sryz by updating 1D extrapolation
                slyz(i,j,k) = sly(i,j,k) &
                     - (dt6/hz)*(wmac(i,j-1,k+1)+wmac(i,j-1,k)) &
                     *(simhz(i,j-1,k+1)-simhz(i,j-1,k))
                sryz(i,j,k) = sry(i,j,k) &
                     - (dt6/hz)*(wmac(i,j  ,k+1)+wmac(i,j  ,k)) &
                     *(simhz(i,j  ,k+1)-simhz(i,j  ,k))
             enddo
          enddo
       enddo
!$omp end parallel do
    end if

    deallocate(simhz)

    ! impose lo side bc's
    if (phys_bc(2,1) .eq. INLET) then
       slyz(is-1:ie+1,js,ks:ke) = s(is-1:ie+1,js-1,ks:ke,comp)
       sryz(is-1:ie+1,js,ks:ke) = s(is-1:ie+1,js-1,ks:ke,comp)
    else if (phys_bc(2,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 2) then
          slyz(is-1:ie+1,js,ks:ke) = ZERO
          sryz(is-1:ie+1,js,ks:ke) = ZERO
       else
          slyz(is-1:ie+1,js,ks:ke) = sryz(is-1:ie+1,js,ks:ke)
       end if
    else if (phys_bc(2,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slyz(is-1:ie+1,js,ks:ke) = ZERO
          sryz(is-1:ie+1,js,ks:ke) = ZERO
       else
          slyz(is-1:ie+1,js,ks:ke) = sryz(is-1:ie+1,js,ks:ke)
       end if
    else if (phys_bc(2,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 2) then
          slyz(is-1:ie+1,js,ks:ke) = min(sryz(is-1:ie+1,js,ks:ke),ZERO)
          sryz(is-1:ie+1,js,ks:ke) = min(sryz(is-1:ie+1,js,ks:ke),ZERO)
       else
          slyz(is-1:ie+1,js,ks:ke) = sryz(is-1:ie+1,js,ks:ke)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(2,2) .eq. INLET) then
       slyz(is-1:ie+1,je+1,ks:ke) = s(is-1:ie+1,je+1,ks:ke,comp)
       sryz(is-1:ie+1,je+1,ks:ke) = s(is-1:ie+1,je+1,ks:ke,comp)
    else if (phys_bc(2,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 2) then
          slyz(is-1:ie+1,je+1,ks:ke) = ZERO
          sryz(is-1:ie+1,je+1,ks:ke) = ZERO
       else
          sryz(is-1:ie+1,je+1,ks:ke) = slyz(is-1:ie+1,je+1,ks:ke)
       end if
    else if (phys_bc(2,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slyz(is-1:ie+1,je+1,ks:ke) = ZERO
          sryz(is-1:ie+1,je+1,ks:ke) = ZERO
       else
          sryz(is-1:ie+1,je+1,ks:ke) = slyz(is-1:ie+1,je+1,ks:ke)
       end if
    else if (phys_bc(2,2) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 2) then
          slyz(is-1:ie+1,je+1,ks:ke) = max(slyz(is-1:ie+1,je+1,ks:ke),ZERO)
          sryz(is-1:ie+1,je+1,ks:ke) = max(slyz(is-1:ie+1,je+1,ks:ke),ZERO)
       else
          sryz(is-1:ie+1,je+1,ks:ke) = slyz(is-1:ie+1,je+1,ks:ke)
       end if
    end if

!$omp parallel do private(i,j,k,savg)
    do k=ks,ke
       do j=js,je+1
          do i=is-1,ie+1
             ! make simhyz by solving Riemann problem
             simhyz(i,j,k) = merge(slyz(i,j,k),sryz(i,j,k),vmac(i,j,k) .gt. ZERO)
             savg = HALF*(slyz(i,j,k)+sryz(i,j,k))
             simhyz(i,j,k) = merge(simhyz(i,j,k),savg,abs(vmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
!$omp end parallel do

    deallocate(slyz,sryz)

    ! loop over appropriate zx faces

    ! These are transverse terms.  The size allocation is tricky.
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    ! lo-1:hi+1 in unused direction
    allocate(slzx  (lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(srzx  (lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))
    allocate(simhzx(lo(1):hi(1),lo(2)-1:hi(2)+1,lo(3):hi(3)+1))

    if (is_conservative) then
!$omp parallel do private(i,j,k)
       do k=ks,ke+1
          do j=js-1,je+1
             do i=is,ie
                ! make slzx, srzx by updating 1D extrapolation
                slzx(i,j,k) = slz(i,j,k) &
                     - (dt3/hx)*(simhx(i+1,j,k-1)*umac(i+1,j,k-1) &
                     - simhx(i,j,k-1)*umac(i,j,k-1))
                srzx(i,j,k) = srz(i,j,k) &
                     - (dt3/hx)*(simhx(i+1,j,k  )*umac(i+1,j,k  ) &
                     - simhx(i,j,k  )*umac(i,j,k  ))
             enddo
          enddo
       enddo
!$omp end parallel do
    else
!$omp parallel do private(i,j,k)
       do k=ks,ke+1
          do j=js-1,je+1
             do i=is,ie
                ! make slzx, srzx by updating 1D extrapolation
                slzx(i,j,k) = slz(i,j,k) &
                     - (dt6/hx)*(umac(i+1,j,k-1)+umac(i,j,k-1)) &
                     *(simhx(i+1,j,k-1)-simhx(i,j,k-1))
                srzx(i,j,k) = srz(i,j,k) &
                     - (dt6/hx)*(umac(i+1,j,k  )+umac(i,j,k  )) &
                     *(simhx(i+1,j,k  )-simhx(i,j,k  ))
             enddo
          enddo
       end do
!$omp end parallel do
    end if

    deallocate(simhx)

    ! impose lo side bc's
    if (phys_bc(3,1) .eq. INLET) then
       slzx(is:ie,js-1:je+1,ks) = s(is:ie,js-1:je+1,ks-1,comp)
       srzx(is:ie,js-1:je+1,ks) = s(is:ie,js-1:je+1,ks-1,comp)
    else if (phys_bc(3,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 3) then
          slzx(is:ie,js-1:je+1,ks) = ZERO
          srzx(is:ie,js-1:je+1,ks) = ZERO
       else
          slzx(is:ie,js-1:je+1,ks) = srzx(is:ie,js-1:je+1,ks)
       end if
    else if (phys_bc(3,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slzx(is:ie,js-1:je+1,ks) = ZERO
          srzx(is:ie,js-1:je+1,ks) = ZERO
       else
          slzx(is:ie,js-1:je+1,ks) = srzx(is:ie,js-1:je+1,ks)
       end if
    else if (phys_bc(3,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 3) then
          slzx(is:ie,js-1:je+1,ks) = min(srzx(is:ie,js-1:je+1,ks),ZERO)
          srzx(is:ie,js-1:je+1,ks) = min(srzx(is:ie,js-1:je+1,ks),ZERO)
       else
          slzx(is:ie,js-1:je+1,ks) = srzx(is:ie,js-1:je+1,ks)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(3,2) .eq. INLET) then
       slzx(is:ie,js-1:je+1,ke+1) = s(is:ie,js-1:je+1,ke+1,comp)
       srzx(is:ie,js-1:je+1,ke+1) = s(is:ie,js-1:je+1,ke+1,comp)
    else if (phys_bc(3,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 3) then
          slzx(is:ie,js-1:je+1,ke+1) = ZERO
          srzx(is:ie,js-1:je+1,ke+1) = ZERO
       else
          srzx(is:ie,js-1:je+1,ke+1) = slzx(is:ie,js-1:je+1,ke+1)
       end if
    else if (phys_bc(3,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slzx(is:ie,js-1:je+1,ke+1) = ZERO
          srzx(is:ie,js-1:je+1,ke+1) = ZERO
       else
          srzx(is:ie,js-1:je+1,ke+1) = slzx(is:ie,js-1:je+1,ke+1)
       end if
    else if (phys_bc(3,2) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 3) then
          slzx(is:ie,js-1:je+1,ke+1) = max(slzx(is:ie,js-1:je+1,ke+1),ZERO)
          srzx(is:ie,js-1:je+1,ke+1) = max(slzx(is:ie,js-1:je+1,ke+1),ZERO)
       else
          srzx(is:ie,js-1:je+1,ke+1) = slzx(is:ie,js-1:je+1,ke+1)
       end if
    end if

!$omp parallel do private(i,j,k,savg)
    do k=ks,ke+1
       do j=js-1,je+1
          do i=is,ie
             ! make simhzx by solving Riemann problem
             simhzx(i,j,k) = merge(slzx(i,j,k),srzx(i,j,k),wmac(i,j,k) .gt. ZERO)
             savg = HALF*(slzx(i,j,k)+srzx(i,j,k))
             simhzx(i,j,k) = merge(simhzx(i,j,k),savg,abs(wmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
!$omp end parallel do

    deallocate(slzx,srzx)

    ! loop over appropriate zy faces

    ! These are transverse terms.  The size allocation is tricky.
    ! lo:hi+1 in normal direction
    ! lo:hi in transverse direction
    ! lo-1:hi+1 in unused direction
    allocate(slzy  (lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))
    allocate(srzy  (lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))
    allocate(simhzy(lo(1)-1:hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1))

    if (is_conservative) then
!$omp parallel do private(i,j,k)
       do k=ks,ke+1
          do j=js,je
             do i=is-1,ie+1
                ! make slzy, srzy by updating 1D extrapolation
                slzy(i,j,k) = slz(i,j,k) &
                     - (dt3/hy)*(simhy(i,j+1,k-1)*vmac(i,j+1,k-1) &
                     - simhy(i,j,k-1)*vmac(i,j,k-1))
                srzy(i,j,k) = srz(i,j,k) &
                     - (dt3/hy)*(simhy(i,j+1,k  )*vmac(i,j+1,k  ) &
                     - simhy(i,j,k  )*vmac(i,j,k  ))
             enddo
          enddo
       enddo
!$omp end parallel do
    else
!$omp parallel do private(i,j,k)
       do k=ks,ke+1
          do j=js,je
             do i=is-1,ie+1
                ! make slzy, srzy by updating 1D extrapolation
                slzy(i,j,k) = slz(i,j,k) &
                     - (dt6/hy)*(vmac(i,j+1,k-1)+vmac(i,j,k-1)) &
                     *(simhy(i,j+1,k-1)-simhy(i,j,k-1))
                srzy(i,j,k) = srz(i,j,k) &
                     - (dt6/hy)*(vmac(i,j+1,k  )+vmac(i,j,k  )) &
                     *(simhy(i,j+1,k  )-simhy(i,j,k  ))
             enddo
          enddo
       enddo
!$omp end parallel do
    end if

    deallocate(simhy)

    ! impose lo side bc's
    if (phys_bc(3,1) .eq. INLET) then
       slzy(is-1:ie+1,js:je,ks) = s(is-1:ie+1,js:je,ks-1,comp)
       srzy(is-1:ie+1,js:je,ks) = s(is-1:ie+1,js:je,ks-1,comp)
    else if (phys_bc(3,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 3) then
          slzy(is-1:ie+1,js:je,ks) = ZERO
          srzy(is-1:ie+1,js:je,ks) = ZERO
       else
          slzy(is-1:ie+1,js:je,ks) = srzy(is-1:ie+1,js:je,ks)
       end if
    else if (phys_bc(3,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slzy(is-1:ie+1,js:je,ks) = ZERO
          srzy(is-1:ie+1,js:je,ks) = ZERO
       else
          slzy(is-1:ie+1,js:je,ks) = srzy(is-1:ie+1,js:je,ks)
       end if
    else if (phys_bc(3,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 3) then
          slzy(is-1:ie+1,js:je,ks) = min(srzy(is-1:ie+1,js:je,ks),ZERO)
          srzy(is-1:ie+1,js:je,ks) = min(srzy(is-1:ie+1,js:je,ks),ZERO)
       else
          slzy(is-1:ie+1,js:je,ks) = srzy(is-1:ie+1,js:je,ks)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(3,2) .eq. INLET) then
       slzy(is-1:ie+1,js:je,ke+1) = s(is-1:ie+1,js:je,ke+1,comp)
       srzy(is-1:ie+1,js:je,ke+1) = s(is-1:ie+1,js:je,ke+1,comp)
    else if (phys_bc(3,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 3) then
          slzy(is-1:ie+1,js:je,ke+1) = ZERO
          srzy(is-1:ie+1,js:je,ke+1) = ZERO
       else
          srzy(is-1:ie+1,js:je,ke+1) = slzy(is-1:ie+1,js:je,ke+1)
       end if
    else if (phys_bc(3,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          slzy(is-1:ie+1,js:je,ke+1) = ZERO
          srzy(is-1:ie+1,js:je,ke+1) = ZERO
       else
          srzy(is-1:ie+1,js:je,ke+1) = slzy(is-1:ie+1,js:je,ke+1)
       end if
    else if (phys_bc(3,2) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 3) then
          slzy(is-1:ie+1,js:je,ke+1) = max(slzy(is-1:ie+1,js:je,ke+1),ZERO)
          srzy(is-1:ie+1,js:je,ke+1) = max(slzy(is-1:ie+1,js:je,ke+1),ZERO)
       else
          srzy(is-1:ie+1,js:je,ke+1) = slzy(is-1:ie+1,js:je,ke+1)
       end if
    end if

!$omp parallel do private(i,j,k,savg)
    do k=ks,ke+1
       do j=js,je
          do i=is-1,ie+1
             ! make simhzy by solving Riemann problem
             simhzy(i,j,k) = merge(slzy(i,j,k),srzy(i,j,k),wmac(i,j,k) .gt. ZERO)
             savg = HALF*(slzy(i,j,k)+srzy(i,j,k))
             simhzy(i,j,k) = merge(simhzy(i,j,k),savg,abs(wmac(i,j,k)) .gt. rel_eps)
          enddo
       enddo
    enddo
!$omp end parallel do

    deallocate(slzy,srzy)

    !******************************************************************
    ! Create sedgelx, etc.
    !******************************************************************

    ! Final edge states.
    ! lo:hi+1 in the normal direction
    ! lo:hi in the transverse directions
    allocate(sedgelx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)))
    allocate(sedgerx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3)))

    ! loop over appropriate x-faces
    if (is_conservative) then
!$omp parallel do private(i,j,k)
       do k=ks,ke
          do j=js,je
             do i=is,ie+1
                ! make sedgelx, sedgerx
                sedgelx(i,j,k) = slx(i,j,k) &
                     - (dt2/hy)*(simhyz(i-1,j+1,k  )*vmac(i-1,j+1,k  ) &
                     - simhyz(i-1,j,k)*vmac(i-1,j,k)) &
                     - (dt2/hz)*(simhzy(i-1,j  ,k+1)*wmac(i-1,j  ,k+1) &
                     - simhzy(i-1,j,k)*wmac(i-1,j,k)) &
                     - (dt2/hx)*s(i-1,j,k,comp)*(umac(i  ,j,k)-umac(i-1,j,k)) &
                     + dt2*force(i-1,j,k,comp)
                sedgerx(i,j,k) = srx(i,j,k) &
                     - (dt2/hy)*(simhyz(i  ,j+1,k  )*vmac(i  ,j+1,  k) &
                     - simhyz(i  ,j,k)*vmac(i  ,j,k)) &
                     - (dt2/hz)*(simhzy(i  ,j  ,k+1)*wmac(i  ,j  ,k+1) &
                     - simhzy(i  ,j,k)*wmac(i  ,j,k)) &
                     - (dt2/hx)*s(i  ,j,k,comp)*(umac(i+1,j,k)-umac(i  ,j,k)) &
                     + dt2*force(i  ,j,k,comp)
             end do
          end do
       end do
!$omp end parallel do
    else
!$omp parallel do private(i,j,k)
       do k=ks,ke
          do j=js,je
             do i=is,ie+1
                ! make sedgelx, sedgerx
                sedgelx(i,j,k) = slx(i,j,k) &
                     - (dt4/hy)*(vmac(i-1,j+1,k  )+vmac(i-1,j,k))* &
                     (simhyz(i-1,j+1,k  )-simhyz(i-1,j,k)) &
                     - (dt4/hz)*(wmac(i-1,j  ,k+1)+wmac(i-1,j,k))* &
                     (simhzy(i-1,j  ,k+1)-simhzy(i-1,j,k)) &
                     + dt2*force(i-1,j,k,comp)
                sedgerx(i,j,k) = srx(i,j,k) &
                     - (dt4/hy)*(vmac(i  ,j+1,k  )+vmac(i  ,j,k))* &
                     (simhyz(i  ,j+1,k  )-simhyz(i  ,j,k)) &
                     - (dt4/hz)*(wmac(i  ,j  ,k+1)+wmac(i  ,j,k))* &
                     (simhzy(i  ,j  ,k+1)-simhzy(i  ,j,k)) &
                     + dt2*force(i  ,j,k,comp)
             end do
          end do
       end do
!$omp end parallel do
    end if

    deallocate(slx,srx,simhyz,simhzy)

    if (is_vel) then
       
       ! add the (Utilde . e_r) d w_0 /dr e_r term here
       if (spherical .eq. 0 .and. comp .eq. 3) then

!$omp parallel do private(i,j,k)
          do k=ks,ke
             do j=js,je
                do i=is,ie+1
                   ! wmac contains w0 so we need to subtract it off
                   sedgelx(i,j,k) = sedgelx(i,j,k) &
                        -(dt4/hz)*(wmac(i-1,j,k)-w0(k)+wmac(i-1,j,k+1)-w0(k+1))* &
                        (w0(k+1)-w0(k))
                   sedgerx(i,j,k) = sedgerx(i,j,k) &
                        -(dt4/hz)*(wmac(i  ,j,k)-w0(k)+wmac(i  ,j,k+1)-w0(k+1))* &
                        (w0(k+1)-w0(k))
                enddo
             enddo
          enddo
!$omp end parallel do

       else if (spherical .eq. 1) then

!$omp parallel do private(i,j,k,Ut_dot_er)
          do k=ks,ke
             do j=js,je
                do i=is,ie+1
                   ! u/v/wmac contain w0, so we need to subtract it off.
                   ! left face centered about (i-1,j,k)
                   Ut_dot_er = (HALF*(umac(i-1,j,k) + umac(i,j,k)) - &
                        HALF*(w0macx(i-1,j,k)+w0macx(i,j,k)))*normal(i-1,j,k,1) + &
                        (HALF*(vmac(i-1,j,k) + vmac(i-1,j+1,k)) - &
                        HALF*(w0macy(i-1,j,k)+w0macy(i-1,j+1,k)))*normal(i-1,j,k,2) + &
                        (HALF*(wmac(i-1,j,k) + wmac(i-1,j,k+1)) - &
                        HALF*(w0macz(i-1,j,k)+w0macz(i-1,j,k+1)))*normal(i-1,j,k,3)

                   sedgelx(i,j,k) = sedgelx(i,j,k) &
                        - dt2*Ut_dot_er*gradw0_cart(i-1,j,k)*normal(i-1,j,k,comp)

                   ! right face centered about (i,j,k)
                   Ut_dot_er = (HALF*(umac(i,j,k) + umac(i+1,j,k)) - &
                        HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)))*normal(i,j,k,1) + &
                        (HALF*(vmac(i,j,k) + vmac(i,j+1,k)) - &
                        HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)))*normal(i,j,k,2) + &
                        (HALF*(wmac(i,j,k) + wmac(i,j,k+1)) - &
                        HALF*(w0macz(i,j,k)+w0macz(i,j,k+1)))*normal(i,j,k,3)

                   sedgerx(i,j,k) = sedgerx(i,j,k) &
                        - dt2*Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,comp)
                enddo
             enddo
          enddo
!$omp end parallel do

       end if

    end if

!$omp parallel do private(i,j,k,savg)
    do k=ks,ke
       do j=js,je
          do i=is,ie+1
             ! make sedgex by solving Riemann problem
             ! boundary conditions enforced outside of i,j,k loop
             sedgex(i,j,k,comp) = merge(sedgelx(i,j,k),sedgerx(i,j,k),umac(i,j,k) .gt. ZERO)
             savg = HALF*(sedgelx(i,j,k)+sedgerx(i,j,k))
             sedgex(i,j,k,comp) = merge(sedgex(i,j,k,comp),savg,abs(umac(i,j,k)).gt.rel_eps)
          enddo
       enddo
    enddo
!$omp end parallel do

    ! impose lo side bc's
    if (phys_bc(1,1) .eq. INLET) then
       sedgex(is,js:je,ks:ke,comp) = s(is-1,js:je,ks:ke,comp)
    else if (phys_bc(1,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 1) then
          sedgex(is,js:je,ks:ke,comp) = ZERO
       else
          sedgex(is,js:je,ks:ke,comp) = sedgerx(is,js:je,ks:ke)
       end if
    else if (phys_bc(1,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          sedgex(is,js:je,ks:ke,comp) = ZERO
       else
          sedgex(is,js:je,ks:ke,comp) = sedgerx(is,js:je,ks:ke)
       end if
    else if (phys_bc(1,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 1) then
          sedgex(is,js:je,ks:ke,comp) = min(sedgerx(is,js:je,ks:ke),ZERO)
       else
          sedgex(is,js:je,ks:ke,comp) = sedgerx(is,js:je,ks:ke)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(1,2) .eq. INLET) then
       sedgex(ie+1,js:je,ks:ke,comp) = s(ie+1,js:je,ks:ke,comp)
    else if (phys_bc(1,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 1) then
          sedgex(ie+1,js:je,ks:ke,comp) = ZERO
       else
          sedgex(ie+1,js:je,ks:ke,comp) = sedgelx(ie+1,js:je,ks:ke)
       end if
    else if (phys_bc(1,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          sedgex(ie+1,js:je,ks:ke,comp) = ZERO
       else
          sedgex(ie+1,js:je,ks:ke,comp) = sedgelx(ie+1,js:je,ks:ke)
       end if
    else if (phys_bc(1,2) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 1) then
          sedgex(ie+1,js:je,ks:ke,comp) = max(sedgelx(ie+1,js:je,ks:ke),ZERO)
       else
          sedgex(ie+1,js:je,ks:ke,comp) = sedgelx(ie+1,js:je,ks:ke)
       end if
    end if

    deallocate(sedgelx,sedgerx)

    allocate(sedgely(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3)))
    allocate(sedgery(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3)))

    ! loop over appropriate y-faces
    if (is_conservative) then
!$omp parallel do private(i,j,k)
       do k=ks,ke
          do j=js,je+1
             do i=is,ie
                ! make sedgely, sedgery
                sedgely(i,j,k) = sly(i,j,k) &
                     - (dt2/hx)*(simhxz(i+1,j-1,k  )*umac(i+1,j-1,k  ) &
                     - simhxz(i,j-1,k)*umac(i,j-1,k)) &
                     - (dt2/hz)*(simhzx(i  ,j-1,k+1)*wmac(i  ,j-1,k+1) &
                     - simhzx(i,j-1,k)*wmac(i,j-1,k)) &
                     - (dt2/hy)*s(i,j-1,k,comp)*(vmac(i,j  ,k)-vmac(i,j-1,k)) &
                     + dt2*force(i,j-1,k,comp)
                sedgery(i,j,k) = sry(i,j,k) &
                     - (dt2/hx)*(simhxz(i+1,j  ,k  )*umac(i+1,j  ,k  ) &
                     - simhxz(i,j  ,k)*umac(i,j  ,k)) &
                     - (dt2/hz)*(simhzx(i  ,j  ,k+1)*wmac(i  ,j  ,k+1) &
                     - simhzx(i,j  ,k)*wmac(i,j  ,k)) &
                     - (dt2/hy)*s(i,j  ,k,comp)*(vmac(i,j+1,k)-vmac(i,j  ,k)) &
                     + dt2*force(i,j  ,k,comp)
             end do
          end do
       end do
!$omp end parallel do
    else
!$omp parallel do private(i,j,k)
       do k=ks,ke
          do j=js,je+1
             do i=is,ie
                ! make sedgely, sedgery
                sedgely(i,j,k) = sly(i,j,k) &
                     - (dt4/hx)*(umac(i+1,j-1,k  )+umac(i,j-1,k))* &
                     (simhxz(i+1,j-1,k  )-simhxz(i,j-1,k)) &
                     - (dt4/hz)*(wmac(i  ,j-1,k+1)+wmac(i,j-1,k))* &
                     (simhzx(i  ,j-1,k+1)-simhzx(i,j-1,k)) &
                     + dt2*force(i,j-1,k,comp)
                sedgery(i,j,k) = sry(i,j,k) &
                     - (dt4/hx)*(umac(i+1,j  ,k  )+umac(i,j  ,k))* &
                     (simhxz(i+1,j  ,k  )-simhxz(i,j  ,k)) &
                     - (dt4/hz)*(wmac(i  ,j  ,k+1)+wmac(i,j  ,k))* &
                     (simhzx(i  ,j  ,k+1)-simhzx(i,j  ,k)) &
                     + dt2*force(i,j  ,k,comp)
             end do
          end do
       end do
!$omp end parallel do
    end if

    deallocate(sly,sry,simhxz,simhzx)

    if (is_vel) then
       
       ! add the (Utilde . e_r) d w_0 /dr e_r term here
       if (spherical .eq. 0 .and. comp .eq. 3) then

!$omp parallel do private(i,j,k)
          do k=ks,ke
             do j=js,je+1
                do i=is,ie
                   ! wmac contains w0 so we need to subtract it off
                   sedgely(i,j,k) = sedgely(i,j,k) &
                        -(dt4/hz)*(wmac(i,j-1,k)-w0(k)+wmac(i,j-1,k+1)-w0(k+1)) &
                        *(w0(k+1)-w0(k))
                   sedgery(i,j,k) = sedgery(i,j,k) &
                        -(dt4/hz)*(wmac(i,j  ,k)-w0(k)+wmac(i,j  ,k+1)-w0(k+1)) &
                        *(w0(k+1)-w0(k))
                enddo
             enddo
          enddo
!$omp end parallel do

       else if (spherical .eq. 1) then

!$omp parallel do private(i,j,k,Ut_dot_er)
          do k=ks,ke
             do j=js,je+1
                do i=is,ie
                   ! u/v/wmac contain w0, so we need to subtract it off.  
                   ! left face centered about (i,j-1,k)
                   Ut_dot_er = (HALF*(umac(i,j-1,k) + umac(i+1,j-1,k)) - &
                        HALF*(w0macx(i,j-1,k)+w0macx(i+1,j-1,k)))*normal(i,j-1,k,1) + &
                        (HALF*(vmac(i,j-1,k) + vmac(i,j,k)) - &
                        HALF*(w0macy(i,j-1,k)+w0macy(i,j,k)))*normal(i,j-1,k,2) + &
                        (HALF*(wmac(i,j-1,k) + wmac(i,j-1,k+1)) - &
                        HALF*(w0macz(i,j-1,k)+w0macz(i,j-1,k+1)))*normal(i,j-1,k,3)

                   sedgely(i,j,k) = sedgely(i,j,k) &
                        - dt2*Ut_dot_er*gradw0_cart(i,j-1,k)*normal(i,j-1,k,comp)

                   ! right face centered about (i,j,k)
                   Ut_dot_er = (HALF*(umac(i,j,k) + umac(i+1,j,k)) - &
                        HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)))*normal(i,j,k,1) + &
                        (HALF*(vmac(i,j,k) + vmac(i,j+1,k)) - &
                        HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)))*normal(i,j,k,2) + &
                        (HALF*(wmac(i,j,k) + wmac(i,j,k+1)) - &
                        HALF*(w0macz(i,j,k)+w0macz(i,j,k+1)))*normal(i,j,k,3)

                   sedgery(i,j,k) = sedgery(i,j,k) &
                        - dt2*Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,comp)
                enddo
             enddo
          enddo
!$omp end parallel do

       end if

    end if

!$omp parallel do private(i,j,k,savg)
    do k=ks,ke
       do j=js,je+1
          do i=is,ie
             ! make sedgey by solving Riemann problem
             ! boundary conditions enforced outside of i,j,k loop
             sedgey(i,j,k,comp) = merge(sedgely(i,j,k),sedgery(i,j,k),vmac(i,j,k) .gt. ZERO)
             savg = HALF*(sedgely(i,j,k)+sedgery(i,j,k))
             sedgey(i,j,k,comp) = merge(sedgey(i,j,k,comp),savg,abs(vmac(i,j,k)).gt.rel_eps)
          enddo
       enddo
    enddo
!$omp end parallel do

    ! impose lo side bc's
    if (phys_bc(2,1) .eq. INLET) then
       sedgey(is:ie,js,ks:ke,comp) = s(is:ie,js-1,ks:ke,comp)
    else if (phys_bc(2,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 2) then
          sedgey(is:ie,js,ks:ke,comp) = ZERO
       else
          sedgey(is:ie,js,ks:ke,comp) = sedgery(is:ie,js,ks:ke)
       end if
    else if (phys_bc(2,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          sedgey(is:ie,js,ks:ke,comp) = ZERO
       else
          sedgey(is:ie,js,ks:ke,comp) = sedgery(is:ie,js,ks:ke)
       end if
    else if (phys_bc(2,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 2) then
          sedgey(is:ie,js,ks:ke,comp) = min(sedgery(is:ie,js,ks:ke),ZERO)
       else
          sedgey(is:ie,js,ks:ke,comp) = sedgery(is:ie,js,ks:ke)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(2,2) .eq. INLET) then
       sedgey(is:ie,je+1,ks:ke,comp) = s(is:ie,je+1,ks:ke,comp)
    else if (phys_bc(2,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 2) then
          sedgey(is:ie,je+1,ks:ke,comp) = ZERO
       else
          sedgey(is:ie,je+1,ks:ke,comp) = sedgely(is:ie,je+1,ks:ke)
       end if
    else if (phys_bc(2,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          sedgey(is:ie,je+1,ks:ke,comp) = ZERO
       else
          sedgey(is:ie,je+1,ks:ke,comp) = sedgely(is:ie,je+1,ks:ke)
       end if
    else if (phys_bc(2,2) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 2) then
          sedgey(is:ie,je+1,ks:ke,comp) = max(sedgely(is:ie,je+1,ks:ke),ZERO)
       else
          sedgey(is:ie,je+1,ks:ke,comp) = sedgely(is:ie,je+1,ks:ke)
       end if
    end if

    deallocate(sedgely,sedgery)

    allocate(sedgelz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1))
    allocate(sedgerz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1))

    ! loop over appropriate z-faces
    if (is_conservative) then
!$omp parallel do private(i,j,k)
       do k=ks,ke+1
          do j=js,je
             do i=is,ie
                ! make sedgelz, sedgerz
                sedgelz(i,j,k) = slz(i,j,k) &
                     - (dt2/hx)*(simhxy(i+1,j  ,k-1)*umac(i+1,j  ,k-1) &
                     - simhxy(i,j,k-1)*umac(i,j,k-1)) &
                     - (dt2/hy)*(simhyx(i  ,j+1,k-1)*vmac(i  ,j+1,k-1) &
                     - simhyx(i,j,k-1)*vmac(i,j,k-1)) &
                     - (dt2/hz)*s(i,j,k-1,comp)*(wmac(i,j,k  )-wmac(i,j,k-1)) &
                     + dt2*force(i,j,k-1,comp)
                sedgerz(i,j,k) = srz(i,j,k) &
                     - (dt2/hx)*(simhxy(i+1,j  ,k  )*umac(i+1,j  ,k  ) &
                     - simhxy(i,j,k  )*umac(i,j,k  )) &
                     - (dt2/hy)*(simhyx(i  ,j+1,k  )*vmac(i  ,j+1,k  ) &
                     - simhyx(i,j,k  )*vmac(i,j,k  )) &
                     - (dt2/hz)*s(i,j,k  ,comp)*(wmac(i,j,k+1)-wmac(i,j,k  )) &
                     + dt2*force(i,j,k  ,comp)
             end do
          end do
       end do
!$omp end parallel do
    else
!$omp parallel do private(i,j,k)
       do k=ks,ke+1
          do j=js,je
             do i=is,ie
                ! make sedgelz, sedgerz
                sedgelz(i,j,k) = slz(i,j,k) &
                     - (dt4/hx)*(umac(i+1,j  ,k-1)+umac(i,j,k-1)) &
                     *(simhxy(i+1,j  ,k-1)-simhxy(i,j,k-1)) &
                     - (dt4/hy)*(vmac(i  ,j+1,k-1)+vmac(i,j,k-1)) &
                     *(simhyx(i  ,j+1,k-1)-simhyx(i,j,k-1)) &
                     + dt2*force(i,j,k-1,comp)
                sedgerz(i,j,k) = srz(i,j,k) &
                     - (dt4/hx)*(umac(i+1,j  ,k  )+umac(i,j,k  )) &
                     *(simhxy(i+1,j  ,k  )-simhxy(i,j,k  )) &
                     - (dt4/hy)*(vmac(i  ,j+1,k  )+vmac(i,j,k  )) &
                     *(simhyx(i  ,j+1,k  )-simhyx(i,j,k  )) &
                     + dt2*force(i,j,k  ,comp)
             end do
          end do
       end do
!$omp end parallel do
    end if

    deallocate(slz,srz,simhxy,simhyx)

    if (is_vel) then
       
       ! add the (Utilde . e_r) d w_0 /dr e_r term here
       if (spherical .eq. 0 .and. comp .eq. 3) then

!$omp parallel do private(i,j,k)
          do k=ks,ke+1
             do j=js,je
                do i=is,ie
                   ! wmac contains w0 so we need to subtract it off
                   if (k .eq. 0) then
                      ! sedgelz unchanged since dw_0 / dr = 0
                      sedgerz(i,j,k) = sedgerz(i,j,k) &
                           -(dt4/hz)*(wmac(i,j,k)-w0(k)+wmac(i,j,k+1)-w0(k+1)) &
                           *(w0(k+1)-w0(k  ))                      
                   else if (k .eq. nr(n)) then
                      sedgelz(i,j,k) = sedgelz(i,j,k) &
                           -(dt4/hz)*(wmac(i,j,k)-w0(k)+wmac(i,j,k-1)-w0(k-1)) &
                           *(w0(k  )-w0(k-1))
                      ! sedgerz unchanged since dw_0 / dr = 0
                   else
                      sedgelz(i,j,k) = sedgelz(i,j,k) &
                           -(dt4/hz)*(wmac(i,j,k)-w0(k)+wmac(i,j,k-1)-w0(k-1)) &
                           *(w0(k  )-w0(k-1))
                      sedgerz(i,j,k) = sedgerz(i,j,k) &
                           -(dt4/hz)*(wmac(i,j,k)-w0(k)+wmac(i,j,k+1)-w0(k+1)) &
                           *(w0(k+1)-w0(k  ))
                   end if
                enddo
             enddo
          enddo
!$omp end parallel do

       else if (spherical .eq. 1) then

!$omp parallel do private(i,j,k,Ut_dot_er)
          do k=ks,ke+1
             do j=js,je
                do i=is,ie
                   ! u/v/wmac contain w0, so we need to subtract it off.  
                   ! left face centered about (i,j,k-1)
                   Ut_dot_er = (HALF*(umac(i,j,k-1) + umac(i+1,j,k-1)) - &
                        HALF*(w0macx(i,j,k-1)+w0macx(i+1,j,k-1)))*normal(i,j,k-1,1) + &
                        (HALF*(vmac(i,j,k-1) + vmac(i,j+1,k-1)) - &
                        HALF*(w0macy(i,j,k-1)+w0macy(i,j+1,k-1)))*normal(i,j,k-1,2) + &
                        (HALF*(wmac(i,j,k-1) + wmac(i,j,k)) - &
                        HALF*(w0macz(i,j,k-1)+w0macz(i,j,k)))*normal(i,j,k-1,3)

                   sedgelz(i,j,k) = sedgelz(i,j,k) &
                        - dt2*Ut_dot_er*gradw0_cart(i,j,k-1)*normal(i,j,k-1,comp)

                   ! right face centered about (i,j,k)
                   Ut_dot_er = (HALF*(umac(i,j,k) + umac(i+1,j,k)) - &
                        HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)))*normal(i,j,k,1) + &
                        (HALF*(vmac(i,j,k) + vmac(i,j+1,k)) - &
                        HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)))*normal(i,j,k,2) + &
                        (HALF*(wmac(i,j,k) + wmac(i,j,k+1)) - &
                        HALF*(w0macz(i,j,k)+w0macz(i,j,k+1)))*normal(i,j,k,3)

                   sedgerz(i,j,k) = sedgerz(i,j,k) &
                        - dt2*Ut_dot_er*gradw0_cart(i,j,k)*normal(i,j,k,comp)
                enddo
             enddo
          enddo
!$omp end parallel do

       end if

    end if

!$omp parallel do private(i,j,k,savg)
    do k=ks,ke+1
       do j=js,je
          do i=is,ie
             ! make sedgez by solving Riemann problem
             ! boundary conditions enforced outside of i,j,k loop
             sedgez(i,j,k,comp) = merge(sedgelz(i,j,k),sedgerz(i,j,k),wmac(i,j,k) .gt. ZERO)
             savg = HALF*(sedgelz(i,j,k)+sedgerz(i,j,k))
             sedgez(i,j,k,comp) = merge(sedgez(i,j,k,comp),savg,abs(wmac(i,j,k)).gt.rel_eps)
          enddo
       enddo
    enddo
!$omp end parallel do

    ! impose lo side bc's
    if (phys_bc(3,1) .eq. INLET) then
       sedgez(is:ie,js:je,ks,comp) = s(is:ie,js:je,ks-1,comp)
    else if (phys_bc(3,1) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 3) then
          sedgez(is:ie,js:je,ks,comp) = ZERO
       else
          sedgez(is:ie,js:je,ks,comp) = sedgerz(is:ie,js:je,ks)
       end if
    else if (phys_bc(3,1) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          sedgez(is:ie,js:je,ks,comp) = ZERO
       else
          sedgez(is:ie,js:je,ks,comp) = sedgerz(is:ie,js:je,ks)
       end if
    else if (phys_bc(3,1) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 3) then
          sedgez(is:ie,js:je,ks,comp) = min(sedgerz(is:ie,js:je,ks),ZERO)
       else
          sedgez(is:ie,js:je,ks,comp) = sedgerz(is:ie,js:je,ks)
       end if
    end if

    ! impose hi side bc's
    if (phys_bc(3,2) .eq. INLET) then
       sedgez(is:ie,js:je,ke+1,comp) = s(is:ie,js:je,ke+1,comp)
    else if (phys_bc(3,2) .eq. SLIP_WALL) then
       if (is_vel .and. comp .eq. 3) then
          sedgez(is:ie,js:je,ke+1,comp) = ZERO
       else
          sedgez(is:ie,js:je,ke+1,comp) = sedgelz(is:ie,js:je,ke+1)
       end if
    else if (phys_bc(3,2) .eq. NO_SLIP_WALL) then
       if (is_vel) then
          sedgez(is:ie,js:je,ke+1,comp) = ZERO
       else
          sedgez(is:ie,js:je,ke+1,comp) = sedgelz(is:ie,js:je,ke+1)
       end if
    else if (phys_bc(3,2) .eq. OUTLET) then
       if (is_vel .and. comp .eq. 3) then
          sedgez(is:ie,js:je,ke+1,comp) = max(sedgelz(is:ie,js:je,ke+1),ZERO)
       else
          sedgez(is:ie,js:je,ke+1,comp) = sedgelz(is:ie,js:je,ke+1)
       end if
    end if

    deallocate(sedgelz,sedgerz)

  end subroutine make_edge_scal_3d

end module make_edge_scal_module
