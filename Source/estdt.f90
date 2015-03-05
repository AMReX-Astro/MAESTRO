! Compute the timestep.  Here we use several different estimations and take
! the most restrictive one.  They are (see paper III): 
!    1. the advective constraint
!    2. the force constraint
!    3. a div{U} constraint
!    4. a dS/dt constraint
!
! the differences between firstdt and estdt are as follows:
!   firstdt does not use the base state velocity w0 since it's supposed to be 0
!   firstdt uses the sound speed time step constraint if the velocity is 0
!

module estdt_module
  
  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: estdt

contains

  subroutine estdt(mla,the_bc_tower,u,s,gpi,divU,dSdt,w0,rho0,p0,gamma1bar, &
                   grav,dx,cflfac,dt)

    use bl_prof_module
    use geometry, only: spherical, nlevs_radial, nr_fine
    use variables, only: rel_eps, rho_comp
    use bl_constants_module
    use mk_vel_force_module
    use fill_3d_module
    use probin_module, only: evolve_base_state, small_dt

    type(ml_layout), intent(inout) :: mla
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: gpi(:)
    type(multifab) , intent(in   ) :: divU(:)
    type(multifab) , intent(in   ) :: dSdt(:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: grav(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: cflfac
    real(kind=dp_t), intent(inout) :: dt
    
    type(multifab) ::      force(mla%nlevel)
    type(multifab) ::      w0mac(mla%nlevel,mla%dim)
    type(multifab) :: umac_dummy(mla%nlevel,mla%dim)
    type(multifab) :: w0_force_cart_dummy(mla%nlevel)
    type(multifab) :: normal_dummy(mla%nlevel)

    logical :: is_final_update

    real(kind=dp_t), allocatable :: w0_force_dummy(:,:)

    real(kind=dp_t), pointer::   uop(:,:,:,:)
    real(kind=dp_t), pointer::   sop(:,:,:,:)
    real(kind=dp_t), pointer::    fp(:,:,:,:)
    real(kind=dp_t), pointer::   dUp(:,:,:,:)
    real(kind=dp_t), pointer:: dSdtp(:,:,:,:)
    real(kind=dp_t), pointer::   wxp(:,:,:,:)
    real(kind=dp_t), pointer::   wyp(:,:,:,:)
    real(kind=dp_t), pointer::   wzp(:,:,:,:)
    
    integer :: lo(mla%dim),hi(mla%dim),i,n,comp,dm,nlevs
    integer :: ng_s,ng_u,ng_f,ng_dU,ng_dS,ng_w
    real(kind=dp_t) :: dt_adv,dt_adv_grid,dt_adv_proc,dt_start,dt_lev
    real(kind=dp_t) :: dt_divu,dt_divu_grid,dt_divu_proc
    real(kind=dp_t) :: umax,umax_grid,umax_proc,umax_lev
    real(kind=dp_t) :: dts_local(2), dts_global(2)
    
    real(kind=dp_t), parameter :: rho_min = 1.d-20

    type(bl_prof_timer), save :: bpt

    call build(bpt, "estdt")

    dm = mla%dim
    nlevs = mla%nlevel
    
    allocate(w0_force_dummy(nlevs_radial,0:nr_fine-1))
    w0_force_dummy = 0.d0

    umax = 0.d0

    do n=1,nlevs
       call multifab_build(force(n), mla%la(n), dm, 1)
       call multifab_build(w0_force_cart_dummy(n),mla%la(n),dm,1)
       call setval(w0_force_cart_dummy(n),0.d0,all=.true.)

       ! create an empty umac so we can call the vel force routine --
       ! this will not be used
       do comp=1,dm
          call multifab_build_edge(umac_dummy(n,comp), mla%la(n),1,1,comp)
          call setval(umac_dummy(n,comp), ZERO, all=.true.)
       end do
 
       ! create a dummy normal also for the vel routine -- this is not
       ! used if we don't add the utilde force
       call multifab_build(normal_dummy(n), mla%la(n), dm, 1)

    end do

    if (spherical .eq. 1) then

       do n=1,nlevs
          do comp=1,dm
             call multifab_build_edge(w0mac(n,comp),mla%la(n),1,1,comp)
             call setval(w0mac(n,comp), ZERO, all=.true.)
          end do
       end do

       if (evolve_base_state) then
          call make_w0mac(mla,w0,w0mac,dx,the_bc_tower%bc_tower_array)
       end if

    end if
    
    is_final_update = .false.
    call mk_vel_force(force,is_final_update, &
                      u,umac_dummy,w0,w0mac,gpi,s,rho_comp,normal_dummy, &
                      rho0,grav,dx,w0_force_dummy,w0_force_cart_dummy, &
                      the_bc_tower%bc_tower_array,mla,.false.)

    do n=1,nlevs
       call destroy(w0_force_cart_dummy(n))
       do comp=1,dm
          call destroy(umac_dummy(n,comp))
       end do
       call destroy(normal_dummy(n))
    end do

    ng_u  = nghost(u(1))
    ng_s  = nghost(s(1))
    ng_f  = nghost(force(1))
    ng_dU = nghost(divU(1))
    ng_dS = nghost(dSdt(1))

    do n=1,nlevs
       
       dt_adv_proc  = 1.d99
       dt_divu_proc = 1.d99
       dt_start     = 1.d99
       umax_grid    = 0.d0
       umax_proc    = 0.d0

       do i = 1, nfabs(u(n))
          uop   => dataptr(u(n), i)
          sop   => dataptr(s(n), i)
          fp    => dataptr(force(n), i)
          dUp   => dataptr(divU(n), i)
          dSdtp => dataptr(dSdt(n), i)
          lo =  lwb(get_box(u(n), i))
          hi =  upb(get_box(u(n), i))

          dt_adv_grid   = HUGE(dt_adv_grid)
          dt_divu_grid  = HUGE(dt_divu_grid)

          select case (dm)
          case (1)
             call estdt_1d(n, uop(:,1,1,1), ng_u, sop(:,1,1,:), ng_s, &
                           fp(:,1,1,1), ng_f, dUp(:,1,1,1), ng_dU, &
                           dSdtp(:,1,1,1), ng_dS, &
                           w0(n,:), p0(n,:), gamma1bar(n,:), lo, hi, &
                           dx(n,:), rho_min, dt_adv_grid, dt_divu_grid, umax_grid, cflfac)
          case (2)
             call estdt_2d(n, uop(:,:,1,:), ng_u, sop(:,:,1,:), ng_s, &
                           fp(:,:,1,:), ng_f, dUp(:,:,1,1), ng_dU, &
                           dSdtp(:,:,1,1), ng_dS, &
                           w0(n,:), p0(n,:), gamma1bar(n,:), lo, hi, &
                           dx(n,:), rho_min, dt_adv_grid, dt_divu_grid, umax_grid, cflfac)
          case (3)
             if (spherical .eq. 1) then
                ng_w = nghost(w0mac(1,1))
                wxp => dataptr(w0mac(n,1), i)
                wyp => dataptr(w0mac(n,2), i)
                wzp => dataptr(w0mac(n,3), i)
                call estdt_3d_sphr(uop(:,:,:,:), ng_u, sop(:,:,:,:), ng_s, &
                                   fp(:,:,:,:), ng_f, dUp(:,:,:,1), ng_dU, &
                                   dSdtp(:,:,:,1), ng_dS, &
                                   w0(1,:),wxp(:,:,:,1),wyp(:,:,:,1),wzp(:,:,:,1),ng_w, &
                                   p0(1,:), gamma1bar(1,:), lo, hi, dx(n,:), &
                                   rho_min, dt_adv_grid, dt_divu_grid, umax_grid, cflfac)
             else
                call estdt_3d_cart(n, uop(:,:,:,:), ng_u, sop(:,:,:,:), ng_s, &
                                   fp(:,:,:,:), ng_f, dUp(:,:,:,1), ng_dU, &
                                   dSdtp(:,:,:,1), ng_dS, &
                                   w0(n,:), p0(n,:), gamma1bar(n,:), lo, hi, dx(n,:), &
                                   rho_min, dt_adv_grid, dt_divu_grid, umax_grid, cflfac)
             end if
          end select

          dt_adv_proc  = min( dt_adv_proc,  dt_adv_grid)
          dt_divu_proc = min(dt_divu_proc, dt_divu_grid)
          umax_proc    = max(   umax_proc,    umax_grid)

       end do

       ! This sets dt to be the min of dt_proc over all processors.
       dts_local(1) = dt_adv_proc
       dts_local(2) = dt_divu_proc
       call parallel_reduce( dts_global,  dts_local, MPI_MIN)
       dt_adv = dts_global(1)
       dt_divu = dts_global(2)

       ! we could pack this in the MIN reduce by looking at 1/umax_proc
       call parallel_reduce(umax_lev,    umax_proc, MPI_MAX)

       umax = max(umax,umax_lev)

       dt_lev = min(dt_adv,dt_divu)

       if (dt_lev .eq. dt_start) then
          if (dm .eq. 1) then 
             dt_lev = dx(n,1)
          else
             dt_lev = min(dx(n,1),dx(n,2))
             if (dm .eq. 3) dt_lev = min(dt_lev,dx(n,3))
          end if
       end if

       dt = min(dt,dt_lev)

    end do   ! end loop over levels

    if (dt < small_dt) then
       call bl_error("ERROR: timestep < small_dt")
    endif

    rel_eps = 1.d-8*umax

     do n=1,nlevs
        call destroy(force(n))
     end do

     if (spherical .eq. 1) then
        do n=1,nlevs
           do comp=1,dm
              call destroy(w0mac(n,comp))
           end do
        end do
     end if

    deallocate(w0_force_dummy)

    call destroy(bpt)

  end subroutine estdt

  subroutine estdt_1d(n, u, ng_u, s, ng_s, force, ng_f, &
                      divU, ng_dU, dSdt, ng_dS, &
                      w0, p0, gamma1bar, lo, hi, &
                      dx, rho_min, dt_adv, dt_divu, umax, cfl)

    use geometry,  only: nr
    use variables, only: rho_comp

    integer, intent(in) :: n, lo(:), hi(:), ng_u, ng_s, ng_f, ng_dU, ng_dS
    real (kind = dp_t), intent(in   ) ::     u(lo(1)-ng_u :)  
    real (kind = dp_t), intent(in   ) ::     s(lo(1)-ng_s :,:)  
    real (kind = dp_t), intent(in   ) :: force(lo(1)-ng_f :)  
    real (kind = dp_t), intent(in   ) ::  divU(lo(1)-ng_dU:)
    real (kind = dp_t), intent(in   ) ::  dSdt(lo(1)-ng_dS:)
    real (kind = dp_t), intent(in   ) :: w0(0:), p0(0:), gamma1bar(0:)
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(in   ) :: rho_min,cfl
    real (kind = dp_t), intent(inout) :: dt_adv,dt_divu,umax
    
    real (kind = dp_t)  :: spdx, spdr
    real (kind = dp_t)  :: fx
    real (kind = dp_t)  :: eps
    real (kind = dp_t)  :: denom, gradp0
    real (kind = dp_t)  :: a, b, c
    integer             :: i
    
    eps = 1.0d-8
    
    ! advective constraints
    spdx = ZERO
    spdr = ZERO
    umax = ZERO
    
    do i = lo(1), hi(1)
       spdx = max(spdx ,abs(u(i)))
       spdr = max(spdr ,abs(w0(i)))
    enddo

    umax = max(spdr,spdx)

    if (spdx > eps) dt_adv = min(dt_adv, dx(1)/spdx)
    if (spdr > eps) dt_adv = min(dt_adv, dx(1)/spdr)

    dt_adv = dt_adv * cfl
    
    ! force constraints
    fx = ZERO
    
    do i = lo(1), hi(1)
       fx = max(fx,abs(force(i)))
    enddo
    
    if (fx > eps) &
       dt_adv = min(dt_adv,sqrt(2.0D0 *dx(1)/fx))
    
    ! divU constraint
    do i = lo(1), hi(1)
       
       if (i .eq. 0) then
          gradp0 = (p0(i+1) - p0(i))/dx(1)
       else if (i .eq. nr(n)-1) then
          gradp0 = (p0(i) - p0(i-1))/dx(1)
       else
          gradp0 = HALF*(p0(i+1) - p0(i-1))/dx(1)
       endif
          
       denom = divU(i) - u(i)*gradp0/(gamma1bar(i)*p0(i))
          
       if (denom > ZERO) then
          
          dt_divu = min(dt_divu, &
               0.4d0*(ONE - rho_min/s(i,rho_comp))/denom)
       endif
          
    enddo
    
    do i = lo(1), hi(1)           
          !
          ! An additional dS/dt timestep constraint originally
          ! used in nova
          ! solve the quadratic equation
          ! (rho - rho_min)/(rho dt) = S + (dt/2)*(dS/dt)
          ! which is equivalent to
          ! (rho/2)*dS/dt*dt^2 + rho*S*dt + (rho_min-rho) = 0
          ! which has solution dt = 2.0d0*c/(-b-sqrt(b**2-4.0d0*a*c))
          !
          if (dSdt(i) .gt. 1.d-20) then
             a = HALF*s(i,rho_comp)*dSdt(i)
             b = s(i,rho_comp)*divU(i)
             c = rho_min - s(i,rho_comp)
             dt_divu = min(dt_divu, 0.4d0*2.0d0*c/(-b-sqrt(b**2-4.0d0*a*c)))
          endif
          
    enddo
    
  end subroutine estdt_1d
  
  subroutine estdt_2d(n, u, ng_u, s, ng_s, force, ng_f, &
                      divU, ng_dU, dSdt, ng_dS, &
                      w0, p0, gamma1bar, lo, hi, &
                      dx, rho_min, dt_adv, dt_divu, umax, cfl)

    use geometry,  only: nr
    use variables, only: rho_comp

    integer, intent(in) :: n, lo(:), hi(:), ng_u, ng_s, ng_f, ng_dU, ng_dS
    real (kind = dp_t), intent(in   ) ::     u(lo(1)-ng_u :,lo(2)-ng_u :,:)  
    real (kind = dp_t), intent(in   ) ::     s(lo(1)-ng_s :,lo(2)-ng_s :,:)  
    real (kind = dp_t), intent(in   ) :: force(lo(1)-ng_f :,lo(2)-ng_f :,:)  
    real (kind = dp_t), intent(in   ) ::  divU(lo(1)-ng_dU:,lo(2)-ng_dU:)
    real (kind = dp_t), intent(in   ) ::  dSdt(lo(1)-ng_dS:,lo(2)-ng_dS:)
    real (kind = dp_t), intent(in   ) :: w0(0:), p0(0:), gamma1bar(0:)
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(in   ) :: rho_min,cfl
    real (kind = dp_t), intent(inout) :: dt_adv,dt_divu,umax
    
    real (kind = dp_t)  :: spdx, spdy, spdr
    real (kind = dp_t)  :: fx, fy
    real (kind = dp_t)  :: eps
    real (kind = dp_t)  :: denom, gradp0
    real (kind = dp_t)  :: a, b, c
    integer             :: i,j
    
    eps = 1.0d-8
    
    ! advective constraints
    spdx = ZERO
    spdy = ZERO
    spdr = ZERO
    umax = ZERO
    
    do j = lo(2), hi(2); do i = lo(1), hi(1)
       spdx = max(spdx ,abs(u(i,j,1)))
    enddo; enddo

    do j = lo(2), hi(2); do i = lo(1), hi(1)
       spdy = max(spdy ,abs(u(i,j,2)+HALF*(w0(j)+w0(j+1))))
    enddo; enddo
    
    do j = lo(2),hi(2)
       spdr = max(spdr ,abs(w0(j)))
    enddo

    umax = max(umax,spdx,spdy,spdr)

    if (spdx > eps) dt_adv = min(dt_adv, dx(1)/spdx)
    if (spdy > eps) dt_adv = min(dt_adv, dx(2)/spdy)
    if (spdr > eps) dt_adv = min(dt_adv, dx(2)/spdr)

    dt_adv = dt_adv * cfl
    
    ! force constraints
    fx = ZERO
    fy = ZERO
    
    do j = lo(2), hi(2); do i = lo(1), hi(1)
       fx = max(fx,abs(force(i,j,1)))
    enddo; enddo

    do j = lo(2), hi(2); do i = lo(1), hi(1)
       fy = max(fy,abs(force(i,j,2)))
    enddo; enddo
    
    if (fx > eps) &
       dt_adv = min(dt_adv,sqrt(2.0D0 *dx(1)/fx))
    
    if (fy > eps) &
       dt_adv = min(dt_adv,sqrt(2.0D0 *dx(2)/fy))
    
    ! divU constraint
    do j = lo(2), hi(2)
       
       if (j .eq. 0) then
          gradp0 = (p0(j+1) - p0(j))/dx(2)
       else if (j .eq. nr(n)-1) then
          gradp0 = (p0(j) - p0(j-1))/dx(2)
       else
          gradp0 = HALF*(p0(j+1) - p0(j-1))/dx(2)
       endif
       
       do i = lo(1), hi(1)
          
          denom = divU(i,j) - u(i,j,2)*gradp0/(gamma1bar(j)*p0(j))
          
          if (denom > ZERO) then
             
             dt_divu = min(dt_divu, &
                  0.4d0*(ONE - rho_min/s(i,j,rho_comp))/denom)
          endif
          
       enddo
    enddo
    
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)           
          !
          ! An additional dS/dt timestep constraint originally
          ! used in nova
          ! solve the quadratic equation
          ! (rho - rho_min)/(rho dt) = S + (dt/2)*(dS/dt)
          ! which is equivalent to
          ! (rho/2)*dS/dt*dt^2 + rho*S*dt + (rho_min-rho) = 0
          ! which has solution dt = 2.0d0*c/(-b-sqrt(b**2-4.0d0*a*c))
          !
          if (dSdt(i,j) .gt. 1.d-20) then
             a = HALF*s(i,j,rho_comp)*dSdt(i,j)
             b = s(i,j,rho_comp)*divU(i,j)
             c = rho_min - s(i,j,rho_comp)
             dt_divu = min(dt_divu, 0.4d0*2.0d0*c/(-b-sqrt(b**2-4.0d0*a*c)))
          endif
          
       enddo
    enddo
    
  end subroutine estdt_2d
  
  subroutine estdt_3d_cart(n, u, ng_u, s, ng_s, force, ng_f, &
                           divU, ng_dU, dSdt, ng_dS, &
                           w0, p0, gamma1bar, lo, hi, &
                           dx, rho_min, dt_adv, dt_divu, umax, cfl)

    use geometry,  only: nr
    use variables, only: rho_comp

    integer, intent(in) :: n, lo(:), hi(:), ng_u, ng_s, ng_f, ng_dU, ng_dS
    real (kind = dp_t), intent(in   ) ::     u(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :,:)  
    real (kind = dp_t), intent(in   ) ::     s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)  
    real (kind = dp_t), intent(in   ) :: force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    real (kind = dp_t), intent(in   ) ::  divU(lo(1)-ng_dU:,lo(2)-ng_dU:,lo(3)-ng_dU:)
    real (kind = dp_t), intent(in   ) ::  dSdt(lo(1)-ng_dS:,lo(2)-ng_dS:,lo(3)-ng_dS:)
    real (kind = dp_t), intent(in   ) :: w0(0:), p0(0:), gamma1bar(0:)
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(in   ) :: rho_min,cfl
    real (kind = dp_t), intent(inout) :: dt_adv, dt_divu, umax
    
    real (kind = dp_t)  :: spdx, spdy, spdz, spdr
    real (kind = dp_t)  :: fx, fy, fz
    real (kind = dp_t)  :: eps,denom,gradp0
    real (kind = dp_t)  :: a, b, c
    integer             :: i,j,k
    
    eps = 1.0d-8
    
    spdx = ZERO
    spdy = ZERO 
    spdz = ZERO 
    spdr = ZERO
    umax = ZERO
    !
    ! Limit dt based on velocity terms.
    !
    !$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(MAX : spdx)
    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       spdx = max(spdx ,abs(u(i,j,k,1)))
    enddo; enddo; enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(MAX : spdy)
    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       spdy = max(spdy ,abs(u(i,j,k,2)))
    enddo; enddo; enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(MAX : spdz)
    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       spdz = max(spdz ,abs(u(i,j,k,3)+HALF*(w0(k)+w0(k+1))))
    enddo; enddo; enddo
    !$OMP END PARALLEL DO
    
    do k = lo(3),hi(3)
       spdr = max(spdr ,abs(w0(k)))
    enddo

    umax = max(umax,spdx,spdy,spdz,spdr)

    if (spdx > eps) dt_adv = min(dt_adv, dx(1)/spdx)
    if (spdy > eps) dt_adv = min(dt_adv, dx(2)/spdy)
    if (spdz > eps) dt_adv = min(dt_adv, dx(3)/spdz)
    if (spdr > eps) dt_adv = min(dt_adv, dx(3)/spdr)

    dt_adv = dt_adv * cfl
    !
    ! Limit dt based on forcing terms
    !
    fx = ZERO 
    fy = ZERO 
    fz = ZERO 

    !$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(MAX : fx)
    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       fx = max(fx,abs(force(i,j,k,1)))
    enddo; enddo; enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(MAX : fy)
    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       fy = max(fy,abs(force(i,j,k,2)))
    enddo; enddo; enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(MAX : fz)
    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       fz = max(fz,abs(force(i,j,k,3)))
    enddo; enddo; enddo
    !$OMP END PARALLEL DO
    
    if (fx > eps) &
       dt_adv = min(dt_adv,sqrt(2.0D0*dx(1)/fx))
    
    if (fy > eps) &
       dt_adv = min(dt_adv,sqrt(2.0D0*dx(2)/fy))
    
    if (fz > eps) &
       dt_adv = min(dt_adv,sqrt(2.0D0*dx(3)/fz))
    !
    ! divU constraint
    !
    !$OMP PARALLEL DO PRIVATE(i,j,k,gradp0,denom) REDUCTION(MIN : dt_divu)
    do k = lo(3), hi(3)
       
       if (k .eq. 0) then
          gradp0 = (p0(k+1) - p0(k))/dx(3)
       else if (k .eq. nr(n)-1) then
          gradp0 = (p0(k) - p0(k-1))/dx(3)
       else
          gradp0 = HALF*(p0(k+1) - p0(k-1))/dx(3)
       endif
       
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             
             denom = divU(i,j,k) - u(i,j,k,3)*gradp0/(gamma1bar(k)*p0(k))
             
             if (denom > ZERO) then
                dt_divu = min(dt_divu, &
                     0.4d0*(ONE - rho_min/s(i,j,k,rho_comp))/denom)
             endif
             
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO
    
    !$OMP PARALLEL DO PRIVATE(i,j,k,a,b,c) REDUCTION(MIN : dt_divu)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)           
             !
             ! An additional dS/dt timestep constraint originally
             ! used in nova
             ! solve the quadratic equation
             ! (rho - rho_min)/(rho dt) = S + (dt/2)*(dS/dt)
             ! which is equivalent to
             ! (rho/2)*dS/dt*dt^2 + rho*S*dt + (rho_min-rho) = 0
             ! which has solution dt = 2.0d0*c/(-b-sqrt(b**2-4.0d0*a*c))
             !
             if (dSdt(i,j,k) .gt. 1.d-20) then
                a = HALF*s(i,j,k,rho_comp)*dSdt(i,j,k)
                b = s(i,j,k,rho_comp)*divU(i,j,k)
                c = rho_min - s(i,j,k,rho_comp)
                dt_divu = min(dt_divu,0.4d0*2.0d0*c/(-b-sqrt(b**2-4.0d0*a*c)))
             endif
             
          enddo
       enddo
    enddo
    
  end subroutine estdt_3d_cart
  
  subroutine estdt_3d_sphr(u, ng_u, s, ng_s, force, ng_f, &
                           divU, ng_dU, dSdt, ng_dS, &
                           w0,w0macx,w0macy,w0macz,ng_w, p0, gamma1bar, &
                           lo, hi, dx, rho_min, dt_adv, dt_divu, umax, cfl)

    use geometry,  only: dr, nr_fine
    use variables, only: rho_comp
    use fill_3d_module
    
    integer           , intent(in   ) :: lo(:),hi(:),ng_u,ng_s,ng_f,ng_dU,ng_dS,ng_w
    real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :,:)  
    real (kind = dp_t), intent(in   ) ::      s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)  
    real (kind = dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)  
    real (kind = dp_t), intent(in   ) ::   divU(lo(1)-ng_dU:,lo(2)-ng_dU:,lo(3)-ng_dU:)
    real (kind = dp_t), intent(in   ) ::   dSdt(lo(1)-ng_dS:,lo(2)-ng_dS:,lo(3)-ng_dS:)
    real (kind = dp_t), intent(in   ) :: w0macx(lo(1)-ng_w :,lo(2)-ng_w :,lo(3)-ng_w :)
    real (kind = dp_t), intent(in   ) :: w0macy(lo(1)-ng_w :,lo(2)-ng_w :,lo(3)-ng_w :)
    real (kind = dp_t), intent(in   ) :: w0macz(lo(1)-ng_w :,lo(2)-ng_w :,lo(3)-ng_w :)
    real (kind = dp_t), intent(in   ) :: w0(0:), p0(0:), gamma1bar(0:)
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(in   ) :: rho_min, cfl
    real (kind = dp_t), intent(inout) :: dt_adv, dt_divu, umax
    
    real (kind = dp_t), allocatable :: gp0_cart(:,:,:,:)

    real (kind = dp_t) :: gp0(0:nr_fine)

    real (kind = dp_t) :: spdx, spdy, spdz, spdr, gp_dot_u, gamma1bar_p_avg
    real (kind = dp_t) :: fx, fy, fz, eps, denom, a, b, c
    integer            :: i,j,k,r

    allocate(gp0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3))
    
    eps = 1.0d-8
    
    spdx = ZERO
    spdy = ZERO 
    spdz = ZERO 
    spdr = ZERO
    umax = ZERO
    !
    ! Limit dt based on velocity terms
    !
    !$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(MAX : spdx)
    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       spdx = max(spdx ,abs(u(i,j,k,1)+HALF*(w0macx(i,j,k)+w0macx(i+1,j,k))))
    enddo; enddo; enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(MAX : spdy)
    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       spdy = max(spdy ,abs(u(i,j,k,2)+HALF*(w0macy(i,j,k)+w0macy(i,j+1,k))))
    enddo; enddo; enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(MAX : spdz)
    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       spdz = max(spdz ,abs(u(i,j,k,3)+HALF*(w0macz(i,j,k)+w0macz(i,j,k+1))))
    enddo; enddo; enddo
    !$OMP END PARALLEL DO
    
    do k=0,nr_fine
       spdr = max(spdr ,abs(w0(k)))
    enddo

    umax = max(umax,spdx,spdy,spdz,spdr)

    if (spdx > eps) dt_adv = min(dt_adv, dx(1)/spdx)
    if (spdy > eps) dt_adv = min(dt_adv, dx(2)/spdy)
    if (spdz > eps) dt_adv = min(dt_adv, dx(3)/spdz)
    if (spdr > eps) dt_adv = min(dt_adv, dr(1)/spdr)

    dt_adv = dt_adv * cfl
    
    ! Limit dt based on forcing terms
    fx = ZERO 
    fy = ZERO 
    fz = ZERO 

    !$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(MAX : fx)
    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       fx = max(fx,abs(force(i,j,k,1)))
    enddo; enddo; enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(MAX : fy)
    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       fy = max(fy,abs(force(i,j,k,2)))
    enddo; enddo; enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(MAX : fz)
    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       fz = max(fz,abs(force(i,j,k,3)))
    enddo; enddo; enddo
    !$OMP END PARALLEL DO
    
    if (fx > eps) &
       dt_adv = min(dt_adv,sqrt(2.0D0*dx(1)/fx))
    
    if (fy > eps) &
       dt_adv = min(dt_adv,sqrt(2.0D0*dx(2)/fy))
    
    if (fz > eps) &
       dt_adv = min(dt_adv,sqrt(2.0D0*dx(3)/fz))
    
    ! divU constraint
    do r=1,nr_fine-1
       gamma1bar_p_avg = HALF * (gamma1bar(r)*p0(r) + gamma1bar(r-1)*p0(r-1))
       gp0(r) = ( (p0(r) - p0(r-1))/dr(1) ) / gamma1bar_p_avg
    end do
    gp0(nr_fine) = gp0(nr_fine-1)
    gp0(      0) = gp0(        1)

    call put_1d_array_on_cart_3d_sphr(.true.,.true.,gp0,gp0_cart,lo,hi,dx,0)
    
    !$OMP PARALLEL DO PRIVATE(i,j,k,a,b,c,gp_dot_u,denom) REDUCTION(MIN : dt_divu)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             
             gp_dot_u = u(i,j,k,1) * gp0_cart(i,j,k,1) + &
                        u(i,j,k,2) * gp0_cart(i,j,k,2) + &
                        u(i,j,k,3) * gp0_cart(i,j,k,3)
             
             denom = divU(i,j,k) - gp_dot_u 
             
             if (denom > ZERO) then
                dt_divu = min(dt_divu,0.4d0*(ONE - rho_min/s(i,j,k,rho_comp))/denom)
             endif
             !
             ! An additional dS/dt timestep constraint originally
             ! used in nova
             ! solve the quadratic equation
             ! (rho - rho_min)/(rho dt) = S + (dt/2)*(dS/dt)
             ! which is equivalent to
             ! (rho/2)*dS/dt*dt^2 + rho*S*dt + (rho_min-rho) = 0
             ! which has solution dt = 2.0d0*c/(-b-sqrt(b**2-4.0d0*a*c))
             !
             if (dSdt(i,j,k) .gt. 1.d-20) then
                a = HALF*s(i,j,k,rho_comp)*dSdt(i,j,k)
                b = s(i,j,k,rho_comp)*divU(i,j,k)
                c = rho_min - s(i,j,k,rho_comp)
                dt_divu = min(dt_divu,0.4d0*2.0d0*c/(-b-sqrt(b**2-4.0d0*a*c)))
             endif
                          
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(gp0_cart)

  end subroutine estdt_3d_sphr
  
end module estdt_module
