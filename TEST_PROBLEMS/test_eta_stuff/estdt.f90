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

  implicit none

  private

  public :: estdt

contains

  subroutine estdt(n, u, s, force, divU, dSdt, normal, w0, p0, gamma1bar, dx, cflfac, dt)

    use bl_prof_module
    use geometry, only: spherical
    
    integer        , intent(in ) :: n
    type(multifab) , intent(in ) :: u
    type(multifab) , intent(in ) :: s
    type(multifab) , intent(in ) :: force
    type(multifab) , intent(in ) :: divU
    type(multifab) , intent(in ) :: dSdt
    type(multifab) , intent(in ) :: normal
    real(kind=dp_t), intent(in ) :: w0(0:), p0(0:), gamma1bar(0:)
    real(kind=dp_t), intent(in ) :: dx(:)
    real(kind=dp_t), intent(in ) :: cflfac
    real(kind=dp_t), intent(out) :: dt
    
    real(kind=dp_t), pointer:: uop(:,:,:,:)
    real(kind=dp_t), pointer:: sop(:,:,:,:)
    real(kind=dp_t), pointer:: fp(:,:,:,:)
    real(kind=dp_t), pointer:: nop(:,:,:,:)
    real(kind=dp_t), pointer:: dUp(:,:,:,:)
    real(kind=dp_t), pointer:: dSdtp(:,:,:,:)
    
    integer :: lo(mla%dim),hi(mla%dim),i,dm,nlevs
    integer :: ng_s,ng_u,ng_f,ng_dU,ng_dS,ng_n
    real(kind=dp_t) :: dt_adv,dt_adv_grid,dt_adv_proc,dt_start
    real(kind=dp_t) :: dt_divu,dt_divu_grid,dt_divu_proc
    
    real(kind=dp_t), parameter :: rho_min = 1.d-20

    type(bl_prof_timer), save :: bpt

    call build(bpt, "estdt")

    dm = mla%dim
    nlevs = mla%nlevel
    
    ng_u = u%ng
    ng_s = s%ng
    ng_f = force%ng
    ng_dU = divU%ng
    ng_dS = dSdt%ng

    dt_adv_proc   = HUGE(dt_adv_proc)
    dt_divu_proc  = HUGE(dt_divu_proc)
    dt_start      = HUGE(dt_start)
    
    do i = 1, u%nboxes
       if ( multifab_remote(u, i) ) cycle
       uop   => dataptr(u, i)
       sop   => dataptr(s, i)
       fp    => dataptr(force, i)
       dUp   => dataptr(divU, i)
       dSdtp => dataptr(dSdt, i)
       lo =  lwb(get_box(u, i))
       hi =  upb(get_box(u, i))

       dt_adv_grid   = HUGE(dt_adv_grid)
       dt_divu_grid  = HUGE(dt_divu_grid)

       select case (dm)
       case (2)
          call estdt_2d(n, uop(:,:,1,:), ng_u, sop(:,:,1,:), ng_s, &
                        fp(:,:,1,:), ng_f, dUp(:,:,1,1), ng_dU, &
                        dSdtp(:,:,1,1), ng_dS, &
                        w0, p0, gamma1bar, lo, hi, &
                        dx, rho_min, dt_adv_grid, dt_divu_grid, cflfac)
       case (3)
          if (spherical .eq. 1) then
             nop => dataptr(normal, i)
             ng_n = normal%ng
             call estdt_3d_sphr(n, uop(:,:,:,:), ng_u, sop(:,:,:,:), ng_s, &
                                fp(:,:,:,:), ng_f, dUp(:,:,:,1), ng_dU, &
                                dSdtp(:,:,:,1), ng_dS, nop(:,:,:,:), ng_n, &
                                w0, p0, gamma1bar, lo, hi, dx, &
                                rho_min, dt_adv_grid, dt_divu_grid, cflfac)
          else
             call estdt_3d_cart(n, uop(:,:,:,:), ng_u, sop(:,:,:,:), ng_s, &
                                fp(:,:,:,:), ng_f, dUp(:,:,:,1), ng_dU, &
                                dSdtp(:,:,:,1), ng_dS, &
                                w0, p0, gamma1bar, lo, hi, dx, rho_min, &
                                dt_adv_grid, dt_divu_grid, cflfac)
          end if
       end select
       
       dt_adv_proc  = min(dt_adv_proc ,dt_adv_grid)
       dt_divu_proc = min(dt_divu_proc,dt_divu_grid)

    end do
    
    ! This sets dt to be the min of dt_proc over all processors.
    call parallel_reduce(dt_adv ,dt_adv_proc ,MPI_MIN)
    call parallel_reduce(dt_divu,dt_divu_proc,MPI_MIN)
    
    dt = min(dt_adv,dt_divu)

    if (dt .eq. dt_start) then
       dt = min(dx(1),dx(2))
       if (dm .eq. 3) dt = min(dt,dx(3))
    end if

    dt = min(dt,0.05d0)

    call destroy(bpt)
    
  end subroutine estdt
  
  
  subroutine estdt_2d(n, u, ng_u, s, ng_s, force, ng_f, &
                      divU, ng_dU, dSdt, ng_dS, &
                      w0, p0, gamma1bar, lo, hi, &
                      dx, rho_min, dt_adv, dt_divu, cfl)

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
    real (kind = dp_t), intent(inout) :: dt_adv,dt_divu
    
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
    
    do j = lo(2), hi(2); do i = lo(1), hi(1)
       spdx = max(spdx ,abs(u(i,j,1)))
    enddo; enddo

    do j = lo(2), hi(2); do i = lo(1), hi(1)
       spdy = max(spdy ,abs(u(i,j,2)+w0(j)))
    enddo; enddo
    
    do j = lo(2),hi(2)
       spdr = max(spdr ,abs(w0(j)))
    enddo

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
                           dx, rho_min, dt_adv, dt_divu, cfl)

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
    real (kind = dp_t), intent(inout) :: dt_adv, dt_divu
    
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
    
    ! Limit dt based on velocity terms
    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       spdx = max(spdx ,abs(u(i,j,k,1)))
    enddo; enddo; enddo

    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       spdy = max(spdy ,abs(u(i,j,k,2)))
    enddo; enddo; enddo

    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       spdz = max(spdz ,abs(u(i,j,k,3)+w0(k)))
    enddo; enddo; enddo
    
    do k = lo(3),hi(3)
       spdr = max(spdr ,abs(w0(k)))
    enddo

    if (spdx > eps) dt_adv = min(dt_adv, dx(1)/spdx)
    if (spdy > eps) dt_adv = min(dt_adv, dx(2)/spdy)
    if (spdz > eps) dt_adv = min(dt_adv, dx(3)/spdz)
    if (spdr > eps) dt_adv = min(dt_adv, dx(3)/spdr)

    dt_adv = dt_adv * cfl
    
    ! Limit dt based on forcing terms
    fx = ZERO 
    fy = ZERO 
    fz = ZERO 
    
    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       fx = max(fx,abs(force(i,j,k,1)))
    enddo; enddo; enddo

    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       fy = max(fy,abs(force(i,j,k,2)))
    enddo; enddo; enddo

    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       fz = max(fz,abs(force(i,j,k,3)))
    enddo; enddo; enddo
    
    if (fx > eps) &
       dt_adv = min(dt_adv,sqrt(2.0D0*dx(1)/fx))
    
    if (fy > eps) &
       dt_adv = min(dt_adv,sqrt(2.0D0*dx(2)/fy))
    
    if (fz > eps) &
       dt_adv = min(dt_adv,sqrt(2.0D0*dx(3)/fz))
    
    ! divU constraint
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
  
  subroutine estdt_3d_sphr(n, u, ng_u, s, ng_s, force, ng_f, &
                           divU, ng_dU, dSdt, ng_dS, normal, ng_n, &
                           w0, p0, gamma1bar, &
                           lo, hi, dx, rho_min, dt_adv, dt_divu, cfl)

    use geometry,  only: dr, nr_fine
    use variables, only: rho_comp
    use fill_3d_module
    
    integer, intent(in) :: n, lo(:), hi(:), ng_u, ng_s, ng_f, ng_dU, ng_dS, ng_n
    real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :,:)  
    real (kind = dp_t), intent(in   ) ::      s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)  
    real (kind = dp_t), intent(in   ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)  
    real (kind = dp_t), intent(in   ) ::   divU(lo(1)-ng_dU:,lo(2)-ng_dU:,lo(3)-ng_dU:)
    real (kind = dp_t), intent(in   ) ::   dSdt(lo(1)-ng_dS:,lo(2)-ng_dS:,lo(3)-ng_dS:)
    real (kind = dp_t), intent(in   ) :: normal(lo(1)-ng_n :,lo(2)-ng_n :,lo(3)-ng_n :,:)
    real (kind = dp_t), intent(in   ) :: w0(0:), p0(0:), gamma1bar(0:)
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(in   ) :: rho_min, cfl
    real (kind = dp_t), intent(inout) :: dt_adv, dt_divu
    
    real (kind = dp_t), allocatable ::  w0_cart(:,:,:,:)
    real (kind = dp_t), allocatable :: gp0_cart(:,:,:,:)
    real (kind = dp_t), allocatable :: gp0(:)
    real (kind = dp_t)  :: spdx, spdy, spdz, spdr, gp_dot_u, gamma1bar_p_avg
    real (kind = dp_t)  :: fx, fy, fz, eps, denom, a, b, c
    integer             :: i,j,k,r
    
    eps = 1.0d-8
    
    spdx = ZERO
    spdy = ZERO 
    spdz = ZERO 
    spdr = ZERO 
    
    allocate( w0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3))
    call put_1d_array_on_cart_3d_sphr(.true.,.true.,w0,w0_cart,lo,hi,dx,0,ng_n,normal)

    ! Limit dt based on velocity terms
    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       spdx = max(spdx ,abs(u(i,j,k,1)+w0_cart(i,j,k,1)))
    enddo; enddo; enddo

    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       spdy = max(spdy ,abs(u(i,j,k,2)+w0_cart(i,j,k,2)))
    enddo; enddo; enddo

    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       spdz = max(spdz ,abs(u(i,j,k,3)+w0_cart(i,j,k,3)))
    enddo; enddo; enddo
    
    deallocate(w0_cart)

    do k = 0,size(w0,dim=1)-1
       spdr = max(spdr ,abs(w0(k)))
    enddo

    if (spdx > eps) dt_adv = min(dt_adv, dx(1)/spdx)
    if (spdy > eps) dt_adv = min(dt_adv, dx(2)/spdy)
    if (spdz > eps) dt_adv = min(dt_adv, dx(3)/spdz)
    if (spdr > eps) dt_adv = min(dt_adv, dr(n)/spdr)

    dt_adv = dt_adv * cfl
    
    ! Limit dt based on forcing terms
    fx = ZERO 
    fy = ZERO 
    fz = ZERO 
    
    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       fx = max(fx,abs(force(i,j,k,1)))
    enddo; enddo; enddo

    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       fy = max(fy,abs(force(i,j,k,2)))
    enddo; enddo; enddo

    do k = lo(3), hi(3); do j = lo(2), hi(2); do i = lo(1), hi(1)
       fz = max(fz,abs(force(i,j,k,3)))
    enddo; enddo; enddo
    
    if (fx > eps) &
       dt_adv = min(dt_adv,sqrt(2.0D0*dx(1)/fx))
    
    if (fy > eps) &
       dt_adv = min(dt_adv,sqrt(2.0D0*dx(2)/fy))
    
    if (fz > eps) &
       dt_adv = min(dt_adv,sqrt(2.0D0*dx(3)/fz))
    
    ! divU constraint
    allocate(gp0(0:nr_fine))
    do r=1,nr_fine-1
       gamma1bar_p_avg = HALF * (gamma1bar(r)*p0(r) + gamma1bar(r-1)*p0(r-1))
       gp0(r) = ( (p0(r) - p0(r-1))/dr(n) ) / gamma1bar_p_avg
    end do
    gp0(nr_fine) = gp0(nr_fine-1)
    gp0(      0) = gp0(        1)
    allocate(gp0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3))
    call put_1d_array_on_cart_3d_sphr(.true.,.true.,gp0,gp0_cart,lo,hi,dx,0,ng_n,normal)
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             
             gp_dot_u = u(i,j,k,1) * gp0_cart(i,j,k,1) + &
                        u(i,j,k,2) * gp0_cart(i,j,k,2) + &
                        u(i,j,k,3) * gp0_cart(i,j,k,3)
             
             denom = divU(i,j,k) - gp_dot_u 
             
             if (denom > ZERO) then
                dt_divu = &
                     min(dt_divu,0.4d0*(ONE - rho_min/s(i,j,k,rho_comp))/denom)
             endif
             
          enddo
       enddo
    enddo
    
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

  end subroutine estdt_3d_sphr
  
end module estdt_module
