! Compute the initial timestep
!
! The differences between firstdt and estdt are as follows:
!   firstdt does not use the base state velocity w0 since it's supposed to be 0
!   firstdt uses the sound speed time step constraint if the velocity is 0
!
! After the initial projection, we should have a source term that gives
! rise to a velocity field, and we can use the normal estdt.

module firstdt_module

  use bl_types
  use multifab_module

  implicit none

  private

  public :: firstdt

contains

  subroutine firstdt(n,u,s,force,divU,p0,gamma1,dx,cflfac,dt)

    integer        , intent(in   ) :: n
    type(multifab) , intent(in   ) :: u,s,force,divU
    real(kind=dp_t), intent(in   ) :: p0(0:), cflfac, gamma1(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(  out) :: dt
    
    real(kind=dp_t), pointer:: uop(:,:,:,:)
    real(kind=dp_t), pointer:: sop(:,:,:,:)
    real(kind=dp_t), pointer:: fp(:,:,:,:)
    real(kind=dp_t), pointer:: divup(:,:,:,:)
    integer :: lo(u%dim),hi(u%dim),ng,dm,i
    real(kind=dp_t) :: dt_hold_proc,dt_grid
    
    ng = u%ng
    dm = u%dim
    
    dt_hold_proc = HUGE(dt_hold_proc)
    dt_grid      = HUGE(dt_grid)
    
    do i = 1, u%nboxes
       if ( multifab_remote(u, i) ) cycle
       uop   => dataptr(u, i)
       sop   => dataptr(s, i)
       fp    => dataptr(force, i)
       divup => dataptr(divU,i)
       lo    =  lwb(get_box(u, i))
       hi    =  upb(get_box(u, i))
       select case (dm)
       case (2)
          call firstdt_2d(n,uop(:,:,1,:), sop(:,:,1,:), fp(:,:,1,:),&
                          divup(:,:,1,1), p0, gamma1, lo, hi, ng, dx, &
                          dt_grid, cflfac)
       case (3)
          call firstdt_3d(n,uop(:,:,:,:), sop(:,:,:,:), fp(:,:,:,:),&
                          divup(:,:,:,1), p0, gamma1, lo, hi, ng, dx, &
                          dt_grid, cflfac)
       end select
       dt_hold_proc = min(dt_hold_proc,dt_grid)
    end do
    
    call parallel_reduce(dt, dt_hold_proc ,MPI_MIN)
    
  end subroutine firstdt
  
  subroutine firstdt_2d(n,u,s,force,divu,p0,gamma1,lo,hi,ng,dx,dt,cfl)

    use eos_module
    use variables, only: rho_comp, temp_comp, spec_comp
    use geometry,  only: nr
    use bl_constants_module
    
    integer, intent(in)             :: n, lo(:), hi(:), ng
    real (kind = dp_t), intent(in ) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in ) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in ) :: force(lo(1)- 1:,lo(2)- 1:,:)
    real (kind = dp_t), intent(in ) :: divu(lo(1):,lo(2):)
    real (kind = dp_t), intent(in ) :: p0(0:), gamma1(0:)
    real (kind = dp_t), intent(in ) :: dx(:)
    real (kind = dp_t), intent(out) :: dt
    real (kind = dp_t), intent(in ) :: cfl
    
    real (kind = dp_t)  :: spdx, spdy
    real (kind = dp_t)  :: pforcex, pforcey
    real (kind = dp_t)  :: ux, uy
    real (kind = dp_t)  :: eps, dt_sound, dt_divu, rho_min
    integer             :: i,j,gradp0,denom
    
    rho_min = 1.d-20
    
    eps = 1.0d-8
    
    spdx    = ZERO
    spdy    = ZERO
    pforcex = ZERO
    pforcey = ZERO
    ux      = ZERO
    uy      = ZERO

    dt = 1.d99
    
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          
          ! compute the sound speed from rho and temp
          den_eos(1)  = s(i,j,rho_comp)
          temp_eos(1) = s(i,j,temp_comp)
          xn_eos(1,:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)
          
          ! dens, temp, and xmass are inputs
          call eos(eos_input_rt, den_eos, temp_eos, &
                   npts, nspec, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, & 
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   do_diag)
          
          spdx    = max(spdx,cs_eos(1))
          spdy    = max(spdy,cs_eos(1))
          pforcex = max(pforcex,abs(force(i,j,1)))
          pforcey = max(pforcey,abs(force(i,j,2)))
          ux      = max(ux,abs(u(i,j,1)))
          uy      = max(uy,abs(u(i,j,2)))
       enddo
    enddo
    
    ux = ux / dx(1)
    uy = uy / dx(2)

    if (ux .ne. ZERO .or. uy .ne. ZERO) then
       dt = cfl / max(ux,uy)
    end if
    
    spdx = spdx / dx(1)
    spdy = spdy / dx(2)
    
    if(max(ux,uy) > max(spdx,spdy)) then
       if (parallel_IOProcessor()) then  
          print*,"Error: initial velocity greater than sound speed"
       end if
       stop
    end if

    if (spdx < eps .and. spdy < eps) then
       dt_sound = min(dx(1),dx(2))
    else
       dt_sound = cfl / max(spdx,spdy)
    end if

    if(dt_sound < dt) then
       dt = min(dt,dt_sound)
    end if
    
    if (pforcex > eps) then
       if(sqrt(2.0D0*dx(1)/pforcex) < dt) then
          dt = sqrt(2.0D0*dx(1)/pforcex)
       end if
    end if

    if (pforcey > eps) then
       if(sqrt(2.0D0*dx(2)/pforcey) < dt) then
          dt = sqrt(2.0D0*dx(2)/pforcey)
       end if
    end if
    
    ! divU constraint
    dt_divu = HUGE(dt_divu)
    
    do j = lo(2), hi(2)
       if (j .eq. 0) then
          gradp0 = (p0(j+1) - p0(j))/dx(2)
       else if (j .eq. nr(n)-1) then
          gradp0 = (p0(j) - p0(j-1))/dx(2)
       else
          gradp0 = HALF*(p0(j+1) - p0(j-1))/dx(2)
       endif
       
       do i = lo(1), hi(1)
          denom = divU(i,j) - u(i,j,2)*gradp0/(gamma1(j)*p0(j))
          if (denom > ZERO) then
             dt_divu = min(dt_divu,0.4d0*(ONE - rho_min/s(i,j,rho_comp))/denom)
          endif
       enddo
    enddo

    if(dt_divu < dt) then
       dt = dt_divu
    end if
    
  end subroutine firstdt_2d
  
  subroutine firstdt_3d(n,u,s,force,divU,p0,gamma1,lo,hi,ng,dx,dt,cfl)

    use geometry,  only: spherical, nr
    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use fill_3d_module
    use bl_constants_module

    integer, intent(in)             :: n,lo(:), hi(:), ng
    real (kind = dp_t), intent(in ) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in ) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in ) :: force(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)
    real (kind = dp_t), intent(in ) :: divU(lo(1):,lo(2):,lo(3):)  
    real (kind = dp_t), intent(in ) :: p0(0:), gamma1(0:)
    real (kind = dp_t), intent(in ) :: dx(:)
    real (kind = dp_t), intent(out) :: dt
    real (kind = dp_t), intent(in ) :: cfl
    
    real (kind = dp_t)  :: spdx, spdy, spdz
    real (kind = dp_t)  :: pforcex, pforcey, pforcez
    real (kind = dp_t)  :: ux, uy, uz
    real (kind = dp_t)  :: eps, dt_sound, dt_divu, gradp0, denom, rho_min
    integer             :: i,j,k
    
    eps = 1.0d-8
    
    rho_min = 1.d-20
    
    spdx    = ZERO
    spdy    = ZERO 
    spdz    = ZERO 
    pforcex = ZERO 
    pforcey = ZERO 
    pforcez = ZERO 
    ux      = ZERO
    uy      = ZERO
    uz      = ZERO

    dt = 1.d99
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             
             ! compute the sound speed from rho and temp
             den_eos(1) = s(i,j,k,rho_comp)
             temp_eos(1) = s(i,j,k,temp_comp)
             xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)
             
             ! dens, temp, and xmass are inputs
             call eos(eos_input_rt, den_eos, temp_eos, &
                      npts, nspec, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, & 
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      do_diag)
             
             spdx    = max(spdx,cs_eos(1))
             spdy    = max(spdy,cs_eos(1))
             spdz    = max(spdz,cs_eos(1))
             pforcex = max(pforcex,abs(force(i,j,k,1)))
             pforcey = max(pforcey,abs(force(i,j,k,2)))
             pforcez = max(pforcez,abs(force(i,j,k,3)))
             ux      = max(ux,abs(u(i,j,k,1)))
             uy      = max(uy,abs(u(i,j,k,2)))
             uz      = max(uz,abs(u(i,j,k,3)))
          enddo
       enddo
    enddo

    ux = ux / dx(1)
    uy = uy / dx(2)
    uz = uz / dx(3)

    if (ux .ne. ZERO .or. uy .ne. ZERO .or. uz .ne. ZERO) then
       dt = cfl / max(ux,uy,uz)
    end if

    spdx = spdx / dx(1)
    spdy = spdy / dx(2)
    spdz = spdz / dx(3)

    if(max(ux,uy,uz) > max(spdx,spdy,spdz)) then
       if (parallel_IOProcessor()) then
          print*,"Error: initial velocity greater than sound speed"
       end if
       stop
    end if
    
    if (spdx < eps .and. spdy < eps .and. spdz < eps) then
       dt_sound = min(dx(1),dx(2),dx(3))
    else
       dt_sound = cfl / max(spdx,spdy,spdz)
    end if

    if(dt_sound < dt) then
       dt = min(dt,dt_sound)
    end if

    if (pforcex > eps) then
       if(sqrt(2.0D0*dx(1)/pforcex) < dt) then
          dt = sqrt(2.0D0*dx(1)/pforcex)
       end if
    end if

    if (pforcey > eps) then
       if(sqrt(2.0D0*dx(2)/pforcey) < dt) then
          dt = sqrt(2.0D0*dx(2)/pforcey)
       end if
    end if

    if (pforcez > eps) then
       if(sqrt(2.0D0*dx(3)/pforcez) < dt) then
          dt = sqrt(2.0D0*dx(3)/pforcez)
       end if
    end if
    
    ! divU constraint
    dt_divu = HUGE(dt_divu)
    
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
             denom = divU(i,j,k) - u(i,j,k,3)*gradp0/(gamma1(k)*p0(k))
             if (denom > ZERO) then
                dt_divu = min(dt_divu,0.4d0*(ONE - rho_min/s(i,j,k,rho_comp))/denom)
             endif
          enddo
       enddo
    enddo

    if(dt_divu < dt) then
       dt = dt_divu
    end if
    
  end subroutine firstdt_3d
  
end module firstdt_module
