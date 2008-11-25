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

  subroutine firstdt(n,u,s,force,divU,normal,p0,gamma1bar,dx,cflfac,dt)

    use geometry, only: dm
    use variables, only: rel_eps

    integer        , intent(in   ) :: n
    type(multifab) , intent(in   ) :: u,s,force,divU,normal
    real(kind=dp_t), intent(in   ) :: p0(0:), cflfac, gamma1bar(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(  out) :: dt
    
    real(kind=dp_t), pointer:: uop(:,:,:,:)
    real(kind=dp_t), pointer:: sop(:,:,:,:)
    real(kind=dp_t), pointer:: fp(:,:,:,:)
    real(kind=dp_t), pointer:: divup(:,:,:,:)
    real(kind=dp_t), pointer:: np(:,:,:,:)

    integer :: lo(dm),hi(dm),ng_u,ng_s,ng_f,ng_dU,ng_n,i
    real(kind=dp_t) :: dt_proc,dt_grid
    real(kind=dp_t) :: umax,umax_proc,umax_grid
    
    ng_u = u%ng
    ng_s = s%ng
    ng_f = force%ng
    ng_dU = divU%ng
    
    dt_proc   = 1.d99
    dt_grid   = 1.d99
    umax_proc = 0.d0
    umax_grid = 0.d0
    
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
          call firstdt_2d(n, uop(:,:,1,:), ng_u, sop(:,:,1,:), ng_s, &
                          fp(:,:,1,:), ng_f, divup(:,:,1,1), ng_dU, &
                          p0, gamma1bar, lo, hi, dx, dt_grid, umax_grid, cflfac)
       case (3)
          np => dataptr(normal, i)
          ng_n = normal%ng
          call firstdt_3d(n, uop(:,:,:,:), ng_u, sop(:,:,:,:), ng_s, &
                          fp(:,:,:,:), ng_f, divup(:,:,:,1), ng_dU, np(:,:,:,:), ng_n, &
                          p0, gamma1bar, lo, hi, dx, dt_grid, umax_grid, cflfac)
       end select

       dt_proc   = min(  dt_proc,   dt_grid)
       umax_proc = max(umax_proc, umax_grid)

    end do
    
    call parallel_reduce(  dt,   dt_proc, MPI_MIN)
    call parallel_reduce(umax, umax_proc, MPI_MAX)

    rel_eps = 1.d-8*umax
    
  end subroutine firstdt
  
  subroutine firstdt_2d(n,u,ng_u,s,ng_s,force,ng_f,divu,ng_dU, &
                        p0,gamma1bar,lo,hi,dx,dt,umax,cfl)

    use eos_module
    use variables, only: rho_comp, temp_comp, spec_comp
    use geometry,  only: nr
    use bl_constants_module
    use probin_module, only: use_soundspeed_firstdt, use_divu_firstdt
    
    integer, intent(in)             :: n, lo(:), hi(:), ng_u, ng_s, ng_f, ng_dU
    real (kind = dp_t), intent(in ) ::     u(lo(1)-ng_u :,lo(2)-ng_u :,:)  
    real (kind = dp_t), intent(in ) ::     s(lo(1)-ng_s :,lo(2)-ng_s :,:)  
    real (kind = dp_t), intent(in ) :: force(lo(1)-ng_f :,lo(2)-ng_f :,:)
    real (kind = dp_t), intent(in ) ::  divu(lo(1)-ng_dU:,lo(2)-ng_dU:)
    real (kind = dp_t), intent(in ) :: p0(0:), gamma1bar(0:)
    real (kind = dp_t), intent(in ) :: dx(:)
    real (kind = dp_t), intent(out) :: dt, umax
    real (kind = dp_t), intent(in ) :: cfl
    
    real (kind = dp_t)  :: spdx, spdy
    real (kind = dp_t)  :: pforcex, pforcey
    real (kind = dp_t)  :: ux, uy
    real (kind = dp_t)  :: eps, dt_divu, dt_sound, rho_min
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
    umax = ZERO
    
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
    
    umax = max(umax,ux)
    umax = max(umax,uy)

    ux = ux / dx(1)
    uy = uy / dx(2)
    
    spdx = spdx / dx(1)
    spdy = spdy / dx(2)

    ! advective constraint
    if (ux .ne. ZERO .or. uy .ne. ZERO) then
       dt = cfl / max(ux,uy)
    else if (spdx .ne. ZERO .and. spdy .ne. ZERO) then
       dt = cfl / max(spdx,spdy)
    end if

    ! sound speed constraint
    if (use_soundspeed_firstdt) then
       if (spdx .eq. ZERO .and. spdy .eq. ZERO) then
          dt_sound = 1.d99
       else
          dt_sound = cfl / max(spdx,spdy)
       end if
       dt = min(dt,dt_sound)
    end if
    
    ! force constraints
    if (pforcex > eps) dt = min(dt,sqrt(2.0D0*dx(1)/pforcex))
    if (pforcey > eps) dt = min(dt,sqrt(2.0D0*dx(2)/pforcey))
    
    ! divU constraint
    if (use_divu_firstdt) then
       
       dt_divu = 1.d99
       
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
                dt_divu = min(dt_divu,0.4d0*(ONE - rho_min/s(i,j,rho_comp))/denom)
             endif
          enddo
       enddo
    
       dt = min(dt,dt_divu)

    end if

  end subroutine firstdt_2d
  
  subroutine firstdt_3d(n,u,ng_u,s,ng_s,force,ng_f,divU,ng_dU,normal,ng_n, &
                        p0,gamma1bar,lo,hi,dx,dt,umax,cfl)

    use geometry,  only: spherical, nr, dr, nr_fine
    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use bl_constants_module
    use probin_module, only: use_soundspeed_firstdt, use_divu_firstdt
    use fill_3d_module

    integer, intent(in)             :: n,lo(:), hi(:), ng_u, ng_s, ng_f, ng_dU, ng_n
    real (kind = dp_t), intent(in ) ::      u(lo(1)-ng_u :,lo(2)-ng_u :,lo(3)-ng_u :,:)
    real (kind = dp_t), intent(in ) ::      s(lo(1)-ng_s :,lo(2)-ng_s :,lo(3)-ng_s :,:)
    real (kind = dp_t), intent(in ) ::  force(lo(1)-ng_f :,lo(2)-ng_f :,lo(3)-ng_f :,:)
    real (kind = dp_t), intent(in ) ::   divU(lo(1)-ng_dU:,lo(2)-ng_dU:,lo(3)-ng_dU:) 
    real (kind = dp_t), intent(in ) :: normal(lo(1)-ng_n :,lo(2)-ng_n :,lo(3)-ng_n :,:) 
    real (kind = dp_t), intent(in ) :: p0(0:), gamma1bar(0:)
    real (kind = dp_t), intent(in ) :: dx(:)
    real (kind = dp_t), intent(out) :: dt, umax
    real (kind = dp_t), intent(in ) :: cfl
    
    real (kind = dp_t)  :: spdx, spdy, spdz
    real (kind = dp_t)  :: pforcex, pforcey, pforcez
    real (kind = dp_t)  :: ux, uy, uz, gp_dot_u, gamma1bar_p_avg
    real (kind = dp_t)  :: eps, dt_divu, dt_sound, gradp0, denom, rho_min
    integer             :: i,j,k,r

    real (kind = dp_t) :: gp0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3)
    real (kind = dp_t) ::      gp0(0:nr_fine)

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
    umax = ZERO

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

    umax = max(umax,ux)
    umax = max(umax,uy)
    umax = max(umax,uz)

    ux = ux / dx(1)
    uy = uy / dx(2)
    uz = uz / dx(3)
    
    spdx = spdx / dx(1)
    spdy = spdy / dx(2)
    spdz = spdz / dx(3)
    
    ! advective constraint
    if (ux .ne. ZERO .or. uy .ne. ZERO .or. uz .ne. ZERO) then
       dt = cfl / max(ux,uy,uz)
    else if (spdx .ne. ZERO .and. spdy .ne. ZERO .and. spdz .ne. ZERO) then
       dt = cfl / max(spdx,spdy,spdz)
    end if

    ! sound speed constraint
    if (use_soundspeed_firstdt) then
       if (spdx .eq. ZERO .and. spdy .eq. ZERO .and. spdz .eq. ZERO) then
          dt_sound = 1.d99
       else
          dt_sound = cfl / max(spdx,spdy,spdz)
       end if
       dt = min(dt,dt_sound)
    end if

    ! force constraints
    if (pforcex > eps) dt = min(dt,sqrt(2.0D0*dx(1)/pforcex))
    if (pforcey > eps) dt = min(dt,sqrt(2.0D0*dx(2)/pforcey))
    if (pforcez > eps) dt = min(dt,sqrt(2.0D0*dx(3)/pforcez))
    
    ! divU constraint
    if (use_divu_firstdt) then

       dt_divu = 1.d99
       
       if (spherical .eq. 0) then

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
                      dt_divu = min(dt_divu,0.4d0*(ONE - rho_min/s(i,j,k,rho_comp))/denom)
                   endif
                enddo
             enddo
          enddo

       else

          ! spherical divU constraint
          do r=1,nr_fine-1
             gamma1bar_p_avg = HALF * (gamma1bar(r)*p0(r) + gamma1bar(r-1)*p0(r-1))
             gp0(r) = ( (p0(r) - p0(r-1))/dr(1) ) / gamma1bar_p_avg
          end do
          gp0(nr_fine) = gp0(nr_fine-1)
          gp0(      0) = gp0(        1)

          call put_1d_array_on_cart_3d_sphr(.true.,.true.,gp0,gp0_cart,lo,hi,dx,0, &
                                            ng_n,normal)
          
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
                   
                enddo
             enddo
          enddo
          
       end if

       dt = min(dt,dt_divu)

    end if

  end subroutine firstdt_3d
  
end module firstdt_module
