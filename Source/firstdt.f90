! in the initial projection, if the velocity is zero everywhere and there
! is no buoyancy terms, then estdt will give 0 as the initial timestep.
! Here we also consider dx/c, where c is the sound speed.  After the
! first iteration, we should have a source term tat gave rise to a 
! velocity field, and we can use the normal estdt.

module firstdt_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use eos_module
  use variables
  use geometry
  use fill_3d_module

  implicit none

contains

   subroutine firstdt (istep, u, s, force, divU, p0, gam1, t0, dx, cflfac, dt, verbose)

      integer        , intent(in   ) :: istep
      type(multifab) , intent(in   ) :: u,s,force,divU
      real(kind=dp_t), intent(in   ) :: p0(0:), cflfac, t0(0:), gam1(0:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      real(kind=dp_t), intent(  out) :: dt
      integer        , intent(in   ) :: verbose

      real(kind=dp_t), pointer:: uop(:,:,:,:)
      real(kind=dp_t), pointer:: sop(:,:,:,:)
      real(kind=dp_t), pointer:: fp(:,:,:,:)
      real(kind=dp_t), pointer:: divup(:,:,:,:)
      integer :: lo(u%dim),hi(u%dim),ng,dm
      real(kind=dp_t) :: dt_hold_proc,dt_grid
      integer         :: i

      ng = u%ng
      dm = u%dim

      dt_hold_proc  = 1.d20
      dt_grid       = 1.d20

      do i = 1, u%nboxes
         if ( multifab_remote(u, i) ) cycle
         uop => dataptr(u, i)
         sop => dataptr(s, i)
         fp => dataptr(force, i)
         divup => dataptr(divU,i)
         lo =  lwb(get_box(u, i))
         hi =  upb(get_box(u, i))
         select case (dm)
            case (2)
              call firstdt_2d(uop(:,:,1,:), sop(:,:,1,:), fp(:,:,1,:),&
                              divup(:,:,1,1), p0, gam1, t0, lo, hi, ng, dx, &
                              dt_grid, cflfac, verbose)
            case (3)
              call firstdt_3d(uop(:,:,:,:), sop(:,:,:,:), fp(:,:,:,:),&
                              divup(:,:,:,1), p0, gam1, t0, lo, hi, ng, dx, &
                              dt_grid, cflfac, verbose)
         end select
         dt_hold_proc = min(dt_hold_proc,dt_grid)
      end do

      call parallel_reduce(dt, dt_hold_proc ,MPI_MIN)

      if (parallel_IOProcessor() .and. verbose .ge. 1) &
        print *,'Using firstdt, at istep',istep,', dt =',dt

    end subroutine firstdt

   subroutine firstdt_2d (u, s, force, divu, p0, gam1, t0, lo, hi, ng, dx, &
                          dt, cfl, verbose)

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in ) :: u(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in ) :: s(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in ) :: force(lo(1)- 1:,lo(2)- 1:,:)
      real (kind = dp_t), intent(in ) :: divu(lo(1):,lo(2):)
      real (kind = dp_t), intent(in ) :: p0(0:), gam1(0:), t0(0:)
      real (kind = dp_t), intent(in ) :: dx(:)
      real (kind = dp_t), intent(out) :: dt
      real (kind = dp_t), intent(in ) :: cfl
      integer           , intent(in ) :: verbose

!     Local variables
      real (kind = dp_t)  :: spdx, spdy
      real (kind = dp_t)  :: pforcex, pforcey
      real (kind = dp_t)  :: ux, uy
      real (kind = dp_t)  :: eps, dt_divu, rho_min
      integer             :: i,j,nr,gradp0,denom

      nr = size(p0,dim=1)

      rho_min = 1.d-20

      eps = 1.0e-8

      spdx    = 0.0D0 
      spdy    = 0.0D0
      pforcex = 0.0D0 
      pforcey = 0.0D0 
      ux      = 0.0D0
      uy      = 0.0D0

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)

            ! compute the sound speed from rho and p0
            den_row(1) = s(i,j,rho_comp)
            temp_row(1) = t0(j)
            p_row(1) = p0(j)
            xn_zone(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

            ! (rho,P) --> T,h
            input_flag = 4

            call eos(input_flag, den_row, temp_row, &
                     npts, nspec, &
                     xn_zone, aion, zion, &
                     p_row, h_row, e_row, & 
                     cv_row, cp_row, xne_row, eta_row, pele_row, &
                     dpdt_row, dpdr_row, dedt_row, dedr_row, &
                     dpdX_row, dhdX_row, &
                     gam1_row, cs_row, s_row, &
                     dsdt_row, dsdr_row, &
                     do_diag)

            spdx = max(spdx,cs_row(1))
            spdy = max(spdy,cs_row(1))
            pforcex = max(pforcex,abs(force(i,j,1)))
            pforcey = max(pforcey,abs(force(i,j,2)))
            ux = max(ux,abs(u(i,j,1)))
            uy = max(uy,abs(u(i,j,2)))
         enddo
      enddo

      ux = ux / dx(1)
      uy = uy / dx(2)

      spdx = spdx / dx(1)
      spdy = spdy / dx(2)

      ! if ux or uy is non-zero use it
      ! otherwise, use soundspeed
      if(ux .ne. ZERO .or. uy .ne. ZERO) then
         dt = cfl / max(ux,uy)
         if ( parallel_IOProcessor() .and. verbose .ge. 1) then
            print*, ''
            print*, 'advective dt =',dt
         endif         
      else if (spdx < eps .and. spdy < eps) then
         dt = min(dx(1),dx(2))
         if ( parallel_IOProcessor() .and. verbose .ge. 1) then
            print*, 'sound speed < eps; dt =',dt
         endif      
      else
         dt = cfl / max(spdx,spdy)
         if ( parallel_IOProcessor() .and. verbose .ge. 1) then
            print*, 'sound speed dt =',dt
         endif
      endif

      if (pforcex > eps) then
         dt = min(dt,sqrt(2.0D0*dx(1)/pforcex))
         if ( parallel_IOProcessor() .and. verbose .ge. 1) then
            print*, 'pforcex dt =',sqrt(2.0D0*dx(1)/pforcex)
         endif      
      endif

      if (pforcey > eps) then
         dt = min(dt,sqrt(2.0D0*dx(2)/pforcey))
         if ( parallel_IOProcessor() .and. verbose .ge. 1) then
            print*, 'pforcey dt =',sqrt(2.0D0*dx(2)/pforcey)
         endif
      endif

     ! divU constraint
     dt_divu = 1.d30
     
     do j = lo(2), hi(2)
        if (j .eq. 0) then
           gradp0 = (p0(j+1) - p0(j))/dx(2)
        else if (j .eq. nr-1) then
           gradp0 = (p0(j) - p0(j-1))/dx(2)
        else
           gradp0 = HALF*(p0(j+1) - p0(j-1))/dx(2)
        endif

        do i = lo(1), hi(1)
           denom = divU(i,j) - u(i,j,2)*gradp0/(gam1(j)*p0(j))
           if (denom > ZERO) then
              dt_divu = min(dt_divu, &
                            0.4d0*(ONE - rho_min/s(i,j,rho_comp))/denom)
           endif
        enddo
     enddo

     if ( parallel_IOProcessor() .and. verbose .ge. 1) then
        print*, 'divu_dt =',dt_divu
     endif     

    end subroutine firstdt_2d

    subroutine firstdt_3d (u, s, force, divU, p0, gam1, t0, lo, hi, ng, dx, dt, cfl, &
                           verbose)

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in ) ::     u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(in ) ::     s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(in ) :: force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent(in ) :: divU(lo(1):,lo(2):,lo(3):)  
      real (kind = dp_t), intent(in ) :: p0(0:), t0(0:), gam1(0:)
      real (kind = dp_t), intent(in ) :: dx(:)
      real (kind = dp_t), intent(out) :: dt
      real (kind = dp_t), intent(in ) :: cfl
      integer           , intent(in ) :: verbose

!     Local variables
      real (kind = dp_t)  :: spdx, spdy, spdz
      real (kind = dp_t)  :: pforcex, pforcey, pforcez
      real (kind = dp_t)  :: ux, uy, uz
      real (kind = dp_t)  :: eps, dt_divu, gradp0, denom, rho_min
      integer             :: i,j,k,nr

      real (kind=dp_t), allocatable :: t0_cart(:,:,:)
      real (kind=dp_t), allocatable :: p0_cart(:,:,:)

      eps = 1.0e-8

      rho_min = 1.d-20

      nr = size(p0,dim=1)

      spdx    = ZERO
      spdy    = ZERO 
      spdz    = ZERO 
      pforcex = ZERO 
      pforcey = ZERO 
      pforcez = ZERO 
      ux      = ZERO
      uy      = ZERO
      uz      = ZERO

      if (spherical == 1) then
         allocate(t0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
         call fill_3d_data(t0_cart,t0,lo,hi,dx,0)
         
         allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
         call fill_3d_data(p0_cart,p0,lo,hi,dx,0)
      endif

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               den_row(1) = s(i,j,k,rho_comp)
               xn_zone(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

               if (spherical == 1) then
                  p_row(1) = p0_cart(i,j,k)
                  temp_row(1) = t0_cart(i,j,k)
               else
                  p_row(1) = p0(k)
                  temp_row(1) = t0(k)
               endif

               ! (rho,P) --> T,h
               input_flag = 4

               call eos(input_flag, den_row, temp_row, &
                        npts, nspec, &
                        xn_zone, aion, zion, &
                        p_row, h_row, e_row, & 
                        cv_row, cp_row, xne_row, eta_row, pele_row, &
                        dpdt_row, dpdr_row, dedt_row, dedr_row, &
                        dpdX_row, dhdX_row, &
                        gam1_row, cs_row, s_row, &
                        dsdt_row, dsdr_row, &
                        do_diag)
               
               spdx    = max(spdx,cs_row(1))
               spdy    = max(spdy,cs_row(1))
               spdz    = max(spdz,cs_row(1))
               pforcex = max(pforcex,abs(force(i,j,k,1)))
               pforcey = max(pforcey,abs(force(i,j,k,2)))
               pforcez = max(pforcez,abs(force(i,j,k,3)))
               ux      = max(ux,abs(u(i,j,k,1)))
               uy      = max(uy,abs(u(i,j,k,2)))
               uz      = max(uz,abs(u(i,j,k,3)))
            enddo
         enddo
      enddo
      
      spdx = spdx / dx(1)
      spdy = spdy / dx(2)
      spdz = spdz / dx(3)

      ux = ux / dx(1)
      uy = uy / dx(2)
      uz = uz / dx(3)

      ! if ux, uy, or uz is non-zero use it
      ! otherwise, use soundspeed
      if(ux .ne. ZERO .or. uy .ne. ZERO .or. uz .ne. ZERO) then
         dt = cfl / max(ux,uy,uz)
         if ( parallel_IOProcessor() .and. verbose .ge. 1) then
            print*, 'advective dt =',dt
         endif        
      else if (spdx < eps .and. spdy < eps .and. spdz < eps) then
         dt = min(dx(1),dx(2),dx(3))
         if ( parallel_IOProcessor() .and. verbose .ge. 1) then
            print*, 'sound speed < eps; dt =',dt
         endif
      else
         dt = cfl / max(spdx,spdy,spdz)
         if ( parallel_IOProcessor() .and. verbose .ge. 1) then
            print*, 'sound speed dt =',dt
         endif
      endif

      if (pforcex > eps) then
         dt = min(dt,sqrt(2.0D0*dx(1)/pforcex))
         if ( parallel_IOProcessor() .and. verbose .ge. 1) then
            print*, 'pforcex dt =',sqrt(2.0D0*dx(1)/pforcex)
         endif     
      endif

      if (pforcey > eps) then
         dt = min(dt,sqrt(2.0D0*dx(2)/pforcey))
         if ( parallel_IOProcessor() .and. verbose .ge. 1) then
            print*, 'pforcey dt =',sqrt(2.0D0*dx(2)/pforcey)
         endif     
      endif

      if (pforcez > eps) then
         dt = min(dt,sqrt(2.0D0*dx(3)/pforcez))
         if ( parallel_IOProcessor() .and. verbose .ge. 1) then
            print*, 'pforcez dt =',sqrt(2.0D0*dx(3)/pforcez)
         endif     
      endif

     ! divU constraint
     dt_divu = 1.d30

     do k = lo(3), hi(3)
        if (k .eq. 0) then
           gradp0 = (p0(k+1) - p0(k))/dx(3)
        else if (k .eq. nr-1) then
           gradp0 = (p0(k) - p0(k-1))/dx(3)
        else
           gradp0 = HALF*(p0(k+1) - p0(k-1))/dx(3)
        endif
        
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              denom = divU(i,j,k) - u(i,j,k,3)*gradp0/(gam1(k)*p0(k))
              if (denom > ZERO) then
                dt_divu = min(dt_divu, &
                              0.4d0*(ONE - rho_min/s(i,j,k,rho_comp))/denom)
              endif
           enddo
        enddo
     enddo

     if ( parallel_IOProcessor() .and. verbose .ge. 1) then
        print*, 'divu_dt =',dt_divu
     endif    

     if (spherical == 1) &
          deallocate(t0_cart,p0_cart)

    end subroutine firstdt_3d

  end module firstdt_module
