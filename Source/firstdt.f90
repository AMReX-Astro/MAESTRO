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

   subroutine firstdt (istep, u, s, force, p0, t0, dx, dt)

      integer        , intent(in   ) :: istep
      type(multifab) , intent(inout) :: u,s,force
      real(kind=dp_t), intent(in   ) :: p0(:), t0(:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      real(kind=dp_t), intent(  out) :: dt

      real(kind=dp_t), pointer:: uop(:,:,:,:)
      real(kind=dp_t), pointer:: sop(:,:,:,:)
      real(kind=dp_t), pointer:: fp(:,:,:,:)
      integer :: lo(u%dim),hi(u%dim),ng,dm
      real(kind=dp_t) :: dt_hold
      integer         :: i

      ng = u%ng
      dm = u%dim

      dt_hold  = 1.d20

      do i = 1, u%nboxes
         if ( multifab_remote(u, i) ) cycle
         uop => dataptr(u, i)
         sop => dataptr(s, i)
          fp => dataptr(force, i)
         lo =  lwb(get_box(u, i))
         hi =  upb(get_box(u, i))
         select case (dm)
            case (2)
              call firstdt_2d(uop(:,:,1,:), sop(:,:,1,:), fp(:,:,1,:),&
                              p0, t0, lo, hi, ng, dx, dt)
            case (3)
              call firstdt_3d(uop(:,:,:,:), sop(:,:,:,:), fp(:,:,:,:),&
                              p0, t0, lo, hi, ng, dx, dt)
         end select
         dt_hold = min(dt_hold,dt)
      end do

      dt = dt_hold

      print *,'Computing dt at istep ',istep,' to be ',dt

    end subroutine firstdt

   subroutine firstdt_2d (u, s, force, p0, t0, lo, hi, ng, dx, dt)

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in ) ::     u(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in ) ::     s(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in ) :: force(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent(in ) :: p0(0:), t0(0:)
      real (kind = dp_t), intent(in ) :: dx(:)
      real (kind = dp_t), intent(out) :: dt

!     Local variables
      real (kind = dp_t)  :: spdx, spdy
      real (kind = dp_t)  :: pforcex, pforcey
      real (kind = dp_t)  :: eps
      integer             :: i,j

      eps = 1.0e-8

      spdx  = 0.0D0 
      spdy  = 0.0D0 
      pforcex = 0.0D0 
      pforcey = 0.0D0 

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

            spdx    = max(spdx ,max(abs(u(i,j,1)), cs_row(1)))
            spdy    = max(spdy ,max(abs(u(i,j,2)), cs_row(1)))
            pforcex = max(pforcex,abs(force(i,j,1)))
            pforcey = max(pforcey,abs(force(i,j,2)))
         enddo
      enddo

      spdx = spdx / dx(1)
      spdy = spdy / dx(2)

      if (spdx < eps .and. spdy < eps) then

        dt = min(dx(1),dx(2))

      else

        dt = 1.0D0  / max(spdx,spdy)

      endif

      if (pforcex > eps) then
        dt = min(dt,sqrt(2.0D0 *dx(1)/pforcex))
      endif

      if (pforcey > eps) then
        dt = min(dt,sqrt(2.0D0 *dx(2)/pforcey))
      endif

    end subroutine firstdt_2d

    subroutine firstdt_3d (u, s, force, p0, t0, lo, hi, ng, dx, dt)

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in ) ::     u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(in ) ::     s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
      real (kind = dp_t), intent(in ) :: force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
      real (kind = dp_t), intent(in ) :: p0(0:), t0(0:)
      real (kind = dp_t), intent(in ) :: dx(:)
      real (kind = dp_t), intent(out) :: dt

!     Local variables
      real (kind = dp_t)  :: spdx, spdy, spdz
      real (kind = dp_t)  :: pforcex, pforcey, pforcez
      real (kind = dp_t)  :: eps
      integer             :: i,j,k

      real (kind=dp_t), allocatable :: t0_cart(:,:,:)
      real (kind=dp_t), allocatable :: p0_cart(:,:,:)

      eps = 1.0e-8

      spdx    = ZERO
      spdy    = ZERO 
      spdz    = ZERO 
      pforcex = ZERO 
      pforcey = ZERO 
      pforcez = ZERO 

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
               
               spdx    = max(spdx ,max(abs(u(i,j,k,1)),cs_row(1)))
               spdy    = max(spdy ,max(abs(u(i,j,k,2)),cs_row(1)))
               spdz    = max(spdz ,max(abs(u(i,j,k,3)),cs_row(1)))
               pforcex = max(pforcex,abs(force(i,j,k,1)))
               pforcey = max(pforcey,abs(force(i,j,k,2)))
               pforcez = max(pforcez,abs(force(i,j,k,3)))
            enddo
         enddo
      enddo
      
      spdx = spdx / dx(1)
      spdy = spdy / dx(2)
      spdz = spdz / dx(3)

      if (spdx < eps .and. spdy < eps .and. spdz < eps) then

        dt = min(dx(1),dx(2),dx(3))

      else

        dt = 1.0D0  / max(spdx,spdy,spdz)

      endif

      if (pforcex > eps) then
        dt = min(dt,sqrt(2.0D0*dx(1)/pforcex))
      endif

      if (pforcey > eps) then
        dt = min(dt,sqrt(2.0D0*dx(2)/pforcey))
      endif

      if (pforcez > eps) then
        dt = min(dt,sqrt(2.0D0*dx(3)/pforcez))
      endif

      if (spherical == 1) &
        deallocate(t0_cart,p0_cart)

    end subroutine firstdt_3d

  end module firstdt_module
