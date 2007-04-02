module estdt_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use fill_3d_module
  use geometry
  use variables

  implicit none

contains

   subroutine estdt (istep, u, s, force, divU, normal, w0, p0, gam1, dx, cflfac, dtold, dt, verbose)

      integer        , intent(in ) :: istep
      type(multifab) , intent(in ) :: u
      type(multifab) , intent(in ) :: s
      type(multifab) , intent(in ) :: force
      type(multifab) , intent(in ) :: divU
      type(multifab) , intent(in ) :: normal
      real(kind=dp_t), intent(in ) :: w0(:), p0(:), gam1(:)
      real(kind=dp_t), intent(in ) :: dx(:)
      real(kind=dp_t), intent(in ) :: cflfac, dtold
      real(kind=dp_t), intent(out) :: dt
      integer        , intent(in ) :: verbose

      real(kind=dp_t), pointer:: uop(:,:,:,:)
      real(kind=dp_t), pointer:: sop(:,:,:,:)
      real(kind=dp_t), pointer:: fp(:,:,:,:)
      real(kind=dp_t), pointer:: np(:,:,:,:)
      real(kind=dp_t), pointer:: dUp(:,:,:,:)     
      integer :: lo(u%dim),hi(u%dim),ng,dm
      real(kind=dp_t) :: dt_adv,  dt_adv_grid,dt_adv_proc
      real(kind=dp_t) :: dt_divu,dt_divu_grid,dt_divu_proc
      real(kind=dp_t) :: dtchange
      integer         :: i

      real(kind=dp_t), parameter :: rho_min = 1.d-20

      ng = u%ng
      dm = u%dim

      dtchange = 1.1d0

      dt_adv_grid   = 1.d20
      dt_adv_proc   = 1.d20
      dt_divu_grid  = 1.d20
      dt_divu_proc  = 1.d20

      do i = 1, u%nboxes
         if ( multifab_remote(u, i) ) cycle
         uop => dataptr(u, i)
         sop => dataptr(s, i)
          fp => dataptr(force, i)
         dUp => dataptr(divU, i)

         lo =  lwb(get_box(u, i))
         hi =  upb(get_box(u, i))

         select case (dm)
            case (2)
              call estdt_2d(uop(:,:,1,:), sop(:,:,1,:), fp(:,:,1,:), dUp(:,:,1,1), &
                            w0, p0, gam1, lo, hi, ng, dx, rho_min, dt_adv_grid, dt_divu_grid)
            case (3)
              if (spherical .eq. 1) then
                np => dataptr(normal, i)
                call estdt_3d_sphr(uop(:,:,:,:), sop(:,:,:,:), fp(:,:,:,:), dUp(:,:,:,1), np(:,:,:,:), &
                                   w0, p0, gam1, lo, hi, ng, dx, rho_min, dt_adv_grid, dt_divu_grid)
              else
                call estdt_3d_cart(uop(:,:,:,:), sop(:,:,:,:), fp(:,:,:,:), dUp(:,:,:,1), &
                                   w0, p0, gam1, lo, hi, ng, dx, rho_min, dt_adv_grid, dt_divu_grid)
              end if
         end select

         dt_adv_proc  = min(dt_adv_proc ,dt_adv_grid)
         dt_divu_proc = min(dt_divu_proc,dt_divu_grid)
      end do

      ! This sets dt to be the min of dt_proc over all processors.
      call parallel_reduce(dt_adv ,dt_adv_proc ,MPI_MIN)
      call parallel_reduce(dt_divu,dt_divu_proc,MPI_MIN)

      if (parallel_IOProcessor() .and. verbose .ge. 1) &
        write(6,*) '%%% timesteps (source, advective): ', dt_divu, dt_adv

      dt = min(dt_adv,dt_divu)

      dt = dt * cflfac

      if (dtold .gt. 0.0D0 ) dt = min(dt,dtchange*dtold)

      if (parallel_IOProcessor() .and. verbose .ge. 1) &
        write(6,*) 'Computing dt at istep ',istep,' to be ',dt

   end subroutine estdt


   subroutine estdt_2d (u, s, force, divU, w0, p0, gam1, lo, hi, ng, dx, rho_min, dt_adv, dt_divu)

     integer, intent(in) :: lo(:), hi(:), ng
     real (kind = dp_t), intent(in ) ::     u(lo(1)-ng:,lo(2)-ng:,:)  
     real (kind = dp_t), intent(in ) ::     s(lo(1)-ng:,lo(2)-ng:,:)  
     real (kind = dp_t), intent(in ) :: force(lo(1)- 1:,lo(2)- 1:,:)  
     real (kind = dp_t), intent(in ) ::  divU(lo(1):,lo(2):)
     real (kind = dp_t), intent( in) ::   w0(0:), p0(0:), gam1(0:)
     real (kind = dp_t), intent(in ) :: dx(:)
     real (kind = dp_t), intent(in ) :: rho_min
     real (kind = dp_t), intent(out) :: dt_adv,dt_divu

!    Local variables
     real (kind = dp_t)  :: spdx, spdy, spdr
     real (kind = dp_t)  :: pforcex, pforcey
     real (kind = dp_t)  :: eps
     real (kind = dp_t)  :: denom, gradp0
     integer             :: i,j,nr

     nr = size(p0,dim=1)

     eps = 1.0e-8

     ! advective constraints
     spdx  = 0.0D0 
     spdy  = 0.0D0 

     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           spdx = max(spdx ,abs(u(i,j,1)))
           spdy = max(spdy ,abs(u(i,j,2)+w0(j)))
        enddo
     enddo

     spdx = spdx / dx(1)
     spdy = spdy / dx(2)

     spdr = ZERO 
     do j = lo(2),hi(2)
        spdr = max(spdr ,abs(w0(j)))
     enddo
     spdr = spdr / dx(2)

     if (spdx < eps .and. spdy < eps .and. spdr < eps) then
        dt_adv = min(dx(1),dx(2))
     else
        dt_adv = 1.0D0  / max(spdx,spdy,spdr)
     endif


     ! force constraints
     pforcex = 0.0D0 
     pforcey = 0.0D0 
     
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           pforcex = max(pforcex,abs(force(i,j,1)))
           pforcey = max(pforcey,abs(force(i,j,2)))
        enddo
     enddo

     if (pforcex > eps) then
        dt_adv = min(dt_adv,sqrt(2.0D0 *dx(1)/pforcex))
     endif

     if (pforcey > eps) then
        dt_adv = min(dt_adv,sqrt(2.0D0 *dx(2)/pforcey))
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
                            HALF*(ONE - rho_min/s(i,j,rho_comp))/denom)
           endif

        enddo
     enddo

   end subroutine estdt_2d

   subroutine estdt_3d_cart (u, s, force, divU, w0, p0, gam1, lo, hi, ng, dx, rho_min, dt_adv, dt_divu)

     integer, intent(in) :: lo(:), hi(:), ng
     real (kind = dp_t), intent(in ) ::      u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
     real (kind = dp_t), intent(in ) ::      s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
     real (kind = dp_t), intent(in ) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
     real (kind = dp_t), intent(in ) ::   divU(lo(1):,lo(2):,lo(3):)
     real (kind = dp_t), intent( in) ::   w0(0:), p0(0:), gam1(0:)
     real (kind = dp_t), intent(in ) :: dx(:)
     real (kind = dp_t), intent(in ) :: rho_min
     real (kind = dp_t), intent(out) :: dt_adv, dt_divu

!    Local variables
     real (kind = dp_t), allocatable :: w0_cart(:,:,:,:)
     real (kind = dp_t)  :: spdx, spdy, spdz, spdr
     real (kind = dp_t)  :: pforcex, pforcey, pforcez
     real (kind = dp_t)  :: eps,denom,gradp0
     integer             :: i,j,k,nr

     eps = 1.0e-8

     nr = size(p0,dim=1)

     spdx    = ZERO
     spdy    = ZERO 
     spdz    = ZERO 

     ! Limit dt based on velocity terms
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              spdx = max(spdx ,abs(u(i,j,k,1)))
              spdy = max(spdy ,abs(u(i,j,k,2)))
              spdz = max(spdz ,abs(u(i,j,k,3)+w0(k)))
           enddo
        enddo
     enddo
     
     spdr = ZERO 
     do k = lo(3),hi(3)
        spdr = max(spdr ,abs(w0(k)))
     enddo
     spdr = spdr / dx(3)

     spdx = spdx / dx(1)
     spdy = spdy / dx(2)
     spdz = spdz / dx(3)
      
     if (spdx < eps .and. spdy < eps .and. spdz < eps .and. spdr < eps) then
        dt_adv = min(dx(1),dx(2),dx(3))
     else
        dt_adv = 1.0D0  / max(spdx,spdy,spdz,spdr)
     endif

     ! Limit dt based on forcing terms
     pforcex = ZERO 
     pforcey = ZERO 
     pforcez = ZERO 

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              pforcex = max(pforcex,abs(force(i,j,k,1)))
              pforcey = max(pforcey,abs(force(i,j,k,2)))
              pforcez = max(pforcez,abs(force(i,j,k,3)))
           enddo
        enddo
     enddo
     
     if (pforcex > eps) then
        dt_adv = min(dt_adv,sqrt(2.0D0*dx(1)/pforcex))
     endif
     
     if (pforcey > eps) then
        dt_adv = min(dt_adv,sqrt(2.0D0*dx(2)/pforcey))
     endif
     
     if (pforcez > eps) then
        dt_adv = min(dt_adv,sqrt(2.0D0*dx(3)/pforcez))
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
                              HALF*(ONE - rho_min/s(i,j,k,rho_comp))/denom)
              endif

           enddo
        enddo
     enddo

   end subroutine estdt_3d_cart

   subroutine estdt_3d_sphr (u, s, force, divU, normal, w0, p0, gam1, lo, hi, ng, dx, rho_min, dt_adv, dt_divu)

     integer, intent(in) :: lo(:), hi(:), ng
     real (kind = dp_t), intent(in ) ::      u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
     real (kind = dp_t), intent(in ) ::      s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
     real (kind = dp_t), intent(in ) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)  
     real (kind = dp_t), intent(in ) ::   divU(lo(1):,lo(2):,lo(3):)
     real (kind = dp_t), intent(in ) :: normal(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
     real (kind = dp_t), intent( in) ::   w0(0:), p0(0:), gam1(0:)
     real (kind = dp_t), intent(in ) :: dx(:)
     real (kind = dp_t), intent(in ) :: rho_min
     real (kind = dp_t), intent(out) :: dt_adv, dt_divu

!    Local variables
     real (kind = dp_t), allocatable :: w0_cart(:,:,:,:)
     real (kind = dp_t)  :: spdx, spdy, spdz, spdr
     real (kind = dp_t)  :: pforcex, pforcey, pforcez
     real (kind = dp_t)  :: eps,denom,gradp0
     integer             :: i,j,k,nr

     eps = 1.0e-8

     nr = size(p0,dim=1)

     spdx    = ZERO
     spdy    = ZERO 
     spdz    = ZERO 

     allocate(w0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),3))
     call put_w0_on_3d_cells_sphr(w0(0:),w0_cart,normal,lo,hi,dx,0)

     ! Limit dt based on velocity terms
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              spdx = max(spdx ,abs(u(i,j,k,1)+w0_cart(i,j,k,1)))
              spdy = max(spdy ,abs(u(i,j,k,2)+w0_cart(i,j,k,2)))
              spdz = max(spdz ,abs(u(i,j,k,3)+w0_cart(i,j,k,3)))
           enddo
        enddo
     enddo

     deallocate(w0_cart)

     spdr = ZERO 
     do k = 0,size(w0,dim=1)-1
        spdr = max(spdr ,abs(w0(k)))
     enddo
     spdr = spdr / dr

     spdx = spdx / dx(1)
     spdy = spdy / dx(2)
     spdz = spdz / dx(3)
      
     if (spdx < eps .and. spdy < eps .and. spdz < eps .and. spdr < eps) then
        dt_adv = min(dx(1),dx(2),dx(3))
     else
        dt_adv = 1.0D0  / max(spdx,spdy,spdz,spdr)
     endif
     
     ! Limit dt based on forcing terms
     pforcex = ZERO 
     pforcey = ZERO 
     pforcez = ZERO 

     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              pforcex = max(pforcex,abs(force(i,j,k,1)))
              pforcey = max(pforcey,abs(force(i,j,k,2)))
              pforcez = max(pforcez,abs(force(i,j,k,3)))
           enddo
        enddo
     enddo
     
     if (pforcex > eps) then
        dt_adv = min(dt_adv,sqrt(2.0D0*dx(1)/pforcex))
     endif
     
     if (pforcey > eps) then
        dt_adv = min(dt_adv,sqrt(2.0D0*dx(2)/pforcey))
     endif
     
     if (pforcez > eps) then
        dt_adv = min(dt_adv,sqrt(2.0D0*dx(3)/pforcez))
     endif


     ! FIX ME FIX ME FIX ME 
     print *,'NOT FIXED IN ESTDT_3D_SPHR '
     stop
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
                              HALF*(ONE - rho_min/s(i,j,k,rho_comp))/denom)
              endif

           enddo
        enddo
     enddo

   end subroutine estdt_3d_sphr

end module estdt_module
