module mkutrans_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use slope_module

  implicit none

contains

      subroutine mkutrans_2d(vel,utrans,vtrans,force,lo,dx,dt,ng_cell,adv_bc,phys_bc)

      integer, intent(in) :: lo(2),ng_cell

      real(kind=dp_t), intent(in   ) ::     vel(lo(1)-ng_cell:,lo(2)-ng_cell:,:)
      real(kind=dp_t), intent(inout) ::  utrans(lo(1)-1:,lo(2)-1:)
      real(kind=dp_t), intent(inout) ::  vtrans(lo(1)-1:,lo(2)-1:)
      real(kind=dp_t), intent(in   ) ::   force(lo(1)-1:,lo(2)-1:,:)

      real(kind=dp_t),intent(in) :: dt,dx(:)
      integer        ,intent(in) ::  adv_bc(:,:,:)
      integer        ,intent(in) :: phys_bc(:,:  )

      real(kind=dp_t), allocatable::  velx(:,:,:)
      real(kind=dp_t), allocatable::  vely(:,:,:)

!     Local variables
      real(kind=dp_t) hx, hy, dth
      real(kind=dp_t) ulft,urgt,vbot,vtop

      real(kind=dp_t) :: eps

      integer :: hi(2)
      integer :: i,j,is,js,ie,je
      integer :: slope_order = 4

      logical :: test

      hi(1) = lo(1) + size(vel,dim=1) - (2*ng_cell+1)
      hi(2) = lo(2) + size(vel,dim=2) - (2*ng_cell+1)

      is = lo(1)
      js = lo(2)
      ie = hi(1)
      je = hi(2)

      eps = 1.0e-8              ! FIXME what should EPS really be?

      dth = HALF * dt

      hx = dx(1)
      hy = dx(2)

      allocate(velx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))
      allocate(vely(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,2))

      call slopex_2d(vel,velx,lo,ng_cell,2,adv_bc,slope_order)
      call slopey_2d(vel,vely,lo,ng_cell,2,adv_bc,slope_order)

!     Create the x-velocity to be used for transverse derivatives.
      do j = js-1,je+1 
        do i = is,ie+1 

          urgt = vel(i  ,j,1) - (HALF + dth*vel(i  ,j,1)/hx) * velx(i  ,j,1)
!    $           + dth * force(i  ,j,1)
          ulft = vel(i-1,j,1) + (HALF - dth*vel(i-1,j,1)/hx) * velx(i-1,j,1)
!    $           + dth * force(i-1,j,1)

          urgt = merge(vel(is-1,j,1),urgt,i.eq.is   .and. phys_bc(1,1) .eq. INLET)
          urgt = merge(vel(ie+1,j,1),urgt,i.eq.ie+1 .and. phys_bc(1,2) .eq. INLET)
          urgt = merge(ZERO     ,urgt,i.eq.is   .and. &
                       (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL))
          urgt = merge(ZERO     ,urgt,i.eq.ie+1 .and. &
                       (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL))

          ulft = merge(vel(is-1,j,1),ulft,i.eq.is   .and. phys_bc(1,1) .eq. INLET)
          ulft = merge(vel(ie+1,j,1),ulft,i.eq.ie+1 .and. phys_bc(1,2) .eq. INLET)
          ulft = merge(ZERO     ,ulft,i.eq.is   .and. &
                       (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL))
          ulft = merge(ZERO     ,ulft,i.eq.ie+1 .and. &
                       (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL))

          utrans(i,j) = merge(ulft,urgt,(ulft+urgt).gt.ZERO)
          test = ( (ulft .le. ZERO  .and.  urgt .ge. ZERO)  .or.  &
                   (abs(ulft+urgt) .lt. eps) )
          utrans(i,j) = merge(ZERO,utrans(i,j),test)

        enddo
      enddo

!     Create the y-velocity to be used for transverse derivatives.
      do j = js,je+1 
        do i = is-1,ie+1 

          vtop = vel(i,j  ,2) - (HALF + dth*vel(i,j  ,2)/hy) * vely(i,j  ,2)
!    $           + dth * force(i,j  ,2)
          vbot = vel(i,j-1,2) + (HALF - dth*vel(i,j-1,2)/hy) * vely(i,j-1,2)
!    $           + dth * force(i,j-1,2)

          vtop = merge(vel(i,js-1,2),vtop,j.eq.js   .and. phys_bc(2,1) .eq. INLET)
          vtop = merge(vel(i,je+1,2),vtop,j.eq.je+1 .and. phys_bc(2,2) .eq. INLET)
          vtop = merge(ZERO     ,vtop,j.eq.js   .and. &
                       (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL))
          vtop = merge(ZERO     ,vtop,j.eq.je+1 .and. &
                       (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL))

          vbot = merge(vel(i,js-1,2),vbot,j.eq.js   .and. phys_bc(2,1) .eq. INLET)
          vbot = merge(vel(i,je+1,2),vbot,j.eq.je+1 .and. phys_bc(2,2) .eq. INLET)
          vbot = merge(ZERO     ,vbot,j.eq.js   .and. &
                       (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL))
          vbot = merge(ZERO     ,vbot,j.eq.je+1 .and. &
                       (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL))

          vtrans(i,j)=merge(vbot,vtop,(vbot+vtop).gt.ZERO)
          test = ( (vbot .le. ZERO  .and.  vtop .ge. ZERO)  .or.  &
                   (abs(vbot+vtop) .lt. eps))
          vtrans(i,j) = merge(ZERO,vtrans(i,j),test)
        enddo
      enddo

      end subroutine mkutrans_2d

      subroutine mkutrans_3d(vel,utrans,vtrans,wtrans,force,lo,dx,dt,ng_cell,adv_bc,phys_bc)

      integer, intent(in) :: lo(3),ng_cell

      real(kind=dp_t), intent(in   ) ::    vel(lo(1)-ng_cell:,lo(2)-ng_cell:,lo(3)-ng_cell:,:)
      real(kind=dp_t), intent(in   ) ::  force(lo(1)-      1:,lo(2)-      1:,lo(3)-      1:,:)
      real(kind=dp_t), intent(inout) :: utrans(lo(1)-      1:,lo(2)-      1:,lo(3)-      1:)
      real(kind=dp_t), intent(inout) :: vtrans(lo(1)-      1:,lo(2)-      1:,lo(3)-      1:)
      real(kind=dp_t), intent(inout) :: wtrans(lo(1)-      1:,lo(2)-      1:,lo(3)-      1:)

      real(kind=dp_t),intent(in) :: dt,dx(:)
      integer        ,intent(in) ::  adv_bc(:,:,:)
      integer        ,intent(in) :: phys_bc(:,:  )

      real(kind=dp_t), allocatable::  velx(:,:,:,:),vely(:,:,:,:),velz(:,:,:,:)

!     Local variables
      real(kind=dp_t) ulft,urgt,vbot,vtop,wbot,wtop
      real(kind=dp_t) hx, hy, hz, dth

      real(kind=dp_t) :: eps

      logical :: test

      integer :: hi(3)
      integer :: i,j,k,is,js,ks,ie,je,ke

      integer :: slope_order = 4

      hi(1) = lo(1) + size(vel,dim=1) - (2*ng_cell+1)
      hi(2) = lo(2) + size(vel,dim=2) - (2*ng_cell+1)
      hi(3) = lo(3) + size(vel,dim=3) - (2*ng_cell+1)
 
      is = lo(1)
      js = lo(2)
      ks = lo(3)
      ie = hi(1)
      je = hi(2)
      ke = hi(3)

      eps = 1.0e-8              ! FIXME what should EPS really be?

      dth = HALF * dt

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      allocate(velx(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
      allocate(vely(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))
      allocate(velz(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,3))

      do k = lo(3)-1,hi(3)+1
         call slopex_2d(vel(:,:,k,:),velx(:,:,k,:),lo,ng_cell,3,adv_bc,slope_order)
         call slopey_2d(vel(:,:,k,:),vely(:,:,k,:),lo,ng_cell,3,adv_bc,slope_order)
      end do
      call slopez_3d(vel,velz,lo,ng_cell,  3,adv_bc,    slope_order)

!     Create the x-velocity to be used for transverse derivatives.
      do k = ks-1,ke+1
      do j = js-1,je+1
        do i = is,ie+1

          urgt = vel(i,j,k  ,1) - (HALF + dth*vel(i  ,j,k,1)/hx) * velx(i  ,j,k,1)
!    $           + dth * force(i  ,j,k,1)
          ulft = vel(i-1,j,k,1) + (HALF - dth*vel(i-1,j,k,1)/hx) * velx(i-1,j,k,1)
!    $           + dth * force(i-1,j,k,1)

          urgt = merge(vel(is-1,j,k,1),urgt,i.eq.is   .and. phys_bc(1,1) .eq. INLET)
          urgt = merge(vel(ie+1,j,k,1),urgt,i.eq.ie+1 .and. phys_bc(1,2) .eq. INLET)
          urgt = merge(ZERO           ,urgt,i.eq.is   .and. &
                       (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL))
          urgt = merge(ZERO           ,urgt,i.eq.ie+1 .and. &
                       (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL))

          ulft = merge(vel(is-1,j,k,1),ulft,i.eq.is   .and. phys_bc(1,1) .eq. INLET)
          ulft = merge(vel(ie+1,j,k,1),ulft,i.eq.ie+1 .and. phys_bc(1,2) .eq. INLET)
          ulft = merge(ZERO           ,ulft,i.eq.is   .and. &
                       (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL))
          ulft = merge(ZERO           ,ulft,i.eq.ie+1 .and. &
                       (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL))

          utrans(i,j,k) = merge(ulft,urgt,(ulft+urgt).gt.ZERO)
          test=( (ulft .le. ZERO  .and.  urgt .ge. ZERO)  .or. &
                (abs(ulft+urgt) .lt. eps) )
          utrans(i,j,k) = merge(ZERO,utrans(i,j,k),test)
        enddo
        enddo
      enddo

!     Create the y-velocity to be used for transverse derivatives.
      do j = js,je+1
        do k = ks-1,ke+1
        do i = is-1,ie+1

          vtop = vel(i,j  ,k,2) - (HALF + dth*vel(i,j  ,k,2)/hy) * vely(i,j  ,k,2)
!    $           + dth * force(i,j  ,k,2)
          vbot = vel(i,j-1,k,2) + (HALF - dth*vel(i,j-1,k,2)/hy) * vely(i,j-1,k,2)
!    $           + dth * force(i,j-1,k,2)

          vtop = merge(vel(i,js-1,k,2),vtop,j.eq.js   .and. phys_bc(2,1) .eq. INLET)
          vtop = merge(vel(i,je+1,k,2),vtop,j.eq.je+1 .and. phys_bc(2,2) .eq. INLET)
          vtop = merge(ZERO           ,vtop,j.eq.js   .and. &
                       (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL))
          vtop = merge(ZERO           ,vtop,j.eq.je+1 .and. &
                       (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL))

          vbot = merge(vel(i,js-1,k,2),vbot,j.eq.js   .and. phys_bc(2,1) .eq. INLET)
          vbot = merge(vel(i,je+1,k,2),vbot,j.eq.je+1 .and. phys_bc(2,2) .eq. INLET)
          vbot = merge(ZERO           ,vbot,j.eq.js   .and. &
                       (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL))
          vbot = merge(ZERO           ,vbot,j.eq.je+1 .and. &
                       (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL))

          vtrans(i,j,k)=merge(vbot,vtop,(vbot+vtop).gt.ZERO)
          test = ( (vbot .le. ZERO  .and.  vtop .ge. ZERO)  .or. &
                   (abs(vbot+vtop) .lt. eps))
          vtrans(i,j,k) = merge(ZERO,vtrans(i,j,k),test)

        enddo
        enddo
      enddo

!     Create the z-velocity to be used for transverse derivatives.
      do k = ks,ke+1
        do j = js-1,je+1
        do i = is-1,ie+1

          wtop = vel(i,j,k  ,3) - (HALF + dth*vel(i,j,k  ,3)/hz) * velz(i,j,k  ,3)
!    $           + dth * force(i,j,k  ,3)
          wbot = vel(i,j,k-1,3) + (HALF - dth*vel(i,j,k-1,3)/hz) * velz(i,j,k-1,3)
!    $           + dth * force(i,j,k-1,3)

          wtop = merge(vel(i,j,ks-1,3),wtop,k.eq.ks   .and. phys_bc(3,1) .eq. INLET)
          wtop = merge(vel(i,j,ke+1,3),wtop,k.eq.ke+1 .and. phys_bc(3,2) .eq. INLET)
          wtop = merge(ZERO           ,wtop,k.eq.ks   .and. &
                       (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL))
          wtop = merge(ZERO           ,wtop,k.eq.ke+1 .and. &
                       (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL))

          wbot = merge(vel(i,j,ks-1,3),wbot,k.eq.ks   .and. phys_bc(3,1) .eq. INLET)
          wbot = merge(vel(i,j,ke+1,3),wbot,k.eq.ke+1 .and. phys_bc(3,2) .eq. INLET)
          wbot = merge(ZERO           ,wbot,k.eq.ks   .and. &
                       (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL))
          wbot = merge(ZERO           ,wbot,k.eq.ke+1 .and. &
                       (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL))

          wtrans(i,j,k)=merge(wbot,wtop,(wbot+wtop).gt.ZERO)
          test = ( (wbot .le. ZERO  .and.  wtop .ge. ZERO)  .or. &
                   (abs(wbot+wtop) .lt. eps))
          wtrans(i,j,k) = merge(ZERO,wtrans(i,j,k),test)

        enddo
        enddo
      enddo

      end subroutine mkutrans_3d

end module mkutrans_module
