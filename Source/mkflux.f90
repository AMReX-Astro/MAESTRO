module mkflux_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use slope_module
  use fill_3d_module
  use geometry

  implicit none

contains

      subroutine mkflux_2d(s,u,sedgex,sedgey,uadv,vadv,utrans,vtrans,&
                           force,w0,lo,dx,dt,is_vel,is_cons,&
                           phys_bc,adv_bc,velpred,ng,base, &
                           advect_in_pert_form,n)

      integer, intent(in) :: lo(:),ng,n

      real(kind=dp_t), intent(inout) ::      s(lo(1)-ng:,lo(2)-ng:,:)
      real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,:)
      real(kind=dp_t), intent(inout) :: sedgex(lo(1)   :,lo(2)   :,:)
      real(kind=dp_t), intent(inout) :: sedgey(lo(1)   :,lo(2)   :,:)
      real(kind=dp_t), intent(inout) ::   uadv(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(inout) ::   vadv(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(in   ) :: utrans(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(in   ) :: vtrans(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t), intent(inout) ::  force(lo(1)- 1:,lo(2)- 1:,:)
      real(kind=dp_t), intent(in   ) ::     w0(0:)

      real(kind=dp_t),intent(in) :: dt,dx(:),base(0:)
      integer        ,intent(in) :: velpred
      integer        ,intent(in) :: phys_bc(:,:)
      integer        ,intent(in) ::  adv_bc(:,:,:)
      logical        ,intent(in) :: is_vel
      logical        ,intent(in) :: is_cons(:)
      logical        ,intent(in) :: advect_in_pert_form

      real(kind=dp_t), allocatable::  slopex(:,:,:),slopey(:,:,:)
      real(kind=dp_t), allocatable::  s_l(:),s_r(:),s_b(:),s_t(:)

!     Local variables
      real(kind=dp_t) ubardth, vbardth
      real(kind=dp_t) hx, hy, dth
      real(kind=dp_t) splus,sminus
      real(kind=dp_t) savg,st
      real(kind=dp_t) vlo,vhi
      real(kind=dp_t) sptop,spbot,smtop,smbot,splft,sprgt,smlft,smrgt

      integer :: hi(2)
      integer :: slope_order = 4
      logical :: test

      real(kind=dp_t) :: abs_eps, eps, umax
      real(kind=dp_t) :: vadv_max

      integer :: i,j,is,js,ie,je,g,nr

      hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
      hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)

      is = lo(1)
      ie = hi(1)
      js = lo(2)
      je = hi(2)

      nr = size(w0,dim=1)-1

      allocate(s_l(lo(1)-1:hi(1)+2))
      allocate(s_r(lo(1)-1:hi(1)+2))
      allocate(s_b(lo(2)-1:hi(2)+2))
      allocate(s_t(lo(2)-1:hi(2)+2))

      allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1))
      allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,1))

      if (.not. is_vel .and. advect_in_pert_form) then
         do j = js,je
           do i = is-ng,ie+ng
             s(i,j,n) = s(i,j,n) - base(j)
           end do
           do g = 1,ng
            do i = is-ng,ie+ng
             s(i,js-g,n) = s(i,js,n)
             s(i,je+g,n) = s(i,je,n)
           end do
          end do
         end do
      end if


      call slopex_2d(s(:,:,n:),slopex,lo,ng,1,adv_bc,slope_order)
      call slopey_2d(s(:,:,n:),slopey,lo,ng,1,adv_bc,slope_order)

      abs_eps = 1.0e-8

      dth = HALF*dt

      hx = dx(1)
      hy = dx(2)

      if (velpred .eq. 1) then

        umax = abs(utrans(is,js))
        do j = js,je
           do i = is,ie+1
             umax = max(umax,abs(utrans(i,j)))
           end do
        end do
        do j = js,je+1
           do i = is,ie
             umax = max(umax,abs(vtrans(i,j)))
           end do
        end do

      else 

        umax = abs(uadv(is,js))
        do j = js,je
           do i = is,ie+1
             umax = max(umax,abs(uadv(i,j)))
           end do
        end do
        do j = js,je+1
           do i = is,ie
             umax = max(umax,abs(vadv(i,j)))
           end do
        end do
      end if

      eps = abs_eps * umax

      ! HACK -- NOTE THAT EPS >= 1.0 FOR VELPRED == 1.
      if (velpred .eq. 1) then
        eps = max(1.d0,eps)
      end if

!
!     Loop for fluxes on x-edges.
!
      vadv_max = ZERO

       do j = js,je 
        if (velpred .eq. 0 .or. n .eq. 1) then
        do i = is-1,ie+1 
 
          vlo = u(i,j  ,2) + HALF * (w0(j  )+w0(j+1))
          if ((j+2).le.nr) then
            vhi = u(i,j+1,2) + HALF * (w0(j+1)+w0(j+2))
          else
            vhi = u(i,j+1,2) + w0(j+1)
          end if

          spbot = s(i,j  ,n) + (HALF - dth*vlo/hy) * slopey(i,j  ,1)
          sptop = s(i,j+1,n) - (HALF + dth*vhi/hy) * slopey(i,j+1,1)

          sptop = merge(s(i,je+1,n),sptop,j.eq.je .and. phys_bc(2,2) .eq. INLET)
          spbot = merge(s(i,je+1,n),spbot,j.eq.je .and. phys_bc(2,2) .eq. INLET)

          if (j .eq. je .and. (phys_bc(2,2).eq.SLIP_WALL.or.phys_bc(2,2).eq.NO_SLIP_WALL)) then
            if (is_vel .and. n .eq. 2) then
              sptop = ZERO
              spbot = ZERO
            elseif (is_vel .and. n .eq. 1) then
              sptop = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
              spbot = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
            else
              sptop = spbot
            endif
          endif

          splus = merge(spbot,sptop,vtrans(i,j+1).gt.ZERO)
          savg  = HALF * (spbot + sptop)
          splus = merge(splus, savg, abs(vtrans(i,j+1)) .gt. eps)

          if (j.ge.1) then
            vlo = u(i,j-1,2) + HALF * (w0(j-1)+w0(j  ))
          else
            vlo = u(i,j-1,2) + w0(j)
          end if
          vhi = u(i,j  ,2) + HALF * (w0(j  )+w0(j+1))

          smtop = s(i,j  ,n) - (HALF + dth*vhi/hy) * slopey(i,j  ,1)
          smbot = s(i,j-1,n) + (HALF - dth*vlo/hy) * slopey(i,j-1,1)

          smtop = merge(s(i,js-1,n),smtop,j.eq.js .and. phys_bc(2,1) .eq. INLET)
          smbot = merge(s(i,js-1,n),smbot,j.eq.js .and. phys_bc(2,1) .eq. INLET)

          if (j .eq. js .and. (phys_bc(2,1).eq.SLIP_WALL.or.phys_bc(2,1).eq.NO_SLIP_WALL)) then
            if (is_vel .and. (n .eq. 2)) then
              smtop = ZERO
              smbot = ZERO
            elseif (is_vel .and. (n .ne. 2)) then
              smbot = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
              smtop = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
            else
              smbot = smtop
            endif
          endif

          sminus = merge(smbot,smtop,vtrans(i,j).gt.ZERO)
          savg   = HALF * (smbot + smtop)
          sminus = merge(sminus, savg, abs(vtrans(i,j)) .gt. eps)

          st = force(i,j,n) - &
                HALF * (vtrans(i,j)+vtrans(i,j+1))*(splus - sminus) / hy

          if (is_vel .and. n.eq.2) then
            st = st - HALF * (vtrans(i,j)+vtrans(i,j+1))*(w0(j+1)-w0(j))/hy
          end if

          ubardth = dth*u(i,j,1)/hx

          s_l(i+1)= s(i,j,n) + (HALF-ubardth)*slopex(i,j,1) + dth*st
          s_r(i  )= s(i,j,n) - (HALF+ubardth)*slopex(i,j,1) + dth*st

         enddo

         if (velpred .eq. 1) then
           do i = is, ie+1 
             savg = HALF*(s_r(i) + s_l(i))
             test = ( (s_l(i) .le. ZERO  .and. &
                       s_r(i) .ge. ZERO)  .or. &
                     (abs(s_l(i) + s_r(i)) .lt. eps) )
             sedgex(i,j,n)=merge(s_l(i),s_r(i),savg.gt.ZERO)
             sedgex(i,j,n)=merge(savg,sedgex(i,j,n),test)
           enddo
         else
           do i = is, ie+1 
             sedgex(i,j,n)=merge(s_l(i),s_r(i),uadv(i,j).gt.ZERO)
             savg = HALF*(s_r(i) + s_l(i))
             sedgex(i,j,n)=merge(savg,sedgex(i,j,n),abs(uadv(i,j)) .lt. eps)
           enddo
         endif

         if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
           if (is_vel .and. n .eq. 1) then
             sedgex(is,j,n) = ZERO
           elseif (is_vel .and. n .ne. 1) then
             sedgex(is,j,n) = merge(ZERO,s_r(is),phys_bc(1,1).eq.NO_SLIP_WALL)
           else 
             sedgex(is,j,n) = s_r(is)
           endif
         elseif (phys_bc(1,1) .eq. INLET) then
           sedgex(is,j,n) = s(is-1,j,n)
         elseif (phys_bc(1,1) .eq. OUTLET) then
           sedgex(is,j,n) = s_r(is)
         endif
         if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
           if (is_vel .and. n .eq. 1) then
             sedgex(ie+1,j,n) = ZERO
           else if (is_vel .and. n .ne. 1) then
             sedgex(ie+1,j,n) = merge(ZERO,s_l(ie+1),phys_bc(1,2).eq.NO_SLIP_WALL)
           else 
             sedgex(ie+1,j,n) = s_l(ie+1)
           endif
         elseif (phys_bc(1,2) .eq. INLET) then
           sedgex(ie+1,j,n) = s(ie+1,j,n)
         elseif (phys_bc(1,2) .eq. OUTLET) then
           sedgex(ie+1,j,n) = s_l(ie+1)
         endif

         if (velpred .eq. 1) then
           do i = is, ie+1 
             uadv(i,j) = sedgex(i,j,1)
           enddo
         endif
         endif
       enddo

!
!     Loop for fluxes on y-edges.
!
       do i = is, ie 
        if (velpred .eq. 0 .or. n .eq. 2) then
        do j = js-1, je+1 

          splft = s(i,j  ,n) + (HALF - dth*u(i  ,j,1)/hx) * slopex(i  ,j,1)
          sprgt = s(i+1,j,n) - (HALF + dth*u(i+1,j,1)/hx) * slopex(i+1,j,1)

          sprgt = merge(s(ie+1,j,n),sprgt,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
          splft = merge(s(ie+1,j,n),splft,i.eq.ie .and. phys_bc(1,2) .eq. INLET)

          if (i .eq. ie .and. (phys_bc(1,2).eq.SLIP_WALL.or.phys_bc(1,2).eq.NO_SLIP_WALL)) then
            if (is_vel .and. n .eq. 1) then
              splft = ZERO
              sprgt = ZERO
            elseif (is_vel .and. n .ne. 1) then
              sprgt = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
              splft = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
            else
              sprgt = splft
            endif
          endif

          splus = merge(splft,sprgt,utrans(i+1,j).gt.ZERO)
          savg  = HALF * (splft + sprgt)
          splus = merge(splus, savg, abs(utrans(i+1,j)) .gt. eps)

          smrgt = s(i  ,j,n) - (HALF + dth*u(i  ,j,1)/hx) * slopex(i  ,j,1)
          smlft = s(i-1,j,n) + (HALF - dth*u(i-1,j,1)/hx) * slopex(i-1,j,1)

          smrgt = merge(s(is-1,j,n),smrgt,i.eq.is .and. phys_bc(1,1) .eq. INLET)
          smlft = merge(s(is-1,j,n),smlft,i.eq.is .and. phys_bc(1,1) .eq. INLET)

          if (i .eq. is .and. (phys_bc(1,1).eq.SLIP_WALL.or.phys_bc(1,1).eq.NO_SLIP_WALL)) then
            if (is_vel .and. n .eq. 1) then
              smlft = ZERO
              smrgt = ZERO
            elseif (is_vel .and. n .ne. 1) then
              smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
              smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
            else
              smlft = smrgt
            endif
          endif

          sminus = merge(smlft,smrgt,utrans(i,j).gt.ZERO)
          savg   = HALF * (smlft + smrgt)
          sminus = merge(sminus, savg, abs(utrans(i,j)) .gt. eps)

          st = force(i,j,n) - &
               HALF * (utrans(i,j)+utrans(i+1,j))*(splus - sminus) / hx

          if (is_vel .and. n.eq.2 .and. j.ge.0 .and. j.lt.nr) then
            st = st - HALF * (vtrans(i,j)+vtrans(i,j+1))*(w0(j+1)-w0(j))/hy
          end if

          if (j .ge. 0 .and. j.lt.nr) then
            vbardth = dth / hy * ( u(i,j,2) + HALF * (w0(j)+w0(j+1)) )
          else
            vbardth = dth / hy * u(i,j,2) 
          end if

          s_b(j+1)= s(i,j,n) + (HALF-vbardth)*slopey(i,j,1) + dth*st
          s_t(j  )= s(i,j,n) - (HALF+vbardth)*slopey(i,j,1) + dth*st
        enddo

        if (velpred .eq. 1) then
          do j = js, je+1 
            savg = HALF*(s_b(j) + s_t(j))
            test = ( (s_b(j) .le. ZERO  .and. &
                      s_t(j) .ge. ZERO)  .or. &
                   (abs(s_b(j) + s_t(j)) .lt. eps) )
            sedgey(i,j,n)=merge(s_b(j),s_t(j),savg.gt.ZERO)
            sedgey(i,j,n)=merge(savg,sedgey(i,j,n),test)
          enddo

        else

          do j = js, je+1 
            sedgey(i,j,n)=merge(s_b(j),s_t(j),vadv(i,j).gt.ZERO)
            savg = HALF*(s_b(j) + s_t(j))
            sedgey(i,j,n)=merge(savg,sedgey(i,j,n),abs(vadv(i,j)) .lt. eps)
          enddo
          ! OUTFLOW HACK 
          if (phys_bc(2,2) .eq. OUTLET .and. vadv(i,je+1).lt.ZERO ) sedgey(i,je+1,n) = s_b(je+1)
        endif

        if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
          if (is_vel .and. n .eq. 2) then
            sedgey(i,js,n) = ZERO
          elseif (is_vel .and. n .ne. 2) then
            sedgey(i,js,n) = merge(ZERO,s_t(js),phys_bc(2,1).eq.NO_SLIP_WALL)
          else 
            sedgey(i,js,n) = s_t(js)
          endif
        elseif (phys_bc(2,1) .eq. INLET) then
          sedgey(i,js,n) = s(i,js-1,n)
        elseif (phys_bc(2,1) .eq. OUTLET) then
          sedgey(i,js,n) = s_t(js)
        endif

        if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
          if (is_vel .and. n .eq. 2) then
            sedgey(i,je+1,n) = ZERO
          elseif (is_vel .and. n .ne. 2) then
            sedgey(i,je+1,n) = merge(ZERO,s_b(je+1),phys_bc(2,2).eq.NO_SLIP_WALL)
          else 
            sedgey(i,je+1,n) = s_b(je+1)
          endif
        elseif (phys_bc(2,2) .eq. INLET) then
          sedgey(i,je+1,n) = s(i,je+1,n)
        elseif (phys_bc(2,2) .eq. OUTLET) then
          sedgey(i,je+1,n) = s_b(je+1)
        endif

        if (velpred .eq. 1) then
          do j = js, je+1 
            vadv(i,j) = sedgey(i,j,2)
            vadv_max = max(vadv_max,abs(sedgey(i,j,2)))
          enddo
        end if

        endif

       enddo

      if (.not. is_vel .and. advect_in_pert_form) then
         do j = js,je
           do i = is-ng,ie+ng
             s(i,j,n) = s(i,j,n) + base(j)
           end do
           do g = 1,ng
            do i = is-ng,ie+ng
             s(i,js-g,n) = s(i,js,n)
           end do
          end do
         end do

!        do j = js,je
!          do i = is,ie+1 
!            sedgex(i,j,n) = sedgex(i,j,n) + base(j)
!          enddo
!        enddo 

!        do i = is,ie
!          do j = js+1,je
!            sedgey(i,j,n) = sedgey(i,j,n) + HALF*(base(j)+base(j-1))
!          enddo
!          sedgey(i,js  ,n) = sedgey(i,js  ,n) + base(js)
!          sedgey(i,je+1,n) = sedgey(i,je+1,n) + base(je)
!        enddo

      end if

      deallocate(s_l)
      deallocate(s_r)
      deallocate(s_b)
      deallocate(s_t)

      deallocate(slopex)
      deallocate(slopey)

      end subroutine mkflux_2d


      subroutine mkflux_3d(s,u,sedgex,sedgey,sedgez,uadv,vadv,wadv,utrans,vtrans,wtrans,&
                           force,w0,w0_cart_vec,lo,dx,dt,is_vel,is_cons,&
                           phys_bc,adv_bc,velpred,ng,base, &
                           advect_in_pert_form,n)

      integer, intent(in) :: lo(:),ng,n

      real(kind=dp_t), intent(inout) ::      s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real(kind=dp_t), intent(in   ) ::      u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real(kind=dp_t), intent(inout) :: sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real(kind=dp_t), intent(inout) :: sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real(kind=dp_t), intent(inout) :: sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real(kind=dp_t), intent(inout) ::   uadv(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
      real(kind=dp_t), intent(inout) ::   vadv(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
      real(kind=dp_t), intent(inout) ::   wadv(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
      real(kind=dp_t), intent(in   ) :: utrans(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
      real(kind=dp_t), intent(in   ) :: vtrans(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
      real(kind=dp_t), intent(in   ) :: wtrans(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
      real(kind=dp_t), intent(inout) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
      real(kind=dp_t), intent(in   ) ::     w0(0:)
      real(kind=dp_t), intent(in   ) :: w0_cart_vec(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)

      real(kind=dp_t),intent(in) :: dt,dx(:),base(0:)
      integer        ,intent(in) :: velpred
      integer        ,intent(in) :: phys_bc(:,:)
      integer        ,intent(in) ::  adv_bc(:,:,:)
      logical        ,intent(in) :: is_vel
      logical        ,intent(in) :: is_cons(:)
      logical        ,intent(in) :: advect_in_pert_form

      real(kind=dp_t), allocatable :: slopex(:,:,:,:),slopey(:,:,:,:),slopez(:,:,:,:)
      real(kind=dp_t), allocatable :: s_l(:),s_r(:),s_b(:),s_t(:),s_u(:),s_d(:)
      real(kind=dp_t), allocatable :: base_cart(:,:,:)

!     Local variables
      real(kind=dp_t) ubardth, vbardth, wbardth
      real(kind=dp_t) hx, hy, hz, dth
      real(kind=dp_t) splus,sminus
      real(kind=dp_t) savg,st
      real(kind=dp_t) ulo,uhi,vlo,vhi,wlo,whi
      real(kind=dp_t) sptop,spbot,smtop,smbot,splft,sprgt,smlft,smrgt

      integer :: hi(3)
      integer :: slope_order = 4
      logical :: test

      real(kind=dp_t) :: abs_eps, eps, umax, w0cell

      integer :: i,j,k,is,js,ie,je,ks,ke,g

      hi(1) = lo(1) + size(s,dim=1) - (2*ng+1)
      hi(2) = lo(2) + size(s,dim=2) - (2*ng+1)
      hi(3) = lo(3) + size(s,dim=3) - (2*ng+1)

      is = lo(1)
      ie = hi(1)
      js = lo(2)
      je = hi(2)
      ks = lo(3)
      ke = hi(3)

      allocate(s_l(lo(1)-1:hi(1)+2))
      allocate(s_r(lo(1)-1:hi(1)+2))
      allocate(s_b(lo(2)-1:hi(2)+2))
      allocate(s_t(lo(2)-1:hi(2)+2))
      allocate(s_d(lo(3)-1:hi(3)+2))
      allocate(s_u(lo(3)-1:hi(3)+2))

      allocate(slopex(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))
      allocate(slopey(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))
      allocate(slopez(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,1))

      if (.not. is_vel .and. advect_in_pert_form) then
         if (spherical .eq. 1) then
           allocate(base_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
           call fill_3d_data(base_cart,base,lo,hi,dx,0)
           do k = ks,ke
             do j = js,je
             do i = is,ie
               s(i,j,k,n) = s(i,j,k,n) - base_cart(i,j,k)
             end do
             end do
             do g = 1,ng
               do j = js,je
                 s(is-g,j,k,n) = s(is,j,k,n)
                 s(ie+g,j,k,n) = s(ie,j,k,n)
               end do
               do i = is-1,ie+1
                 s(i,js-g,k,n) = s(i,j,k,n)
                 s(i,je+g,k,n) = s(i,je,k,n)
               end do
             end do
           end do
         else
           do k = ks,ke
             do j = js-ng,je+ng
             do i = is-ng,ie+ng
               s(i,j,k,n) = s(i,j,k,n) - base(k)
             end do
             end do
           end do
         end if

        do k = ks,ke
          do g = 1,ng
           do j = js-ng,je+ng
           do i = is-ng,ie+ng
            s(i,j,ks-g,n) = s(i,j,ks,n)
            s(i,j,ke+g,n) = s(i,j,ke,n)
          end do
          end do
         end do
        end do
      end if

      do k = lo(3)-1,hi(3)+1
         call slopex_2d(s(:,:,k,n:),slopex(:,:,k,:),lo,ng,1,adv_bc,slope_order)
         call slopey_2d(s(:,:,k,n:),slopey(:,:,k,:),lo,ng,1,adv_bc,slope_order)
      end do
      call slopez_3d(s(:,:,:,n:),slopez,lo,ng,1,adv_bc,slope_order)

      abs_eps = 1.0e-8

      dth = HALF*dt

      hx = dx(1)
      hy = dx(2)
      hz = dx(3)

      if (velpred .eq. 1) then

        umax = abs(utrans(is,js,ks))
        do k = ks,ke
        do j = js,je
        do i = is,ie+1
           umax = max(umax,abs(utrans(i,j,k)))
        end do
        end do
        end do
        do k = ks,ke
        do j = js,je+1
        do i = is,ie
          umax = max(umax,abs(vtrans(i,j,k)))
        end do
        end do
        end do
        do k = ks,ke+1
        do j = js,je
        do i = is,ie
          umax = max(umax,abs(wtrans(i,j,k)))
        end do
        end do
        end do

      else 

        umax = abs(uadv(is,js,ks))
        do k = ks,ke
        do j = js,je
        do i = is,ie+1
          umax = max(umax,abs(uadv(i,j,k)))
        end do
        end do
        end do
        do k = ks,ke
        do j = js,je+1
        do i = is,ie
          umax = max(umax,abs(vadv(i,j,k)))
        end do
        end do
        end do
        do k = ks,ke+1
        do j = js,je
        do i = is,ie
          umax = max(umax,abs(wadv(i,j,k)))
        end do
        end do
        end do
      end if

      eps = abs_eps * umax

      ! HACK -- NOTE THAT EPS >= 1.0 FOR VELPRED == 1.
      if (velpred .eq. 1) then
        eps = max(1.d0,eps)
      end if
!
!     Loop for fluxes on x-edges.
!
      if (velpred .eq. 0 .or. n .eq. 1) then
       do k = ks,ke 
       do j = js,je 
        do i = is-1,ie+1 

          ! Do transverse in j direction

          vlo = u(i,j  ,k,2) + w0_cart_vec(i,j  ,k,2)
          vhi = u(i,j+1,k,2) + w0_cart_vec(i,j+1,k,2)

          spbot = s(i,j  ,k,n) + (HALF - dth*vlo/hy) * slopey(i,j  ,k,1)
          sptop = s(i,j+1,k,n) - (HALF + dth*vhi/hy) * slopey(i,j+1,k,1)

          sptop = merge(s(i,je+1,k,n),sptop,j.eq.je .and. phys_bc(2,2) .eq. INLET)
          spbot = merge(s(i,je+1,k,n),spbot,j.eq.je .and. phys_bc(2,2) .eq. INLET)

          if (j .eq. je .and. (phys_bc(2,2).eq.SLIP_WALL.or.phys_bc(2,2).eq.NO_SLIP_WALL)) then
            if (is_vel .and. n.eq.2) then
              sptop = ZERO
              spbot = ZERO
            elseif (is_vel .and. n.ne.2) then
              sptop = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
              spbot = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
            else
              sptop = spbot
            endif
          endif

          splus = merge(spbot,sptop,vtrans(i,j+1,k).gt.ZERO)
          savg  = HALF * (spbot + sptop)
          splus = merge(splus, savg, abs(vtrans(i,j+1,k)) .gt. eps)

          vlo = u(i,j-1,k,2) + w0_cart_vec(i,j-1,k,2)
          vhi = u(i,j  ,k,2) + w0_cart_vec(i,j  ,k,2)

          smbot = s(i,j-1,k,n) + (HALF - dth*vlo/hy) * slopey(i,j-1,k,1)
          smtop = s(i,j  ,k,n) - (HALF + dth*vhi/hy) * slopey(i,j  ,k,1)

          smtop = merge(s(i,js-1,k,n),smtop,j.eq.js .and. phys_bc(2,1) .eq. INLET)
          smbot = merge(s(i,js-1,k,n),smbot,j.eq.js .and. phys_bc(2,1) .eq. INLET)

          if (j .eq. js .and. (phys_bc(2,1).eq.SLIP_WALL.or.phys_bc(2,1).eq.NO_SLIP_WALL)) then
            if (is_vel .and. (n .eq. 2)) then
              smtop = ZERO
              smbot = ZERO
            elseif (is_vel .and. (n .ne. 2)) then
              smbot = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
              smtop = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
            else
              smbot = smtop
            endif
          endif

          sminus = merge(smbot,smtop,vtrans(i,j,k).gt.ZERO)
          savg   = HALF * (smbot + smtop)
          sminus = merge(sminus, savg, abs(vtrans(i,j,k)) .gt. eps)

          st = force(i,j,k,n) - &
                HALF * (vtrans(i,j,k)+vtrans(i,j+1,k))*(splus - sminus) / hy

          ! Do transverse in k direction

          wlo = u(i,j,k  ,3) + w0_cart_vec(i,j,k  ,3)
          whi = u(i,j,k+1,3) + w0_cart_vec(i,j,k+1,3)

          spbot = s(i,j,k  ,n) + (HALF - dth*wlo/hz) * slopez(i,j,k  ,1)
          sptop = s(i,j,k+1,n) - (HALF + dth*whi/hz) * slopez(i,j,k+1,1)

          sptop = merge(s(i,j,ke+1,n),sptop,k.eq.ke .and. phys_bc(3,2) .eq. INLET)
          spbot = merge(s(i,j,ke+1,n),spbot,k.eq.ke .and. phys_bc(3,2) .eq. INLET)

          if (k .eq. ke .and. (phys_bc(3,2).eq.SLIP_WALL.or.phys_bc(3,2).eq.NO_SLIP_WALL)) then
            if (is_vel .and. n .eq. 3) then
              sptop = ZERO
              spbot = ZERO
            elseif (is_vel .and. (n.ne.3)) then
              sptop = merge(ZERO,spbot,phys_bc(3,2).eq.NO_SLIP_WALL)
              spbot = merge(ZERO,spbot,phys_bc(3,2).eq.NO_SLIP_WALL)
            else
              sptop = spbot
            endif
          endif

          splus = merge(spbot,sptop,wtrans(i,j,k+1).gt.ZERO)
          savg  = HALF * (spbot + sptop)
          splus = merge(splus, savg, abs(wtrans(i,j,k+1)) .gt. eps)

          wlo = u(i,j,k-1,3) + w0_cart_vec(i,j,k-1,3)
          whi = u(i,j,k  ,3) + w0_cart_vec(i,j,k  ,3)

          smtop = s(i,j,k  ,n) - (HALF + dth*whi/hz) * slopez(i,j,k  ,1)
          smbot = s(i,j,k-1,n) + (HALF - dth*wlo/hz) * slopez(i,j,k-1,1)

          smtop = merge(s(i,j,ks-1,n),smtop,k.eq.ks .and. phys_bc(3,1) .eq. INLET)
          smbot = merge(s(i,j,ks-1,n),smbot,k.eq.ks .and. phys_bc(3,1) .eq. INLET)

          if (k .eq. ks .and. (phys_bc(3,1).eq.SLIP_WALL.or.phys_bc(3,1).eq.NO_SLIP_WALL)) then
            if (is_vel .and. (n.eq.3)) then
              smtop = ZERO
              smbot = ZERO
            elseif (is_vel .and. (n.ne.3)) then
              smbot = merge(ZERO,smtop,phys_bc(3,1).eq.NO_SLIP_WALL)
              smtop = merge(ZERO,smtop,phys_bc(3,1).eq.NO_SLIP_WALL)
            else
              smbot = smtop
            endif
          endif

          sminus = merge(smbot,smtop,wtrans(i,j,k).gt.ZERO)
          savg   = HALF * (smbot + smtop)
          sminus = merge(sminus, savg, abs(wtrans(i,j,k)) .gt. eps)

          st = st - HALF * (wtrans(i,j,k)+wtrans(i,j+1,k))*(splus - sminus) / hz

          ! NOTE NOTE : THIS IS WRONG FOR SPHERICAL !!
          if (spherical .eq. 0 .and. is_vel .and. n.eq.3) then
            st = st - HALF * (wtrans(i,j,k)+wtrans(i,j,k+1))*(w0(k+1)-w0(k))/hz
          end if

          ubardth = dth/hx * ( u(i,j,k,1) + w0_cart_vec(i,j,k,1))

          s_l(i+1)= s(i,j,k,n) + (HALF-ubardth)*slopex(i,j,k,1) + dth*st
          s_r(i  )= s(i,j,k,n) - (HALF+ubardth)*slopex(i,j,k,1) + dth*st

         enddo

         if (velpred .eq. 1) then
           do i = is, ie+1 
             savg = HALF*(s_r(i) + s_l(i))
             test = ( (s_l(i) .le. ZERO  .and. &
                       s_r(i) .ge. ZERO)  .or. &
                     (abs(s_l(i) + s_r(i)) .lt. eps) )
             sedgex(i,j,k,n)=merge(s_l(i),s_r(i),savg.gt.ZERO)
             sedgex(i,j,k,n)=merge(savg,sedgex(i,j,k,n),test)
           enddo
         else
           do i = is, ie+1 
             sedgex(i,j,k,n)=merge(s_l(i),s_r(i),uadv(i,j,k).gt.ZERO)
             savg = HALF*(s_r(i) + s_l(i))
             sedgex(i,j,k,n)=merge(savg,sedgex(i,j,k,n),abs(uadv(i,j,k)) .lt. eps)
           enddo
         endif

         if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_WALL) then
           if (is_vel .and. n .eq. 1) then
             sedgex(is,j,k,n) = ZERO
           elseif (is_vel .and. n .ne. 1) then
             sedgex(is,j,k,n) = merge(ZERO,s_r(is),phys_bc(1,1).eq.NO_SLIP_WALL)
           else 
             sedgex(is,j,k,n) = s_r(is)
           endif
         elseif (phys_bc(1,1) .eq. INLET) then
           sedgex(is,j,k,n) = s(is-1,j,k,n)
         elseif (phys_bc(1,1) .eq. OUTLET) then
           sedgex(is,j,k,n) = s_r(is)
         endif
         if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_WALL) then
           if (is_vel .and. n .eq. 1) then
             sedgex(ie+1,j,k,n) = ZERO
           else if (is_vel .and. n .ne. 1) then
             sedgex(ie+1,j,k,n) = merge(ZERO,s_l(ie+1),phys_bc(1,2).eq.NO_SLIP_WALL)
           else 
             sedgex(ie+1,j,k,n) = s_l(ie+1)
           endif
         elseif (phys_bc(1,2) .eq. INLET) then
           sedgex(ie+1,j,k,n) = s(ie+1,j,k,n)
         elseif (phys_bc(1,2) .eq. OUTLET) then
           sedgex(ie+1,j,k,n) = s_l(ie+1)
         endif

         if (velpred .eq. 1) then
           do i = is, ie+1 
             uadv(i,j,k) = sedgex(i,j,k,1)
           enddo
         endif

         enddo
         enddo
       endif
!
!     Loop for fluxes on y-edges.
!
       if (velpred .eq. 0 .or. n .eq. 2) then
       do k = ks, ke 
       do i = is, ie 
        do j = js-1, je+1 

          ! Do transverse in i direction

          ulo = u(i  ,j,k,1) + w0_cart_vec(i  ,j,k,1)
          uhi = u(i+1,j,k,1) + w0_cart_vec(i+1,j,k,1)

          splft = s(i  ,j,k,n) + (HALF - dth*ulo/hx) * slopex(i  ,j,k,1)
          sprgt = s(i+1,j,k,n) - (HALF + dth*uhi/hx) * slopex(i+1,j,k,1)

          sprgt = merge(s(ie+1,j,k,n),sprgt,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
          splft = merge(s(ie+1,j,k,n),splft,i.eq.ie .and. phys_bc(1,2) .eq. INLET)

          if (i .eq. ie .and. (phys_bc(1,2).eq.SLIP_WALL.or.phys_bc(1,2).eq.NO_SLIP_WALL)) then
            if (is_vel .and. n .eq. 1) then
              splft = ZERO
              sprgt = ZERO
            elseif (is_vel .and. n .ne. 1) then
              sprgt = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
              splft = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
            else
              sprgt = splft
            endif
          endif

          splus = merge(splft,sprgt,utrans(i+1,j,k).gt.ZERO)
          savg  = HALF * (splft + sprgt)
          splus = merge(splus, savg, abs(utrans(i+1,j,k)) .gt. eps)

          ulo = u(i-1,j,k,1) + w0_cart_vec(i-1,j,k,1)
          uhi = u(i  ,j,k,1) + w0_cart_vec(i  ,j,k,1)

          smlft = s(i-1,j,k,n) + (HALF - dth*ulo/hx) * slopex(i-1,j,k,1)
          smrgt = s(i  ,j,k,n) - (HALF + dth*uhi/hx) * slopex(i  ,j,k,1)

          smrgt = merge(s(is-1,j,k,n),smrgt,i.eq.is .and. phys_bc(1,1) .eq. INLET)
          smlft = merge(s(is-1,j,k,n),smlft,i.eq.is .and. phys_bc(1,1) .eq. INLET)

          if (i .eq. is .and. (phys_bc(1,1).eq.SLIP_WALL.or.phys_bc(1,1).eq.NO_SLIP_WALL)) then
            if (is_vel .and. n .eq. 1) then
              smlft = ZERO
              smrgt = ZERO
            elseif (is_vel .and. n .ne. 1) then
              smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
              smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
            else
              smlft = smrgt
            endif
          endif

          sminus = merge(smlft,smrgt,utrans(i,j,k).gt.ZERO)
          savg   = HALF * (smlft + smrgt)
          sminus = merge(sminus, savg, abs(utrans(i,j,k)) .gt. eps)

          st = force(i,j,k,n) - &
               HALF * (utrans(i,j,k)+utrans(i+1,j,k))*(splus - sminus) / hx

          ! Do transverse in k direction

          wlo = u(i,j,k  ,3) + w0_cart_vec(i,j,k  ,3)
          whi = u(i,j,k+1,3) + w0_cart_vec(i,j,k+1,3)

          splft = s(i,j,k  ,n) + (HALF - dth*wlo/hz) * slopex(i  ,j,k,1)
          sprgt = s(i,j,k+1,n) - (HALF + dth*whi/hz) * slopex(i+1,j,k,1)

          sprgt = merge(s(i,j,ke+1,n),sprgt,k.eq.ke .and. phys_bc(3,2) .eq. INLET)
          splft = merge(s(i,j,ke+1,n),splft,k.eq.ke .and. phys_bc(3,2) .eq. INLET)

          if (k .eq. ke .and. (phys_bc(3,2).eq.SLIP_WALL.or.phys_bc(3,2).eq.NO_SLIP_WALL)) then
            if (is_vel .and. n .eq. 3) then
              splft = ZERO
              sprgt = ZERO
            elseif (is_vel .and. n .ne. 3) then
              sprgt = merge(ZERO,splft,phys_bc(3,2).eq.NO_SLIP_WALL)
              splft = merge(ZERO,splft,phys_bc(3,2).eq.NO_SLIP_WALL)
            else
              sprgt = splft
            endif
          endif

          splus = merge(splft,sprgt,wtrans(i,j,k+1).gt.ZERO)
          savg  = HALF * (splft + sprgt)
          splus = merge(splus, savg, abs(wtrans(i,j,k+1)) .gt. eps)

          wlo = u(i,j,k-1,3) + w0_cart_vec(i,j,k-1,3)
          whi = u(i,j,k  ,3) + w0_cart_vec(i,j,k  ,3)

          smrgt = s(i,j,k  ,n) - (HALF + dth*whi/hz) * slopez(i,j,k  ,1)
          smlft = s(i,j,k-1,n) + (HALF - dth*wlo/hz) * slopez(i,j,k-1,1)

          smrgt = merge(s(i,j,ks-1,n),smrgt,k.eq.ks .and. phys_bc(3,1) .eq. INLET)
          smlft = merge(s(i,j,ks-1,n),smlft,k.eq.ks .and. phys_bc(3,1) .eq. INLET)

          if (k .eq. ks .and. (phys_bc(3,1).eq.SLIP_WALL.or.phys_bc(3,1).eq.NO_SLIP_WALL)) then
            if (is_vel .and. n .eq. 3) then
              smlft = ZERO
              smrgt = ZERO
            elseif (is_vel .and. n .ne. 3) then
              smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
              smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
            else
              smlft = smrgt
            endif
          endif

          sminus = merge(smlft,smrgt,wtrans(i,j,k).gt.ZERO)
          savg   = HALF * (smlft + smrgt)
          sminus = merge(sminus, savg, abs(wtrans(i,j,k)) .gt. eps)

          st = st - HALF * (wtrans(i,j,k)+wtrans(i,j,k+1))*(splus - sminus) / hz

          ! NOTE NOTE : THIS IS WRONG FOR SPHERICAL !!
          if (spherical .eq. 0 .and. is_vel .and. n.eq.3) then
            st = st - HALF * (wtrans(i,j,k)+wtrans(i,j,k+1))*(w0(k+1)-w0(k))/hz
          end if

          vbardth = dth/hy * ( u(i,j,k,2) + w0_cart_vec(i,j,k,2))

          s_b(j+1)= s(i,j,k,n) + (HALF-vbardth)*slopey(i,j,k,1) + dth*st
          s_t(j  )= s(i,j,k,n) - (HALF+vbardth)*slopey(i,j,k,1) + dth*st
        enddo

        if (velpred .eq. 1) then
          do j = js, je+1 
            savg = HALF*(s_b(j) + s_t(j))
            test = ( (s_b(j) .le. ZERO  .and. &
                      s_t(j) .ge. ZERO)  .or. &
                   (abs(s_b(j) + s_t(j)) .lt. eps) )
            sedgey(i,j,k,n)=merge(s_b(j),s_t(j),savg.gt.ZERO)
            sedgey(i,j,k,n)=merge(savg,sedgey(i,j,k,n),test)
          enddo
        else
          do j = js, je+1 
            sedgey(i,j,k,n)=merge(s_b(j),s_t(j),vadv(i,j,k).gt.ZERO)
            savg = HALF*(s_b(j) + s_t(j))
            sedgey(i,j,k,n)=merge(savg,sedgey(i,j,k,n),abs(vadv(i,j,k)) .lt. eps)
          enddo
        endif

        if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_WALL) then
          if (is_vel .and. n .eq. 2) then
            sedgey(i,js,k,n) = ZERO
          elseif (is_vel .and. n .ne. 2) then
            sedgey(i,js,k,n) = merge(ZERO,s_t(js),phys_bc(2,1).eq.NO_SLIP_WALL)
          else 
            sedgey(i,js,k,n) = s_t(js)
          endif
        elseif (phys_bc(2,1) .eq. INLET) then
          sedgey(i,js,k,n) = s(i,js-1,k,n)
        elseif (phys_bc(2,1) .eq. OUTLET) then
          sedgey(i,js,k,n) = s_t(js)
        endif

        if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_WALL) then
          if (is_vel .and. n .eq. 2) then
            sedgey(i,je+1,k,n) = ZERO
          elseif (is_vel .and. n .ne. 2) then
            sedgey(i,je+1,k,n) = merge(ZERO,s_b(je+1),phys_bc(2,2).eq.NO_SLIP_WALL)
          else 
            sedgey(i,je+1,k,n) = s_b(je+1)
          endif
        elseif (phys_bc(2,2) .eq. INLET) then
          sedgey(i,je+1,k,n) = s(i,je+1,k,n)
        elseif (phys_bc(2,2) .eq. OUTLET) then
          sedgey(i,je+1,k,n) = s_b(je+1)
        endif

        if (velpred .eq. 1) then
          do j = js, je+1 
            vadv(i,j,k) = sedgey(i,j,k,2)
          enddo
        endif

        enddo
        enddo
       endif
!
!     Loop for fluxes on z-edges.
!
       if (velpred .eq. 0 .or. n .eq. 3) then
       do j = js, je 
       do i = is, ie 
        do k = ks-1,ke+1

          ! Do transverse in i direction

          ulo = u(i  ,j,k,1) + w0_cart_vec(i  ,j,k,1)
          uhi = u(i+1,j,k,1) + w0_cart_vec(i+1,j,k,1)

          splft = s(i  ,j,k,n) + (HALF - dth*ulo/hx) * slopex(i  ,j,k,1)
          sprgt = s(i+1,j,k,n) - (HALF + dth*uhi/hx) * slopex(i+1,j,k,1)

          sprgt = merge(s(ie+1,j,k,n),sprgt,i.eq.ie .and. phys_bc(1,2) .eq. INLET)
          splft = merge(s(ie+1,j,k,n),splft,i.eq.ie .and. phys_bc(1,2) .eq. INLET)

          if (i .eq. ie .and. (phys_bc(1,2).eq.SLIP_WALL.or.phys_bc(1,2).eq.NO_SLIP_WALL)) then
            if (is_vel .and. n .eq. 1) then
              splft = ZERO
              sprgt = ZERO
            elseif (is_vel .and. n .ne. 1) then
              sprgt = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
              splft = merge(ZERO,splft,phys_bc(1,2).eq.NO_SLIP_WALL)
            else
              sprgt = splft
            endif
          endif

          splus = merge(splft,sprgt,utrans(i+1,j,k).gt.ZERO)
          savg  = HALF * (splft + sprgt)
          splus = merge(splus, savg, abs(utrans(i+1,j,k)) .gt. eps)

          ulo = u(i-1,j,k,1) + w0_cart_vec(i-1,j,k,1)
          uhi = u(i  ,j,k,1) + w0_cart_vec(i  ,j,k,1)

          smlft = s(i-1,j,k,n) + (HALF - dth*ulo/hx) * slopex(i-1,j,k,1)
          smrgt = s(i  ,j,k,n) - (HALF + dth*uhi/hx) * slopex(i  ,j,k,1)

          smrgt = merge(s(is-1,j,k,n),smrgt,i.eq.is .and. phys_bc(1,1) .eq. INLET)
          smlft = merge(s(is-1,j,k,n),smlft,i.eq.is .and. phys_bc(1,1) .eq. INLET)

          if (i .eq. is .and. (phys_bc(1,1).eq.SLIP_WALL.or.phys_bc(1,1).eq.NO_SLIP_WALL)) then
            if (is_vel .and. n .eq. 1) then
              smlft = ZERO
              smrgt = ZERO
            elseif (is_vel .and. n .ne. 1) then
              smlft = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
              smrgt = merge(ZERO,smrgt,phys_bc(1,1).eq.NO_SLIP_WALL)
            else
              smlft = smrgt
            endif
          endif

          sminus = merge(smlft,smrgt,utrans(i,j,k).gt.ZERO)
          savg   = HALF * (smlft + smrgt)
          sminus = merge(sminus, savg, abs(utrans(i,j,k)) .gt. eps)

          st = force(i,j,k,n) - &
               HALF * (utrans(i,j,k)+utrans(i+1,j,k))*(splus - sminus) / hx

          ! Do transverse in j direction

          vlo = u(i,j  ,k,2) + w0_cart_vec(i,j  ,k,2)
          vhi = u(i,j+1,k,2) + w0_cart_vec(i,j+1,k,2)

          spbot = s(i,j  ,k,n) + (HALF - dth*vlo/hy) * slopey(i,j  ,k,1)
          sptop = s(i,j+1,k,n) - (HALF + dth*vhi/hy) * slopey(i,j+1,k,1)

          sptop = merge(s(i,je+1,k,n),sptop,j.eq.je .and. phys_bc(2,2) .eq. INLET)
          spbot = merge(s(i,je+1,k,n),spbot,j.eq.je .and. phys_bc(2,2) .eq. INLET)

          if (j .eq. je .and. (phys_bc(2,2).eq.SLIP_WALL.or.phys_bc(2,2).eq.NO_SLIP_WALL)) then
            if (is_vel .and. n.eq.2) then
              sptop = ZERO
              spbot = ZERO
            elseif (is_vel .and. n.ne.2) then
              sptop = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
              spbot = merge(ZERO,spbot,phys_bc(2,2).eq.NO_SLIP_WALL)
            else
              sptop = spbot
            endif
          endif

          splus = merge(spbot,sptop,vtrans(i,j+1,k).gt.ZERO)
          savg  = HALF * (spbot + sptop)
          splus = merge(splus, savg, abs(vtrans(i,j+1,k)) .gt. eps)

          vlo = u(i,j-1,k,2) + w0_cart_vec(i,j-1,k,2)
          vhi = u(i,j  ,k,2) + w0_cart_vec(i,j  ,k,2)

          smbot = s(i,j-1,k,n) + (HALF - dth*vlo/hy) * slopey(i,j-1,k,1)
          smtop = s(i,j  ,k,n) - (HALF + dth*vhi/hy) * slopey(i,j  ,k,1)

          smtop = merge(s(i,js-1,k,n),smtop,j.eq.js .and. phys_bc(2,1) .eq. INLET)
          smbot = merge(s(i,js-1,k,n),smbot,j.eq.js .and. phys_bc(2,1) .eq. INLET)

          if (j .eq. js .and. (phys_bc(2,1).eq.SLIP_WALL.or.phys_bc(2,1).eq.NO_SLIP_WALL)) then
            if (is_vel .and. (n .eq. 2)) then
              smtop = ZERO
              smbot = ZERO
            elseif (is_vel .and. (n .ne. 2)) then
              smbot = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
              smtop = merge(ZERO,smtop,phys_bc(2,1).eq.NO_SLIP_WALL)
            else
              smbot = smtop
            endif
          endif

          sminus = merge(smbot,smtop,vtrans(i,j,k).gt.ZERO)
          savg   = HALF * (smbot + smtop)
          sminus = merge(sminus, savg, abs(vtrans(i,j,k)) .gt. eps)

          st = st - HALF * (vtrans(i,j,k)+vtrans(i,j+1,k))*(splus - sminus) / hy

          ! NOTE NOTE : THIS IS WRONG FOR SPHERICAL !!
          if (spherical .eq. 0 .and. is_vel .and. n.eq.3) then
            st = st - HALF * (wtrans(i,j,k)+wtrans(i,j,k+1))*(w0(k+1)-w0(k))/hz
          end if

          wbardth = dth/hz * ( u(i,j,k,3) + w0_cart_vec(i,j,k,3))

          s_d(k+1)= s(i,j,k,n) + (HALF-wbardth)*slopez(i,j,k,1) + dth*st
          s_u(k  )= s(i,j,k,n) - (HALF+wbardth)*slopez(i,j,k,1) + dth*st
        enddo

        if (velpred .eq. 1) then
          do k = ks, ke+1 
            savg = HALF*(s_d(k) + s_u(k))
            test = ( (s_d(k) .le. ZERO  .and. &
                      s_u(k) .ge. ZERO)  .or. &
                   (abs(s_d(k) + s_u(k)) .lt. eps) )
            sedgez(i,j,k,n)=merge(s_d(k),s_u(k),savg.gt.ZERO)
            sedgez(i,j,k,n)=merge(savg,sedgez(i,j,k,n),test)
          enddo
        else
          do k = ks, ke+1 
            sedgez(i,j,k,n)=merge(s_d(k),s_u(k),wadv(i,j,k).gt.ZERO)
            savg = HALF*(s_d(k) + s_u(k))
            sedgez(i,j,k,n)=merge(savg,sedgez(i,j,k,n),abs(wadv(i,j,k)) .lt. eps)
          enddo
          ! OUTFLOW HACK 
          if (spherical .eq. 0 .and. phys_bc(3,2) .eq. OUTLET) then
            if (wadv(i,j,ke+1).lt.ZERO) sedgez(i,j,ke+1,n) = s_d(ke+1)
          end if
        endif

        if (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_WALL) then
          if (is_vel .and. n .eq. 2) then
            sedgez(i,j,ks,n) = ZERO
          elseif (is_vel .and. n .ne. 2) then
            sedgez(i,j,ks,n) = merge(ZERO,s_u(ks),phys_bc(3,1).eq.NO_SLIP_WALL)
          else 
            sedgez(i,j,ks,n) = s_u(ks)
          endif
        elseif (phys_bc(3,1) .eq. INLET) then
          sedgez(i,j,ks,n) = s(i,j,ks-1,n)
        elseif (phys_bc(3,1) .eq. OUTLET) then
          sedgez(i,j,ks,n) = s_u(ks)
        endif

        if (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_WALL) then
          if (is_vel .and. n .eq. 2) then
            sedgez(i,j,ke+1,n) = ZERO
          elseif (is_vel .and. n .ne. 2) then
            sedgez(i,j,ke+1,n) = merge(ZERO,s_d(ke+1),phys_bc(3,2).eq.NO_SLIP_WALL)
          else 
            sedgez(i,j,ke+1,n) = s_d(ke+1)
          endif
        elseif (phys_bc(3,2) .eq. INLET) then
          sedgez(i,j,ke+1,n) = s(i,j,ke+1,n)
        elseif (phys_bc(3,2) .eq. OUTLET) then
          sedgez(i,j,ke+1,n) = s_d(ke+1)
        endif

        if (velpred .eq. 1) then
          do k = ks, ke+1 
            wadv(i,j,k) = sedgez(i,j,k,3)
          enddo
        endif

        enddo
        enddo
       endif

      if (.not. is_vel .and. advect_in_pert_form) then
         if (spherical .eq. 1) then
           do k = ks,ke
             do j = js,je
             do i = is,ie
               s(i,j,k,n) = s(i,j,k,n) + base_cart(i,j,k)
             end do
             end do
             do g = 1,ng
               do j = js,je
                 s(is-g,j,k,n) = s(is,j,k,n)
                 s(ie+g,j,k,n) = s(ie,j,k,n)
               end do
               do i = is-1,ie+1
                 s(i,js-g,k,n) = s(i,j,k,n)
                 s(i,je+g,k,n) = s(i,je,k,n)
               end do
             end do
           end do
         else
           do k = ks,ke
             do j = js-ng,je+ng
             do i = is-ng,ie+ng
               s(i,j,k,n) = s(i,j,k,n) + base(k)
             end do
             end do
           end do
         end if

        do k = ks,ke
          do g = 1,ng
           do j = js-ng,je+ng
           do i = is-ng,ie+ng
            s(i,j,ks-g,n) = s(i,j,ks,n)
            s(i,j,ke+g,n) = s(i,j,ke,n)
          end do
          end do
         end do
        end do
      end if

      deallocate(s_l)
      deallocate(s_r)
      deallocate(s_b)
      deallocate(s_t)
      deallocate(s_d)
      deallocate(s_u)

      deallocate(slopex)
      deallocate(slopey)
      deallocate(slopez)

      if (.not. is_vel .and. advect_in_pert_form .and. (spherical .eq. 1) ) &
        deallocate(base_cart)

      end subroutine mkflux_3d

      subroutine mkflux_1d(s,sedgex,uadv,force,lo,dx,dt)

      integer        , intent(in   ) :: lo
      real(kind=dp_t), intent(in   ) ::      s(lo:)
      real(kind=dp_t), intent(inout) :: sedgex(lo:)
      real(kind=dp_t), intent(in   ) ::   uadv(lo:)
      real(kind=dp_t), intent(in   ) ::  force(lo:)
      real(kind=dp_t), intent(in   ) :: dx,dt

!     Local variables
      real(kind=dp_t), allocatable::  slopex(:)
      real(kind=dp_t), allocatable::  s_l(:),s_r(:)
      real(kind=dp_t), allocatable:: dxscr(:,:)
      real(kind=dp_t) :: dmin,dpls,ds
      real(kind=dp_t) :: ubardth, dth, savg
      real(kind=dp_t) :: abs_eps, eps, umax, u
      real(kind=dp_t) :: fourthirds,sixth

      integer :: i,is,ie
      integer :: hi,cen,lim,flag,fromm
      integer :: slope_order = 4
      logical :: test

      parameter( cen = 1 )
      parameter( lim = 2 )
      parameter( flag = 3 )
      parameter( fromm = 4 )

      hi = lo + size(s,dim=1) - 1

      allocate(s_l(lo-1:hi+2),s_r(lo-1:hi+2))
      allocate(slopex(lo:hi))

      allocate(dxscr(lo:hi,4))

      abs_eps = 1.0e-8

      dth = HALF*dt

      is = lo
      ie = lo + size(s,dim=1) - 1

      umax = ZERO
      do i = is,ie+1
        umax = max(umax,abs(uadv(i)))
      end do

      eps = abs_eps * umax

      ! Compute fourth-order slopes
      do i = is+1,ie-1
        dxscr(i,cen) = half*(s(i+1)-s(i-1))
        dmin = two*(s(i  )-s(i-1))
        dpls = two*(s(i+1)-s(i  ))
        dxscr(i,lim)= min(abs(dmin),abs(dpls))
        dxscr(i,lim) = merge(dxscr(i,lim),zero,dpls*dmin.gt.ZERO)
        dxscr(i,flag) = sign(one,dxscr(i,cen))
        dxscr(i,fromm)= dxscr(i,flag)*min(dxscr(i,lim), &
                        abs(dxscr(i,cen)))
      enddo

      dxscr(is,fromm) = ZERO
      dxscr(ie,fromm) = ZERO

      fourthirds = 4.0_dp_t / 3.0_dp_t
      sixth      = 1.0_dp_t / 6.0_dp_t

      do i = is+1,ie-1
         ds = fourthirds * dxscr(i,cen) - &
              sixth * (dxscr(i+1,fromm) + dxscr(i-1,fromm))
         slopex(i) = dxscr(i,flag)*min(abs(ds),dxscr(i,lim))
      enddo

      slopex(is) = ZERO
      slopex(ie) = ZERO

      ! Use fourth-order slopes to compute edge values
      do i = is,ie

         u = HALF * (uadv(i) + uadv(i+1))
         ubardth = dth*u/dx

         s_l(i+1)= s(i) + (HALF-ubardth)*slopex(i) + dth * force(i)
         s_r(i  )= s(i) - (HALF+ubardth)*slopex(i) + dth * force(i)

      enddo

      sedgex(is  ) = s_r(is  )
      sedgex(ie+1) = s_l(ie+1)

      do i = is+1, ie 
        sedgex(i)=merge(s_l(i),s_r(i),uadv(i).gt.ZERO)
        savg = HALF*(s_r(i) + s_l(i))
        sedgex(i)=merge(savg,sedgex(i),abs(uadv(i)) .lt. eps)
      enddo

      deallocate(s_l)
      deallocate(s_r)
      deallocate(slopex)
      deallocate(dxscr)

      end subroutine mkflux_1d

end module mkflux_module
