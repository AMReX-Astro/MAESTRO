
program oned 

      implicit none

      integer, parameter :: dp_t = selected_real_kind(15,307)

      integer, parameter :: nz = 768
      real(kind=dp_t) :: z(nz),rho0(nz),rhoh0(nz),p0(nz),temp0(nz)
      real(kind=dp_t) :: temp_new(nz),den_new(nz),sx(nz)
      real(kind=dp_t) :: rho0_edge(nz+1), rhoh0_edge(nz+1), p0_edge(nz+1)
      real(kind=dp_t) :: rho0_new,rhoh0_new
      real(kind=dp_t) :: vel(nz+1),maxvel
      real(kind=dp_t) :: base_state(7,nz)
      real(kind=dp_t) :: dx(2),dt

      integer, parameter :: IONMAX = 2
      integer, parameter :: NP = 1
      real(kind=dp_t) :: xmass(IONMAX), aion(IONMAX), zion(IONMAX)
      real(kind=dp_t) :: temp_row(NP)
      real(kind=dp_t) :: den_row(NP)
      real(kind=dp_t) :: abar_row(NP)
      real(kind=dp_t) :: zbar_row(NP)
      real(kind=dp_t) :: e_row(NP)
      real(kind=dp_t) :: p_row(NP)
      real(kind=dp_t) :: h_row(NP)
      real(kind=dp_t) :: cv_row(NP)
      real(kind=dp_t) :: cp_row(NP)
      real(kind=dp_t) :: xne_row(NP)
      real(kind=dp_t) :: eta_row(NP)
      real(kind=dp_t) :: pele_row(NP)
      real(kind=dp_t) :: dpdt_row(NP)
      real(kind=dp_t) :: dpdr_row(NP)
      real(kind=dp_t) :: dedr_row(NP)
      real(kind=dp_t) :: dedt_row(NP)
      real(kind=dp_t) :: gam1_row(NP)
      real(kind=dp_t) ::   cs_row(NP)
      real(kind=dp_t) ::    s_row(NP)

      real(kind=dp_t) :: y,y0,Hbar,coeff,disp
      real(kind=dp_t) :: rho0_lo,rho0_hi,p0_lo,p0_hi
      real(kind=dp_t) :: log_rho_new, log_rho_orig
      integer :: input_flag,npts,nspecies
      integer :: it,nt,i,j
      logical :: do_diag
      
      do_diag = .false.

      y0 = 0.40e8

!     dt = 0.1_dp_t
      dt = 0.01_dp_t

!     nt = 5.0000001_dp_t / dt
      nt = 2.0000001_dp_t / dt
      print *,'NT ',nt

      !..set the mass fractions, z's and a's of the composition

      !..carbon 12
      aion(1)  = 12.0d0
      zion(1)  = 6.0d0

      !..oxygen 16
      aion(2)  = 16.0d0
      zion(2)  = 8.0d0

      !..call the eos
      npts     = 1
      nspecies = 2

      xmass(1) = 0.2999999997
      xmass(2) = 0.6999999993

      open(99,file="model.hse")
      do i = 1,nz
        read(99,*) base_state(1,i),base_state(2,i),base_state(3,i),base_state(4,i), &
                   base_state(5,i),base_state(6,i), base_state(7,i)
            z(i) = base_state(1,i)
         rho0(i) = base_state(2,i)
        temp0(i) = base_state(3,i)
           p0(i) = base_state(4,i)
           write(11,999) z(i), rho0(i)
           write(12,999) z(i), temp0(i)
           write(13,999) z(i), p0(i)
      end do
      close(99)
 
      dx(2) = base_state(1,2) - base_state(1,1)

      call helmeos_init

!     do j = 1,nz

           ! (rho, p) --> T,h, etc
!          input_flag = 4
!           den_row(1) =  rho0(j)
!             p_row(1) =    p0(j)
!          temp_row(1) = temp0(j)

!          call eos(input_flag, den_row, temp_row, npts, nspecies, &
!               xmass, aion, zion, &
!               p_row, h_row, e_row, &
!               cv_row, cp_row, xne_row, eta_row, &
!               pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
!               s_row, do_diag)
!          write(14,*) z(j),temp0(j),temp_row(1)

!     end do
!     stop

      do it = 1, nt

        print *,'DOING LOOP ',it, ' TO TIME ', it*dt
      ! coeff = p_T / (rho * c_p * p_rho)
      vel(1) = 0.0_dp_t
      maxvel = 0.0_dp_t
      do j = 2,nz+1

        ! Compute the horizontal average of the heating term
        Hbar = 1.e17 * exp(-((z(j-1)-y0)**2) / 1.e14 )

           ! (rho, T) --> p,h, etc
           input_flag = 1
            den_row(1) = rho0(j-1)
           temp_row(1) = temp0(j-1)

           call eos(input_flag, den_row, temp_row, npts, nspecies, &
                xmass, aion, zion, &
                p_row, h_row, e_row, &
                cv_row, cp_row, xne_row, eta_row, &
                pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
                s_row, do_diag)
           rhoh0(j-1) = den_row(1)*h_row(1)

           ! Compute the coefficient of heating in the divu expression
           coeff = dpdt_row(1) / (rho0(j-1) * cp_row(1) * dpdr_row(1))

           ! Compute the displacement velocity
           vel(j) = vel(j-1) + coeff * Hbar * dx(2)
           if (j.le.510) then
             maxvel = max(maxvel,abs(vel(j)))
!            print *,'W DOT GRAD B0 ', j,coeff * Hbar * dx(2), &
!              0.25_dp_t * (vel(j)+vel(j-1)) * (rho0(j+1) - rho0(j-1)) / rho0(j)
           end if

           ! Compute the p0 corresponding to T and rho
           if (it.eq.1) write(14,*) z(j-1),p0(j-1),p_row(1)
           p0(j-1) = p_row(1)

      end do
      do j = 2,nz
         write(55,*) 0.5_dp_t*(z(j)+z(j-1)),log(vel(j))/log(10.)
      end do
!     stop

      print *,'MAXVEL ',it,maxvel

      disp = vel(nz+1) * dt / dx(2)

!     print *,'MAX CFL FRAC OF DISPL ',disp 

      ! Compute the new base state.

      ! Update p0
      call mkflux_1d(   p0,   p0_edge,vel,1,dx(2),dt)
      do j = 1,nz
         p0(j) = p0(j) - dt / dx(2) * 0.5_dp_t * (vel(j)+vel(j+1)) *  (p0_edge(j+1) - p0_edge(j))
      end do

      ! Update rho*h
      call mkflux_1d(rhoh0,rhoh0_edge,vel,1,dx(2),dt)
      do j = 1,nz
         rhoh0(j)= rhoh0(j) - dt / dx(2) * ( rhoh0_edge(j+1) * vel(j+1) -  rhoh0_edge(j) * vel(j))
      end do

      call mkflux_1d( rho0, rho0_edge,vel,1,dx(2),dt)
      do j = 1,nz

         rho0_new = rho0(j) - dt / dx(2) * (rho0_edge(j+1) * vel(j+1) - rho0_edge(j) * vel(j))
         rho0(j) = max(rho0_new, rho0(nz))
  
         ! Initial guess
           den_row(1) = rho0(j)
             p_row(1) = p0(j)
          temp_row(1) = temp0(j)

         ! (rho, p) --> T
         input_flag = 4

         call eos(input_flag, den_row, temp_row, npts, nspecies, &
              xmass, aion, zion, &
              p_row, h_row, e_row, &
              cv_row, cp_row, xne_row, eta_row, &
              pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
              s_row, do_diag)

         temp0(j) = temp_row(1)

      end do

      end do

      do j = 1,nz
         write(41,999) z(j), rho0(j)
         write(42,999) z(j),temp0(j)
         write(43,999) z(j),   p0(j)
      end do

 999  format(e18.10,x,e18.10)
1000  format(e18.10,x,e18.10,x,e18.10)
1001  format(e18.10,x,e18.10,x,e18.10,x,e18.10)

   contains

      subroutine mkflux_1d(s,sedgex,uadv,lo,dx,dt)

      use bl_constants_module

      integer        , intent(in   ) :: lo
      real(kind=dp_t), intent(inout) ::      s(lo:)
      real(kind=dp_t), intent(inout) :: sedgex(lo:)
      real(kind=dp_t), intent(inout) ::   uadv(lo:)
      real(kind=dp_t), intent(in   ) :: dx,dt

!     Local variables
      real(kind=dp_t), allocatable::  slopex(:)
      real(kind=dp_t), allocatable::  s_l(:),s_r(:)
      real(kind=dp_t), allocatable:: dxscr(:,:)
      real(kind=dp_t) :: dmin,dpls,ds
      real(kind=dp_t) :: ubardth, dth, savg
      real(kind=dp_t) :: abs_eps, eps, umax, u
      real(kind=dp_t) :: fourthirds

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

      do i = is,ie
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

      slopex(is) = ZERO
      slopex(ie) = ZERO

      fourthirds = 4.0_dp_t / 3.0_dp_t
!     sixth      = 1.0_dp_t / 6.0_dp_t

      do i = is+1,ie-1
         ds = fourthirds * dxscr(i,cen) - &
              sixth * (dxscr(i+1,fromm) + dxscr(i-1,fromm))
         slopex(i) = dxscr(i,flag)*min(abs(ds),dxscr(i,lim))
      enddo

      ! Use fourth-order slopes to compute edge values
      do i = is,ie

         u = HALF * (uadv(i) + uadv(i+1))
         ubardth = dth*u/dx

         s_l(i+1)= s(i) + (HALF-ubardth)*slopex(i)
         s_r(i  )= s(i) - (HALF+ubardth)*slopex(i)

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

   end program oned
