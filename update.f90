module update_module

  use bl_types
  use multifab_module
  use variables
  use network

  implicit none

  real (kind = dp_t), private, parameter :: ZERO  = 0.0_dp_t
  real (kind = dp_t), private, parameter :: HALF  = 0.5_dp_t

  contains

   subroutine update_scal_2d (nstart,nstop,sold,snew,umac,vmac,w0,sedgex,sedgey,force, &
                              base_old,base_new,lo,hi,ng,dx,dt,verbose)

     ! update each scalar in time.  Here, it is assumed that the edge
     ! states (sedgex and sedgey) are for the perturbational quantities.

      implicit none

      integer              , intent(in) :: nstart, nstop, lo(:), hi(:), ng, verbose
      real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,:)
      real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,:)
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:)
      real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:)
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,:)
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,:)
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,:)
      real (kind = dp_t), intent(in   ) ::   base_old(lo(2)   :,:)
      real (kind = dp_t), intent(in   ) ::   base_new(lo(2)   :,:)
      real (kind = dp_t), intent(in   ) :: w0(lo(2):)
      real (kind = dp_t), intent(in   ) :: dt,dx(:)

      integer :: i, j, n
      real (kind = dp_t) :: divsu,divbaseu,smax_rho,smin_rho
      real (kind = dp_t), allocatable :: smin(:)
      real (kind = dp_t), allocatable :: smax(:)
      real (kind = dp_t), allocatable :: base_edge(:)

      allocate(base_edge(lo(2):hi(2)+1))
      allocate(smax(nstart:nstop))
      allocate(smin(nstart:nstop))

      if (0.eq.1 .and. verbose .ge. 1) then
        smax(:) = -1.d20
        smin(:) =  1.d20
        do n = nstart, nstop
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            if (n .ge. spec_comp .and. n .le. spec_comp+nspec-1) then
              smax(n) = max(smax(n),sold(i,j,n)/sold(i,j,rho_comp))
              smin(n) = min(smin(n),sold(i,j,n)/sold(i,j,rho_comp))
            else
              smax(n) = max(smax(n),sold(i,j,n))
              smin(n) = min(smin(n),sold(i,j,n))
            endif

          enddo
          enddo

          if (n.eq.rhoh_comp) write(6,1001) smin(n),smax(n)
          if (n.gt.rhoh_comp .and. n.lt.trac_comp) write(6,1002) spec_names(n-rhoh_comp),smin(n),smax(n)
          if (n.ge.trac_comp) write(6,1003) smin(n),smax(n)
        enddo

        if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
          n = rho_comp
          smax_rho = -1.d20
          smin_rho =  1.d20
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
              smax_rho = max(smax_rho,sold(i,j,n))
              smin_rho = min(smin_rho,sold(i,j,n))
          end do
          end do
          write(6,1000) smin_rho,smax_rho
        end if

      end if

1000  format('OLD MIN/MAX : density           ',e15.10,2x,e15.10)
1001  format('OLD MIN/MAX : rho * H           ',e15.10,2x,e15.10)
1002  format('OLD MIN/MAX : ',a16,2x,e15.10,2x,e15.10)
1003  format('OLD MIN/MAX :           tracer',2x,e15.10,2x,e15.10)


      do n = nstart, nstop

        base_edge(lo(2)  ) = base_old(lo(2),n)
        base_edge(hi(2)+1) = base_old(hi(2),n)
      
        base_edge(lo(2)+1) = HALF*(base_old(lo(2),n)+base_old(lo(2)+1,n))
        base_edge(hi(2)  ) = HALF*(base_old(hi(2),n)+base_old(hi(2)-1,n))

        do j = lo(2)+2,hi(2)-1
           base_edge(j) = 7.d0/12.d0 * (base_old(j  ,n) + base_old(j-1,n)) &
                         -1.d0/12.d0 * (base_old(j+1,n) + base_old(j-2,n))
        end do

        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
  
          divsu = (umac(i+1,j) * sedgex(i+1,j,n) &
                  -umac(i  ,j) * sedgex(i  ,j,n) ) / dx(1) + &
                 ((vmac(i,j+1)+w0(j+1)) * sedgey(i,j+1,n) &
                 -(vmac(i,j  )+w0(j  )) * sedgey(i,j  ,n) ) / dx(2)

          divbaseu = (umac(i+1,j) - umac(i,j) ) * base_old(j,n) / dx(1) &
                    +(vmac(i,j+1) * base_edge(j+1) - vmac(i,j) * base_edge(j) ) / dx(2)

          snew(i,j,n) = sold(i,j,n) + (base_new(j,n) - base_old(j,n)) &
                      - dt * (divsu + divbaseu) + dt * force(i,j,n)
  
        enddo
        enddo
      enddo

      if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
        snew(:,:,rho_comp) = ZERO
        do n = nstart, nstop
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           snew(i,j,rho_comp) = snew(i,j,rho_comp) + snew(i,j,n)
        enddo
        enddo
        enddo
      end if
  
      if (verbose .ge. 1) then

        smax(:) = -1.d20
        smin(:) =  1.d20

        do n = nstart, nstop
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
  
            if (n .ge. spec_comp .and. n .le. spec_comp+nspec-1) then
              smax(n) = max(smax(n),snew(i,j,n)/snew(i,j,rho_comp))
              smin(n) = min(smin(n),snew(i,j,n)/snew(i,j,rho_comp))
            else
              smax(n) = max(smax(n),snew(i,j,n))
              smin(n) = min(smin(n),snew(i,j,n))
            endif
    
          enddo
          enddo
  
          if (n.eq.rhoh_comp) write(6,2001) smin(n),smax(n)
          if (n.gt.rhoh_comp .and. n.lt.trac_comp) write(6,2002) spec_names(n-rhoh_comp),smin(n),smax(n)
          if (n.ge.trac_comp) write(6,2003) smin(n),smax(n)
          if (n.eq.rhoh_comp) write(6,2004) 
        enddo

        if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
          n = rho_comp
          smax_rho = -1.d20
          smin_rho =  1.d20
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
              smax_rho = max(smax_rho,snew(i,j,n))
              smin_rho = min(smin_rho,snew(i,j,n))
          end do
          end do
          write(6,2000) smin_rho,smax_rho
        end if

      end if

2000  format('... new min/max : density           ',e15.10,2x,e15.10)
2001  format('... new min/max : rho * H           ',e15.10,2x,e15.10)
2002  format('... new min/max : ',a16,2x,e15.10,2x,e15.10)
2003  format('... new min/max :           tracer',2x,e15.10,2x,e15.10)
2004  format(' ')

      deallocate(smin,smax)
      deallocate(base_edge)

   end subroutine update_scal_2d

   subroutine update_velocity_2d (uold,unew,umac,vmac,sedgex,sedgey,force,w0, &
                                  lo,hi,ng,dx,time,dt,verbose)


     ! update the velocity in time (to get the provisional velocity that
     ! does not yet satisfy the divergence constraint).  This is the first
     ! part of step 5 in ABRZ2.

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng, verbose
      real (kind = dp_t), intent(in   ) ::    uold(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(  out) ::    unew(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,:)  
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,:)  
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent(in   ) ::      w0(          lo(2)   :)  
      real (kind = dp_t), intent(in   ) :: dx(:)
      real (kind = dp_t), intent(in   ) :: time,dt

      integer :: i, j, n
      real (kind = dp_t) ubar,vbar
      real (kind = dp_t) ugradu,ugradv,ugrads
      real (kind = dp_t) :: divsu
      real (kind = dp_t) :: smin,smax,umin,umax,vmin,vmax
      real (kind = dp_t) :: fac

      print *,'<<< updating velocity ',n,' >>> '

      do j = lo(2), hi(2)
      do i = lo(1), hi(1)

           ubar = HALF*(umac(i,j) + umac(i+1,j))
           vbar = HALF*(vmac(i,j) + vmac(i,j+1))

           ugradu = ubar*(sedgex(i+1,j,1) - sedgex(i,j,1))/dx(1) + &
                    vbar*(sedgey(i,j+1,1) - sedgey(i,j,1))/dx(2)

           ugradv = ubar*(sedgex(i+1,j,2) - sedgex(i,j,2))/dx(1) + &
                    vbar*(sedgey(i,j+1,2) - sedgey(i,j,2))/dx(2)

           unew(i,j,1) = uold(i,j,1) - dt * ugradu + dt * force(i,j,1)
           unew(i,j,2) = uold(i,j,2) - dt * ugradv + dt * force(i,j,2)

           ! Add w dot grad w0 term to w.
           unew(i,j,2) = unew(i,j,2) - dt * vbar*(w0(j+1) - w0(j))/dx(2)

           ! Add w0 dot grad u term to u and w.
           vbar = HALF*(w0(j) + w0(j+1))
           unew(i,j,:) = unew(i,j,:) - dt * vbar*(sedgey(i,j+1,:) - sedgey(i,j,:))/dx(2)

      enddo
      enddo

      umax = unew(lo(1),lo(2),1) 
      umin = unew(lo(1),lo(2),1) 
      vmax = unew(lo(1),lo(2),2) 
      vmin = unew(lo(1),lo(2),2) 
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          umax = max(umax,unew(i,j,1))
          umin = min(umin,unew(i,j,1))
          vmax = max(vmax,unew(i,j,2))
          vmin = min(vmin,unew(i,j,2))
        enddo
      enddo
      if (verbose .ge. 1) then
        write(6,1000) umin,umax
        write(6,1001) vmin,vmax
        write(6,1002)
      end if

1000  format('... new min/max : x-velocity       ',e17.10,2x,e17.10)
1001  format('... new min/max : y-velocity       ',e17.10,2x,e17.10)
1002  format(' ')

   end subroutine update_velocity_2d

   subroutine update_scal_3d (nstart,nstop,sold,snew,umac,vmac,wmac,w0,sedgex,sedgey,sedgez,force, &
                              base_old,base_new,lo,hi,ng,dx,dt,verbose)

      implicit none

      integer, intent(in) :: nstart,nstop, lo(:), hi(:), ng, verbose
      real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
      real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
      real (kind = dp_t), intent(in   ) ::    wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real (kind = dp_t), intent(in   ) ::  sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
      real (kind = dp_t), intent(in   ) ::   base_old(lo(3)   :,:)
      real (kind = dp_t), intent(in   ) ::   base_new(lo(3)   :,:)
      real (kind = dp_t), intent(in   ) :: w0(lo(3):)
      real (kind = dp_t), intent(in   ) :: dt,dx(:)

      integer :: i, j, k, n
      real (kind = dp_t) :: divsu,divbaseu,smax_rho,smin_rho
      real (kind = dp_t), allocatable :: smin(:)
      real (kind = dp_t), allocatable :: smax(:)
      real (kind = dp_t), allocatable :: base_edge(:)

      allocate(base_edge(lo(3):hi(3)+1))
      allocate(smax(nstart:nstop))
      allocate(smin(nstart:nstop))

      smax(:) = -1.d20
      smin(:) =  1.d20

      do n = nstart, nstop

        base_edge(lo(3)  ) = base_old(lo(3),n)
        base_edge(hi(3)+1) = base_old(hi(3),n)
        
        base_edge(lo(3)+1) = HALF*(base_old(lo(3),n)+base_old(lo(3)+1,n))
        base_edge(hi(3)  ) = HALF*(base_old(hi(3),n)+base_old(hi(3)-1,n))
  
        do k = lo(3)+2,hi(3)-1
           base_edge(k) = 7.d0/12.d0 * (base_old(k  ,n) + base_old(k-1,n)) &
                         -1.d0/12.d0 * (base_old(k+1,n) + base_old(k-2,n))
        end do

        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
  
          divsu = (umac(i+1,j,k) * sedgex(i+1,j,k,n) &
                  -umac(i  ,j,k) * sedgex(i  ,j,k,n) ) / dx(1) + &
                  (vmac(i,j+1,k) * sedgey(i,j+1,k,n) &
                  -vmac(i,j  ,k) * sedgey(i,j  ,k,n) ) / dx(2) + &
                 ((wmac(i,j,k+1)+w0(k+1)) * sedgez(i,j,k+1,n) &
                 -(wmac(i,j,k  )+w0(k  )) * sedgez(i,j,k  ,n) ) / dx(3)

          divbaseu = (umac(i+1,j,k) - umac(i,j,k) ) * base_old(k,n) / dx(1) &
                    +(vmac(i,j+1,k) - vmac(i,j,k) ) * base_old(k,n) / dx(2) &
                    +(wmac(i,j,k+1) * base_edge(k+1) - wmac(i,j,k) * base_edge(k) ) / dx(3)

          snew(i,j,k,n) = sold(i,j,k,n) + (base_new(k,n) - base_old(k,n)) &
                          - dt * (divsu + divbaseu) + dt * force(i,j,k,n)
  
        enddo
        enddo
        enddo

      end do

      if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
        snew(:,:,:,rho_comp) = ZERO
        do n = nstart, nstop
        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           snew(i,j,k,rho_comp) = snew(i,j,k,rho_comp) + snew(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
      end if

      if (verbose .ge. 1) then
        do n = nstart, nstop
          do k = lo(3), hi(3)
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
  
            if (n .ge. spec_comp .and. n .le. spec_comp+nspec-1) then
              smax(n) = max(smax(n),snew(i,j,k,n)/snew(i,j,k,rho_comp))
              smin(n) = min(smin(n),snew(i,j,k,n)/snew(i,j,k,rho_comp))
            else
              smax(n) = max(smax(n),snew(i,j,k,n))
              smin(n) = min(smin(n),snew(i,j,k,n))
            endif
  
          enddo
          enddo
          enddo

          if (n.eq.rhoh_comp) write(6,1001) smin(n),smax(n)
          if (n.gt.rhoh_comp .and. n.lt.trac_comp) write(6,1002) spec_names(n-rhoh_comp),smin(n),smax(n)
          if (n.ge.trac_comp) write(6,1003) smin(n),smax(n)
          if (n.eq.rhoh_comp) write(6,1004)
        enddo

        if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
          n = rho_comp
          smax_rho = -1.d20
          smin_rho =  1.d20
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
              smax_rho = max(smax_rho,snew(i,j,k,n))
              smin_rho = min(smin_rho,snew(i,j,k,n))
          end do
          end do
          write(6,1000) smin_rho,smax_rho
        end if

      end if

1000  format('... new min/max : density           ',e15.10,2x,e15.10)
1001  format('... new min/max : rho * H           ',e15.10,2x,e15.10)
1002  format('... new min/max : ',a16,2x,e15.10,2x,e15.10)
1003  format('... new min/max :           tracer',2x,e15.10,2x,e15.10)
1004  format(' ')

      deallocate(smin,smax)
      deallocate(base_edge)

   end subroutine update_scal_3d

   subroutine update_velocity_3d (uold,unew,umac,vmac,wmac,sedgex,sedgey,sedgez, &
                                  force,w0,lo,hi,ng,dx,time,dt,verbose)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng, verbose
      real (kind = dp_t), intent(in   ) ::    uold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(  out) ::    unew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:  )
      real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:  )
      real (kind = dp_t), intent(in   ) ::    wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:  )
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real (kind = dp_t), intent(in   ) ::  sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
      real (kind = dp_t), intent(in   ) ::      w0(          lo(3)   :)  
      real (kind = dp_t), intent(in   ) :: dx(:)
      real (kind = dp_t), intent(in   ) :: time,dt

      integer :: i, j, k, n
      real (kind = dp_t) ubar,vbar,wbar
      real (kind = dp_t) ugradu,ugradv,ugradw,ugrads
      real (kind = dp_t) :: divsu
      real (kind = dp_t) :: smin,smax,umin,umax,vmin,vmax,wmin,wmax
      real (kind = dp_t) :: fac

      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)

           ubar = HALF*(umac(i,j,k) + umac(i+1,j,k))
           vbar = HALF*(vmac(i,j,k) + vmac(i,j+1,k))
           wbar = HALF*(wmac(i,j,k) + wmac(i,j,k+1))

           ugradu = ubar*(sedgex(i+1,j,k,1) - sedgex(i,j,k,1))/dx(1) + &
                    vbar*(sedgey(i,j+1,k,1) - sedgey(i,j,k,1))/dx(2) + &
                    wbar*(sedgez(i,j,k+1,1) - sedgez(i,j,k,1))/dx(3)

           ugradv = ubar*(sedgex(i+1,j,k,2) - sedgex(i,j,k,2))/dx(1) + &
                    vbar*(sedgey(i,j+1,k,2) - sedgey(i,j,k,2))/dx(2) + &
                    wbar*(sedgez(i,j,k+1,2) - sedgez(i,j,k,2))/dx(3)

           ugradw = ubar*(sedgex(i+1,j,k,3) - sedgex(i,j,k,3))/dx(1) + &
                    vbar*(sedgey(i,j+1,k,3) - sedgey(i,j,k,3))/dx(2) + &
                    wbar*(sedgez(i,j,k+1,3) - sedgez(i,j,k,3))/dx(3)

           unew(i,j,k,1) = uold(i,j,k,1) - dt * ugradu + dt * force(i,j,k,1)
           unew(i,j,k,2) = uold(i,j,k,2) - dt * ugradv + dt * force(i,j,k,2)
           unew(i,j,k,3) = uold(i,j,k,3) - dt * ugradw + dt * force(i,j,k,3)

           ! Add w dot grad w0 term to w.
           unew(i,j,k,3) = unew(i,j,k,3) - dt * wbar*(w0(k+1) - w0(k))/dx(3)

           ! Add w0 dot grad u term to u and w.
           wbar = HALF*(w0(k) + w0(k+1))
           unew(i,j,k,:) = unew(i,j,k,:) - dt * wbar*(sedgez(i,j,k+1,:) - sedgez(i,j,k,:))/dx(3)

      enddo
      enddo
      enddo

      umax = unew(lo(1),lo(2),lo(3),1) 
      umin = unew(lo(1),lo(2),lo(3),1) 
      vmax = unew(lo(1),lo(2),lo(3),2) 
      vmin = unew(lo(1),lo(2),lo(3),2) 
      wmax = unew(lo(1),lo(2),lo(3),3) 
      wmin = unew(lo(1),lo(2),lo(3),3) 
      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          umax = max(umax,unew(i,j,k,1))
          umin = min(umin,unew(i,j,k,1))
          vmax = max(vmax,unew(i,j,k,2))
          vmin = min(vmin,unew(i,j,k,2))
          wmax = max(wmax,unew(i,j,k,3))
          wmin = min(wmin,unew(i,j,k,3))
        enddo
      enddo
      enddo
      if (verbose .ge. 1) then
        write(6,1000) umin,umax
        write(6,1001) vmin,vmax
        write(6,1002) wmin,wmax
        write(6,1003)
      end if

1000  format('... new min/max : x-velocity       ',e17.10,2x,e17.10)
1001  format('... new min/max : y-velocity       ',e17.10,2x,e17.10)
1002  format('... new min/max : z-velocity       ',e17.10,2x,e17.10)
1003  format(' ')

   end subroutine update_velocity_3d

   subroutine mk_shalf_2d (sold,snew,ng,shalf,ng_shalf,lo,hi)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng, ng_shalf
      real (kind = dp_t), intent(in   ) ::   sold(lo(1)-ng      :,lo(2)-ng      :)  
      real (kind = dp_t), intent(in   ) ::   snew(lo(1)-ng      :,lo(2)-ng      :)  
      real (kind = dp_t), intent(  out) ::  shalf(lo(1)-ng_shalf:,lo(2)-ng_shalf:)

      integer :: i, j

      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        shalf(i,j) = HALF * (sold(i,j) + snew(i,j))
      end do
      end do

   end subroutine mk_shalf_2d

   subroutine mk_shalf_3d (sold,snew,ng,shalf,ng_shalf,lo,hi)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng, ng_shalf
      real (kind = dp_t), intent(in   ) ::   sold(lo(1)-ng      :,lo(2)-ng      :,lo(3)-ng      :)  
      real (kind = dp_t), intent(in   ) ::   snew(lo(1)-ng      :,lo(2)-ng      :,lo(3)-ng      :)  
      real (kind = dp_t), intent(  out) ::  shalf(lo(1)-ng_shalf:,lo(2)-ng_shalf:,lo(3)-ng_shalf:)

      integer :: i, j, k

      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        shalf(i,j,k) = HALF * (sold(i,j,k) + snew(i,j,k))
      end do
      end do
      end do

   end subroutine mk_shalf_3d

   subroutine mk_density_fromrhoX_2d (rho,rhoX,lo,hi,ng,verbose)

      implicit none

      integer              , intent(in) :: lo(:), hi(:), ng, verbose
      real (kind = dp_t), intent(inout) ::    rho(lo(1)-ng:,lo(2)-ng:)
      real (kind = dp_t), intent(in   ) ::   rhoX(lo(1)-ng:,lo(2)-ng:,:)  

      integer :: i, j, n
      real(dp_t) :: smin,smax

      smax = -1.d20
      smin =  1.d20

      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        rho(i,j) = ZERO
        do n = 1,nspec
          rho(i,j) = rho(i,j) + rhoX(i,j,n)
        end do
        smax = max(smax,rho(i,j))
        smin = min(smin,rho(i,j))
      end do
      end do

      if (verbose .ge. 1) then
        write(6,1000) smin,smax
      end if

1000  format('... new min/max : density           ',e15.10,2x,e15.10)

   end subroutine mk_density_fromrhoX_2d

   subroutine mk_density_fromrhoX_3d (rho,rhoX,lo,hi,ng,verbose)

      implicit none

      integer              , intent(in) :: lo(:), hi(:), ng, verbose
      real (kind = dp_t), intent(  out) ::    rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
      real (kind = dp_t), intent(in   ) ::   rhoX(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  

      integer :: i, j, k, n
      real(dp_t) :: smin,smax

      smax = -1.d20
      smin =  1.d20

      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        rho(i,j,k) = ZERO
        do n = 1,nspec
          rho(i,j,k) = rho(i,j,k) + rhoX(i,j,k,n)
        end do
        smax = max(smax,rho(i,j,k))
        smin = min(smin,rho(i,j,k))
      end do
      end do
      end do

      if (verbose .ge. 1) then
        write(6,1000) smin,smax
      end if

1000  format('... new min/max : density           ',e15.10,2x,e15.10)

   end subroutine mk_density_fromrhoX_3d

end module update_module
