module update_module

  use bl_types
  use multifab_module
  use bl_constants_module
  use fill_3d_module
  use addw0_module
  use geometry
  use variables
  use network

  implicit none

  contains

   subroutine update_scal_2d (nstart,nstop,sold,snew,umac,vmac,w0,sedgex,sedgey,force, &
                              base_old,base_new,lo,hi,ng,dx,dt)

     ! update each scalar in time.  Here, it is assumed that the edge
     ! states (sedgex and sedgey) are for the perturbational quantities.

      implicit none

      integer              , intent(in) :: nstart, nstop, lo(:), hi(:), ng
      real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,:)
      real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,:)
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:)
      real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:)
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,:)
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,:)
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,:)
      real (kind = dp_t), intent(in   ) ::   base_old(0:,:)
      real (kind = dp_t), intent(in   ) ::   base_new(0:,:)
      real (kind = dp_t), intent(in   ) :: w0(0:)
      real (kind = dp_t), intent(in   ) :: dt,dx(:)

      integer :: i, j, n, nr
      real (kind = dp_t) :: divsu,divbaseu
      real (kind = dp_t), allocatable :: base_edge(:)

      allocate(base_edge(lo(2):hi(2)+1))
      nr = size(base_old,dim=1)

      do n = nstart, nstop

        ! In case lo(2) = 0
        base_edge(lo(2)  ) = base_old(lo(2),n)
        base_edge(lo(2)+1) = HALF*(base_old(lo(2),n)+base_old(lo(2)+1,n))

        ! In case hi(2) = nr
        base_edge(hi(2)+1) = base_old(hi(2),n)
        base_edge(hi(2)  ) = HALF*(base_old(hi(2),n)+base_old(hi(2)-1,n))

        do j = max(2,lo(2)),min(nr-2,hi(2)+1)
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
        snew(:,:,rho_comp) = sold(:,:,rho_comp)
        do n = nstart, nstop
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           snew(i,j,rho_comp) = snew(i,j,rho_comp) + (snew(i,j,n)-sold(i,j,n))
        enddo
        enddo
        enddo
      end if
  
      deallocate(base_edge)

   end subroutine update_scal_2d

   subroutine update_velocity_2d (uold,unew,umac,vmac,sedgex,sedgey,force,w0, &
                                  lo,hi,ng,dx,time,dt)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in   ) ::    uold(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(  out) ::    unew(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,:)  
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,:)  
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent(in   ) ::      w0(0:)
      real (kind = dp_t), intent(in   ) :: dx(:)
      real (kind = dp_t), intent(in   ) :: time,dt

      integer :: i, j, n
      real (kind = dp_t) ubar,vbar
      real (kind = dp_t) ugradu,ugradv,ugrads
      real (kind = dp_t) :: divsu
      real (kind = dp_t) :: fac

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

   end subroutine update_velocity_2d

   subroutine update_scal_3d (nstart,nstop,sold,snew,umac,vmac,wmac,w0,w0_cart,sedgex,sedgey,sedgez,&
                              force,base_old,base_new,lo,hi,ng,dx,dt)

      implicit none

      integer, intent(in) :: nstart,nstop, lo(:), hi(:), ng
      real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(inout) ::    umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
      real (kind = dp_t), intent(inout) ::    vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
      real (kind = dp_t), intent(inout) ::    wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real (kind = dp_t), intent(in   ) ::  sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
      real (kind = dp_t), intent(in   ) ::   base_old(0:,:)
      real (kind = dp_t), intent(in   ) ::   base_new(0:,:)
      real (kind = dp_t), intent(in   ) :: w0(0:)
      real (kind = dp_t), intent(in   ) :: w0_cart(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
      real (kind = dp_t), intent(in   ) :: dt,dx(:)

      integer :: i, j, k, n, nr
      real (kind = dp_t) :: divsu,divbaseu,mult
      real (kind = dp_t), allocatable :: delta_base(:),delta_base_cart(:,:,:)
      real (kind = dp_t), allocatable :: base_cart(:,:,:)
      real (kind = dp_t), allocatable :: base_edge(:)

      allocate(base_edge(lo(3):hi(3)+1))

      nr = size(base_old,dim=1)

      if (spherical .eq. 1) then

        allocate(delta_base(nr))
        allocate(delta_base_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
        allocate(      base_cart(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))

        do n = nstart, nstop
          do k = 1,nr
            delta_base(k) = base_new(lo(3)+k-1,n) - base_old(lo(3)+k-1,n)
          end do
          call fill_3d_data(delta_base_cart,delta_base,lo,hi,dx,0)
          call fill_3d_data(base_cart,base_old(:,n),lo,hi,dx,1)
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            base_cart(i,j,lo(3)-1) = base_cart(i,j,lo(3))
            base_cart(i,j,hi(3)+1) = base_cart(i,j,hi(3))
          end do
          end do
          do k = lo(3), hi(3)
          do i = lo(1), hi(1)
            base_cart(i,lo(2)-1,k) = base_cart(i,lo(2),k)
            base_cart(i,hi(2)+1,k) = base_cart(i,hi(2),k)
          end do
          end do
          do k = lo(3), hi(3)
          do j = lo(2), hi(2)
            base_cart(lo(1)-1,j,k) = base_cart(lo(1),j,k)
            base_cart(hi(1)+1,j,k) = base_cart(hi(1),j,k)
          end do
          end do

          ! Note the umac here does NOT have w0 in it
          do k = lo(3), hi(3)
            do j = lo(2), hi(2)
            do i = lo(1), hi(1)
    
              divbaseu = HALF * (  &
                   (umac(i+1,j,k)*(base_cart(i,j,k)+base_cart(i+1,j,k)) &
                   -umac(i  ,j,k)*(base_cart(i,j,k)+base_cart(i-1,j,k)) ) / dx(1) &
                  +(vmac(i,j+1,k)*(base_cart(i,j,k)+base_cart(i,j+1,k)) &
                   -vmac(i,j  ,k)*(base_cart(i,j,k)+base_cart(i,j-1,k)) ) / dx(2) &
                  +(wmac(i,j,k+1)*(base_cart(i,j,k)+base_cart(i,j,k+1)) &
                   -wmac(i,j,k  )*(base_cart(i,j,k)+base_cart(i,j,k-1)) ) ) / dx(3) 

              snew(i,j,k,n) = sold(i,j,k,n) + delta_base_cart(i,j,k) &
                              - dt * divbaseu + dt * force(i,j,k,n)
      
            enddo
            enddo
          enddo
        end do

        deallocate(delta_base,delta_base_cart,base_cart)

        mult = ONE
        call addw0_3d_sphr(umac,vmac,wmac,w0_cart,lo,hi,dx,mult)

        do n = nstart, nstop

          ! Note the umac here DOES have w0 in it
          do k = lo(3), hi(3)
            do j = lo(2), hi(2)
            do i = lo(1), hi(1)
    
              divsu = (umac(i+1,j,k) * sedgex(i+1,j,k,n) &
                      -umac(i  ,j,k) * sedgex(i  ,j,k,n) ) / dx(1) + &
                      (vmac(i,j+1,k) * sedgey(i,j+1,k,n) &
                      -vmac(i,j  ,k) * sedgey(i,j  ,k,n) ) / dx(2) + &
                      (wmac(i,j,k+1) * sedgez(i,j,k+1,n) &
                      -wmac(i,j,k  ) * sedgez(i,j,k  ,n) ) / dx(3)

              snew(i,j,k,n) = snew(i,j,k,n) - dt * divsu

            enddo
            enddo
          enddo
        enddo

        mult = -ONE
        call addw0_3d_sphr(umac,vmac,wmac,w0_cart,lo,hi,dx,mult)

      ! not spherical
      else 

        allocate(delta_base(lo(3):hi(3)))
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
            delta_base(k) = base_new(k,n) - base_old(k,n)
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
    
              snew(i,j,k,n) = sold(i,j,k,n) + delta_base(k) &
                              - dt * (divsu + divbaseu) + dt * force(i,j,k,n)
      
            enddo
            enddo
          enddo
        end do
        deallocate(delta_base)
 
      end if


      if (nstart .eq. spec_comp .and. nstop .eq. (spec_comp+nspec-1)) then
        snew(:,:,:,rho_comp) = sold(:,:,:,rho_comp)
        do n = nstart, nstop
        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           snew(i,j,k,rho_comp) = snew(i,j,k,rho_comp) + (snew(i,j,k,n)-sold(i,j,k,n))
        enddo
        enddo
        enddo
        enddo
      end if

      deallocate(base_edge)

   end subroutine update_scal_3d

   subroutine update_velocity_3d (uold,unew,umac,vmac,wmac,sedgex,sedgey,sedgez, &
                                  force,w0,w0_cart,lo,hi,ng,dx,time,dt)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in   ) ::    uold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(  out) ::    unew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:  )
      real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:  )
      real (kind = dp_t), intent(in   ) ::    wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:  )
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real (kind = dp_t), intent(in   ) ::  sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
      real (kind = dp_t), intent(in   ) ::      w0(0:)
      real (kind = dp_t), intent(in   ) ::  w0_cart(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
      real (kind = dp_t), intent(in   ) :: dx(:)
      real (kind = dp_t), intent(in   ) :: time,dt

      integer :: i, j, k, n, nr
      real (kind = dp_t) ubar,vbar,wbar
      real (kind = dp_t) ugradu,ugradv,ugradw,ugrads
      real (kind = dp_t) :: divsu
      real (kind = dp_t) :: gradux,graduy,graduz
      real (kind = dp_t) :: gradvx,gradvy,gradvz
      real (kind = dp_t) :: gradwx,gradwy,gradwz
      real (kind = dp_t) :: w0_gradur,w0_gradvr,w0_gradwr
      real (kind = dp_t) :: gradw0
      real (kind = dp_t), allocatable :: divw0(:)
      real (kind = dp_t), allocatable :: divw0_cart(:,:,:)

      ! 1) Subtract (Utilde dot grad) Utilde term from old Utilde
      ! 2) Add forcing term to new Utilde
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


      enddo
      enddo
      enddo

      ! A) Subtract (Utilde dot er) dot grad w0 term from new Utilde.
      ! B) Subtract w0 dot grad U term from new Utilde
      if (spherical .eq. 0) then

        do k = lo(3), hi(3)

          gradw0 = (w0(k+1) - w0(k)) /dx(3)
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            wbar = HALF*(wmac(i,j,k) + wmac(i,j,k+1))
            unew(i,j,k,3) = unew(i,j,k,3) - dt * wbar * gradw0
          enddo
          enddo

          wbar = HALF*(w0(k) + w0(k+1))
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            unew(i,j,k,:) = unew(i,j,k,:) - dt * wbar*(sedgez(i,j,k+1,:) - sedgez(i,j,k,:))/dx(3)
          enddo
          enddo

        enddo

      else

      ! A) Subtract (Utilde dot er) dot grad w0 term from new Utilde.

        nr = size(w0,dim=1)-1

!       THIS IS CURRENTLY WRONG !!
!       allocate(divw0(nr))
!       do k = 1,nr
!         i = k + (lo(3)-1)
!         divw0(k) = (w0(i+1)-w0(i))/dr
!       end do
!       allocate(divw0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
!       call fill_3d_data(divw0_cart,divw0,lo,hi,dx,0)

!       do k = lo(3), hi(3)
!       do j = lo(2), hi(2)
!       do i = lo(1), hi(1)
!          ubar = HALF*(umac(i,j,k) + umac(i+1,j,k))
!          vbar = HALF*(vmac(i,j,k) + vmac(i,j+1,k))
!          wbar = HALF*(wmac(i,j,k) + wmac(i,j,k+1))
!          unew(i,j,k,1) = unew(i,j,k,1) - dt * ubar*ddx_w0dotex + vbar*ddy_w0dotex + wbar*ddz_w0dotex
!          unew(i,j,k,2) = unew(i,j,k,2) - dt * vbar*normal(i,j,k,2)*divw0_cart(i,j,k)
!          unew(i,j,k,3) = unew(i,j,k,3) - dt * wbar*normal(i,j,k,3)*divw0_cart(i,j,k)
!       enddo
!       enddo
!       enddo
!       deallocate(divw0,divw0_cart)

        ! B) Subtract (w0 dot grad) U term from new Utilde

        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           gradux = (sedgex(i+1,j,k,1) - sedgex(i,j,k,1))/dx(1)
           gradvx = (sedgex(i+1,j,k,2) - sedgex(i,j,k,2))/dx(1)
           gradwx = (sedgex(i+1,j,k,3) - sedgex(i,j,k,3))/dx(1)
           graduy = (sedgey(i,j+1,k,1) - sedgey(i,j,k,1))/dx(2)
           gradvy = (sedgey(i,j+1,k,2) - sedgey(i,j,k,2))/dx(2)
           gradwy = (sedgey(i,j+1,k,3) - sedgey(i,j,k,3))/dx(2)
           graduz = (sedgez(i,j,k+1,1) - sedgez(i,j,k,1))/dx(3)
           gradvz = (sedgez(i,j,k+1,2) - sedgez(i,j,k,2))/dx(3)
           gradwz = (sedgez(i,j,k+1,3) - sedgez(i,j,k,3))/dx(3)

           w0_gradur = gradux * w0_cart(i,j,k,1) + graduy * w0_cart(i,j,k,2) + graduz * w0_cart(i,j,k,3)
           w0_gradvr = gradvx * w0_cart(i,j,k,1) + gradvy * w0_cart(i,j,k,2) + gradvz * w0_cart(i,j,k,3)
           w0_gradwr = gradwx * w0_cart(i,j,k,1) + gradwy * w0_cart(i,j,k,2) + gradwz * w0_cart(i,j,k,3)

           unew(i,j,k,1) = unew(i,j,k,1) - dt * w0_gradur
           unew(i,j,k,2) = unew(i,j,k,2) - dt * w0_gradvr
           unew(i,j,k,2) = unew(i,j,k,3) - dt * w0_gradwr

        enddo
        enddo
        enddo

      end if
   
   end subroutine update_velocity_3d

end module update_module
