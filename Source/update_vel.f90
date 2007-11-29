module update_vel_module

  use bl_types
  use multifab_module
  use bl_constants_module
  use geometry
  use sponge_module
  use variables
  use network

  implicit none

  private
  public :: update_velocity_2d, update_velocity_3d

  contains

   subroutine update_velocity_2d (uold,unew,umac,vmac,sedgex,sedgey,force,w0,w0_force, &
                                  lo,hi,ng,dx,dt,sponge,do_sponge)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in   ) ::     uold(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(  out) ::     unew(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in   ) ::     umac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::     vmac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::   sedgex(lo(1)   :,lo(2)   :,:)  
      real (kind = dp_t), intent(in   ) ::   sedgey(lo(1)   :,lo(2)   :,:)  
      real (kind = dp_t), intent(in   ) ::    force(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent(in   ) ::   sponge(lo(1)   :,lo(2)   :  )
      real (kind = dp_t), intent(in   ) ::       w0(0:)
      real (kind = dp_t), intent(in   ) :: w0_force(0:)
      real (kind = dp_t), intent(in   ) :: dx(:)
      real (kind = dp_t), intent(in   ) :: dt
      logical           , intent(in   ) :: do_sponge

      integer :: i, j
      real (kind = dp_t) ubar,vbar
      real (kind = dp_t) ugradu,ugradv

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

           ! Add in the pi0 term.
           unew(i,j,2) = unew(i,j,2) - dt * w0_force(j)

           ! Add the sponge
           if (do_sponge) unew(i,j,:) = unew(i,j,:) * sponge(i,j)

      enddo
      enddo

   end subroutine update_velocity_2d

   subroutine update_velocity_3d (uold,unew,umac,vmac,wmac,sedgex,sedgey,sedgez, &
                                  force,w0,w0_cart,w0_force,w0_force_cart,lo,hi,ng,dx,dt, &
                                  sponge,do_sponge)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(in   ) ::     uold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(  out) ::     unew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(in   ) ::     umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:  )
      real (kind = dp_t), intent(in   ) ::     vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:  )
      real (kind = dp_t), intent(in   ) ::     wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:  )
      real (kind = dp_t), intent(in   ) ::   sedgex(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real (kind = dp_t), intent(in   ) ::   sedgey(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real (kind = dp_t), intent(in   ) ::   sedgez(lo(1)   :,lo(2)   :,lo(3)   :,:)
      real (kind = dp_t), intent(in   ) ::    force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
      real (kind = dp_t), intent(in   ) ::   sponge(lo(1)   :,lo(2)   :,lo(3)   :  ) 
      real (kind = dp_t), intent(in   ) ::       w0(0:)
      real (kind = dp_t), intent(in   ) ::  w0_cart(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
      real (kind = dp_t), intent(in   ) :: w0_force(0:)
      real (kind = dp_t), intent(in   ) ::  w0_force_cart(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
      real (kind = dp_t), intent(in   ) :: dx(:)
      real (kind = dp_t), intent(in   ) :: dt
      logical           , intent(in   ) :: do_sponge

      integer :: i, j, k
      real (kind = dp_t) ubar,vbar,wbar
      real (kind = dp_t) ugradu,ugradv,ugradw
      real (kind = dp_t) :: gradux,graduy,graduz
      real (kind = dp_t) :: gradvx,gradvy,gradvz
      real (kind = dp_t) :: gradwx,gradwy,gradwz
      real (kind = dp_t) :: w0_gradur,w0_gradvr,w0_gradwr
      real (kind = dp_t) :: gradw0

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

            ! Add in the pi0 term.
            unew(i,j,k,3) = unew(i,j,k,3) - dt * w0_force(k)

            ! Add the sponge
            if (do_sponge) unew(i,j,k,:) = unew(i,j,k,:) * sponge(i,j,k)

          enddo
          enddo

        enddo

      else

      ! A) Subtract (Utilde dot er) dot grad w0 term from new Utilde.

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
           unew(i,j,k,3) = unew(i,j,k,3) - dt * w0_gradwr

           ! Add in the pi0 term.
           unew(i,j,k,:) = unew(i,j,k,:) - dt * w0_force_cart(i,j,k,:)

           ! Add the sponge
           if (do_sponge) unew(i,j,k,:) = unew(i,j,k,:) * sponge(i,j,k)

          enddo
         enddo
        enddo

      end if
   
   end subroutine update_velocity_3d

end module update_vel_module
