module update_module

  ! do the conservative updating of the scalars and velocity
  ! This is used both in step 2 of Almgren et al. 2006, paper II
  ! (ABRZ2) (when pred_vs_corr = 1) and in step 4 (pred_vs_corr = 2).

  use bl_types
  use multifab_module
  use eos_module

  implicit none

  real (kind = dp_t), private, parameter :: HALF  = 0.5_dp_t

  contains

   subroutine update_density_2d (sold,snew,umac,vmac,w0,sedgex,sedgey,rhohalf, &
                                 rho0_old,rho0_new,lo,hi,ng,dx,dt,pred_vs_corr,div_coeff,div_coeff_half,&
                                 verbose)

     ! in the density update, the edge states (sedgex and sedgey) are
     ! taken to be the perturbational density (rhopert = rho - rho0)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng, pred_vs_corr, verbose
      real (kind = dp_t), intent(in   ) ::     sold(lo(1)-ng:,lo(2)-ng:)  
      real (kind = dp_t), intent(  out) ::     snew(lo(1)-ng:,lo(2)-ng:)  
      real (kind = dp_t), intent(in   ) ::     umac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::     vmac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::   sedgex(lo(1)   :,lo(2)   :)  
      real (kind = dp_t), intent(in   ) ::   sedgey(lo(1)   :,lo(2)   :)  
      real (kind = dp_t), intent(inout) ::  rhohalf(lo(1)- 1:,lo(2)- 1:)
      real (kind = dp_t), intent(in   ) :: rho0_old(lo(2)   :)
      real (kind = dp_t), intent(in   ) :: rho0_new(lo(2)   :)
      real (kind = dp_t), intent(in   ) :: div_coeff(lo(2)   :)
      real (kind = dp_t), intent(in   ) :: div_coeff_half(lo(2)   :)
      real (kind = dp_t), intent(in   ) :: w0(lo(2):)
      real (kind = dp_t), intent(in   ) :: dx(:)
      real (kind = dp_t), intent(in   ) :: dt

      integer :: i, j
      real (kind = dp_t) :: divsu,divbaseu
      real (kind = dp_t) :: smin,smax
      real (kind = dp_t), allocatable :: rho0_edge(:)

      allocate(rho0_edge(lo(2):hi(2)+1))

      smax = -1.e20
      smin =  1.e20
           
      rho0_edge(lo(2)  ) = rho0_old(lo(2))
      rho0_edge(hi(2)+1) = rho0_old(hi(2))
      
      rho0_edge(lo(2)+1) = HALF*(rho0_old(lo(2))+rho0_old(lo(2)+1))
      rho0_edge(hi(2)  ) = HALF*(rho0_old(hi(2))+rho0_old(hi(2)-1))

      do j = lo(2)+2,hi(2)-1
         rho0_edge(j) = 7.d0/12.d0 * (rho0_old(j  ) + rho0_old(j-1)) &
                       -1.d0/12.d0 * (rho0_old(j+1) + rho0_old(j-2))
      end do

      if (pred_vs_corr .eq. 1) then

         ! see ABRZ2, step 2

        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
  
          divsu = (umac(i+1,j) * sedgex(i+1,j) &
                  -umac(i  ,j) * sedgex(i  ,j) ) / dx(1) + &
                  (vmac(i,j+1) * sedgey(i,j+1) &
                  -vmac(i,j  ) * sedgey(i,j  ) ) / dx(2)

          divbaseu = (umac(i+1,j) - umac(i,j) ) * rho0_old(j) / dx(1) &
                    +(vmac(i,j+1) * rho0_edge(j+1) - vmac(i,j) * rho0_edge(j) ) / dx(2)

          snew(i,j) = sold(i,j) - dt * (divsu + divbaseu)
          snew(i,j) = max(snew(i,j),rho0_old(hi(2)))
  
          smax = max(smax,snew(i,j))
          smin = min(smin,snew(i,j))

        enddo
        enddo

      else if (pred_vs_corr .eq. 2) then

         ! see ABRZ2, step 4

        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
  
          divsu = (umac(i+1,j) * sedgex(i+1,j) &
                  -umac(i  ,j) * sedgex(i  ,j) ) / dx(1) + &
                  ((vmac(i,j+1)+w0(j+1)) * sedgey(i,j+1) &
                  -(vmac(i,j  )+w0(j  )) * sedgey(i,j  ) ) / dx(2)

          divbaseu = (umac(i+1,j) - umac(i,j) ) * rho0_old(j) / dx(1) &
                    +(vmac(i,j+1) * rho0_edge(j+1) - vmac(i,j) * rho0_edge(j) ) / dx(2)

          snew(i,j) = sold(i,j) + (rho0_new(j) - rho0_old(j)) - dt * (divsu + divbaseu)
          snew(i,j) = max(snew(i,j),rho0_new(hi(2)))
  
          smax = max(smax,snew(i,j))
          smin = min(smin,snew(i,j))

        enddo
        enddo

      end if

      if (verbose .ge. 1) then
        print *,'NEW MIN/MAX OF RHO ',smin,smax
        print *,' '
      end if

      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        rhohalf(i,j) = HALF * (sold(i,j) + snew(i,j))
      end do
      end do

      deallocate(rho0_edge)

   end subroutine update_density_2d

   subroutine update_rhoh_2d (sold,snew,umac,vmac,w0,sedgex,sedgey,force, &
                              rhoh0_old,rhoh0_new,lo,hi,ng,dx,dt,pred_vs_corr,&
                              verbose)

     ! update the enthalpy in time.  Here, it is assumed that the edge
     ! states (sedgex and sedgey) are for the enthalpy perturbation
     ! (rhohpert = rhoh - rhoh0).  

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng, pred_vs_corr, verbose
      real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:)  
      real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:)  
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :)  
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :)  
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::   rhoh0_old(lo(2)   :)
      real (kind = dp_t), intent(in   ) ::   rhoh0_new(lo(2)   :)
      real (kind = dp_t), intent(in   ) :: w0(lo(2):)
      real (kind = dp_t), intent(in   ) :: dt,dx(:)

      integer :: i, j, n
      real (kind = dp_t) :: divsu,divbaseu
      real (kind = dp_t) :: smin,smax
      real (kind = dp_t), allocatable :: rhoh0_edge(:)

      allocate(rhoh0_edge(lo(2):hi(2)+1))

      do_diag = .false.

      smax = -1.d20
      smin =  1.d20

      rhoh0_edge(lo(2)  ) = rhoh0_old(lo(2))
      rhoh0_edge(hi(2)+1) = rhoh0_old(hi(2))
      
      rhoh0_edge(lo(2)+1) = HALF*(rhoh0_old(lo(2))+rhoh0_old(lo(2)+1))
      rhoh0_edge(hi(2)  ) = HALF*(rhoh0_old(hi(2))+rhoh0_old(hi(2)-1))

      do j = lo(2)+2,hi(2)-1
         rhoh0_edge(j) = 7.d0/12.d0 * (rhoh0_old(j  ) + rhoh0_old(j-1)) &
                        -1.d0/12.d0 * (rhoh0_old(j+1) + rhoh0_old(j-2))
      end do
      
           
      if (pred_vs_corr .eq. 1) then

         ! see ABRZ2, step 2

        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
  
          divsu = (umac(i+1,j) * sedgex(i+1,j) &
                  -umac(i  ,j) * sedgex(i  ,j) ) / dx(1) + &
                  (vmac(i,j+1) * sedgey(i,j+1) &
                  -vmac(i,j  ) * sedgey(i,j  ) ) / dx(2)

          divbaseu = (umac(i+1,j) - umac(i,j) ) * rhoh0_old(j) / dx(1) &
                    +(vmac(i,j+1) * rhoh0_edge(j+1) - vmac(i,j) * rhoh0_edge(j) ) / dx(2)
  
          snew(i,j) = sold(i,j) - dt * (divsu + divbaseu) + dt * force(i,j)
          snew(i,j) = max(snew(i,j),rhoh0_old(hi(2)))
  
          smax = max(smax,snew(i,j))
          smin = min(smin,snew(i,j))
  
        enddo
        enddo

      else if (pred_vs_corr .eq. 2) then

         ! see ABRZ2, step 4

        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
  
          divsu = (umac(i+1,j) * sedgex(i+1,j) &
                  -umac(i  ,j) * sedgex(i  ,j) ) / dx(1) + &
                 ((vmac(i,j+1)+w0(j+1)) * sedgey(i,j+1) &
                 -(vmac(i,j  )+w0(j  )) * sedgey(i,j  ) ) / dx(2)

          divbaseu = (umac(i+1,j) - umac(i,j) ) * rhoh0_old(j) / dx(1) &
                    +(vmac(i,j+1) * rhoh0_edge(j+1) - vmac(i,j) * rhoh0_edge(j) ) / dx(2)

          snew(i,j) = sold(i,j) + (rhoh0_new(j) - rhoh0_old(j)) - dt * (divsu + divbaseu) + dt * force(i,j)

          snew(i,j) = max(snew(i,j),rhoh0_new(hi(2)))
  
          smax = max(smax,snew(i,j))
          smin = min(smin,snew(i,j))
  
        enddo
        enddo

      end if

      if (verbose .ge. 1) then
        print *,'NEW MIN/MAX OF (RHO H) ',smin,smax
        print *,' '
      end if

      deallocate(rhoh0_edge)

   end subroutine update_rhoh_2d

   subroutine update_species_2d (sold,snew,umac,vmac,w0,sedgex,sedgey,force, &
                                 lo,hi,ng,dx,dt,pred_vs_corr,verbose)

     ! Update the species in time.  

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng, pred_vs_corr, verbose
      real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:)  
      real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:)  
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :)  
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :)  
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) :: w0(lo(2):)
      real (kind = dp_t), intent(in   ) :: dt,dx(:)

      integer :: i, j, n
      real (kind = dp_t) :: divsu,divbaseu
      real (kind = dp_t) :: somin,somax,snmin,snmax

      somax = -1.d20
      somin =  1.d20

      snmax = -1.d20
      snmin =  1.d20

      if (pred_vs_corr .eq. 1) then

         ! see ABRZ2, step 2

        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
  
          divsu = (umac(i+1,j) * sedgex(i+1,j) &
                  -umac(i  ,j) * sedgex(i  ,j) ) / dx(1) + &
                  (vmac(i,j+1) * sedgey(i,j+1) &
                  -vmac(i,j  ) * sedgey(i,j  ) ) / dx(2)
  
          snew(i,j) = sold(i,j) - dt * divsu + dt * force(i,j)
          if (abs(sold(i,j)) .gt. 1000.) print *,'SPEC ',i,j,sold(i,j)
  
          somax = max(somax,sold(i,j))
          somin = min(somin,sold(i,j))
  
          snmax = max(snmax,snew(i,j))
          snmin = min(snmin,snew(i,j))
  
        enddo
        enddo

      else if (pred_vs_corr .eq. 2) then

         ! see ABRZ2, step 4

        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
  
          divsu = (umac(i+1,j) * sedgex(i+1,j) &
                  -umac(i  ,j) * sedgex(i  ,j) ) / dx(1) + &
                 ((vmac(i,j+1)+w0(j+1)) * sedgey(i,j+1) &
                 -(vmac(i,j  )+w0(j  )) * sedgey(i,j  ) ) / dx(2)

          snew(i,j) = sold(i,j) - dt * divsu + dt * force(i,j)
  
          somax = max(somax,sold(i,j))
          somin = min(somin,sold(i,j))
  
          snmax = max(snmax,snew(i,j))
          snmin = min(snmin,snew(i,j))
  
        enddo
        enddo

      end if

      if (verbose .ge. 1) then
        print *,'OLD MIN/MAX OF SPECIES ',somin,somax
        print *,'NEW MIN/MAX OF SPECIES ',snmin,snmax
        print *,' '
      end if

   end subroutine update_species_2d

   subroutine update_velocity_2d (uold,unew,rhoold,rhonew,umac,vmac,sedgex,sedgey,force,w0, &
                                  lo,hi,ng,dx,time,dt,do_mom,verbose)


     ! update the velocity in time (to get the provisional velocity that
     ! does not yet satisfy the divergence constraint).  This is the first
     ! part of step 5 in ABRZ2.

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng, verbose
      real (kind = dp_t), intent(in   ) ::    uold(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(  out) ::    unew(lo(1)-ng:,lo(2)-ng:,:)  
      real (kind = dp_t), intent(in   ) ::  rhoold(lo(1)-ng:,lo(2)-ng:)  
      real (kind = dp_t), intent(in   ) ::  rhonew(lo(1)-ng:,lo(2)-ng:)  
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:)  
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,:)  
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,:)  
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,:)  
      real (kind = dp_t), intent(in   ) ::      w0(          lo(2)   :)  
      real (kind = dp_t), intent(in   ) :: dx(:)
      real (kind = dp_t), intent(in   ) :: time,dt
      logical           , intent(in   ) :: do_mom

      integer :: i, j, n
      real (kind = dp_t) ubar,vbar
      real (kind = dp_t) ugradu,ugradv,ugrads
      real (kind = dp_t) :: divsu
      real (kind = dp_t) :: smin,smax,umin,umax,vmin,vmax
      real (kind = dp_t) :: fac

      if (do_mom) then
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)

             divsu = (umac(i+1,j) * sedgex(i+1,j,1) &
                     -umac(i  ,j) * sedgex(i  ,j,1) ) / dx(1) + &
                     (vmac(i,j+1) * sedgey(i,j+1,1) &
                     -vmac(i,j  ) * sedgey(i,j  ,1) ) / dx(2)
             unew(i,j,1) = rhoold(i,j)*uold(i,j,1) - dt * divsu + dt * force(i,j,1)

             divsu = (umac(i+1,j) * sedgex(i+1,j,2) &
                     -umac(i  ,j) * sedgex(i  ,j,2) ) / dx(1) + &
                     (vmac(i,j+1) * sedgey(i,j+1,2) &
                     -vmac(i,j  ) * sedgey(i,j  ,2) ) / dx(2)
             unew(i,j,2) = rhoold(i,j)*uold(i,j,2) - dt * divsu + dt * force(i,j,2)

             unew(i,j,1) = unew(i,j,1) / rhonew(i,j)
             unew(i,j,2) = unew(i,j,2) / rhonew(i,j)

        enddo
        enddo
      else
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
      end if

      umax = uold(lo(1),lo(2),1) 
      umin = uold(lo(1),lo(2),1) 
      vmax = uold(lo(1),lo(2),2) 
      vmin = uold(lo(1),lo(2),2) 
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          umax = max(umax,uold(i,j,1))
          umin = min(umin,uold(i,j,1))
          vmax = max(vmax,uold(i,j,2))
          vmin = min(vmin,uold(i,j,2))
        enddo
      enddo
      print *,'MIN/MAX OF UOLD ',umin,umax
      print *,'MIN/MAX OF VOLD ',vmin,vmax

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
        print *,'MIN/MAX OF UNEW ',umin,umax
        print *,'MIN/MAX OF VNEW ',vmin,vmax
        print *,' '
      end if

   end subroutine update_velocity_2d

   subroutine update_density_3d (sold,snew,umac,vmac,wmac,w0,sedgex,sedgey,sedgez,rhohalf, &
                                 rho0_old,rho0_new,lo,hi,ng,dx,dt,pred_vs_corr,div_coeff,div_coeff_half,&
                                 verbose)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng, pred_vs_corr, verbose
      real (kind = dp_t), intent(in   ) ::     sold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
      real (kind = dp_t), intent(  out) ::     snew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
      real (kind = dp_t), intent(in   ) ::     umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
      real (kind = dp_t), intent(in   ) ::     vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
      real (kind = dp_t), intent(in   ) ::     wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
      real (kind = dp_t), intent(in   ) ::   sedgex(lo(1)   :,lo(2)   :,lo(3)   :)  
      real (kind = dp_t), intent(in   ) ::   sedgey(lo(1)   :,lo(2)   :,lo(3)   :)  
      real (kind = dp_t), intent(in   ) ::   sedgez(lo(1)   :,lo(2)   :,lo(3)   :)  
      real (kind = dp_t), intent(inout) ::  rhohalf(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
      real (kind = dp_t), intent(in   ) :: rho0_old(lo(3)   :)
      real (kind = dp_t), intent(in   ) :: rho0_new(lo(3)   :)
      real (kind = dp_t), intent(in   ) :: div_coeff(lo(3)   :)
      real (kind = dp_t), intent(in   ) :: div_coeff_half(lo(3)   :)
      real (kind = dp_t), intent(in   ) :: w0(lo(3):)
      real (kind = dp_t), intent(in   ) :: dx(:)
      real (kind = dp_t), intent(in   ) :: dt

      integer :: i, j, k
      real (kind = dp_t) :: divsu,divbaseu
      real (kind = dp_t) :: smin,smax
      real (kind = dp_t), allocatable :: rho0_edge(:)

      allocate(rho0_edge(lo(3):hi(3)+1))

      smax = -1.e20
      smin =  1.e20
           
      rho0_edge(lo(3)  ) = rho0_old(lo(3))
      rho0_edge(hi(3)+1) = rho0_old(hi(3))
      
      rho0_edge(lo(3)+1) = HALF*(rho0_old(lo(3))+rho0_old(lo(3)+1))
      rho0_edge(hi(3)  ) = HALF*(rho0_old(hi(3))+rho0_old(hi(3)-1))

      do k = lo(3)+2,hi(3)-1
         rho0_edge(k) = 7.d0/12.d0 * (rho0_old(k  ) + rho0_old(k-1)) &
                       -1.d0/12.d0 * (rho0_old(k+1) + rho0_old(k-2))
      end do

      if (pred_vs_corr .eq. 1) then

        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
  
          divsu = (umac(i+1,j,k) * sedgex(i+1,j,k) &
                  -umac(i  ,j,k) * sedgex(i  ,j,k) ) / dx(1) + &
                  (vmac(i,j+1,k) * sedgey(i,j+1,k) &
                  -vmac(i,j  ,k) * sedgey(i,j  ,k) ) / dx(2) + &
                  (wmac(i,j,k+1) * sedgez(i,j,k+1) &
                  -wmac(i,j,k  ) * sedgez(i,j,k  ) ) / dx(3)

          divbaseu = (umac(i+1,j,k) - umac(i,j,k) ) * rho0_old(k) / dx(1) &
                    +(vmac(i,j+1,k) - vmac(i,j,k) ) * rho0_old(k) / dx(2) &
                    +(wmac(i,j,k+1) * rho0_edge(k+1) - wmac(i,j,k) * rho0_edge(k) ) / dx(3)

          snew(i,j,k) = sold(i,j,k) - dt * (divsu + divbaseu)

          snew(i,j,k) = max(snew(i,j,k),rho0_old(hi(3)))
  
          smax = max(smax,snew(i,j,k))
          smin = min(smin,snew(i,j,k))

        enddo
        enddo
        enddo

      else if (pred_vs_corr .eq. 2) then

        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
  
          divsu = (umac(i+1,j,k) * sedgex(i+1,j,k) &
                  -umac(i  ,j,k) * sedgex(i  ,j,k) ) / dx(1) + &
                  (vmac(i,j+1,k) * sedgey(i,j+1,k) &
                  -vmac(i,j  ,k) * sedgey(i,j  ,k) ) / dx(2) + &
                  ((wmac(i,j,k+1)+w0(k+1)) * sedgez(i,j,k+1) &
                  -(wmac(i,j,k  )+w0(k  )) * sedgez(i,j,k  ) ) / dx(3)

          divbaseu = (umac(i+1,j,k) - umac(i,j,k) ) * rho0_old(k) / dx(1) &
                    +(vmac(i,j+1,k) - vmac(i,j,k) ) * rho0_old(k) / dx(2) &
                    +(wmac(i,j,k+1) * rho0_edge(k+1) - wmac(i,j,k) * rho0_edge(k) ) / dx(3)

          snew(i,j,k) = sold(i,j,k) + (rho0_new(k) - rho0_old(k)) - dt * (divsu + divbaseu)
          snew(i,j,k) = max(snew(i,j,k),rho0_new(hi(3)))
  
          smax = max(smax,snew(i,j,k))
          smin = min(smin,snew(i,j,k))
          

        enddo
        enddo
        enddo

      end if

      if (verbose .ge. 1) then
        print *,'NEW MIN/MAX OF RHO ',smin,smax
        print *,' '
      end if

      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
        rhohalf(i,j,k) = HALF * (sold(i,j,k) + snew(i,j,k))
      end do
      end do
      end do

      deallocate(rho0_edge)

   end subroutine update_density_3d

   subroutine update_rhoh_3d (sold,snew,umac,vmac,wmac,w0,sedgex,sedgey,sedgez,force, &
                              rhoh0_old,rhoh0_new,lo,hi,ng,dx,dt,pred_vs_corr,verbose)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng, pred_vs_corr, verbose
      real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)  
      real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)  
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
      real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
      real (kind = dp_t), intent(in   ) ::    wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,lo(3)   :)  
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,lo(3)   :)  
      real (kind = dp_t), intent(in   ) ::  sedgez(lo(1)   :,lo(2)   :,lo(3)   :)  
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
      real (kind = dp_t), intent(in   ) ::   rhoh0_old(lo(3)   :)
      real (kind = dp_t), intent(in   ) ::   rhoh0_new(lo(3)   :)
      real (kind = dp_t), intent(in   ) :: w0(lo(3):)
      real (kind = dp_t), intent(in   ) :: dt,dx(:)

      integer :: i, j, k, n
      real (kind = dp_t) :: divsu,divbaseu
      real (kind = dp_t) :: smin,smax
      real (kind = dp_t), allocatable :: rhoh0_edge(:)

      allocate(rhoh0_edge(lo(3):hi(3)+1))

      do_diag = .false.

      smax = -1.d20
      smin =  1.d20

      rhoh0_edge(lo(3)  ) = rhoh0_old(lo(3))
      rhoh0_edge(hi(3)+1) = rhoh0_old(hi(3))
      
      rhoh0_edge(lo(3)+1) = HALF*(rhoh0_old(lo(3))+rhoh0_old(lo(3)+1))
      rhoh0_edge(hi(3)  ) = HALF*(rhoh0_old(hi(3))+rhoh0_old(hi(3)-1))

      do k = lo(3)+2,hi(3)-1
         rhoh0_edge(k) = 7.d0/12.d0 * (rhoh0_old(k  ) + rhoh0_old(k-1)) &
                        -1.d0/12.d0 * (rhoh0_old(k+1) + rhoh0_old(k-2))
      end do
           
      if (pred_vs_corr .eq. 1) then

        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
  
          divsu = (umac(i+1,j,k) * sedgex(i+1,j,k) &
                  -umac(i  ,j,k) * sedgex(i  ,j,k) ) / dx(1) + &
                  (vmac(i,j+1,k) * sedgey(i,j+1,k) &
                  -vmac(i,j  ,k) * sedgey(i,j  ,k) ) / dx(2) + &
                  (wmac(i,j,k+1) * sedgez(i,j,k+1) &
                  -wmac(i,j,k  ) * sedgez(i,j,k  ) ) / dx(3)

          divbaseu = (umac(i+1,j,k) - umac(i,j,k) ) * rhoh0_old(k) / dx(1) &
                    +(vmac(i,j+1,k) - vmac(i,j,k) ) * rhoh0_old(k) / dx(2) &
                    +(wmac(i,j,k+1) * rhoh0_edge(k+1) - wmac(i,j,k) * rhoh0_edge(k) ) / dx(3)

          snew(i,j,k) = sold(i,j,k) - dt * (divsu + divbaseu) + dt * force(i,j,k)
          snew(i,j,k) = max(snew(i,j,k),rhoh0_old(hi(3)))
  
          smax = max(smax,snew(i,j,k))
          smin = min(smin,snew(i,j,k))
  
        enddo
        enddo
        enddo

      else if (pred_vs_corr .eq. 2) then

        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
  
          divsu = (umac(i+1,j,k) * sedgex(i+1,j,k) &
                  -umac(i  ,j,k) * sedgex(i  ,j,k) ) / dx(1) + &
                  (vmac(i,j+1,k) * sedgey(i,j+1,k) &
                  -vmac(i,j  ,k) * sedgey(i,j  ,k) ) / dx(2) + &
                 ((wmac(i,j,k+1)+w0(k+1)) * sedgez(i,j,k+1) &
                 -(wmac(i,j,k  )+w0(k  )) * sedgez(i,j,k  ) ) / dx(3)

          divbaseu = (umac(i+1,j,k) - umac(i,j,k) ) * rhoh0_old(k) / dx(1) &
                    +(vmac(i,j+1,k) - vmac(i,j,k) ) * rhoh0_old(k) / dx(2) &
                    +(wmac(i,j,k+1) * rhoh0_edge(k+1) - wmac(i,j,k) * rhoh0_edge(k) ) / dx(3)

          snew(i,j,k) = sold(i,j,k) + (rhoh0_new(k) - rhoh0_old(k)) - dt * (divsu + divbaseu) + dt * force(i,j,k)

          snew(i,j,k) = max(snew(i,j,k),rhoh0_new(hi(3)))
  
          smax = max(smax,snew(i,j,k))
          smin = min(smin,snew(i,j,k))
  
        enddo
        enddo
        enddo

      end if

      if (verbose .ge. 1) then
        print *,'NEW MIN/MAX OF (RHO H) ',smin,smax
        print *,' '
      end if

      deallocate(rhoh0_edge)

   end subroutine update_rhoh_3d

   subroutine update_species_3d (sold,snew,umac,vmac,wmac,w0,sedgex,sedgey,sedgez,force, &
                                 lo,hi,ng,dx,dt,pred_vs_corr,verbose)

     ! Update the species in time.  

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng, pred_vs_corr, verbose
      real (kind = dp_t), intent(in   ) ::    sold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)  
      real (kind = dp_t), intent(  out) ::    snew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)  
      real (kind = dp_t), intent(in   ) ::    umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
      real (kind = dp_t), intent(in   ) ::    vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
      real (kind = dp_t), intent(in   ) ::    wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
      real (kind = dp_t), intent(in   ) ::  sedgex(lo(1)   :,lo(2)   :,lo(3)   :)  
      real (kind = dp_t), intent(in   ) ::  sedgey(lo(1)   :,lo(2)   :,lo(3)   :)  
      real (kind = dp_t), intent(in   ) ::  sedgez(lo(1)   :,lo(2)   :,lo(3)   :)  
      real (kind = dp_t), intent(in   ) ::   force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)  
      real (kind = dp_t), intent(in   ) :: w0(lo(3):)
      real (kind = dp_t), intent(in   ) :: dt,dx(:)

      integer :: i, j, k, n
      real (kind = dp_t) :: divsu
      real (kind = dp_t) :: smin,smax

      smax = -1.d20
      smin =  1.d20

      if (pred_vs_corr .eq. 1) then

         ! see ABRZ2, step 2

        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
  
          divsu = (umac(i+1,j,k) * sedgex(i+1,j,k) &
                  -umac(i  ,j,k) * sedgex(i  ,j,k) ) / dx(1) + &
                  (vmac(i,j+1,k) * sedgey(i,j+1,k) &
                  -vmac(i,j  ,k) * sedgey(i,j  ,k) ) / dx(2) + &
                  (wmac(i,j,k+1) * sedgez(i,j,k+1) &
                  -wmac(i,j,k  ) * sedgez(i,j,k  ) ) / dx(3)
  
          snew(i,j,k) = sold(i,j,k) - dt * divsu + dt * force(i,j,k)
  
          smax = max(smax,snew(i,j,k))
          smin = min(smin,snew(i,j,k))
  
        enddo
        enddo
        enddo

      else if (pred_vs_corr .eq. 2) then

         ! see ABRZ2, step 4

        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)
  
          divsu = (umac(i+1,j,k) * sedgex(i+1,j,k) &
                  -umac(i  ,j,k) * sedgex(i  ,j,k) ) / dx(1) + &
                  (vmac(i,j+1,k) * sedgey(i,j+1,k) &
                  -vmac(i,j  ,k) * sedgey(i,j  ,k) ) / dx(2) + &
                 ((wmac(i,j,k+1)+w0(k+1)) * sedgez(i,j,k+1) &
                 -(wmac(i,j,k  )+w0(k  )) * sedgez(i,j,k  ) ) / dx(3)

          snew(i,j,k) = sold(i,j,k) - dt * divsu + dt * force(i,j,k)
  
          smax = max(smax,snew(i,j,k))
          smin = min(smin,snew(i,j,k))
  
        enddo
        enddo
        enddo

      end if

      if (verbose .ge. 1) then
        print *,'NEW MIN/MAX OF SPECIES ',smin,smax
        print *,' '
      end if

   end subroutine update_species_3d

   subroutine update_velocity_3d (uold,unew,rhoold,rhonew,umac,vmac,wmac,sedgex,sedgey,sedgez, &
                                  force,w0,lo,hi,ng,dx,time,dt,do_mom,verbose)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng, verbose
      real (kind = dp_t), intent(in   ) ::    uold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(  out) ::    unew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(in   ) ::  rhoold(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:  )
      real (kind = dp_t), intent(in   ) ::  rhonew(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:  )
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
      logical           , intent(in   ) :: do_mom

      integer :: i, j, k, n
      real (kind = dp_t) ubar,vbar,wbar
      real (kind = dp_t) ugradu,ugradv,ugradw,ugrads
      real (kind = dp_t) :: divsu
      real (kind = dp_t) :: smin,smax,umin,umax,vmin,vmax,wmin,wmax
      real (kind = dp_t) :: fac

      if (do_mom) then
        do n = 1,3
        do k = lo(3), hi(3)
        do j = lo(2), hi(2)
        do i = lo(1), hi(1)

             divsu = (umac(i+1,j,k) * sedgex(i+1,j,k,n) &
                     -umac(i  ,j,k) * sedgex(i  ,j,k,n) ) / dx(1) + &
                     (vmac(i,j+1,k) * sedgey(i,j+1,k,n) &
                     -vmac(i,j  ,k) * sedgey(i,j  ,k,n) ) / dx(2) + &
                     (wmac(i,j,k+1) * sedgez(i,j,k+1,n) &
                     -wmac(i,j,k  ) * sedgez(i,j,k  ,n) ) / dx(3)
             unew(i,j,k,n) = rhoold(i,j,k)*uold(i,j,k,n) - dt * divsu + dt * force(i,j,k,n)

             unew(i,j,k,n) = unew(i,j,k,n) / rhonew(i,j,k)

        enddo
        enddo
        enddo
        enddo
      else
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
      end if

      umax = uold(lo(1),lo(2),lo(3),1) 
      umin = uold(lo(1),lo(2),lo(3),1) 
      vmax = uold(lo(1),lo(2),lo(3),2) 
      vmin = uold(lo(1),lo(2),lo(3),2) 
      wmax = uold(lo(1),lo(2),lo(3),3) 
      wmin = uold(lo(1),lo(2),lo(3),3) 
      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
        do i = lo(1), hi(1)
          umax = max(umax,uold(i,j,k,1))
          umin = min(umin,uold(i,j,k,1))
          vmax = max(vmax,uold(i,j,k,2))
          vmin = min(vmin,uold(i,j,k,2))
          wmax = max(wmax,uold(i,j,k,3))
          wmin = min(wmin,uold(i,j,k,3))
        enddo
      enddo
      enddo
      print *,'MIN/MAX OF UOLD ',umin,umax
      print *,'MIN/MAX OF VOLD ',vmin,vmax
      print *,'MIN/MAX OF WOLD ',wmin,wmax

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
        print *,'MIN/MAX OF UNEW ',umin,umax
        print *,'MIN/MAX OF VNEW ',vmin,vmax
        print *,'MIN/MAX OF WNEW ',wmin,wmax
        print *,' '
      end if

   end subroutine update_velocity_3d

end module update_module
