module update_scal_module

  use bl_types
  use multifab_module
  use bl_constants_module
  use fill_3d_module
  use addw0_module
  use geometry
  use sponge
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

   subroutine update_scal_3d_cart (nstart,nstop,sold,snew,umac,vmac,wmac,w0,w0_cart,sedgex,sedgey,sedgez,&
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

      ! not spherical

        allocate(delta_base(lo(3):hi(3)))
        do n = nstart, nstop

          if (lo(3) .eq. 0) then
             base_edge(lo(3)  ) = base_old(lo(3),n)
             base_edge(lo(3)+1) = HALF*(base_old(lo(3),n)+base_old(lo(3)+1,n))
          else
             k = lo(3)
             base_edge(k) = 7.d0/12.d0 * (base_old(k  ,n) + base_old(k-1,n)) &
                           -1.d0/12.d0 * (base_old(k+1,n) + base_old(k-2,n))
             k = lo(3)+1
             base_edge(k) = 7.d0/12.d0 * (base_old(k  ,n) + base_old(k-1,n)) &
                           -1.d0/12.d0 * (base_old(k+1,n) + base_old(k-2,n))
          end if

          if (hi(3) .eq. nr) then
             base_edge(hi(3)+1) = base_old(hi(3),n)
             base_edge(hi(3)  ) = HALF*(base_old(hi(3),n)+base_old(hi(3)-1,n))
          else
             k = hi(3)+1
             base_edge(k) = 7.d0/12.d0 * (base_old(k  ,n) + base_old(k-1,n)) &
                           -1.d0/12.d0 * (base_old(k+1,n) + base_old(k-2,n))
             k = hi(3)
             base_edge(k) = 7.d0/12.d0 * (base_old(k  ,n) + base_old(k-1,n)) &
                           -1.d0/12.d0 * (base_old(k+1,n) + base_old(k-2,n))
          end if
          
    
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

   end subroutine update_scal_3d_cart

   subroutine update_scal_3d_sphr (nstart,nstop,sold,snew,umac,vmac,wmac,w0,w0_cart,sedgex,sedgey,sedgez,&
                                   force,base_old,base_new,base_cart,lo,hi,ng,dx,dt)

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
      real (kind = dp_t), intent(in   ) ::   base_cart(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
      real (kind = dp_t), intent(in   ) :: w0(0:)
      real (kind = dp_t), intent(in   ) :: w0_cart(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:,:)
      real (kind = dp_t), intent(in   ) :: dt,dx(:)

      integer :: i, j, k, n, nr
      real (kind = dp_t) :: divsu,divbaseu,mult
      real (kind = dp_t), allocatable :: delta_base(:),delta_base_cart(:,:,:)

      real (kind = dp_t) :: divu

      nr = size(base_old,dim=1)

      ! is spherical

      allocate(delta_base(0:nr-1))
      allocate(delta_base_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

      do n = nstart, nstop
          do k = 0,nr-1
            delta_base(k) = base_new(k,n) - base_old(k,n)
          end do
          call fill_3d_data(delta_base_cart,delta_base,lo,hi,dx,0)

          ! Note the umac here does NOT have w0 in it
          do k = lo(3), hi(3)
            do j = lo(2), hi(2)
            do i = lo(1), hi(1)
    
              divbaseu = HALF * (  &
                   (umac(i+1,j,k)*(base_cart(i,j,k,n)+base_cart(i+1,j,k,n)) &
                   -umac(i  ,j,k)*(base_cart(i,j,k,n)+base_cart(i-1,j,k,n)) ) / dx(1) &
                  +(vmac(i,j+1,k)*(base_cart(i,j,k,n)+base_cart(i,j+1,k,n)) &
                   -vmac(i,j  ,k)*(base_cart(i,j,k,n)+base_cart(i,j-1,k,n)) ) / dx(2) &
                  +(wmac(i,j,k+1)*(base_cart(i,j,k,n)+base_cart(i,j,k+1,n)) &
                   -wmac(i,j,k  )*(base_cart(i,j,k,n)+base_cart(i,j,k-1,n)) ) / dx(3)  )

              divu = (umac(i+1,j,k) - umac(i  ,j,k) ) / dx(1) &
                    +(vmac(i,j+1,k) - vmac(i,j  ,k) ) / dx(2) &
                    +(wmac(i,j,k+1) - wmac(i,j,k  ) ) / dx(3) 

              snew(i,j,k,n) = sold(i,j,k,n) + delta_base_cart(i,j,k) &
                              - dt * divbaseu + dt * force(i,j,k,n)

            enddo
            enddo
          enddo
        end do

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

      deallocate(delta_base,delta_base_cart)

      mult = -ONE
      call addw0_3d_sphr(umac,vmac,wmac,w0_cart,lo,hi,dx,mult)

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

   end subroutine update_scal_3d_sphr

end module update_scal_module
