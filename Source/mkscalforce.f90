module mkscalforce_module

  ! this module contains the 2d and 3d routines that make the 
  ! forcing term, w dp/dr,  for rho*h .  

  use eos_module
  use bl_constants_module
  use fill_3d_module
  use variables
  use geometry

  implicit none

contains


  subroutine mkrhohforce_2d(force, wmac, lo, hi, &
                            s_old, s_new, ng, dx, time, &
                            p0_old, p0_new, s0_old, s0_new, temp0, dr)

    ! compute the source terms for the non-reactive part of the enthalpy equation { w dp0/dr }
    
    ! note, in the prediction of the interface states, we will set
    ! both p0_old and p0_new to the same old value.  In the computation
    ! of the force for the update, they will be used to time-center.

    integer,         intent(in   ) :: lo(:), hi(:), ng
    real(kind=dp_t), intent(  out) ::  force(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(in   ) ::   wmac(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(in   ) ::  s_old(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(in   ) ::  s_new(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: time
    real(kind=dp_t), intent(in   ) :: p0_old(lo(2):)
    real(kind=dp_t), intent(in   ) :: p0_new(lo(2):)
    real(kind=dp_t), intent(in   ) :: s0_old(lo(2):,:)
    real(kind=dp_t), intent(in   ) :: s0_new(lo(2):,:)
    real(kind=dp_t), intent(in   ) :: temp0(lo(2):)
    real(kind=dp_t), intent(in   ) ::     dr

    real(kind=dp_t) :: gradp0, wadv, denom
    real(kind=dp_t) :: coeff_old, coeff_new, sigma_H, sigma0_old, sigma0_new
    integer :: i,j

    denom = ONE/dble(hi(1)-lo(1)+1)

    force = ZERO

!   Add w d(p0)/dz 
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          wadv = HALF*(wmac(i,j)+wmac(i,j+1))
          if (j.eq.lo(2)) then
             gradp0 = HALF * ( p0_old(j+1) + p0_new(j+1) &
                              -p0_old(j  ) - p0_new(j  ) ) / dr
          else if (j.eq.hi(2)) then
             gradp0 = HALF * ( p0_old(j  ) + p0_new(j  ) &
                              -p0_old(j-1) - p0_new(j-1) ) / dr
          else
             gradp0 = FOURTH * ( p0_old(j+1) + p0_new(j+1) &
                                -p0_old(j-1) - p0_new(j-1) ) / dr
          end if
          force(i,j) =  wadv * gradp0 
       end do
    end do

  end subroutine mkrhohforce_2d

  subroutine mkrhohforce_3d(force, wmac, lo, hi, &
                            s_old, s_new, ng, dx, time, &
                            p0_old, p0_new, s0_old, s0_new, temp0, dr)

    ! compute the source terms for the non-reactive part of the enthalpy equation { w dp0/dr }

    integer,         intent(in   ) :: lo(:), hi(:), ng
    real(kind=dp_t), intent(  out) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real(kind=dp_t), intent(in   ) ::   wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real(kind=dp_t), intent(in   ) ::  s_old(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) ::  s_new(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: time
    real(kind=dp_t), intent(in   ) :: p0_old(lo(2):)
    real(kind=dp_t), intent(in   ) :: p0_new(lo(2):)
    real(kind=dp_t), intent(in   ) :: s0_old(lo(2):,:)
    real(kind=dp_t), intent(in   ) :: s0_new(lo(2):,:)
    real(kind=dp_t), intent(in   ) :: temp0(lo(2):)
    real(kind=dp_t), intent(in   ) ::     dr

    real(kind=dp_t) :: gradp0,wadv
    integer :: i,j,k
 
    force = ZERO

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             wadv = HALF*(wmac(i,j,k)+wmac(i,j,k+1))

             if (k.eq.lo(3)) then
                gradp0 = HALF * ( p0_old(k+1) + p0_new(k+1) &
                                 -p0_old(k  ) - p0_new(k  ) ) / dr
             else if (k.eq.hi(3)) then
                gradp0 = HALF * ( p0_old(k  ) + p0_new(k  ) &
                                 -p0_old(k-1) - p0_new(k-1) ) / dr
             else
                gradp0 = FOURTH * ( p0_old(k+1) + p0_new(k+1) &
                                   -p0_old(k-1) - p0_new(k-1) ) / dr
             end if
             force(i,j,k) = wadv * gradp0 
          end do
       end do
    end do

  end subroutine mkrhohforce_3d

  subroutine mkrhohforce_3d_sphr(force, umac, vmac, wmac, lo, hi, &
                                 s_old, s_new, ng, dx, time, normal, &
                                 p0_old, p0_new, s0_old, s0_new, temp0)

    ! compute the source terms for the non-reactive part of the enthalpy equation { w dp0/dr }

    integer,         intent(in   ) :: lo(:), hi(:), ng
    real(kind=dp_t), intent(  out) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)- 1:,lo(2)- 1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::   vmac(lo(1)- 1:,lo(2)- 1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::   wmac(lo(1)- 1:,lo(2)- 1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: normal(lo(1)- 1:,lo(2)- 1:,lo(3)-1:,:)
    real(kind=dp_t), intent(in   ) ::  s_old(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) ::  s_new(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: time
    real(kind=dp_t), intent(in   ) :: p0_old(0:)
    real(kind=dp_t), intent(in   ) :: p0_new(0:)
    real(kind=dp_t), intent(in   ) :: s0_old(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(0:,:)
    real(kind=dp_t), intent(in   ) :: temp0(0:)

    real(kind=dp_t) :: uadv,vadv,wadv,normal_vel
    real(kind=dp_t), allocatable :: gradp_rad(:)
    real(kind=dp_t), allocatable :: gradp_cart(:,:,:)
    integer :: i,j,k,nr

    nr = size(p0_old,dim=1)

    allocate(gradp_rad(0:nr-1))
    allocate(gradp_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
 
    force = ZERO

    do k = 0, nr-1
       
       if (k.eq.0) then
          gradp_rad(k) = HALF * ( p0_old(k+1) + p0_new(k+1) &
                                 -p0_old(k  ) - p0_new(k  ) ) / dr
       else if (k.eq.nr-1) then 
          gradp_rad(k) = HALF * ( p0_old(k  ) + p0_new(k  ) &
                                 -p0_old(k-1) - p0_new(k-1) ) / dr
       else
          gradp_rad(k) = FOURTH * ( p0_old(k+1) + p0_new(k+1) &
                                   -p0_old(k-1) - p0_new(k-1) ) / dr
       end if
    end do

    call fill_3d_data(gradp_cart,gradp_rad,lo,hi,dx,0)

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             uadv = HALF*(umac(i,j,k)+umac(i+1,j,k))
             vadv = HALF*(vmac(i,j,k)+vmac(i,j+1,k))
             wadv = HALF*(wmac(i,j,k)+wmac(i,j,k+1))
             normal_vel = uadv*normal(i,j,k,1)+vadv*normal(i,j,k,2)+wadv*normal(i,j,k,3)
             force(i,j,k) = gradp_cart(i,j,k) * normal_vel
          end do
       end do
    end do

    deallocate(gradp_rad, gradp_cart)

  end subroutine mkrhohforce_3d_sphr

end module mkscalforce_module
