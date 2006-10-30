module mkscalforce_module

  ! this module contains the 2d and 3d routines that make the 
  ! forcing term, w dp/dr,  for rho*h .  

  use bl_constants_module
  use variables

  implicit none

contains


  subroutine mkrhohforce_2d(force, wmac, p0_old, p0_new, dr)

    ! compute the source terms for the non-reactive part of the enthalpy equation { w dp0/dr }

    real(kind=dp_t), intent(  out) ::  force(0:,0:)
    real(kind=dp_t), intent(in   ) ::   wmac(0:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:)
    real(kind=dp_t), intent(in   ) :: p0_new(:)
    real(kind=dp_t), intent(in   ) ::     dr

    real(kind=dp_t) :: gradp0,wadv
    integer :: i,j,nx,ny

    nx = size(force,dim=1) - 2
    ny = size(force,dim=2) - 2
 
    force = ZERO
 
!   Add w d(p0)/dz and full heating term to forcing term for (rho h)
    do j = 1,ny
       do i = 1,nx
          wadv = HALF*(wmac(i,j)+wmac(i,j+1))
          if (j.eq.1) then
             gradp0 = HALF * ( p0_old(j+1) + p0_new(j+1) &
                              -p0_old(j  ) - p0_new(j  ) ) / dr
          else if (j.eq.ny) then
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

  subroutine mkrhohforce_3d(force, wmac, p0_old, p0_new, dr)

    ! compute the source terms for the non-reactive part of the enthalpy equation { w dp0/dr }

    real(kind=dp_t), intent(  out) ::  force(0:,0:,0:)
    real(kind=dp_t), intent(in   ) ::   wmac(0:,0:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:)
    real(kind=dp_t), intent(in   ) :: p0_new(:)
    real(kind=dp_t), intent(in   ) ::     dr

    real(kind=dp_t) :: gradp0,wadv
    integer :: i,j,k,nx,ny,nz

    nx = size(force,dim=1) - 2
    ny = size(force,dim=2) - 2
    nz = size(force,dim=3) - 2
 
    force = ZERO

    do k = 1,nz
       do j = 1,ny
       do i = 1,nx
          wadv = HALF*(wmac(i,j,k)+wmac(i,j,k+1))
          if (k.eq.1) then
             gradp0 = HALF * ( p0_old(k+1) + p0_new(k+1) &
                              -p0_old(k  ) - p0_new(k  ) ) / dr
          else if (k.eq.nz) then
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

end module mkscalforce_module
