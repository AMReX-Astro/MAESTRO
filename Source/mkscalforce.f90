module mkscalforce_module

  ! this module contains the 2d and 3d routines that make the 
  ! forcing term, w dp/dr,  for rho*h .  

  use eos_module
  use bl_constants_module
  use fill_3d_module
  use variables
  use geometry
  use heating_module

  implicit none

contains


  subroutine mkrhohforce_2d(force, wmac, lo, hi, &
                            s, ng, dx, time, &
                            p0_old, p0_new, s0_old, s0_new, temp0, dr)

    ! compute the source terms for the non-reactive part of the enthalpy equation { w dp0/dr }

    real(kind=dp_t), intent(  out) ::  force(0:,0:)
    real(kind=dp_t), intent(in   ) ::   wmac(0:,0:)
    integer,         intent(in   ) :: lo(:), hi(:), ng
    real(kind=dp_t), intent(in   ) ::      s(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: time
    real(kind=dp_t), intent(in   ) :: p0_old(:)
    real(kind=dp_t), intent(in   ) :: p0_new(:)
    real(kind=dp_t), intent(in   ) :: s0_old(:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(:,:)
    real(kind=dp_t), intent(in   ) :: temp0(:)
    real(kind=dp_t), intent(in   ) ::     dr

    real(kind=dp_t), allocatable :: H(:,:)
    real(kind=dp_t) :: gradp0,wadv,denom,coeff,sigma_H, sigma0
    integer :: i,j,nx,ny

    nx = size(force,dim=1) - 2
    ny = size(force,dim=2) - 2
 
    allocate(H(lo(1):hi(1),lo(2):hi(2)))

    denom = ONE/dble(hi(1)-lo(1)+1)

    force = ZERO

    call get_H_2d(H,lo,hi,dx,time)

    do j = lo(2),hi(2)

       den_row(1) = HALF*(s0_old(j,rho_comp) + s0_new(j,rho_comp))
       temp_row(1) = temp0(j)
       p_row(1) = HALF*(p0_old(j) + p0_new(j))

       xn_zone(:) = HALF*(s0_old(j,spec_comp:spec_comp-1+nspec)/s0_old(j,rho_comp) + &
                          s0_new(j,spec_comp:spec_comp-1+nspec)/s0_new(j,rho_comp))

       ! (rho,P) --> h, etc
       input_flag = 4

       call eos(input_flag, den_row, temp_row, &
                npts, nspec, &
                xn_zone, aion, zion, &
                p_row, h_row, e_row, &
                cv_row, cp_row, xne_row, eta_row, pele_row, &
                dpdt_row, dpdr_row, dedt_row, dedr_row, &
                dpdX_row, dhdX_row, &
                gam1_row, cs_row, s_row, &
                do_diag)

       sigma0 = dpdt_row(1) / (den_row(1) * cp_row(1) * dpdr_row(1))

       sigma_H = ZERO

       do i = lo(1), hi(1)
          den_row(1) = s(i,j,rho_comp)
          temp_row(1) = temp0(j)
          p_row(1) = HALF*(p0_old(j) + p0_new(j))

          xn_zone(:) = s(i,j,spec_comp:spec_comp-1+nspec)/s(i,j,rho_comp) 

          ! (rho,P) --> h, etc
          input_flag = 4

          call eos(input_flag, den_row, temp_row, &
                   npts, nspec, &
                   xn_zone, aion, zion, &
                   p_row, h_row, e_row, &
                   cv_row, cp_row, xne_row, eta_row, pele_row, &
                   dpdt_row, dpdr_row, dedt_row, dedr_row, &
                   dpdX_row, dhdX_row, &
                   gam1_row, cs_row, s_row, &
                   do_diag)

          coeff = dpdt_row(1) / (den_row(1) * cp_row(1) * dpdr_row(1))
          sigma_H = sigma_H + coeff * H(i,j)
       end do

       sigma_H = sigma_H * denom

       do i = lo(1),hi(1)
          H(i,j) = s(i,j,rho_comp) * H(i,j) - (HALF*(s0_old(j,rho_comp)+s0_new(j,rho_comp)) / sigma0) * sigma_H
       end do
    end do

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
          force(i,j) =  wadv * gradp0 + H(i,j)
       end do
    end do

  end subroutine mkrhohforce_2d

  subroutine mkrhohforce_3d(force, wmac, lo, hi, &
                            s, ng, dx, time, &
                            p0_old, p0_new, s0_old, s0_new, temp0, dr)

    ! compute the source terms for the non-reactive part of the enthalpy equation { w dp0/dr }

    real(kind=dp_t), intent(  out) ::  force(0:,0:,0:)
    real(kind=dp_t), intent(in   ) ::   wmac(0:,0:,0:)
    integer,         intent(in   ) :: lo(:), hi(:), ng
    real(kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: time
    real(kind=dp_t), intent(in   ) :: p0_old(:)
    real(kind=dp_t), intent(in   ) :: p0_new(:)
    real(kind=dp_t), intent(in   ) :: s0_old(:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(:,:)
    real(kind=dp_t), intent(in   ) :: temp0(:)
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

  subroutine mkrhohforce_3d_sphr(force, umac, vmac, wmac, lo, hi, &
                                 s, ng, dx, time, normal, &
                                 p0_old, p0_new, s0_old, s0_new, temp0)

    ! compute the source terms for the non-reactive part of the enthalpy equation { w dp0/dr }

    real(kind=dp_t), intent(  out) ::  force(0:,0:,0:)
    real(kind=dp_t), intent(in   ) ::   umac(0:,0:,0:)
    real(kind=dp_t), intent(in   ) ::   vmac(0:,0:,0:)
    real(kind=dp_t), intent(in   ) ::   wmac(0:,0:,0:)
    real(kind=dp_t), intent(in   ) :: normal(0:,0:,0:,:)
    integer,         intent(in   ) :: lo(:), hi(:), ng
    real(kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: time
    real(kind=dp_t), intent(in   ) :: p0_old(:)
    real(kind=dp_t), intent(in   ) :: p0_new(:)
    real(kind=dp_t), intent(in   ) :: s0_old(:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(:,:)
    real(kind=dp_t), intent(in   ) :: temp0(:)

    real(kind=dp_t) :: uadv,vadv,wadv,normal_vel
    real(kind=dp_t), allocatable :: gradp_rad(:)
    real(kind=dp_t), allocatable :: gradp_cart(:,:,:)
    integer :: i,j,k,nx,ny,nz,nr

    nx = size(force,dim=1) - 2
    ny = size(force,dim=2) - 2
    nz = size(force,dim=3) - 2

    nr = size(p0_old,dim=1)

    allocate(gradp_rad(nr))
    allocate(gradp_cart(nx,ny,nz))
 
    force = ZERO

    do k = 1,nr
      if (k.eq.1) then
         gradp_rad(k) = HALF * ( p0_old(k+1) + p0_new(k+1) &
                                -p0_old(k  ) - p0_new(k  ) ) / dr
      else if (k.eq.nr) then
         gradp_rad(k) = HALF * ( p0_old(k  ) + p0_new(k  ) &
                                -p0_old(k-1) - p0_new(k-1) ) / dr
      else
         gradp_rad(k) = FOURTH * ( p0_old(k+1) + p0_new(k+1) &
                                  -p0_old(k-1) - p0_new(k-1) ) / dr
      end if
    end do

    call fill_3d_data(gradp_cart,gradp_rad,dx,0)

    do k = 1,nz
       do j = 1,ny
       do i = 1,nx
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
