module mkscalforce_module

  ! this module contains the routines that make the force terms for
  ! the scalar quantities (rho and rho h).  There are two sets of 
  ! routines.  mkrhohforce makes the forcing terms that appear as
  ! source terms in the enthalpy equation (i.e. the heating).  
  ! modify_force adds to the forces the advective quantities that
  ! result from writing the equations in convective and perturbational 
  ! form for the edge-state prediction.

  use bl_constants_module
  use laplac_module
  use heating_module
  use eos_module

  implicit none

contains


  subroutine mkrhohforce_2d(force, s, ng, rho, ng_r, wmac, dx, bc, &
                            diff_coef, diff_fac, p0, rho0, temp0, time, pred_vs_corr)

    ! compute the source terms for the perturbational form of the
    ! enthalpy equation { w dp0/dr + (rho H - (rho0/sigma0) <sigma H>) }

    integer        , intent(in   ) :: ng,ng_r,pred_vs_corr
    real(kind=dp_t), intent(  out) :: force(0:,0:)
    real(kind=dp_t), intent(in   ) ::     s(1-ng  :,1-ng  :)
    real(kind=dp_t), intent(in   ) ::   rho(1-ng_r:,1-ng_r:)
    real(kind=dp_t), intent(in   ) ::  wmac(0     :,0     :)
    real(kind=dp_t), intent(in   ) ::    dx(:)
    integer        , intent(in   ) ::    bc(:,:) 
    real(kind=dp_t), intent(in   ) :: diff_coef, diff_fac
    real(kind=dp_t), intent(in   ) ::    p0(:), rho0(:), temp0(:)
    real(kind=dp_t), intent(in   ) :: time

    real(kind=dp_t) :: lapu
    real(kind=dp_t) :: gradp0,denom,Hbar,wadv,sigma0,coeff,sigma_H
    real(kind=dp_t), allocatable :: diff(:,:)
    real(kind=dp_t), allocatable ::    H(:,:)
    integer :: i,j,n,nx,ny,lo(2),hi(2)

    nx = size(force,dim=1) - 2
    ny = size(force,dim=2) - 2
    lo(:) = 1
    hi(1) = nx
    hi(2) = ny

    do_diag = .false.

    allocate(   H(lo(1):hi(1),lo(2):hi(2)))
    allocate(diff(lo(1):hi(1),lo(2):hi(2)))
 
    force = ZERO
    diff = ZERO

    denom = ONE / dble(hi(1)-lo(1)+1)

    if (diff_coef * diff_fac > 0.0_dp_t) then
       call laplac_2d(s(:,:),diff,dx,ng,bc(:,:))
       do j = 1,size(force,dim=2)-2
          do i = 1,size(force,dim=1)-2
             force(i,j) = diff_coef * diff_fac * diff(i,j)
          end do
       end do
    end if

    call get_H_2d(H,lo,hi,dx,time)

    if (pred_vs_corr .eq. 1) then

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             H(i,j) = rho(i,j) * H(i,j)
          end do
       end do

    else if (pred_vs_corr .eq. 2) then

       do j = lo(2),hi(2)
          den_row(1) = rho0(j)
          temp_row(1) = temp0(j)
          p_row(1) = p0(j)
          ! (rho,P) --> h, etc
          input_flag = 4
          call eos(input_flag, den_row, temp_row, npts, nspec, &
               xmass, aion, zion, &
               p_row, h_row, e_row, &
               cv_row, cp_row, xne_row, eta_row, &
               pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
               s_row, do_diag)
          sigma0 = dpdt_row(1) / (den_row(1) * cp_row(1) * dpdr_row(1))

          sigma_H = ZERO

          do i = lo(1), hi(1)
             den_row(1) = rho(i,j)
             temp_row(1) = temp0(j)
             p_row(1) = p0(j)
             ! (rho,P) --> h, etc
             input_flag = 4
             call eos(input_flag, den_row, temp_row, npts, nspec, &
                  xmass, aion, zion, &
                  p_row, h_row, e_row, &
                  cv_row, cp_row, xne_row, eta_row, &
                  pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
                  s_row, do_diag)
             coeff = dpdt_row(1) / (den_row(1) * cp_row(1) * dpdr_row(1))
             sigma_H = sigma_H + coeff * H(i,j)
          end do
          sigma_H = sigma_H * denom

          do i = lo(1),hi(1)
             H(i,j) = rho(i,j) * H(i,j) - (rho0(j) / sigma0) * sigma_H
          end do
       end do
       
    end if
 
!   Add w d(p0)/dz and full heating term to forcing term for (rho h)
    do j = 1,ny
       do i = 1,nx
          wadv = HALF*(wmac(i,j)+wmac(i,j+1))
          if (j.eq.1) then
             gradp0 =        (p0(j+1) - p0(j  )) / dx(2)
          else if (j.eq.ny) then
             gradp0 =        (p0(j  ) - p0(j-1)) / dx(2)
          else
             gradp0 = HALF * (p0(j+1) - p0(j-1)) / dx(2)
          end if
          force(i,j) = rho(i,j)*force(i,j) + wadv * gradp0 + H(i,j)
       end do
    end do

    deallocate(H)
    deallocate(diff)

  end subroutine mkrhohforce_2d

  subroutine mkrhohforce_3d(force, s, ng, rho, ng_r, wmac, dx, bc, &
                            diff_coef, diff_fac, p0, rho0, temp0, time, pred_vs_corr)

    ! compute the source terms for the perturbational form of the
    ! enthalpy equation { w dp0/dr + (rho H - (rho0/sigma0) <sigma H>) }

    integer        , intent(in   ) :: ng,ng_r,pred_vs_corr
    real(kind=dp_t), intent(  out) :: force(0:,0:,0:)
    real(kind=dp_t), intent(in   ) ::     s(1-ng  :,1-ng  :,1-ng:)
    real(kind=dp_t), intent(in   ) ::   rho(1-ng_r:,1-ng_r:,1-ng_r:)
    real(kind=dp_t), intent(in   ) ::  wmac(0     :,0     :,0:)
    real(kind=dp_t), intent(in   ) ::    dx(:)
    integer        , intent(in   ) ::    bc(:,:) 
    real(kind=dp_t), intent(in   ) :: diff_coef, diff_fac
    real(kind=dp_t), intent(in   ) ::    p0(:), rho0(:), temp0(:)
    real(kind=dp_t), intent(in   ) :: time

    real(kind=dp_t) :: lapu
    real(kind=dp_t) :: gradp0,denom,Hbar,wadv,sigma0,coeff,sigma_H
    real(kind=dp_t), allocatable ::    H(:,:,:)
    real(kind=dp_t), allocatable :: diff(:,:,:)
    integer :: i,j,k,n,nx,ny,nz,lo(3),hi(3)

    nx = size(force,dim=1) - 2
    ny = size(force,dim=2) - 2
    nz = size(force,dim=3) - 2
    lo(:) = 1
    hi(1) = nx
    hi(2) = ny
    hi(3) = nz

    do_diag = .false.

    allocate(   H(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    allocate(diff(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
 
    force = ZERO
    diff = ZERO

    denom = ONE / dble(hi(1)-lo(1)+1)

    if (diff_coef * diff_fac > 0.0_dp_t) then
       call laplac_3d(s(:,:,:),diff,dx,ng,bc(:,:))
       do k = 1,nz
          do j = 1,ny
          do i = 1,nx
             force(i,j,k) = diff_coef * diff_fac * diff(i,j,k)
          end do
          end do
       end do
    end if

    call get_H_3d(H,lo,hi,dx,time)

    if (pred_vs_corr .eq. 1) then

       do k = lo(3),hi(3)
       do j = lo(2),hi(2)
       do i = lo(1),hi(1)
          H(i,j,k) = rho(i,j,k) * H(i,j,k)
       end do
       end do
       end do

    else if (pred_vs_corr .eq. 2) then

       do k = lo(3),hi(3)
          den_row(1) = rho0(k)
          temp_row(1) = temp0(k)
          p_row(1) = p0(k)
          ! (rho,P) --> h, etc
          input_flag = 4
          call eos(input_flag, den_row, temp_row, npts, nspec, &
               xmass, aion, zion, &
               p_row, h_row, e_row, &
               cv_row, cp_row, xne_row, eta_row, &
               pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
               s_row, do_diag)
          sigma0 = dpdt_row(1) / (den_row(1) * cp_row(1) * dpdr_row(1))

          sigma_H = ZERO

          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             den_row(1) = rho(i,j,k)
             temp_row(1) = temp0(k)
             p_row(1) = p0(k)
             ! (rho,P) --> h, etc
             input_flag = 4
             call eos(input_flag, den_row, temp_row, npts, nspec, &
                  xmass, aion, zion, &
                  p_row, h_row, e_row, &
                  cv_row, cp_row, xne_row, eta_row, &
                  pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
                  s_row, do_diag)
             coeff = dpdt_row(1) / (den_row(1) * cp_row(1) * dpdr_row(1))
             sigma_H = sigma_H + coeff * H(i,j,k)
          end do
          end do
          sigma_H = sigma_H * denom

          do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             H(i,j,k) = rho(i,j,k) * H(i,j,k) - (rho0(k) / sigma0) * sigma_H
          end do
          end do
       end do
       
    end if
 
!   Add w d(p0)/dz and full heating term to forcing term for (rho h)
    do k = 1,nz
       do j = 1,ny
       do i = 1,nx
          wadv = HALF*(wmac(i,j,k)+wmac(i,j,k+1))
          if (k.eq.1) then
             gradp0 =        (p0(k+1) - p0(k  )) / dx(3)
          else if (k.eq.nz) then
          else
             if (wadv .gt. ZERO) then
                gradp0 = (p0(k) - p0(k-1)) / dx(2)
             else if (wadv .lt. ZERO) then
                gradp0 = (p0(k+1) - p0(k)) / dx(2)
             else 
                gradp0 = HALF * (p0(k+1) - p0(k-1)) / dx(2)
             end if
          end if
          force(i,j,k) = rho(i,j,k)*force(i,j,k) + wadv * gradp0 + H(i,j,k)
       end do
       end do
    end do

    deallocate(H)
    deallocate(diff)

  end subroutine mkrhohforce_3d

end module mkscalforce_module
