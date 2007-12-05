module mkscalforce_module

  ! this module contains the 2d and 3d routines that make the 
  ! forcing term, w dp/dr,  for rho*h .  

  use bl_constants_module
  use fill_3d_module
  use variables
  use geometry
  use eos_module
  use multifab_module

  implicit none

  private
  public :: mkrhohforce
  public :: mktempforce
contains

  subroutine mkrhohforce(n,force,comp,umac,p0_old,p0_new,normal,dx)

    integer        , intent(in   ) :: n
    type(multifab) , intent(inout) :: force
    integer        , intent(in   ) :: comp
    type(multifab) , intent(in   ) :: umac(:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:)
    real(kind=dp_t), intent(in   ) :: p0_new(0:)
    type(multifab) , intent(in   ) :: normal
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    integer                  :: i,dm
    integer                  :: lo(force%dim),hi(force%dim)    
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)

    dm = force%dim
      
    do i=1,force%nboxes
       if ( multifab_remote(force,i) ) cycle
       fp => dataptr(force, i)
       ump => dataptr(umac(1),i)
       vmp => dataptr(umac(2),i)
       lo = lwb(get_box(force,i))
       hi = upb(get_box(force,i))
       select case (dm)
       case (2)
          call  mkrhohforce_2d(n,fp(:,:,1,comp), vmp(:,:,1,1), lo, hi, p0_old, p0_new)
       case(3)
          wmp  => dataptr(umac(3), i)
          if (spherical .eq. 0) then
             call mkrhohforce_3d(n,fp(:,:,:,comp), wmp(:,:,:,1), lo, hi, p0_old, p0_new)
          else
             np => dataptr(normal, i)
             call mkrhohforce_3d_sphr(n,fp(:,:,:,comp), &
                                      ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                      lo, hi, dx, np(:,:,:,:), p0_old, p0_new)
          end if
       end select
    end do
    
  end subroutine mkrhohforce

  subroutine mkrhohforce_2d(n,force,wmac,lo,hi,p0_old,p0_new)

    ! compute the source terms for the non-reactive part of the enthalpy equation {w dp0/dr}
    
    ! note, in the prediction of the interface states, we will set
    ! both p0_old and p0_new to the same old value.  In the computation
    ! of the force for the update, they will be used to time-center.

    integer,         intent(in   ) :: n,lo(:),hi(:)
    real(kind=dp_t), intent(  out) ::  force(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(in   ) ::   wmac(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:)
    real(kind=dp_t), intent(in   ) :: p0_new(0:)

    real(kind=dp_t) :: gradp0, wadv
    integer :: i,j

    force = ZERO

!   Add w d(p0)/dz 
    do j = lo(2),hi(2)
       if (j.eq.0) then
          gradp0 = HALF * ( p0_old(j+1) + p0_new(j+1) &
                           -p0_old(j  ) - p0_new(j  ) ) / dr(1)
       else if (j.eq.nr(n)-1) then
          gradp0 = HALF * ( p0_old(j  ) + p0_new(j  ) &
                           -p0_old(j-1) - p0_new(j-1) ) / dr(1)
       else
          gradp0 = FOURTH * ( p0_old(j+1) + p0_new(j+1) &
                             -p0_old(j-1) - p0_new(j-1) ) / dr(1)
       end if
       do i = lo(1),hi(1)
          wadv = HALF*(wmac(i,j)+wmac(i,j+1))
          force(i,j) =  wadv * gradp0 
       end do
    end do

  end subroutine mkrhohforce_2d

  subroutine mkrhohforce_3d(n,force,wmac,lo,hi,p0_old,p0_new)

    ! compute the source terms for the non-reactive part of the enthalpy equation {w dp0/dr}

    integer,         intent(in   ) :: n,lo(:),hi(:)
    real(kind=dp_t), intent(  out) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real(kind=dp_t), intent(in   ) ::   wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:)
    real(kind=dp_t), intent(in   ) :: p0_new(0:)

    real(kind=dp_t) :: gradp0,wadv
    integer :: i,j,k

    force = ZERO
 
    do k = lo(3),hi(3)

       if (k.eq.0) then
          gradp0 = HALF * ( p0_old(k+1) + p0_new(k+1) &
                           -p0_old(k  ) - p0_new(k  ) ) / dr(1)
       else if (k.eq.nr(n)-1) then
          gradp0 = HALF * ( p0_old(k  ) + p0_new(k  ) &
                           -p0_old(k-1) - p0_new(k-1) ) / dr(1)
       else
          gradp0 = FOURTH * ( p0_old(k+1) + p0_new(k+1) &
                             -p0_old(k-1) - p0_new(k-1) ) / dr(1)
       end if

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             wadv = HALF*(wmac(i,j,k)+wmac(i,j,k+1))
             force(i,j,k) = wadv * gradp0 
          end do
       end do

    end do

  end subroutine mkrhohforce_3d

  subroutine mkrhohforce_3d_sphr(n,force,umac,vmac,wmac,lo,hi,dx,normal,p0_old,p0_new)

    ! compute the source terms for the non-reactive part of the enthalpy equation {w dp0/dr}

    integer,         intent(in   ) :: n,lo(:),hi(:)
    real(kind=dp_t), intent(  out) ::  force(lo(1)- 1:,lo(2)- 1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::   umac(lo(1)- 1:,lo(2)- 1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::   vmac(lo(1)- 1:,lo(2)- 1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::   wmac(lo(1)- 1:,lo(2)- 1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: normal(lo(1)- 1:,lo(2)- 1:,lo(3)-1:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:)
    real(kind=dp_t), intent(in   ) :: p0_new(0:)

    real(kind=dp_t) :: uadv,vadv,wadv,normal_vel
    real(kind=dp_t), allocatable :: gradp_rad(:)
    real(kind=dp_t), allocatable :: gradp_cart(:,:,:)
    integer :: i,j,k

    allocate(gradp_rad(0:nr(n)-1))
    allocate(gradp_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
 
    force = ZERO

    do k = 0, nr(n)-1
       
       if (k.eq.0) then
          gradp_rad(k) = HALF * ( p0_old(k+1) + p0_new(k+1) &
                                 -p0_old(k  ) - p0_new(k  ) ) / dr(1)
       else if (k.eq.nr(n)-1) then 
          gradp_rad(k) = HALF * ( p0_old(k  ) + p0_new(k  ) &
                                 -p0_old(k-1) - p0_new(k-1) ) / dr(1)
       else
          gradp_rad(k) = FOURTH * ( p0_old(k+1) + p0_new(k+1) &
                                   -p0_old(k-1) - p0_new(k-1) ) / dr(1)
       end if
    end do

    call fill_3d_data(n,gradp_cart,gradp_rad,lo,hi,dx,0)

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

  subroutine mktempforce(n,force,comp,s,thermal,p0_old,dx)

    integer        , intent(in   ) :: n
    type(multifab) , intent(inout) :: force
    integer        , intent(in   ) :: comp
    type(multifab) , intent(in   ) :: s
    type(multifab) , intent(in   ) :: thermal
    real(kind=dp_t), intent(in   ) :: p0_old(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    integer                  :: i,dm,ng
    integer                  :: lo(force%dim),hi(force%dim)
    real(kind=dp_t), pointer :: tp(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)

    dm = force%dim
    ng = s%ng

    do i=1,force%nboxes
       if ( multifab_remote(force,i) ) cycle
       fp => dataptr(force,i)
       sp => dataptr(s, i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       select case (dm)
       case (2)
          tp => dataptr(thermal,i)
          call mktempforce_2d(fp(:,:,1,comp), sp(:,:,1,:), tp(:,:,1,1), lo, hi, &
                              ng, p0_old)
       case(3)
          if (spherical .eq. 1) then
             tp => dataptr(thermal, i)
             call mktempforce_3d_sphr(n,fp(:,:,:,comp), sp(:,:,:,:), tp(:,:,:,1), &
                                      lo, hi, ng, p0_old, dx)
          else
             tp => dataptr(thermal, i)
             call mktempforce_3d(fp(:,:,:,comp), sp(:,:,:,:), tp(:,:,:,1), lo, hi, &
                                 ng, p0_old)
          end if
       end select
    end do
    
  end subroutine mktempforce

  subroutine mktempforce_2d(force, s, thermal, lo, hi, ng, p0)

    ! compute the source terms for temperature

    ! note, in the prediction of the interface states, we will set
    ! both p0_old and p0_new to the same old value.  In the computation
    ! of the force for the update, they will be used to time-center.

    integer,         intent(in   ) :: lo(:), hi(:), ng
    real(kind=dp_t), intent(  out) ::   force(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(in   ) :: thermal(lo(1)-1 :,lo(2)-1 :)
    real(kind=dp_t), intent(in   ) ::      p0(0:)

    integer :: i,j

    force = ZERO

!   HACK HACK HACK 
!   We ignore the w d(p0) / dz term since p0 is essentially constant

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)

          temp_eos(1) = s(i,j,temp_comp)
          den_eos(1) = s(i,j,rho_comp)
          xn_eos(1,:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)

          ! dens, temp, xmass inputs
         call eos(eos_input_rt, den_eos, temp_eos, &
                  npts, nspec, &
                  xn_eos, &
                  p_eos, h_eos, e_eos, &
                  cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                  dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                  dpdX_eos, dhdX_eos, &
                  gam1_eos, cs_eos, s_eos, &
                  dsdt_eos, dsdr_eos, &
                  do_diag)

          force(i,j) =  thermal(i,j) / (s(i,j,rho_comp) * cp_eos(1))

       end do
    end do

  end subroutine mktempforce_2d

  subroutine mktempforce_3d(force, s, thermal, lo, hi, ng, p0)

    ! compute the source terms for temperature

    ! note, in the prediction of the interface states, we will set
    ! both p0_old and p0_new to the same old value.  In the computation
    ! of the force for the update, they will be used to time-center.

    integer,         intent(in   ) :: lo(:), hi(:), ng
    real(kind=dp_t), intent(  out) ::   force(lo(1)-1 :,lo(2)-1 :,lo(3)-1 :)
    real(kind=dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: thermal(lo(1)-1 :,lo(2)-1 :,lo(3)-1 :)
    real(kind=dp_t), intent(in   ) ::      p0(0:)

    integer :: i,j,k

    force = ZERO

!   HACK HACK HACK 
!   We ignore the w d(p0) / dz term since p0 is essentially constant

    do k = lo(3),hi(3)
     do j = lo(2),hi(2)
       do i = lo(1),hi(1)

          temp_eos(1) = s(i,j,k,temp_comp)
          den_eos(1) = s(i,j,k,rho_comp)
          xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)
          
          ! dens, temp, xmass inputs
         call eos(eos_input_rt, den_eos, temp_eos, &
                  npts, nspec, &
                  xn_eos, &
                  p_eos, h_eos, e_eos, &
                  cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                  dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                  dpdX_eos, dhdX_eos, &
                  gam1_eos, cs_eos, s_eos, &
                  dsdt_eos, dsdr_eos, &
                  do_diag)

          force(i,j,k) =  thermal(i,j,k) / (s(i,j,k,rho_comp) * cp_eos(1))

       end do
     end do
    end do

  end subroutine mktempforce_3d

  subroutine mktempforce_3d_sphr(n,force, s, thermal, lo, hi, ng, p0, dx)

    ! compute the source terms for temperature

    ! note, in the prediction of the interface states, we will set
    ! both p0_old and p0_new to the same old value.  In the computation
    ! of the force for the update, they will be used to time-center.

    integer,         intent(in   ) :: n,lo(:), hi(:), ng
    real(kind=dp_t), intent(  out) ::   force(lo(1)-1 :,lo(2)-1 :,lo(3)-1 :)
    real(kind=dp_t), intent(in   ) ::       s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: thermal(lo(1)-1 :,lo(2)-1 :,lo(3)-1 :)
    real(kind=dp_t), intent(in   ) ::      p0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), allocatable   :: p0_cart(:,:,:)

    integer :: i,j,k

    allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(n,p0_cart,p0,lo,hi,dx,0)

    force = ZERO

!   HACK HACK HACK 
!   We ignore the w d(p0) / dz term since p0 is essentially constant

    do k = lo(3),hi(3)
     do j = lo(2),hi(2)
       do i = lo(1),hi(1)

          temp_eos(1) = s(i,j,k,temp_comp)
          den_eos(1) = s(i,j,k,rho_comp)
          xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)
          
          ! dens, temp, xmass inputs
         call eos(eos_input_rt, den_eos, temp_eos, &
                  npts, nspec, &
                  xn_eos, &
                  p_eos, h_eos, e_eos, &
                  cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                  dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                  dpdX_eos, dhdX_eos, &
                  gam1_eos, cs_eos, s_eos, &
                  dsdt_eos, dsdr_eos, &
                  do_diag)

          force(i,j,k) =  thermal(i,j,k) / (s(i,j,k,rho_comp) * cp_eos(1))

       end do
     end do
    end do

    deallocate(p0_cart)

  end subroutine mktempforce_3d_sphr

end module mkscalforce_module
