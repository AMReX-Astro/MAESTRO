module mkscalforce_module
  !
  ! this module contains the 2d and 3d routines that make the forcing
  ! terms for the scalar equations.
  ! 
  ! mkrhohforce computes the  w dp/dr  for rho*h evolution
  !
  ! mktempforce computes the source terms that appear in the
  ! temperature evolution equation.  Note, this formulation is only
  ! used in the prediction of the temperature on the edges, when
  ! enthalpy_pred_type = predict_T_then_rhohprime or predict_T_then_h.  
  ! The edge temperatures are then converted to enthalpy for the full update.
  !

  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: mkrhohforce, mktempforce

contains

  subroutine mkrhohforce(nlevs,scal_force,thermal,umac,p0_old,p0_new,psi,normal,dx, &
                         add_thermal,mla,the_bc_level)

    use bl_prof_module
    use variables, only: foextrap_comp, rhoh_comp
    use geometry, only: spherical
    use ml_restriction_module, only: ml_cc_restriction_c
    use multifab_fill_ghost_module
    use multifab_physbc_module

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: scal_force(:)
    type(multifab) , intent(in   ) :: thermal(:)
    type(multifab) , intent(in   ) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    logical        , intent(in   ) :: add_thermal
    type(ml_layout), intent(inout) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    integer                  :: i,dm,n
    integer                  :: lo(scal_force(1)%dim),hi(scal_force(1)%dim)    
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: tp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mkrhohforce")

    dm = scal_force(1)%dim
      
    do n=1,nlevs
       do i=1,scal_force(n)%nboxes
          if ( multifab_remote(scal_force(n),i) ) cycle
          fp => dataptr(scal_force(n), i)
          ump => dataptr(umac(n,1),i)
          vmp => dataptr(umac(n,2),i)
          tp  => dataptr(thermal(n),i)
          lo = lwb(get_box(scal_force(n),i))
          hi = upb(get_box(scal_force(n),i))
          select case (dm)
          case (2)
             call mkrhohforce_2d(n,fp(:,:,1,rhoh_comp), vmp(:,:,1,1), tp(:,:,1,1), lo, hi, &
                                 p0_old(n,:), p0_new(n,:), psi(n,:), add_thermal)
          case(3)
             wmp  => dataptr(umac(n,3), i)
             if (spherical .eq. 0) then
                call mkrhohforce_3d(n,fp(:,:,:,rhoh_comp), wmp(:,:,:,1), tp(:,:,:,1), &
                                    lo, hi, p0_old(n,:), p0_new(n,:), psi(n,:), add_thermal)
             else
                np => dataptr(normal(n), i)
                call mkrhohforce_3d_sphr(n,fp(:,:,:,rhoh_comp), &
                                         ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                         tp(:,:,:,1), lo, hi, dx(n,:), np(:,:,:,:), &
                                         p0_old(n,:), p0_new(n,:), psi(n,:), add_thermal)
             end if
          end select
       end do
    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(scal_force(nlevs),rhoh_comp,1)

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(scal_force(nlevs),rhoh_comp,foextrap_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(scal_force(n-1),rhoh_comp,scal_force(n),rhoh_comp, &
                                   mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(scal_force(n),scal_force(n-1), &
                                         scal_force(n)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         rhoh_comp,foextrap_comp,1)      
       end do

    end if

    call destroy(bpt)
    
  end subroutine mkrhohforce

  subroutine mkrhohforce_2d(n,rhoh_force,wmac,thermal,lo,hi,p0_old,p0_new,psi,add_thermal)

    use geometry, only: dr, nr
    use probin_module, only: enthalpy_pred_type
    use pred_parameters

    ! compute the source terms for the non-reactive part of the enthalpy equation {w dp0/dr}
    
    ! note, in the prediction of the interface states, we will set
    ! both p0_old and p0_new to the same old value.  In the computation
    ! of the rhoh_force for the update, they will be used to time-center.

    integer,         intent(in   ) :: n,lo(:),hi(:)
    real(kind=dp_t), intent(  out) :: rhoh_force(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) :: wmac(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) :: thermal(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:)
    real(kind=dp_t), intent(in   ) :: p0_new(0:)
    real(kind=dp_t), intent(in   ) :: psi(0:)
    logical        , intent(in   ) :: add_thermal

    real(kind=dp_t) :: gradp0, wadv
    integer :: i,j

!   Add w d(p0)/dz 
    do j = lo(2),hi(2)
       if (j.eq.0) then
          gradp0 = HALF * ( p0_old(j+1) + p0_new(j+1) &
                           -p0_old(j  ) - p0_new(j  ) ) / dr(n)
       else if (j.eq.nr(n)-1) then
          gradp0 = HALF * ( p0_old(j  ) + p0_new(j  ) &
                           -p0_old(j-1) - p0_new(j-1) ) / dr(n)
       else
          gradp0 = FOURTH * ( p0_old(j+1) + p0_new(j+1) &
                             -p0_old(j-1) - p0_new(j-1) ) / dr(n)
       end if
       do i = lo(1),hi(1)
          wadv = HALF*(wmac(i,j)+wmac(i,j+1))
          rhoh_force(i,j) =  wadv * gradp0           
       end do

       ! add psi as a force only if we are dealing with the full (rho h) or h,
       ! not the perturbational form.
       if ((enthalpy_pred_type .eq. predict_h) .or. &
           (enthalpy_pred_type .eq. predict_T_then_h)) then
          do i = lo(1),hi(1)
             rhoh_force(i,j) =  rhoh_force(i,j)  + psi(j)
          end do
       endif

    end do

    if (add_thermal) then
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             rhoh_force(i,j) = rhoh_force(i,j) + thermal(i,j)
          end do
       end do
    end if

  end subroutine mkrhohforce_2d

  subroutine mkrhohforce_3d(n,rhoh_force,wmac,thermal,lo,hi,p0_old,p0_new,psi,add_thermal)

    use geometry, only: dr, nr
    use probin_module, only: enthalpy_pred_type
    use pred_parameters

    ! compute the source terms for the non-reactive part of the enthalpy equation {w dp0/dr}

    integer,         intent(in   ) :: n,lo(:),hi(:)
    real(kind=dp_t), intent(  out) :: rhoh_force(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: wmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: thermal(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:)
    real(kind=dp_t), intent(in   ) :: p0_new(0:)
    real(kind=dp_t), intent(in   ) :: psi(0:)
    logical        , intent(in   ) :: add_thermal

    real(kind=dp_t) :: gradp0,wadv
    integer :: i,j,k

    do k = lo(3),hi(3)

       if (k.eq.0) then
          gradp0 = HALF * ( p0_old(k+1) + p0_new(k+1) &
                           -p0_old(k  ) - p0_new(k  ) ) / dr(n)
       else if (k.eq.nr(n)-1) then
          gradp0 = HALF * ( p0_old(k  ) + p0_new(k  ) &
                           -p0_old(k-1) - p0_new(k-1) ) / dr(n)
       else
          gradp0 = FOURTH * ( p0_old(k+1) + p0_new(k+1) &
                             -p0_old(k-1) - p0_new(k-1) ) / dr(n)
       end if

       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             wadv = HALF*(wmac(i,j,k)+wmac(i,j,k+1))
             rhoh_force(i,j,k) = wadv * gradp0 
          end do
       end do

       ! add psi as a force only if we are dealing with the full (rho h) or h,
       ! not the perturbational form.
       if ((enthalpy_pred_type .eq. predict_h) .or. &
           (enthalpy_pred_type .eq. predict_T_then_h)) then
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                rhoh_force(i,j,k) = rhoh_force(i,j,k) + psi(k)
             end do
          end do
       endif
       
    end do

    if (add_thermal) then
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhoh_force(i,j,k) = rhoh_force(i,j,k) + thermal(i,j,k)
             end do
          end do
       end do
    end if

  end subroutine mkrhohforce_3d

  subroutine mkrhohforce_3d_sphr(n,rhoh_force,umac,vmac,wmac,thermal,lo,hi,dx,normal, &
                                 p0_old,p0_new,psi,add_thermal)

    use fill_3d_module
    use geometry, only: nr, dr
    use probin_module, only: enthalpy_pred_type
    use pred_parameters

    ! compute the source terms for the non-reactive part of the enthalpy equation {w dp0/dr}

    integer,         intent(in   ) :: n,lo(:),hi(:)
    real(kind=dp_t), intent(  out) :: rhoh_force(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: umac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: vmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: wmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: thermal(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: normal(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:)
    real(kind=dp_t), intent(in   ) :: p0_new(0:)
    real(kind=dp_t), intent(in   ) :: psi(0:)
    logical        , intent(in   ) :: add_thermal

    real(kind=dp_t) :: uadv,vadv,wadv,normal_vel
    real(kind=dp_t), allocatable :: gradp_rad(:)
    real(kind=dp_t), allocatable :: gradp_cart(:,:,:,:)
    real(kind=dp_t), allocatable :: psi_cart(:,:,:,:)
    integer :: i,j,k,r

    allocate(gradp_rad(0:nr(n)-1))
    allocate(gradp_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    allocate(psi_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,1,psi,psi_cart, &
                                      lo,hi,dx,0)

    do r = 0, nr(n)-1
       
       if (r.eq.0) then
          gradp_rad(r) = HALF * ( p0_old(r+1) + p0_new(r+1) &
                                 -p0_old(r  ) - p0_new(r  ) ) / dr(n)
       else if (r.eq.nr(n)-1) then 
          gradp_rad(r) = HALF * ( p0_old(r  ) + p0_new(r  ) &
                                 -p0_old(r-1) - p0_new(r-1) ) / dr(n)
       else
          gradp_rad(r) = FOURTH * ( p0_old(r+1) + p0_new(r+1) &
                                   -p0_old(r-1) - p0_new(r-1) ) / dr(n)
       end if
       gradp_rad(r) = gradp_rad(r)
    end do

    call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,1,gradp_rad,gradp_cart, &
                                      lo,hi,dx,0)

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             uadv = HALF*(umac(i,j,k)+umac(i+1,j,k))
             vadv = HALF*(vmac(i,j,k)+vmac(i,j+1,k))
             wadv = HALF*(wmac(i,j,k)+wmac(i,j,k+1))

             normal_vel = uadv*normal(i,j,k,1)+vadv*normal(i,j,k,2)+wadv*normal(i,j,k,3)

             rhoh_force(i,j,k) = gradp_cart(i,j,k,1) * normal_vel

          end do
       end do
    end do

    ! add psi as a force only if we are dealing with the full (rho h) or h,
    ! not the perturbational form.
    if ((enthalpy_pred_type .eq. predict_h) .or. &
        (enthalpy_pred_type .eq. predict_T_then_h)) then 

       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                rhoh_force(i,j,k) = rhoh_force(i,j,k) + psi_cart(i,j,k,1)
             enddo
          enddo
       enddo
       
    endif

    if (add_thermal) then
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                rhoh_force(i,j,k) = rhoh_force(i,j,k) + thermal(i,j,k)
             end do
          end do
       end do
    end if

    deallocate(gradp_rad, gradp_cart, psi_cart)

  end subroutine mkrhohforce_3d_sphr



  subroutine mktempforce(nlevs,temp_force,umac,s,thermal,p0_old,p0_new,psi,normal, &
                         dx,mla,the_bc_level)

    use bl_prof_module
    use variables, only: foextrap_comp, temp_comp
    use geometry, only: spherical
    use ml_restriction_module, only: ml_cc_restriction_c
    use multifab_fill_ghost_module
    use multifab_physbc_module

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: temp_force(:)
    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: thermal(:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(ml_layout), intent(inout) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    integer                  :: i,dm,ng,n
    integer                  :: lo(temp_force(1)%dim),hi(temp_force(1)%dim)
    real(kind=dp_t), pointer :: tp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "mktempforce")

    dm = temp_force(1)%dim
    ng = s(1)%ng

    do n=1,nlevs

       do i=1,temp_force(n)%nboxes
          if ( multifab_remote(temp_force(n),i) ) cycle
          fp  => dataptr(temp_force(n),i)
          ump => dataptr(umac(n,1),i)
          vmp => dataptr(umac(n,2),i)
          sp  => dataptr(s(n),i)
          lo  =  lwb(get_box(s(n),i))
          hi  =  upb(get_box(s(n),i))
          tp  => dataptr(thermal(n),i)
          select case (dm)
          case (2)
             call mktempforce_2d(n, fp(:,:,1,temp_comp), sp(:,:,1,:), vmp(:,:,1,1), &
                                 tp(:,:,1,1), lo, hi, ng, p0_old(n,:), p0_new(n,:), psi(n,:))
          case(3)
             wmp => dataptr(umac(n,3),i)
             if (spherical .eq. 1) then
                np => dataptr(normal(n),i)
                call mktempforce_3d_sphr(n,fp(:,:,:,temp_comp), sp(:,:,:,:), &
                                         ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                         tp(:,:,:,1), lo, hi, ng, &
                                         p0_old(n,:), p0_new(n,:), psi(n,:), &
                                         np(:,:,:,:), dx(n,:))
             else
                call mktempforce_3d(n, fp(:,:,:,temp_comp), sp(:,:,:,:), wmp(:,:,:,1), &
                                    tp(:,:,:,1), lo, hi, ng, p0_old(n,:), p0_new(n,:), psi(n,:))
             end if
          end select
       end do
   
    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary_c(temp_force(nlevs),temp_comp,1)

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(temp_force(nlevs),temp_comp,foextrap_comp,1,the_bc_level(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction_c(temp_force(n-1),temp_comp,temp_force(n),temp_comp, &
                                   mla%mba%rr(n-1,:),1)

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(temp_force(n),temp_force(n-1), &
                                         temp_force(n)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_level(n-1),the_bc_level(n), &
                                         temp_comp,foextrap_comp,1)
       enddo

    end if

    call destroy(bpt)

  end subroutine mktempforce

  subroutine mktempforce_2d(n, temp_force, s, wmac, thermal, lo, hi, ng, p0_old, p0_new, psi)

    use geometry, only: dr, nr
    use variables, only: temp_comp, rho_comp, spec_comp
    use eos_module

    ! compute the source terms for temperature

    ! note, in the prediction of the interface states, we will set
    ! both p0_old and p0_new to the same old value.  In the computation
    ! of the temp_force for the update, they will be used to time-center.

    integer,         intent(in   ) :: lo(:),hi(:),ng
    integer,         intent(in   ) :: n
    real(kind=dp_t), intent(  out) :: temp_force(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t), intent(in   ) :: wmac(lo(1)-1:, lo(2)-1:)
    real(kind=dp_t), intent(in   ) :: thermal(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:), p0_new(0:), psi(0:)

    integer :: i,j

    real(kind=dp_t) :: gradp0, wadv, dhdp


    do j = lo(2),hi(2)

       if (j.eq.0) then
          gradp0 = HALF * ( p0_old(j+1) + p0_new(j+1) &
                           -p0_old(j  ) - p0_new(j  ) ) / dr(n)
       else if (j.eq.nr(n)-1) then
          gradp0 = HALF * ( p0_old(j  ) + p0_new(j  ) &
                           -p0_old(j-1) - p0_new(j-1) ) / dr(n)
       else
          gradp0 = FOURTH * ( p0_old(j+1) + p0_new(j+1) &
                             -p0_old(j-1) - p0_new(j-1) ) / dr(n)
       end if

       do i = lo(1),hi(1)

          temp_eos(1) = s(i,j,temp_comp)
          den_eos(1) = s(i,j,rho_comp)
          xn_eos(1,:) = s(i,j,spec_comp:spec_comp+nspec-1)

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

         dhdp = ONE / s(i,j,rho_comp) + ( s(i,j,rho_comp) * dedr_eos(1) -                  &
                                          p_eos(1) / s(i,j,rho_comp) )                     &
                                        / ( s(i,j,rho_comp) * dpdr_eos(1) )

         wadv = HALF*(wmac(i,j)+wmac(i,j+1))

         temp_force(i,j) =  thermal(i,j) + (ONE - s(i,j,rho_comp) * dhdp) * &
                            (wadv * gradp0 + psi(j))
         temp_force(i,j) = temp_force(i,j) / (cp_eos(1) * s(i,j,rho_comp))

       end do
    end do

  end subroutine mktempforce_2d

  subroutine mktempforce_3d(n, temp_force, s, wmac, thermal, lo, hi, ng, p0_old, p0_new, psi)

    use geometry,  only: dr, nr
    use variables, only: temp_comp, rho_comp, spec_comp
    use eos_module

    ! compute the source terms for temperature

    ! note, in the prediction of the interface states, we will set
    ! both p0_old and p0_new to the same old value.  In the computation
    ! of the temp_force for the update, they will be used to time-center.

    integer,         intent(in   ) :: lo(:),hi(:),ng, n
    real(kind=dp_t), intent(  out) :: temp_force(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: wmac(lo(1)-1:, lo(2)-1:, lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: thermal(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:), p0_new(0:), psi(0:)

    integer         :: i,j,k
    real(kind=dp_t) :: dhdp, gradp0, wadv


    do k = lo(3),hi(3)
       if (k.eq.0) then
          gradp0 = HALF * ( p0_old(k+1) + p0_new(k+1) &
                           -p0_old(k  ) - p0_new(k  ) ) / dr(n)
       else if (k.eq.nr(n)-1) then
          gradp0 = HALF * ( p0_old(k  ) + p0_new(k  ) &
                           -p0_old(k-1) - p0_new(k-1) ) / dr(n)
       else
          gradp0 = FOURTH * ( p0_old(k+1) + p0_new(k+1) &
                             -p0_old(k-1) - p0_new(k-1) ) / dr(n)
       end if

     do j = lo(2),hi(2)
       do i = lo(1),hi(1)

          temp_eos(1) = s(i,j,k,temp_comp)
          den_eos(1) = s(i,j,k,rho_comp)
          xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)
             
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

         dhdp = ONE / s(i,j,k,rho_comp) + ( s(i,j,k,rho_comp) * dedr_eos(1) -              &
                                            p_eos(1) / s(i,j,k,rho_comp) )                 &
                                          / ( s(i,j,k,rho_comp) * dpdr_eos(1) )

         wadv = HALF * (wmac(i,j,k+1) + wmac(i,j,k))

         temp_force(i,j,k) =  thermal(i,j,k) + &
                              (ONE - s(i,j,k,rho_comp) * dhdp) * &
                              (wadv * gradp0 + psi(k))

         temp_force(i,j,k) = temp_force(i,j,k) / (cp_eos(1) * s(i,j,k,rho_comp))

       end do
     end do
    end do

  end subroutine mktempforce_3d

  subroutine mktempforce_3d_sphr(n,temp_force, s, umac, vmac, wmac, thermal, &
                                 lo, hi, ng, p0_old, p0_new, psi, normal, dx)

    use fill_3d_module
    use variables, only: temp_comp, rho_comp, spec_comp
    use eos_module
    use geometry,  only: dr, nr

    ! compute the source terms for temperature

    ! note, in the prediction of the interface states, we will set
    ! both p0_old and p0_new to the same old value.  In the computation
    ! of the temp_force for the update, they will be used to time-center.

    integer,         intent(in   ) :: n,lo(:),hi(:),ng
    real(kind=dp_t), intent(  out) :: temp_force(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t), intent(in   ) :: umac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: vmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: wmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: thermal(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:), p0_new(0:), psi(0:)
    real(kind=dp_t), intent(in   ) :: normal(lo(1)-1:,lo(2)-1:,lo(3)-1:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    integer :: i,j,k,r
    real(kind=dp_t) :: uadv,vadv,wadv,normal_vel,dhdp
    real(kind=dp_t), allocatable :: gradp_rad(:)
    real(kind=dp_t), allocatable :: gradp_cart(:,:,:,:)
    real(kind=dp_t), allocatable :: psi_cart(:,:,:,:)

    allocate(gradp_rad(0:nr(n)-1))
    allocate(gradp_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    allocate(psi_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,1,psi,psi_cart, &
                                      lo,hi,dx,0)

    do r = 0, nr(n)-1
       
       if (r.eq.0) then
          gradp_rad(r) = HALF * ( p0_old(r+1) + p0_new(r+1) &
                                 -p0_old(r  ) - p0_new(r  ) ) / dr(n)
       else if (r.eq.nr(n)-1) then 
          gradp_rad(r) = HALF * ( p0_old(r  ) + p0_new(r  ) &
                                 -p0_old(r-1) - p0_new(r-1) ) / dr(n)
       else
          gradp_rad(r) = FOURTH * ( p0_old(r+1) + p0_new(r+1) &
                                   -p0_old(r-1) - p0_new(r-1) ) / dr(n)
       end if
    end do

    call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,1,gradp_rad,gradp_cart, &
                                      lo,hi,dx,0)

    do k = lo(3),hi(3)
     do j = lo(2),hi(2)
       do i = lo(1),hi(1)
        
          temp_eos(1) = s(i,j,k,temp_comp)
          den_eos(1) = s(i,j,k,rho_comp)
          xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)

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

         dhdp = ONE / s(i,j,k,rho_comp) + ( s(i,j,k,rho_comp) * dedr_eos(1) -              &
                                            p_eos(1) / s(i,j,k,rho_comp) )                 &
                                          / ( s(i,j,k,rho_comp) * dpdr_eos(1) )

         uadv = HALF * ( umac(i,j,k) + umac(i+1,j,k) )
         vadv = HALF * ( vmac(i,j,k) + vmac(i,j+1,k) )
         wadv = HALF * ( wmac(i,j,k) + wmac(i,j,k+1) )

         normal_vel = uadv*normal(i,j,k,1) + vadv*normal(i,j,k,2) + wadv*normal(i,j,k,3)

         temp_force(i,j,k) =  thermal(i,j,k) + &
                              (ONE - s(i,j,k,rho_comp) * dhdp) * &
                              (normal_vel * gradp_cart(i,j,k,1) + psi_cart(i,j,k,1))

         temp_force(i,j,k) = temp_force(i,j,k) / (cp_eos(1) * s(i,j,k,rho_comp))

       end do
     end do
    end do

    deallocate(gradp_rad, gradp_cart, psi_cart)

  end subroutine mktempforce_3d_sphr

end module mkscalforce_module
