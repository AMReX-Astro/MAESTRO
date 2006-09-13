module scalar_advance_module

  ! IF pred_vs_corr = 1:
  ! This routine does the advancement of rho and (rho h), treating the
  ! base state as if it were constant.  This is step 2 in Almgren et
  ! al. 2006 (paper II).

  ! IF pred_vs_corr = 2:
  ! This routine does the advancement of rho and (rho h), after the
  ! base state adjustmented.  This is step 4 in Almgren et
  ! al. 2006 (paper II).

  use bl_types
  use bl_constants_module
  use multifab_module
  use viscous_module
  use mkflux_module
  use mkscalforce_module
  use update_module
  use addw0_module
  use setbc_module
  use variables

  implicit none

contains

   subroutine scalar_advance (uold, sold, snew, rhohalf,&
                              umac, w0, sedge, utrans, ext_scal_force, &
                              rho0_old , rho0_new , rho0_nph, &
                              rhoh0_old, rhoh0_new, &
                              rhoX0_old, rhoX0_new, &
                              p0_old, p0_new, temp0, &
                              dx,time, dt, the_bc_level,diff_coef,&
                              verbose, div_coeff, div_coeff_half, &
                              pred_vs_corr)
 
      type(multifab) , intent(inout) :: uold
      type(multifab) , intent(inout) :: sold
      type(multifab) , intent(inout) :: snew
      type(multifab) , intent(inout) :: rhohalf
      type(multifab) , intent(inout) :: umac(:)
      type(multifab) , intent(inout) :: sedge(:)
      type(multifab) , intent(inout) :: utrans(:)
      type(multifab) , intent(inout) :: ext_scal_force
!
      real(kind=dp_t), intent(inout) :: w0(:)
      real(kind=dp_t), intent(inout) :: dx(:),time,dt
      type(bc_level) , intent(in   ) :: the_bc_level
      real(kind=dp_t), intent(in   ) :: diff_coef
! 
      integer        , intent(in   ) :: verbose, pred_vs_corr
! 
      type(multifab) :: force,scal_force
      real(kind=dp_t), intent(in   ) ::  rho0_old(:)
      real(kind=dp_t), intent(in   ) ::  rho0_new(:)
      real(kind=dp_t), intent(in   ) ::  rho0_nph(:)
      real(kind=dp_t), intent(in   ) :: rhoh0_old(:)
      real(kind=dp_t), intent(in   ) :: rhoh0_new(:)
      real(kind=dp_t), intent(in   ) :: rhoX0_old(:,:)
      real(kind=dp_t), intent(in   ) :: rhoX0_new(:,:)
      real(kind=dp_t), intent(in   ) ::    p0_old(:)
      real(kind=dp_t), intent(in   ) ::    p0_new(:)
      real(kind=dp_t), intent(in   ) ::     temp0(:)
      real(kind=dp_t), intent(in   ) :: div_coeff(:)
      real(kind=dp_t), intent(in   ) :: div_coeff_half(:)
! 
      real(kind=dp_t), pointer:: uop(:,:,:,:)
      real(kind=dp_t), pointer:: ump(:,:,:,:)
      real(kind=dp_t), pointer:: vmp(:,:,:,:)
      real(kind=dp_t), pointer:: wmp(:,:,:,:)
      real(kind=dp_t), pointer:: utp(:,:,:,:)
      real(kind=dp_t), pointer:: vtp(:,:,:,:)
      real(kind=dp_t), pointer:: wtp(:,:,:,:)
      real(kind=dp_t), pointer::  ep(:,:,:,:)
      real(kind=dp_t), pointer::  fp(:,:,:,:)
!
      real(kind=dp_t), pointer:: sop(:,:,:,:)
      real(kind=dp_t), pointer:: snp(:,:,:,:)
      real(kind=dp_t), pointer::  rp(:,:,:,:) 
      real(kind=dp_t), pointer::  dp(:,:,:,:) 
      real(kind=dp_t), pointer:: sepx(:,:,:,:)
      real(kind=dp_t), pointer:: sepy(:,:,:,:)
      real(kind=dp_t), pointer:: sepz(:,:,:,:)
!
      type(multifab) :: divu
      real(dp_t)     :: mult

      integer :: nscal,nspec,velpred
      integer :: lo(uold%dim),hi(uold%dim)
      integer :: i,n,bc_comp,dm,ng_cell,ng_rho
      logical :: is_vel, make_divu, advect_in_pert_form
      logical :: do_mom
      logical, allocatable :: is_conservative(:)
      real(kind=dp_t) :: visc_fac, diff_fac
      real(kind=dp_t) :: half_time
      real(kind=dp_t), allocatable :: p0_nph(:)

      velpred = 0

      do_mom = .false.

      ng_cell = uold%ng
      ng_rho  = rhohalf%ng
      dm      = uold%dim

      half_time = time + HALF*dt

      nscal   = ncomp(ext_scal_force)
      nspec   = nscal - 2
      is_vel  = .false.

      allocate(is_conservative(nscal))
      is_conservative(1) = .true.
      is_conservative(2) = .true.

      call build(scal_force, ext_scal_force%la, nscal, 1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Scalar force at time n for density 
!
!        This is used in the prediction of the interface states.
!        Here build the interface states of the density perturbation
!        (rhopert) in convective form, so the density `force' =
!        -rhopert divu - div (rho0 u) .
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call setval(scal_force,ZERO)
      do i = 1, scal_force%nboxes
        if ( multifab_remote(scal_force, i) ) cycle
        fp  => dataptr(scal_force, i)
        sop => dataptr(sold      , i)
        ump  => dataptr(umac(1), i)
        vmp  => dataptr(umac(2), i)
        select case(dm)
        case(2)
          call modify_force_2d(fp(:,:,1,1),sop(:,:,1,1),ng_cell,rho0_old,&
                               ump(:,:,1,1),vmp(:,:,1,1),w0,dx,pred_vs_corr)
        case(3)
          wmp  => dataptr(umac(3), i)
          call modify_force_3d(fp(:,:,:,1),sop(:,:,:,1),ng_cell,rho0_old,&
                               ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),w0,dx,pred_vs_corr)
        end select
      end do

      call multifab_fill_boundary(scal_force)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create edge states of the density perturbation (rhopert) using 
!     the MAC velocity 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (pred_vs_corr .eq. 2) then
         mult = ONE
         do i = 1, umac(dm)%nboxes
            if ( multifab_remote(umac(dm), i) ) cycle
            wmp  => dataptr(umac(dm), i)
            lo =  lwb(get_box(uold, i))
            hi =  upb(get_box(uold, i))
            select case (dm)
               case (2)
                 call addw0_2d(wmp(:,:,1,1),w0,lo,hi,mult)
               case (3)
                 call addw0_3d(wmp(:,:,:,1),w0,lo,hi,mult)
            end select
         end do
      end if

      n = 1
      advect_in_pert_form = .true.
      do i = 1, sold%nboxes
         if ( multifab_remote(sold, i) ) cycle
         sop  => dataptr(sold, i)
         uop  => dataptr(uold, i)
         sepx => dataptr(sedge(1), i)
         sepy => dataptr(sedge(2), i)
         ump  => dataptr(umac(1), i)
         vmp  => dataptr(umac(2), i)
         utp  => dataptr(utrans(1), i)
         vtp  => dataptr(utrans(2), i)
          fp  => dataptr(scal_force , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call mkflux_2d(sop(:,:,1,:), uop(:,:,1,:), sop(:,:,1,1), &
                             sepx(:,:,1,:), sepy(:,:,1,:), &
                             ump(:,:,1,1), vmp(:,:,1,1), &
                             utp(:,:,1,1), vtp(:,:,1,1), fp(:,:,1,:), &
                             lo, dx, dt, is_vel, is_conservative, &
                             the_bc_level%phys_bc_level_array(i,:,:), &
                             the_bc_level%adv_bc_level_array(i,:,:,rho_comp+dm:), &
                             velpred, ng_cell, rho0_old, advect_in_pert_form, do_mom, n)
           case (3)
              sepz => dataptr( sedge(3), i)
              wmp  => dataptr(  umac(3), i)
              wtp  => dataptr(utrans(3), i)
              call mkflux_3d(sop(:,:,:,:), uop(:,:,:,:), sop(:,:,:,1), &
                             sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                             ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                             utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), fp(:,:,:,:), &
                             lo, dx, dt, is_vel, is_conservative, &
                             the_bc_level%phys_bc_level_array(i,:,:), &
                             the_bc_level%adv_bc_level_array(i,:,:,rho_comp+dm:), &
                             velpred, ng_cell, rho0_old, advect_in_pert_form, do_mom, n)
         end select
      end do

      if (pred_vs_corr .eq. 2) then
         mult = -ONE
         do i = 1, umac(dm)%nboxes
            if ( multifab_remote(umac(dm), i) ) cycle
            wmp  => dataptr(umac(dm), i)
            lo =  lwb(get_box(uold, i))
            hi =  upb(get_box(uold, i))
            select case(dm)
            case(2)
              call addw0_2d(wmp(:,:,1,1),w0,lo,hi,mult)
            case(3)
              call addw0_3d(wmp(:,:,:,1),w0,lo,hi,mult)
            end select
         end do
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Scalar force in update for density at time n+1/2 is 0.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call setval(scal_force,ZERO)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Update density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i = 1, sold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         sop => dataptr(sold, i)
         ump => dataptr(umac(1), i)
         vmp => dataptr(umac(2), i)
         sepx => dataptr(sedge(1), i)
         sepy => dataptr(sedge(2), i)
         snp => dataptr(snew, i)
          rp => dataptr(rhohalf, i)
          fp => dataptr(scal_force , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              call update_scal_2d( &
                             sop(:,:,1,1), snp(:,:,1,1), ump(:,:,1,1), vmp(:,:,1,1), w0, &
                             sepx(:,:,1,1), sepy(:,:,1,1), rp(:,:,1,1) , rho0_old, rho0_new, &
                             lo, hi, ng_cell, dx, dt, pred_vs_corr, verbose)
              call setbc_2d(snp(:,:,1,1), lo, ng_cell, &
                            the_bc_level%adv_bc_level_array(i,:,:,rho_comp+dm),dx,rho_comp+dm)
              call setbc_2d(rp(:,:,1,1), lo, ng_rho, &
                            the_bc_level%adv_bc_level_array(i,:,:,rho_comp+dm),dx,rho_comp+dm)
              call mk_rhohalf(sop(:,:,1,1),snp(:,:,1,1),rp(:,:,1,1),lo,hi,ng_cell)
           case (3)
              wmp => dataptr(umac(3), i)
              sepz => dataptr(sedge(3), i)
              call update_scal_3d( &
                             sop(:,:,:,1), snp(:,:,:,1), ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), w0, &
                             sepx(:,:,:,1), sepy(:,:,:,1), sepz(:,:,:,1), rp(:,:,:,1) , rho0_old, rho0_new, &
                             lo, hi, ng_cell, dx, dt, pred_vs_corr, verbose)
              call setbc_3d(snp(:,:,:,1), lo, ng_cell, &
                            the_bc_level%adv_bc_level_array(i,:,:,rho_comp+dm),dx,rho_comp+dm)
              call setbc_3d(rp(:,:,:,1), lo, ng_rho, &
                            the_bc_level%adv_bc_level_array(i,:,:,rho_comp+dm),dx,rho_comp+dm)
              call mk_rhohalf(sop(:,:,1,1),snp(:,:,1,1),rp(:,:,1,1),lo,hi,ng_cell)
         end select
      end do
      call multifab_fill_boundary(rhohalf)
      call multifab_fill_boundary(snew)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create scalar force at time n for (rho h) only.  First we
!     add the explicit forces (i.e. heating) through mkrhohforce and
!     then we add those advective terms that appear as forces when 
!     we write it in convective/perturbational form (through 
!     modify_force)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      visc_fac = ONE
      n = rhoh_comp
      do i = 1, scal_force%nboxes
         if ( multifab_remote(scal_force, i) ) cycle
          fp => dataptr(scal_force, i)
         sop => dataptr(sold , i)
         ump  => dataptr(umac(1), i)
         vmp  => dataptr(umac(2), i)
         lo =  lwb(get_box(sold, i))
         hi =  upb(get_box(sold, i))
         select case (dm)
            case (2)
              call  mkrhohforce_2d(fp(:,:,1,n), sop(:,:,1,n), ng_cell, & 
                                   sop(:,:,1,1), ng_cell, vmp(:,:,1,1), &
                                   dx, the_bc_level%ell_bc_level_array(i,:,:,rhoh_comp+dm), &
                                   diff_coef, visc_fac, p0_old, rho0_old, temp0, time, pred_vs_corr)
              call modify_force_2d(fp(:,:,1,n),sop(:,:,1,n),ng_cell,rhoh0_old, &
                                   ump(:,:,1,1),vmp(:,:,1,1),w0,dx,pred_vs_corr)
            case(3)
              wmp  => dataptr(umac(3), i)
              call  mkrhohforce_3d(fp(:,:,:,n), sop(:,:,:,n), ng_cell, & 
                                   sop(:,:,:,1), ng_cell, wmp(:,:,:,1), &
                                   dx, the_bc_level%ell_bc_level_array(i,:,:,rhoh_comp+dm), &
                                   diff_coef, visc_fac, p0_old, rho0_old, temp0, time, pred_vs_corr)
              call modify_force_3d(fp(:,:,:,n),sop(:,:,:,n),ng_cell,rhoh0_old, &
                                   ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),w0,dx,pred_vs_corr)
         end select
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create scalar source term at time n for (rho X)_i.  The source term
!     for (rho X) looks like (rho wdot_i).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i = 1, scal_force%nboxes
         if ( multifab_remote(scal_force, i) ) cycle
          fp => dataptr(scal_force, i)
         sop => dataptr(sold , i)
         ump => dataptr(umac(1), i)
         vmp => dataptr(umac(2), i)
         lo =  lwb(get_box(sold, i))
         hi =  upb(get_box(sold, i))
         select case (dm)
            case (2)
              n = spec_comp
              call mkspecforce_2d(fp(:,:,1,n:), sop(:,:,1,n:), ng_cell, & 
                                  sop(:,:,1,1), ng_cell, &
                                  dx, the_bc_level%ell_bc_level_array(i,:,:,n+dm:), &
                                  time, pred_vs_corr, nspec)
              do n = spec_comp,spec_comp+nspec-1
                call modify_force_2d(fp(:,:,1,n),sop(:,:,1,n),ng_cell,rhoX0_old(:,n), &
                                     ump(:,:,1,1),vmp(:,:,1,1),w0,dx,pred_vs_corr)
              end do
            case(3)
              wmp  => dataptr(umac(3), i)
              n = spec_comp
              call mkspecforce_3d(fp(:,:,:,n:), sop(:,:,:,n:), ng_cell, & 
                                  sop(:,:,:,1), ng_cell, &
                                  dx, the_bc_level%ell_bc_level_array(i,:,:,n+dm:), &
                                  time, pred_vs_corr, nspec)
              do n = spec_comp,spec_comp+nspec-1
                call modify_force_3d(fp(:,:,:,n),sop(:,:,:,n),ng_cell,rhoX0_old(:,n), &
                                     ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),w0,dx,pred_vs_corr)
              end do
         end select
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Do fill_boundary for all components of scal_force
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call multifab_fill_boundary(scal_force)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the edge states of (rho h)' and (rho X)_i using the MAC velocity 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (pred_vs_corr .eq. 2) then
         mult = ONE
         do i = 1, umac(dm)%nboxes
            if ( multifab_remote(umac(dm), i) ) cycle
            wmp  => dataptr(umac(dm), i)
            lo =  lwb(get_box(uold, i))
            hi =  upb(get_box(uold, i))
            select case(dm)
            case(2)
              call addw0_2d(wmp(:,:,1,1),w0,lo,hi,mult)
            case(3)
              call addw0_3d(wmp(:,:,:,1),w0,lo,hi,mult)
            end select
         end do
      end if

      n = nscal
      advect_in_pert_form = .true.
      do i = 1, sold%nboxes
         if ( multifab_remote(sold, i) ) cycle
         sop  => dataptr(sold, i)
         uop  => dataptr(uold, i)
         sepx => dataptr(sedge(1), i)
         sepy => dataptr(sedge(2), i)
         ump  => dataptr(umac(1), i)
         vmp  => dataptr(umac(2), i)
         utp  => dataptr(utrans(1), i)
         vtp  => dataptr(utrans(2), i)
          fp  => dataptr(scal_force , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         select case (dm)
            case (2)
              n = rhoh_comp
                bc_comp = dm+n
                call mkflux_2d(sop(:,:,1,:), uop(:,:,1,:), sop(:,:,1,1), &
                               sepx(:,:,1,:), sepy(:,:,1,:), &
                               ump(:,:,1,1), vmp(:,:,1,1), &
                               utp(:,:,1,1), vtp(:,:,1,1), fp(:,:,1,:), &
                               lo, dx, dt, is_vel, is_conservative, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%adv_bc_level_array(i,:,:,bc_comp:), &
                               velpred, ng_cell, rhoh0_old, advect_in_pert_form, do_mom, n)
              do n = spec_comp,nscal
                bc_comp = dm+n
                call mkflux_2d(sop(:,:,1,:), uop(:,:,1,:), sop(:,:,1,1), &
                               sepx(:,:,1,:), sepy(:,:,1,:), &
                               ump(:,:,1,1), vmp(:,:,1,1), &
                               utp(:,:,1,1), vtp(:,:,1,1), fp(:,:,1,:), &
                               lo, dx, dt, is_vel, is_conservative, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%adv_bc_level_array(i,:,:,bc_comp:), &
                               velpred, ng_cell, rhoX0_old(:,n), advect_in_pert_form, do_mom, n)
              end do
            case (3)
              wmp  => dataptr(  umac(3), i)
              wtp  => dataptr(utrans(3), i)
              sepz => dataptr( sedge(3), i)
              n = rhoh_comp
                bc_comp = dm+n
                call mkflux_3d(sop(:,:,:,:), uop(:,:,:,:), sop(:,:,:,1), &
                               sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                               ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                               utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), fp(:,:,:,:), &
                               lo, dx, dt, is_vel, is_conservative, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%adv_bc_level_array(i,:,:,bc_comp:), &
                               velpred, ng_cell, rhoh0_old, advect_in_pert_form, do_mom, n)
              do n = spec_comp,nscal
                bc_comp = dm+n
                call mkflux_3d(sop(:,:,:,:), uop(:,:,:,:), sop(:,:,:,1), &
                               sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                               ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                               utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), fp(:,:,:,:), &
                               lo, dx, dt, is_vel, is_conservative, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%adv_bc_level_array(i,:,:,bc_comp:), &
                               velpred, ng_cell, rhoX0_old(:,n), advect_in_pert_form, do_mom, n)
              end do
         end select
      end do

      if (pred_vs_corr .eq. 2) then
         mult = -ONE
         do i = 1, umac(dm)%nboxes
            if ( multifab_remote(umac(dm), i) ) cycle
            wmp  => dataptr(umac(dm), i)
            lo =  lwb(get_box(uold, i))
            hi =  upb(get_box(uold, i))
            select case(dm)
            case(2)
              call addw0_2d(wmp(:,:,1,1),w0,lo,hi,mult)
            case(3)
              call addw0_3d(wmp(:,:,:,1),w0,lo,hi,mult)
            end select
         end do
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create (rhoh)' force at time n+1/2.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      diff_fac = HALF
      allocate(p0_nph(lo(dm):hi(dm)))
      p0_nph = HALF * (p0_old + p0_new)
      do i = 1, uold%nboxes
         if ( multifab_remote(uold, i) ) cycle
          fp => dataptr(scal_force, i)
          rp => dataptr(rhohalf, i)
          sop => dataptr(sold , i)
          lo =  lwb(get_box(sold, i))
          hi =  upb(get_box(sold, i))
         select case (dm)
            case (2)
              vmp  => dataptr(umac(2), i)
              n = rhoh_comp
              call mkrhohforce_2d(fp(:,:,1,n), sop(:,:,1,n), ng_cell, &
                                  rp(:,:,1,1), ng_rho, vmp(:,:,1,1), &
                                  dx, the_bc_level%ell_bc_level_array(i,:,:,rhoh_comp+dm), &
                                  diff_coef, diff_fac, p0_nph, rho0_nph, temp0, half_time, pred_vs_corr)
            case (3)
              wmp  => dataptr(umac(3), i)
              n = rhoh_comp
              call mkrhohforce_3d(fp(:,:,:,n), sop(:,:,:,n), ng_cell, &
                                  rp(:,:,:,1), ng_rho, wmp(:,:,:,1), &
                                  dx, the_bc_level%ell_bc_level_array(i,:,:,rhoh_comp+dm), &
                                  diff_coef, diff_fac, p0_nph, rho0_nph, temp0, half_time, pred_vs_corr)
         end select
      end do
      deallocate(p0_nph)
      call multifab_fill_boundary(scal_force)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create scalar source term at time n+1/2 for (rho X)_i.  The source term
!     for (rho X) looks like (rho wdot_i).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      n = spec_comp
      do i = 1, scal_force%nboxes
         if ( multifab_remote(scal_force, i) ) cycle
          fp => dataptr(scal_force, i)
         sop => dataptr(sold , i)
         ump => dataptr(umac(1), i)
         vmp => dataptr(umac(2), i)
          rp => dataptr(rhohalf, i)
         lo =  lwb(get_box(sold, i))
         hi =  upb(get_box(sold, i))
         select case (dm)
            case (2)
              call mkspecforce_2d(fp(:,:,1,n:), sop(:,:,1,n:), ng_cell, & 
                                  rp(:,:,1,1), ng_cell, &
                                  dx, the_bc_level%ell_bc_level_array(i,:,:,n+dm:), &
                                  time, pred_vs_corr, nspec)
            case(3)
              wmp  => dataptr(umac(3), i)
              call mkspecforce_3d(fp(:,:,:,n:), sop(:,:,:,n:), ng_cell, & 
                                  rp(:,:,:,1), ng_cell, &
                                  dx, the_bc_level%ell_bc_level_array(i,:,:,n+dm:), &
                                  time, pred_vs_corr, nspec)
         end select
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Update (rho h)' and (rho X)'_i with conservative differencing.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i = 1, sold%nboxes
         if ( multifab_remote(uold, i) ) cycle
         sop => dataptr(sold, i)
         snp => dataptr(snew, i)
         ump => dataptr(umac(1), i)
         vmp => dataptr(umac(2), i)
         sepx => dataptr(sedge(1), i)
         sepy => dataptr(sedge(2), i)
          fp => dataptr(scal_force , i)
         lo =  lwb(get_box(uold, i))
         hi =  upb(get_box(uold, i))
         n = rhoh_comp
         select case (dm)
            case (2)
              call update_scal_2d( &
                             sop(:,:,1,n), snp(:,:,1,n), ump(:,:,1,1), vmp(:,:,1,1), w0, &
                             sepx(:,:,1,n), sepy(:,:,1,n), fp(:,:,1,n), &
                             rhoh0_old, rhoh0_new, &
                             lo, hi, ng_cell, dx, dt, pred_vs_corr, verbose)
              call setbc_2d(snp(:,:,1,n), lo, ng_cell, &
                            the_bc_level%adv_bc_level_array(i,:,:,rhoh_comp),dx,n+dm)
              do n = spec_comp,nscal
                call update_scal_2d( &
                               sop(:,:,1,n), snp(:,:,1,n), &
                               ump(:,:,1,1), vmp(:,:,1,1), w0, &
                               sepx(:,:,1,n), sepy(:,:,1,n), fp(:,:,1,n), &
                               rhoX0_old(:,n), rhoX0_new(:,n), &
                               lo, hi, ng_cell, dx, dt, pred_vs_corr, verbose)
                call setbc_2d(snp(:,:,1,n), lo, ng_cell, &
                              the_bc_level%adv_bc_level_array(i,:,:,n+dm),dx,n+dm)
              end do
            case (3)
              wmp => dataptr(umac(3), i)
              sepz => dataptr(sedge(3), i)
              call update_scal_3d( &
                             sop(:,:,:,n), snp(:,:,:,n), ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), w0, &
                             sepx(:,:,:,n), sepy(:,:,:,n), sepz(:,:,:,n), fp(:,:,:,n), &
                             rhoh0_old, rhoh0_new, &
                             lo, hi, ng_cell, dx, dt, pred_vs_corr, verbose)
              call setbc_3d(snp(:,:,:,2), lo, ng_cell, & 
                            the_bc_level%adv_bc_level_array(i,:,:,rhoh_comp),dx,n+dm)
              do n = spec_comp,nscal
                call update_scal_3d( &
                               sop(:,:,:,n), snp(:,:,:,n), &
                               ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), w0, &
                               sepx(:,:,:,n), sepy(:,:,:,n), sepz(:,:,:,n), fp(:,:,:,n), &
                               rhoX0_old(:,n), rhoX0_new(:,n), &
                               lo, hi, ng_cell, dx, dt, pred_vs_corr, verbose)
                call setbc_3d(snp(:,:,:,n), lo, ng_cell, &
                              the_bc_level%adv_bc_level_array(i,:,:,n+dm),dx,n+dm)
              end do
         end select
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Call fill_boundary for all components of snew
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call multifab_fill_boundary(snew)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      deallocate(is_conservative)
      call multifab_destroy(scal_force)

   end subroutine scalar_advance

   subroutine modify_force_2d(force,s,ng,base,umac,vmac,w0,dx,pred_vs_corr)

    ! When we write the scalar equation in perturbational and convective
    ! form, the terms other than s'_t + U.grad s' act as source terms.  Add
    ! them to the forces here.

    integer        , intent(in   ) :: ng,pred_vs_corr
    real(kind=dp_t), intent(  out) :: force(0:,0:)
    real(kind=dp_t), intent(in   ) :: s(1-ng:,1-ng:)
    real(kind=dp_t), intent(in   ) :: base(:)
    real(kind=dp_t), intent(in   ) :: umac(0:,0:)
    real(kind=dp_t), intent(in   ) :: vmac(0:,0:)
    real(kind=dp_t), intent(in   ) :: w0(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    integer :: i,j,nx,ny
    real(kind=dp_t) :: divu,divbaseu
    real(kind=dp_t) :: base_half_lo,base_half_hi
    
    nx = size(force,dim=1)-2
    ny = size(force,dim=2)-2

    do j = 1,ny
       if (j.eq.1) then
          base_half_lo = base(j)
          base_half_hi = HALF*(base(j)+base(j+1))
       else if (j.eq.ny) then
          base_half_lo = HALF*(base(j)+base(j-1))
          base_half_hi = base(j)
       else if (j.eq.2) then
          base_half_lo = HALF*(base(j)+base(j-1))
          base_half_hi = 7.d0/12.d0 * (base(j  )+base(j+1)) &
               -1.d0/12.d0 * (base(j-1)+base(j+2))
       else if (j.eq.ny-1) then
          base_half_lo = 7.d0/12.d0 * (base(j  )+base(j-1)) &
               -1.d0/12.d0 * (base(j+1)+base(j-2))
          base_half_hi = HALF*(base(j)+base(j+1))
       else 
          base_half_lo = 7.d0/12.d0 * (base(j  )+base(j-1)) &
               -1.d0/12.d0 * (base(j+1)+base(j-2))
          base_half_hi = 7.d0/12.d0 * (base(j  )+base(j+1)) &
               -1.d0/12.d0 * (base(j-1)+base(j+2))
       end if
        do i = 1,nx
           divu = (umac(i+1,j) - umac(i,j)) / dx(1) &
                 +(vmac(i,j+1) - vmac(i,j)) / dx(2)
           divbaseu = base(j)*(umac(i+1,j) - umac(i,j))/dx(1) &
                             +(vmac(i,j+1) * base_half_hi &
                             - vmac(i,j  ) * base_half_lo)/ dx(2)
           if (pred_vs_corr .eq. 2) then
              divu = divu + (w0(j+1)-w0(j))/dx(2)
           end if
           force(i,j) = force(i,j) - (s(i,j)-base(j))*divu - divbaseu
        end do
     end do
     
   end subroutine modify_force_2d

   subroutine modify_force_3d(force,s,ng,base,umac,vmac,wmac,w0,dx,pred_vs_corr)

    ! When we write the scalar equation in perturbational and convective
    ! form, the terms other than s'_t + U.grad s' act as source terms.  Add
    ! them to the forces here.

    integer        , intent(in   ) :: ng,pred_vs_corr
    real(kind=dp_t), intent(  out) :: force(0:,0:,0:)
    real(kind=dp_t), intent(in   ) :: s(1-ng:,1-ng:,1-ng:)
    real(kind=dp_t), intent(in   ) :: base(:)
    real(kind=dp_t), intent(in   ) :: umac(0:,0:,0:)
    real(kind=dp_t), intent(in   ) :: vmac(0:,0:,0:)
    real(kind=dp_t), intent(in   ) :: wmac(0:,0:,0:)
    real(kind=dp_t), intent(in   ) :: w0(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    integer :: i,j,k,nx,ny,nz
    real(kind=dp_t) :: divu,divbaseu
    real(kind=dp_t) :: base_half_lo,base_half_hi
    
    nx = size(force,dim=1)-2
    ny = size(force,dim=2)-2
    nz = size(force,dim=3)-2

    do k = 1,nz
       if (k.eq.1) then
          base_half_lo = base(k)
          base_half_hi = HALF*(base(k)+base(k+1))
       else if (k.eq.nz) then
          base_half_lo = HALF*(base(k)+base(k-1))
          base_half_hi = base(k)
       else if (k.eq.2) then
          base_half_lo = HALF*(base(k)+base(k-1))
          base_half_hi = 7.d0/12.d0 * (base(k  )+base(k+1)) &
               -1.d0/12.d0 * (base(k-1)+base(k+2))
       else if (k.eq.nz-1) then
          base_half_lo = 7.d0/12.d0 * (base(k  )+base(k-1)) &
               -1.d0/12.d0 * (base(k+1)+base(k-2))
          base_half_hi = HALF*(base(k)+base(k+1))
       else 
          base_half_lo = 7.d0/12.d0 * (base(k  )+base(k-1)) &
               -1.d0/12.d0 * (base(k+1)+base(k-2))
          base_half_hi = 7.d0/12.d0 * (base(k  )+base(k+1)) &
               -1.d0/12.d0 * (base(k-1)+base(k+2))
       end if
        do j = 1,ny
        do i = 1,nx
           divu = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) &
                 +(vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) &
                 +(wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)
           divbaseu = base(k)*( (umac(i+1,j,k) - umac(i,j,k))/dx(1) &
                               +(vmac(i,j+1,k) - vmac(i,j,k))/dx(2) ) &
                               +(wmac(i,j,k+1) * base_half_hi &
                               - wmac(i,j,k  ) * base_half_lo)/ dx(3)
           if (pred_vs_corr .eq. 2) then
              divu = divu + (w0(k+1)-w0(k))/dx(3)
           end if
           force(i,j,k) = force(i,j,k) - (s(i,j,k)-base(k))*divu - divbaseu
        end do
        end do
     end do
     
   end subroutine modify_force_3d

end module scalar_advance_module
