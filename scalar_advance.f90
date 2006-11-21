module scalar_advance_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use mkflux_module
  use mkscalforce_module
  use update_module
  use addw0_module
  use define_bc_module
  use setbc_module
  use fill_3d_module
  use variables
  use geometry
  use network

  implicit none

contains

   subroutine scalar_advance (uold, sold, snew, &
                              umac, w0, sedge, utrans, ext_scal_force, normal, &
                              s0_old , s0_new , &
                              p0_old, p0_new, &
                              dx,time, dt, the_bc_level, &
                              verbose)
 
      type(multifab) , intent(inout) :: uold
      type(multifab) , intent(inout) :: sold
      type(multifab) , intent(inout) :: snew
      type(multifab) , intent(inout) :: umac(:)
      type(multifab) , intent(inout) :: sedge(:)
      type(multifab) , intent(inout) :: utrans(:)
      type(multifab) , intent(inout) :: ext_scal_force
      type(multifab) , intent(in   ) :: normal
!
      real(kind=dp_t), intent(inout) :: w0(:)
      real(kind=dp_t), intent(in   ) :: dx(:),time,dt
      type(bc_level) , intent(in   ) :: the_bc_level
! 
      integer        , intent(in   ) :: verbose
! 
      type(multifab) :: force,scal_force
      real(kind=dp_t), intent(in   ) ::  s0_old(:,:)
      real(kind=dp_t), intent(in   ) ::  s0_new(:,:)
      real(kind=dp_t), intent(in   ) ::    p0_old(:)
      real(kind=dp_t), intent(in   ) ::    p0_new(:)
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
      real(kind=dp_t), pointer::  np(:,:,:,:)
!
      real(kind=dp_t), pointer:: sop(:,:,:,:)
      real(kind=dp_t), pointer:: snp(:,:,:,:)
      real(kind=dp_t), pointer::  dp(:,:,:,:) 
      real(kind=dp_t), pointer:: sepx(:,:,:,:)
      real(kind=dp_t), pointer:: sepy(:,:,:,:)
      real(kind=dp_t), pointer:: sepz(:,:,:,:)
!
      type(multifab) :: divu
      real(dp_t)     :: mult

      integer :: nscal,ntrac,velpred
      integer :: lo(uold%dim),hi(uold%dim)
      integer :: i,n,bc_comp,dm,ng_cell
      logical :: is_vel, make_divu, advect_in_pert_form
      logical, allocatable :: is_conservative(:)
      real(dp_t), allocatable :: s0_cart(:,:,:)
      real(kind=dp_t) :: half_time

      print *,'<<< advect state >>> '

      velpred = 0

      ng_cell = sold%ng
      dm      = sold%dim

      half_time = time + HALF*dt

      nscal  = ncomp(ext_scal_force)
      ntrac  = nscal - nspec - 2
      is_vel = .false.

      allocate(is_conservative(nscal))
      is_conservative(1) = .true.
      is_conservative(2) = .true.
      is_conservative(spec_comp:spec_comp+nspec-1) = .true.
      if (ntrac .ge. 1) &
        is_conservative(trac_comp:trac_comp+ntrac-1) = .false.

      call build(scal_force, ext_scal_force%la, nscal, 1)
      call setval(scal_force,ZERO)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Add w0 to radial velocity.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      mult = ONE
      call addw0(umac,normal,w0,dx,mult)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create scalar source term at time n for (rho X)_i and (rho H).  
!     The source term for (rho X) is zero.
!     The source term for (rho h) has only the w dp0/dr term.

!     The call to modify_force is used to add those advective terms 
!     that appear as forces when we write it in convective/perturbational form.
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

              do n = spec_comp,spec_comp+nspec-1
                call modify_force_2d(fp(:,:,1,n),sop(:,:,1,n),ng_cell,&
                                     s0_old(:,n), &
                                     ump(:,:,1,1),vmp(:,:,1,1),dx)
              end do

              n = rhoh_comp
              call  mkrhohforce_2d(fp(:,:,1,n), vmp(:,:,1,1), p0_old, p0_new, dx(dm))
              call modify_force_2d(fp(:,:,1,n),sop(:,:,1,n),ng_cell,s0_old(:,rhoh_comp), &
                                   ump(:,:,1,1),vmp(:,:,1,1),dx)

            case(3)
              wmp  => dataptr(umac(3), i)

              if (spherical .eq. 1) then
                allocate(s0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
                do n = spec_comp,spec_comp+nspec-1
                  call fill_3d_data(s0_cart,s0_old(:,n),dx,ng_cell)
                  call modify_force_3d_sphr(fp(:,:,:,n),sop(:,:,:,n),ng_cell,&
                                            s0_cart, &
                                            ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),w0,dx)
                end do

                n = rhoh_comp
                np => dataptr(normal, i)
                call  mkrhohforce_3d_sphr(fp(:,:,:,n), &
                                          ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                          np(:,:,:,:), p0_old, p0_new, dx)

                call fill_3d_data(s0_cart,s0_old(:,n),dx,ng_cell)
                call modify_force_3d_sphr(fp(:,:,:,n),sop(:,:,:,n),ng_cell,s0_cart, &
                                          ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),w0,dx)
                deallocate(s0_cart)
              else
                do n = spec_comp,spec_comp+nspec-1
                  call modify_force_3d_cart(fp(:,:,:,n),sop(:,:,:,n),ng_cell,&
                                            s0_old(:,n), &
                                            ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),w0,dx)
                end do

                n = rhoh_comp
                call  mkrhohforce_3d(fp(:,:,:,n), wmp(:,:,:,1), p0_old, p0_new, dx(dm))

                call modify_force_3d_cart(fp(:,:,:,n),sop(:,:,:,n),ng_cell,s0_old(:,n), &
                                          ump(:,:,:,1),vmp(:,:,:,1),wmp(:,:,:,1),w0,dx)
              end if
         end select
      end do

      call multifab_fill_boundary(scal_force)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the edge states of (rho h)' and (rho X)_i using the MAC velocity 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
                call mkflux_2d(sop(:,:,1,:), uop(:,:,1,:), &
                               sepx(:,:,1,:), sepy(:,:,1,:), &
                               ump(:,:,1,1), vmp(:,:,1,1), &
                               utp(:,:,1,1), vtp(:,:,1,1), fp(:,:,1,:), &
                               lo, dx, dt, is_vel, is_conservative, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%adv_bc_level_array(i,:,:,bc_comp:), &
                               velpred, ng_cell, s0_old(:,n), &
                               advect_in_pert_form, n)

              do n = spec_comp,spec_comp+nspec-1
                bc_comp = dm+n
                call mkflux_2d(sop(:,:,1,:), uop(:,:,1,:), &
                               sepx(:,:,1,:), sepy(:,:,1,:), &
                               ump(:,:,1,1), vmp(:,:,1,1), &
                               utp(:,:,1,1), vtp(:,:,1,1), fp(:,:,1,:), &
                               lo, dx, dt, is_vel, is_conservative, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%adv_bc_level_array(i,:,:,bc_comp:), &
                               velpred, ng_cell, s0_old(:,n), &
                               advect_in_pert_form, n)
              end do
            case (3)
              wmp  => dataptr(  umac(3), i)
              wtp  => dataptr(utrans(3), i)
              sepz => dataptr( sedge(3), i)
              n = rhoh_comp
                bc_comp = dm+n
                call mkflux_3d(sop(:,:,:,:), uop(:,:,:,:), &
                               sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                               ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                               utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), fp(:,:,:,:), &
                               lo, dx, dt, is_vel, is_conservative, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%adv_bc_level_array(i,:,:,bc_comp:), &
                               velpred, ng_cell, s0_old(:,n), &
                                advect_in_pert_form, n)

              do n = spec_comp,spec_comp+nspec-1
                bc_comp = dm+n
                call mkflux_3d(sop(:,:,:,:), uop(:,:,:,:), &
                               sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                               ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                               utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), fp(:,:,:,:), &
                               lo, dx, dt, is_vel, is_conservative, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%adv_bc_level_array(i,:,:,bc_comp:), &
                               velpred, ng_cell, s0_old(:,n), &
                               advect_in_pert_form, n)
              end do
         end select
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the edge states of tracers using the MAC velocity 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (ntrac .ge. 1) then
      advect_in_pert_form = .false.
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
              do n = trac_comp,trac_comp+ntrac-1
                bc_comp = dm+n
                call mkflux_2d(sop(:,:,1,:), uop(:,:,1,:), &
                               sepx(:,:,1,:), sepy(:,:,1,:), &
                               ump(:,:,1,1), vmp(:,:,1,1), &
                               utp(:,:,1,1), vtp(:,:,1,1), fp(:,:,1,:), &
                               lo, dx, dt, is_vel, is_conservative, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%adv_bc_level_array(i,:,:,bc_comp:), &
                               velpred, ng_cell, s0_old(:,n), &
                               advect_in_pert_form, n)
              end do
            case (3)
              wmp  => dataptr(  umac(3), i)
              wtp  => dataptr(utrans(3), i)
              sepz => dataptr( sedge(3), i)
              do n = trac_comp,trac_comp+ntrac-1
                bc_comp = dm+n
                call mkflux_3d(sop(:,:,:,:), uop(:,:,:,:), &
                               sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                               ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                               utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), fp(:,:,:,:), &
                               lo, dx, dt, is_vel, is_conservative, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%adv_bc_level_array(i,:,:,bc_comp:), &
                               velpred, ng_cell, s0_old(:,n), &
                               advect_in_pert_form, n)
              end do
         end select
      end do
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Subtract w0 from radial velocity.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      mult = -ONE
      call addw0(umac,normal,w0,dx,mult)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     1) Set force for (rho X)_i at time n+1/2 = 0.
!     2) Update (rho X)'_i with conservative differencing.
!     3) Define density as the sum of the (rho X)_i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call setval(scal_force,ZERO)

      do i = 1, sold%nboxes
         if ( multifab_remote(sold, i) ) cycle
         sop => dataptr(sold, i)
         snp => dataptr(snew, i)
         ump => dataptr(umac(1), i)
         vmp => dataptr(umac(2), i)
         sepx => dataptr(sedge(1), i)
         sepy => dataptr(sedge(2), i)
          fp => dataptr(scal_force , i)
         lo =  lwb(get_box(sold, i))
         hi =  upb(get_box(sold, i))
         select case (dm)
            case (2)
                call update_scal_2d(spec_comp, spec_comp+nspec-1, &
                                    sop(:,:,1,:), snp(:,:,1,:), &
                                    ump(:,:,1,1), vmp(:,:,1,1), w0, &
                                    sepx(:,:,1,:), sepy(:,:,1,:), fp(:,:,1,:), &
                                    s0_old(:,:), s0_new(:,:), &
                                    lo, hi, ng_cell, dx, dt, verbose)
              do n = spec_comp,spec_comp+nspec-1
                bc_comp = dm+n
                call setbc_2d(snp(:,:,1,n), lo, ng_cell, &
                              the_bc_level%adv_bc_level_array(i,:,:,bc_comp),dx,bc_comp)
              end do
              ! Dont forget to call setbc for density also
              n = rho_comp
              bc_comp = dm+n
              call setbc_2d(snp(:,:,1,n), lo, ng_cell, &
                            the_bc_level%adv_bc_level_array(i,:,:,bc_comp),dx,bc_comp)
            case (3)
               wmp => dataptr(umac(3), i)
              sepz => dataptr(sedge(3), i)
                np => dataptr(normal , i)
              call update_scal_3d(spec_comp, spec_comp+nspec-1, &
                                  sop(:,:,:,:), snp(:,:,:,:), &
                                  ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), w0, &
                                  sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                  np(:,:,:,:), fp(:,:,:,:), &
                                  s0_old(:,:), s0_new(:,:), &
                                  lo, hi, ng_cell, dx, dt, verbose)
              do n = spec_comp,spec_comp+nspec-1
                bc_comp = dm+n
                call setbc_3d(snp(:,:,:,n), lo, ng_cell, &
                              the_bc_level%adv_bc_level_array(i,:,:,bc_comp),dx,bc_comp)
              end do
              ! Dont forget to call setbc for density also
              n = rho_comp
              bc_comp = dm+n
              call setbc_3d(snp(:,:,:,n), lo, ng_cell, &
                            the_bc_level%adv_bc_level_array(i,:,:,bc_comp),dx,bc_comp)
         end select
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     2) Update tracers with convective differencing.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (ntrac .ge. 1) then
      do i = 1, sold%nboxes
         if ( multifab_remote(sold, i) ) cycle
         sop => dataptr(sold, i)
         snp => dataptr(snew, i)
         ump => dataptr(umac(1), i)
         vmp => dataptr(umac(2), i)
         sepx => dataptr(sedge(1), i)
         sepy => dataptr(sedge(2), i)
          fp => dataptr(scal_force , i)
         lo =  lwb(get_box(sold, i))
         hi =  upb(get_box(sold, i))
         select case (dm)
            case (2)
              call update_scal_2d(trac_comp,trac_comp+ntrac-1, &
                             sop(:,:,1,:), snp(:,:,1,:), &
                             ump(:,:,1,1), vmp(:,:,1,1), w0, &
                             sepx(:,:,1,:), sepy(:,:,1,:), fp(:,:,1,:), &
                             s0_old(:,:), s0_new(:,:), &
                             lo, hi, ng_cell, dx, dt, verbose)
              do n = trac_comp,trac_comp+ntrac-1
                bc_comp = dm+n
                call setbc_2d(snp(:,:,1,n), lo, ng_cell, &
                              the_bc_level%adv_bc_level_array(i,:,:,bc_comp),dx,bc_comp)
              end do
            case (3)
               wmp => dataptr(umac(3), i)
              sepz => dataptr(sedge(3), i)
                np => dataptr(normal , i)
              call update_scal_3d(trac_comp,trac_comp+ntrac-1, &
                                  sop(:,:,:,:), snp(:,:,:,:), &
                                  ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), w0, &
                                  sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                  np(:,:,:,:), fp(:,:,:,:), &
                                  s0_old(:,:), s0_new(:,:), &
                                  lo, hi, ng_cell, dx, dt, verbose)
              do n = trac_comp,trac_comp+ntrac-1
                bc_comp = dm+n
                call setbc_3d(snp(:,:,:,n), lo, ng_cell, &
                              the_bc_level%adv_bc_level_array(i,:,:,bc_comp),dx,bc_comp)
              end do
         end select
      end do
      end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     1) Create (rhoh)' force at time n+1/2.
!          (NOTE: we don't worry about filling ghost cells of the scal_force
!                 because we only need them in valid regions...)     
!     2) Update (rho h)' with conservative differencing.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      n = rhoh_comp
      bc_comp = dm+n

      do i = 1, sold%nboxes
         if ( multifab_remote(sold, i) ) cycle
         sop => dataptr(sold, i)
         snp => dataptr(snew, i)
         ump => dataptr(umac(1), i)
         vmp => dataptr(umac(2), i)
         sepx => dataptr(sedge(1), i)
         sepy => dataptr(sedge(2), i)
          fp => dataptr(scal_force , i)
         lo =  lwb(get_box(sold, i))
         hi =  upb(get_box(sold, i))
         select case (dm)
            case (2)
              call  mkrhohforce_2d(fp(:,:,1,n), vmp(:,:,1,1), p0_old, p0_new, dx(dm))

              call update_scal_2d(rhoh_comp, rhoh_comp, &
                             sop(:,:,1,:), snp(:,:,1,:), &
                             ump(:,:,1,1), vmp(:,:,1,1), w0, &
                             sepx(:,:,1,:), sepy(:,:,1,:), fp(:,:,1,:), &
                             s0_old(:,:), s0_new(:,:), &
                             lo, hi, ng_cell, dx, dt, verbose)

              call setbc_2d(snp(:,:,1,n), lo, ng_cell, &
                            the_bc_level%adv_bc_level_array(i,:,:,bc_comp),dx,bc_comp)

            case(3)
              wmp  => dataptr(umac(3), i)
              sepz => dataptr(sedge(3), i)

              if (spherical .eq. 1) then

                np => dataptr(normal, i)
                call  mkrhohforce_3d_sphr(fp(:,:,:,n), &
                                          ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                          np(:,:,:,:), p0_old, p0_new, dx)

              else

                call  mkrhohforce_3d(fp(:,:,:,n), wmp(:,:,:,1), p0_old, p0_new, dx(dm))

              end if

               np => dataptr(normal , i)
              call update_scal_3d(rhoh_comp, rhoh_comp, &
                             sop(:,:,:,:), snp(:,:,:,:), &
                             ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), w0, &
                             sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                             np(:,:,:,:), fp(:,:,:,:), &
                             s0_old(:,:), s0_new(:,:), &
                             lo, hi, ng_cell, dx, dt, verbose)

              call setbc_3d(snp(:,:,:,n), lo, ng_cell, & 
                            the_bc_level%adv_bc_level_array(i,:,:,bc_comp),dx,bc_comp)
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

   subroutine modify_force_2d(force,s,ng,base,umac,vmac,dx)

    ! When we write the scalar equation in perturbational and convective
    ! form, the terms other than s'_t + U.grad s' act as source terms.  Add
    ! them to the forces here.

    ! Note that the MAC velocity has w0 already added to it here.

    integer        , intent(in   ) :: ng
    real(kind=dp_t), intent(  out) :: force(0:,0:)
    real(kind=dp_t), intent(in   ) :: s(1-ng:,1-ng:)
    real(kind=dp_t), intent(in   ) :: base(:)
    real(kind=dp_t), intent(in   ) :: umac(0:,0:)
    real(kind=dp_t), intent(in   ) :: vmac(0:,0:)
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
           force(i,j) = force(i,j) - (s(i,j)-base(j))*divu - divbaseu
        end do
     end do
     
   end subroutine modify_force_2d

   subroutine modify_force_3d_cart(force,s,ng,base,umac,vmac,wmac,w0,dx)

    ! When we write the scalar equation in perturbational and convective
    ! form, the terms other than s'_t + U.grad s' act as source terms.  Add
    ! them to the forces here.

    integer        , intent(in   ) :: ng
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
           divu = divu + (w0(k+1)-w0(k))/dx(3)
           force(i,j,k) = force(i,j,k) - (s(i,j,k)-base(k))*divu - divbaseu
        end do
        end do
     end do
     
   end subroutine modify_force_3d_cart

   subroutine modify_force_3d_sphr(force,s,ng,base_cart,umac,vmac,wmac,w0,dx)

    ! When we write the scalar equation in perturbational and convective
    ! form, the terms other than s'_t + U.grad s' act as source terms.  Add
    ! them to the forces here.

    integer        , intent(in   ) :: ng
    real(kind=dp_t), intent(  out) :: force(0:,0:,0:)
    real(kind=dp_t), intent(in   ) :: s(1-ng:,1-ng:,1-ng:)
    real(kind=dp_t), intent(in   ) :: base_cart(:,:,:)
    real(kind=dp_t), intent(in   ) :: umac(0:,0:,0:)
    real(kind=dp_t), intent(in   ) :: vmac(0:,0:,0:)
    real(kind=dp_t), intent(in   ) :: wmac(0:,0:,0:)
    real(kind=dp_t), intent(in   ) :: w0(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    ! Local variables
    integer :: i,j,k,nx,ny,nz,nr
    real(kind=dp_t) :: divumac,divbaseu
    real(kind=dp_t) :: base_xlo,base_xhi
    real(kind=dp_t) :: base_ylo,base_yhi
    real(kind=dp_t) :: base_zlo,base_zhi

    real(kind=dp_t), allocatable :: divu(:),divu_cart(:,:,:)
    
    nx = size(force,dim=1)-2
    ny = size(force,dim=2)-2
    nz = size(force,dim=3)-2

    nr = size(w0,dim=1)-1
    allocate(divu(nr))
    allocate(divu_cart(nx,ny,nz))

    do k = 1,nr
      divu(k) = (zl(k+1)**2 * w0(k+1)- zl(k)**2 * w0(k))/(dr*z(k)**2)
    end do
    call fill_3d_data(divu_cart,divu,dx,0)

    do k = 1,nz
        do j = 1,ny
        do i = 1,nx

           divumac = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) &
                    +(vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) &
                    +(wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)

           if (i.lt.nx) then
             base_xhi = HALF * (base_cart(i,j,k) + base_cart(i+1,j,k))
           else
             base_xhi = base_cart(i,j,k)
           end if
           if (i.gt.1) then
             base_xlo = HALF * (base_cart(i,j,k) + base_cart(i-1,j,k))
           else
             base_xlo = base_cart(i,j,k)
           end if
           if (j.lt.ny) then
             base_yhi = HALF * (base_cart(i,j,k) + base_cart(i,j+1,k))
           else
             base_yhi = base_cart(i,j,k)
           end if
           if (j.gt.1) then
             base_ylo = HALF * (base_cart(i,j,k) + base_cart(i,j-1,k))
           else
             base_ylo = base_cart(i,j,k)
           end if
           if (k.lt.nz) then
             base_zhi = HALF * (base_cart(i,j,k) + base_cart(i,j,k+1))
           else
             base_zhi = base_cart(i,j,k)
           end if
           if (k.gt.1) then
             base_zlo = HALF * (base_cart(i,j,k) + base_cart(i,j,k-1))
           else
             base_zlo = base_cart(i,j,k)
           end if

           divbaseu =  (umac(i+1,j,k) * base_xhi &
                      - umac(i  ,j,k) * base_xlo)/ dx(3) &
                      +(vmac(i,j+1,k) * base_yhi &
                       -vmac(i,j  ,k) * base_ylo)/ dx(3) &
                      +(wmac(i,j,k+1) * base_zhi &
                       -wmac(i,j,k  ) * base_zlo)/ dx(3)

           force(i,j,k) = force(i,j,k) - divbaseu &
                          -(s(i,j,k)-base_cart(i,j,k))*(divumac+divu_cart(i,j,k)) 
        end do
        end do
     end do

     deallocate(divu,divu_cart)
     
   end subroutine modify_force_3d_sphr

end module scalar_advance_module
