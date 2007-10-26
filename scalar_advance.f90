module scalar_advance_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use mkflux_module
  use mkscalforce_module
  use update_scal_module
  use addw0_module
  use define_bc_module
  use setbc_module
  use fill_3d_module
  use pert_form_module
  use cell_to_edge_module
  use rhoh_vs_t_module
  use variables
  use geometry
  use network
  use probin_module

  implicit none

contains

   subroutine scalar_advance (which_step, uold, sold, snew, thermal, &
                              umac, w0, w0_cart_vec, sedge, utrans, &
                              ext_scal_force, normal, &
                              s0_old , s0_new , &
                              p0_old, p0_new, &
                              dx, dt, the_bc_level, &
                              verbose)
 
      integer        , intent(in   ) :: which_step
      type(multifab) , intent(inout) :: uold
      type(multifab) , intent(inout) :: sold
      type(multifab) , intent(inout) :: snew
      type(multifab) , intent(inout) :: thermal
      type(multifab) , intent(inout) :: umac(:)
      type(multifab) , intent(inout) :: sedge(:)
      type(multifab) , intent(inout) :: utrans(:)
      type(multifab) , intent(inout) :: ext_scal_force
      type(multifab) , intent(in   ) :: normal
!
      real(kind=dp_t), intent(inout) :: w0(0:)
      type(multifab) , intent(in   ) :: w0_cart_vec
      real(kind=dp_t), intent(in   ) :: dx(:),dt
      type(bc_level) , intent(in   ) :: the_bc_level
! 
      integer        , intent(in   ) :: verbose
!
      type(multifab) :: scal_force,s0_cart

      real(kind=dp_t), intent(in   ) ::  s0_old(0:,:)
      real(kind=dp_t), intent(in   ) ::  s0_new(0:,:)
      real(kind=dp_t), intent(in   ) ::    p0_old(0:)
      real(kind=dp_t), intent(in   ) ::    p0_new(0:)

      real(kind=dp_t), allocatable :: s0_edge_old(:,:)
      real(kind=dp_t), allocatable :: s0_edge_new(:,:)
! 
      real(kind=dp_t), pointer:: uop(:,:,:,:)
      real(kind=dp_t), pointer:: ump(:,:,:,:)
      real(kind=dp_t), pointer:: vmp(:,:,:,:)
      real(kind=dp_t), pointer:: wmp(:,:,:,:)
      real(kind=dp_t), pointer:: utp(:,:,:,:)
      real(kind=dp_t), pointer:: vtp(:,:,:,:)
      real(kind=dp_t), pointer:: wtp(:,:,:,:)
      real(kind=dp_t), pointer:: w0p(:,:,:,:)
      real(kind=dp_t), pointer::  fp(:,:,:,:)
      real(kind=dp_t), pointer::  dp(:,:,:,:)
      real(kind=dp_t), pointer::  tp(:,:,:,:)
      real(kind=dp_t), pointer::  np(:,:,:,:)
      real(kind=dp_t), pointer:: s0p(:,:,:,:)
!
      real(kind=dp_t), pointer:: sop(:,:,:,:)
      real(kind=dp_t), pointer:: snp(:,:,:,:)
      real(kind=dp_t), pointer:: sepx(:,:,:,:)
      real(kind=dp_t), pointer:: sepy(:,:,:,:)
      real(kind=dp_t), pointer:: sepz(:,:,:,:)
!
      real(dp_t)     :: mult
      real(dp_t)     :: smin,smax

      integer :: velpred
      integer :: lo(uold%dim),hi(uold%dim)
      integer :: i,n,bc_comp,dm,ng_cell,nr
      logical :: is_vel
      type(box)       :: domain
      integer         :: domlo(uold%dim), domhi(uold%dim)

      nr = size(s0_old,dim=1)
      allocate(s0_edge_old(0:nr,nscal))
      allocate(s0_edge_new(0:nr,nscal))

      call cell_to_edge_n(s0_old,s0_edge_old)
      call cell_to_edge_n(s0_new,s0_edge_new)

      domain = layout_get_pd(uold%la)
      domlo = lwb(domain)
      domhi = upb(domain)

      velpred = 0

      ng_cell = sold%ng
      dm      = sold%dim

      is_vel = .false.

      call build(scal_force, ext_scal_force%la, nscal, 1)
      call setval(scal_force,ZERO)

      call build(s0_cart   , ext_scal_force%la, nscal, 1)
      call setval(s0_cart   ,ZERO,all=.true.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create scalar source term at time n for (rho X)_i and (rho H).  
!     The source term for (rho X) is zero.
!     The source term for (rho h) has only the w dp0/dr term.

!     The call to modify_scal_force is used to add those advective terms 
!     that appear as forces when we write it in convective/perturbational form.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Define s0_cart 
      if (spherical .eq. 1) then
        do i = 1, s0_cart%nboxes
           if ( multifab_remote(s0_cart, i) ) cycle
           s0p => dataptr(s0_cart, i)
           lo =  lwb(get_box(s0_cart, i))
           hi =  upb(get_box(s0_cart, i))
           do n = spec_comp,spec_comp+nspec-1
             call fill_3d_data(s0p(:,:,:,n),s0_old(:,n),lo,hi,dx,1)
           end do
           n = rhoh_comp
           call fill_3d_data(s0p(:,:,:,n),s0_old(:,n),lo,hi,dx,1)
        end do
        call multifab_fill_boundary_c(s0_cart,spec_comp,nspec)
        call multifab_fill_boundary_c(s0_cart,rhoh_comp,1)
      end if

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
               call modify_scal_force_2d(fp(:,:,1,n),sop(:,:,1,n), lo, hi, &
                                         ng_cell,ump(:,:,1,1),vmp(:,:,1,1), &
                                         s0_old(:,n), s0_edge_old(:,n), w0, dx)
            end do

            if (use_temp_in_mkflux) then
              n = temp_comp
              tp => dataptr(thermal, i)

              ! This can be uncommented if you wish to compute T           
              ! call makeTfromRhoH_2d(sop(:,:,1,:), lo, hi, ng_cell, s0_old(:,temp_comp))

              call mktempforce_2d(fp(:,:,1,n), sop(:,:,1,:), tp(:,:,1,1), lo, hi, ng_cell, p0_old)

            else
              n = rhoh_comp
              call mkrhohforce_2d(fp(:,:,1,n), vmp(:,:,1,1), lo, hi, &
                                  p0_old, p0_old)

              call modify_scal_force_2d(fp(:,:,1,n),sop(:,:,1,n), lo, hi, &
                                        ng_cell, ump(:,:,1,1),vmp(:,:,1,1), &
                                        s0_old(:,rhoh_comp),s0_edge_old(:,rhoh_comp), w0,dx)            
            end if

         case(3)
            wmp  => dataptr(umac(3), i)

            if (spherical .eq. 1) then
               s0p => dataptr(s0_cart, i)
               do n = spec_comp,spec_comp+nspec-1
                  call modify_scal_force_3d_sphr(fp(:,:,:,n),sop(:,:,:,n), &
                                                 lo,hi,domlo,domhi,ng_cell, &
                                                 ump(:,:,:,1),vmp(:,:,:,1), &
                                                 wmp(:,:,:,1),s0p(:,:,:,n), &
                                                 w0,dx)
               end do
                
               if (use_temp_in_mkflux) then

                 n = temp_comp
                 tp => dataptr(thermal, i)

                 ! This can be uncommented if you wish to compute T      
                 ! call makeTfromRhoH_3d(sop(:,:,:,:), lo, hi, ng_cell, s0_old(:,temp_comp))
   
                 call mktempforce_3d_sphr(fp(:,:,:,n), sop(:,:,:,:), tp(:,:,:,1), lo, hi, ng_cell, p0_old, dx)

               else 

                  n = rhoh_comp
                  np => dataptr(normal, i)
                  call  mkrhohforce_3d_sphr(fp(:,:,:,n), &
                                            ump(:,:,:,1), vmp(:,:,:,1), &
                                            wmp(:,:,:,1), lo, hi, &
                                            dx, np(:,:,:,:), p0_old, p0_old)

                  call modify_scal_force_3d_sphr(fp(:,:,:,n),sop(:,:,:,n),lo,hi, &
                                                 domlo,domhi,ng_cell,&
                                                 ump(:,:,:,1),vmp(:,:,:,1), &
                                                 wmp(:,:,:,1), s0p(:,:,:,n),w0,dx)
               end if

            else
               do n = spec_comp,spec_comp+nspec-1
                  call modify_scal_force_3d_cart(fp(:,:,:,n),sop(:,:,:,n), &
                                                 lo,hi,ng_cell,ump(:,:,:,1), &
                                                 vmp(:,:,:,1),wmp(:,:,:,1), &
                                                 s0_old(:,n),s0_edge_old(:,n),w0,dx)
               end do

               if (use_temp_in_mkflux) then

                 n = temp_comp
                 tp => dataptr(thermal, i)

                 ! This can be uncommented if you wish to compute T      
                 ! call makeTfromRhoH_3d(sop(:,:,:,:), lo, hi, ng_cell, s0_old(:,temp_comp))

                 call mktempforce_3d(fp(:,:,:,n), sop(:,:,:,:), tp(:,:,:,1), lo, hi, ng_cell, p0_old)

               else

                  n = rhoh_comp
                  call  mkrhohforce_3d(fp(:,:,:,n), wmp(:,:,:,1), lo, hi, &
                                       p0_old, p0_old)
   
                  call modify_scal_force_3d_cart(fp(:,:,:,n),sop(:,:,:,n),lo,hi, &
                                                 ng_cell,ump(:,:,:,1), &
                                                 vmp(:,:,:,1),wmp(:,:,:,1), &
                                                 s0_old(:,n),s0_edge_old(:,n),w0,dx)
               end if
            end if
         end select

         ! add to the rhoh component of force if NOT use_temp_in_mkflux
         if ( (.not. use_temp_in_mkflux) .and. use_thermal_diffusion) &
            call multifab_plus_plus_c(scal_force,rhoh_comp,thermal,1,1)

      end do

      ! Do this because we have just defined temperature in the valid region
      if (use_temp_in_mkflux) &
        call multifab_fill_boundary_c(sold,temp_comp,1)

      call multifab_fill_boundary(scal_force)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Add w0 to MAC velocities (trans velocities already have w0).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      mult = ONE
      call addw0(umac,w0,w0_cart_vec,dx,mult)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the edge states of (rho h)' and (rho X)_i.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (.not. use_temp_in_mkflux) &
        call put_in_pert_form(sold,s0_old,dx,rhoh_comp,    1,.true.)
      call put_in_pert_form(sold,s0_old,dx,spec_comp,nspec,.true.)

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
              if (use_temp_in_mkflux) then
                n = temp_comp
              else
                n = rhoh_comp
              end if
              bc_comp = dm+n

              call mkflux_2d(sop(:,:,1,:), uop(:,:,1,:), &
                             sepx(:,:,1,:), sepy(:,:,1,:), &
                             ump(:,:,1,1), vmp(:,:,1,1), &
                             utp(:,:,1,1), vtp(:,:,1,1), fp(:,:,1,:), w0, &
                             lo, dx, dt, is_vel, &
                             the_bc_level%phys_bc_level_array(i,:,:), &
                             the_bc_level%adv_bc_level_array(i,:,:,bc_comp:), &
                             velpred, ng_cell, n)

              do n = spec_comp,spec_comp+nspec-1
                bc_comp = dm+n
                call mkflux_2d(sop(:,:,1,:), uop(:,:,1,:), &
                               sepx(:,:,1,:), sepy(:,:,1,:), &
                               ump(:,:,1,1), vmp(:,:,1,1), &
                               utp(:,:,1,1), vtp(:,:,1,1), fp(:,:,1,:), w0, &
                               lo, dx, dt, is_vel, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%adv_bc_level_array(i,:,:,bc_comp:), &
                               velpred, ng_cell, n)
              end do

              if (use_temp_in_mkflux) &
                call makeRhoHfromT_2d(sepx(:,:,1,:), sepy(:,:,1,:), &
                                      s0_old, s0_edge_old, &
                                      s0_new, s0_edge_new, lo, hi)

            case (3)
              wmp  => dataptr(  umac(3), i)
              wtp  => dataptr(utrans(3), i)
              sepz => dataptr( sedge(3), i)
               w0p => dataptr(w0_cart_vec, i)

                if (use_temp_in_mkflux) then
                  n = temp_comp
                else
                  n = rhoh_comp
                end if
                bc_comp = dm+n
                call mkflux_3d(sop(:,:,:,:), uop(:,:,:,:), &
                               sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                               ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                               utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), fp(:,:,:,:), w0, w0p(:,:,:,:), &
                               lo, dx, dt, is_vel, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%adv_bc_level_array(i,:,:,bc_comp:), &
                               velpred, ng_cell, n)

              do n = spec_comp,spec_comp+nspec-1
                bc_comp = dm+n
                call mkflux_3d(sop(:,:,:,:), uop(:,:,:,:), &
                               sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                               ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                               utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), fp(:,:,:,:), w0, w0p(:,:,:,:), &
                               lo, dx, dt, is_vel, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%adv_bc_level_array(i,:,:,bc_comp:), &
                               velpred, ng_cell, n)
              end do

              if (use_temp_in_mkflux) &
                call makeRhoHfromT_3d(sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                      s0_old, s0_edge_old, &
                                      s0_new, s0_edge_new, lo, hi)

         end select
      end do
      if (.not. use_temp_in_mkflux) &
        call put_in_pert_form(sold,s0_old,dx,rhoh_comp,    1,.false.)
      call put_in_pert_form(sold,s0_old,dx,spec_comp,nspec,.false.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the edge states of tracers.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (ntrac .ge. 1) then
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
                               utp(:,:,1,1), vtp(:,:,1,1), fp(:,:,1,:), w0, &
                               lo, dx, dt, is_vel, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%adv_bc_level_array(i,:,:,bc_comp:), &
                               velpred, ng_cell, n)
              end do
            case (3)
              wmp  => dataptr(  umac(3), i)
              wtp  => dataptr(utrans(3), i)
              sepz => dataptr( sedge(3), i)
              w0p  => dataptr(w0_cart_vec, i)
              do n = trac_comp,trac_comp+ntrac-1
                bc_comp = dm+n
                call mkflux_3d(sop(:,:,:,:), uop(:,:,:,:), &
                               sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                               ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                               utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), fp(:,:,:,:), w0, w0p(:,:,:,:), &
                               lo, dx, dt, is_vel, &
                               the_bc_level%phys_bc_level_array(i,:,:), &
                               the_bc_level%adv_bc_level_array(i,:,:,bc_comp:), &
                               velpred, ng_cell, n)
              end do
         end select
      end do
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Subtract w0 from MAC velocities.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      mult = -ONE
      call addw0(umac,w0,w0_cart_vec,dx,mult)

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
            call update_scal_2d(which_step, spec_comp, spec_comp+nspec-1, &
                                sop(:,:,1,:), snp(:,:,1,:), &
                                ump(:,:,1,1), vmp(:,:,1,1), w0, &
                                sepx(:,:,1,:), sepy(:,:,1,:), fp(:,:,1,:), &
                                s0_old(:,:), s0_edge_old(:,:), s0_new(:,:), &
                                lo, hi, ng_cell, dx, dt)
         case (3)
            wmp => dataptr(umac(3), i)
            sepz => dataptr(sedge(3), i)
            w0p => dataptr(w0_cart_vec, i)
            if (spherical .eq. 0) then
               call update_scal_3d_cart(which_step, spec_comp, spec_comp+nspec-1, &
                                        sop(:,:,:,:), snp(:,:,:,:), &
                                        ump(:,:,:,1), vmp(:,:,:,1), &
                                        wmp(:,:,:,1), w0, &
                                        sepx(:,:,:,:), sepy(:,:,:,:), &
                                        sepz(:,:,:,:), fp(:,:,:,:), &
                                        s0_old(:,:), s0_edge_old(:,:), s0_new(:,:), &
                                        lo, hi, ng_cell, dx, dt)
            else
               s0p => dataptr(s0_cart, i)
               call update_scal_3d_sphr(which_step, spec_comp, spec_comp+nspec-1, &
                                        sop(:,:,:,:), snp(:,:,:,:), &
                                        ump(:,:,:,1), vmp(:,:,:,1), &
                                        wmp(:,:,:,1), w0p(:,:,:,:), &
                                        sepx(:,:,:,:), sepy(:,:,:,:), &
                                        sepz(:,:,:,:), fp(:,:,:,:), &
                                        s0_old(:,:), s0_new(:,:), &
                                        s0p(:,:,:,:), lo, hi, domlo, domhi, &
                                        ng_cell, dx, dt)
            end if
         end select
      end do
            
      ! Make sure we do this before the calls to setbc.
      call multifab_fill_boundary(snew)

      do i = 1, snew%nboxes
         if ( multifab_remote(snew, i) ) cycle
         snp => dataptr(snew, i)
         select case (dm)
         case (2)
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

      if (verbose .ge. 1) then
         do n = spec_comp,spec_comp+nspec-1
            call multifab_div_div_c(snew,n,snew,rho_comp,1)

            smin = multifab_min_c(snew,n) 
            smax = multifab_max_c(snew,n)

            if (parallel_IOProcessor()) &
                 write(6,2002) spec_names(n-spec_comp+1), smin,smax
            call multifab_mult_mult_c(snew,n,snew,rho_comp,1)
         end do

         smin = multifab_min_c(snew,rho_comp) 
         smax = multifab_max_c(snew,rho_comp)

         if (parallel_IOProcessor()) &
              write(6,2000) smin,smax
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     2) Update tracers with convective differencing.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Define s0_cart 
      if (spherical .eq. 1) then
        do i = 1, sold%nboxes
           if ( multifab_remote(s0_cart, i) ) cycle
           s0p => dataptr(s0_cart, i)
           lo =  lwb(get_box(s0_cart, i))
           hi =  upb(get_box(s0_cart, i))
           do n = trac_comp,trac_comp+ntrac-1
             call fill_3d_data(s0p(:,:,:,n),s0_old(:,n),lo,hi,dx,1)
           end do
        end do
        do n = trac_comp,trac_comp+ntrac-1
           call multifab_fill_boundary_c(s0_cart,n,1)
        end do
      end if

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
            call update_scal_2d(which_step, trac_comp,trac_comp+ntrac-1, &
                                sop(:,:,1,:), snp(:,:,1,:), &
                                ump(:,:,1,1), vmp(:,:,1,1), w0, &
                                sepx(:,:,1,:), sepy(:,:,1,:), fp(:,:,1,:), &
                                s0_old(:,:), s0_edge_old(:,:), s0_new(:,:), &
                                lo, hi, ng_cell, dx, dt)

         case (3)
            wmp => dataptr(umac(3), i)
            sepz => dataptr(sedge(3), i)
            w0p => dataptr(w0_cart_vec, i)
            if (spherical .eq. 0) then
               call update_scal_3d_cart(which_step, trac_comp,trac_comp+ntrac-1, &
                                        sop(:,:,:,:), snp(:,:,:,:), &
                                        ump(:,:,:,1), vmp(:,:,:,1), &
                                        wmp(:,:,:,1), w0, &
                                        sepx(:,:,:,:), sepy(:,:,:,:), &
                                        sepz(:,:,:,:),   fp(:,:,:,:), &
                                        s0_old(:,:), s0_edge_old(:,:), s0_new(:,:), &
                                        lo, hi, ng_cell, dx, dt)
              else
                s0p => dataptr(s0_cart, i)
                call update_scal_3d_sphr(which_step, trac_comp,trac_comp+ntrac-1, &
                                         sop(:,:,:,:), snp(:,:,:,:), &
                                         ump(:,:,:,1), vmp(:,:,:,1), &
                                         wmp(:,:,:,1), w0p(:,:,:,:), &
                                         sepx(:,:,:,:), sepy(:,:,:,:), &
                                         sepz(:,:,:,:), fp(:,:,:,:), &
                                         s0_old(:,:), s0_new(:,:), &
                                         s0p(:,:,:,:), &
                                         lo, hi, domlo, domhi, ng_cell, dx, dt)
              end if
         end select
      end do

      ! Make sure we do this before the calls to setbc.
      call multifab_fill_boundary(snew)

      do i = 1, snew%nboxes
         if ( multifab_remote(snew, i) ) cycle
         snp => dataptr(snew, i)
         select case (dm)
         case (2)
            do n = trac_comp,trac_comp+ntrac-1
               bc_comp = dm+n
               call setbc_2d(snp(:,:,1,n), lo, ng_cell, &
                             the_bc_level%adv_bc_level_array(i,:,:,bc_comp), &
                             dx,bc_comp)
            end do
         case (3)
              do n = trac_comp,trac_comp+ntrac-1
                bc_comp = dm+n
                call setbc_3d(snp(:,:,:,n), lo, ng_cell, &
                              the_bc_level%adv_bc_level_array(i,:,:,bc_comp), &
                              dx,bc_comp)
              end do
         end select
      end do

      if (verbose .eq. 1) then
        smin = multifab_min_c(snew,trac_comp) 
        smax = multifab_max_c(snew,trac_comp)
        if (parallel_IOProcessor()) &
          write(6,2003) smin,smax
      end if

      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     1) Create (rhoh)' force at time n+1/2.
!          (NOTE: we don't worry about filling ghost cells of the scal_force
!                 because we only need them in valid regions...)     
!     2) Update (rho h)' with conservative differencing.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      n = rhoh_comp
      bc_comp = dm+n

      ! Define s0_cart 
      if (spherical .eq. 1) then
        do i = 1, sold%nboxes
           if ( multifab_remote(s0_cart, i) ) cycle
           s0p => dataptr(s0_cart, i)
           lo =  lwb(get_box(s0_cart, i))
           hi =  upb(get_box(s0_cart, i))
           n = rhoh_comp
           call fill_3d_data(s0p(:,:,:,n),s0_old(:,n),lo,hi,dx,1)
        end do
        call multifab_fill_boundary_c(s0_cart,rhoh_comp,1)
      end if

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
              call  mkrhohforce_2d(fp(:,:,1,n), vmp(:,:,1,1), lo, hi, &
                                   p0_old, p0_new)

              call update_scal_2d(which_step, rhoh_comp, rhoh_comp, &
                             sop(:,:,1,:), snp(:,:,1,:), &
                             ump(:,:,1,1), vmp(:,:,1,1), w0, &
                             sepx(:,:,1,:), sepy(:,:,1,:), fp(:,:,1,:), &
                             s0_old(:,:), s0_edge_old(:,:), s0_new(:,:), &
                             lo, hi, ng_cell, dx, dt)
              
            case(3)
              wmp  => dataptr(umac(3), i)
              w0p => dataptr(w0_cart_vec, i)
              sepz => dataptr(sedge(3), i)

              if (spherical .eq. 1) then

                np => dataptr(normal, i)
                call mkrhohforce_3d_sphr(fp(:,:,:,n), &
                                         ump(:,:,:,1), vmp(:,:,:,1), &
                                         wmp(:,:,:,1), lo, hi, &
                                         dx, np(:,:,:,:), p0_old, p0_new)
                s0p => dataptr(s0_cart, i)
                call update_scal_3d_sphr(which_step, rhoh_comp, rhoh_comp, &
                                         sop(:,:,:,:), snp(:,:,:,:), &
                                         ump(:,:,:,1), vmp(:,:,:,1), &
                                         wmp(:,:,:,1), w0p(:,:,:,:), &
                                         sepx(:,:,:,:), sepy(:,:,:,:), &
                                         sepz(:,:,:,:), fp(:,:,:,:), &
                                         s0_old(:,:), s0_new(:,:), &
                                         s0p(:,:,:,:), lo, hi, &
                                         domlo, domhi, ng_cell, dx, dt)

              else

                call mkrhohforce_3d(fp(:,:,:,n), wmp(:,:,:,1), lo, hi, &
                                    p0_old, p0_new)

                call update_scal_3d_cart(which_step, rhoh_comp, rhoh_comp, &
                                         sop(:,:,:,:), snp(:,:,:,:), &
                                         ump(:,:,:,1), vmp(:,:,:,1), &
                                         wmp(:,:,:,1), w0, &
                                         sepx(:,:,:,:), sepy(:,:,:,:), &
                                         sepz(:,:,:,:), fp(:,:,:,:), &
                                         s0_old(:,:), s0_edge_old(:,:), s0_new(:,:), &
                                         lo, hi, ng_cell, dx, dt)
               end if
         end select
      end do

      if(.not. use_thermal_diffusion) then
         ! compute updated temperature
         do i=1,snew%nboxes
            if (multifab_remote(snew,i)) cycle
            snp => dataptr(snew,i)
            lo = lwb(get_box(snew, i))
            hi = upb(get_box(snew, i))
            select case (dm)
            case (2)
               call makeTfromRhoH_2d(snp(:,:,1,:), lo, hi, 3, s0_new(:,temp_comp))
            case (3)
               call makeTfromRhoH_3d(snp(:,:,:,:), lo, hi, 3, s0_new(:,temp_comp))
            end select
         end do
      endif

      ! Make sure we do this before the calls to setbc.
      call multifab_fill_boundary(snew)

      do i = 1, snew%nboxes
         if ( multifab_remote(snew, i) ) cycle
         snp => dataptr(snew, i)
         select case (dm)
            case (2)
               do n=rho_comp,rho_comp+nscal-1
                  bc_comp = dm+n
                  call setbc_2d(snp(:,:,1,n), lo, ng_cell, &
                                the_bc_level%adv_bc_level_array(i,:,:,bc_comp), &
                                dx,bc_comp)
               enddo
            case (3)
               do n=rho_comp,rho_comp+nscal-1
                  bc_comp = dm+n
                  call setbc_3d(snp(:,:,:,n), lo, ng_cell, & 
                                the_bc_level%adv_bc_level_array(i,:,:,bc_comp), &
                                dx,bc_comp)
               enddo
         end select
      end do

      if (verbose .eq. 1) then
        smin = multifab_min_c(snew,rhoh_comp) 
        smax = multifab_max_c(snew,rhoh_comp)
        if (parallel_IOProcessor()) then
          write(6,2001) smin,smax
          write(6,2004) 
        end if
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Call fill_boundary for all components of snew
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call multifab_fill_boundary(snew)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call multifab_destroy(scal_force)
      call multifab_destroy(s0_cart)

      deallocate(s0_edge_old,s0_edge_new)

2000  format('... new min/max : density           ',e17.10,2x,e17.10)
2001  format('... new min/max : rho * H           ',e17.10,2x,e17.10)
2002  format('... new min/max : ',a16,2x,e17.10,2x,e17.10)
2003  format('... new min/max : tracer            ',e17.10,2x,e17.10)
2004  format(' ')

   end subroutine scalar_advance

   subroutine modify_scal_force_2d(force,s,lo,hi,ng,umac,vmac,base,base_edge,w0,dx)

    ! When we write the scalar equation in perturbational and convective
    ! form, the terms other than s'_t + U.grad s' act as source terms.  Add
    ! them to the forces here.

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(  out) :: force(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) ::     s(lo(1)-ng:,lo(2)-ng:)
    real(kind=dp_t), intent(in   ) ::  umac(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) ::  vmac(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(in   ) ::  base(0:), base_edge(0:)
    real(kind=dp_t), intent(in   ) ::    w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    integer :: i,j,nr
    real(kind=dp_t) :: divu,divbaseu

    nr = size(base,dim=1)

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
           divu =  (umac(i+1,j) - umac(i,j)) / dx(1) &
                +( (vmac(i,j+1) + w0(j+1)) &
                  -(vmac(i,j)   + w0(j)  ) ) / dx(2)
           divbaseu = base(j)*(umac(i+1,j) - umac(i,j))/dx(1) &
                             +(vmac(i,j+1) * base_edge(j+1) &
                             - vmac(i,j  ) * base_edge(j  ) )/ dx(2)
           force(i,j) = force(i,j) - (s(i,j)-base(j))*divu - divbaseu
       end do
     end do
     
   end subroutine modify_scal_force_2d

   subroutine modify_scal_force_3d_cart(force,s,lo,hi,ng,umac,vmac,wmac,base,base_edge,w0,dx)

    ! When we write the scalar equation in perturbational and convective
    ! form, the terms other than s'_t + U.grad s' act as source terms.  Add
    ! them to the forces here.

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(  out) :: force(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(kind=dp_t), intent(in   ) ::  umac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::  vmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) ::  wmac(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: base(0:), base_edge(0:)
    real(kind=dp_t), intent(in   ) ::   w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    integer :: i,j,k
    real(kind=dp_t) :: divu,divbaseu

    do k = lo(3),hi(3)
     do j = lo(2),hi(2)
       do i = lo(1),hi(1)
           divu = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) &
                 +(vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) &
                 +(wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)
           divbaseu = base(k)*( (umac(i+1,j,k) - umac(i,j,k))/dx(1) &
                               +(vmac(i,j+1,k) - vmac(i,j,k))/dx(2) ) &
                               +(wmac(i,j,k+1) * base_edge(k+1) &
                               - wmac(i,j,k  ) * base_edge(k  ))/ dx(3)
           divu = divu + (w0(k+1)-w0(k))/dx(3)
           force(i,j,k) = force(i,j,k) - (s(i,j,k)-base(k))*divu - divbaseu
       end do
     end do
    end do
     
   end subroutine modify_scal_force_3d_cart

   subroutine modify_scal_force_3d_sphr(force,s,lo,hi,domlo,domhi,ng, &
                                        umac,vmac,wmac,base_cart,w0,dx)

    ! When we write the scalar equation in perturbational and convective
    ! form, the terms other than s'_t + U.grad s' act as source terms.  Add
    ! them to the forces here.

    integer        , intent(in   ) :: lo(:),hi(:),domlo(:),domhi(:),ng
    real(kind=dp_t), intent(  out) :: force(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real(kind=dp_t), intent(in   ) ::     s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real(kind=dp_t), intent(in   ) ::  umac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real(kind=dp_t), intent(in   ) ::  vmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)
    real(kind=dp_t), intent(in   ) ::  wmac(lo(1)- 1:,lo(2)- 1:,lo(3)- 1:)

    real(kind=dp_t), intent(in   ) :: base_cart(lo(1)-1:,lo(2)-1:,lo(3)-1:)
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    
    ! Local variables
    integer :: i,j,k,nr
    real(kind=dp_t) :: divumac,divbaseu
    real(kind=dp_t) :: base_xlo,base_xhi
    real(kind=dp_t) :: base_ylo,base_yhi
    real(kind=dp_t) :: base_zlo,base_zhi

    real(kind=dp_t), allocatable :: divu(:),divu_cart(:,:,:)

    nr = size(w0,dim=1)-1
    allocate(divu(0:nr-1))
    allocate(divu_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))

    do k = 0,nr-1
      divu(k) = (zl(k+1)**2 * w0(k+1)- zl(k)**2 * w0(k))/(dr*z(k)**2)
    end do
    call fill_3d_data(divu_cart,divu,lo,hi,dx,0)

    do k = lo(3),hi(3)
        do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           divumac = (umac(i+1,j,k) - umac(i,j,k)) / dx(1) &
                    +(vmac(i,j+1,k) - vmac(i,j,k)) / dx(2) &
                    +(wmac(i,j,k+1) - wmac(i,j,k)) / dx(3)

           if (i.lt.domhi(1)) then
             base_xhi = HALF * (base_cart(i,j,k) + base_cart(i+1,j,k))
           else
             base_xhi = base_cart(i,j,k)
           end if
           if (i.gt.domlo(1)) then
             base_xlo = HALF * (base_cart(i,j,k) + base_cart(i-1,j,k))
           else
             base_xlo = base_cart(i,j,k)
           end if
           if (j.lt.domhi(2)) then
             base_yhi = HALF * (base_cart(i,j,k) + base_cart(i,j+1,k))
           else
             base_yhi = base_cart(i,j,k)
           end if
           if (j.gt.domlo(2)) then
             base_ylo = HALF * (base_cart(i,j,k) + base_cart(i,j-1,k))
           else
             base_ylo = base_cart(i,j,k)
           end if
           if (k.lt.domhi(3)) then
             base_zhi = HALF * (base_cart(i,j,k) + base_cart(i,j,k+1))
           else
             base_zhi = base_cart(i,j,k)
           end if
           if (k.gt.domlo(3)) then
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
     
   end subroutine modify_scal_force_3d_sphr

end module scalar_advance_module
