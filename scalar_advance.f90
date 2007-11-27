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
  use ml_restriction_module
  use multifab_fill_ghost_module

  implicit none
  
contains

  subroutine scalar_advance(nlevs, mla, which_step, uold, sold, snew, thermal, &
                            umac, w0, w0_cart_vec, eta, sedge, utrans, &
                            ext_scal_force, normal, &
                            s0_old, s0_new , &
                            p0_old, p0_new, &
                            dx, dt, the_bc_level, &
                            verbose)
 
    integer        , intent(in   ) :: nlevs
    type(ml_layout), intent(inout) :: mla
    integer        , intent(in   ) :: which_step
    type(multifab) , intent(inout) :: uold(:)
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: thermal(:)
    type(multifab) , intent(inout) :: umac(:,:)

    real(kind=dp_t), intent(inout) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0_cart_vec(:)
    real(kind=dp_t), intent(inout) :: eta(0:,:)

    type(multifab) , intent(inout) :: sedge(:,:)
    type(multifab) , intent(inout) :: utrans(:,:)
    type(multifab) , intent(inout) :: ext_scal_force(:)
    type(multifab) , intent(in   ) :: normal(:)

    real(kind=dp_t), intent(in   ) :: s0_old(0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(0:,:)
    real(kind=dp_t), intent(in   ) :: p0_old(0:)
    real(kind=dp_t), intent(in   ) :: p0_new(0:)

    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    integer        , intent(in   ) :: verbose

    ! local
    type(multifab), allocatable :: scal_force(:)
    type(multifab), allocatable :: s0_old_cart(:)
    type(multifab), allocatable :: s0_new_cart(:)
    
    real(kind=dp_t), allocatable :: s0_edge_old(:,:)
    real(kind=dp_t), allocatable :: s0_edge_new(:,:)

    real(kind=dp_t), pointer :: uop(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: utp(:,:,:,:)
    real(kind=dp_t), pointer :: vtp(:,:,:,:)
    real(kind=dp_t), pointer :: wtp(:,:,:,:)
    real(kind=dp_t), pointer :: w0p(:,:,:,:)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: dp(:,:,:,:)
    real(kind=dp_t), pointer :: tp(:,:,:,:)
    real(kind=dp_t), pointer :: np(:,:,:,:)
    real(kind=dp_t), pointer :: s0op(:,:,:,:)
    real(kind=dp_t), pointer :: s0np(:,:,:,:)

    real(kind=dp_t), pointer :: sop(:,:,:,:)
    real(kind=dp_t), pointer :: snp(:,:,:,:)
    real(kind=dp_t), pointer :: sepx(:,:,:,:)
    real(kind=dp_t), pointer :: sepy(:,:,:,:)
    real(kind=dp_t), pointer :: sepz(:,:,:,:)

    real(dp_t) :: mult
    real(dp_t) :: smin,smax

    type(box) :: domain

    integer :: velpred
    integer :: lo(uold(1)%dim),hi(uold(1)%dim)
    integer :: i,comp,n,bc_comp,dm,ng_cell,nr
    logical :: is_vel
    integer :: domlo(uold(1)%dim),domhi(uold(1)%dim)
    
    allocate(scal_force(nlevs))
    allocate(s0_old_cart(nlevs))
    allocate(s0_new_cart(nlevs))

    nr = size(s0_old,dim=1)
    allocate(s0_edge_old(0:nr,nscal))
    allocate(s0_edge_new(0:nr,nscal))
    
    velpred = 0    
    is_vel = .false.

    call cell_to_edge_n(s0_old,s0_edge_old)
    call cell_to_edge_n(s0_new,s0_edge_new)

    do n = 1, nlevs
    
    domain = layout_get_pd(uold(n)%la)
    domlo = lwb(domain)
    domhi = upb(domain)
        
    ng_cell = sold(n)%ng
    dm      = sold(n)%dim

    call build(scal_force(n),  ext_scal_force(n)%la, nscal, 1)
    call build(s0_old_cart(n), ext_scal_force(n)%la, nscal, 1)
    call build(s0_new_cart(n), ext_scal_force(n)%la, nscal, 1)

    call setval(scal_force(n), ZERO,all=.true.)
    call setval(s0_old_cart(n),ZERO,all=.true.)
    call setval(s0_new_cart(n),ZERO,all=.true.)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create scalar source term at time n for (rho X)_i and (rho H).  
!     The source term for (rho X) is zero.
!     The source term for (rho h) has only the w dp0/dr term.

!     The call to modify_scal_force is used to add those advective terms 
!     that appear as forces when we write it in convective/perturbational form.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Define s0_old_cart and s0_new_cart
    if (spherical .eq. 1) then
       do i = 1, s0_old_cart(n)%nboxes
          if ( multifab_remote(s0_old_cart(n), i) ) cycle
          s0op => dataptr(s0_old_cart(n), i)
          s0np => dataptr(s0_new_cart(n), i)
          lo =  lwb(get_box(s0_old_cart(n), i))
          hi =  upb(get_box(s0_old_cart(n), i))
          do comp = spec_comp,spec_comp+nspec-1
             call fill_3d_data(s0op(:,:,:,comp),s0_old(:,comp),lo,hi,dx(n,:),1)
             call fill_3d_data(s0np(:,:,:,comp),s0_new(:,comp),lo,hi,dx(n,:),1)
          end do
          comp = rhoh_comp
          call fill_3d_data(s0op(:,:,:,comp),s0_old(:,comp),lo,hi,dx(n,:),1)
          call fill_3d_data(s0np(:,:,:,comp),s0_new(:,comp),lo,hi,dx(n,:),1)
       end do
       call multifab_fill_boundary_c(s0_old_cart(n),spec_comp,nspec)
       call multifab_fill_boundary_c(s0_old_cart(n),rhoh_comp,1)
       call multifab_fill_boundary_c(s0_new_cart(n),spec_comp,nspec)
       call multifab_fill_boundary_c(s0_new_cart(n),rhoh_comp,1)
    end if
    
    do i = 1, scal_force(n)%nboxes
       if ( multifab_remote(scal_force(n), i) ) cycle
       fp => dataptr(scal_force(n), i)
       sop => dataptr(sold(n) , i)
       ump => dataptr(umac(n,1), i)
       vmp => dataptr(umac(n,2), i)
       lo =  lwb(get_box(sold(n), i))
       hi =  upb(get_box(sold(n), i))
       select case (dm)
       case (2)
          
          do comp = spec_comp,spec_comp+nspec-1
             call modify_scal_force_2d(fp(:,:,1,comp),sop(:,:,1,comp), lo, hi, &
                                       ng_cell,ump(:,:,1,1),vmp(:,:,1,1), &
                                       s0_old(:,comp), s0_edge_old(:,comp), w0(n,:), dx(n,:))
          end do
          
          if (use_temp_in_mkflux) then
             comp = temp_comp
             tp => dataptr(thermal(n), i)
             
             ! This can be uncommented if you wish to compute T           
             ! call makeTfromRhoH_2d(sop(:,:,1,:), lo, hi, ng_cell, s0_old(:,temp_comp))
             
             call mktempforce_2d(fp(:,:,1,comp), sop(:,:,1,:), tp(:,:,1,1), lo, hi, &
                                 ng_cell, p0_old)
          else
             comp = rhoh_comp
             call mkrhohforce_2d(fp(:,:,1,comp), vmp(:,:,1,1), lo, hi, &
                                 p0_old, p0_old)
             
             call modify_scal_force_2d(fp(:,:,1,comp),sop(:,:,1,comp), lo, hi, &
                                       ng_cell, ump(:,:,1,1),vmp(:,:,1,1), &
                                       s0_old(:,rhoh_comp),s0_edge_old(:,rhoh_comp), &
                                       w0(n,:),dx(n,:))
          end if

       case(3)
          wmp  => dataptr(umac(n,3), i)
          
          if (spherical .eq. 1) then
             s0op => dataptr(s0_old_cart(n), i)
             do comp = spec_comp,spec_comp+nspec-1
                call modify_scal_force_3d_sphr(fp(:,:,:,comp),sop(:,:,:,comp), &
                                               lo,hi,domlo,domhi,ng_cell, &
                                               ump(:,:,:,1),vmp(:,:,:,1), &
                                               wmp(:,:,:,1),s0op(:,:,:,comp), &
                                               w0(n,:),dx(n,:))
             end do
             
             if (use_temp_in_mkflux) then
                comp = temp_comp
                tp => dataptr(thermal(n), i)
                
                ! This can be uncommented if you wish to compute T      
                ! call makeTfromRhoH_3d(sop(:,:,:,:), lo, hi, ng_cell, s0_old(:,temp_comp))
                
                call mktempforce_3d_sphr(fp(:,:,:,comp), sop(:,:,:,:), tp(:,:,:,1), &
                                         lo, hi, ng_cell, p0_old, dx(n,:))
               else 
                  comp = rhoh_comp
                  np => dataptr(normal(n), i)
                  call  mkrhohforce_3d_sphr(fp(:,:,:,comp), &
                                            ump(:,:,:,1), vmp(:,:,:,1), &
                                            wmp(:,:,:,1), lo, hi, &
                                            dx(n,:), np(:,:,:,:), p0_old, p0_old)

                  call modify_scal_force_3d_sphr(fp(:,:,:,comp),sop(:,:,:,comp),lo,hi, &
                                                 domlo,domhi,ng_cell,&
                                                 ump(:,:,:,1),vmp(:,:,:,1), &
                                                 wmp(:,:,:,1), s0op(:,:,:,comp),w0(n,:), &
                                                 dx(n,:))
               end if
               
            else
               do comp = spec_comp,spec_comp+nspec-1
                  call modify_scal_force_3d_cart(fp(:,:,:,comp),sop(:,:,:,comp), &
                                                 lo,hi,ng_cell,ump(:,:,:,1), &
                                                 vmp(:,:,:,1),wmp(:,:,:,1), &
                                                 s0_old(:,comp),s0_edge_old(:,comp), &
                                                 w0(n,:),dx(n,:))
               end do
               
               if (use_temp_in_mkflux) then

                  comp = temp_comp
                  tp => dataptr(thermal(n), i)
                  
                  ! This can be uncommented if you wish to compute T      
                  ! call makeTfromRhoH_3d(sop(:,:,:,:), lo, hi, ng_cell, s0_old(:,temp_comp))
                  
                  call mktempforce_3d(fp(:,:,:,comp), sop(:,:,:,:), tp(:,:,:,1), lo, hi, &
                                      ng_cell, p0_old)
                  
               else
                  
                  comp = rhoh_comp
                  call  mkrhohforce_3d(fp(:,:,:,comp), wmp(:,:,:,1), lo, hi, &
                                       p0_old, p0_old)
                  
                  call modify_scal_force_3d_cart(fp(:,:,:,comp),sop(:,:,:,comp),lo,hi, &
                                                 ng_cell,ump(:,:,:,1), &
                                                 vmp(:,:,:,1),wmp(:,:,:,1), &
                                                 s0_old(:,comp),s0_edge_old(:,comp), &
                                                 w0(n,:),dx(n,:))
               end if
            end if
         end select
         
         ! add to the rhoh component of force if NOT use_temp_in_mkflux
         if ( (.not. use_temp_in_mkflux) .and. use_thermal_diffusion) &
              call multifab_plus_plus_c(scal_force(n),rhoh_comp,thermal(n),1,1)
         
      end do
      
      ! Do this because we have just defined temperature in the valid region
      if (use_temp_in_mkflux) &
           call multifab_fill_boundary_c(sold(n),temp_comp,1)
      
      call multifab_fill_boundary(scal_force(n))
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Add w0 to MAC velocities (trans velocities already have w0).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      mult = ONE
      call addw0(umac(n,:),w0(n,:),w0_cart_vec(n),dx(n,:),mult)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the edge states of (rho h)' and (rho X)_i.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (.not. use_temp_in_mkflux) &
           call put_in_pert_form(sold(n),s0_old,dx(n,:),rhoh_comp,    1,.true.)
      
      call put_in_pert_form(sold(n),s0_old,dx(n,:),spec_comp,nspec,.true.)
      
      do i = 1, sold(n)%nboxes
         if ( multifab_remote(sold(n), i) ) cycle
         sop  => dataptr(sold(n), i)
         uop  => dataptr(uold(n), i)
         sepx => dataptr(sedge(n,1), i)
         sepy => dataptr(sedge(n,2), i)
         ump  => dataptr(umac(n,1), i)
         vmp  => dataptr(umac(n,2), i)
         utp  => dataptr(utrans(n,1), i)
         vtp  => dataptr(utrans(n,2), i)
         fp  => dataptr(scal_force(n) , i)
         lo =  lwb(get_box(uold(n), i))
         hi =  upb(get_box(uold(n), i))
         select case (dm)
         case (2)
            if (use_temp_in_mkflux) then
               comp = temp_comp
            else
               comp = rhoh_comp
            end if
            bc_comp = dm+comp
            
            call mkflux_2d(sop(:,:,1,:), uop(:,:,1,:), &
                           sepx(:,:,1,:), sepy(:,:,1,:), &
                           ump(:,:,1,1), vmp(:,:,1,1), &
                           utp(:,:,1,1), vtp(:,:,1,1), fp(:,:,1,:), w0(1,:), &
                           lo, dx(n,:), dt, is_vel, &
                           the_bc_level(n)%phys_bc_level_array(i,:,:), &
                           the_bc_level(n)%adv_bc_level_array(i,:,:,bc_comp:), &
                           velpred, ng_cell, comp)
            
            do comp = spec_comp,spec_comp+nspec-1
               bc_comp = dm+comp
               call mkflux_2d(sop(:,:,1,:), uop(:,:,1,:), &
                              sepx(:,:,1,:), sepy(:,:,1,:), &
                              ump(:,:,1,1), vmp(:,:,1,1), &
                              utp(:,:,1,1), vtp(:,:,1,1), fp(:,:,1,:), w0(1,:), &
                              lo, dx(n,:), dt, is_vel, &
                              the_bc_level(n)%phys_bc_level_array(i,:,:), &
                              the_bc_level(n)%adv_bc_level_array(i,:,:,bc_comp:), &
                              velpred, ng_cell, comp)
            end do

            if (use_temp_in_mkflux) &
                 call makeRhoHfromT_2d(sepx(:,:,1,:), sepy(:,:,1,:), &
                                       s0_old, s0_edge_old, &
                                       s0_new, s0_edge_new, lo, hi)
            
         case (3)
            wmp  => dataptr(  umac(n,3), i)
            wtp  => dataptr(utrans(n,3), i)
            sepz => dataptr( sedge(n,3), i)
            w0p => dataptr(w0_cart_vec(n), i)
            
            if (use_temp_in_mkflux) then
               comp = temp_comp
            else
               comp = rhoh_comp
            end if
            bc_comp = dm+comp
            call mkflux_3d(sop(:,:,:,:), uop(:,:,:,:), &
                           sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                           ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                           utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), fp(:,:,:,:), &
                           w0(1,:), w0p(:,:,:,:), &
                           lo, dx(n,:), dt, is_vel, &
                           the_bc_level(n)%phys_bc_level_array(i,:,:), &
                           the_bc_level(n)%adv_bc_level_array(i,:,:,bc_comp:), &
                           velpred, ng_cell, comp)
            
            do comp = spec_comp,spec_comp+nspec-1
               bc_comp = dm+comp
               call mkflux_3d(sop(:,:,:,:), uop(:,:,:,:), &
                              sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                              ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                              utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), fp(:,:,:,:), &
                              w0(1,:), w0p(:,:,:,:), &
                              lo, dx(n,:), dt, is_vel, &
                              the_bc_level(n)%phys_bc_level_array(i,:,:), &
                              the_bc_level(n)%adv_bc_level_array(i,:,:,bc_comp:), &
                              velpred, ng_cell, comp)
            end do
            
            if (use_temp_in_mkflux) &
                 call makeRhoHfromT_3d(sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                       s0_old, s0_edge_old, &
                                       s0_new, s0_edge_new, lo, hi)
            
         end select
      end do
      if (.not. use_temp_in_mkflux) &
           call put_in_pert_form(sold(n),s0_old,dx(n,:),rhoh_comp,    1,.false.)

      call put_in_pert_form(sold(n),s0_old,dx(n,:),spec_comp,nspec,.false.)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the edge states of tracers.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (ntrac .ge. 1) then
         do i = 1, sold(n)%nboxes
            if ( multifab_remote(sold(n), i) ) cycle
            sop  => dataptr(sold(n), i)
            uop  => dataptr(uold(n), i)
            sepx => dataptr(sedge(n,1), i)
            sepy => dataptr(sedge(n,2), i)
            ump  => dataptr(umac(n,1), i)
            vmp  => dataptr(umac(n,2), i)
            utp  => dataptr(utrans(n,1), i)
            vtp  => dataptr(utrans(n,2), i)
            fp  => dataptr(scal_force(n) , i)
            lo =  lwb(get_box(uold(n), i))
            hi =  upb(get_box(uold(n), i))
            select case (dm)
            case (2)
               do comp = trac_comp,trac_comp+ntrac-1
                  bc_comp = dm+comp
                  call mkflux_2d(sop(:,:,1,:), uop(:,:,1,:), &
                                 sepx(:,:,1,:), sepy(:,:,1,:), &
                                 ump(:,:,1,1), vmp(:,:,1,1), &
                                 utp(:,:,1,1), vtp(:,:,1,1), fp(:,:,1,:), w0(1,:), &
                                 lo, dx(n,:), dt, is_vel, &
                                 the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                 the_bc_level(n)%adv_bc_level_array(i,:,:,bc_comp:), &
                                 velpred, ng_cell, comp)
               end do
            case (3)
               wmp  => dataptr(  umac(n,3), i)
               wtp  => dataptr(utrans(n,3), i)
               sepz => dataptr( sedge(n,3), i)
               w0p  => dataptr(w0_cart_vec(n), i)
               do comp = trac_comp,trac_comp+ntrac-1
                  bc_comp = dm+comp
                  call mkflux_3d(sop(:,:,:,:), uop(:,:,:,:), &
                                 sepx(:,:,:,:), sepy(:,:,:,:), sepz(:,:,:,:), &
                                 ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                 utp(:,:,:,1), vtp(:,:,:,1), wtp(:,:,:,1), fp(:,:,:,:), &
                                 w0(1,:), w0p(:,:,:,:), &
                                 lo, dx(n,:), dt, is_vel, &
                                 the_bc_level(n)%phys_bc_level_array(i,:,:), &
                                 the_bc_level(n)%adv_bc_level_array(i,:,:,bc_comp:), &
                                 velpred, ng_cell, comp)
               end do
            end select
         end do
      end if
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Subtract w0 from MAC velocities.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      mult = -ONE
      call addw0(umac(n,:),w0(n,:),w0_cart_vec(n),dx(n,:),mult)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     1) Set force for (rho X)_i at time n+1/2 = 0.
!     2) Update (rho X)'_i with conservative differencing.
!     3) Define density as the sum of the (rho X)_i
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call setval(scal_force(n),ZERO)
      
      do i = 1, sold(n)%nboxes
         if ( multifab_remote(sold(n), i) ) cycle
         sop => dataptr(sold(n), i)
         snp => dataptr(snew(n), i)
         ump => dataptr(umac(n,1), i)
         vmp => dataptr(umac(n,2), i)
         sepx => dataptr(sedge(n,1), i)
         sepy => dataptr(sedge(n,2), i)
         fp => dataptr(scal_force(n) , i)
         lo =  lwb(get_box(sold(n), i))
         hi =  upb(get_box(sold(n), i))
         select case (dm)
         case (2)
            call update_scal_2d(which_step, spec_comp, spec_comp+nspec-1, &
                                sop(:,:,1,:), snp(:,:,1,:), &
                                ump(:,:,1,1), vmp(:,:,1,1), w0(1,:), eta, &
                                sepx(:,:,1,:), sepy(:,:,1,:), fp(:,:,1,:), &
                                s0_old(:,:), s0_edge_old(:,:), &
                                s0_new(:,:), s0_edge_new(:,:), &
                                lo, hi, ng_cell, dx(n,:), dt)
         case (3)
            wmp => dataptr(umac(n,3), i)
            sepz => dataptr(sedge(n,3), i)
            w0p => dataptr(w0_cart_vec(n), i)
            if (spherical .eq. 0) then
               call update_scal_3d_cart(which_step, spec_comp, spec_comp+nspec-1, &
                                        sop(:,:,:,:), snp(:,:,:,:), &
                                        ump(:,:,:,1), vmp(:,:,:,1), &
                                        wmp(:,:,:,1), w0(1,:), eta, &
                                        sepx(:,:,:,:), sepy(:,:,:,:), &
                                        sepz(:,:,:,:), fp(:,:,:,:), &
                                        s0_old(:,:), s0_edge_old(:,:), &
                                        s0_new(:,:), s0_edge_new(:,:), &
                                        lo, hi, ng_cell, dx(n,:), dt)
            else
               s0op => dataptr(s0_old_cart(n), i)
               s0np => dataptr(s0_new_cart(n), i)
               call update_scal_3d_sphr(which_step, spec_comp, spec_comp+nspec-1, &
                                        sop(:,:,:,:), snp(:,:,:,:), &
                                        ump(:,:,:,1), vmp(:,:,:,1), &
                                        wmp(:,:,:,1), w0p(:,:,:,:), &
                                        sepx(:,:,:,:), sepy(:,:,:,:), &
                                        sepz(:,:,:,:), fp(:,:,:,:), &
                                        s0_old(:,:), s0_new(:,:), &
                                        s0op(:,:,:,:), s0np(:,:,:,:), &
                                        lo, hi, domlo, domhi, ng_cell, dx(n,:), dt)
            end if
         end select
      end do
      
      ! Make sure we do this before the calls to setbc.
      call multifab_fill_boundary(snew(n))
      
      do i = 1, snew(n)%nboxes
         if ( multifab_remote(snew(n), i) ) cycle
         snp => dataptr(snew(n), i)
         select case (dm)
         case (2)
            do comp = spec_comp,spec_comp+nspec-1
               bc_comp = dm+comp
               call setbc_2d(snp(:,:,1,comp), lo, ng_cell, &
                             the_bc_level(n)%adv_bc_level_array(i,:,:,bc_comp),dx(n,:),bc_comp)
            end do
            ! Dont forget to call setbc for density also
            comp = rho_comp
            bc_comp = dm+comp
            call setbc_2d(snp(:,:,1,comp), lo, ng_cell, &
                          the_bc_level(n)%adv_bc_level_array(i,:,:,bc_comp),dx(n,:),bc_comp)
         case (3)
            do comp = spec_comp,spec_comp+nspec-1
               bc_comp = dm+comp
               call setbc_3d(snp(:,:,:,comp), lo, ng_cell, &
                             the_bc_level(n)%adv_bc_level_array(i,:,:,bc_comp),dx(n,:),bc_comp)
            end do
            ! Dont forget to call setbc for density also
            comp = rho_comp
            bc_comp = dm+comp
            call setbc_3d(snp(:,:,:,comp), lo, ng_cell, &
                          the_bc_level(n)%adv_bc_level_array(i,:,:,bc_comp),dx(n,:),bc_comp)
         end select
      end do
      
      if (verbose .ge. 1) then
         do comp = spec_comp,spec_comp+nspec-1
            call multifab_div_div_c(snew(n),comp,snew(n),rho_comp,1)
            
            smin = multifab_min_c(snew(n),comp) 
            smax = multifab_max_c(snew(n),comp)
            
            if (parallel_IOProcessor()) &
                 write(6,2002) spec_names(comp-spec_comp+1), smin,smax
            call multifab_mult_mult_c(snew(n),comp,snew(n),rho_comp,1)
         end do
         
         smin = multifab_min_c(snew(n),rho_comp) 
         smax = multifab_max_c(snew(n),rho_comp)
         
         if (parallel_IOProcessor()) &
              write(6,2000) smin,smax
      end if
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     2) Update tracers with convective differencing.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Define s0_old_cart and s0_new_cart
      if (spherical .eq. 1) then
         do i = 1, sold(n)%nboxes
            if ( multifab_remote(s0_old_cart(n), i) ) cycle
            s0op => dataptr(s0_old_cart(n), i)
            s0np => dataptr(s0_new_cart(n), i)
            lo =  lwb(get_box(s0_old_cart(n), i))
            hi =  upb(get_box(s0_old_cart(n), i))
            do comp = trac_comp,trac_comp+ntrac-1
               call fill_3d_data(s0op(:,:,:,comp),s0_old(:,comp),lo,hi,dx(n,:),1)
               call fill_3d_data(s0np(:,:,:,comp),s0_new(:,comp),lo,hi,dx(n,:),1)
            end do
         end do
         do comp = trac_comp,trac_comp+ntrac-1
            call multifab_fill_boundary_c(s0_old_cart(n),comp,1)
            call multifab_fill_boundary_c(s0_new_cart(n),comp,1)
         end do
      end if
      
      if (ntrac .ge. 1) then
         do i = 1, sold(n)%nboxes
            if ( multifab_remote(sold(n), i) ) cycle
            sop => dataptr(sold(n), i)
            snp => dataptr(snew(n), i)
            ump => dataptr(umac(n,1), i)
            vmp => dataptr(umac(n,2), i)
            sepx => dataptr(sedge(n,1), i)
            sepy => dataptr(sedge(n,2), i)
            fp => dataptr(scal_force(n) , i)
            lo =  lwb(get_box(sold(n), i))
            hi =  upb(get_box(sold(n), i))
            select case (dm)
            case (2)
               call update_scal_2d(which_step, trac_comp,trac_comp+ntrac-1, &
                                   sop(:,:,1,:), snp(:,:,1,:), &
                                   ump(:,:,1,1), vmp(:,:,1,1), w0(1,:), eta, &
                                   sepx(:,:,1,:), sepy(:,:,1,:), fp(:,:,1,:), &
                                   s0_old(:,:), s0_edge_old(:,:), &
                                   s0_new(:,:), s0_edge_new(:,:), &
                                   lo, hi, ng_cell, dx(n,:), dt)
               
            case (3)
               wmp => dataptr(umac(n,3), i)
               sepz => dataptr(sedge(n,3), i)
               w0p => dataptr(w0_cart_vec(n), i)
               if (spherical .eq. 0) then
                  call update_scal_3d_cart(which_step, trac_comp,trac_comp+ntrac-1, &
                                           sop(:,:,:,:), snp(:,:,:,:), &
                                           ump(:,:,:,1), vmp(:,:,:,1), &
                                           wmp(:,:,:,1), w0(1,:), eta, &
                                           sepx(:,:,:,:), sepy(:,:,:,:), &
                                           sepz(:,:,:,:),   fp(:,:,:,:), &
                                           s0_old(:,:), s0_edge_old(:,:), &
                                           s0_new(:,:), s0_edge_new(:,:), &
                                           lo, hi, ng_cell, dx(n,:), dt)
               else
                  s0op => dataptr(s0_old_cart(n), i)
                  s0np => dataptr(s0_new_cart(n), i)
                  call update_scal_3d_sphr(which_step, trac_comp,trac_comp+ntrac-1, &
                                           sop(:,:,:,:), snp(:,:,:,:), &
                                           ump(:,:,:,1), vmp(:,:,:,1), &
                                           wmp(:,:,:,1), w0p(:,:,:,:), &
                                           sepx(:,:,:,:), sepy(:,:,:,:), &
                                           sepz(:,:,:,:), fp(:,:,:,:), &
                                           s0_old(:,:), s0_new(:,:), &
                                           s0op(:,:,:,:), s0np(:,:,:,:), &
                                           lo, hi, domlo, domhi, ng_cell, dx(n,:), dt)
               end if
            end select
         end do
         
         ! Make sure we do this before the calls to setbc.
         call multifab_fill_boundary(snew(n))
         
         do i = 1, snew(n)%nboxes
            if ( multifab_remote(snew(n), i) ) cycle
            snp => dataptr(snew(n), i)
            select case (dm)
            case (2)
               do comp = trac_comp,trac_comp+ntrac-1
                  bc_comp = dm+comp
                  call setbc_2d(snp(:,:,1,comp), lo, ng_cell, &
                                the_bc_level(n)%adv_bc_level_array(i,:,:,bc_comp), &
                                dx(n,:),bc_comp)
               end do
            case (3)
               do comp = trac_comp,trac_comp+ntrac-1
                  bc_comp = dm+comp
                  call setbc_3d(snp(:,:,:,comp), lo, ng_cell, &
                                the_bc_level(n)%adv_bc_level_array(i,:,:,bc_comp), &
                                dx(n,:),bc_comp)
               end do
            end select
         end do
         
         if (verbose .eq. 1) then
            smin = multifab_min_c(snew(n),trac_comp) 
            smax = multifab_max_c(snew(n),trac_comp)
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

      comp = rhoh_comp
      bc_comp = dm+comp
      
      ! Define s0_old_cart and s0_new_cart
      if (spherical .eq. 1) then
         do i = 1, sold(n)%nboxes
            if ( multifab_remote(s0_old_cart(n), i) ) cycle
            s0op => dataptr(s0_old_cart(n), i)
            s0np => dataptr(s0_new_cart(n), i)
            lo =  lwb(get_box(s0_old_cart(n), i))
            hi =  upb(get_box(s0_old_cart(n), i))
            comp = rhoh_comp
            call fill_3d_data(s0op(:,:,:,comp),s0_old(:,comp),lo,hi,dx(n,:),1)
            call fill_3d_data(s0np(:,:,:,comp),s0_new(:,comp),lo,hi,dx(n,:),1)
         end do
         call multifab_fill_boundary_c(s0_old_cart(n),rhoh_comp,1)
         call multifab_fill_boundary_c(s0_new_cart(n),rhoh_comp,1)
      end if
      
      do i = 1, sold(n)%nboxes
         if ( multifab_remote(sold(n), i) ) cycle
         sop => dataptr(sold(n), i)
         snp => dataptr(snew(n), i)
         ump => dataptr(umac(n,1), i)
         vmp => dataptr(umac(n,2), i)
         sepx => dataptr(sedge(n,1), i)
         sepy => dataptr(sedge(n,2), i)
         fp => dataptr(scal_force(n) , i)
         lo =  lwb(get_box(sold(n), i))
         hi =  upb(get_box(sold(n), i))
         select case (dm)
         case (2)
            call  mkrhohforce_2d(fp(:,:,1,comp), vmp(:,:,1,1), lo, hi, &
                                 p0_old, p0_new)
            
            call update_scal_2d(which_step, rhoh_comp, rhoh_comp, &
                                sop(:,:,1,:), snp(:,:,1,:), &
                                ump(:,:,1,1), vmp(:,:,1,1), w0(1,:), eta, &
                                sepx(:,:,1,:), sepy(:,:,1,:), fp(:,:,1,:), &
                                s0_old(:,:), s0_edge_old(:,:), &
                                s0_new(:,:), s0_edge_new(:,:), &
                                lo, hi, ng_cell, dx(n,:), dt)
            
         case(3)
            wmp  => dataptr(umac(n,3), i)
            w0p => dataptr(w0_cart_vec(n), i)
            sepz => dataptr(sedge(n,3), i)
            
            if (spherical .eq. 0) then
               
               call mkrhohforce_3d(fp(:,:,:,comp), wmp(:,:,:,1), lo, hi, &
                                   p0_old, p0_new)

               call update_scal_3d_cart(which_step, rhoh_comp, rhoh_comp, &
                                        sop(:,:,:,:), snp(:,:,:,:), &
                                        ump(:,:,:,1), vmp(:,:,:,1), &
                                        wmp(:,:,:,1), w0(1,:), eta, &
                                        sepx(:,:,:,:), sepy(:,:,:,:), &
                                        sepz(:,:,:,:), fp(:,:,:,:), &
                                        s0_old(:,:), s0_edge_old(:,:), &
                                        s0_new(:,:), s0_edge_new(:,:), &
                                        lo, hi, ng_cell, dx(n,:), dt)
              
            else
               
               np => dataptr(normal(n), i)
               call mkrhohforce_3d_sphr(fp(:,:,:,comp), &
                                        ump(:,:,:,1), vmp(:,:,:,1), &
                                        wmp(:,:,:,1), lo, hi, &
                                        dx(n,:), np(:,:,:,:), p0_old, p0_new)
               s0op => dataptr(s0_old_cart(n), i)
               s0np => dataptr(s0_new_cart(n), i)
               call update_scal_3d_sphr(which_step, rhoh_comp, rhoh_comp, &
                                        sop(:,:,:,:), snp(:,:,:,:), &
                                        ump(:,:,:,1), vmp(:,:,:,1), &
                                        wmp(:,:,:,1), w0p(:,:,:,:), &
                                        sepx(:,:,:,:), sepy(:,:,:,:), &
                                        sepz(:,:,:,:), fp(:,:,:,:), &
                                        s0_old(:,:), s0_new(:,:), &
                                        s0op(:,:,:,:), s0np(:,:,:,:), &
                                        lo, hi, domlo, domhi, ng_cell, dx(n,:), dt)
               
            end if
         end select
      end do
      
      if(.not. use_thermal_diffusion) then
         ! compute updated temperature
         do i=1,snew(n)%nboxes
            if (multifab_remote(snew(n),i)) cycle
            snp => dataptr(snew(n),i)
            lo = lwb(get_box(snew(n), i))
            hi = upb(get_box(snew(n), i))
            select case (dm)
            case (2)
               call makeTfromRhoH_2d(snp(:,:,1,:), lo, hi, 3, s0_new(:,temp_comp))
            case (3)
               call makeTfromRhoH_3d(snp(:,:,:,:), lo, hi, 3, s0_new(:,temp_comp))
            end select
         end do
      endif
      
      ! Make sure we do this before the calls to setbc.
      call multifab_fill_boundary(snew(n))
      
      do i = 1, snew(n)%nboxes
         if ( multifab_remote(snew(n), i) ) cycle
         snp => dataptr(snew(n), i)
         select case (dm)
         case (2)
            do comp = rho_comp,rho_comp+nscal-1
               bc_comp = dm+comp
               call setbc_2d(snp(:,:,1,comp), lo, ng_cell, &
                             the_bc_level(n)%adv_bc_level_array(i,:,:,bc_comp), &
                             dx(n,:),bc_comp)
            enddo
         case (3)
            do comp = rho_comp,rho_comp+nscal-1
               bc_comp = dm+comp
               call setbc_3d(snp(:,:,:,comp), lo, ng_cell, & 
                             the_bc_level(n)%adv_bc_level_array(i,:,:,bc_comp), &
                             dx(n,:),bc_comp)
            enddo
         end select
      end do
      
      if (verbose .eq. 1) then
         smin = multifab_min_c(snew(n),rhoh_comp) 
         smax = multifab_max_c(snew(n),rhoh_comp)
         if (parallel_IOProcessor()) then
            write(6,2001) smin,smax
            write(6,2004) 
         end if
      end if
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Call fill_boundary for all components of snew
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      enddo ! do n = 1, nlevs

      do n = 2, nlevs
         call ml_cc_restriction(snew(n-1),snew(n),mla%mba%rr(n-1,:))

         domain = layout_get_pd(mla%la(n))
         call multifab_fill_ghost_cells(snew(n),snew(n-1),domain, &
                                        ng_cell,mla%mba%rr(n-1,:), &
                                        the_bc_level(n-1)%adv_bc_level_array(0,:,:,:), &
                                        1,dm+rho_comp,nscal)
      end do
      
      do n = 1,nlevs
         call multifab_destroy(scal_force(n))
         call multifab_destroy(s0_old_cart(n))
         call multifab_destroy(s0_new_cart(n))
      enddo

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
     
     integer :: i,j
     real(kind=dp_t) :: divu,divbaseu
     
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
 
