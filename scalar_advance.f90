module scalar_advance_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use mkflux_module
  use mkscalforce_module
  use update_scal_module
  use addw0_module
  use define_bc_module
  use multifab_physbc_module
  use fill_3d_module
  use pert_form_module
  use cell_to_edge_module
  use rhoh_vs_t_module
  use variables
  use geometry
  use network
  use probin_module, ONLY: use_temp_in_mkflux, use_thermal_diffusion, evolve_base_state
  use ml_restriction_module
  use ml_layout_module
  use multifab_fill_ghost_module

  implicit none

  private
  public :: scalar_advance

contains

  subroutine scalar_advance(nlevs,mla,which_step,uold,sold,snew,thermal,umac,w0, &
                            w0_cart_vec,eta,sedge,utrans,ext_scal_force,normal, &
                            s0_old,s0_new,p0_old,p0_new,dx,dt,the_bc_level,verbose)
 
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
    real(kind=dp_t), intent(inout) :: eta(:,0:,:)
    type(multifab) , intent(inout) :: sedge(:,:)
    type(multifab) , intent(inout) :: utrans(:,:)
    type(multifab) , intent(inout) :: ext_scal_force(:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: s0_old(:,0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    integer        , intent(in   ) :: verbose

    ! local
    type(multifab), allocatable :: scal_force(:)
    type(multifab), allocatable :: s0_old_cart(:)
    type(multifab), allocatable :: s0_new_cart(:)
    
    real(kind=dp_t), allocatable :: s0_edge_old(:,:)
    real(kind=dp_t), allocatable :: s0_edge_new(:,:)

    real(kind=dp_t), pointer :: s0op(:,:,:,:)
    real(kind=dp_t), pointer :: s0np(:,:,:,:)

    real(dp_t) :: mult
    real(dp_t) :: smin,smax

    type(box) :: domain

    integer :: lo(uold(1)%dim),hi(uold(1)%dim)
    integer :: domlo(uold(1)%dim),domhi(uold(1)%dim)
    integer :: velpred,i,comp,n,bc_comp,dm,ng_cell,nr

    logical :: is_vel
    
    allocate(scal_force(nlevs))
    allocate(s0_old_cart(nlevs))
    allocate(s0_new_cart(nlevs))

    nr = size(s0_old,dim=2)
    allocate(s0_edge_old(0:nr,nscal))
    allocate(s0_edge_new(0:nr,nscal))
    
    velpred = 0    
    is_vel = .false.
    ng_cell = sold(1)%ng
    dm = sold(1)%dim

    do n = 1, nlevs
       call cell_to_edge_n(s0_old(n,:,:),s0_edge_old)
       call cell_to_edge_n(s0_new(n,:,:),s0_edge_new)
    enddo

    do n = 1, nlevs
    
       domain = layout_get_pd(uold(n)%la)
       domlo = lwb(domain)
       domhi = upb(domain)
              
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
                call fill_3d_data(s0op(:,:,:,comp),s0_old(n,:,comp),lo,hi,dx(n,:),1)
                call fill_3d_data(s0np(:,:,:,comp),s0_new(n,:,comp),lo,hi,dx(n,:),1)
             end do
             comp = rhoh_comp
             call fill_3d_data(s0op(:,:,:,comp),s0_old(n,:,comp),lo,hi,dx(n,:),1)
             call fill_3d_data(s0np(:,:,:,comp),s0_new(n,:,comp),lo,hi,dx(n,:),1)
          end do
          call multifab_fill_boundary_c(s0_old_cart(n),spec_comp,nspec)
          call multifab_fill_boundary_c(s0_old_cart(n),rhoh_comp,1)
          call multifab_fill_boundary_c(s0_new_cart(n),spec_comp,nspec)
          call multifab_fill_boundary_c(s0_new_cart(n),rhoh_comp,1)
       end if

       ! This can be uncommented if you wish to compute T
       ! call makeTfromRhoH(sold(n),s0_old(n,:,temp_comp))
       
       call modify_scal_force(scal_force(n),sold(n),umac(n,:),s0_old(n,:,:),s0_edge_old, &
                              w0(n,:),dx(n,:),domlo,domhi,s0_old_cart(n),spec_comp,nspec)
       
       if(use_temp_in_mkflux) then
          call mktempforce(scal_force(n),temp_comp,sold(n),thermal(n),p0_old(n,:),dx(n,:))
       else
          call mkrhohforce(scal_force(n),rhoh_comp,umac(n,:),p0_old(n,:),p0_new(n,:), &
                           normal(n),dx(n,:))
          call modify_scal_force(scal_force(n),sold(n),umac(n,:),s0_old(n,:,:),s0_edge_old, &
                                 w0(n,:),dx(n,:),domlo,domhi,s0_old_cart(n),rhoh_comp,1)
       endif
       
       ! add to the rhoh component of force if NOT use_temp_in_mkflux
       if ( (.not. use_temp_in_mkflux) .and. use_thermal_diffusion) &
            call multifab_plus_plus_c(scal_force(n),rhoh_comp,thermal(n),1,1)
       
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
            call put_in_pert_form(sold(n),s0_old(n,:,:),dx(n,:),rhoh_comp,1,.true.)
       
       call put_in_pert_form(sold(n),s0_old(n,:,:),dx(n,:),spec_comp,nspec,.true.)
       
       if (use_temp_in_mkflux) then
          comp = temp_comp
       else
          comp = rhoh_comp
       end if
       
       ! create temperature or enthalpy edge states
       call mkflux(sold(n),uold(n),sedge(n,:),umac(n,:),utrans(n,:),scal_force(n),w0(n,:), &
                   w0_cart_vec(n),dx(n,:),dt,is_vel,the_bc_level(n),velpred,comp,dm+comp,1)
       
       ! create species edge states
       call mkflux(sold(n),uold(n),sedge(n,:),umac(n,:),utrans(n,:),scal_force(n),w0(n,:), &
                   w0_cart_vec(n),dx(n,:),dt,is_vel,the_bc_level(n),velpred, &
                   spec_comp,dm+spec_comp,nspec)
       
       if(use_temp_in_mkflux) then
          call makeRhoHfromT(uold(n),sedge(n,:),s0_old(n,:,:),s0_edge_old,s0_new(n,:,:), &
                             s0_edge_new)
       endif
       
       if (.not. use_temp_in_mkflux) &
            call put_in_pert_form(sold(n),s0_old(n,:,:),dx(n,:),rhoh_comp,1,.false.)
       
       call put_in_pert_form(sold(n),s0_old(n,:,:),dx(n,:),spec_comp,nspec,.false.)
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the edge states of tracers.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (ntrac .ge. 1) then
          call mkflux(sold(n),uold(n),sedge(n,:),umac(n,:),utrans(n,:),scal_force(n), &
                      w0(n,:),w0_cart_vec(n),dx(n,:),dt,is_vel,the_bc_level(n),velpred,&
                      trac_comp,dm+trac_comp,ntrac)
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
       
       call update_scal(which_step,spec_comp,spec_comp+nspec-1,sold(n),snew(n),umac(n,:), &
                        w0(n,:),w0_cart_vec(n),eta(n,:,:),sedge(n,:),scal_force(n), &
                        s0_old(n,:,:),s0_edge_old,s0_new(n,:,:),s0_edge_new,s0_old_cart(n), &
                        s0_new_cart(n),domlo,domhi,dx(n,:),dt,evolve_base_state)
       
       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary conditions
       call multifab_fill_boundary(snew(n))
       
       ! fill ghost cells for physical boundary conditions at domain boundaries
       call multifab_physbc(snew(n),rho_comp ,dm+rho_comp ,1    ,dx(n,:),the_bc_level(n))
       call multifab_physbc(snew(n),spec_comp,dm+spec_comp,nspec,dx(n,:),the_bc_level(n))
       
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
                call fill_3d_data(s0op(:,:,:,comp),s0_old(n,:,comp),lo,hi,dx(n,:),1)
                call fill_3d_data(s0np(:,:,:,comp),s0_new(n,:,comp),lo,hi,dx(n,:),1)
             end do
          end do
          do comp = trac_comp,trac_comp+ntrac-1
             call multifab_fill_boundary_c(s0_old_cart(n),comp,1)
             call multifab_fill_boundary_c(s0_new_cart(n),comp,1)
          end do
       end if
       
       if (ntrac .ge. 1) then
       
          call update_scal(which_step,trac_comp,trac_comp+ntrac-1,sold(n),snew(n), &
                           umac(n,:),w0(n,:),w0_cart_vec(n),eta(n,:,:),sedge(n,:), &
                           scal_force(n),s0_old(n,:,:),s0_edge_old,s0_new(n,:,:), &
                           s0_edge_new,s0_old_cart(n),s0_new_cart(n),domlo,domhi,dx(n,:),dt, &
                           evolve_base_state)
          
          ! fill ghost cells for two adjacent grids at the same level
          ! this includes periodic domain boundary conditions
          call multifab_fill_boundary(snew(n))
          
          ! fill ghost cells for physical boundary conditions at domain boundaries
          call multifab_physbc(snew(n),trac_comp,dm+trac_comp,ntrac,dx(n,:),the_bc_level(n))
          
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

       ! Define s0_old_cart and s0_new_cart
       if (spherical .eq. 1) then
          do i = 1, sold(n)%nboxes
             if ( multifab_remote(s0_old_cart(n), i) ) cycle
             s0op => dataptr(s0_old_cart(n), i)
             s0np => dataptr(s0_new_cart(n), i)
             lo =  lwb(get_box(s0_old_cart(n), i))
             hi =  upb(get_box(s0_old_cart(n), i))
             call fill_3d_data(s0op(:,:,:,rhoh_comp),s0_old(n,:,rhoh_comp),lo,hi,dx(n,:),1)
             call fill_3d_data(s0np(:,:,:,rhoh_comp),s0_new(n,:,rhoh_comp),lo,hi,dx(n,:),1)
          end do
          call multifab_fill_boundary_c(s0_old_cart(n),rhoh_comp,1)
          call multifab_fill_boundary_c(s0_new_cart(n),rhoh_comp,1)
       end if
       
       call mkrhohforce(scal_force(n),rhoh_comp,umac(n,:),p0_old(n,:),p0_new(n,:), &
                        normal(n),dx(n,:))
       
       call update_scal(which_step,rhoh_comp,rhoh_comp,sold(n),snew(n), &
                        umac(n,:),w0(n,:),w0_cart_vec(n),eta(n,:,:),sedge(n,:), &
                        scal_force(n),s0_old(n,:,:),s0_edge_old,s0_new(n,:,:),s0_edge_new, &
                        s0_old_cart(n),s0_new_cart(n),domlo,domhi,dx(n,:),dt, &
                        evolve_base_state)
       
       if(.not. use_thermal_diffusion) then
          call makeTfromRhoH(snew(n),s0_new(n,:,temp_comp))
       endif
       
       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary conditions
       call multifab_fill_boundary(snew(n))
       
       ! fill ghost cells for physical boundary conditions at domain boundaries
       call multifab_physbc(snew(n),rho_comp,dm+rho_comp,nscal,dx(n,:),the_bc_level(n))
       
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

    do n = nlevs, 2, -1
       call ml_cc_restriction(snew(n-1),snew(n),mla%mba%rr(n-1,:))
       
       domain = layout_get_pd(mla%la(n))
       call multifab_fill_ghost_cells(snew(n),snew(n-1),domain, &
                                      ng_cell,mla%mba%rr(n-1,:), &
                                      the_bc_level(n-1), the_bc_level(n), &
                                      1,dm+rho_comp,nscal)
    end do
    
    do n = 1,nlevs
       call multifab_destroy(scal_force(n))
       call multifab_destroy(s0_old_cart(n))
       call multifab_destroy(s0_new_cart(n))
    enddo
    
    deallocate(s0_edge_old,s0_edge_new)
    
2000 format('... new min/max : density           ',e17.10,2x,e17.10)
2001 format('... new min/max : rho * H           ',e17.10,2x,e17.10)
2002 format('... new min/max : ',a16,2x,e17.10,2x,e17.10)
2003 format('... new min/max : tracer            ',e17.10,2x,e17.10)
2004 format(' ')
    
  end subroutine scalar_advance

  subroutine modify_scal_force(force,s,umac,base,base_edge,w0,dx,domlo,domhi,base_cart, &
                               start_comp,num_comp)

    ! When we write the scalar equation in perturbational and convective
    ! form, the terms other than s'_t + U.grad s' act as source terms.  Add
    ! them to the forces here.
    
    type(multifab) , intent(inout) :: force
    type(multifab) , intent(in   ) :: s
    type(multifab) , intent(in   ) :: umac(:)
    real(kind=dp_t), intent(in   ) :: base(0:,:), base_edge(0:,:), w0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: domlo(:),domhi(:)
    type(multifab) , intent(in   ) :: base_cart
    
    ! local
    integer :: i,ng,dm
    integer :: comp,start_comp,num_comp
    integer :: lo(s%dim),hi(s%dim)
    real(kind=dp_t), pointer :: fp(:,:,:,:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: bcp(:,:,:,:)
    
    dm = s%dim
    ng = s%ng
    
    do i=1,force%nboxes
       if ( multifab_remote(force,i) ) cycle
       fp => dataptr(force,i)
       sp => dataptr(s,i)
       ump => dataptr(umac(1),i)
       vmp => dataptr(umac(2),i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       select case (dm)
       case (2)
          do comp = start_comp, start_comp+num_comp-1
             call modify_scal_force_2d(fp(:,:,1,comp),sp(:,:,1,comp), lo, hi, &
                                       ng,ump(:,:,1,1),vmp(:,:,1,1), &
                                       base(:,comp), base_edge(:,comp), w0, dx)
          end do
       case(3)
          wmp  => dataptr(umac(3), i)
          if (spherical .eq. 1) then
             bcp => dataptr(base_cart, i)
             do comp = start_comp, start_comp+num_comp-1
                call modify_scal_force_3d_sphr(fp(:,:,:,comp),sp(:,:,:,comp), &
                                               lo,hi,domlo,domhi,ng, &
                                               ump(:,:,:,1),vmp(:,:,:,1), &
                                               wmp(:,:,:,1),bcp(:,:,:,comp), &
                                               w0,dx)
             end do
          else
             do comp = start_comp, start_comp+num_comp-1
                call modify_scal_force_3d_cart(fp(:,:,:,comp),sp(:,:,:,comp), &
                                               lo,hi,ng,ump(:,:,:,1), &
                                               vmp(:,:,:,1),wmp(:,:,:,1), &
                                               base(:,comp),base_edge(:,comp), &
                                               w0,dx)
             end do
          end if
       end select
    end do
    
  end subroutine modify_scal_force
  
  subroutine modify_scal_force_2d(force,s,lo,hi,ng,umac,vmac,base,base_edge,w0,dx)

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
       divu(k) = (zl(k+1)**2 * w0(k+1)- zl(k)**2 * w0(k))/(dr(1)*z(k)**2)
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
 
