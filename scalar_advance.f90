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
  use modify_scal_force_module

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
    integer :: velpred,i,comp,n,bc_comp,dm,ng_cell

    logical :: is_vel
    
    allocate(scal_force(nlevs))
    allocate(s0_old_cart(nlevs))
    allocate(s0_new_cart(nlevs))

    allocate(s0_edge_old(0:nr(nlevs),nscal))
    allocate(s0_edge_new(0:nr(nlevs),nscal))
    
    velpred = 0    
    is_vel = .false.
    ng_cell = sold(1)%ng
    dm = sold(1)%dim

    do n = 1, nlevs
       call cell_to_edge_allcomps(n,s0_old(n,:,:),s0_edge_old)
       call cell_to_edge_allcomps(n,s0_new(n,:,:),s0_edge_new)
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
          do comp = spec_comp, spec_comp+nspec-1
             call fill_3d_data_wrapper(n,s0_old_cart(n),s0_old(n,:,comp),dx(n,:),comp)
             call fill_3d_data_wrapper(n,s0_new_cart(n),s0_new(n,:,comp),dx(n,:),comp)
          enddo
          call fill_3d_data_wrapper(n,s0_old_cart(n),s0_old(n,:,rhoh_comp),dx(n,:),rhoh_comp)
          call fill_3d_data_wrapper(n,s0_new_cart(n),s0_new(n,:,rhoh_comp),dx(n,:),rhoh_comp)
       end if

       ! This can be uncommented if you wish to compute T
       ! call makeTfromRhoH(sold(n),s0_old(n,:,temp_comp))
       ! call multifab_fill_boundary_c(sold(n),temp_comp,1)
       ! call multifab_physbc(sold(n),temp_comp,dm+temp_comp,1,dx(n,:),the_bc_level(n))
       
       ! make force for species
       call modify_scal_force(n,scal_force(n),sold(n),umac(n,:),s0_old(n,:,:),s0_edge_old, &
                              w0(n,:),dx(n,:),domlo,domhi,s0_old_cart(n),spec_comp,nspec)

       call multifab_fill_boundary_c(scal_force(n),spec_comp,nspec)

       do comp = spec_comp, spec_comp+nspec-1
          call multifab_physbc(scal_force(n),comp,neumann_comp,1,dx(n,:),the_bc_level(n))
       enddo
       
       if(use_temp_in_mkflux) then

          ! make force for temperature
          call mktempforce(n,scal_force(n),temp_comp,sold(n),thermal(n),p0_old(n,:),dx(n,:))

          call multifab_fill_boundary_c(scal_force(n),temp_comp,1)

          call multifab_physbc(scal_force(n),temp_comp,neumann_comp,1,dx(n,:), &
                               the_bc_level(n))
       else

          ! make force for rhoh
          call mkrhohforce(n,scal_force(n),rhoh_comp,umac(n,:),p0_old(n,:),p0_old(n,:), &
                           normal(n),dx(n,:))

          call modify_scal_force(n,scal_force(n),sold(n),umac(n,:),s0_old(n,:,:), &
                                 s0_edge_old,w0(n,:),dx(n,:),domlo,domhi,s0_old_cart(n), &
                                 rhoh_comp,1)

          if(use_thermal_diffusion) then
            call multifab_plus_plus_c(scal_force(n),rhoh_comp,thermal(n),1,1)
          endif

          call multifab_fill_boundary_c(scal_force(n),rhoh_comp,1)

          call multifab_physbc(scal_force(n),rhoh_comp,neumann_comp,1,dx(n,:), &
                               the_bc_level(n))
       endif
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Add w0 to MAC velocities (trans velocities already have w0).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       mult = ONE
       call addw0(umac(n,:),w0(n,:),w0_cart_vec(n),dx(n,:),mult)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the edge states of (rho h)' and (rho X)_i.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (.not. use_temp_in_mkflux) &
            call put_in_pert_form(n,sold(n),s0_old(n,:,:),dx(n,:),rhoh_comp,1,.true.)
       
       call put_in_pert_form(n,sold(n),s0_old(n,:,:),dx(n,:),spec_comp,nspec,.true.)
       
       if (use_temp_in_mkflux) then
          comp = temp_comp
       else
          comp = rhoh_comp
       end if
       
       ! create temperature or enthalpy edge states
       call mkflux(n,sold(n),uold(n),sedge(n,:),umac(n,:),utrans(n,:),scal_force(n), &
                   w0(n,:),w0_cart_vec(n),dx(n,:),dt,is_vel,the_bc_level(n),velpred,comp, &
                   dm+comp,1)
       
       ! create species edge states
       call mkflux(n,sold(n),uold(n),sedge(n,:),umac(n,:),utrans(n,:),scal_force(n), &
                   w0(n,:),w0_cart_vec(n),dx(n,:),dt,is_vel,the_bc_level(n),velpred, &
                   spec_comp,dm+spec_comp,nspec)
       
       if(use_temp_in_mkflux) then
          call makeRhoHfromT(uold(n),sedge(n,:),s0_old(n,:,:),s0_edge_old,s0_new(n,:,:), &
                             s0_edge_new)
       endif
       
       if (.not. use_temp_in_mkflux) &
            call put_in_pert_form(n,sold(n),s0_old(n,:,:),dx(n,:),rhoh_comp,1,.false.)
       
       call put_in_pert_form(n,sold(n),s0_old(n,:,:),dx(n,:),spec_comp,nspec,.false.)
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Create the edge states of tracers.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (ntrac .ge. 1) then
          call mkflux(n,sold(n),uold(n),sedge(n,:),umac(n,:),utrans(n,:),scal_force(n), &
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
          do comp = trac_comp, trac_comp+ntrac-1
             call fill_3d_data_wrapper(n,s0_old_cart(n),s0_old(n,:,comp),dx(n,:),comp)
             call fill_3d_data_wrapper(n,s0_new_cart(n),s0_new(n,:,comp),dx(n,:),comp)
          enddo
       end if
       
       if (ntrac .ge. 1) then
          call update_scal(which_step,trac_comp,trac_comp+ntrac-1,sold(n),snew(n), &
                           umac(n,:),w0(n,:),w0_cart_vec(n),eta(n,:,:),sedge(n,:), &
                           scal_force(n),s0_old(n,:,:),s0_edge_old,s0_new(n,:,:), &
                           s0_edge_new,s0_old_cart(n),s0_new_cart(n),domlo,domhi,dx(n,:), &
                           dt,evolve_base_state)
          
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
          call fill_3d_data_wrapper(n,s0_old_cart(n),s0_old(n,:,rhoh_comp),dx(n,:),rhoh_comp)
          call fill_3d_data_wrapper(n,s0_new_cart(n),s0_new(n,:,rhoh_comp),dx(n,:),rhoh_comp)
       end if
       
       call mkrhohforce(n,scal_force(n),rhoh_comp,umac(n,:),p0_old(n,:),p0_new(n,:), &
                        normal(n),dx(n,:))
       
       call update_scal(which_step,rhoh_comp,rhoh_comp,sold(n),snew(n), &
                        umac(n,:),w0(n,:),w0_cart_vec(n),eta(n,:,:),sedge(n,:), &
                        scal_force(n),s0_old(n,:,:),s0_edge_old,s0_new(n,:,:),s0_edge_new, &
                        s0_old_cart(n),s0_new_cart(n),domlo,domhi,dx(n,:),dt, &
                        evolve_base_state)
       
       if(.not. use_thermal_diffusion) then
          call makeTfromRhoH(snew(n),s0_new(n,:,temp_comp))
       endif
       
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
       
       call multifab_fill_boundary(snew(n))
       call multifab_physbc(snew(n),rho_comp,dm+rho_comp,nscal,dx(n,:),the_bc_level(n))

    enddo ! do n = 1, nlevs

    do n = nlevs, 2, -1
       call ml_cc_restriction(snew(n-1),snew(n),mla%mba%rr(n-1,:))
       
       call multifab_fill_ghost_cells(snew(n),snew(n-1), &
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

end module scalar_advance_module
