module scalar_advance_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: scalar_advance

contains

  subroutine scalar_advance(nlevs,mla,which_step,uold,sold,snew,thermal,umac,w0, &
                            w0_cart_vec,eta,utrans,normal, &
                            s0_old,s0_new,p0_old,p0_new,dx,dt,the_bc_level)

    use bl_prof_module
    use bl_constants_module
    use make_edge_state_module
    use mkflux_module
    use mkscalforce_module
    use update_scal_module
    use addw0_module
    use define_bc_module
    use fill_3d_module
    use pert_form_module
    use cell_to_edge_module
    use rhoh_vs_t_module
    use network,       only: nspec, spec_names
    use geometry,      only: spherical, nr
    use variables,     only: nscal, ntrac, spec_comp, trac_comp, temp_comp, rho_comp, &
                             rhoh_comp
    use probin_module, only: predict_temp_at_edges, use_thermal_diffusion, verbose, &
                             evolve_base_state, use_eta
    use modify_scal_force_module
    use make_eta_module
    use add_thermal_to_force_module

    integer        , intent(in   ) :: nlevs
    type(ml_layout), intent(inout) :: mla
    integer        , intent(in   ) :: which_step
    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(in   ) :: thermal(:)
    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0_cart_vec(:)
    real(kind=dp_t), intent(inout) :: eta(:,0:,:)
    type(multifab) , intent(in   ) :: utrans(:,:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: s0_old(:,0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    type(multifab) :: scal_force(nlevs)
    type(multifab) :: s0_old_cart(nlevs)
    type(multifab) :: s0_new_cart(nlevs)
    type(multifab) :: sedge(nlevs,mla%dim)
    type(multifab) :: sflux(nlevs,mla%dim)
    type(multifab) :: etaflux(nlevs)

    integer    :: velpred,comp,pred_comp,n,dm
    logical    :: umac_nodal_flag(sold(1)%dim), is_vel
    real(dp_t) :: smin,smax

    real(kind=dp_t), allocatable :: s0_edge_old(:,:,:), s0_edge_new(:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "scalar_advance")

    allocate(s0_edge_old(nlevs,0:nr(nlevs),nscal))
    allocate(s0_edge_new(nlevs,0:nr(nlevs),nscal))

    dm      = sold(1)%dim
    is_vel  = .false.
    velpred = 0    

    do n = 1, nlevs
       call cell_to_edge_allcomps(n,s0_old(n,:,:),s0_edge_old(n,:,:))
       call cell_to_edge_allcomps(n,s0_new(n,:,:),s0_edge_new(n,:,:))
    end do

    !**************************************************************************
    !     Create scalar source term at time n for (rho X)_i and (rho H).  
    !     The source term for (rho X) is zero.
    !     The source term for (rho h) has only the w dp0/dr term.
    !
    !     The call to modify_scal_force is used to add those advective terms 
    !     that appear as forces when we write it in convective/perturbational form.
    !**************************************************************************

    ! Define s0_old_cart and s0_new_cart
    if (spherical .eq. 1) then
       do n=1,nlevs
          call build(s0_old_cart(n), sold(n)%la, nscal, 1)
          call build(s0_new_cart(n), sold(n)%la, nscal, 1)
       end do
       call fill_3d_data_c(nlevs,dx,the_bc_level,mla, &
                           s0_old_cart,s0_old(:,:,rhoh_comp),rhoh_comp,dm+rhoh_comp)
       call fill_3d_data_c(nlevs,dx,the_bc_level,mla, &
                           s0_new_cart,s0_new(:,:,rhoh_comp),rhoh_comp,dm+rhoh_comp)
       do comp = spec_comp, spec_comp+nspec-1
          call fill_3d_data_c(nlevs,dx,the_bc_level,mla, &
                              s0_old_cart,s0_old(:,:,comp),comp,dm+comp)
          call fill_3d_data_c(nlevs,dx,the_bc_level,mla, &
                              s0_new_cart,s0_new(:,:,comp),comp,dm+comp)
       end do
       do comp = trac_comp, trac_comp+ntrac-1
          call fill_3d_data_c(nlevs,dx,the_bc_level,mla, &
                              s0_old_cart,s0_old(:,:,comp),comp,dm+comp)
          call fill_3d_data_c(nlevs,dx,the_bc_level,mla, &
                              s0_new_cart,s0_new(:,:,comp),comp,dm+comp)
       end do
    end if

    ! This can be uncommented if you wish to compute T
    ! call makeTfromRhoH(nlevs,sold,s0_old(:,:,temp_comp),mla,the_bc_level,dx)

    do n = 1, nlevs
       call build(scal_force(n), sold(n)%la, nscal, 1)       
       call setval(scal_force(n),ZERO,all=.true.)
    end do

    ! make force for species
    call modify_scal_force(which_step,nlevs,scal_force,sold,umac,s0_old,s0_edge_old,w0,eta,dx, &
                           s0_old_cart,spec_comp,nspec,mla,the_bc_level)

    if(predict_temp_at_edges) then

       ! make force for temperature
       call mktempforce(nlevs,scal_force,temp_comp,umac,sold,thermal,p0_old,p0_old,normal, &
                        dx,mla,the_bc_level)

    else

       ! Make force for rhoh -- just use p0_old
       call mkrhohforce(nlevs,scal_force,rhoh_comp,umac,p0_old,p0_old,normal,dx, &
                        mla,the_bc_level)
       
       call modify_scal_force(which_step,nlevs,scal_force,sold,umac,s0_old,s0_edge_old,w0,eta,dx, &
                              s0_old_cart,rhoh_comp,1,mla,the_bc_level)
        
       if(use_thermal_diffusion) then
          call add_thermal_to_force(nlevs,scal_force,thermal,the_bc_level,mla)
       end if

    end if
      
    !**************************************************************************
    !     Add w0 to MAC velocities (trans velocities already have w0).
    !**************************************************************************

    call addw0(nlevs,umac,w0,w0_cart_vec,mult=ONE)

    !**************************************************************************
    !     Create the edge states of (rho h)' (or T) and (rho X)_i.
    !**************************************************************************

    if (.not. predict_temp_at_edges) then
       call put_in_pert_form(nlevs,sold,s0_old,dx,rhoh_comp,1,.true.,mla,the_bc_level)
    end if

    call put_in_pert_form(nlevs,sold,s0_old,dx,spec_comp,nspec,.true.,mla,the_bc_level)

    do n=1,nlevs
       do comp = 1,dm
          umac_nodal_flag = .false.
          umac_nodal_flag(comp) = .true.
          call multifab_build( sedge(n,comp), mla%la(n), nscal, 0, nodal = umac_nodal_flag)
       end do
    end do

    if (predict_temp_at_edges) then
       pred_comp = temp_comp
    else
       pred_comp = rhoh_comp
    end if
    
    ! create temperature or enthalpy edge states
    call make_edge_state(nlevs,sold,uold,sedge,umac,utrans,scal_force,w0, &
                         w0_cart_vec,dx,dt,is_vel,the_bc_level,velpred,pred_comp, &
                         dm+pred_comp,1,mla)

    ! create species edge states
    call make_edge_state(nlevs,sold,uold,sedge,umac,utrans,scal_force,w0, &
                         w0_cart_vec,dx,dt,is_vel,the_bc_level,velpred,spec_comp, &
                         dm+spec_comp,nspec,mla)

    ! compute enthalpy edge states
    if(predict_temp_at_edges) then
       call makeRhoHfromT(nlevs,uold,sedge,s0_old,s0_edge_old,s0_new,s0_edge_new,the_bc_level,dx)
    end if

    if (.not. predict_temp_at_edges) then
       call put_in_pert_form(nlevs,sold,s0_old,dx,rhoh_comp,1,.false.,mla,the_bc_level)
    end if
    
    call put_in_pert_form(nlevs,sold,s0_old,dx,spec_comp,nspec,.false.,mla,the_bc_level)

    !**************************************************************************
    !     Create the edge states of tracers.
    !**************************************************************************

    if (ntrac .ge. 1) then
       call make_edge_state(nlevs,sold,uold,sedge,umac,utrans,scal_force,w0, &
                            w0_cart_vec,dx,dt,is_vel,the_bc_level,velpred,trac_comp, &
                            dm+trac_comp,ntrac,mla)
    end if

    !**************************************************************************
    !     Subtract w0 from MAC velocities.
    !**************************************************************************

    call addw0(nlevs,umac,w0,w0_cart_vec,mult=-ONE)

    !**************************************************************************
    !     Compute fluxes
    !**************************************************************************

    do n=1,nlevs
       do comp = 1,dm
          umac_nodal_flag = .false.
          umac_nodal_flag(comp) = .true.
          call multifab_build(sflux(n,comp), mla%la(n), nscal, 0, nodal = umac_nodal_flag)
       end do

       umac_nodal_flag = .false.
       umac_nodal_flag(dm) = .true.
       call multifab_build(etaflux(n), mla%la(n), nscal, 0, nodal = umac_nodal_flag)
    end do

    ! compute enthalpy fluxes
    call mkflux(nlevs,sflux,etaflux,sold,sedge,umac,w0,w0_cart_vec,s0_old,s0_edge_old, &
                s0_old_cart,s0_new,s0_edge_new,s0_new_cart,rhoh_comp,rhoh_comp, &
                which_step,mla)

    ! compute species fluxes
    call mkflux(nlevs,sflux,etaflux,sold,sedge,umac,w0,w0_cart_vec,s0_old,s0_edge_old, &
                s0_old_cart,s0_new,s0_edge_new,s0_new_cart,spec_comp, &
                spec_comp+nspec-1,which_step,mla)

    if (ntrac .ge. 1) then
       ! compute tracer fluxes
       call mkflux(nlevs,sflux,etaflux,sold,sedge,umac,w0,w0_cart_vec,s0_old,s0_edge_old, &
                   s0_old_cart,s0_new,s0_edge_new,s0_new_cart,trac_comp, &
                   trac_comp+ntrac-1,which_step,mla)
    end if

    !**************************************************************************
    !     1) Set force for (rho X)_i at time n+1/2 = 0.
    !     2) Update (rho X)'_i with conservative differencing.
    !     3) Define density as the sum of the (rho X)_i
    !**************************************************************************
    
    do n=1,nlevs
       call setval(scal_force(n),ZERO,all=.true.)
    end do

    call update_scal(nlevs,spec_comp,spec_comp+nspec-1,sold,snew,umac,w0, &
                     w0_cart_vec,eta,sedge,sflux,scal_force,s0_old,s0_edge_old,s0_new, &
                     s0_edge_new,s0_old_cart,s0_new_cart,dx,dt,the_bc_level,mla)
    
    if ( verbose .ge. 1 ) then
       do n=1, nlevs
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
       end do
    end if

    !**************************************************************************
    !     2) Update tracers with convective differencing.
    !**************************************************************************
    
    if ( ntrac .ge. 1 ) then
       call update_scal(nlevs,trac_comp,trac_comp+ntrac-1,sold,snew,umac,w0, &
                        w0_cart_vec,eta,sedge,sflux,scal_force,s0_old,s0_edge_old,s0_new, &
                        s0_edge_new,s0_old_cart,s0_new_cart,dx,dt,the_bc_level,mla)

       if ( verbose .ge. 1 ) then
          do n=1,nlevs
             smin = multifab_min_c(snew(n),trac_comp) 
             smax = multifab_max_c(snew(n),trac_comp)
             if (parallel_IOProcessor()) &
                  write(6,2003) smin,smax
          end do
       end if
    end if

    !**************************************************************************
    !     1) Create (rhoh)' force at time n+1/2.
    !          (NOTE: we don't worry about filling ghost cells of the scal_force
    !                 because we only need them in valid regions...)     
    !     2) Update (rho h)' with conservative differencing.
    !**************************************************************************
       
    if (which_step .eq. 1) then
      ! Here just send p0_old and p0_old
      call mkrhohforce(nlevs,scal_force,rhoh_comp,umac,p0_old,p0_old,normal,dx, &
                       mla,the_bc_level)
    else
      ! Here send p0_old and p0_new
      call mkrhohforce(nlevs,scal_force,rhoh_comp,umac,p0_old,p0_new,normal,dx, &
                       mla,the_bc_level)
    end if

    call update_scal(nlevs,rhoh_comp,rhoh_comp,sold,snew,umac,w0,w0_cart_vec, &
                     eta,sedge,sflux,scal_force,s0_old,s0_edge_old,s0_new,s0_edge_new, &
                     s0_old_cart,s0_new_cart,dx,dt,the_bc_level,mla)

    if (use_eta .and. evolve_base_state) then
       call make_eta(nlevs,eta,sold,etaflux,mla)
    end if

    if (spherical .eq. 1) then
       do n=1,nlevs
          call destroy(s0_old_cart(n))
          call destroy(s0_new_cart(n))
       end do
    end if

    do n = 1, nlevs
       call destroy(scal_force(n))
       call destroy(etaflux(n))
       do comp = 1,dm
          call destroy(sedge(n,comp))
          call destroy(sflux(n,comp))
       end do
    end do

    if (.not. use_thermal_diffusion) then
       call makeTfromRhoH(nlevs,snew,s0_new(:,:,temp_comp),mla,the_bc_level,dx)
    end if

    if ( verbose .ge. 1 ) then
       do n=1,nlevs
          smin = multifab_min_c(snew(n),rhoh_comp) 
          smax = multifab_max_c(snew(n),rhoh_comp)
          if (parallel_IOProcessor()) then
             write(6,2001) smin,smax
             write(6,2004) 
          end if
       end do
    end if

    call destroy(bpt)

2000 format('... new min/max : density           ',e17.10,2x,e17.10)
2001 format('... new min/max : rho * H           ',e17.10,2x,e17.10)
2002 format('... new min/max : ',a16,2x,e17.10,2x,e17.10)
2003 format('... new min/max : tracer            ',e17.10,2x,e17.10)
2004 format(' ')

  end subroutine scalar_advance

end module scalar_advance_module
