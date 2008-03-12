module scalar_advance_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: scalar_advance

contains

  subroutine scalar_advance(nlevs,mla,which_step,uold,sold,snew,thermal, &
                            umac,w0,w0_cart_vec,eta,utrans,normal, &
                            s0_old,s0_new,p0_old,p0_new, &
                            s0_predicted_edge, &
                            dx,dt,the_bc_level)

    use bl_prof_module
    use bl_constants_module
    use make_edge_state_module
    use make_edge_scal_module
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
    use variables,     only: nscal, ntrac, spec_comp, trac_comp, temp_comp, &
                             rho_comp, rhoh_comp
    use probin_module, only: predict_temp_at_edges, predict_X_at_edges, predict_h_at_edges, &
                             use_thermal_diffusion, verbose, evolve_base_state, use_eta
    use modify_scal_force_module
    use make_eta_module
    use convert_rhoX_to_X_module

    integer        , intent(in   ) :: nlevs
    type(ml_layout), intent(inout) :: mla
    integer        , intent(in   ) :: which_step
    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: thermal(:)
    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0_cart_vec(:)
    real(kind=dp_t), intent(inout) :: eta(:,0:,:)
    type(multifab) , intent(in   ) :: utrans(:,:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(inout) :: s0_old(:,0:,:)
    real(kind=dp_t), intent(in   ) :: s0_new(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: s0_predicted_edge(:,0:,:)
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

    real(kind=dp_t), allocatable :: s0_edge_old(:,:,:)
    real(kind=dp_t), allocatable :: s0_edge_new(:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "scalar_advance")

    dm      = sold(1)%dim
    is_vel  = .false.
    velpred = 0    

    ! create edge-centered base state quantities.  Note: s0_edge_{old,new} 
    ! contain edge-centered quantities created via spatial interpolation.
    ! This is to be contrasted to s0_predicted_edge which is the half-time
    ! edge state created in advect_base.
    allocate(s0_edge_old(nlevs,0:nr(nlevs),nscal))
    allocate(s0_edge_new(nlevs,0:nr(nlevs),nscal))

    do n = 1, nlevs
       call cell_to_edge_allcomps(n,s0_old(n,:,:),s0_edge_old(n,:,:))
       call cell_to_edge_allcomps(n,s0_new(n,:,:),s0_edge_new(n,:,:))
    end do

    ! Define s0_old_cart and s0_new_cart
    if (spherical .eq. 1) then
       do n=1,nlevs
          call build(s0_old_cart(n), sold(n)%la, nscal, 1)
          call build(s0_new_cart(n), sold(n)%la, nscal, 1)
       end do

       call fill_3d_data_c(nlevs,dx,the_bc_level,mla,s0_old_cart,s0_old(:,:,rhoh_comp), &
                           rhoh_comp,dm+rhoh_comp)
       call fill_3d_data_c(nlevs,dx,the_bc_level,mla,s0_new_cart,s0_new(:,:,rhoh_comp), &
                           rhoh_comp,dm+rhoh_comp)

       do comp = spec_comp, spec_comp+nspec-1
          call fill_3d_data_c(nlevs,dx,the_bc_level,mla,s0_old_cart,s0_old(:,:,comp), &
                              comp,dm+comp)
          call fill_3d_data_c(nlevs,dx,the_bc_level,mla,s0_new_cart,s0_new(:,:,comp), &
                              comp,dm+comp)
       end do

       do comp = trac_comp, trac_comp+ntrac-1
          call fill_3d_data_c(nlevs,dx,the_bc_level,mla,s0_old_cart,s0_old(:,:,comp), &
                              comp,dm+comp)
          call fill_3d_data_c(nlevs,dx,the_bc_level,mla,s0_new_cart,s0_new(:,:,comp), &
                              comp,dm+comp)
       end do
    end if

    ! This can be uncommented if you wish to compute T
    ! call makeTfromRhoH(nlevs,sold,s0_old(:,:,temp_comp),mla,the_bc_level,dx)

    ! if we are predicting X on the edges, then convert the state arrays
    ! (and base state) from (rho X) to X.  Note, only the time-level n
    ! stuff need be converted, since that's all the prediction uses
    if (predict_X_at_edges) then

       call convert_rhoX_to_X(nlevs,sold,dx,.true.,mla,the_bc_level)

       if (spherical .eq. 1) then
          call convert_rhoX_to_X(nlevs,s0_old_cart,dx,.true.,mla,the_bc_level)
       end if

       do comp = spec_comp, spec_comp + nspec - 1
          s0_old(:,:,comp) = s0_old(:,:,comp)/s0_old(:,:,rho_comp)
          s0_edge_old(:,:,comp) = s0_edge_old(:,:,comp)/s0_edge_old(:,:,rho_comp)
       enddo

    endif

    ! if we are predicting h on the edges, then convert the state arrays
    ! (and base state) from (rho h) to h.  Note, only the time-level n
    ! stuff need be converted, since that's all the prediction uses
    if (predict_h_at_edges) then

       call convert_rhoh_to_h(nlevs,sold,dx,.true.,mla,the_bc_level)

       if (spherical .eq. 1) then
          call convert_rhoh_to_h(nlevs,s0_old_cart,dx,.true.,mla,the_bc_level)
       end if

       s0_old(:,:,rhoh_comp) = s0_old(:,:,rhoh_comp)/s0_old(:,:,rho_comp)
       s0_edge_old(:,:,rhoh_comp) = s0_edge_old(:,:,rhoh_comp)/s0_edge_old(:,:,rho_comp)

    end if

    !**************************************************************************
    ! Create scalar source term at time n
    ! if predict_X_at_edges is true, we compute a source term for X
    ! if predict_X_at_edges is false, we compute a source term for (rho X)'
    ! if predict_h_at_edges is true, we compute a source term for h
    ! if predict_temp_at_edges is true, we compute a source term for temperature
    ! if predict_h_at_edges and predict_temp_at_edges are both false, we compute a
    !    source term for (rho h)'
    ! if either predict_X_at_edges or predic_h_at_edges is true, we compute a source
    !    term for rho'
    ! 
    ! The call to modify_scal_force is used to add those advective terms 
    ! that appear as forces when we write it in convective/perturbational form.
    !**************************************************************************

    do n = 1, nlevs
       call build(scal_force(n), sold(n)%la, nscal, 1)       
       call setval(scal_force(n),ZERO,all=.true.)
    end do

    ! make force for rho', X, and/or (rho X)'
    if (predict_X_at_edges) then

       ! X force is zero - do nothing

       ! make force for rho'
       if (predict_X_at_edges) then
          call modify_scal_force(which_step,nlevs,scal_force,sold,umac,s0_old,s0_edge_old, &
                                 w0,dx,s0_old_cart,rho_comp,1,mla,the_bc_level)
       end if

    else

       ! make force for (rho X)'
       call modify_scal_force(which_step,nlevs,scal_force,sold,umac,s0_old,s0_edge_old,w0,&
                              dx,s0_old_cart,spec_comp,nspec,mla,the_bc_level)

    endif

    ! make force for either h, T, or (rho h)'
    if (predict_h_at_edges) then

       ! make force for h by calling mkrhohforce then dividing by rho
       call mkrhohforce(nlevs,scal_force,thermal,umac,p0_old,p0_old,normal,dx,.true., &
                        mla,the_bc_level)
       do n=1,nlevs
          call multifab_div_div_c(scal_force(n),rhoh_comp,sold(n),rho_comp,1,1)
       end do

    else if (predict_temp_at_edges) then

       ! make force for temperature
       call mktempforce(nlevs,scal_force,umac,sold,thermal,p0_old,p0_old,normal,dx,mla, &
                        the_bc_level)

    else

       ! make force for (rho h)'
       call mkrhohforce(nlevs,scal_force,thermal,umac,p0_old,p0_old,normal,dx,.true., &
                        mla,the_bc_level)

       call modify_scal_force(which_step,nlevs,scal_force,sold,umac,s0_old,s0_edge_old,w0,&
                              dx,s0_old_cart,rhoh_comp,1,mla,the_bc_level)
        
    end if
      
    !**************************************************************************
    !     Add w0 to MAC velocities (trans velocities already have w0).
    !**************************************************************************

    call addw0(nlevs,umac,w0,w0_cart_vec,mult=ONE)

    !**************************************************************************
    !     Create the edge states of (rho h)' or h or T and (rho X)' or X and rho'
    !**************************************************************************

    if (.not. predict_temp_at_edges .and. .not. predict_h_at_edges) then
       ! convert (rho h) -> (rho h)'
       call put_in_pert_form(nlevs,sold,s0_old,dx,rhoh_comp,1,.true.,mla,the_bc_level)
    end if

    if (predict_X_at_edges) then
       ! convert rho -> rho'
       call put_in_pert_form(nlevs,sold,s0_old,dx,rho_comp,1,.true.,mla,the_bc_level)
    else
       ! convert (rho X) -> (rho X)'
       call put_in_pert_form(nlevs,sold,s0_old,dx,spec_comp,nspec,.true.,mla,the_bc_level)
    end if

    do n=1,nlevs
       do comp = 1,dm
          umac_nodal_flag = .false.
          umac_nodal_flag(comp) = .true.
          call multifab_build(sedge(n,comp), mla%la(n), nscal, 0, nodal = umac_nodal_flag)
       end do
    end do

    ! predict either T, h, or (rho h)' at the edges
    if (predict_temp_at_edges) then
       pred_comp = temp_comp
    else
       pred_comp = rhoh_comp
    end if
!   call make_edge_state(nlevs,sold,uold,sedge,umac,utrans,scal_force,w0,w0_cart_vec,dx,dt, &
!                        is_vel,the_bc_level,velpred,pred_comp,dm+pred_comp,1,mla)
    call make_edge_scal(nlevs,sold,sedge,umac,scal_force,w0,w0_cart_vec,dx,dt,is_vel, &
                        the_bc_level,pred_comp,dm+pred_comp,1,mla)

    ! predict either X or (rho X)' at the edges
!      call make_edge_state(nlevs,sold,uold,sedge,umac,utrans,scal_force,w0,w0_cart_vec,dx, &
!                           dt,is_vel,the_bc_level,velpred,spec_comp,dm+spec_comp,nspec,mla)
       call make_edge_scal(nlevs,sold,sedge,umac,scal_force,w0,w0_cart_vec,dx,dt,is_vel, &
                           the_bc_level,spec_comp,dm+spec_comp,nspec,mla)

    if (predict_X_at_edges) then
       ! predict rho' at the edges
!      call make_edge_state(nlevs,sold,uold,sedge,umac,utrans,scal_force,w0,w0_cart_vec,dx, &
!                           dt,is_vel,the_bc_level,velpred,rho_comp,dm+rho_comp,1,mla)
       call make_edge_scal(nlevs,sold,sedge,umac,scal_force,w0,w0_cart_vec,dx,dt,is_vel, &
                           the_bc_level,rho_comp,dm+rho_comp,1,mla)
    end if

    if (.not. predict_temp_at_edges .and. .not. predict_h_at_edges) then
       ! convert (rho h)' -> (rho h)
       call put_in_pert_form(nlevs,sold,s0_old,dx,rhoh_comp,1,.false.,mla,the_bc_level)
    end if

    if (predict_X_at_edges) then
       ! convert rho' -> rho
       call put_in_pert_form(nlevs,sold,s0_old,dx,rho_comp,1,.false.,mla,the_bc_level)
    else
       ! convert (rho X)' -> (rho X)
       call put_in_pert_form(nlevs,sold,s0_old,dx,spec_comp,nspec,.false.,mla,the_bc_level)
    end if

    ! if we were predicting X at the edges, then restore the state arrays 
    ! (and base state) from X to (rho X)
    if (predict_X_at_edges) then

       call convert_rhoX_to_X(nlevs,sold,dx,.false.,mla,the_bc_level)

       if (spherical .eq. 1) then
          call convert_rhoX_to_X(nlevs,s0_old_cart,dx,.false.,mla,the_bc_level)
       end if

       do comp = spec_comp, spec_comp + nspec - 1
          s0_old(:,:,comp) = s0_old(:,:,rho_comp)*s0_old(:,:,comp)
          s0_edge_old(:,:,comp) = s0_edge_old(:,:,rho_comp)*s0_edge_old(:,:,comp)
       enddo

    endif

    ! if we are predicing h at the edges, then restore the state arrays
    ! (and base state) from h to (rho h)
    if (predict_h_at_edges) then

       call convert_rhoh_to_h(nlevs,sold,dx,.false.,mla,the_bc_level)

       if (spherical .eq. 1) then
          call convert_rhoh_to_h(nlevs,s0_old_cart,dx,.false.,mla,the_bc_level)
       end if

       s0_old(:,:,rhoh_comp) = s0_old(:,:,rho_comp)*s0_old(:,:,rhoh_comp)
       s0_edge_old(:,:,rhoh_comp) = s0_edge_old(:,:,rho_comp)*s0_edge_old(:,:,rhoh_comp)

    end if

    ! Compute enthalpy edge states if we were predicting temperature.  This
    ! needs to be done after the state was returned to the full state, and 
    ! the species are back to (rho X) instead of X.
    if (predict_temp_at_edges) then
       call makeRhoHfromT(nlevs,uold,sedge,s0_old,s0_edge_old,s0_new,s0_edge_new, &
                          the_bc_level,dx)
    end if

    !**************************************************************************
    !     Create the edge states of tracers.
    !**************************************************************************

    if (ntrac .ge. 1) then
       call make_edge_state(nlevs,sold,uold,sedge,umac,utrans,scal_force,w0, &
                            w0_cart_vec,dx,dt,is_vel,the_bc_level,velpred, &
                            trac_comp,dm+trac_comp,ntrac,mla)
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

    ! for which_step .eq. 1, we pass in only the old base state quantities
    ! (s0_old, s0_edge_old, s0_old_cart)
    ! for which_step .eq. 2, we pass in the old and new for averaging within mkflux
    if (which_step .eq. 1) then

    ! compute enthalpy fluxes
       call mkflux(nlevs,sflux,etaflux,sold,sedge,umac,w0,w0_cart_vec, &
                   s0_old,s0_edge_old,s0_old_cart,s0_old,s0_edge_old,s0_old_cart, &
                   s0_predicted_edge,rhoh_comp,rhoh_comp,mla,dx,dt)

       ! compute species fluxes
       call mkflux(nlevs,sflux,etaflux,sold,sedge,umac,w0,w0_cart_vec, &
                   s0_old,s0_edge_old,s0_old_cart,s0_old,s0_edge_old,s0_old_cart, &
                   s0_predicted_edge,spec_comp,spec_comp+nspec-1,mla,dx,dt)

       if (ntrac .ge. 1) then
          ! compute tracer fluxes
          call mkflux(nlevs,sflux,etaflux,sold,sedge,umac,w0,w0_cart_vec, &
                      s0_old,s0_edge_old,s0_old_cart,s0_old,s0_edge_old,s0_old_cart, &
                      s0_predicted_edge,trac_comp,trac_comp+ntrac-1,mla,dx,dt)
       end if

    else if (which_step .eq. 2) then

       ! compute enthalpy fluxes
       call mkflux(nlevs,sflux,etaflux,sold,sedge,umac,w0,w0_cart_vec, &
                   s0_old,s0_edge_old,s0_old_cart,s0_new,s0_edge_new,s0_new_cart, &
                   s0_predicted_edge,rhoh_comp,rhoh_comp,mla,dx,dt)

       ! compute species fluxes
       call mkflux(nlevs,sflux,etaflux,sold,sedge,umac,w0,w0_cart_vec, &
                   s0_old,s0_edge_old,s0_old_cart,s0_new,s0_edge_new,s0_new_cart, &
                   s0_predicted_edge,spec_comp,spec_comp+nspec-1,mla,dx,dt)

       if (ntrac .ge. 1) then
          ! compute tracer fluxes
          call mkflux(nlevs,sflux,etaflux,sold,sedge,umac,w0,w0_cart_vec, &
                      s0_old,s0_edge_old,s0_old_cart,s0_new,s0_edge_new,s0_new_cart, &
                      s0_predicted_edge,trac_comp,trac_comp+ntrac-1,mla,dx,dt)
       end if

    end if


    !**************************************************************************
    !     1) Set force for (rho X)'_i at time n+1/2 = 0.
    !     2) Update (rho X)_i with conservative differencing.
    !     3) Define density as the sum of the (rho X)_i
    !**************************************************************************
    
    do n=1,nlevs
       call setval(scal_force(n),ZERO,all=.true.)
    end do

    call update_scal(nlevs,spec_comp,spec_comp+nspec-1,sold,snew,umac,w0, &
                     w0_cart_vec,sedge,sflux,scal_force,s0_old,s0_edge_old,s0_new, &
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
    !     Update tracers with convective differencing.
    !**************************************************************************
    
    if ( ntrac .ge. 1 ) then
       call update_scal(nlevs,trac_comp,trac_comp+ntrac-1,sold,snew,umac,w0, &
                        w0_cart_vec,sedge,sflux,scal_force,s0_old,s0_edge_old,s0_new, &
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
    !     1) Create (rho h)' force at time n+1/2.
    !          (NOTE: we don't worry about filling ghost cells of the scal_force
    !                 because we only need them in valid regions...)     
    !     2) Update (rho h) with conservative differencing.
    !**************************************************************************
       
    if (which_step .eq. 1) then
      ! Here just send p0_old and p0_old
      call mkrhohforce(nlevs,scal_force,thermal,umac,p0_old,p0_old,normal,dx,.false., &
                       mla,the_bc_level)
    else
      ! Here send p0_old and p0_new
      call mkrhohforce(nlevs,scal_force,thermal,umac,p0_old,p0_new,normal,dx,.false., &
                       mla,the_bc_level)
    end if

    call update_scal(nlevs,rhoh_comp,rhoh_comp,sold,snew,umac,w0,w0_cart_vec, &
                     sedge,sflux,scal_force,s0_old,s0_edge_old,s0_new,s0_edge_new, &
                     s0_old_cart,s0_new_cart,dx,dt,the_bc_level,mla)

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


    !**************************************************************************
    !     Create the new eta
    !**************************************************************************

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

    call destroy(bpt)

2000 format('... new min/max : density           ',e17.10,2x,e17.10)
2001 format('... new min/max : rho * H           ',e17.10,2x,e17.10)
2002 format('... new min/max : ',a16,2x,e17.10,2x,e17.10)
2003 format('... new min/max : tracer            ',e17.10,2x,e17.10)
2004 format(' ')

  end subroutine scalar_advance

end module scalar_advance_module
