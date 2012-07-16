module density_advance_module

  use parallel, only: parallel_IOProcessor
  use bl_types, only: dp_t
  use multifab_module, only: multifab, multifab_max_c, multifab_min_c, get_layout, &
                             build, destroy, setval, multifab_build_edge, &
                             multifab_div_div_c, multifab_mult_mult_c, multifab_copy_c, &
                             multifab_plus_plus_c
  use ml_layout_module, only: ml_layout
  use define_bc_module, only: bc_level

  implicit none

  private

  public :: density_advance

contains

  subroutine density_advance(mla,which_step,sold,snew,sedge,sflux,scal_force,intra,&
                             umac,w0,w0mac,etarhoflux, &
                             rho0_old,rho0_new,p0_dummy, &
                             rho0_predicted_edge,dx,dt,the_bc_level)

    use bl_prof_module, only: bl_prof_timer, build, destroy
    use bl_constants_module, only: ZERO, ONE
    use make_edge_scal_module, only: make_edge_scal
    use bds_module, only: bds
    use mkflux_module, only: mk_rhoX_flux
    use update_scal_module, only: update_scal
    use addw0_module, only: addw0
    use fill_3d_module, only: put_1d_array_on_cart, make_s0mac
    use pert_form_module, only: put_in_pert_form
    use cell_to_edge_module, only: cell_to_edge
    use network,       only: nspec, spec_names
    use geometry,      only: spherical, nr_fine, nlevs_radial
    use variables,     only: nscal, ntrac, spec_comp, rho_comp, trac_comp, foextrap_comp
    use probin_module, only: verbose, bds_type, species_pred_type
    use modify_scal_force_module, only: modify_scal_force
    use convert_rhoX_to_X_module, only: convert_rhoX_to_X
    use pred_parameters

    type(ml_layout), intent(inout) :: mla
    integer        , intent(in   ) :: which_step
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: sedge(:,:)
    type(multifab) , intent(inout) :: sflux(:,:)
    type(multifab) , intent(inout) :: scal_force(:)
    type(multifab) , intent(inout) :: intra(:)
    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    type(multifab) , intent(inout) :: etarhoflux(:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_dummy(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    type(multifab) :: rho0_old_cart(mla%nlevel)
    type(multifab) :: p0_dummy_cart(mla%nlevel)

    type(multifab) :: rho0mac_old(mla%nlevel,mla%dim)
    type(multifab) :: rho0mac_new(mla%nlevel,mla%dim)

    integer    :: comp,i,n,dm,nlevs
    logical    :: is_vel
    real(dp_t) :: smin,smax

    real(kind=dp_t) :: rho0_edge_old(nlevs_radial,0:nr_fine)
    real(kind=dp_t) :: rho0_edge_new(nlevs_radial,0:nr_fine)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "density_advance")

    is_vel  = .false.

    dm = mla%dim
    nlevs = mla%nlevel

    if (spherical .eq. 0) then

       ! create edge-centered base state quantities.
       ! Note: rho0_edge_{old,new} 
       ! contains edge-centered quantities created via spatial interpolation.
       ! This is to be contrasted to rho0_predicted_edge which is the half-time
       ! edge state created in advect_base.       
       call cell_to_edge(rho0_old,rho0_edge_old)
       call cell_to_edge(rho0_new,rho0_edge_new)

    end if

    !**************************************************************************
    ! Create source terms at time n
    !**************************************************************************

    ! Source terms for X and for tracers are zero - do nothing
    do n = 1, nlevs
       call setval(scal_force(n),ZERO,all=.true.)
    end do

    if (spherical .eq. 1) then
       do n=1,nlevs
          call build(rho0_old_cart(n), get_layout(sold(n)), 1, 1)
       end do
       call put_1d_array_on_cart(rho0_old,rho0_old_cart,dm+rho_comp,.false., &
                                 .false.,dx,the_bc_level,mla)
    end if

    ! ** density source term **

    ! Make source term for rho or rho' 
    if (species_pred_type == predict_rhoprime_and_X) then
       ! rho' source term
       !   . this is needed for pred_rhoprime_and_X
       call modify_scal_force(scal_force,sold,umac,rho0_old, &
                              rho0_edge_old,w0,dx,rho0_old_cart,rho_comp, &
                              mla,the_bc_level)

    else if (species_pred_type == predict_rho_and_X) then
       ! rho source term
       call modify_scal_force(scal_force,sold,umac,rho0_old, &
                              rho0_edge_old,w0,dx,rho0_old_cart,rho_comp, &
                              mla,the_bc_level,.true.)
    endif

    ! ** species source term **

    ! for species_pred_types predict_rhoprime_and_X and
    ! predict_rho_and_X, there is no force for X.

    ! for predict_rhoX, we are predicting (rho X)
    ! as a conservative equation, and there is no force.

    if (spherical .eq. 1) then
       do n=1,nlevs
          call destroy(rho0_old_cart(n))
       end do
    end if

    ! SDC HACK
    do n=1,nlevs
       call multifab_plus_plus_c(scal_force(n), spec_comp, intra(n), spec_comp, nspec, 1)
    end do

    !**************************************************************************
    !     Add w0 to MAC velocities (trans velocities already have w0).
    !**************************************************************************

    call addw0(umac,the_bc_level,mla,w0,w0mac,mult=ONE)

    !**************************************************************************
    !     Create the edge states of (rho X)' or X and rho'
    !**************************************************************************

    if ((species_pred_type == predict_rhoprime_and_X) .or. &
        (species_pred_type == predict_rho_and_X)) then
       ! we are predicting X to the edges, so convert the scalar
       ! data to those quantities

       ! convert (rho X) --> X in sold 
       call convert_rhoX_to_X(sold,.true.,mla,the_bc_level)
    endif

    if (species_pred_type == predict_rhoprime_and_X) then
       ! convert rho -> rho' in sold
       !   . this is needed for predict_rhoprime_and_X
       call put_in_pert_form(mla,sold,rho0_old,dx,rho_comp,foextrap_comp,.true.,the_bc_level)
    endif


    ! predict species at the edges -- note, either X or (rho X) will be
    ! predicted here, depending on species_pred_type

    if ((species_pred_type == predict_rhoprime_and_X) .or. &
        (species_pred_type == predict_rho_and_X)) then

       ! we are predicting X to the edges, using the advective form of
       ! the prediction
       if (bds_type .eq. 0) then
          call make_edge_scal(sold,sedge,umac,scal_force, &
                              dx,dt,is_vel,the_bc_level, &
                              spec_comp,dm+spec_comp,nspec,.false.,mla)
       else if (bds_type .eq. 1) then
          call bds(sold,sedge,umac,scal_force, &
                   dx,dt,is_vel,the_bc_level, &
                   spec_comp,dm+spec_comp,nspec,.false.,mla)
       end if

    else if (species_pred_type == predict_rhoX) then

       ! we are predicting (rho X) to the edges, using the
       ! conservative form of the prediction
       if (bds_type .eq. 0) then
          call make_edge_scal(sold,sedge,umac,scal_force, &
                              dx,dt,is_vel,the_bc_level, &
                              spec_comp,dm+spec_comp,nspec,.true.,mla)
       else if (bds_type .eq. 1) then
          call bds(sold,sedge,umac,scal_force, &
                   dx,dt,is_vel,the_bc_level, &
                   spec_comp,dm+spec_comp,nspec,.true.,mla)
       end if

    endif

    if (species_pred_type .eq. predict_rhoX) then
       ! compute rho = sum(rhoX) at the edges
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(sedge(n,i),rho_comp,sedge(n,i),spec_comp,1,0)
             do comp=2,nspec
                call multifab_plus_plus_c(sedge(n,i),rho_comp,sedge(n,i),spec_comp+comp-1,1,0)
             end do
          end do
       end do
    else
       ! predict rho or rho' at the edges (depending on species_pred_type)
       if (bds_type .eq. 0) then
          call make_edge_scal(sold,sedge,umac,scal_force, &
                              dx,dt,is_vel,the_bc_level, &
                              rho_comp,dm+rho_comp,1,.false.,mla)
       else if (bds_type .eq. 1) then
          call bds(sold,sedge,umac,scal_force, &
                   dx,dt,is_vel,the_bc_level, &
                   rho_comp,dm+rho_comp,1,.false.,mla)
       end if
    end if

    if (species_pred_type == predict_rhoprime_and_X) then
       ! convert rho' -> rho in sold 
       call put_in_pert_form(mla,sold,rho0_old,dx,rho_comp,dm+rho_comp,.false.,the_bc_level)
    endif

    if ((species_pred_type == predict_rhoprime_and_X) .or. &
        (species_pred_type == predict_rho_and_X)) then
       ! convert X --> (rho X) in sold 
       call convert_rhoX_to_X(sold,.false.,mla,the_bc_level)
    endif


    !**************************************************************************
    !     Create edge states of tracers
    !**************************************************************************
    if (ntrac .ge. 1) then
       if (bds_type .eq. 0) then
          call make_edge_scal(sold,sedge,umac,scal_force, &
                              dx,dt,is_vel,the_bc_level, &
                              trac_comp,dm+trac_comp,ntrac,.false.,mla)
       else if (bds_type .eq. 1) then
          call bds(sold,sedge,umac,scal_force, &
                   dx,dt,is_vel,the_bc_level, &
                   trac_comp,dm+trac_comp,ntrac,.false.,mla)
       end if
    end if

    !**************************************************************************
    !     Subtract w0 from MAC velocities.
    !**************************************************************************

    call addw0(umac,the_bc_level,mla,w0,w0mac,mult=-ONE)

    !**************************************************************************
    !     Compute fluxes
    !**************************************************************************

    ! for which_step .eq. 1, we pass in only the old base state quantities
    ! for which_step .eq. 2, we pass in the old and new for averaging within mkflux
    if (which_step .eq. 1) then

       if (spherical .eq. 1) then
          do n=1,nlevs
             do comp=1,dm
                call multifab_build_edge(rho0mac_old(n,comp),mla%la(n),1,1,comp) 
             end do
          end do

          call make_s0mac(mla,rho0_old,rho0mac_old,dx,dm+rho_comp,the_bc_level)

       end if

       ! compute species fluxes
       call mk_rhoX_flux(mla,sflux,etarhoflux,sold,sedge,umac,w0,w0mac, &
                         rho0_old,rho0_edge_old,rho0mac_old, &
                         rho0_old,rho0_edge_old,rho0mac_old, &
                         rho0_predicted_edge,spec_comp,spec_comp+nspec-1)

       ! compute tracer fluxes
       if (ntrac .ge. 1) then
          call mk_rhoX_flux(mla,sflux,etarhoflux,sold,sedge,umac,w0,w0mac, &
                            rho0_old,rho0_edge_old,rho0mac_old, &
                            rho0_old,rho0_edge_old,rho0mac_old, &
                            rho0_predicted_edge,trac_comp,trac_comp+ntrac-1)
       end if


       if (spherical .eq. 1) then
          do n=1,nlevs
             do comp=1,dm
                call destroy(rho0mac_old(n,comp))
             end do
          end do
       end if

    else if (which_step .eq. 2) then

       if (spherical .eq. 1) then
          do n=1,nlevs
             do comp=1,dm
                call multifab_build_edge(rho0mac_old(n,comp),mla%la(n),1,1,comp)
                call multifab_build_edge(rho0mac_new(n,comp),mla%la(n),1,1,comp)
             end do
          end do

          call make_s0mac(mla,rho0_old,rho0mac_old,dx,dm+rho_comp,the_bc_level)
          call make_s0mac(mla,rho0_new,rho0mac_new,dx,dm+rho_comp,the_bc_level)

       end if

       ! compute species fluxes
       call mk_rhoX_flux(mla,sflux,etarhoflux,sold,sedge,umac,w0,w0mac, &
                         rho0_old,rho0_edge_old,rho0mac_old, &
                         rho0_new,rho0_edge_new,rho0mac_new, &
                         rho0_predicted_edge,spec_comp,spec_comp+nspec-1)

       ! compute tracer fluxes
       if (ntrac .ge. 1) then
          call mk_rhoX_flux(mla,sflux,etarhoflux,sold,sedge,umac,w0,w0mac, &
                            rho0_old,rho0_edge_old,rho0mac_old, &
                            rho0_new,rho0_edge_new,rho0mac_new, &
                            rho0_predicted_edge,trac_comp,trac_comp+ntrac-1)
       end if

       if (spherical .eq. 1) then
          do n=1,nlevs
             do comp=1,dm
                call destroy(rho0mac_old(n,comp))
                call destroy(rho0mac_new(n,comp))
             end do
          end do
       end if

    end if

    !**************************************************************************
    !     1) Set force for (rho X)'_i at time n+1/2 = 0.
    !     2) Update (rho X)_i with conservative differencing.
    !     3) Define density as the sum of the (rho X)_i
    !     4) Update tracer with conservative differencing as well.
    !**************************************************************************
    
    do n=1,nlevs
       call setval(scal_force(n),ZERO,all=.true.)
    end do

    ! p0 only used in rhoh update so we just pass in a dummy version
    if (spherical .eq. 1) then
       do n=1,nlevs
          call build(p0_dummy_cart(n), get_layout(sold(n)), 1, 1)          
       end do
    end if

    call update_scal(mla,spec_comp,spec_comp+nspec-1,sold,snew,sflux,scal_force, &
                     p0_dummy,p0_dummy_cart,dx,dt,the_bc_level)

    if (ntrac .ge. 1) then
       call update_scal(mla,trac_comp,trac_comp+ntrac-1,sold,snew,sflux,scal_force, &
                        p0_dummy,p0_dummy_cart,dx,dt,the_bc_level)
    end if

    if (spherical .eq. 1) then
       do n=1,nlevs
          call destroy(p0_dummy_cart(n))
       end do
    end if
    
    if (verbose .ge. 1) then
       do n=1, nlevs
          if (parallel_IOProcessor()) write(6,1999) n

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

          if (ntrac .ge. 1) then
             smin = multifab_min_c(snew(n),trac_comp) 
             smax = multifab_max_c(snew(n),trac_comp)
             if (parallel_IOProcessor()) &
                  write(6,2003) smin,smax
          end if
       end do
    end if

    call destroy(bpt)

1999 format('... Level ', i1, ' update:')
2000 format('... new min/max : density           ',e17.10,2x,e17.10)
2002 format('... new min/max : ',a16,2x,e17.10,2x,e17.10)
2003 format('... new min/max : tracer            ',e17.10,2x,e17.10)

  end subroutine density_advance

end module density_advance_module
