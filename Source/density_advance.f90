module density_advance_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: density_advance

contains

  subroutine density_advance(nlevs,mla,which_step,sold,snew,sedge,&
                            umac,w0,w0mac,etarhoflux,normal, &
                            rho0_old,rho0_new,p0_new, &
                            rho0_predicted_edge,dx,dt,the_bc_level)

    use bl_prof_module
    use bl_constants_module
    use make_edge_scal_module
    use mk_rhoX_flux_module
    use update_scal_module
    use addw0_module
    use define_bc_module
    use fill_3d_module
    use pert_form_module
    use cell_to_edge_module
    use network,       only: nspec, spec_names
    use geometry,      only: spherical, nr_fine
    use variables,     only: nscal, spec_comp, rho_comp, foextrap_comp
    use probin_module, only: verbose
    use modify_scal_force_module
    use convert_rhoX_to_X_module

    integer        , intent(in   ) :: nlevs
    type(ml_layout), intent(inout) :: mla
    integer        , intent(in   ) :: which_step
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: sedge(:,:)
    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    type(multifab) , intent(inout) :: etarhoflux(:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_predicted_edge(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    type(multifab) :: scal_force(nlevs)
    type(multifab) :: rho0_old_cart(nlevs)
    type(multifab) :: rho0_new_cart(nlevs)
    type(multifab) :: p0_new_cart(nlevs)
    type(multifab) :: sflux(nlevs,mla%dim)

    integer    :: comp,n,dm
    logical    :: umac_nodal_flag(sold(1)%dim), is_vel
    real(dp_t) :: smin,smax

    real(kind=dp_t), allocatable :: rho0_edge_old(:,:)
    real(kind=dp_t), allocatable :: rho0_edge_new(:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "density_advance")

    dm      = sold(1)%dim
    is_vel  = .false.

    ! create edge-centered base state quantities.
    ! Note: rho0_edge_{old,new} 
    ! contains edge-centered quantities created via spatial interpolation.
    ! This is to be contrasted to rho0_predicted_edge which is the half-time
    ! edge state created in advect_base.
    allocate(rho0_edge_old (nlevs,0:nr_fine))
    allocate(rho0_edge_new (nlevs,0:nr_fine))

    do n = 1, nlevs
       call cell_to_edge(n,rho0_old(n,:),rho0_edge_old(n,:))
       call cell_to_edge(n,rho0_new(n,:),rho0_edge_new(n,:))
    end do

    ! Define rho0_old_cart and rho0_new_cart
    if (spherical .eq. 1) then
       do n=1,nlevs
          call build(rho0_old_cart(n), sold(n)%la, 1, 1)
          call build(rho0_new_cart(n), sold(n)%la, 1, 1)
          call build(p0_new_cart(n), sold(n)%la, 1, 1)          
       end do

       call put_1d_array_on_cart(nlevs,rho0_old,rho0_old_cart,dm+rho_comp,.false., &
                                 .false.,dx,the_bc_level,mla)
       call put_1d_array_on_cart(nlevs,rho0_new,rho0_new_cart,dm+rho_comp,.false., &
                                 .false.,dx,the_bc_level,mla)
       call put_1d_array_on_cart(nlevs,p0_new,p0_new_cart,foextrap_comp,.false., &
                                 .false.,dx,the_bc_level,mla)
    end if

    !**************************************************************************
    ! Create source terms at time n
    !**************************************************************************

    do n = 1, nlevs
       call build(scal_force(n), sold(n)%la, nscal, 1)       
       call setval(scal_force(n),ZERO,all=.true.)
    end do

    ! X force is zero - do nothing

    ! make force for rho'
    call modify_scal_force(nlevs,scal_force,sold,umac,rho0_old, &
                           rho0_edge_old,w0,dx,rho0_old_cart,rho_comp,mla,the_bc_level)

    !**************************************************************************
    !     Add w0 to MAC velocities (trans velocities already have w0).
    !**************************************************************************

    call addw0(nlevs,umac,w0,w0mac,mult=ONE)

    !**************************************************************************
    !     Create the edge states of (rho X)' or X and rho'
    !**************************************************************************

    ! convert (rho X) --> X in sold 
    call convert_rhoX_to_X(nlevs,sold,.true.,mla,the_bc_level)

    ! convert rho -> rho' in sold 
    call put_in_pert_form(nlevs,sold,rho0_old,dx,rho_comp,.true.,mla,the_bc_level)

    ! predict X at the edges
    call make_edge_scal(nlevs,sold,sedge,umac,scal_force,normal, &
                        w0,w0mac,dx,dt,is_vel,the_bc_level, &
                        spec_comp,dm+spec_comp,nspec,.false.,mla)

    ! predict rho' at the edges
    call make_edge_scal(nlevs,sold,sedge,umac,scal_force,normal, &
                        w0,w0mac,dx,dt,is_vel,the_bc_level, &
                        rho_comp,dm+rho_comp,1,.false.,mla)

    ! convert rho' -> rho in sold
    call put_in_pert_form(nlevs,sold,rho0_old,dx,rho_comp,.false.,mla,the_bc_level)

    ! convert X --> (rho X) in sold 
    call convert_rhoX_to_X(nlevs,sold,.false.,mla,the_bc_level)

    !**************************************************************************
    !     Subtract w0 from MAC velocities.
    !**************************************************************************

    call addw0(nlevs,umac,w0,w0mac,mult=-ONE)

    !**************************************************************************
    !     Compute fluxes
    !**************************************************************************

    do comp = 1,dm
       umac_nodal_flag = .false.
       umac_nodal_flag(comp) = .true.
       do n=1,nlevs
          call multifab_build(sflux(n,comp), mla%la(n), nscal, 0, nodal = umac_nodal_flag)
       end do
    end do

    ! for which_step .eq. 1, we pass in only the old base state quantities
    ! for which_step .eq. 2, we pass in the old and new for averaging within mkflux
    if (which_step .eq. 1) then

       ! compute species fluxes
       call mk_rhoX_flux(mla,sflux,etarhoflux,sold,sedge,umac,w0,w0mac, &
                         rho0_old,rho0_edge_old,rho0_old_cart, &
                         rho0_old,rho0_edge_old,rho0_old_cart, &
                         rho0_predicted_edge,spec_comp,spec_comp+nspec-1)

    else if (which_step .eq. 2) then

       ! compute species fluxes
       call mk_rhoX_flux(mla,sflux,etarhoflux,sold,sedge,umac,w0,w0mac, &
                         rho0_old,rho0_edge_old,rho0_old_cart, &
                         rho0_new,rho0_edge_new,rho0_new_cart, &
                         rho0_predicted_edge,spec_comp,spec_comp+nspec-1)

    end if

    !**************************************************************************
    !     1) Set force for (rho X)'_i at time n+1/2 = 0.
    !     2) Update (rho X)_i with conservative differencing.
    !     3) Define density as the sum of the (rho X)_i
    !**************************************************************************
    
    do n=1,nlevs
       call setval(scal_force(n),ZERO,all=.true.)
    end do

    call update_scal(mla,spec_comp,spec_comp+nspec-1,sold,snew,sflux,scal_force, &
                     p0_new,p0_new_cart,dx,dt,the_bc_level)
    
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

    if (parallel_IOProcessor()) write(6,2004) 

    if (spherical .eq. 1) then
       do n=1,nlevs
          call destroy(rho0_old_cart(n))
          call destroy(rho0_new_cart(n))
          call destroy(p0_new_cart(n))
       end do
    end if

    do n = 1, nlevs
       call destroy(scal_force(n))
       do comp = 1,dm
          call destroy(sflux(n,comp))
       end do
    end do

    call destroy(bpt)

2000 format('... new min/max : density           ',e17.10,2x,e17.10)
2002 format('... new min/max : ',a16,2x,e17.10,2x,e17.10)
2004 format(' ')

  end subroutine density_advance

end module density_advance_module
