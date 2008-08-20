module scalar_advance_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: scalar_advance

contains

  subroutine scalar_advance(mla,which_step,sold,snew,sedge,sflux,&
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
    use cell_to_edge_module
    use geometry,      only: spherical, nr_fine
    use variables,     only: nscal, ntrac, trac_comp, temp_comp, &
                             rho_comp, rhoh_comp, foextrap_comp
    use probin_module, only: verbose

    type(ml_layout), intent(inout) :: mla
    integer        , intent(in   ) :: which_step
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: sedge(:,:)
    type(multifab) , intent(inout) :: sflux(:,:)
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

    type(multifab) :: scal_force(mla%nlevel)
    type(multifab) :: rho0_old_cart(mla%nlevel)
    type(multifab) :: rho0_new_cart(mla%nlevel)
    type(multifab) :: p0_new_cart(mla%nlevel)

    integer    :: velpred,comp,pred_comp,n,dm,nlevs
    logical    :: is_vel
    real(dp_t) :: smin,smax
    logical    :: is_prediction

    real(kind=dp_t), allocatable :: rho0_edge_old(:,:)
    real(kind=dp_t), allocatable :: rho0_edge_new(:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "scalar_advance")

    dm    = mla%dim
    nlevs = mla%nlevel
    is_vel  = .false.
    velpred = 0    

    ! Create edge-centered base state quantities.
    ! Note: rho0_edge_{old,new}
    ! contain edge-centered quantities created via spatial interpolation.
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
    ! Create and zero scalar source term.
    !**************************************************************************

    do n = 1, nlevs
       call build(scal_force(n), sold(n)%la, nscal, 1)       
       call setval(scal_force(n),ZERO,all=.true.)
    end do

    !**************************************************************************
    !     Add w0 to MAC velocities (trans velocities already have w0).
    !**************************************************************************

    call addw0(nlevs,umac,w0,w0mac,mult=ONE)

    !**************************************************************************
    !     Create the edge states of tracers (is_conservative = false).
    !**************************************************************************

    call make_edge_scal(nlevs,sold,sedge,umac,scal_force,normal, &
                        w0,w0mac,dx,dt,is_vel,the_bc_level, &
                        trac_comp,dm+trac_comp,ntrac,.false.,mla)

    !**************************************************************************
    !     Subtract w0 from MAC velocities.
    !**************************************************************************

    call addw0(nlevs,umac,w0,w0mac,mult=-ONE)

    !**************************************************************************
    !     Compute tracer fluxes
    !**************************************************************************

    ! for which_step .eq. 1, we pass in only the old base state quantities
    ! for which_step .eq. 2, we pass in the old and new for averaging within mkflux
    if (which_step .eq. 1) then

       call mk_rhoX_flux(mla,sflux,etarhoflux,sold,sedge,umac,w0,w0mac, &
                         rho0_old,rho0_edge_old,rho0_old_cart, &
                         rho0_old,rho0_edge_old,rho0_old_cart, &
                         rho0_predicted_edge,trac_comp,trac_comp+ntrac-1)

    else if (which_step .eq. 2) then

       call mk_rhoX_flux(mla,sflux,etarhoflux,sold,sedge,umac,w0,w0mac, &
                         rho0_old,rho0_edge_old,rho0_old_cart, &
                         rho0_new,rho0_edge_new,rho0_new_cart, &
                         rho0_predicted_edge,trac_comp,trac_comp+ntrac-1)

    end if
    
    !**************************************************************************
    !     Update tracers with convective differencing.
    !**************************************************************************
    
    call update_scal(mla,trac_comp,trac_comp+ntrac-1,sold,snew,sflux,scal_force, &
                     p0_new,p0_new_cart,dx,dt,the_bc_level)

    if ( verbose .ge. 1 ) then
       do n=1,nlevs
          smin = multifab_min_c(snew(n),trac_comp) 
          smax = multifab_max_c(snew(n),trac_comp)
          if (parallel_IOProcessor()) &
               write(6,2003) smin,smax
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
    end do

    call destroy(bpt)

2003 format('... new min/max : tracer            ',e17.10,2x,e17.10)
2004 format(' ')

  end subroutine scalar_advance

end module scalar_advance_module
