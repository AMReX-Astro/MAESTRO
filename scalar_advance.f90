module scalar_advance_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: scalar_advance

contains

  subroutine scalar_advance(nlevs,mla,which_step,uold,sold,snew,sedge,thermal, &
                            umac,w0,w0mac,utrans,normal, &
                            rho0_old,rhoh0_old,rho0_new,rhoh0_new,p0_old,p0_new, &
                            tempbar,psi,dx,dt,the_bc_level)

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
    use geometry,      only: spherical, nr_fine
    use variables,     only: nscal, ntrac, trac_comp, temp_comp, &
                             rho_comp, rhoh_comp, foextrap_comp
    use probin_module, only: enthalpy_pred_type, use_thermal_diffusion, verbose
    use pred_parameters
    use modify_scal_force_module
    use convert_rhoX_to_X_module

    integer        , intent(in   ) :: nlevs
    type(ml_layout), intent(inout) :: mla
    integer        , intent(in   ) :: which_step
    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: sedge(:,:)
    type(multifab) , intent(in   ) :: thermal(:)
    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
    type(multifab) , intent(in   ) :: utrans(:,:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) :: rho0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: rhoh0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_old(:,0:)
    real(kind=dp_t), intent(in   ) :: p0_new(:,0:)
    real(kind=dp_t), intent(in   ) :: tempbar(:,0:)
    real(kind=dp_t), intent(in   ) :: psi(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    type(multifab) :: scal_force(nlevs)
    type(multifab) :: rho0_old_cart(nlevs)
    type(multifab) :: rho0_new_cart(nlevs)
    type(multifab) :: rhoh0_old_cart(nlevs)
    type(multifab) :: rhoh0_new_cart(nlevs)
    type(multifab) :: p0_new_cart(nlevs)
    type(multifab) :: sflux(nlevs,mla%dim)

    integer    :: velpred,comp,pred_comp,n,dm
    logical    :: umac_nodal_flag(sold(1)%dim), is_vel
    real(dp_t) :: smin,smax
    logical    :: is_prediction

    real(kind=dp_t), allocatable :: rho0_edge_old(:,:)
    real(kind=dp_t), allocatable :: rho0_edge_new(:,:)
    real(kind=dp_t), allocatable :: rhoh0_edge_old(:,:)
    real(kind=dp_t), allocatable :: rhoh0_edge_new(:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "scalar_advance")

    dm      = sold(1)%dim
    is_vel  = .false.
    velpred = 0    

    ! Create edge-centered base state quantities.
    ! Note: rho0_edge_{old,new} and rhoh0_edge_{old,new}
    ! contain edge-centered quantities created via spatial interpolation.
    allocate(rho0_edge_old (nlevs,0:nr_fine))
    allocate(rho0_edge_new (nlevs,0:nr_fine))
    allocate(rhoh0_edge_old(nlevs,0:nr_fine))
    allocate(rhoh0_edge_new(nlevs,0:nr_fine))

    do n = 1, nlevs
       call cell_to_edge(n,rho0_old(n,:),rho0_edge_old(n,:))
       call cell_to_edge(n,rho0_new(n,:),rho0_edge_new(n,:))
       call cell_to_edge(n,rhoh0_old(n,:),rhoh0_edge_old(n,:))
       call cell_to_edge(n,rhoh0_new(n,:),rhoh0_edge_new(n,:))
    end do

    ! Define rho0_old_cart and rho0_new_cart
    ! Define rhoh0_old_cart and rhoh0_new_cart
    if (spherical .eq. 1) then
       do n=1,nlevs
          call build(rho0_old_cart(n), sold(n)%la, 1, 1)
          call build(rho0_new_cart(n), sold(n)%la, 1, 1)
          call build(rhoh0_old_cart(n), sold(n)%la, 1, 1)
          call build(rhoh0_new_cart(n), sold(n)%la, 1, 1)
          call build(p0_new_cart(n), sold(n)%la, 1, 1)          
       end do

       call put_1d_array_on_cart(nlevs,rho0_old,rho0_old_cart,dm+rho_comp,.false., &
                                 .false.,dx,the_bc_level,mla)
       call put_1d_array_on_cart(nlevs,rho0_new,rho0_new_cart,dm+rho_comp,.false., &
                                 .false.,dx,the_bc_level,mla)
       call put_1d_array_on_cart(nlevs,p0_new,p0_new_cart,foextrap_comp,.false., &
                                 .false.,dx,the_bc_level,mla)


       if (enthalpy_pred_type .eq. predict_T_then_rhohprime .or. &
           enthalpy_pred_type .eq. predict_rhohprime) then
          call put_1d_array_on_cart(nlevs,rhoh0_old,rhoh0_old_cart,dm+rhoh_comp,.false., &
                                    .false.,dx,the_bc_level,mla)
          call put_1d_array_on_cart(nlevs,rhoh0_new,rhoh0_new_cart,dm+rhoh_comp,.false., &
                                    .false.,dx,the_bc_level,mla)
       end if
    end if

    ! This can be uncommented if you wish to compute T
    ! note -- not sure if p0_old or p0_new should be used here
    ! call makeTfromRhoH(nlevs,sold,p0_old(:,:),tempbar(:,:),mla,the_bc_level,dx)

    ! if we are predicting h on the edges, then convert the state arrays
    ! (and base state) from (rho h) to h.  Note, only the time-level n
    ! stuff need be converted, since that's all the prediction uses

    if (enthalpy_pred_type .eq. predict_h) then
       call convert_rhoh_to_h(nlevs,sold,.true.,mla,the_bc_level)
    end if

    !**************************************************************************
    ! Create scalar source term at time n
    !**************************************************************************

    do n = 1, nlevs
       call build(scal_force(n), sold(n)%la, nscal, 1)       
       call setval(scal_force(n),ZERO,all=.true.)
    end do

    ! make force for either h, T, or (rho h)'
    if (enthalpy_pred_type .eq. predict_rhohprime) then
       
       ! make force for (rho h)'
       is_prediction = .true.
       call mkrhohforce(nlevs,scal_force,is_prediction, &
                        thermal,umac,p0_old,p0_old,rho0_old,rho0_old,&
                        psi,normal,dx,.true.,mla,the_bc_level)

       call modify_scal_force(nlevs,scal_force,sold,umac,rhoh0_old, &
                              rhoh0_edge_old,w0,dx,rhoh0_old_cart,rhoh_comp,mla,the_bc_level)

    else if (enthalpy_pred_type .eq. predict_h) then

       ! make force for h by calling mkrhohforce then dividing by rho
       is_prediction = .true.
       call mkrhohforce(nlevs,scal_force,is_prediction, &
                        thermal,umac,p0_old,p0_old,rho0_old,rho0_old,&
                        psi,normal,dx,.true.,mla,the_bc_level)
       do n=1,nlevs
          call multifab_div_div_c(scal_force(n),rhoh_comp,sold(n),rho_comp,1,1)
       end do

    else if ( (enthalpy_pred_type .eq. predict_T_then_rhohprime) .or. &
              (enthalpy_pred_type .eq. predict_T_then_h        ) ) then

       ! make force for temperature
       call mktempforce(nlevs,scal_force,umac,sold,thermal,p0_old,p0_old,psi, &
                        normal,dx,mla,the_bc_level)

    end if
        
      
    !**************************************************************************
    !     Add w0 to MAC velocities (trans velocities already have w0).
    !**************************************************************************

    call addw0(nlevs,umac,w0,w0mac,mult=ONE)

    !**************************************************************************
    !     Create the edge states of (rho h)' or h or T and (rho X)' or X and rho'
    !**************************************************************************

    if (enthalpy_pred_type .eq. predict_rhohprime) then
       ! convert (rho h) -> (rho h)'
       call put_in_pert_form(nlevs,sold,rhoh0_old,dx,rhoh_comp, &
                             .true.,mla,the_bc_level)
    end if

    ! predict either T, h, or (rho h)' at the edges
    if ( (enthalpy_pred_type .eq. predict_T_then_rhohprime) .or. &
         (enthalpy_pred_type .eq. predict_T_then_h        )  ) then
       pred_comp = temp_comp
    else
       pred_comp = rhoh_comp
    end if

    call make_edge_scal(nlevs,sold,sedge,umac,scal_force,normal, &
                        w0,w0mac, &
                        dx,dt,is_vel,the_bc_level, &
                        pred_comp,dm+pred_comp,1,.false.,mla)

    if (enthalpy_pred_type .eq. predict_rhohprime) then
       ! convert (rho h)' -> (rho h)
       call put_in_pert_form(nlevs,sold,rhoh0_old,dx,rhoh_comp, &
                             .false.,mla,the_bc_level)
    end if

    ! if we are predicting h at the edges, then restore the state arrays
    ! (and base state) from h to (rho h)
    if (enthalpy_pred_type .eq. predict_h) then
       call convert_rhoh_to_h(nlevs,sold,.false.,mla,the_bc_level)
    end if

    ! Compute enthalpy edge states if we were predicting temperature.  This
    ! needs to be done after the state was returned to the full state.

    if ( (enthalpy_pred_type .eq. predict_T_then_rhohprime) .or. &
         (enthalpy_pred_type .eq. predict_T_then_h        ) ) then
       call makeRhoHfromT(nlevs,uold,sedge,rho0_old,rhoh0_old, &
                          rho0_edge_old,rhoh0_edge_old, &
                          rho0_new,rhoh0_new, &
                          rho0_edge_new,rhoh0_edge_new,the_bc_level,dx)
    end if

    !**************************************************************************
    !     Create the edge states of tracers.
    !**************************************************************************

    if (ntrac .ge. 1) then
       call make_edge_state(nlevs,sold,uold,sedge,umac,utrans,scal_force, &
                            normal,w0,w0mac, &
                            dx,dt,is_vel,the_bc_level,velpred, &
                            trac_comp,dm+trac_comp,ntrac,mla)
    end if

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

       ! compute enthalpy fluxes
       call mkflux(nlevs,sflux,sold,sedge,umac,w0,w0mac, &
                   rho0_old,rho0_edge_old,rho0_old_cart, &
                   rho0_old,rho0_edge_old,rho0_old_cart, &
                   rhoh0_old,rhoh0_edge_old,rhoh0_old_cart, &
                   rhoh0_old,rhoh0_edge_old,rhoh0_old_cart, &
                   rhoh_comp,rhoh_comp,mla)

       if (ntrac .ge. 1) then
          ! compute tracer fluxes
          call mkflux(nlevs,sflux,sold,sedge,umac,w0,w0mac, &
                      rho0_old,rho0_edge_old,rho0_old_cart, &
                      rho0_old,rho0_edge_old,rho0_old_cart, &
                      rhoh0_old,rhoh0_edge_old,rhoh0_old_cart, &
                      rhoh0_old,rhoh0_edge_old,rhoh0_old_cart, &
                      trac_comp,trac_comp+ntrac-1,mla)
       end if

    else if (which_step .eq. 2) then

       ! compute enthalpy fluxes
       call mkflux(nlevs,sflux,sold,sedge,umac,w0,w0mac, &
                   rho0_old,rho0_edge_old,rho0_old_cart, &
                   rho0_new,rho0_edge_new,rho0_new_cart, &
                   rhoh0_old,rhoh0_edge_old,rhoh0_old_cart, &
                   rhoh0_new,rhoh0_edge_new,rhoh0_new_cart, &
                   rhoh_comp,rhoh_comp,mla)

       if (ntrac .ge. 1) then
          ! compute tracer fluxes
          call mkflux(nlevs,sflux,sold,sedge,umac,w0,w0mac, &
                      rho0_old,rho0_edge_old,rho0_old_cart, &
                      rho0_new,rho0_edge_new,rho0_new_cart, &
                      rhoh0_old,rhoh0_edge_old,rhoh0_old_cart, &
                      rhoh0_new,rhoh0_edge_new,rhoh0_new_cart, &
                      trac_comp,trac_comp+ntrac-1,mla)
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
    
    !**************************************************************************
    !     Update tracers with convective differencing.
    !**************************************************************************
    
    if ( ntrac .ge. 1 ) then
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
    end if

    !**************************************************************************
    !     1) Create (rho h)' force at time n+1/2.
    !          (NOTE: we don't worry about filling ghost cells of the scal_force
    !                 because we only need them in valid regions...)     
    !     2) Update (rho h) with conservative differencing.
    !**************************************************************************
       
    if (which_step .eq. 1) then
      ! Here just send p0_old and p0_old
       is_prediction = .false.
       call mkrhohforce(nlevs,scal_force,is_prediction, &
                        thermal,umac,p0_old,p0_old,rho0_old,rho0_old,&
                        psi,normal,dx,.false.,mla,the_bc_level)
    else
      ! Here send p0_old and p0_new
       is_prediction = .false.
       call mkrhohforce(nlevs,scal_force,is_prediction, &
                        thermal,umac,p0_old,p0_new,rho0_old,rho0_new,&
                        psi,normal,dx,.false.,mla,the_bc_level)
    end if

    call update_scal(mla,rhoh_comp,rhoh_comp,sold,snew,sflux,scal_force, &
                     p0_new,p0_new_cart,dx,dt,the_bc_level)

    if ( verbose .ge. 1 ) then
       do n=1,nlevs
          smin = multifab_min_c(snew(n),rhoh_comp) 
          smax = multifab_max_c(snew(n),rhoh_comp)
          if (parallel_IOProcessor()) then
             write(6,2001) smin,smax
          end if
       end do
    end if

    if (parallel_IOProcessor()) write(6,2004) 

    if (spherical .eq. 1) then
       do n=1,nlevs
          call destroy(rho0_old_cart(n))
          call destroy(rho0_new_cart(n))
          call destroy(rhoh0_old_cart(n))
          call destroy(rhoh0_new_cart(n))
          call destroy(p0_new_cart(n))
       end do
    end if

    do n = 1, nlevs
       call destroy(scal_force(n))
       do comp = 1,dm
          call destroy(sflux(n,comp))
       end do
    end do

    if (.not. use_thermal_diffusion) then
       call makeTfromRhoH(nlevs,snew,p0_new(:,:),tempbar(:,:),mla,the_bc_level,dx)
    else
       do n=1,nlevs
          call multifab_copy_c(snew(n),temp_comp,sold(n),temp_comp,1)
       end do
    end if

    call destroy(bpt)

2001 format('... new min/max : rho * H           ',e17.10,2x,e17.10)
2003 format('... new min/max : tracer            ',e17.10,2x,e17.10)
2004 format(' ')

  end subroutine scalar_advance

end module scalar_advance_module
