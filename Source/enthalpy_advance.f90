module enthalpy_advance_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: enthalpy_advance

contains

  subroutine enthalpy_advance(mla,which_step,uold,sold,snew,sedge,sflux,scal_force,&
                              thermal,umac,w0,w0mac,normal, &
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
    use geometry,      only: spherical, nr_fine, dm, r_start_coord, r_end_coord, &
         numdisjointchunks
    use variables,     only: nscal, temp_comp, rho_comp, rhoh_comp, foextrap_comp
    use probin_module, only: enthalpy_pred_type, use_thermal_diffusion, edge_nodal_flag, &
         verbose, use_tfromp, nlevs
    use pred_parameters
    use modify_scal_force_module
    use convert_rhoX_to_X_module

    type(ml_layout), intent(inout) :: mla
    integer        , intent(in   ) :: which_step
    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: sedge(:,:)
    type(multifab) , intent(inout) :: sflux(:,:)
    type(multifab) , intent(inout) :: scal_force(:)
    type(multifab) , intent(in   ) :: thermal(:)
    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
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

    type(multifab) :: rho0_old_cart(mla%nlevel)
    type(multifab) :: rho0_new_cart(mla%nlevel)
    type(multifab) :: rhoh0_old_cart(mla%nlevel)
    type(multifab) :: rhoh0_new_cart(mla%nlevel)
    type(multifab) :: t0_old_cart(mla%nlevel)
    type(multifab) :: t0_new_cart(mla%nlevel)
    type(multifab) :: p0_new_cart(mla%nlevel)

    integer    :: pred_comp,n,r,i
    logical    :: is_vel
    real(dp_t) :: smin,smax
    logical    :: is_prediction

    real(kind=dp_t) :: h0_old(nlevs,0:nr_fine-1)

    ! Create edge-centered base state quantities.
    ! Note: rho0_edge_{old,new} and rhoh0_edge_{old,new}
    ! contain edge-centered quantities created via spatial interpolation.
    real(kind=dp_t) :: rho0_edge_old(nlevs,0:nr_fine)
    real(kind=dp_t) :: rho0_edge_new(nlevs,0:nr_fine)
    real(kind=dp_t) :: rhoh0_edge_old(nlevs,0:nr_fine)
    real(kind=dp_t) :: rhoh0_edge_new(nlevs,0:nr_fine)
    real(kind=dp_t) :: t0_edge_old(nlevs,0:nr_fine)
    real(kind=dp_t) :: t0_edge_new(nlevs,0:nr_fine)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "enthalpy_advance")

    is_vel  = .false.

    do n = 1, nlevs
       call cell_to_edge(n,rho0_old(n,:),rho0_edge_old(n,:))
       call cell_to_edge(n,rho0_new(n,:),rho0_edge_new(n,:))
       call cell_to_edge(n,rhoh0_old(n,:),rhoh0_edge_old(n,:))
       call cell_to_edge(n,rhoh0_new(n,:),rhoh0_edge_new(n,:))
       call cell_to_edge(n,tempbar(n,:),t0_edge_old(n,:))
       call cell_to_edge(n,tempbar(n,:),t0_edge_new(n,:))
    end do

    ! Define rho0_old_cart and rho0_new_cart
    ! Define rhoh0_old_cart and rhoh0_new_cart
    if (spherical .eq. 1) then
       do n=1,nlevs
          call build(rho0_old_cart(n), sold(n)%la, 1, 1)
          call build(rho0_new_cart(n), sold(n)%la, 1, 1)
          call build(rhoh0_old_cart(n), sold(n)%la, 1, 1)
          call build(rhoh0_new_cart(n), sold(n)%la, 1, 1)
          call build(t0_old_cart(n), sold(n)%la, 1, 1)
          call build(t0_new_cart(n), sold(n)%la, 1, 1)
          call build(p0_new_cart(n), sold(n)%la, 1, 1)          
       end do

       call put_1d_array_on_cart(rho0_old,rho0_old_cart,dm+rho_comp,.false., &
                                 .false.,dx,the_bc_level,mla)
       call put_1d_array_on_cart(rho0_new,rho0_new_cart,dm+rho_comp,.false., &
                                 .false.,dx,the_bc_level,mla)
       call put_1d_array_on_cart(p0_new,p0_new_cart,foextrap_comp,.false., &
                                 .false.,dx,the_bc_level,mla)

       if (enthalpy_pred_type .eq. predict_T_then_rhohprime .or. &
           enthalpy_pred_type .eq. predict_rhohprime .or. &
           enthalpy_pred_type .eq. predict_hprime) then
          call put_1d_array_on_cart(rhoh0_old,rhoh0_old_cart,dm+rhoh_comp,.false., &
                                    .false.,dx,the_bc_level,mla)
          call put_1d_array_on_cart(rhoh0_new,rhoh0_new_cart,dm+rhoh_comp,.false., &
                                    .false.,dx,the_bc_level,mla)
       end if

       if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
          call put_1d_array_on_cart(tempbar,t0_old_cart,dm+temp_comp,.false., &
                                    .false.,dx,the_bc_level,mla)
          call put_1d_array_on_cart(tempbar,t0_new_cart,dm+temp_comp,.false., &
                                    .false.,dx,the_bc_level,mla)
       end if

    end if

    ! This can be uncommented if you wish to compute T
    ! note -- not sure if p0_old or p0_new should be used here
    ! if (use_tfromp) then
    !    call makeTfromRhoP(sold,p0_old,tempbar,mla,the_bc_level,dx)
    ! else
    !    call makeTfromRhoH(sold,p0_old,tempbar,mla,the_bc_level,dx)
    ! end if

    if (enthalpy_pred_type .eq. predict_h .or. &
        enthalpy_pred_type .eq. predict_hprime) then
       ! convert (rho h) -> h
       call convert_rhoh_to_h(sold,.true.,mla,the_bc_level)
    end if

    !**************************************************************************
    ! Create scalar source term at time n
    !**************************************************************************

    do n = 1, nlevs
       call setval(scal_force(n),ZERO,all=.true.)
    end do

    ! compute forcing terms
    if (enthalpy_pred_type .eq. predict_rhohprime) then
       
       ! make force for (rho h)'
       is_prediction = .true.
       call mkrhohforce(mla,scal_force,is_prediction, &
                        thermal,umac,p0_old,p0_old,rho0_old,rho0_old,&
                        psi,dx,.true.,the_bc_level)

       call modify_scal_force(scal_force,sold,umac,rhoh0_old, &
                              rhoh0_edge_old,w0,dx,rhoh0_old_cart,rhoh_comp,mla,the_bc_level)

    else if (enthalpy_pred_type .eq. predict_h) then

       ! make force for h by calling mkrhohforce then dividing by rho
       is_prediction = .true.
       call mkrhohforce(mla,scal_force,is_prediction,&
                        thermal,umac,p0_old,p0_old,rho0_old,rho0_old,&
                        psi,dx,.true.,the_bc_level)
       do n=1,nlevs
          call multifab_div_div_c(scal_force(n),rhoh_comp,sold(n),rho_comp,1,1)
       end do

    else if (enthalpy_pred_type .eq. predict_hprime) then

       ! first compute h0_old
       do n=1,nlevs
          do i=1,numdisjointchunks(n)
             do r=r_start_coord(n,i),r_end_coord(n,i)
                h0_old(n,r) = rhoh0_old(n,r) / rho0_old(n,r)
             end do
          end do
       end do

       ! make force for hprime
       is_prediction = .true.
       call mkhprimeforce(mla,sold,sold,scal_force,is_prediction,&
                          thermal,umac,p0_old,p0_old,rho0_old,rho0_old,&
                          h0_old,h0_old,psi,dx,.true.,the_bc_level)

    else if ( (enthalpy_pred_type .eq. predict_T_then_rhohprime) .or. &
              (enthalpy_pred_type .eq. predict_T_then_h        ) .or. &
              (enthalpy_pred_type .eq. predict_Tprime_then_h) ) then

       ! make force for temperature
       call mktempforce(mla,scal_force,umac,sold,thermal,p0_old,p0_old, &
                        tempbar,tempbar,psi,dx,the_bc_level)

    end if        
      
    !**************************************************************************
    !     Add w0 to MAC velocities (trans velocities already have w0).
    !**************************************************************************

    call addw0(umac,w0,w0mac,mult=ONE)

    !**************************************************************************
    !     Create the edge states of (rho h)' or h or T and (rho X)' or X and rho'
    !**************************************************************************

    if (enthalpy_pred_type .eq. predict_rhohprime) then
       ! convert (rho h) -> (rho h)' or
       call put_in_pert_form(mla,sold,rhoh0_old,dx,rhoh_comp,foextrap_comp,.true., &
                             the_bc_level)
    end if

    if (enthalpy_pred_type .eq. predict_hprime) then
       ! convert h -> h'
       call put_in_pert_form(mla,sold,h0_old,dx,rhoh_comp,foextrap_comp,.true.,the_bc_level)
    end if

    if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
       ! convert T -> T'
       call put_in_pert_form(mla,sold,tempbar,dx,temp_comp,foextrap_comp,.true.,the_bc_level)
    end if

    ! predict either T, h, or (rho h)' at the edges
    if ( (enthalpy_pred_type .eq. predict_T_then_rhohprime) .or. &
         (enthalpy_pred_type .eq. predict_T_then_h        ) .or. &
         (enthalpy_pred_type .eq. predict_Tprime_then_h) ) then
       pred_comp = temp_comp
    else
       pred_comp = rhoh_comp
    end if

    call make_edge_scal(sold,sedge,umac,scal_force,normal, &
                        w0,w0mac, &
                        dx,dt,is_vel,the_bc_level, &
                        pred_comp,dm+pred_comp,1,.false.,mla)

    if (enthalpy_pred_type .eq. predict_rhohprime) then
       ! convert (rho h)' -> (rho h)
       call put_in_pert_form(mla,sold,rhoh0_old,dx,rhoh_comp,dm+rhoh_comp,.false., &
                             the_bc_level)
    end if

    if (enthalpy_pred_type .eq. predict_hprime) then
       ! convert h' -> h
       call put_in_pert_form(mla,sold,h0_old,dx,rhoh_comp,foextrap_comp,.false.,the_bc_level)
    end if

    if (enthalpy_pred_type .eq. predict_Tprime_then_h) then
       ! convert T' -> T
       call put_in_pert_form(mla,sold,tempbar,dx,temp_comp,dm+temp_comp,.false., &
                             the_bc_level)
    end if

    if (enthalpy_pred_type .eq. predict_h .or. &
        enthalpy_pred_type .eq. predict_hprime) then
       ! convert h -> (rho h)
       call convert_rhoh_to_h(sold,.false.,mla,the_bc_level)
    end if

    ! Compute enthalpy edge states if we were predicting temperature.  This
    ! needs to be done after the state was returned to the full state.
    if ( (enthalpy_pred_type .eq. predict_T_then_rhohprime) .or. &
         (enthalpy_pred_type .eq. predict_T_then_h        ) .or. &
         (enthalpy_pred_type .eq. predict_Tprime_then_h) ) then
       call makeHfromRhoT_edge(uold,sedge,rho0_old,rhoh0_old,tempbar, &
                               rho0_edge_old,rhoh0_edge_old,t0_edge_old, &
                               rho0_new,rhoh0_new,tempbar, &
                               rho0_edge_new,rhoh0_edge_new,t0_edge_new, &
                               the_bc_level,dx)
    end if

    !**************************************************************************
    !     Subtract w0 from MAC velocities.
    !**************************************************************************

    call addw0(umac,w0,w0mac,mult=-ONE)

    !**************************************************************************
    !     Compute fluxes
    !**************************************************************************

    ! for which_step .eq. 1, we pass in only the old base state quantities
    ! for which_step .eq. 2, we pass in the old and new for averaging within mkflux
    if (which_step .eq. 1) then

       ! compute enthalpy fluxes
       call mkflux(sflux,sold,sedge,umac,w0,w0mac, &
                   rho0_old,rho0_edge_old,rho0_old_cart, &
                   rho0_old,rho0_edge_old,rho0_old_cart, &
                   rhoh0_old,rhoh0_edge_old,rhoh0_old_cart, &
                   rhoh0_old,rhoh0_edge_old,rhoh0_old_cart, &
                   rhoh_comp,rhoh_comp,mla)

    else if (which_step .eq. 2) then

       ! compute enthalpy fluxes
       call mkflux(sflux,sold,sedge,umac,w0,w0mac, &
                   rho0_old,rho0_edge_old,rho0_old_cart, &
                   rho0_new,rho0_edge_new,rho0_new_cart, &
                   rhoh0_old,rhoh0_edge_old,rhoh0_old_cart, &
                   rhoh0_new,rhoh0_edge_new,rhoh0_new_cart, &
                   rhoh_comp,rhoh_comp,mla)

    end if

    do n=1,nlevs
       call setval(scal_force(n),ZERO,all=.true.)
    end do

    !**************************************************************************
    !     1) Create (rho h)' force at time n+1/2.
    !          (NOTE: we don't worry about filling ghost cells of the scal_force
    !                 because we only need them in valid regions...)     
    !     2) Update (rho h) with conservative differencing.
    !**************************************************************************
       
    if (which_step .eq. 1) then
      ! Here just send p0_old and p0_old
       is_prediction = .false.
       call mkrhohforce(mla,scal_force,is_prediction, &
                        thermal,umac,p0_old,p0_old,rho0_old,rho0_old,&
                        psi,dx,.false.,the_bc_level)
    else
      ! Here send p0_old and p0_new
       is_prediction = .false.
       call mkrhohforce(mla,scal_force,is_prediction, &
                        thermal,umac,p0_old,p0_new,rho0_old,rho0_new,&
                        psi,dx,.false.,the_bc_level)
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

    if (.not. use_thermal_diffusion) then
       if (use_tfromp) then
          call makeTfromRhoP(snew,p0_new,tempbar,mla,the_bc_level,dx)
       else
          call makeTfromRhoH(snew,p0_new,tempbar,mla,the_bc_level,dx)
       end if
    else
       do n=1,nlevs
          call multifab_copy_c(snew(n),temp_comp,sold(n),temp_comp,1)
       end do
    end if

    call destroy(bpt)

2001 format('... new min/max : rho * H           ',e17.10,2x,e17.10)
2004 format(' ')

  end subroutine enthalpy_advance

end module enthalpy_advance_module
