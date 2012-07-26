module enthalpy_advance_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: enthalpy_advance

contains

  subroutine enthalpy_advance(mla,which_step,uold,sold,snew,sedge,sflux,scal_force,intra,&
                              thermal,umac,w0,w0mac, &
                              rho0_old,rhoh0_old,rho0_new,rhoh0_new,p0_old,p0_new, &
                              tempbar,psi,dx,dt,the_bc_level)

    use bl_prof_module
    use bl_constants_module
    use make_edge_scal_module
    use bds_module
    use mkflux_module
    use mkscalforce_module
    use update_scal_module
    use addw0_module
    use define_bc_module
    use fill_3d_module
    use pert_form_module
    use cell_to_edge_module
    use rhoh_vs_t_module
    use geometry,      only: spherical, nr_fine, r_start_coord, r_end_coord, &
         numdisjointchunks, nlevs_radial
    use variables,     only: nscal, temp_comp, rho_comp, rhoh_comp, foextrap_comp
    use probin_module, only: enthalpy_pred_type, use_thermal_diffusion, verbose, bds_type
    use pred_parameters
    use modify_scal_force_module
    use convert_rhoX_to_X_module
    use multifab_physbc_module

    type(ml_layout), intent(inout) :: mla
    integer        , intent(in   ) :: which_step
    type(multifab) , intent(in   ) :: uold(:)
    type(multifab) , intent(inout) :: sold(:)
    type(multifab) , intent(inout) :: snew(:)
    type(multifab) , intent(inout) :: sedge(:,:)
    type(multifab) , intent(inout) :: sflux(:,:)
    type(multifab) , intent(inout) :: scal_force(:)
    type(multifab) , intent(inout) :: intra(:)
    type(multifab) , intent(in   ) :: thermal(:)
    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    type(multifab) , intent(in   ) :: w0mac(:,:)
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

    type(multifab) :: rhoh0_old_cart(mla%nlevel)
    type(multifab) :: p0_new_cart(mla%nlevel)

    type(multifab) :: rho0mac_old(mla%nlevel,mla%dim)
    type(multifab) :: rho0mac_new(mla%nlevel,mla%dim)
    type(multifab) :: rhoh0mac_old(mla%nlevel,mla%dim)
    type(multifab) :: rhoh0mac_new(mla%nlevel,mla%dim)
    type(multifab) :: h0mac_old(mla%nlevel,mla%dim)
    type(multifab) :: h0mac_new(mla%nlevel,mla%dim)

    integer    :: pred_comp,n,r,i,comp,dm,nlevs
    logical    :: is_vel
    real(dp_t) :: smin,smax
    logical    :: is_prediction

    ! Create cell-centered base state quantity
    real(kind=dp_t), allocatable :: h0_old(:,:)
    real(kind=dp_t), allocatable :: h0_new(:,:)

    ! Create edge-centered base state quantities.
    ! Note: rho0_edge_{old,new} and rhoh0_edge_{old,new}
    ! contain edge-centered quantities created via spatial interpolation.
    real(kind=dp_t), allocatable ::  rho0_edge_old(:,:)
    real(kind=dp_t), allocatable ::  rho0_edge_new(:,:)
    real(kind=dp_t), allocatable :: rhoh0_edge_old(:,:)
    real(kind=dp_t), allocatable :: rhoh0_edge_new(:,:)
    real(kind=dp_t), allocatable ::    t0_edge_old(:,:)
    real(kind=dp_t), allocatable ::    t0_edge_new(:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "enthalpy_advance")

    allocate( h0_old(nlevs_radial,0:nr_fine-1))
    allocate( h0_new(nlevs_radial,0:nr_fine-1))
    allocate(  rho0_edge_old(nlevs_radial,0:nr_fine))
    allocate(  rho0_edge_new(nlevs_radial,0:nr_fine))
    allocate( rhoh0_edge_old(nlevs_radial,0:nr_fine))
    allocate( rhoh0_edge_new(nlevs_radial,0:nr_fine))
    allocate(    t0_edge_old(nlevs_radial,0:nr_fine))
    allocate(    t0_edge_new(nlevs_radial,0:nr_fine))

    is_vel  = .false.
    dm = mla%dim
    nlevs = mla%nlevel

    if (spherical .eq. 0) then
       call cell_to_edge(rho0_old,rho0_edge_old)
       call cell_to_edge(rho0_new,rho0_edge_new)
       call cell_to_edge(rhoh0_old,rhoh0_edge_old)
       call cell_to_edge(rhoh0_new,rhoh0_edge_new)
       call cell_to_edge(tempbar,t0_edge_old)
       call cell_to_edge(tempbar,t0_edge_new)
    end if

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

       do n=1,nlevs
          call build(rhoh0_old_cart(n), get_layout(sold(n)), 1, 1)
       end do

       call put_1d_array_on_cart(rhoh0_old,rhoh0_old_cart,dm+rhoh_comp,.false., &
                                 .false.,dx,the_bc_level,mla)

       call modify_scal_force(scal_force,sold,umac,rhoh0_old, &
                              rhoh0_edge_old,w0,dx,rhoh0_old_cart,rhoh_comp,mla,the_bc_level)

       do n=1,nlevs
          call destroy(rhoh0_old_cart(n))
       end do

    else if (enthalpy_pred_type .eq. predict_h .or. &
             enthalpy_pred_type .eq. predict_rhoh) then

       is_prediction = .true.
       call mkrhohforce(mla,scal_force,is_prediction,&
                        thermal,umac,p0_old,p0_old,rho0_old,rho0_old,&
                        psi,dx,.true.,the_bc_level)

       if (enthalpy_pred_type .eq. predict_h) then
          ! make force for h by calling mkrhohforce then dividing by rho
          do n=1,nlevs
             call multifab_div_div_c(scal_force(n),rhoh_comp,sold(n),rho_comp,1,1)
          end do
       end if

    else if (enthalpy_pred_type .eq. predict_hprime) then

       ! first compute h0_old
       do n=1,nlevs_radial
          do i=1,numdisjointchunks(n)
             do r=r_start_coord(n,i),r_end_coord(n,i)
                h0_old(n,r) = rhoh0_old(n,r) / rho0_old(n,r)
             end do
          end do
       end do

       ! make force for hprime
       is_prediction = .true.
       call mkhprimeforce(mla,sold,sold,scal_force,is_prediction,thermal,umac,p0_old, &
                          p0_old,h0_old,h0_old,psi,dx,.true.,the_bc_level)

    else if ( (enthalpy_pred_type .eq. predict_T_then_rhohprime) .or. &
              (enthalpy_pred_type .eq. predict_T_then_h        ) .or. &
              (enthalpy_pred_type .eq. predict_Tprime_then_h) ) then

       ! make force for temperature
       call mktempforce(mla,scal_force,umac,sold,thermal,p0_old,p0_old,psi,dx,the_bc_level)

    end if        

    ! reaction forcing terms
    do n=1,nlevs
       call multifab_plus_plus_c(scal_force(n), rhoh_comp, intra(n), rhoh_comp, 1, 0)
       call multifab_fill_boundary(scal_force(n))
       call multifab_physbc(scal_force(n), rhoh_comp, foextrap_comp, 1, the_bc_level(n))
    end do
      
    !**************************************************************************
    !     Add w0 to MAC velocities (trans velocities already have w0).
    !**************************************************************************

    call addw0(umac,the_bc_level,mla,w0,w0mac,mult=ONE)

    !**************************************************************************
    !     Create the edge states of (rho h)' or h or T    
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

    if (enthalpy_pred_type .eq. predict_rhoh) then
       if (bds_type .eq. 0) then
          call make_edge_scal(sold,sedge,umac,scal_force, &
                              dx,dt,is_vel,the_bc_level, &
                              pred_comp,dm+pred_comp,1,.true.,mla)
       else if (bds_type .eq. 1) then
          call bds(sold,sedge,umac,scal_force, &
                   dx,dt,is_vel,the_bc_level, &
                   pred_comp,dm+pred_comp,1,.true.,mla)
       end if
    else
       if (bds_type .eq. 0) then
          call make_edge_scal(sold,sedge,umac,scal_force, &
               dx,dt,is_vel,the_bc_level, &
               pred_comp,dm+pred_comp,1,.false.,mla)
       else if (bds_type .eq. 1) then
          call bds(sold,sedge,umac,scal_force, &
               dx,dt,is_vel,the_bc_level, &
               pred_comp,dm+pred_comp,1,.false.,mla)
       end if

    end if

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
       call makeHfromRhoT_edge(uold,sedge,rho0_old,rhoh0_old,tempbar,rho0_edge_old, &
                               rhoh0_edge_old,t0_edge_old,rho0_new,rhoh0_new,tempbar, &
                               rho0_edge_new,rhoh0_edge_new,t0_edge_new,the_bc_level,dx,mla)
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
                call multifab_build_edge(rhoh0mac_old(n,comp),mla%la(n),1,1,comp)
                call multifab_build_edge(h0mac_old(n,comp),mla%la(n),1,1,comp)
             end do
          end do

          do n=1,nlevs_radial
             do i=1,numdisjointchunks(n)
                do r=r_start_coord(n,i),r_end_coord(n,i)
                   h0_old(n,r) = rhoh0_old(n,r) / rho0_old(n,r)
                end do
             end do
          end do

          call make_s0mac(mla,rho0_old,rho0mac_old,dx,dm+rho_comp,the_bc_level)
          call make_s0mac(mla,rhoh0_old,rhoh0mac_old,dx,dm+rhoh_comp,the_bc_level)
          call make_s0mac(mla,h0_old,h0mac_old,dx,foextrap_comp,the_bc_level)
       end if

       ! compute enthalpy fluxes
       call mk_rhoh_flux(mla,sflux,sold,sedge,umac,w0,w0mac, &
                         rho0_old,rho0_edge_old,rho0mac_old, &
                         rho0_old,rho0_edge_old,rho0mac_old, &
                         rhoh0_old,rhoh0_edge_old,rhoh0mac_old, &
                         rhoh0_old,rhoh0_edge_old,rhoh0mac_old, &
                         h0mac_old,h0mac_old)

      if (spherical .eq. 1) then
          do n=1,nlevs
             do comp=1,dm
                call destroy(rho0mac_old(n,comp))
                call destroy(rhoh0mac_old(n,comp))
                call destroy(h0mac_old(n,comp))
             end do
          end do
       end if

    else if (which_step .eq. 2) then

       if (spherical .eq. 1) then
          do n=1,nlevs
             do comp=1,dm
                call multifab_build_edge( rho0mac_old(n,comp),mla%la(n),1,1,comp)
                call multifab_build_edge(rhoh0mac_old(n,comp),mla%la(n),1,1,comp)
                call multifab_build_edge(   h0mac_old(n,comp),mla%la(n),1,1,comp)
                call multifab_build_edge( rho0mac_new(n,comp),mla%la(n),1,1,comp)
                call multifab_build_edge(rhoh0mac_new(n,comp),mla%la(n),1,1,comp)
                call multifab_build_edge(   h0mac_new(n,comp),mla%la(n),1,1,comp)
             end do
          end do

          do n=1,nlevs_radial
             do i=1,numdisjointchunks(n)
                do r=r_start_coord(n,i),r_end_coord(n,i)
                   h0_old(n,r) = rhoh0_old(n,r) / rho0_old(n,r)
                   h0_new(n,r) = rhoh0_new(n,r) / rho0_new(n,r)
                end do
             end do
          end do

          call make_s0mac(mla,rho0_old,rho0mac_old,dx,dm+rho_comp,the_bc_level)
          call make_s0mac(mla,rhoh0_old,rhoh0mac_old,dx,dm+rhoh_comp,the_bc_level)
          call make_s0mac(mla,h0_old,h0mac_old,dx,foextrap_comp,the_bc_level)
          call make_s0mac(mla,rho0_new,rho0mac_new,dx,dm+rho_comp,the_bc_level)
          call make_s0mac(mla,rhoh0_new,rhoh0mac_new,dx,dm+rhoh_comp,the_bc_level)
          call make_s0mac(mla,h0_new,h0mac_new,dx,foextrap_comp,the_bc_level)
       end if

       ! compute enthalpy fluxes
       call mk_rhoh_flux(mla,sflux,sold,sedge,umac,w0,w0mac, &
                         rho0_old,rho0_edge_old,rho0mac_old, &
                         rho0_new,rho0_edge_new,rho0mac_new, &
                         rhoh0_old,rhoh0_edge_old,rhoh0mac_old, &
                         rhoh0_new,rhoh0_edge_new,rhoh0mac_new, &
                         h0mac_old,h0mac_new)

      if (spherical .eq. 1) then
          do n=1,nlevs
             do comp=1,dm
                call destroy(rho0mac_old(n,comp))
                call destroy(rhoh0mac_old(n,comp))
                call destroy(h0mac_old(n,comp))
                call destroy(rho0mac_new(n,comp))
                call destroy(rhoh0mac_new(n,comp))
                call destroy(h0mac_new(n,comp))
             end do
          end do
       end if

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

    ! reaction forcing terms
    do n=1,nlevs
       call multifab_plus_plus_c(scal_force(n), rhoh_comp, intra(n), rhoh_comp, 1, 0)
    end do

    if (spherical .eq. 1) then
       do n=1,nlevs
          call build(p0_new_cart(n), get_layout(sold(n)), 1, 1)
       end do

       call put_1d_array_on_cart(p0_new,p0_new_cart,foextrap_comp,.false., &
                                 .false.,dx,the_bc_level,mla)
    end if

    call update_scal(mla,rhoh_comp,rhoh_comp,sold,snew,sflux,scal_force, &
                     p0_new,p0_new_cart,dx,dt,the_bc_level)

    if (spherical .eq. 1) then
       do n=1,nlevs
          call destroy(p0_new_cart(n))
       end do
    end if

    if ( verbose .ge. 1 ) then
       do n=1,nlevs
          smin = multifab_min_c(snew(n),rhoh_comp) 
          smax = multifab_max_c(snew(n),rhoh_comp)
          if (parallel_IOProcessor()) then
             write (6,1999) n
             write (6,2001) smin,smax
          end if
       end do
    end if

    if (parallel_IOProcessor()) write(6,2004)

    call destroy(bpt)

1999 format('... Level ', i1, ' update:')
2001 format('... new min/max : rho * H           ',e17.10,2x,e17.10)
2004 format(' ')

  end subroutine enthalpy_advance

end module enthalpy_advance_module
