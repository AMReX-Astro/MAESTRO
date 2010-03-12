module make_plotfile_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_boxarray_module

  implicit none

  private
  public :: get_plot_names, make_plotfile

contains

  subroutine get_plot_names(plot_names)

    use plot_variables_module
    use variables
    use network, only: nspec, short_spec_names
    use probin_module, only: plot_spec, plot_trac, plot_base, use_thermal_diffusion
    use geometry, only: spherical, dm

    character(len=20), intent(inout) :: plot_names(:)

    ! Local variables
    integer :: comp

    plot_names(icomp_vel  ) = "x_vel"
    if (dm > 1) then
       plot_names(icomp_vel+1) = "y_vel"
    end if
    if (dm > 2) then
       plot_names(icomp_vel+2) = "z_vel"
    end if
    plot_names(icomp_rho)  = "density"
    plot_names(icomp_rhoh) = "rhoh"

    if (plot_spec) then
       do comp = 1, nspec
          plot_names(icomp_spec+comp-1) = "X(" // trim(short_spec_names(comp)) // ")"
       end do
    end if

    if (plot_trac) then
       do comp = 1, ntrac
          plot_names(icomp_trac+comp-1) = "tracer"
       end do
    end if

    if (plot_base) then
       plot_names(icomp_w0)   = "w0_x"
       if (dm > 1) plot_names(icomp_w0+1) = "w0_y"
       if (dm > 2) plot_names(icomp_w0+2) = "w0_z"
       plot_names(icomp_divw0) = "divw0"
       plot_names(icomp_rho0)  = "rho0"
       plot_names(icomp_rhoh0) = "rhoh0"
       plot_names(icomp_h0)    = "h0"
       plot_names(icomp_p0)    = "p0"
    end if

    if (spherical .eq. 1) then
       plot_names(icomp_velr) = "radial_velocity"
    endif

    plot_names(icomp_magvel)      = "magvel"
    plot_names(icomp_mom)         = "momentum"
    plot_names(icomp_vort)        = "vort"
    plot_names(icomp_divu)        = "divu"
    plot_names(icomp_enthalpy)    = "enthalpy"
    plot_names(icomp_rhopert)     = "rhopert"
    plot_names(icomp_rhohpert)    = "rhohpert"
    plot_names(icomp_tfromp)      = "tfromp"
    plot_names(icomp_tfromH)      = "tfromh"
    plot_names(icomp_tpert)       = "tpert"
    plot_names(icomp_machno)      = "Machnumber"
    plot_names(icomp_dp)          = "deltap"
    plot_names(icomp_dg)          = "deltagamma"
    plot_names(icomp_entropy)     = "entropy"
    plot_names(icomp_entropypert) = "entropypert"
    plot_names(icomp_dT)          = "deltaT"
    plot_names(icomp_sponge)      = "sponge"
    plot_names(icomp_gp)          = "gpx"
    if (dm > 1) plot_names(icomp_gp+1) = "gpy"
    if (dm > 2) plot_names(icomp_gp+2) = "gpz"

    if (plot_spec) then
       do comp = 1, nspec
          plot_names(icomp_omegadot+comp-1) = &
               "omegadot(" // trim(short_spec_names(comp)) // ")"
       end do
       plot_names(icomp_enuc) = "enucdot"
    end if

    if (use_thermal_diffusion) then
       plot_names(icomp_thermal) = "thermal"
       plot_names(icomp_conductivity) = "conductivity"
    endif

  end subroutine get_plot_names

  subroutine make_plotfile(dirname,mla,u,s,gpres,rho_omegadot,rho_Hnuc,thermal,Source,sponge,&
                           mba,plot_names,time,dx,the_bc_tower,w0,rho0,rhoh0,p0,tempbar, &
                           gamma1bar,normal)

    use bl_prof_module
    use fabio_module
    use variables
    use plot_variables_module
    use fill_3d_module
    use probin_module, only: nOutFiles, lUsingNFiles, plot_spec, plot_trac, plot_base, &
         single_prec_plotfiles, edge_nodal_flag, do_smallscale, use_thermal_diffusion, &
         evolve_base_state, prob_lo, prob_hi
    use geometry, only: spherical, nr_fine, dm, nlevs, nlevs_radial
    use average_module
    use ml_restriction_module
    use multifab_physbc_module
    use multifab_fill_ghost_module
    use bl_constants_module
    use network, only: nspec

    character(len=*) , intent(in   ) :: dirname
    type(ml_layout)  , intent(in   ) :: mla
    type(multifab)   , intent(in   ) :: u(:)
    type(multifab)   , intent(in   ) :: s(:)
    type(multifab)   , intent(in   ) :: gpres(:)
    type(multifab)   , intent(in   ) :: rho_omegadot(:)
    type(multifab)   , intent(in   ) :: rho_Hnuc(:)
    type(multifab)   , intent(in   ) :: thermal(:)
    type(multifab)   , intent(in   ) :: Source(:)
    type(multifab)   , intent(in   ) :: sponge(:)
    type(ml_boxarray), intent(in   ) :: mba
    character(len=20), intent(in   ) :: plot_names(:)
    real(dp_t)       , intent(in   ) :: time,dx(:,:)
    type(bc_tower)   , intent(in   ) :: the_bc_tower
    real(dp_t)       , intent(in   ) :: w0(:,0:)
    real(dp_t)       , intent(in   ) :: rho0(:,0:)
    real(dp_t)       , intent(in   ) :: rhoh0(:,0:)
    real(dp_t)       , intent(in   ) :: p0(:,0:)
    real(dp_t)       , intent(in   ) :: tempbar(:,0:)
    real(dp_t)       , intent(in   ) :: gamma1bar(:,0:)
    type(multifab)   , intent(in   ) :: normal(:)
    
    type(multifab) :: plotdata(nlevs)
    type(multifab) ::  tempfab(nlevs)
    type(multifab) ::    w0mac(nlevs,dm)
    type(multifab) :: w0r_cart(nlevs)

    real(dp_t) :: entropybar(nlevs_radial,0:nr_fine-1)
    real(dp_t) ::         h0(nlevs_radial,0:nr_fine-1)

    integer :: n,prec,comp

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_plotfile")

    if (single_prec_plotfiles) then
       prec = FABIO_SINGLE
    else
       prec = FABIO_DOUBLE
    endif

    do n = 1,nlevs

       call multifab_build(plotdata(n), mla%la(n), n_plot_comps, 0)
       call multifab_build(tempfab(n),  mla%la(n), dm,           1)
              
       ! VELOCITY 
       call multifab_copy_c(plotdata(n),icomp_vel,u(n),1,dm)

       ! DENSITY AND (RHO H) 
       call multifab_copy_c(plotdata(n),icomp_rho,s(n),rho_comp,2)
       
       if (plot_spec) then
          
          ! SPECIES
          call multifab_copy_c(plotdata(n),icomp_spec,s(n),spec_comp,nspec)
          do comp=1,nspec
             call multifab_div_div_c(plotdata(n),icomp_spec+comp-1,s(n),rho_comp,1)
          end do

          ! OMEGADOT
          call multifab_copy_c(plotdata(n),icomp_omegadot,rho_omegadot(n),1,nspec)
          do comp=1,nspec
             call multifab_div_div_c(plotdata(n),icomp_omegadot+comp-1,s(n),rho_comp,1)
          end do

          ! ENUCDOT
          call multifab_copy_c(plotdata(n),icomp_enuc,rho_Hnuc(n),1)
          call multifab_div_div_c(plotdata(n),icomp_enuc,s(n),rho_comp,1)
         
       end if

       if (use_thermal_diffusion) then
          call multifab_copy_c(plotdata(n),icomp_thermal,thermal(n),1)
       endif


       ! TRACER
       if (plot_trac .and. ntrac .ge. 1) then
          call multifab_copy_c(plotdata(n),icomp_trac,s(n),trac_comp,ntrac)
       end if

    end do

    if (spherical .eq. 1) then

       do n=1,nlevs

          do comp=1,dm
             ! w0mac will contain an edge-centered w0 on a Cartesian grid,
             ! for use in computing divergences.
             call multifab_build(w0mac(n,comp), mla%la(n),1,1,nodal=edge_nodal_flag(comp,:))
             call setval(w0mac(n,comp), ZERO, all=.true.)
          enddo

          ! w0r_cart is w0 but onto a Cartesian grid in cell-centered as
          ! a scalar.  Since w0 is the radial expansion velocity, w0r_cart
          ! is the radial w0 in a zone
          call multifab_build(w0r_cart(n), mla%la(n),1,0)
          call setval(w0r_cart(n), ZERO, all=.true.)

       end do

       if (evolve_base_state) then
          ! put w0 on Cartesian edges as a vector
          call make_w0mac(mla,w0,w0mac,dx,the_bc_tower%bc_tower_array)

          ! put w0 in Cartesian cell-centers as a scalar (the radial expansion velocity)
          call put_1d_array_on_cart(w0,w0r_cart,foextrap_comp,.true.,.false.,dx, &
                                    the_bc_tower%bc_tower_array,mla)
       end if

    end if

    if (plot_base) then

       ! w0
       if (evolve_base_state) then
          call put_1d_array_on_cart(w0,tempfab,1,.true.,.true.,dx, &
                                    the_bc_tower%bc_tower_array,mla)
       else
          do n=1,nlevs
             call setval(tempfab(n), ZERO, all=.true.)
          end do
       end if

       do n=1,nlevs
          call multifab_copy_c(plotdata(n),icomp_w0,tempfab(n),1,dm)
       end do

       ! divw0
       if (spherical .eq. 1) then
          do n=1,nlevs
             call make_divw0(plotdata(n),icomp_divw0,w0(1,:),w0mac(n,:),dx(n,:))
          end do
       else
          do n=1,nlevs
             call make_divw0(plotdata(n),icomp_divw0,w0(n,:),w0mac(n,:),dx(n,:))
          end do
       end if

       ! rho0
       call put_1d_array_on_cart(rho0,tempfab,dm+rho_comp,.false.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla)

       do n=1,nlevs
          call multifab_copy_c(plotdata(n),icomp_rho0,tempfab(n),1,1)
       end do

       ! rhoh0
       if (do_smallscale) then
          h0 = ZERO
       else
          h0 = rhoh0 / rho0
       end if
       call put_1d_array_on_cart(rhoh0,tempfab,dm+rhoh_comp,.false.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla)

       do n=1,nlevs
          call multifab_copy_c(plotdata(n),icomp_rhoh0,tempfab(n),1,1)
       end do

       ! h0
       call put_1d_array_on_cart(h0,tempfab,foextrap_comp,.false.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla)

       do n=1,nlevs
          call multifab_copy_c(plotdata(n),icomp_h0,tempfab(n),1,1)
       end do

       ! p0
       call put_1d_array_on_cart(p0,tempfab,foextrap_comp,.false.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla)
       do n=1,nlevs
          call multifab_copy_c(plotdata(n),icomp_p0,tempfab(n),1,1)
       end do

    end if

    do n = 1,nlevs

       ! RADIAL VELOCITY (spherical only)
       if (spherical .eq. 1) then
          call make_velr(plotdata(n),icomp_velr,u(n),w0r_cart(n),normal(n))
       endif

       ! MAGVEL = |U + w0|
       if (spherical .eq. 1) then
          call make_magvel(plotdata(n),icomp_magvel,icomp_mom,s(n),u(n),w0(1,:),w0mac(n,:))
       else
          call make_magvel(plotdata(n),icomp_magvel,icomp_mom,s(n),u(n),w0(n,:),w0mac(n,:))
       end if

       ! VORTICITY
       call make_vorticity(plotdata(n),icomp_vort,u(n),dx(n,:), &
                           the_bc_tower%bc_tower_array(n))

       ! DIVU
       call multifab_copy_c(plotdata(n),icomp_divu,Source(n),1,1)

       ! ENTHALPY 
       call multifab_copy_c(plotdata(n),icomp_enthalpy,s(n),rhoh_comp)
       call multifab_div_div_c(plotdata(n),icomp_enthalpy,s(n),rho_comp,1)

    end do

    do n=1,nlevs

       ! make_tfromp -> RHOPERT, TFROMP, TPERT, MACHNUMBER, DELTAGAMMA, ENTROPY, AND RHOPERT
       ! make_tfromH -> TFROMP AND DELTA_P
       ! make_conductivity -> CONDUCTIVITY
       if (spherical .eq. 1) then
          
          call make_tfromp(plotdata(n),icomp_tfromp,icomp_tpert,icomp_rhopert, &
                           icomp_rhohpert,icomp_machno,icomp_dg,icomp_entropy,s(n),u(n), &
                           rho0(1,:),rhoh0(1,:),tempbar(1,:),gamma1bar(1,:),p0(1,:),dx(n,:))

          call make_tfromH(plotdata(n),icomp_tfromH,icomp_tpert,icomp_dp,s(n),p0(1,:), &
                           tempbar(1,:),dx(n,:))

       else

          call make_tfromp(plotdata(n),icomp_tfromp,icomp_tpert,icomp_rhopert, &
                           icomp_rhohpert,icomp_machno,icomp_dg,icomp_entropy,s(n),u(n), &
                           rho0(n,:),rhoh0(n,:),tempbar(n,:),gamma1bar(n,:),p0(n,:),dx(n,:))

          call make_tfromH(plotdata(n),icomp_tfromH,icomp_tpert,icomp_dp,s(n),p0(n,:), &
                           tempbar(n,:),dx(n,:))

       end if
       
    end do

    if (use_thermal_diffusion) then
       do n=1,nlevs
          ! this just uses (rho, T, X_k) ---> conductivity
          ! it doesn't need to do anything fancy for spherical
          call make_conductivity(plotdata(n),icomp_conductivity,s(n))
       end do
    end if

    ! the loop over nlevs must count backwards to make sure the finer grids are done first
    do n=nlevs,2,-1
       ! set level n-1 data to be the average of the level n data covering it
       call ml_cc_restriction_c(plotdata(n-1),icomp_tfromp,plotdata(n),icomp_tfromp, &
                                mla%mba%rr(n-1,:),1)
       call ml_cc_restriction_c(plotdata(n-1),icomp_tpert,plotdata(n),icomp_tpert, &
                                mla%mba%rr(n-1,:),1)
       call ml_cc_restriction_c(plotdata(n-1),icomp_rhopert,plotdata(n),icomp_rhopert, &
                                mla%mba%rr(n-1,:),1)
       call ml_cc_restriction_c(plotdata(n-1),icomp_rhohpert,plotdata(n),icomp_rhohpert, &
                                mla%mba%rr(n-1,:),1)
       call ml_cc_restriction_c(plotdata(n-1),icomp_machno,plotdata(n),icomp_machno, &
                                mla%mba%rr(n-1,:),1)
       call ml_cc_restriction_c(plotdata(n-1),icomp_dg,plotdata(n),icomp_dg, &
                                mla%mba%rr(n-1,:),1)
       call ml_cc_restriction_c(plotdata(n-1),icomp_entropy,plotdata(n),icomp_entropy, &
                                mla%mba%rr(n-1,:),1)
       call ml_cc_restriction_c(plotdata(n-1),icomp_tfromH,plotdata(n),icomp_tfromH, &
                                mla%mba%rr(n-1,:),1)
       call ml_cc_restriction_c(plotdata(n-1),icomp_dp,plotdata(n),icomp_dp, &
                                mla%mba%rr(n-1,:),1)
    end do

    do n=1,nlevs

       ! DIFF BETWEEN TFROMP AND TFROMH
       call make_deltaT(plotdata(n),icomp_dT,icomp_tfromp,icomp_tfromH)

       ! PRESSURE GRADIENT
       call multifab_copy_c(plotdata(n),icomp_gp,gpres(n),1,dm)

       ! SPONGE
       call multifab_copy_c(plotdata(n),icomp_sponge,sponge(n),1,1)

    end do

    ! we just made the entropy above.  To compute s - sbar, we need to average
    ! the entropy first, and then compute that.
    ! an average quantity needs ghostcells, so copy entropy into
    ! tempfab
    do n=1,nlevs
       call multifab_copy_c(tempfab(n),1,plotdata(n),icomp_entropy,1)
    end do

    ! fill the ghostcells of tempfab (entropy) so we can properly average it
    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(tempfab(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(tempfab(nlevs),1,foextrap_comp,1, &
                            the_bc_tower%bc_tower_array(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(tempfab(n-1),tempfab(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(tempfab(n),tempfab(n-1), &
                                         tempfab(n)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), &
                                         the_bc_tower%bc_tower_array(n), &
                                         1,foextrap_comp,1,fill_crse_input=.false.)
       enddo

    end if
    
    call average(mla,tempfab,entropybar,dx,1)

    call make_entropypert(plotdata,icomp_entropy,icomp_entropypert,entropybar,dx)

    ! the loop over nlevs must count backwards to make sure the finer grids are done first
    do n=nlevs,2,-1
       ! set level n-1 data to be the average of the level n data covering it
       call ml_cc_restriction_c(plotdata(n-1),icomp_entropypert,plotdata(n), &
                                icomp_entropypert,mla%mba%rr(n-1,:),1)
    end do

    call fabio_ml_multifab_write_d(plotdata, mba%rr(:,1), dirname, plot_names, &
                                   mba%pd(1), prob_lo, prob_hi, time, dx(1,:), &
                                   nOutFiles = nOutFiles, &
                                   lUsingNFiles = lUsingNFiles, prec = prec)

    do n = 1,nlevs
       call destroy(plotdata(n))
       call destroy(tempfab(n))
    end do

    if (spherical .eq. 1) then
       do n=1,nlevs
          call destroy(w0r_cart(n))
          do comp=1,dm
             call destroy(w0mac(n,comp))
          end do
       end do
    end if

    call destroy(bpt)

  end subroutine make_plotfile

end module make_plotfile_module
