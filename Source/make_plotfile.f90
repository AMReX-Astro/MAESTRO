module make_plotfile_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_boxarray_module

  implicit none

  private
  public :: get_plot_names, make_plotfile

contains

  subroutine get_plot_names(dm,plot_names,plot_spec,plot_trac)

    use plot_variables_module
    use probin_module, ONLY: use_big_h
    use variables
    use network, only: nspec, short_spec_names

    integer          , intent(in   ) :: dm
    logical          , intent(in   ) :: plot_spec,plot_trac
    character(len=20), intent(inout) :: plot_names(:)

    ! Local variables
    integer :: comp

    plot_names(icomp_vel  ) = "x_vel"
    plot_names(icomp_vel+1) = "y_vel"
    if (dm > 2) then
       plot_names(icomp_vel+2) = "z_vel"
    end if
    plot_names(icomp_rho)  = "density"
    if(use_big_h) then
       plot_names(icomp_rhoh) = "rhoH"
    else
       plot_names(icomp_rhoh) = "rhoh"
    end if

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

    plot_names(icomp_magvel)   = "magvel"
    plot_names(icomp_mom)      = "momentum"
    plot_names(icomp_vort)     = "vort"
    plot_names(icomp_divu)     = "divu"
    plot_names(icomp_enthalpy) = "enthalpy"
    plot_names(icomp_rhopert)  = "rhopert"
    plot_names(icomp_tfromrho) = "tfromrho"
    if(use_big_h) then
       plot_names(icomp_tfromH)   = "tfromH"
    else
       plot_names(icomp_tfromH)   = "tfromh"
    end if
    plot_names(icomp_tpert)    = "tpert"
    plot_names(icomp_machno)   = "Machnumber"
    plot_names(icomp_dp)       = "deltap"
    plot_names(icomp_dg)       = "deltagamma"
    plot_names(icomp_spert)    = "spert"
    plot_names(icomp_dT)       = "deltaT"
    plot_names(icomp_sponge)   = "sponge"
    plot_names(icomp_gp)       = "gpx"
    plot_names(icomp_gp+1)     = "gpy"
    if (dm > 2) plot_names(icomp_gp+2) = "gpz"

    if (plot_spec) then
       do comp = 1, nspec
          plot_names(icomp_omegadot+comp-1) = &
               "omegadot(" // trim(short_spec_names(comp)) // ")"
       end do
       plot_names(icomp_enuc) = "enucdot"
    end if

  end subroutine get_plot_names

  subroutine make_plotfile(dirname,mla,u,s,gpres,rho_omegadot,Source,sponge, &
                           mba,plot_names,time,dx,the_bc_tower,s0,p0,tempbar, &
                           plot_spec,plot_trac)

    use bl_prof_module
    use fabio_module
    use vort_module
    use variables
    use plot_variables_module
    use probin_module, ONLY: use_big_h

    character(len=*) , intent(in   ) :: dirname
    type(ml_layout)  , intent(in   ) :: mla
    type(multifab)   , intent(inout) :: u(:)
    type(multifab)   , intent(inout) :: s(:)
    type(multifab)   , intent(in   ) :: gpres(:)
    type(multifab)   , intent(in   ) :: rho_omegadot(:)
    type(multifab)   , intent(in   ) :: Source(:)
    type(multifab)   , intent(in   ) :: sponge(:)
    type(ml_boxarray), intent(in   ) :: mba
    character(len=20), intent(in   ) :: plot_names(:)
    real(dp_t)       , intent(in   ) :: time,dx(:,:)
    type(bc_tower)   , intent(in   ) :: the_bc_tower
    real(dp_t)       , intent(in   ) :: s0(:,0:,:)
    real(dp_t)       , intent(in   ) :: p0(:,0:)
    real(dp_t)       , intent(in   ) :: tempbar(:,0:,:)
    logical          , intent(in   ) :: plot_spec,plot_trac

    type(multifab) :: plotdata(mla%nlevel)

    integer :: n,dm,nlevs

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_plotfile")

    dm = get_dim(mba)
    nlevs = size(u)

    do n = 1,nlevs

       call multifab_build(plotdata(n), mla%la(n), n_plot_comps, 0)

       ! VELOCITY 
       call multifab_copy_c(plotdata(n),icomp_vel,u(n),1,dm)

       ! DENSITY AND (RHO H) 
       call multifab_copy_c(plotdata(n),icomp_rho,s(n),rho_comp,2)

       ! SPECIES
       if (plot_spec) then
          call make_XfromrhoX(plotdata(n),icomp_spec,s(n))
       end if

       ! TRACER
       if (plot_trac .and. ntrac .ge. 1) then
          call multifab_copy_c(plotdata(n),icomp_trac,s(n),trac_comp,ntrac)
       end if

       ! MAGVEL & MOMENTUM
       call make_magvel (plotdata(n),icomp_magvel,icomp_mom,u(n),s(n))

       ! VORTICITY
       call make_vorticity(plotdata(n),icomp_vort,u(n),dx(n,:), &
                           the_bc_tower%bc_tower_array(n))

       ! DIVU
       call multifab_copy_c(plotdata(n),icomp_divu,Source(n),1,1)

       ! ENTHALPY 
       call make_enthalpy(plotdata(n),icomp_enthalpy,s(n))

       ! RHOPERT & TEMP (FROM RHO) & TPERT & MACHNO & (GAM1 - GAM10)
       call make_tfromrho(n,plotdata(n),icomp_tfromrho,icomp_tpert,icomp_rhopert, &
                          icomp_machno,icomp_dg,icomp_spert, &
                          s(n),u(n),s0(n,:,:),tempbar(n,:,1),p0(n,:),dx(n,:))

       ! TEMP (FROM H) & DELTA_P
       call make_tfromH(n,plotdata(n),icomp_tfromH,icomp_dp,s(n),p0(n,:), &
                        tempbar(n,:,1),dx(n,:))
       
       ! DIFF BETWEEN TFROMRHO AND TFROMH
       call make_deltaT (plotdata(n),icomp_dT,icomp_tfromrho,icomp_tfromH)

       ! PRESSURE GRADIENT
       call multifab_copy_c(plotdata(n),icomp_gp,gpres(n),1,dm)

       ! SPONGE
       call multifab_copy_c(plotdata(n),icomp_sponge,sponge(n),1,1)

    end do

    if (plot_spec) then
       ! OMEGADOT
       do n = 1,nlevs
          call make_omegadot(plotdata(n),icomp_omegadot,icomp_enuc,s(n),rho_omegadot(n))
       end do
    end if

    call fabio_ml_multifab_write_d(plotdata, mba%rr(:,1), dirname, plot_names, &
                                   mba%pd(1), time, dx(1,:))
    do n = 1,nlevs
       call destroy(plotdata(n))
    end do

    call destroy(bpt)

  end subroutine make_plotfile

end module make_plotfile_module
