module make_plotfile_module

  use bl_types
  use multifab_module
  use define_bc_module
  use ml_boxarray_module

  implicit none

  private
  public :: get_plot_names, make_plotfile

contains

  subroutine get_plot_names(dm,plot_names)

    use plot_variables_module
    use variables
    use network, only: nspec, short_spec_names
    use probin_module, only: plot_spec, plot_trac, plot_base
    use geometry, only: spherical

    integer          , intent(in   ) :: dm
    character(len=20), intent(inout) :: plot_names(:)

    ! Local variables
    integer :: comp

    plot_names(icomp_vel  ) = "x_vel"
    plot_names(icomp_vel+1) = "y_vel"
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
       plot_names(icomp_w0+1) = "w0_y"
       if (dm > 2) plot_names(icomp_w0+2) = "w0_z"
       plot_names(icomp_divw0) = "divw0"
       plot_names(icomp_rho0)  = "rho0"
       plot_names(icomp_rhoh0) = "rhoh0"
       plot_names(icomp_p0)    = "p0"
    end if

    if (spherical .eq. 1) then
       plot_names(icomp_velr) = "radial_velocity"
    endif

    plot_names(icomp_magvel)      = "magvel"
    plot_names(icomp_velplusw0)   = "velplusw0"
    plot_names(icomp_mom)         = "momentum"
    plot_names(icomp_vort)        = "vort"
    plot_names(icomp_divu)        = "divu"
    plot_names(icomp_enthalpy)    = "enthalpy"
    plot_names(icomp_rhopert)     = "rhopert"
    plot_names(icomp_rhohpert)     = "rhohpert"
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
    plot_names(icomp_gp+1)        = "gpy"
    if (dm > 2) plot_names(icomp_gp+2) = "gpz"

    if (plot_spec) then
       do comp = 1, nspec
          plot_names(icomp_omegadot+comp-1) = &
               "omegadot(" // trim(short_spec_names(comp)) // ")"
       end do
       plot_names(icomp_enuc) = "enucdot"
    end if

  end subroutine get_plot_names

  subroutine make_plotfile(dirname,mla,u,s,gpres,rho_omegadot,Source,sponge,&
                           mba,plot_names,time,dx,the_bc_tower,w0,rho0,rhoh0,p0,tempbar, &
                           gamma1bar,div_coeff,normal)

    use bl_prof_module
    use fabio_module
    use vort_module
    use variables
    use plot_variables_module
    use fill_3d_module
    use probin_module, only: nOutFiles, lUsingNFiles, plot_spec, plot_trac, plot_base
    use probin_module, only: single_prec_plotfiles
    use geometry, only: spherical, nr_fine
    use average_module
    use ml_restriction_module
    use multifab_physbc_module
    use multifab_fill_ghost_module
    use bl_constants_module

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
    real(dp_t)       , intent(in   ) :: w0(:,0:)
    real(dp_t)       , intent(in   ) :: rho0(:,0:)
    real(dp_t)       , intent(in   ) :: rhoh0(:,0:)
    real(dp_t)       , intent(in   ) :: p0(:,0:)
    real(dp_t)       , intent(in   ) :: tempbar(:,0:)
    real(dp_t)       , intent(in   ) :: gamma1bar(:,0:)
    real(dp_t)       , intent(in   ) :: div_coeff(:,0:)
    type(multifab)   , intent(in   ) :: normal(:)
    
    type(multifab) :: plotdata(mla%nlevel)
    type(multifab) :: tempfab(mla%nlevel)
    type(multifab) :: w0mac(mla%nlevel,mla%dim)
    real(dp_t), allocatable :: entropybar(:,:)

    integer :: n,dm,nlevs,prec,comp

    logical :: umac_nodal_flag(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_plotfile")

    dm = get_dim(mba)
    nlevs = size(u)

    if (single_prec_plotfiles) then
       prec = FABIO_SINGLE
    else
       prec = FABIO_DOUBLE
    endif

    do n = 1,nlevs

       call multifab_build(plotdata(n), mla%la(n), n_plot_comps, 0)
       call multifab_build(tempfab(n), mla%la(n), dm, 1)
              
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

    end do

    if (plot_base) then

       ! w0
       call put_1d_array_on_cart(nlevs,w0,tempfab,1,.true.,.true.,dx, &
                                 the_bc_tower%bc_tower_array,mla,normal=normal)

       do n=1,nlevs
          call multifab_copy_c(plotdata(n),icomp_w0,tempfab(n),1,dm)
       end do

       if (spherical .eq. 1) then
          do n=1,nlevs
             do comp=1,dm
                umac_nodal_flag = .false.
                umac_nodal_flag(comp) = .true.
                call multifab_build(w0mac(n,comp), mla%la(n),  1, 1, nodal = umac_nodal_flag)
                call setval(w0mac(n,comp), ZERO, all=.true.)
             end do
          end do
          call put_w0_on_edges(mla,w0,w0mac,dx,div_coeff,the_bc_tower)
       end if

       ! divw0
       do n=1,nlevs
          call make_divw0(tempfab(n),w0(n,:),w0mac(n,:),dx(n,:))
       end do
         
       do n=1,nlevs
          call multifab_copy_c(plotdata(n),icomp_divw0,tempfab(n),1,1)
       end do

       ! rho0
       call put_1d_array_on_cart(nlevs,rho0,tempfab,dm+rho_comp,.false.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla,normal=normal)

       do n=1,nlevs
          call multifab_copy_c(plotdata(n),icomp_rho0,tempfab(n),1,1)
       end do

       ! rhoh0
       call put_1d_array_on_cart(nlevs,rhoh0,tempfab,dm+rhoh_comp,.false.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla,normal=normal)

       do n=1,nlevs
          call multifab_copy_c(plotdata(n),icomp_rhoh0,tempfab(n),1,1)
       end do

       ! p0
       call put_1d_array_on_cart(nlevs,p0,tempfab,foextrap_comp,.false.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla,normal=normal)
       do n=1,nlevs
          call multifab_copy_c(plotdata(n),icomp_p0,tempfab(n),1,1)
       end do

    end if

    do n = 1,nlevs

       ! MAGVEL & MOMENTUM
       call make_magvel (plotdata(n),icomp_magvel,icomp_mom,u(n),s(n))

       ! RADIAL VELOCITY (spherical only)
       if (spherical .eq. 1) then
          call make_velr (n,plotdata(n),icomp_velr,u(n),w0(n,:),normal(n),dx(n,:))
       endif

       ! VEL_PLUS_W0
       call make_velplusw0 (n,plotdata(n),icomp_velplusw0,u(n),w0(n,:),normal(n),dx(n,:))

       ! VORTICITY
       call make_vorticity(plotdata(n),icomp_vort,u(n),dx(n,:), &
                           the_bc_tower%bc_tower_array(n))

       ! DIVU
       call multifab_copy_c(plotdata(n),icomp_divu,Source(n),1,1)

       ! ENTHALPY 
       call make_enthalpy(plotdata(n),icomp_enthalpy,s(n))

       ! RHOPERT & TEMP (FROM RHO) & TPERT & MACHNO & (GAM1 - GAM10) & Entropy & RHOHPERT
       call make_tfromp(n,plotdata(n), &
                        icomp_tfromp,icomp_tpert,icomp_rhopert,icomp_rhohpert, &
                        icomp_machno,icomp_dg,icomp_entropy, &
                        s(n),u(n),rho0(n,:),rhoh0(n,:),tempbar(n,:),gamma1bar(n,:), &
                        p0(n,:),dx(n,:))

       ! TEMP (FROM H) & DELTA_P
       call make_tfromH(n,plotdata(n),icomp_tfromH,icomp_dp,s(n),p0(n,:), &
                        tempbar(n,:),dx(n,:))
       
       ! DIFF BETWEEN TFROMP AND TFROMH
       call make_deltaT (plotdata(n),icomp_dT,icomp_tfromp,icomp_tfromH)

       ! PRESSURE GRADIENT
       call multifab_copy_c(plotdata(n),icomp_gp,gpres(n),1,dm)

       ! SPONGE
       call multifab_copy_c(plotdata(n),icomp_sponge,sponge(n),1,1)

    end do

    ! we just made the entropy above.  To compute s - sbar, we need to average
    ! the entropy first, and then compute that.
    allocate(entropybar(nlevs,0:nr_fine-1))

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
       call multifab_physbc(tempfab(nlevs),1,foextrap_comp,1,the_bc_tower%bc_tower_array(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(tempfab(n-1)    ,tempfab(n)    ,mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(tempfab(n),tempfab(n-1), &
                                         tempfab(n)%ng,mla%mba%rr(n-1,:), &
                                         the_bc_tower%bc_tower_array(n-1), the_bc_tower%bc_tower_array(n), &
                                         1,foextrap_comp,1)
       enddo

    end if
    

    call average(mla,tempfab,entropybar,dx,1)

    do n = 1,nlevs
       call make_entropypert(n,plotdata(n),icomp_entropy,icomp_entropypert,entropybar(n,0:),dx(n,:))
    enddo


    if (plot_spec) then
       ! OMEGADOT
       do n = 1,nlevs
          call make_omegadot(plotdata(n),icomp_omegadot,icomp_enuc,s(n),rho_omegadot(n))
       end do
    end if

    call fabio_ml_multifab_write_d(plotdata, mba%rr(:,1), dirname, plot_names, &
                                   mba%pd(1), time, dx(1,:), nOutFiles = nOutFiles, &
                                   lUsingNFiles = lUsingNFiles, prec = prec)
    do n = 1,nlevs
       call destroy(plotdata(n))
       call destroy(tempfab(n))
    end do

    if (spherical .eq. 1 .and. plot_base) then
       do n=1,nlevs
          do comp=1,dm
             call destroy(w0mac(n,comp))
          end do
       end do
    end if

    call destroy(bpt)

  end subroutine make_plotfile

end module make_plotfile_module
