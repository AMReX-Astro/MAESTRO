module make_plotfile_module

  use bl_error_module
  use bl_string_module
  use bl_IO_module
  use bl_types
  use fab_module
  use fabio_module
  use boxarray_module
  use ml_boxarray_module
  use multifab_module
  use parallel
  use vort_module
  use geometry
  use plot_variables_module

  use variables

  implicit none

contains

  subroutine make_plotfile(istep,plotdata,u,s,gp,rho_omegadot,mba, &
                           plot_names,time,dx, &
                           the_bc_tower, &
                           s0,p0,temp0,ntrac)

    integer          , intent(in   ) :: istep
    integer          , intent(in   ) :: ntrac
    type(multifab)   , intent(inout) :: plotdata(:)
    type(multifab)   , intent(in   ) :: u(:)
    type(multifab)   , intent(in   ) :: s(:)
    type(multifab)   , intent(in   ) :: gp(:)
    type(multifab)   , intent(in   ) :: rho_omegadot(:)
    type(ml_boxarray), intent(in   ) :: mba
    character(len=20), intent(in   ) :: plot_names(:)
    real(dp_t)       , intent(in   ) :: time,dx(:,:)
    type(bc_tower)   , intent(in   ) :: the_bc_tower
    real(dp_t)       , intent(in   ) :: s0(:,:),p0(:)
    real(dp_t)       , intent(inout) :: temp0(:)

    integer :: n,dm,nlevs,nscal,icomp
    integer :: icomp_tfromrho,icomp_tpert,icomp_rhopert
    integer :: icomp_machno,icomp_deltag
    integer :: icomp_tfromH,icomp_dp
    character(len=7) :: sd_name

    dm = get_dim(mba)
    nlevs = size(plotdata)
    nscal = multifab_ncomp(s(nlevs))

    do n = 1,nlevs

       ! VELOCITY 
       call multifab_copy_c(plotdata(n),1            ,    u(n),1,dm)

       ! DENSITY AND (RHO H) 
       icomp = dm+rho_comp
       call multifab_copy_c(plotdata(n),icomp,s(n),1,2)

       ! SPECIES
       icomp = dm+spec_comp
       call make_XfromrhoX(plotdata(n),icomp,s(n))

       ! TRACER
       if (ntrac .ge. 1) then
         icomp = dm+trac_comp
         call multifab_copy_c(plotdata(n),icomp,s(n),trac_comp,ntrac)
       end if

       ! VORTICITY
       icomp = derive_comp
       call make_vorticity (plotdata(n),icomp,u(n),dx(n,:), &
                            the_bc_tower%bc_tower_array(n))

       ! ENTHALPY (RHO H)
       icomp = derive_comp+1
       call make_enthalpy  (plotdata(n),icomp,s(n))

    end do

    if (spherical .eq. 1) then

      do n = 1,nlevs

       ! RHOPERT & TEMP (FROM RHO) & TPERT & MACHNO & (GAM1 - GAM10)
       icomp_rhopert  = derive_comp+2
       icomp_tfromrho = derive_comp+3
       icomp_tpert    = derive_comp+5
       icomp_machno   = derive_comp+6
       icomp_deltag   = derive_comp+8
       call make_tfromrho  (plotdata(n),icomp_tfromrho,icomp_tpert,icomp_rhopert, &
                            icomp_machno,icomp_deltag, &
                            s(n),u(n),s0,temp0,p0,time,dx(n,:))

       ! TEMP (FROM H) & DELTA_P
       icomp_tfromH  = derive_comp+4
       icomp_dp      = derive_comp+7
       call make_tfromH    (plotdata(n),icomp_tfromH,icomp_dp,s(n),p0,temp0,dx(n,:))

      end do

    else

      do n = 1,nlevs

       ! RHOPERT & TEMP (FROM RHO) & TPERT & MACHNO & (GAM1 - GAM10)
       icomp_rhopert  = derive_comp+2
       icomp_tfromrho = derive_comp+3
       icomp_tpert    = derive_comp+5
       icomp_machno   = derive_comp+6
       icomp_deltag   = derive_comp+8
       call make_tfromrho  (plotdata(n),icomp_tfromrho,icomp_tpert,icomp_rhopert, &
                            icomp_machno,icomp_deltag, &
                            s(n),u(n),s0,temp0,p0,time,dx(n,:))

       ! TEMP (FROM H) & DELTA_P
       icomp_tfromH  = derive_comp+4
       icomp_dp      = derive_comp+7
       call make_tfromH    (plotdata(n),icomp_tfromH,icomp_dp,s(n),p0,temp0,dx(n,:))

      end do

    end if

    do n = 1,nlevs

      ! PRESSURE GRADIENT
      icomp = derive_comp+9
      call multifab_copy_c(plotdata(n),icomp,gp(n),1,dm)

    end do


    do n = 1,nlevs

       ! OMEGADOT
       icomp = derive_spec_comp
       call make_omegadot(plotdata(n),icomp,s(n),rho_omegadot(n))

    enddo


    write(unit=sd_name,fmt='("plt",i4.4)') istep
    call fabio_ml_multifab_write_d(plotdata, mba%rr(:,1), sd_name, plot_names, &
                                    mba%pd(1), time, dx(1,:))

  end subroutine make_plotfile

end module make_plotfile_module
