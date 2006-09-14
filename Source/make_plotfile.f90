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
  use rhopert_module
  use tfromrho_module
  use tfromh_module
  use tpert_module
  use enthalpy_module
  use machno_module
  use deltap_module
  use XfromrhoX_module

  use variables

  implicit none

contains

  subroutine make_plotfile(istep,plotdata,u,s,gp,mba,plot_names,time,dx,the_bc_tower, &
                           rho0,p0,temp0)

    integer          , intent(in   ) :: istep
    type(multifab)   , intent(inout) :: plotdata(:)
    type(multifab)   , intent(in   ) :: u(:)
    type(multifab)   , intent(in   ) :: s(:)
    type(multifab)   , intent(in   ) :: gp(:)
    type(ml_boxarray), intent(in   ) :: mba
    character(len=20), intent(in   ) :: plot_names(:)
    real(dp_t)       , intent(in   ) :: time,dx(:,:)
    type(bc_tower)   , intent(in   ) :: the_bc_tower
    real(dp_t)       , intent(in   ) :: rho0(:),p0(:),temp0(:)

    integer :: n,dm,nlevs,nscal,icomp
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

       ! VORTICITY
       icomp = derive_comp
       call make_vorticity (plotdata(n),icomp,u(n),dx(n,:), &
                            the_bc_tower%bc_tower_array(n))

       ! DENSITY PERTURBATION
       icomp = derive_comp+1
       call make_rhopert   (plotdata(n),icomp,s(n),rho0)

       ! ENTHALPY (RHO H)
       icomp = derive_comp+2
       call make_enthalpy  (plotdata(n),icomp,s(n))

       ! TEMP (FROM RHO)
       icomp = derive_comp+3
       call make_tfromrho  (plotdata(n),icomp,s(n),temp0,p0,time,dx(n,:))

       ! TEMP (FROM H)
       icomp = derive_comp+4
       call make_tfromH    (plotdata(n),icomp,s(n),p0,temp0)

       ! TEMPERATURE PERTURBATION
       icomp = derive_comp+5
       call make_tpert     (plotdata(n),icomp,s(n),p0,temp0)

       ! MACH NUMBER
       icomp = derive_comp+6
       call make_machno    (plotdata(n),icomp,u(n),s(n),p0,temp0)

       ! DELTA P (P - P0)
       icomp = derive_comp+7
       call make_deltap    (plotdata(n),icomp,s(n),p0,temp0)

       ! PRESSURE GRADIENT
       icomp = derive_comp+8
       call multifab_copy_c(plotdata(n),icomp,gp(n),1,dm)

     end do

     write(unit=sd_name,fmt='("plt",i4.4)') istep
     call fabio_ml_multifab_write_d(plotdata, mba%rr(:,1), sd_name, plot_names, &
                                    mba%pd(1), time, dx(1,:))

  end subroutine make_plotfile

end module make_plotfile_module
