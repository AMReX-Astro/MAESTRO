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

  use variables

  implicit none

contains

  subroutine make_plotfile(istep,plotdata,uold,sold,gp,mba,plot_names,time,dx,the_bc_tower, &
                           rho0,p0,temp0)

    integer          , intent(in   ) :: istep
    type(multifab)   , intent(inout) :: plotdata(:)
    type(multifab)   , intent(in   ) :: uold(:)
    type(multifab)   , intent(in   ) :: sold(:)
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
    nscal = multifab_ncomp(sold(nlevs))

    do n = 1,nlevs
       call multifab_copy_c(plotdata(n),1            ,    uold(n),1,dm)
       call multifab_copy_c(plotdata(n),rho_comp+dm  ,    sold(n),1,nscal)

       icomp = derive_comp
       call make_vorticity(plotdata(n),icomp,uold(n),dx(n,:), &
                           the_bc_tower%bc_tower_array(n))

       icomp = derive_comp+1
       call make_rhopert (plotdata(n) ,icomp,sold(n),rho0)

       icomp = derive_comp+2
       call make_enthalpy(plotdata(n),icomp,sold(n))

       icomp = derive_comp+3
       call make_tfromrho(plotdata(n) ,icomp,sold(n),temp0,p0,time,dx(n,:))

       icomp = derive_comp+4
       call make_tfromH  (plotdata(n)   ,icomp,sold(n),p0,temp0)

       icomp = derive_comp+5
       call make_tpert   (plotdata(n)   ,icomp,sold(n),p0,temp0)

       icomp = derive_comp+6
       call make_machno  (plotdata(n)  ,icomp,uold(n),sold(n),p0,temp0)

       icomp = derive_comp+7
       call make_deltap  (plotdata(n)  ,icomp,sold(n),p0,temp0)

       icomp = derive_comp+8
       call multifab_copy_c(plotdata(n),icomp,gp(n),1,dm)

     end do

     write(unit=sd_name,fmt='("plt",i4.4)') istep
     call fabio_ml_multifab_write_d(plotdata, mba%rr(:,1), sd_name, plot_names, &
                                    mba%pd(1), time, dx(1,:))

  end subroutine make_plotfile

end module make_plotfile_module
