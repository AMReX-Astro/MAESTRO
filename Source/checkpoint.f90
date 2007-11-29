module checkpoint_module

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

  implicit none

  private
  public :: checkpoint_write, checkpoint_read

contains

  subroutine checkpoint_write(dirname, mfs, mfs_nodal, dSdt, Source_old, &
                              rho_omegadot2, rrs, dx, time_in, dt_in, verbose)
    type(multifab), intent(in) :: mfs(:), mfs_nodal(:)
    type(multifab), intent(in) :: dSdt(:), Source_old(:)
    type(multifab), intent(in) :: rho_omegadot2(:)
    integer        , intent(in) :: rrs(:,:)
    real(kind=dp_t), intent(in) :: dx(:,:)
    character(len=*), intent(in) :: dirname
    real(kind=dp_t), intent(in) :: time_in, dt_in
    integer        , intent(in) :: verbose
    integer :: n, i
    character(len=128) :: header, sd_name, sd_name_nodal
    integer :: un

    integer         :: nlevs, dm
    real(kind=dp_t) :: time, dt

    namelist /chkpoint/ time
    namelist /chkpoint/ dt
    namelist /chkpoint/ nlevs
    namelist /chkpoint/ dm

    if ( parallel_IOProcessor() ) then
       call fabio_mkdir(dirname)
    end if

    call parallel_barrier()

    write(unit=sd_name, fmt='(a,"/State")') trim(dirname)
    call fabio_ml_multifab_write_d(mfs, rrs(:,1), sd_name)

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
      write(6,*) 'Writing    state to checkpoint file ',trim(sd_name)
    end if

    write(unit=sd_name, fmt='(a,"/dSdt")') trim(dirname)
    call fabio_ml_multifab_write_d(dSdt, rrs(:,1), sd_name)

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
      write(6,*) 'Writing    state to checkpoint file ',trim(sd_name)
    end if

    write(unit=sd_name, fmt='(a,"/Source_old")') trim(dirname)
    call fabio_ml_multifab_write_d(Source_old, rrs(:,1), sd_name)

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
      write(6,*) 'Writing    state to checkpoint file ',trim(sd_name)
    end if

    write(unit=sd_name, fmt='(a,"/rho_omegadot2")') trim(dirname)
    call fabio_ml_multifab_write_d(rho_omegadot2, rrs(:,1), sd_name)

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
      write(6,*) 'Writing    state to checkpoint file ',trim(sd_name)
    end if

    write(unit=sd_name_nodal, fmt='(a,"/Pressure")') trim(dirname)
    call fabio_ml_multifab_write_d(mfs_nodal, rrs(:,1), sd_name_nodal)

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
      write(6,*) 'Writing state to checkpoint file ',trim(sd_name_nodal)
    end if

    time = time_in
      dt =   dt_in

      dm = size(dx,dim=2)

    ! Note: parallel fails on Bassi if this is done on all processors
    if (parallel_IOProcessor()) then
       header = "Header"
       un = unit_new()
       open(unit=un, &
            file = trim(dirname) // "/" // trim(header), &
            form = "formatted", access = "sequential", &
            status = "replace", action = "write")
       nlevs = size(mfs)
       write(unit=un, nml = chkpoint)
       do n = 1,nlevs
          write(unit=un,fmt=*) (dx(n,i), i=1,dm)
       end do
       do n = 1,nlevs-1
          write(unit=un,fmt=*) rrs(n,1)
       end do
       close(un)
    end if

  end subroutine checkpoint_write

  subroutine checkpoint_read(mfs, mfs_nodal, dSdt, Source_old, rho_omegadot2, &
       dirname, time_out, dt_out, nlevs_out)
    use bl_IO_module
    type(multifab  ),                pointer :: mfs(:), mfs_nodal(:)
    type(multifab  ),                pointer :: dSdt(:), Source_old(:), rho_omegadot2(:)
    character(len=*), intent(in   )          :: dirname
    integer         , intent(  out)          :: nlevs_out
    real(kind=dp_t) , intent(  out)          :: time_out, dt_out

    integer         ,                pointer :: rrs(:)
    real(kind=dp_t) ,                pointer :: dx(:,:)

    integer :: n, i
    character(len=128) :: header, sd_name
    integer :: un

    integer         :: nlevs, dm
    real(kind=dp_t) :: time, dt

    namelist /chkpoint/ nlevs
    namelist /chkpoint/ time
    namelist /chkpoint/ dt
    namelist /chkpoint/ dm

!   First read the header information
    header = "Header"
    un = unit_new()
    open(unit=un, &
         file = trim(dirname) // "/" // trim(header), &
         status = "old", &
         action = "read")
    read(unit=un, nml = chkpoint)

    allocate( dx(nlevs,dm))
    allocate(rrs(nlevs-1))

    do n = 1,nlevs
       read(unit=un,fmt=*) (dx(n,i), i=1,dm)
    end do
    do n = 1,nlevs-1
       read(unit=un,fmt=*) rrs(n)
    end do
    close(un)

     time_out = time
       dt_out = dt
    nlevs_out = nlevs

!   Read the state data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/State")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs, sd_name)

!   Read the pressure data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/Pressure")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_nodal, sd_name)

!   Read the dSdt data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/dSdt")') trim(dirname)
    call fabio_ml_multifab_read_d(dSdt, sd_name)

!   Read the Source_old data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/Source_old")') trim(dirname)
    call fabio_ml_multifab_read_d(Source_old, sd_name)

!   Read the rho_omegadot2 data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/rho_omegadot2")') trim(dirname)
    call fabio_ml_multifab_read_d(rho_omegadot2, sd_name)

    deallocate(dx,rrs)

  end subroutine checkpoint_read

end module checkpoint_module
