module checkpoint_module

  use bl_types
  use multifab_module

  implicit none

  private

  public :: checkpoint_write, checkpoint_read

contains

  subroutine checkpoint_write(dirname, mfs, mfs_nodal, dSdt, Source_old, Source_new, &
                              rho_omegadot2, rho_Hnuc2, thermal2, rrs, time, dt)

    use parallel
    use bl_IO_module
    use fabio_module
    use bl_prof_module
    use probin_module, only: verbose, nOutFiles, lUsingNFiles, use_thermal_diffusion
    use variables, only: rel_eps

    type(multifab)  , intent(in) :: mfs(:), mfs_nodal(:)
    type(multifab)  , intent(in) :: dSdt(:), Source_old(:), Source_new(:)
    type(multifab)  , intent(in) :: rho_omegadot2(:), rho_Hnuc2(:), thermal2(:)
    integer         , intent(in) :: rrs(:,:)
    character(len=*), intent(in) :: dirname
    real(kind=dp_t) , intent(in) :: time, dt

    ! local
    integer :: un, nlevs
    character(len=256) :: header, sd_name, sd_name_nodal

    namelist /chkpoint/ nlevs

    type(bl_prof_timer), save :: bpt

    call build(bpt, "checkpoint_write")

    if ( parallel_IOProcessor() ) then
       call fabio_mkdir(dirname)
    end if

    call parallel_barrier()

    write(unit=sd_name, fmt='(a,"/State")') trim(dirname)
    call fabio_ml_multifab_write_d(mfs, rrs(:,1), sd_name, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
      write(6,*) 'Writing state to checkpoint file ',trim(sd_name)
    end if

    write(unit=sd_name, fmt='(a,"/dSdt")') trim(dirname)
    call fabio_ml_multifab_write_d(dSdt, rrs(:,1), sd_name, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
      write(6,*) 'Writing state to checkpoint file ',trim(sd_name)
    end if

    write(unit=sd_name, fmt='(a,"/Source_old")') trim(dirname)
    call fabio_ml_multifab_write_d(Source_old, rrs(:,1), sd_name, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
      write(6,*) 'Writing state to checkpoint file ',trim(sd_name)
    end if

    write(unit=sd_name, fmt='(a,"/Source_new")') trim(dirname)
    call fabio_ml_multifab_write_d(Source_new, rrs(:,1), sd_name, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
      write(6,*) 'Writing state to checkpoint file ',trim(sd_name)
    end if

    write(unit=sd_name, fmt='(a,"/rho_omegadot2")') trim(dirname)
    call fabio_ml_multifab_write_d(rho_omegadot2, rrs(:,1), sd_name, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
      write(6,*) 'Writing state to checkpoint file ',trim(sd_name)
    end if

    write(unit=sd_name, fmt='(a,"/rho_Hnuc2")') trim(dirname)
    call fabio_ml_multifab_write_d(rho_Hnuc2, rrs(:,1), sd_name, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)

    if (use_thermal_diffusion) then
       write(unit=sd_name, fmt='(a,"/thermal2")') trim(dirname)
       call fabio_ml_multifab_write_d(thermal2, rrs(:,1), sd_name, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)
    end if

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
      write(6,*) 'Writing state to checkpoint file ',trim(sd_name)
    end if

    write(unit=sd_name_nodal, fmt='(a,"/Pressure")') trim(dirname)
    call fabio_ml_multifab_write_d(mfs_nodal, rrs(:,1), sd_name_nodal, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
      write(6,*) 'Writing state to checkpoint file ',trim(sd_name_nodal)
      write(6,*)
    end if

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
       write(unit=un,fmt=1000) dt
       write(unit=un,fmt=1000) time
       write(unit=un,fmt=1000) rel_eps
       close(un)
    end if

    call destroy(bpt)

1000 format(32(e30.20,1x))

  end subroutine checkpoint_write

  subroutine checkpoint_read(mfs, mfs_nodal, dSdt, Source_old, Source_new, rho_omegadot2, &
                             rho_Hnuc2, thermal2, dirname, time_out, dt_out, nlevs_out)

    use parallel
    use bl_IO_module
    use fabio_module
    use bl_prof_module
    use variables, only: rel_eps
    use probin_module, only: use_thermal_diffusion

    type(multifab  ),                pointer :: mfs(:), mfs_nodal(:)
    type(multifab  ),                pointer :: dSdt(:), Source_old(:), Source_new(:)
    type(multifab  ),                pointer :: rho_omegadot2(:), rho_Hnuc2(:), thermal2(:)
    character(len=*), intent(in   )          :: dirname
    integer         , intent(  out)          :: nlevs_out
    real(kind=dp_t) , intent(  out)          :: time_out, dt_out

    ! local
    integer            :: un
    character(len=256) :: header, sd_name
    integer            :: nlevs
    real(kind=dp_t)    :: time, dt

    namelist /chkpoint/ nlevs

    type(bl_prof_timer), save :: bpt

    call build(bpt, "checkpoint_read")

!   First read the header information
    header = "Header"
    un = unit_new()
    open(unit=un, &
         file = trim(dirname) // "/" // trim(header), &
         status = "old", &
         action = "read")
    read(unit=un, nml = chkpoint)

    read(unit=un,fmt=*) dt
    read(unit=un,fmt=*) time
    read(unit=un,fmt=*) rel_eps
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

!   Read the Source_new data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/Source_new")') trim(dirname)
    call fabio_ml_multifab_read_d(Source_new, sd_name)

!   Read the rho_omegadot2 data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/rho_omegadot2")') trim(dirname)
    call fabio_ml_multifab_read_d(rho_omegadot2, sd_name)

!   Read the rho_Hnuc2 data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/rho_Hnuc2")') trim(dirname)
    call fabio_ml_multifab_read_d(rho_Hnuc2, sd_name)

!   Read the thermal2 data into a multilevel multifab.
    if (use_thermal_diffusion) then
       write(unit=sd_name, fmt='(a,"/thermal2")') trim(dirname)
       call fabio_ml_multifab_read_d(thermal2, sd_name)
    end if

    call destroy(bpt)

  end subroutine checkpoint_read

end module checkpoint_module
