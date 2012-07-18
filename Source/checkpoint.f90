! checkpoint_write and checkpoint_read write/read the full state of
! the simulation to/from a restart file.  A separate set of routines
! (in base_io_module) are used to store the base state).

module checkpoint_module

  use bl_types, only: dp_t
  use multifab_module

  implicit none

  private

  public :: checkpoint_write, checkpoint_read

contains

  subroutine checkpoint_write(dirname, mfs, mfs_nodal, dSdt, Source_old, Source_new, &
                              rho_omegadot2, rho_Hnuc2, rho_Hext, thermal2, &
                              rrs, dt)

    use parallel, only: parallel_IOProcessor, parallel_barrier
    use bl_IO_module, only: unit_new
    use fabio_module, only: fabio_mkdir, fabio_ml_multifab_write_d
    use bl_prof_module, only: bl_prof_timer, build, destroy
    use probin_module, only: verbose, nOutFiles, lUsingNFiles, &
                             use_thermal_diffusion, plot_Hext
    use variables, only: rel_eps
    use time_module, only: time
    use cputime_module, only: get_cputime

    type(multifab)  , intent(in) :: mfs(:), mfs_nodal(:)
    type(multifab)  , intent(in) :: dSdt(:), Source_old(:), Source_new(:)
    type(multifab)  , intent(in) :: rho_omegadot2(:), rho_Hnuc2(:), rho_Hext(:), thermal2(:)
    integer         , intent(in) :: rrs(:,:)
    character(len=*), intent(in) :: dirname
    real(kind=dp_t) , intent(in) :: dt

    ! local
    integer :: un, nlevs
    character(len=256) :: header, sd_name, sd_name_nodal

    real(dp_t) :: writetime1, writetime2

    namelist /chkpoint/ nlevs

    type(bl_prof_timer), save :: bpt

    call build(bpt, "checkpoint_write")

    if ( parallel_IOProcessor() ) then
       call fabio_mkdir(dirname)
    end if

    call parallel_barrier()

    writetime1 = parallel_wtime()

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

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
      write(6,*) 'Writing state to checkpoint file ',trim(sd_name)
    end if

    if (plot_Hext) then
       write(unit=sd_name, fmt='(a,"/rho_Hext")') trim(dirname)
       call fabio_ml_multifab_write_d(rho_Hext, rrs(:,1), sd_name, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)
       
       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) 'Writing state to checkpoint file ',trim(sd_name)
       end if
    endif

    if (use_thermal_diffusion) then
       write(unit=sd_name, fmt='(a,"/thermal2")') trim(dirname)
       call fabio_ml_multifab_write_d(thermal2, rrs(:,1), sd_name, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)

       if (parallel_IOProcessor() .and. verbose .ge. 1) then
          write(6,*) 'Writing state to checkpoint file ',trim(sd_name)
       end if
    end if

    write(unit=sd_name_nodal, fmt='(a,"/Pressure")') trim(dirname)
    call fabio_ml_multifab_write_d(mfs_nodal, rrs(:,1), sd_name_nodal, nOutFiles = nOutFiles, lUsingNFiles = lUsingNFiles)

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
      write(6,*) 'Writing state to checkpoint file ',trim(sd_name_nodal)
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

    if (parallel_IOProcessor()) then
       un = unit_new()
       open(unit=un, file=trim(dirname) // "/CPUtime", &
            form="formatted", action="write", status="replace")
       write(unit=un,fmt=*) get_cputime()
       close(un)
    endif

    writetime2 = parallel_wtime() - writetime1
    call parallel_reduce(writetime1, writetime2, MPI_MAX, proc=parallel_IOProcessorNode())
    if (parallel_IOProcessor()) then
       print*,'Time to write checkpoint: ',writetime1,' seconds'
       print*,''
    end if

    call destroy(bpt)

1000 format(32(e30.20,1x))

  end subroutine checkpoint_write

  subroutine checkpoint_read(mfs, mfs_nodal, dSdt, Source_old, Source_new, rho_omegadot2, &
                             rho_Hnuc2, rho_Hext, thermal2, dirname, dt_out, nlevs_out)

    use bl_IO_module, only: unit_new
    use fabio_module, only: fabio_ml_multifab_read_d
    use bl_prof_module, only: bl_prof_timer, build, destroy
    use variables, only: rel_eps
    use probin_module, only: use_thermal_diffusion, plot_Hext
    use time_module, only: time

    type(multifab  ),                pointer :: mfs(:), mfs_nodal(:)
    type(multifab  ),                pointer :: dSdt(:), Source_old(:), Source_new(:)
    type(multifab  ),                pointer :: rho_omegadot2(:), rho_Hnuc2(:), rho_Hext(:), thermal2(:)
    character(len=*), intent(in   )          :: dirname
    integer         , intent(  out)          :: nlevs_out
    real(kind=dp_t) , intent(  out)          :: dt_out

    ! local
    integer            :: un
    character(len=256) :: header, sd_name
    integer            :: nlevs
    real(kind=dp_t)    :: dt_in

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

    read(unit=un,fmt=*) dt_in
    read(unit=un,fmt=*) time
    read(unit=un,fmt=*) rel_eps
    close(un)

       dt_out = dt_in
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

!   Read the rho_Hnuc2 data into a multilevel multifab.
    if(plot_Hext) then
       write(unit=sd_name, fmt='(a,"/rho_Hext")') trim(dirname)
       call fabio_ml_multifab_read_d(rho_Hext, sd_name)
    endif

!   Read the thermal2 data into a multilevel multifab.
    if (use_thermal_diffusion) then
       write(unit=sd_name, fmt='(a,"/thermal2")') trim(dirname)
       call fabio_ml_multifab_read_d(thermal2, sd_name)
    end if

    call destroy(bpt)

  end subroutine checkpoint_read

end module checkpoint_module
