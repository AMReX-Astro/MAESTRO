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

contains

  subroutine checkpoint_write(dirname, mfs, mfs_nodal, Source_nm1, Source_old, &
                              rrs, dx, time_in, dt_in, verbose)
    type(multifab), intent(in) :: mfs(:), mfs_nodal(:)
    type(multifab), intent(in) :: Source_nm1(:), Source_old(:)
    integer        , intent(in) :: rrs(:,:)
    real(kind=dp_t), intent(in) :: dx(:,:)
    character(len=*), intent(in) :: dirname
    real(kind=dp_t), intent(in) :: time_in, dt_in
    integer        , intent(in) :: verbose
    integer :: i, j, k, n
    character(len=128) :: header, sd_name, sd_name_nodal
    integer :: nc, un, nl, dm
    integer, allocatable ::  lo(:),  hi(:)
    integer :: idummy, rdummy
    type(box) :: lbbox

    integer         :: nlevs
    real(kind=dp_t) :: time, dt

    namelist /chkpoint/ time
    namelist /chkpoint/ dt
    namelist /chkpoint/ nlevs

    if ( parallel_IOProcessor() ) then
       call fabio_mkdir(dirname)
    end if

    write(unit=sd_name, fmt='(a,"/State")') trim(dirname)
    call fabio_ml_multifab_write_d(mfs, rrs(:,1), sd_name)

    write(unit=sd_name, fmt='(a,"/Source_nm1")') trim(dirname)
    call fabio_ml_multifab_write_d(Source_nm1, rrs(:,1), sd_name)

    write(unit=sd_name, fmt='(a,"/Source_old")') trim(dirname)
    call fabio_ml_multifab_write_d(Source_old, rrs(:,1), sd_name)

    write(unit=sd_name_nodal, fmt='(a,"/Pressure")') trim(dirname)
    call fabio_ml_multifab_write_d(mfs_nodal, rrs(:,1), sd_name_nodal)

    if (parallel_IOProcessor() .and. verbose .ge. 1) then
      print *,'Writing    state to checkpoint file ',trim(sd_name)
      print *,'Writing pressure to checkpoint file ',trim(sd_name_nodal)
      print *,' '
    end if
    
    nl = size(mfs)
    nc = ncomp(mfs(1))
    dm = mfs(1)%dim
    allocate(lo(dm),hi(dm))
    lbbox = bbox(get_boxarray(mfs(1)))

    idummy = 0
    rdummy = 0.0_dp_t
    lo = lwb(lbbox); hi = upb(lbbox)

    time = time_in
      dt =   dt_in

    header = "Header"
    un = unit_new()
    open(unit=un, &
         file = trim(dirname) // "/" // trim(header), &
         form = "formatted", access = "sequential", &
         status = "replace", action = "write")
    nlevs = size(mfs)
    write(unit=un, nml = chkpoint)
    do n = 1,nlevs
       write(unit=un,fmt=*) dx(n,1), dx(n,2)
    end do
    do n = 1,nlevs-1
       write(unit=un,fmt=*) rrs(n,1)
    end do
    
  end subroutine checkpoint_write

  subroutine checkpoint_read(mfs, mfs_nodal, Source_nm1, Source_old, dirname, time_out, dt_out, nlevs_out)
    use bl_IO_module
    type(multifab  ),                pointer :: mfs(:), mfs_nodal(:)
    type(multifab  ),                pointer :: Source_nm1(:), Source_old(:)
    character(len=*), intent(in   )          :: dirname
    integer         , intent(  out)          :: nlevs_out
    real(kind=dp_t) , intent(  out)          :: time_out, dt_out

    integer         ,                pointer :: rrs(:)
    real(kind=dp_t) ,                pointer :: dx(:,:)

    integer :: i, j, k, n
    character(len=128) :: header, sd_name
    integer :: nc, un, nl, dm
    integer :: idummy, rdummy
    type(box) :: lbbox

    integer         :: nlevs
    real(kind=dp_t) :: time, dt

    namelist /chkpoint/ nlevs
    namelist /chkpoint/ time
    namelist /chkpoint/ dt

!   First read the header information
    header = "Header"
    un = unit_new()
    open(unit=un, &
         file = trim(dirname) // "/" // trim(header), &
         status = "old", &
         action = "read")
    read(unit=un, nml = chkpoint)
    allocate( dx(nlevs,2))
    allocate(rrs(nlevs-1))
    do n = 1,nlevs
       read(unit=un,fmt=*) dx(n,1), dx(n,2)
    end do
    do n = 1,nlevs-1
       read(unit=un,fmt=*) rrs(n)
    end do
     time_out = time
       dt_out = dt
    nlevs_out = nlevs

!   Read the state data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/State")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs, sd_name)

!   Read the pressure data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/Pressure")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_nodal, sd_name)

!   Read the Source_nm1 data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/Source_nm1")') trim(dirname)
    call fabio_ml_multifab_read_d(Source_nm1, sd_name)

!   Read the Source_old data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/Source_old")') trim(dirname)
    call fabio_ml_multifab_read_d(Source_old, sd_name)
    
    nl = nlevs
    nc = ncomp(mfs(1))
    dm = mfs(1)%dim
    lbbox = bbox(get_boxarray(mfs(1)))

    deallocate(dx,rrs)

  end subroutine checkpoint_read

end module checkpoint_module
