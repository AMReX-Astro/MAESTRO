subroutine write_job_info(dirname)

  ! write out some basic information about the way the job was run
  ! to a file called job_info in the directory dir_name.  Usually
  ! dir_name will be the name of the checkpoint or plotfile toplevel
  ! directory

  use parallel

  implicit none
  character (len=*)   :: dirname
  character (len=256) :: build_date, build_dir, build_machine
  character (len=64) :: out_name

  call build_info(build_date, build_dir, build_machine)

  out_name = dirname // "/job_info"

1001 format(a,a)

  if (parallel_IOProcessor()) then
     open(unit=99,file=out_name,form = "formatted", access = "sequential",action="write")

     write (99,1001) "build date:    ", trim(build_date)
     write (99,1001) "build machine: ", trim(build_machine)
     write (99,1001) "build dir:     ", trim(build_dir)

     close(99)
  endif

end subroutine write_job_info
  
