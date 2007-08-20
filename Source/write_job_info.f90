subroutine write_job_info(dirname)

  ! write out some basic information about the way the job was run
  ! to a file called job_info in the directory dir_name.  Usually
  ! dir_name will be the name of the checkpoint or plotfile toplevel
  ! directory

  use parallel
  use probin_module

  implicit none
  character (len=*)   :: dirname
  character (len=256) :: build_date, build_dir, build_machine
  character (len=256) :: out_name
  character (len=16) :: date, time
  integer, dimension(8) :: values

  call build_info(build_date, build_dir, build_machine)
  call date_and_time(date, time, VALUES=values)

  out_name = trim(dirname) // "/job_info"

1000 format(79('-'))
1001 format(a,a)
1002 format(a,i6)
1003 format(a,i4.4,'-',i2.2,'-',i2.2)
1004 format(a,i2.2,':',i2.2,':',i2.2)

  if (parallel_IOProcessor()) then
     open(unit=99,file=out_name,form = "formatted", access = "sequential",action="write")
     
     write (99,*) "Job Information"
     write (99,1000)
     write (99,1002) "number of processors: ", parallel_nprocs()
     write (99,1003) "output date:          ", values(1), values(2), values(3)
     write (99,1004) "output time:          ", values(5), values(6), values(7)
     write (99,*) " "
     write (99,*) " "

     write (99,*) "Build Information"
     write (99,1000)
     write (99,1001) "build date:    ", trim(build_date)
     write (99,1001) "build machine: ", trim(build_machine)
     write (99,1001) "build dir:     ", trim(build_dir)
     write (99,*) " "
     write (99,*) " "

     write (99,*) "Runtime Parameter Information"
     write (99,1000)
     write (99,nml=probin)
     close(99)
  endif

end subroutine write_job_info
  
