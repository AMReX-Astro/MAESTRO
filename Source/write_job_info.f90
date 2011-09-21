subroutine write_job_info(dirname, mba)

  ! write out some basic information about the way the job was run
  ! to a file called job_info in the directory dir_name.  Usually
  ! dir_name will be the name of the checkpoint or plotfile toplevel
  ! directory

  use parallel
  use probin_module, only: job_name, probin, inputs_file_used
  use bl_system_module, only: BL_CWD_SIZE, get_cwd 
  use ml_boxarray_module
  use build_info_module, only: build_date, build_dir, build_machine, boxlib_dir, &
                               module_list, f90_compile_line, f_compile_line, &
                               C_compile_line, link_line
  use omp_module

  implicit none

  character (len=*), intent(in) :: dirname
  type(ml_boxarray), intent(in) :: mba



  character (len=256) :: out_name
  character (len=16) :: date_in, time_in
  integer, dimension(8) :: values
  character (len=BL_CWD_SIZE) :: cwd

  integer :: i, n

  call date_and_time(date_in, time_in, VALUES=values)
  call get_cwd(cwd)
 
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
     write (99,1001) "job name:    ", trim(job_name)
     write (99,1001) "inputs file: ", trim(inputs_file_used)
     write (99,*) " "     
     write (99,1002) "number of MPI processes ", parallel_nprocs()
     write (99,1002) "number of threads       ", omp_get_max_threads()
     write (99,1003) "output date:            ", values(1), values(2), values(3)
     write (99,1004) "output time:            ", values(5), values(6), values(7)
     write (99,1001) "output dir:             ", trim(cwd)

     write (99,*) " "
     write (99,*) " "


     write (99,*) "Build Information"
     write (99,1000)
     write (99,1001) "build date:    ", trim(build_date)
     write (99,1001) "build machine: ", trim(build_machine)
     write (99,1001) "build dir:     ", trim(build_dir)
     write (99,1001) "BoxLib dir:    ", trim(boxlib_dir)
     write (99,*) " "
     write (99,1001) "modules used:  ", trim(module_list)
     write (99,*) " "
     write (99,1001) "F90 compile line: ", trim(f90_compile_line)
     write (99,*) " "
     write (99,1001) "F77 compile line: ", trim(f_compile_line)
     write (99,*) " "     
     write (99,1001) "C compile line:   ", trim(C_compile_line)
     write (99,*) " "
     write (99,1001) "linker line:      ", trim(link_line)


     write (99,*) " "
     write (99,*) " "


     write (99,*) "Grid Information"
     write (99,1000)
     do n = 1, mba%nlevel
        write (99,*) "level: ", n
        write (99,*) "   number of boxes = ", nboxes(mba, n)
        write (99,*) "   maximum zones   = ", (extent(mba%pd(n),i),i=1,mba%dim)
     end do

     write (99,*) " "
     write (99,*) " "


     write (99,*) "Runtime Parameter Information"
     write (99,1000)
     write (99,nml=probin)
     close(99)
  endif

end subroutine write_job_info
  
