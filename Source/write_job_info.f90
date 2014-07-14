subroutine write_job_info(dirname, mba, the_bc_tower, write_pf_time)

  ! write out some basic information about the way the job was run
  ! to a file called job_info in the directory dir_name.  Usually
  ! dir_name will be the name of the checkpoint or plotfile toplevel
  ! directory

  use bl_types
  use parallel
  use probin_module, only: job_name, inputs_file_used
  use runtime_init_module, only: runtime_pretty_print
  use bl_system_module, only: BL_CWD_SIZE, get_cwd 
  use bc_module
  use define_bc_module
  use ml_boxarray_module
  use build_info_module, only: build_date, build_dir, build_machine, boxlib_dir, &
                               NUM_MODULES, modules, FCOMP, FCOMP_version, &
                               f90_compile_line, f_compile_line, &
                               C_compile_line, link_line, &
                               source_git_hash, boxlib_git_hash, extra_git_hash, &
                               different_build_tree, build_git_hash
  use omp_module
  use network
  use cputime_module, only: get_cputime
  use mg_eps_module
  use simple_log_module

  implicit none

  character (len=*), intent(in) :: dirname
  type(ml_boxarray), intent(in) :: mba
  type(bc_tower)   , intent(in) :: the_bc_tower
  real(kind=dp_t)  , intent(in) :: write_pf_time

  character (len=256) :: out_name
  character (len=16) :: date_in, time_in
  integer, dimension(8) :: values
  character (len=BL_CWD_SIZE) :: cwd

  integer :: i, n

  logical :: is_overview 

  call date_and_time(date_in, time_in, VALUES=values)
  call get_cwd(cwd)

  if (dirname == "") then
     out_name = "maestro-overview.out"
     is_overview = .true.
  else
     out_name = trim(dirname) // "/job_info"
     is_overview = .false.
  endif

999  format(79('='))
1000 format(79('-'))
1001 format(a,a)
1002 format(a,i6)
1003 format(a,i4.4,'-',i2.2,'-',i2.2)
1004 format(a,i2.2,':',i2.2,':',i2.2)
1005 format(a,g20.10)
2001 format(a5,1x,a20,1x,a20,1x,a8,1x,a8)
2002 format(i5,1x,a20,1x,a20,1x,f8.2,1x,f8.2)

  if (parallel_IOProcessor()) then
     open(unit=99,file=out_name,form = "formatted", access = "sequential",action="write")
     
     write (99,999)
     write (99,*) "Job Information"
     write (99,999)
     write (99,1001) "job name:    ", trim(job_name)
     write (99,1001) "inputs file: ", trim(inputs_file_used)
     write (99,*) " "     
     write (99,1002) "number of MPI processes ", parallel_nprocs()
     write (99,1002) "number of threads       ", omp_get_max_threads()
     write (99,*) " "
     write (99,1005) "CPU time used since start of simulation (CPU-hours) ", get_cputime()/3600.0_dp_t

     write (99,*) " "
     write (99,*) " "

     write (99,999)
     write (99,*) "Plotfile Information"
     write (99,999)
     write (99,1003) "output date:              ", values(1), values(2), values(3)
     write (99,1004) "output time:              ", values(5), values(6), values(7)
     write (99,1001) "output dir:               ", trim(cwd)
     write (99,1005) "time to write plotfile (s)", write_pf_time

     write (99,*) " "
     write (99,*) " "


     write (99,999)
     write (99,*) "Build Information"
     write (99,999)
     write (99,1001) "build date:    ", trim(build_date)
     write (99,1001) "build machine: ", trim(build_machine)
     write (99,1001) "build dir:     ", trim(build_dir)
     write (99,1001) "BoxLib dir:    ", trim(boxlib_dir)
     write (99,*) " "
     write (99,1001) "MAESTRO    git hash: ", trim(source_git_hash)
     write (99,1001) "BoxLib     git hash: ", trim(boxlib_git_hash)
     write (99,1001) "AstroDev   git hash: ", trim(extra_git_hash)
     if (different_build_tree) then
        write (99,1001) "build tree git hash: ", trim(build_git_hash)     
     endif

     write (99,*) " "
     write (99,1001) "modules used:  ", " "
     do i=1, NUM_MODULES
        write (99,1001) "  ", trim(modules(i))
     enddo
     write (99,*) " "
     write (99,1001) "FCOMP:            ", trim(FCOMP)
     write (99,1001) "FCOMP version:    ", trim(FCOMP_version)
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


     write (99,999)
     write (99,*) "Grid Information"
     write (99,999)
     do n = 1, mba%nlevel
        write (99,*) "level: ", n
        write (99,*) "   number of boxes = ", nboxes(mba, n)
        write (99,*) "   maximum zones   = ", (extent(mba%pd(n),i),i=1,mba%dim)
     end do

     write (99,*) " "
     write (99,*) "Boundary Conditions"     
     write (99,*) "  -x: ", bc_integer_to_string(the_bc_tower%bc_tower_array(1)%phys_bc_level_array(0,1,1))
     write (99,*) "  +x: ", bc_integer_to_string(the_bc_tower%bc_tower_array(1)%phys_bc_level_array(0,1,2))

     if (mba%dim >= 2) then
        write (99,*) " "
        write (99,*) "  -y: ", bc_integer_to_string(the_bc_tower%bc_tower_array(1)%phys_bc_level_array(0,2,1))
        write (99,*) "  +y: ", bc_integer_to_string(the_bc_tower%bc_tower_array(1)%phys_bc_level_array(0,2,2))
     endif

     if (mba%dim == 3) then
        write (99,*) " "
        write (99,*) "  -z: ", bc_integer_to_string(the_bc_tower%bc_tower_array(1)%phys_bc_level_array(0,3,1))
        write (99,*) "  +z: ", bc_integer_to_string(the_bc_tower%bc_tower_array(1)%phys_bc_level_array(0,3,2))
     endif

     write (99,*) " "
     write (99,*) " "


     write (99,999)
     write (99,*) "Species Information"
     write (99,999)
     write (99,2001) "index", "name", "short name", "A", "Z"
     write (99,1000)
     do n = 1, nspec
        write (99,2002) n, spec_names(n), short_spec_names(n), aion(n), zion(n)
     enddo

     write (99,*) " "
     write (99,*) " "


     write (99,999)
     write (99,*) "MG Tolerance Information"
     write (99,999)
     write (99,*) "  eps_mac:          ", eps_mac
     write (99,*) "  eps_mac_max:      ", eps_mac_max
     write (99,*) "  mac_level_factor: ", mac_level_factor
     write (99,*) "  eps_mac_bottom:   ", eps_mac_bottom

     write (99,*) " "

     write (99,*) "  eps_hg:           ", eps_hg
     write (99,*) "  eps_hg_max:       ", eps_hg_max
     write (99,*) "  hg_level_factor:  ", hg_level_factor
     write (99,*) "  eps_hg_bottom:    ", eps_hg_bottom

     write (99,*) " "
     write (99,*) " "


     write (99,999)
     write (99,*) "Initialization Logs/Warnings"
     write (99,999)

     do i = 1, log_lines
        write (99,*) trim(log_data(i))
     enddo

     write (99,999)
     write (99,*) "Runtime Parameter Information"
     write (99,999)
     call runtime_pretty_print(99)

     if (is_overview) then
        write (99,*) " "
        write (99,*) "Restart information: "
     endif

     close(99)

  endif

end subroutine write_job_info
  
