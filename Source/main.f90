program main

  use BoxLib
  use parallel
  use layout_module
  use bl_prof_module

  implicit none

  real(dp_t) :: r1, r2

  call boxlib_initialize()

  r1 = parallel_wtime()

  call bl_prof_initialize(on = .true.)

  !call layout_set_verbosity(1)

  call varden()
  !
  ! TODO -- add ability to specify filename via inputs file.
  !
  call bl_prof_glean("bl_prof_res")

  call bl_prof_finalize()

  r2 = parallel_wtime() - r1

  call parallel_reduce(r1, r2, MPI_MAX, proc = parallel_IOProcessorNode())

  if (parallel_IOProcessor()) print*, 'Run Time = ', r1

  call boxlib_finalize()

end program main
