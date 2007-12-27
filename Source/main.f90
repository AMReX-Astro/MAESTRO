program main

  use BoxLib
  use layout_module
  use bl_prof_module

  implicit none

  call boxlib_initialize()
  call bl_prof_initialize(on = .true.)

  !call layout_set_verbosity(1)

  call varden()

  call bl_prof_glean("bl_prof_res")
  call bl_prof_finalize()
  call boxlib_finalize()

end program main
