program main
  use BoxLib
  use layout_module

  implicit none

  call boxlib_initialize()

  !call layout_set_verbosity(1)

  call varden()

  call boxlib_finalize()

end program main
