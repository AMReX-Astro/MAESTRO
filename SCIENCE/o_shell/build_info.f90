subroutine build_info(build_date, build_dir, build_machine)

  implicit none

! the size of the strings is set to be the same as the size used in the
! cut in the shell script, so we don't go out of bounds

  character (len=128) :: build_date, build_dir, build_machine

  build_date    = &
   "Tue Apr 20 14:40:02 PDT 2010"

  build_dir     = &
    "/home/gilet/MAESTRO/fParallel/MAESTRO/o_shell"

  build_machine = &
    "Linux orga 2.6.28-18-generic #60-Ubuntu SMP Fri Mar 12 04:26:47 UTC 2010 x86_64 GNU/Linux"

  return
end subroutine build_info

