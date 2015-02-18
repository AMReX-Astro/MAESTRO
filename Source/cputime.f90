! keep track of the CPU time since the simulation began.
!
! Note: we only expect this information to be valid on the IOProcessor

module cputime_module
  
  use bl_types
  use bl_constants_module
  use bl_error_module, only: bl_warn
  use parallel
  use omp_module

  implicit none

  ! previous_elapsed_cputime is zero at the start of a simulation.
  ! When restarting from a checkpoint, previous_elapsed_cputime will
  ! hold the CPU time used before now.
  real (kind=dp_t), save :: previous_elapsed_cputime = ZERO
  real (kind=dp_t), save :: start_cputime
  logical,          save :: initialized = .FALSE.

  private

  public :: initialize_elapsed_cputime, start_cputime_clock, get_cputime
contains

  subroutine initialize_elapsed_cputime(cputime)

    real (kind=dp_t) :: cputime

    previous_elapsed_cputime = cputime

  end subroutine initialize_elapsed_cputime


  subroutine start_cputime_clock()

    start_cputime = parallel_wtime()
    initialized = .TRUE.
    
  end subroutine start_cputime_clock


  function get_cputime() result (time)

    real (kind=dp_t) :: time
    integer :: ncores

    if (.not. initialized) then
      call start_cputime_clock() 
      call bl_warn("WARNING! You called get_cputime() without first calling start_cputime_clock()")
    endif
    ncores = parallel_nprocs() * omp_get_max_threads()
    time = ncores*(parallel_wtime() - start_cputime) + previous_elapsed_cputime

    return

  end function get_cputime

end module cputime_module
