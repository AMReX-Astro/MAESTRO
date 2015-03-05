module tranlib_module

  implicit none

  integer, save :: lmcwork, lmciwork
  double precision, allocatable, save :: mcwork(:)
  integer, allocatable, save :: mciwork(:)

  integer, parameter :: maxtp = 3
  integer, parameter :: LLINKMC = 44, LOUT = 6

!$omp threadprivate(mcwork,mciwork)

  private

  public mcwork, mciwork, lmcwork, lmciwork, tranlib_init, tranlib_close

contains

  subroutine tranlib_init(nspecies)
    integer, intent(in) :: nspecies
    integer :: ierr
    integer :: MAXFIT, NO, NFDIM, NT, NRANGE, NLITEMAX

    MAXFIT=7
    NO=4
    NFDIM=165
    NT=50
    NRANGE = MAXTP-1
    NLITEMAX=2
    lmcwork = nspecies*(19+2*NO+NO*NLITEMAX)+(NO+15)*nspecies**2
    lmciwork = 4*nspecies + NLITEMAX

    !$omp parallel

    allocate(mcwork(lmcwork))
    allocate(mciwork(lmciwork))

    CALL MCINITCD (LLINKMC, LOUT, lmciwork, lmcwork, mciwork, mcwork, ierr)
    if (ierr .gt. 0) then
       WRITE(LOUT,*)' QUITTING BECAUSE MCINIT IFLAG = ', ierr
       stop
    end if

    !$omp end parallel

  end subroutine tranlib_init

  subroutine tranlib_close()
    !$omp parallel
    deallocate(mcwork, mciwork)
    !$omp end parallel
  end subroutine tranlib_close

end module tranlib_module
