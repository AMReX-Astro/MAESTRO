module eglib_module

  implicit none

  integer, save :: EGLIB_NP=-1
  integer, save :: legwork, legiwork
  double precision, allocatable, save :: egwork(:) 
  integer, allocatable, save :: egiwork(:) 

!$omp threadprivate(egwork,egiwork)

  private

  public egwork, egiwork, legwork, legiwork, eglib_init, eglib_close

contains

  subroutine eglib_init(nspecies, NP, ITLS, IFLAG)
    integer, intent(in) :: nspecies, NP, ITLS, IFLAG

    if (NP .ne. EGLIB_NP) then

       call eglib_close()

       EGLIB_NP = NP

       ! for ITLS=1
       if (ITLS .eq. 1) then
          legwork = 23 + 14*nspecies + 32*nspecies**2 + 13*NP &
               + 14*NP*nspecies + NP*nspecies**2
       else if (ITLS .eq. 2) then
          legwork = 23 + 14*nspecies + 32*nspecies**2 + 13*NP &
               + 21*NP*nspecies + NP*(2*nspecies**2+(nspecies*(nspecies+1))/2)
       else
          legwork = 23 + 14*nspecies + 32*nspecies**2 + 13*NP  & 
               + 30*NP*nspecies + 5*NP*nspecies**2
       end if
       
       legiwork = nspecies
       
       !$omp parallel 
       allocate(egwork(legwork))
       allocate(egiwork(legiwork))
       call egini(NP, 6, IFLAG, ITLS, egwork, legwork, egiwork, legiwork)
       !$omp end parallel
       
    end if

  end subroutine eglib_init


  subroutine eglib_close()
    if (EGLIB_NP > 0) then
       !$omp parallel 
       deallocate(egwork)
       deallocate(egiwork)
       !$omp end parallel
       EGLIB_NP = -1
    end if
  end subroutine eglib_close

end module eglib_module
