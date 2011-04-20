! the aux data module is used to store auxillary data (i.e. for
! diagnostic purposes) that is needed to persist through
! checkpoint/restarts.

! NOTE: only the IOProcessor's data is stored and read-back -- it
! is assumed that the other processors have the same data.

module aux_data_module

  use bl_types

  implicit none

  private

  integer, save, public :: naux = 0
  real (kind=dp_t), allocatable, public, save :: aux_data(:)

  public :: init_aux_data, write_aux_data, read_aux_data

contains

  subroutine init_aux_data(naux_pass)

    integer, intent(in) :: naux_pass

    naux = naux_pass
    allocate (aux_data(naux))

  end subroutine init_aux_data


  subroutine write_aux_data(istep,chk_name)

    ! the IOProcessor will write the aux_data to the checkpoint
    ! directory

    use parallel, only: parallel_IOProcessor
    use bl_IO_module

    integer          , intent(in) :: istep
    character(len=*) , intent(in) :: chk_name

    integer :: n, un
    character (len=256) :: aux_name, out_name

    if (naux == 0) return

    ! create the names of the files that will store the output
    if (istep <= 99999) then
       write(unit=aux_name,fmt='("aux_",i5.5)') istep
    else
       write(unit=aux_name,fmt='("aux_",i6.6)') istep
    endif
  

    if (parallel_IOProcessor()) then

       out_name = trim(chk_name) // "/" // trim(aux_name)

       print *, "Writing auxillary data to ", trim(out_name)

       un = unit_new()
       open(unit=un,file=out_name,form = "formatted", access = "sequential",action="write")

       write(un,*) naux
       do n = 1, naux
          write (un,*) aux_data(n)
       enddo

       close(un)

    end if

  end subroutine write_aux_data


  subroutine read_aux_data(restart,chk_name)

    ! each processor will read in the data directly

    use parallel, only: parallel_IOProcessor
    use bl_IO_module

    integer          , intent(in   ) :: restart
    character(len=*) , intent(in   ) :: chk_name    

    integer :: n, un
    logical :: lexist
    character (len=256) :: aux_name, out_name

    if (restart <= 99999) then
       write(unit=aux_name,fmt='("aux_",i5.5)') restart
    else
       write(unit=aux_name,fmt='("aux_",i6.6)') restart
    endif

    out_name = trim(chk_name) // "/" // trim(aux_name)

    ! check to see if the aux file exists.  If it does not, then there
    ! is no auxillary data, and we do nothing.
    inquire(file=trim(out_name), exist=lexist)
    
    if (lexist) then

       if (parallel_IOProcessor()) then
          print *, "Reading auxillary data from ", trim(out_name)
       endif

       un = unit_new()
       open(unit=un,file=out_name)

       read (un,*) naux
       allocate (aux_data(naux))
       
       do n = 1, naux
          read (un,*) aux_data(n)
       enddo

    endif

  end subroutine read_aux_data

end module aux_data_module
