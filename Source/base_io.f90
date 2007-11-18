module base_io_module

  use bl_types
  use bl_constants_module
  use parallel
  use variables
  use network
  use geometry

  implicit none

contains

  subroutine write_base_state(state_name,w0_name,chk_name,s0,p0,gam1,w0,div_coeff)

    character(len=10), intent(in) :: state_name
    character(len=7), intent(in) :: w0_name
    character(len=7), intent(in) :: chk_name
    real(kind=dp_t) , intent(in) :: s0(:,:),p0(:),gam1(:),div_coeff(:), w0(:)
    real(kind=dp_t) :: base_r

    character(len=18) :: out_name
    integer :: i, n, nr

    nr = size(s0,dim=1)


    if (parallel_IOProcessor()) then

       ! write out the base state quantities
       out_name = chk_name // "/" // state_name
       write(6,*) 'Writing base state to ',out_name

       open(unit=99,file=out_name,form = "formatted", access = "sequential",action="write")
       do i = 1, nr
          base_r = (dble(i)-HALF) * dr
          write(99,1000)  base_r,s0(i,rho_comp), p0(i), gam1(i), s0(i,rhoh_comp), &
               (s0(i,n), n=spec_comp,spec_comp+nspec-1), s0(i,temp_comp), div_coeff(i)
       end do
       close(99)

       ! write out w0 (it is nodal, so it gets a separate file)
       out_name = chk_name // "/" // w0_name
       write(6,*) 'Writing w0 state to ',out_name
       write(6,*) ''

       open(unit=99,file=out_name,form = "formatted", access = "sequential",action="write")
       do i = 1, nr+1
          base_r = (dble(i)-1) * dr
          write(99,1000)  base_r,w0(i)
       end do
       close(99)

    endif


1000 format(32(e30.20,1x))

  end subroutine write_base_state


  subroutine read_base_state(state_name,w0_name,chk_name,s0,p0,gam1,w0,div_coeff)
    
    character(len=10), intent(in   ) :: state_name
    character(len=7) , intent(in   ) :: w0_name
    character(len=7) , intent(in   ) :: chk_name    
    real(kind=dp_t) , intent(inout) :: s0(:,:),p0(:),gam1(:),div_coeff(:),w0(:)
    real(kind=dp_t) , allocatable   :: base_r(:)

    real(kind=dp_t) :: r_dummy
    character(len=18) :: out_name
    integer :: i, n, nr

    nr = size(s0,dim=1)
    allocate(base_r(nr))

    ! read in the state variables
    out_name = chk_name // "/" // state_name
    if (parallel_IOProcessor()) then
      print *,'Reading base state from ',out_name
    end if

    open(unit=99,file=out_name)
    do i = 1, size(s0,dim=1)
       read(99,*)  base_r(i),s0(i,rho_comp), p0(i), gam1(i),s0(i,rhoh_comp), &
                   (s0(i,n), n=spec_comp,spec_comp+nspec-1), s0(i,temp_comp), div_coeff(i)
    end do
    close(99)


    ! read in w0
    out_name = chk_name // "/" // w0_name
    if (parallel_IOProcessor()) then
      print *,'Reading w0 state from ',out_name
    end if

    open(unit=99,file=out_name)
    do i = 1, size(w0,dim=1)
       read(99,*)  r_dummy, w0(i)
    end do
    close(99)

    deallocate(base_r)

  end subroutine read_base_state

end module base_io_module

