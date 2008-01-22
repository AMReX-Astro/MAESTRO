module base_io_module

  use bl_types

  implicit none

  private

  public :: write_base_state, read_base_state

contains

  subroutine write_base_state(nlevs,state_name,w0_name,chk_name,s0,p0,gam1,w0, &
                              div_coeff,problo)
    
    use parallel
    use bl_prof_module
    use geometry, only : dr, nr
    use network, only: nspec
    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp
    use bl_constants_module
    use probin_module, only: prob_lo_y, prob_lo_z

    integer          , intent(in) :: nlevs
    character(len=10), intent(in) :: state_name
    character(len=7) , intent(in) :: w0_name
    character(len=7) , intent(in) :: chk_name
    real(kind=dp_t)  , intent(in) :: s0(:,:),p0(:),gam1(:),div_coeff(:), w0(:)

    real(kind=dp_t) :: base_r, problo
    character(len=18) :: out_name
    integer :: i, n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "write_base_state")

    if (parallel_IOProcessor()) then

       ! write out the base state quantities
       out_name = chk_name // "/" // state_name
       write(6,*) 'Writing base state to ',out_name

       open(unit=99,file=out_name,form = "formatted", access = "sequential",action="write")
       do i = 1, nr(1)
          base_r = problo + (dble(i)-HALF) * dr(1)
          write(99,1000)  base_r,s0(i,rho_comp), p0(i), gam1(i), s0(i,rhoh_comp), &
               (s0(i,n), n=spec_comp,spec_comp+nspec-1), s0(i,temp_comp), div_coeff(i)
       end do
       close(99)

       ! write out w0 (it is nodal, so it gets a separate file)
       out_name = chk_name // "/" // w0_name
       write(6,*) 'Writing w0 state to ',out_name
       write(6,*) ''

       open(unit=99,file=out_name,form = "formatted", access = "sequential",action="write")
       do i = 1, nr(1)+1
          base_r = problo + (dble(i)-1) * dr(1)
          write(99,1000)  base_r,w0(i)
       end do
       close(99)

    endif

    call destroy(bpt)

1000 format(32(e30.20,1x))

  end subroutine write_base_state


  subroutine read_base_state(nlevs,state_name,w0_name,chk_name,s0,p0,gam1,w0,div_coeff)


    use parallel
    use bl_prof_module
    use variables, only: rho_comp, rhoh_comp, spec_comp, temp_comp
    use network, only: nspec
    use geometry, only : dr, nr
    use bl_constants_module
    
    integer          , intent(in   ) :: nlevs
    character(len=10), intent(in   ) :: state_name
    character(len=7) , intent(in   ) :: w0_name
    character(len=7) , intent(in   ) :: chk_name    
    real(kind=dp_t)  , intent(inout) :: s0(:,:),p0(:),gam1(:),div_coeff(:),w0(:)
    real(kind=dp_t)  , allocatable   :: base_r(:)

    real(kind=dp_t) :: r_dummy
    character(len=18) :: out_name
    integer :: i, n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "read_base_state")

    allocate(base_r(nr(1)))

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

    call destroy(bpt)

  end subroutine read_base_state

end module base_io_module

