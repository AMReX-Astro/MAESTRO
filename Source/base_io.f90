module base_io_module

  use bl_types

  implicit none

  private

  public :: write_base_state, read_base_state

contains

  subroutine write_base_state(nlevs,state_name,w0_name,eta_name,chk_name, &
                              s0,p0,gam1,w0,eta,div_coeff,problo)
    
    use parallel
    use bl_prof_module
    use geometry, only : dr, nr
    use network, only: nspec
    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp
    use bl_constants_module
    use probin_module, only: prob_lo_y, prob_lo_z

    integer          , intent(in) :: nlevs
    character(len=11), intent(in) :: state_name
    character(len=8) , intent(in) :: w0_name
    character(len=9) , intent(in) :: eta_name
    character(len=8) , intent(in) :: chk_name
    real(kind=dp_t)  , intent(in) :: s0(:,:,:),p0(:,:),gam1(:,:),div_coeff(:,:)
    real(kind=dp_t)  , intent(in) :: w0(:,:),eta(:,:,:)

    real(kind=dp_t) :: base_r, problo
    character(len=20) :: out_name
    integer :: i, comp, n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "write_base_state")

    if (parallel_IOProcessor()) then

       print*,"chk_name",chk_name
       print*,"state_name",state_name

       ! write out the base state quantities
       out_name = chk_name // "/" // state_name
       write(6,*) 'Writing base state to ',out_name

       open(unit=99,file=out_name,form = "formatted", access = "sequential",action="write")

       do n=1,nlevs
          do i=1,nr(n)
             base_r = problo + (dble(i)-HALF) * dr(n)
             write(99,1000)  base_r,s0(n,i,rho_comp), p0(n,i), gam1(n,i), &
                  s0(n,i,rhoh_comp), &
                  (s0(n,i,comp), comp=spec_comp,spec_comp+nspec-1), &
                  s0(n,i,temp_comp), div_coeff(n,i)
          end do
       end do
       close(99)

       ! write out w0 (it is edge-based, so it gets a separate file)
       out_name = chk_name // "/" // w0_name
       write(6,*) 'Writing w0 state to ',out_name
       write(6,*) ''

       open(unit=99,file=out_name,form = "formatted", access = "sequential",action="write")
       do n=1,nlevs
          do i=1,nr(n)+1
             base_r = problo + (dble(i)-1) * dr(n)
             write(99,1000)  base_r,w0(n,i)
          end do
       end do
       close(99)

       ! write out eta (it is edge-based, so it gets a separate file)
       out_name = chk_name // "/" // eta_name
       write(6,*) 'Writing eta on edges to ',out_name
       write(6,*) ''

       open(unit=99,file=out_name,form = "formatted", access = "sequential",action="write")
       do n=1,nlevs
          do i=1,nr(n)+1
             base_r = problo + (dble(i)-1) * dr(n)
             write(99,1000)  base_r,eta(n,i,rho_comp), eta(n,i,rhoh_comp), &
                             (eta(n,i,comp), comp=spec_comp,spec_comp+nspec-1)
                               
          end do
       end do
       close(99)

    endif

    call destroy(bpt)

1000 format(32(e30.20,1x))

  end subroutine write_base_state


  subroutine read_base_state(nlevs,state_name,w0_name,eta_name,chk_name, &
                             s0,p0,gam1,w0,eta,div_coeff)

    use parallel
    use bl_prof_module
    use variables, only: rho_comp, rhoh_comp, spec_comp, temp_comp
    use network, only: nspec
    use geometry, only : dr, nr
    use bl_constants_module
    
    integer          , intent(in   ) :: nlevs
    character(len=11), intent(in   ) :: state_name
    character(len=8) , intent(in   ) :: w0_name
    character(len=9) , intent(in   ) :: eta_name
    character(len=8) , intent(in   ) :: chk_name    
    real(kind=dp_t)  , intent(inout) :: s0(:,:,:),p0(:,:),gam1(:,:),div_coeff(:,:)
    real(kind=dp_t)  , intent(inout) :: w0(:,:),eta(:,:,:)
    real(kind=dp_t)  , allocatable   :: base_r(:,:)

    real(kind=dp_t) :: r_dummy
    character(len=20) :: out_name
    integer :: i, comp, n

    type(bl_prof_timer), save :: bpt

    call build(bpt, "read_base_state")

    allocate(base_r(nlevs,nr(nlevs)))

    ! read in the state variables
    out_name = chk_name // "/" // state_name
    if (parallel_IOProcessor()) then
      print *,'Reading base state from ',out_name
    end if

    open(unit=99,file=out_name)

    do n=1,nlevs
       do i=1,nr(n)
          read(99,*)  base_r(n,i), s0(n,i,rho_comp), p0(n,i), gam1(n,i), s0(n,i,rhoh_comp), &
               (s0(n,i,comp), comp=spec_comp,spec_comp+nspec-1), s0(n,i,temp_comp), &
               div_coeff(n,i)
       end do
    end do
    close(99)

    ! read in w0
    out_name = chk_name // "/" // w0_name
    if (parallel_IOProcessor()) then
      print *,'Reading w0 state from ',out_name
    end if

    open(unit=99,file=out_name)
    do n=1,nlevs
       do i=1,nr(n)+1
          read(99,*)  r_dummy, w0(n,i)
       end do
    end do
    close(99)

    ! read in eta
    out_name = chk_name // "/" // eta_name
    if (parallel_IOProcessor()) then
      print *,'Reading eta state from ',out_name
    end if

    open(unit=99,file=out_name)
    do n=1,nlevs
       do i=1,nr(n)+1
          read(99,*)  r_dummy,eta(n,i,rho_comp), eta(n,i,rhoh_comp), &
                     (eta(n,i,comp), comp=spec_comp,spec_comp+nspec-1)
       end do
    end do
    close(99)

    deallocate(base_r)

    call destroy(bpt)

  end subroutine read_base_state

end module base_io_module

