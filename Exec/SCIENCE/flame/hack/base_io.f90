module base_io_module

  use bl_types

  implicit none

  private

  public :: write_base_state, read_base_state

contains

  subroutine write_base_state(istep,chk_name, &
                              rho0,rhoh0,p0,gamma1bar,w0,etarho_ec,etarho_cc, &
                              div_coeff,psi,tempbar,problo)
    
    use parallel
    use bl_prof_module
    use geometry, only : dr, nr, nlevs_radial
    use network, only: nspec
    use variables, only: rho_comp, rhoh_comp
    use bl_constants_module

    integer          , intent(in) :: istep
    character(len=*) , intent(in) :: chk_name
    real(kind=dp_t)  , intent(in) :: rho0(:,0:),rhoh0(:,0:)
    real(kind=dp_t)  , intent(in) :: p0(:,0:),gamma1bar(:,0:)
    real(kind=dp_t)  , intent(in) :: div_coeff(:,0:), psi(:,0:), tempbar(:,0:)
    real(kind=dp_t)  , intent(in) :: w0(:,0:)
    real(kind=dp_t)  , intent(in) :: etarho_ec(:,0:),etarho_cc(:,0:)

    character(len=256) :: state_name, w0_name
    real(kind=dp_t) :: base_r, problo
    character(len=256) :: out_name
    integer :: n, r, i

    type(bl_prof_timer), save :: bpt

    call build(bpt, "write_base_state")

    ! create the names of the files that will store the output
    if (istep <= 99999) then
       write(unit=state_name,fmt='("model_cc_",i5.5)') istep
       write(unit=w0_name,fmt='("model_ec_",i5.5)') istep
    else
       write(unit=state_name,fmt='("model_cc_",i6.6)') istep
       write(unit=w0_name,fmt='("model_ec_",i6.6)') istep
    endif
    
    if (parallel_IOProcessor()) then

       print*,"Writing 1D data to ", trim(chk_name)

       ! write out cell-centered 1D data
       out_name = trim(chk_name) // "/" // trim(state_name)
       write(6,*) 'Writing cell-centered 1D data to ', trim(out_name)

       open(unit=99,file=out_name,form = "formatted", access = "sequential",action="write")

       do n=1,nlevs_radial
          do r=0,nr(n)-1
             base_r = problo + (dble(r)+HALF) * dr(n)
             write(99,1000)  r, base_r, rho0(n,r), p0(n,r), gamma1bar(n,r), rhoh0(n,r), &
                  div_coeff(n,r), psi(n,r), tempbar(n,r), etarho_cc(n,r)
          end do
       end do
       close(99)

       ! write out edge-centered 1D data
       out_name = trim(chk_name) // "/" // trim(w0_name)
       write(6,*) 'Writing edge-centered 1D data to ', trim(out_name)

       open(unit=99,file=out_name,form = "formatted", access = "sequential",action="write")
       do n=1,nlevs_radial
          do r=0,nr(n)
             base_r = problo + dble(r) * dr(n)
             write(99,1000)  r, base_r, w0(n,r), etarho_ec(n,r)
          end do
       end do
       close(99)

       write(6,*)

    endif

    call destroy(bpt)

1000 format(i8,32(e30.20,1x))

  end subroutine write_base_state


  subroutine read_base_state(restart,chk_name, &
                             rho0,rhoh0,p0,gamma1bar,w0, &
                             etarho_ec,etarho_cc,div_coeff,psi,tempbar)

    use parallel
    use bl_prof_module
    use variables, only: rho_comp, rhoh_comp
    use network, only: nspec
    use geometry, only : dr, nr, nlevs_radial
    use bl_constants_module
    use restrict_base_module, only: fill_ghost_base
    use inlet_bc_module, only: set_initial_inlet_bcs
    
    integer          , intent(in   ) :: restart
    character(len=*) , intent(in   ) :: chk_name    
    real(kind=dp_t)  , intent(inout) :: rho0(:,0:),rhoh0(:,0:)
    real(kind=dp_t)  , intent(inout) :: p0(:,0:),gamma1bar(:,0:)
    real(kind=dp_t)  , intent(inout) :: div_coeff(:,0:), psi(:,0:), tempbar(:,0:)
    real(kind=dp_t)  , intent(inout) :: w0(:,0:)
    real(kind=dp_t)  , intent(inout) :: etarho_ec(:,0:),etarho_cc(:,0:)

    ! local
    character(len=256) :: state_name, w0_name
    real(kind=dp_t) :: r_dummy
    character(len=256) :: out_name
    integer :: r, n, i, r_dummy_int

    type(bl_prof_timer), save :: bpt

    call build(bpt, "read_base_state")

    if (restart <= 99999) then
       write(unit=state_name,fmt='("model_cc_",i5.5)') restart
       write(unit=w0_name,fmt='("model_ec_",i5.5)') restart
    else
       write(unit=state_name,fmt='("model_",i6.6)') restart
       write(unit=w0_name,fmt='("w0_",i6.6)') restart
    endif

    ! read in cell-centered 1D data
    out_name = trim(chk_name) // "/" // trim(state_name)
    if (parallel_IOProcessor()) then
      print *,'Reading cell-centered 1D data from ', trim(out_name)
    end if

    open(unit=99,file=out_name)

    do n=1,nlevs_radial
       do r=0,nr(n)-1
          read(99,*) r_dummy_int, r_dummy, rho0(n,r), p0(n,r), gamma1bar(n,r), &
               rhoh0(n,r), div_coeff(n,r), psi(n,r), tempbar(n,r), etarho_cc(n,r)
       end do
    end do
    close(99)

    ! read in edge-centered 1D data
    out_name = trim(chk_name) // "/" // trim(w0_name)
    if (parallel_IOProcessor()) then
      print *,'Reading edge-centered 1D data from ', trim(out_name)
    end if

    open(unit=99,file=out_name)
    do n=1,nlevs_radial
       do r=0,nr(n)
          read(99,*) r_dummy_int, r_dummy, w0(n,r), etarho_ec(n,r)
       end do
    end do
    close(99)

    ! set the inlet boundary condition parameters
    call set_initial_inlet_bcs()

    call destroy(bpt)

  end subroutine read_base_state

end module base_io_module

