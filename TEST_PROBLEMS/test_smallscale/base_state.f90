module base_state_module

  use bl_types
  use bl_constants_module
  use bc_module
  use setbc_module
  use define_bc_module
  use multifab_module
  use eos_module
  use variables
  use network
  use geometry

  implicit none

contains

  subroutine init_base_state (model_file,n_base,s0,p0,gam1,dx,prob_lo,prob_hi)

    character (len=256), intent(in) :: model_file ! I'm not using this anymore
    integer        ,     intent(in   ) :: n_base
    real(kind=dp_t),     intent(inout) ::    s0(0:,:)
    real(kind=dp_t),     intent(inout) ::    p0(0:)
    real(kind=dp_t),     intent(inout) ::  gam1(0:)
    real(kind=dp_t),     intent(in   ) :: prob_lo(:)
    real(kind=dp_t),     intent(in   ) :: prob_hi(:)
    real(kind=dp_t),     intent(in   ) :: dx(:)

    ! local
    integer ndum, i, j, dm, nspec

    parameter (ndum = 30)
    parameter (nspec = 3)

    character(len=128) :: lamsolfile
    real(kind=dp_t) :: state1d(ndum), Pamb, temporary
    real(kind=dp_t) :: loloc,hiloc,flameloc,qreact
    
    call helmeos_init

    dm = size(dx)

    lamsolfile = 'flame_4.e7_screen_left.out'

    ! first set the inflow boundary condition
    call asin1d(lamsolfile, -.00125d0, 0.d0, state1d, ndum, .false.)

    Pamb = state1d(18)
    p_eos(1) = Pamb

    den_eos(1) = state1d(3)
    temp_eos(1) = state1d(9)
    do j=1,nspec
       if(spec_names(j) .eq. "carbon-12") then
          xn_eos(1,j) = state1d(21)
       else if(spec_names(j) .eq. "magnesium-24") then
          xn_eos(1,j) = state1d(22)
       else if(spec_names(j) .eq. "oxygen-16") then
          xn_eos(1,j) = state1d(23)
       else
          print*,"In initdata, spec_names(",j,") invalid"
       endif
    enddo

    ! given P, T, and X, compute rho
    call eos(eos_input_tp, den_eos, temp_eos, &
             npts, nspec, &
             xn_eos, &
             p_eos, h_eos, e_eos, & 
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             do_diag)

    ! given rho, T, and X, compute h
    call eos(eos_input_rt, den_eos, temp_eos, &
             npts, nspec, &
             xn_eos, &
             p_eos, h_eos, e_eos, & 
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             do_diag)

    INLET_VN = 0.0d0
    INLET_VT = 0.0d0
    INLET_RHO = den_eos(1)
    if(use_big_h) then
       qreact = 0.0d0
       do j=1,nspec
          qreact = qreact + ebin(j)*xn_eos(1,j)
       enddo
       INLET_RHOH = den_eos(1)*(h_eos(1) + qreact)
    else
       INLET_RHOH = den_eos(1)*h_eos(1)
    endif
    do j=1,nspec
       if(spec_names(j) .eq. "carbon-12") then
          INLET_RHOC12 = den_eos(1)*xn_eos(1,j)
       else if(spec_names(j) .eq. "magnesium-24") then
          INLET_RHOMG24 = den_eos(1)*xn_eos(1,j)
       else if(spec_names(j) .eq. "oxygen-16") then
          INLET_RHOO16 = den_eos(1)*xn_eos(1,j)
       endif
    enddo
    INLET_TEMP = temp_eos(1)
    INLET_TRA = 0.0d0

    ! Now do the interior cells
    flameloc = ONE

    do i=0,n_base-1

       loloc = dble(i)*dx(dm) - flameloc
       hiloc = (dble(i) + ONE)*dx(dm) - flameloc

       call asin1d(lamsolfile, loloc, hiloc, state1d, ndum, .false.)

       p_eos(1) = Pamb
       den_eos(1) = state1d(3)
       temp_eos(1) = state1d(9)
       do j=1,nspec
          if(spec_names(j) .eq. "carbon-12") then
             xn_eos(1,j) = state1d(21)
          else if(spec_names(j) .eq. "magnesium-24") then
             xn_eos(1,j) = state1d(22)
          else if(spec_names(j) .eq. "oxygen-16") then
             xn_eos(1,j) = state1d(23)
          endif
       enddo

       ! given P, T, and X, compute rho
       call eos(eos_input_tp, den_eos, temp_eos, &
                npts, nspec, &
                xn_eos, &
                p_eos, h_eos, e_eos, & 
                cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                dpdX_eos, dhdX_eos, &
                gam1_eos, cs_eos, s_eos, &
                dsdt_eos, dsdr_eos, &
                do_diag)

       ! given rho, T, and X, compute h.
       call eos(eos_input_rt, den_eos, temp_eos, &
                npts, nspec, &
                xn_eos, &
                p_eos, h_eos, e_eos, & 
                cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                dpdX_eos, dhdX_eos, &
                gam1_eos, cs_eos, s_eos, &
                dsdt_eos, dsdr_eos, &
                do_diag)

       s0(i,rho_comp) = den_eos(1)
       if(use_big_h) then
          qreact = ZERO
          do j=1,nspec
             qreact = qreact + ebin(j)*xn_eos(1,j)
          enddo
          temporary = h_eos(1) + qreact
          s0(i,rhoh_comp) = den_eos(1)*temporary
       else
          s0(i,rhoh_comp) = den_eos(1)*h_eos(1)
       endif
       do j=1,nspec
          s0(i,spec_comp+j-1) = den_eos(1)*xn_eos(1,j)
       enddo
       s0(i,trac_comp) = 0.0d0
       s0(i,temp_comp) = temp_eos(1)
       p0(i) = pamb
       gam1(i) = gam1_eos(1)

    enddo

  end subroutine init_base_state


  subroutine write_base_state(state_name,w0_name,chk_name,s0,p0,w0,div_coeff)

    character(len=10), intent(in) :: state_name
    character(len=7), intent(in) :: w0_name
    character(len=7), intent(in) :: chk_name
    real(kind=dp_t) , intent(in) :: s0(:,:),p0(:),div_coeff(:), w0(:)
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
          write(99,1000)  base_r,s0(i,rho_comp), p0(i), s0(i,rhoh_comp), &
               (s0(i,n), n=spec_comp,temp_comp), div_coeff(i)
       end do
       close(99)

       ! write out w0 (it is nodal, so it gets a separate file)
       out_name = chk_name // "/" // w0_name
       write(6,*) 'Writing w0 state to ',out_name

       open(unit=99,file=out_name,form = "formatted", access = "sequential",action="write")
       do i = 1, nr+1
          base_r = (dble(i)-1) * dr
          write(99,1000)  base_r,w0(i)
       end do
       close(99)

    endif


1000 format(16(e30.20,1x))

  end subroutine write_base_state


  subroutine read_base_state(state_name,w0_name,chk_name,s0,p0,w0,div_coeff)
    
    character(len=10), intent(in   ) :: state_name
    character(len=7) , intent(in   ) :: w0_name
    character(len=7) , intent(in   ) :: chk_name    
    real(kind=dp_t) , intent(inout) :: s0(:,:),p0(:),div_coeff(:),w0(:)
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
       read(99,*)  base_r(i),s0(i,rho_comp), p0(i), s0(i,rhoh_comp), &
                   (s0(i,n), n=spec_comp,temp_comp), div_coeff(i)
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

end module base_state_module
