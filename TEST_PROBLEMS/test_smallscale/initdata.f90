module init_module

  use bl_types
  use bl_constants_module
  use bc_module
  use setbc_module
  use define_bc_module
  use multifab_module
  use fill_3d_module
  use eos_module
  use variables
  use network
  use geometry
  use probin_module
  use inlet_bc_module

  implicit none

contains

  subroutine initscalardata (s,s0,p0,dx,perturb_model, &
                             prob_lo,prob_hi,bc)

    type(multifab) , intent(inout) :: s
    real(kind=dp_t), intent(inout) ::    s0(0:,:)
    real(kind=dp_t), intent(in   ) ::    p0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    logical,         intent(in   ) :: perturb_model
    real(kind=dp_t), intent(in   ) :: prob_lo(:)
    real(kind=dp_t), intent(in   ) :: prob_hi(:)
    type(bc_level) , intent(in   ) :: bc

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng,dm
    integer :: i,n
    
    ng = s%ng
    dm = s%dim

    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       sop => dataptr(s, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (2)
          call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx, perturb_model, &
                                 prob_lo, prob_hi, s0, p0)
       case (3)
          call initscalardata_3d(sop(:,:,:,:), lo, hi, ng, dx, perturb_model, &
                                 prob_lo, prob_hi, s0, p0)
       end select
    end do

    call multifab_fill_boundary(s)

    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       sop => dataptr(s, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (2)
          do n = 1,nscal
             call setbc_2d(sop(:,:,1,n), lo, ng, &
                           bc%adv_bc_level_array(i,:,:,dm+n),dx,dm+n)
          end do

       case (3)
          do n = 1,nscal
             call setbc_3d(sop(:,:,:,n), lo, ng, &
                           bc%adv_bc_level_array(i,:,:,dm+n),dx,dm+n)
          end do
       end select
    end do

  end subroutine initscalardata

  subroutine initscalardata_2d (s,lo,hi,ng,dx, perturb_model, &
                                prob_lo,prob_hi,s0,p0)

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    logical,            intent(in ) :: perturb_model
    real (kind = dp_t), intent(in ) :: prob_lo(:)
    real (kind = dp_t), intent(in ) :: prob_hi(:)
    real(kind=dp_t), intent(inout) ::    s0(0:,:)
    real(kind=dp_t), intent(in   ) ::    p0(0:)

    !     Local variables
    integer :: i, j, n
    real(kind=dp_t) :: x,y,r,r0,r1,r2,temp
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)

    ! initial the domain with the base state
    s = ZERO

    ! initialize the scalars
    do n = rho_comp,rho_comp+nscal-1
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             s(i,j,n) = s0(j,n)
          enddo
       enddo
    enddo
    
    ! set density in base state to the constant, "bottom domain" value
    do j=lo(2),hi(2)
       s0(j,rho_comp) = s0(j,rho_comp)
    enddo
    
  end subroutine initscalardata_2d

  subroutine initscalardata_3d (s,lo,hi,ng,dx, perturb_model, &
                                prob_lo,prob_hi,s0,p0)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    logical,            intent(in ) :: perturb_model
    real (kind = dp_t), intent(in ) :: prob_lo(:)
    real (kind = dp_t), intent(in ) :: prob_hi(:)
    real(kind=dp_t), intent(inout) ::    s0(0:,:)
    real(kind=dp_t), intent(in   ) ::    p0(0:)

    !     Local variables
    integer :: i, j, k, n
    real(kind=dp_t) :: x,y,z,r,r0,r1,r2,temp
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)

    ! initial the domain with the base state
    s = ZERO
  
    if (spherical .eq. 1) then
       do n = rho_comp,rho_comp+nscal-1
          call fill_3d_data(s(:,:,:,n),s0(:,n),lo,hi,dx,ng)
       end do
    else 
        ! initialize the scalars
       do n = rho_comp,rho_comp+nscal-1
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   s(i,j,k,n) = s0(k,n)
                enddo
             enddo
          enddo
       enddo
       
       ! set density in base state to the constant, "bottom domain" value
       do k=lo(3),hi(3)
          s0(k,rho_comp) = s0(k,rho_comp)
       enddo
    end if
    
  end subroutine initscalardata_3d

  subroutine initveldata (u,s0,p0,dx,prob_lo,prob_hi,bc)

    type(multifab) , intent(inout) :: u
    real(kind=dp_t), intent(in   ) ::    s0(:,:)
    real(kind=dp_t), intent(in   ) ::    p0(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: prob_lo(:)
    real(kind=dp_t), intent(in   ) :: prob_hi(:)
    type(bc_level) , intent(in   ) :: bc

    real(kind=dp_t), pointer:: uop(:,:,:,:)
    integer :: lo(u%dim),hi(u%dim),ng,dm
    integer :: i,n
    
    ng = u%ng
    dm = u%dim

    do i = 1, u%nboxes
       if ( multifab_remote(u, i) ) cycle
       uop => dataptr(u, i)
       lo =  lwb(get_box(u, i))
       hi =  upb(get_box(u, i))
       select case (dm)
       case (2)
          call initveldata_2d(uop(:,:,1,:), lo, hi, ng, dx, &
                              prob_lo, prob_hi, s0, p0)
       case (3)
          call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx, &
                              prob_lo, prob_hi, s0, p0)
       end select
    end do

    call multifab_fill_boundary(u)

    do i = 1, u%nboxes
       if ( multifab_remote(u, i) ) cycle
       uop => dataptr(u, i)
       lo =  lwb(get_box(u, i))
       hi =  upb(get_box(u, i))
       select case (dm)
       case (2)
          do n = 1,dm
             call setbc_2d(uop(:,:,1,n), lo, ng, &
                           bc%adv_bc_level_array(i,:,:,   n),dx,   n)
          end do
       case (3)
          do n = 1, dm
             call setbc_3d(uop(:,:,:,n), lo, ng, &
                           bc%adv_bc_level_array(i,:,:,   n),dx,   n)
          end do
       end select
    end do

  end subroutine initveldata

  subroutine initveldata_2d (u,lo,hi,ng,dx, &
                             prob_lo,prob_hi,s0,p0)

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real (kind = dp_t), intent(in ) :: prob_lo(:)
    real (kind = dp_t), intent(in ) :: prob_hi(:)
    real(kind=dp_t), intent(in   ) ::    s0(0:,:)
    real(kind=dp_t), intent(in   ) ::    p0(0:)

    ! local
    integer ndum, i, dm
    parameter (ndum = 31)

    character(len=128) :: lamsolfile
    real(kind=dp_t) :: state1d(ndum)
    real(kind=dp_t) :: loloc,hiloc,flameloc
    
    dm = size(dx)

    lamsolfile = 'flame_4.e7_screen_left.out'

    flameloc = ONE

    do i=lo(2),hi(2)

       loloc = dble(i)*dx(dm) - flameloc
       hiloc = (dble(i) + ONE)*dx(dm) - flameloc

       call asin1d(lamsolfile, loloc, hiloc, state1d, ndum, .false.)

       u(lo(1):hi(1),i,1) = 0.0d0
       u(lo(1):hi(1),i,2) = state1d(2)

    enddo

  end subroutine initveldata_2d

  subroutine initveldata_3d (u,lo,hi,ng,dx, &
                             prob_lo,prob_hi,s0,p0)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real (kind = dp_t), intent(in ) :: prob_lo(:)
    real (kind = dp_t), intent(in ) :: prob_hi(:)
    real(kind=dp_t), intent(in   ) ::    s0(0:,:)
    real(kind=dp_t), intent(in   ) ::    p0(0:)

    ! local
    integer ndum, i, dm
    parameter (ndum = 31)

    character(len=128) :: lamsolfile
    real(kind=dp_t) :: state1d(ndum)
    real(kind=dp_t) :: loloc,hiloc
    
    dm = size(dx)

    lamsolfile = 'flame_4.e7_screen_left.out'

    do i=lo(3),hi(3)

       loloc = dble(i)*dx(dm)
       hiloc = (dble(i) + ONE)*dx(dm)

       call asin1d(lamsolfile, loloc, hiloc, state1d, ndum, .false.)

       u(lo(1):hi(1),lo(2):hi(2),i,1:2) = 0.0d0
       u(lo(1):hi(1),lo(2):hi(2),i,3) = state1d(2)

    enddo

  end subroutine initveldata_3d

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
    integer input_flag
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
    p_row(1) = Pamb

    den_row(1) = state1d(3)
    temp_row(1) = state1d(9)
    do j=1,nspec
       if(spec_names(j) .eq. "carbon-12") then
          xn_zone(j) = state1d(21)
       else if(spec_names(j) .eq. "magnesium-24") then
          xn_zone(j) = state1d(22)
       else if(spec_names(j) .eq. "oxygen-16") then
          xn_zone(j) = state1d(23)
       else
          print*,"In initdata, spec_names(",j,") invalid"
       endif
    enddo

    ! given P, T, and X, compute rho
    input_flag = 3
    call eos(input_flag, den_row, temp_row, &
         npts, nspec, &
         xn_zone, aion, zion, &
         p_row, h_row, e_row, & 
         cv_row, cp_row, xne_row, eta_row, pele_row, &
         dpdt_row, dpdr_row, dedt_row, dedr_row, &
         dpdX_row, dhdX_row, &
         gam1_row, cs_row, s_row, &
         dsdt_row, dsdr_row, &
         do_diag)

    ! given rho, T, and X, compute h
    input_flag = 1
    call eos(input_flag, den_row, temp_row, &
         npts, nspec, &
         xn_zone, aion, zion, &
         p_row, h_row, e_row, & 
         cv_row, cp_row, xne_row, eta_row, pele_row, &
         dpdt_row, dpdr_row, dedt_row, dedr_row, &
         dpdX_row, dhdX_row, &
         gam1_row, cs_row, s_row, &
         dsdt_row, dsdr_row, &
         do_diag)

    INLET_VN = 0.0d0
    INLET_VT = 0.0d0
    INLET_RHO = den_row(1)
    if(use_big_h) then
       qreact = 0.0d0
       do j=1,nspec
          qreact = qreact + ebin(j)*xn_zone(j)
       enddo
       INLET_RHOH = den_row(1)*(h_row(1) + qreact)
    else
       INLET_RHOH = den_row(1)*h_row(1)
    endif
    do j=1,nspec
       if(spec_names(j) .eq. "carbon-12") then
          INLET_RHOC12 = den_row(1)*xn_zone(j)
       else if(spec_names(j) .eq. "magnesium-24") then
          INLET_RHOMG24 = den_row(1)*xn_zone(j)
       else if(spec_names(j) .eq. "oxygen-16") then
          INLET_RHOO16 = den_row(1)*xn_zone(j)
       endif
    enddo
    INLET_TEMP = temp_row(1)
    INLET_TRA = 0.0d0

    ! Now do the interior cells
    flameloc = ONE

    do i=0,n_base-1

       loloc = dble(i)*dx(dm) - flameloc
       hiloc = (dble(i) + ONE)*dx(dm) - flameloc

       call asin1d(lamsolfile, loloc, hiloc, state1d, ndum, .false.)

       p_row(1) = Pamb
       den_row(1) = state1d(3)
       temp_row(1) = state1d(9)
       do j=1,nspec
          if(spec_names(j) .eq. "carbon-12") then
             xn_zone(j) = state1d(21)
          else if(spec_names(j) .eq. "magnesium-24") then
             xn_zone(j) = state1d(22)
          else if(spec_names(j) .eq. "oxygen-16") then
             xn_zone(j) = state1d(23)
          endif
       enddo

       ! given P, T, and X, compute rho
       input_flag = 3
       call eos(input_flag, den_row, temp_row, &
                npts, nspec, &
                xn_zone, aion, zion, &
                p_row, h_row, e_row, & 
                cv_row, cp_row, xne_row, eta_row, pele_row, &
                dpdt_row, dpdr_row, dedt_row, dedr_row, &
                dpdX_row, dhdX_row, &
                gam1_row, cs_row, s_row, &
                dsdt_row, dsdr_row, &
                do_diag)

       ! given rho, T, and X, compute h.
       input_flag = 1
       call eos(input_flag, den_row, temp_row, &
                npts, nspec, &
                xn_zone, aion, zion, &
                p_row, h_row, e_row, & 
                cv_row, cp_row, xne_row, eta_row, pele_row, &
                dpdt_row, dpdr_row, dedt_row, dedr_row, &
                dpdX_row, dhdX_row, &
                gam1_row, cs_row, s_row, &
                dsdt_row, dsdr_row, &
                do_diag)

       s0(i,rho_comp) = den_row(1)
       if(use_big_h) then
          qreact = ZERO
          do j=1,nspec
             qreact = qreact + ebin(j)*xn_zone(j)
          enddo
          temporary = h_row(1) + qreact
          s0(i,rhoh_comp) = den_row(1)*temporary
       else
          s0(i,rhoh_comp) = den_row(1)*h_row(1)
       endif
       do j=1,nspec
          s0(i,spec_comp+j-1) = den_row(1)*xn_zone(j)
       enddo
       s0(i,trac_comp) = 0.0d0
       s0(i,temp_comp) = temp_row(1)
       p0(i) = pamb
       gam1(i) = gam1_row(1)

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

  subroutine scalar_diags (istep,s,s0,dx)

    integer        , intent(in   ) :: istep
    type(multifab) , intent(inout) :: s
    real(kind=dp_t), intent(in)    :: s0(:,:)
    real(kind=dp_t), intent(in)    :: dx(:)

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng,dm
    integer :: i,n
    
    ng = s%ng
    dm = s%dim

    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       sop => dataptr(s, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))

       select case (dm)
       case (2)
          call scalar_diags_2d(istep, sop(:,:,1,:), lo, hi, ng, dx, s0)
       case (3)
!         call scalar_diags_3d(istep, sop(:,:,:,:), lo, hi, ng, dx, s0)
       end select
    end do

  end subroutine scalar_diags

  subroutine scalar_diags_2d (istep, s,lo,hi,ng,dx,s0)

    integer, intent(in) :: istep, lo(:), hi(:), ng
    real (kind = dp_t), intent(in) ::  s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in) :: dx(:)
    real(kind=dp_t)   , intent(in) :: s0(0:,:)

    ! Local variables
    integer :: i, j, n
    real(kind=dp_t) :: fac, stot, smax
    character(len=11) :: file_name

    write(unit=file_name,fmt='("rhodiag",i4.4)') istep
    open(90,file=file_name)

    fac = ONE / dble(hi(1)-lo(1)+1)
    do j = lo(2), hi(2)
      stot = ZERO
      smax = ZERO
      do i = lo(1), hi(1)
         stot = stot + (s(i,j,rho_comp) - s0(j,rho_comp))
         smax = max(smax,abs(s(i,j,rho_comp) - s0(j,rho_comp)))
      enddo
      write(90,*) j,stot*fac/ s0(j,rho_comp), smax / s0(j,rho_comp)
    enddo
    
  end subroutine scalar_diags_2d

end module init_module
