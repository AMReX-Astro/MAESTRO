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

          do n = 1,nscal
             call setbc_2d(sop(:,:,1,n), lo, ng, &
                           bc%adv_bc_level_array(i,:,:,dm+n),dx,dm+n)
          end do

       case (3)
          call initscalardata_3d(sop(:,:,:,:), lo, hi, ng, dx, perturb_model, &
                                 prob_lo, prob_hi, s0, p0)

          do n = 1, nscal
             call setbc_3d(sop(:,:,:,n), lo, ng, &
                           bc%adv_bc_level_array(i,:,:,dm+n),dx,dm+n)
          end do
       end select
    end do

    call multifab_fill_boundary(s)

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
    do n = rho_comp,temp_comp
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             s(i,j,n) = s0(j,n)
          enddo
       enddo
    enddo
    
    ! add an optional perturbation
    if (perturb_model) then
       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j)+HALF) * dx(2)
       
          do i = lo(1), hi(1)
             x = prob_lo(1) + (dble(i)+HALF) * dx(1)
          
             call perturb_2d(x, y, p0(j), s0(j,:), &
                             dens_pert, rhoh_pert, rhoX_pert, &
                             temp_pert, trac_pert)

             s(i,j,rho_comp) = dens_pert
             s(i,j,rhoh_comp) = rhoh_pert
             s(i,j,spec_comp:spec_comp+nspec-1) = rhoX_pert(1:)
             s(i,j,trac_comp:trac_comp+ntrac-1) = trac_pert(:)
          enddo
       enddo
    endif
    
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

       do n = rho_comp, spec_comp+nspec-1
          call fill_3d_data (s(:,:,:,n),s0(:,n),lo,hi,dx,ng)
       end do

    else 
 
       ! initialize the scalars
       do n = rho_comp,spec_comp+nspec-1
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
       
       if (perturb_model) then

          ! add an optional perturbation
          do k = lo(3), hi(3)
             z = prob_lo(3) + (dble(k)+HALF) * dx(3)
             
             do j = lo(2), hi(2)
                y = prob_lo(2) + (dble(j)+HALF) * dx(2)
                
                do i = lo(1), hi(1)
                   x = prob_lo(1) + (dble(i)+HALF) * dx(1)
                   
                   call perturb_3d(x, y, z, p0(k), s0(k,:), &
                                   dens_pert, rhoh_pert, rhoX_pert, &
                                   temp_pert, trac_pert)

                   s(i,j,k,rho_comp) = dens_pert
                   s(i,j,k,rhoh_comp) = rhoh_pert
                   s(i,j,k,spec_comp:spec_comp+nspec-1) = rhoX_pert(:)
                   s(i,j,k,trac_comp:trac_comp+ntrac-1) = trac_pert(:)
                enddo
             enddo
          enddo
       endif

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
   
          do n = 1,dm
             call setbc_2d(uop(:,:,1,n), lo, ng, &
                           bc%adv_bc_level_array(i,:,:,   n),dx,   n)
          end do

       case (3)
          call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx, &
                              prob_lo, prob_hi, s0, p0)

          do n = 1, dm
             call setbc_3d(uop(:,:,:,n), lo, ng, &
                           bc%adv_bc_level_array(i,:,:,   n),dx,   n)
          end do

       end select
    end do

    call multifab_fill_boundary(u)

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

    do i=0,hi(2)

       loloc = dble(i)*dx(dm) - flameloc
       hiloc = (dble(i) + ONE)*dx(dm) - flameloc

       call asin1d(lamsolfile, loloc, hiloc, state1d, ndum, .false.)

       u(lo(1):hi(1),i,1) = 0.0d0
       u(lo(1):hi(1),i,2) = state1d(2)

    enddo

    INLET_VN = 0.0d0
    INLET_VT = 0.0d0

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

    do i=0,hi(3)

       loloc = dble(i)*dx(dm)
       hiloc = (dble(i) + ONE)*dx(dm)

       call asin1d(lamsolfile, loloc, hiloc, state1d, ndum, .false.)

       u(lo(1):hi(1),lo(2):hi(2),i,1:2) = 0.0d0
       u(lo(1):hi(1),lo(2):hi(2),i,3) = state1d(2)

    enddo

    INLET_VN = 0.0d0
    INLET_VT = 0.0d0

  end subroutine initveldata_3d

  subroutine perturb_2d(x, y, p0, s0, dens_pert, rhoh_pert, rhoX_pert, temp_pert, trac_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y
    real(kind=dp_t), intent(in ) :: p0, s0(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp,t0
    real(kind=dp_t) :: x0, y0, x1, y1, x2, y2
    integer :: i, j
    real(kind=dp_t) :: r0, r1, r2

    t0 = s0(temp_comp)

    x0 = 5.0d7
    y0 = 6.5d7
    
    x1 = 1.2d8
    y1 = 8.5d7
    
    x2 = 2.0d8
    y2 = 7.5d7

    ! Tanh bubbles
    r0 = sqrt( (x-x0)**2 +(y-y0)**2 ) / 2.5e6
    r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / 0.05
    r2 = sqrt( (x-x2)**2 +(y-y2)**2 ) / 2.5e6
    
    ! This case works
    temp = t0 * (ONE + TWO * ( &
         .15_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r0))) + &
         .3_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r1))) + &
         .225_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r2)))  ) )

    ! This case breaks
!   temp = t0 * (ONE + FOUR * ( &
!        .15_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r0))) + &
!        .3_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r1))) + &
!        .225_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r2)))  ) )
          
    ! Use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    temp_row(1) = temp
    p_row(1) = p0
    den_row(1) = s0(rho_comp)
    xn_zone(:) = s0(spec_comp:)/s0(rho_comp)

    input_flag = 3      ! (t, p) -> (rho, h)

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

    dens_pert = den_row(1)
    rhoh_pert = den_row(1)*h_row(1)
    rhoX_pert(:) = dens_pert*xn_zone(:)

    temp_pert = temp
    
!   if ( (r0 .lt. 2.0) .or. (r1 .lt. 2.0) .or. (r2 .lt. 2.0) ) then
!     trac_pert(:) = ONE
!   else
      trac_pert(:) = ZERO
!   end if

  end subroutine perturb_2d

  subroutine perturb_3d(x, y, z, p0, s0, dens_pert, rhoh_pert, rhoX_pert, temp_pert, trac_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y, z
    real(kind=dp_t), intent(in ) :: p0, s0(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp,t0
    real(kind=dp_t) :: x0, y0, z0, x1, y1, z1, x2, y2, z2
    integer :: i, j, k
    real(kind=dp_t) :: r0, r1, r2

    x0 = 5.0d7
    y0 = 5.0d7
    z0 = 6.5d7
    
    x1 = 1.2d8
    y1 = 1.2d8
    z1 = 8.5d7
    
    x2 = 2.0d8
    y2 = 2.0d8
    z2 = 7.5d7

    t0 = s0(temp_comp)

!   temp = t0 * (ONE + TWO * ( &
!        .0625_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r0))) + &
!        .1875_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r1))) + &
!        .1250_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r2)))  ) )

    ! Tanh bubbles from perturb_2d
    r0 = sqrt( (y-y0)**2 +(z-z0)**2 ) / 2.5e6
    r1 = sqrt( (y-y1)**2 +(z-z1)**2 ) / 2.5e6
    r2 = sqrt( (y-y2)**2 +(z-z2)**2 ) / 2.5e6
    
    ! This case works
    temp = t0 * (ONE + TWO * ( &
         .150_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r0))) + &
         .300_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r1))) + &
         .225_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r2)))  ) )

    ! Use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    temp_row(1) = temp
    p_row(1) = p0
    den_row(1) = s0(rho_comp)
    xn_zone(:) = s0(spec_comp:)/s0(rho_comp)

    input_flag = 3      ! (t, p) -> (rho, h)

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

    dens_pert = den_row(1)
    rhoh_pert = den_row(1)*h_row(1)
    rhoX_pert(:) = dens_pert*xn_zone(:)

    temp_pert = temp
    
!   if (r1 .lt. 2.0) then
!     trac_pert(:) = ONE
!   else
      trac_pert(:) = ZERO
!   end if

  end subroutine perturb_3d

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
    parameter (ndum = 31)
    parameter (nspec = 3)

    character(len=128) :: lamsolfile
    real(kind=dp_t) :: state1d(ndum)
    real(kind=dp_t) :: loloc,hiloc,flameloc,qreact
    
    call helmeos_init

    dm = size(dx)

    lamsolfile = 'flame_4.e7_screen_left.out'

    ! ambient pressure from john's problem
    p_row(1) = 6.08741290776d24

    flameloc = ONE

    do i=0,n_base-1

       loloc = dble(i)*dx(dm) - flameloc
       hiloc = (dble(i) + ONE)*dx(dm) - flameloc

       call asin1d(lamsolfile, loloc, hiloc, state1d, ndum, .false.)

       temp_row(1) = state1d(9)

       xn_zone(1) = state1d(21)
       xn_zone(2) = state1d(23)
       xn_zone(3) = state1d(22)
       
       den_row(1) = state1d(3)
       
       ! given P, T, and X, compute rho and h.
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

       s0(i,rho_comp) = den_row(1)

       if(use_big_h) then
          qreact = ZERO
          do j=1,nspec
             qreact = qreact + ebin(j)*xn_zone(j)
          enddo
          s0(i,rhoh_comp) = den_row(1)*(h_row(1) + qreact)
       else
          s0(i,rhoh_comp) = den_row(1)*h_row(1)
       endif
       do j=1,nspec
          s0(i,spec_comp+j-1) = den_row(1)*xn_zone(j)
       enddo
       s0(i,trac_comp) = 0.0d0
       
       s0(i,temp_comp) = temp_row(1)

       p0(i) = p_row(1)
       gam1(i) = gam1_row(1)

    enddo

    ! temporary - when you fix make sure you do 3d case as well
    INLET_RHO = 3.9999999999798551d7
    INLET_RHOH = 3.00943325470894d25
    INLET_RHOC12 = 0.5d0*3.9999999999798551d7
    INLET_RHOO16 = 0.5d0*3.9999999999798551d7
    INLET_RHOMG24 = 0.0d0
    INLET_TEMP = 1.00000000047000d8
    INLET_TRA = 0.0d0
        
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
