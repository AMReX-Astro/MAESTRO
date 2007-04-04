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

  implicit none

contains

  subroutine initscalardata (s,s0,p0,temp0,dx,perturb_model, &
                             prob_lo,prob_hi,bc,nscal,ntrac)

    type(multifab) , intent(inout) :: s
    real(kind=dp_t), intent(in   ) ::    s0(:,:)
    real(kind=dp_t), intent(in   ) ::    p0(:)
    real(kind=dp_t), intent(in   ) :: temp0(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    logical,         intent(in   ) :: perturb_model
    real(kind=dp_t), intent(in   ) :: prob_lo(:)
    real(kind=dp_t), intent(in   ) :: prob_hi(:)
    type(bc_level) , intent(in   ) :: bc
    integer        , intent(in   ) :: nscal,ntrac

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
                                 prob_lo, prob_hi, s0, p0, temp0, ntrac)

          do n = 1,nscal
             call setbc_2d(sop(:,:,1,n), lo, ng, &
                           bc%adv_bc_level_array(i,:,:,dm+n),dx,dm+n)
          end do

       case (3)
          call initscalardata_3d(sop(:,:,:,:), lo, hi, ng, dx, perturb_model, &
                                 prob_lo, prob_hi, s0, p0, temp0, ntrac)

          do n = 1, nscal
             call setbc_3d(sop(:,:,:,n), lo, ng, &
                           bc%adv_bc_level_array(i,:,:,dm+n),dx,dm+n)
          end do
       end select
    end do

    call multifab_fill_boundary(s)

  end subroutine initscalardata

  subroutine initscalardata_2d (s,lo,hi,ng,dx, perturb_model, &
                                prob_lo,prob_hi,s0,p0,temp0,ntrac)

    integer, intent(in) :: lo(:), hi(:), ng, ntrac
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    logical,            intent(in ) :: perturb_model
    real (kind = dp_t), intent(in ) :: prob_lo(:)
    real (kind = dp_t), intent(in ) :: prob_hi(:)
    real(kind=dp_t), intent(in   ) ::    s0(0:,:)
    real(kind=dp_t), intent(in   ) ::    p0(0:)
    real(kind=dp_t), intent(in   ) :: temp0(0:)

    !     Local variables
    integer :: i, j, n
    real(kind=dp_t) :: x,y,r,r0,r1,r2,temp
    real(kind=dp_t) :: dens_pert, rhoh_pert, rhoX_pert(nspec), trac_pert(ntrac)


    ! initial the domain with the base state
    s = ZERO

    ! initialize the scalars
    do n = rho_comp,spec_comp+nspec-1
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
          
             call perturb_2d(x, y, temp0(j), p0(j), s0(j,:), &
                             dens_pert, rhoh_pert, rhoX_pert, trac_pert)

             s(i,j,rho_comp) = dens_pert
             s(i,j,rhoh_comp) = rhoh_pert
             s(i,j,spec_comp:spec_comp+nspec-1) = rhoX_pert(1:)
             s(i,j,trac_comp:trac_comp+ntrac-1) = trac_pert(:)
          enddo
       enddo
    endif
    
  end subroutine initscalardata_2d

  subroutine initscalardata_3d (s,lo,hi,ng,dx, perturb_model, &
                                prob_lo,prob_hi,s0,p0,temp0,ntrac)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng, ntrac
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    logical,            intent(in ) :: perturb_model
    real (kind = dp_t), intent(in ) :: prob_lo(:)
    real (kind = dp_t), intent(in ) :: prob_hi(:)
    real(kind=dp_t), intent(in   ) ::    s0(0:,:)
    real(kind=dp_t), intent(in   ) ::    p0(0:)
    real(kind=dp_t), intent(in   ) :: temp0(0:)

    !     Local variables
    integer :: i, j, k, n
    real(kind=dp_t) :: x,y,z,r,r0,r1,r2,temp
    real(kind=dp_t) :: dens_pert, rhoh_pert, rhoX_pert(nspec), trac_pert(ntrac)

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
       
       if (perturb_model) then

          ! add an optional perturbation
          do k = lo(3), hi(3)
             z = prob_lo(3) + (dble(k)+HALF) * dx(3)
             
             do j = lo(2), hi(2)
                y = prob_lo(2) + (dble(j)+HALF) * dx(2)
                
                do i = lo(1), hi(1)
                   x = prob_lo(1) + (dble(i)+HALF) * dx(1)
                   
                   call perturb_3d(x, y, z, temp0(k), p0(k), s0(k,:), &
                                   dens_pert, rhoh_pert, rhoX_pert, trac_pert)

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

  subroutine initveldata (u,s0,p0,temp0,dx,prob_lo,prob_hi,bc,nscal,ntrac)

    type(multifab) , intent(inout) :: u
    real(kind=dp_t), intent(in   ) ::    s0(:,:)
    real(kind=dp_t), intent(in   ) ::    p0(:)
    real(kind=dp_t), intent(in   ) :: temp0(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: prob_lo(:)
    real(kind=dp_t), intent(in   ) :: prob_hi(:)
    type(bc_level) , intent(in   ) :: bc
    integer        , intent(in   ) :: nscal,ntrac

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
                              prob_lo, prob_hi, s0, p0, temp0, ntrac)
   
          do n = 1,dm
             call setbc_2d(uop(:,:,1,n), lo, ng, &
                           bc%adv_bc_level_array(i,:,:,   n),dx,   n)
          end do

       case (3)
          call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx, &
                              prob_lo, prob_hi, s0, p0, temp0, ntrac)

          do n = 1, dm
             call setbc_3d(uop(:,:,:,n), lo, ng, &
                           bc%adv_bc_level_array(i,:,:,   n),dx,   n)
          end do

       end select
    end do

    call multifab_fill_boundary(u)

  end subroutine initveldata

  subroutine initveldata_2d (u,lo,hi,ng,dx, &
                             prob_lo,prob_hi,s0,p0,temp0,ntrac)

    integer, intent(in) :: lo(:), hi(:), ng, ntrac
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real (kind = dp_t), intent(in ) :: prob_lo(:)
    real (kind = dp_t), intent(in ) :: prob_hi(:)
    real(kind=dp_t), intent(in   ) ::    s0(0:,:)
    real(kind=dp_t), intent(in   ) ::    p0(0:)
    real(kind=dp_t), intent(in   ) :: temp0(0:)

    !     Local variables
    integer :: i, j, n
    real(kind=dp_t) :: x,y,r,r0,r1,r2,temp
    real(kind=dp_t) :: dens_pert, rhoh_pert, rhoX_pert(nspec), trac_pert(ntrac)

    ! initial the velocity
    u = ZERO

  end subroutine initveldata_2d

  subroutine initveldata_3d (u,lo,hi,ng,dx, &
                             prob_lo,prob_hi,s0,p0,temp0,ntrac)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng, ntrac
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real (kind = dp_t), intent(in ) :: prob_lo(:)
    real (kind = dp_t), intent(in ) :: prob_hi(:)
    real(kind=dp_t), intent(in   ) ::    s0(0:,:)
    real(kind=dp_t), intent(in   ) ::    p0(0:)
    real(kind=dp_t), intent(in   ) :: temp0(0:)

    !     Local variables
    integer :: i, j, k, n
    real(kind=dp_t) :: x,y,z,r,r0,r1,r2,temp
    real(kind=dp_t) :: dens_pert, rhoh_pert, rhoX_pert(nspec), trac_pert(ntrac)

    ! initial the velocity
    u = ZERO
    
  end subroutine initveldata_3d


  subroutine perturb_2d(x, y, t0, p0, s0, dens_pert, rhoh_pert, rhoX_pert, trac_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y
    real(kind=dp_t), intent(in ) :: t0, p0, s0(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp
    real(kind=dp_t) :: x0, y0, x1, y1, x2, y2
    integer :: i, j
    real(kind=dp_t) :: r0, r1, r2


    x0 = 5.0d7
    y0 = 6.5d7
    
    x1 = 1.2d8
    y1 = 8.5d7
    
    x2 = 2.0d8
    y2 = 7.5d7

    ! Tanh bubbles
    r0 = sqrt( (x-x0)**2 +(y-y0)**2 ) / 2.5e6
    r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / 2.5e6
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
    
    if ( (r0 .lt. 2.0) .or. (r1 .lt. 2.0) .or. (r2 .lt. 2.0) ) then
      trac_pert(:) = ONE
    else
      trac_pert(:) = ZERO
    end if

  end subroutine perturb_2d

  subroutine perturb_3d(x, y, z, t0, p0, s0, dens_pert, rhoh_pert, rhoX_pert, trac_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y, z
    real(kind=dp_t), intent(in ) :: t0, p0, s0(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp
    real(kind=dp_t) :: x0, y0, z0, x1, y1, z1, x2, y2, z2
    integer :: i, j, k
    real(kind=dp_t) :: r0, r1, r2

    x0 = 5.0d7
    z0 = 6.5d7
    
    x1 = 1.2d8
    z1 = 8.5d7
    
    x2 = 2.0d8
    z2 = 7.5d7

!   temp = t0 * (ONE + TWO * ( &
!        .0625_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r0))) + &
!        .1875_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r1))) + &
!        .1250_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r2)))  ) )

    ! Tanh bubbles from perturb_2d
    r0 = sqrt( (x-x0)**2 +(z-z0)**2 ) / 2.5e6
    r1 = sqrt( (x-x1)**2 +(z-z1)**2 ) / 2.5e6
    r2 = sqrt( (x-x2)**2 +(z-z2)**2 ) / 2.5e6
    
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
    
    if (r1 .lt. 2.0) then
      trac_pert(:) = ONE
    else
      trac_pert(:) = ZERO
    end if

  end subroutine perturb_3d

  subroutine init_base_state (n_base,s0,temp0,p0,gam1,dx,prob_lo,prob_hi)

    integer        , intent(in   ) :: n_base
    real(kind=dp_t), intent(inout) ::    s0(0:,:)
    real(kind=dp_t), intent(inout) :: temp0(0:)
    real(kind=dp_t), intent(inout) ::    p0(0:)
    real(kind=dp_t), intent(inout) ::  gam1(0:)
    real(kind=dp_t), intent(in   ) :: prob_lo(:)
    real(kind=dp_t), intent(in   ) :: prob_hi(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer :: i,j,n,j_cutoff

    real(kind=dp_t) :: r,dr_in
    real(kind=dp_t) :: d_ambient,t_ambient,p_ambient, xn_ambient(nspec)
    real(kind=dp_t) :: integral, temp_term_lo, temp_term_hi
    real(kind=dp_t) :: temp_min,p0_lo,p0_hi

    ! these indices define how the initial model is stored in the 
    ! base_state array
    integer, parameter :: nvars_model = 3 + nspec
    integer, parameter :: idens_model = 1
    integer, parameter :: itemp_model = 2
    integer, parameter :: ipres_model = 3
    integer, parameter :: ispec_model = 4

    integer, parameter :: MAX_VARNAME_LENGTH=80
    integer :: npts_model, nvars_model_file

    real(kind=dp_t), allocatable :: base_state(:,:), base_r(:)
    real(kind=dp_t), allocatable :: vars_stored(:)
    character(len=MAX_VARNAME_LENGTH), allocatable :: varnames_stored(:)
    logical :: found

    integer :: ipos,dm
    character (len=256) :: header_line

    logical :: do_diag

    real(kind=dp_t), parameter :: cutoff_density = 3.0d6

    do_diag = .false.

    ! open the model file and read in the header
    ! the model file is assumed to be of the follow form:
    ! # npts = 896
    ! # num of variables = 6
    ! # density
    ! # temperature
    ! # pressure
    ! # carbon-12
    ! # oxygen-16
    ! # magnesium-24
    ! 195312.5000  5437711139.  8805500.952   .4695704813E+28  0.3  0.7  0
    ! 585937.5000  5410152416.  8816689.836  0.4663923963E+28  0.3  0.7  0

    ! we read in the number of variables and their order and use this to map 
    ! them into the base_state array.  We ignore anything other than density, 
    ! temperature, pressure and composition.  

    ! Presently, we take density, temperature, and composition as the 
    ! independent variables and use them to define the thermodynamic state.

    ! composition is assumed to be in terms of mass fractions
    
    open(99,file="model.hse")

    ! the first line has the number of points in the model
    read (99, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) npts_model

    if ( parallel_IOProcessor() ) &
      print *, npts_model, '    points found in the initial model file'

    ! now read in the number of variables
    read (99, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) nvars_model_file

    if ( parallel_IOProcessor() ) &
      print *, nvars_model_file, ' variables found in the initial model file'

    allocate (vars_stored(nvars_model_file))
    allocate (varnames_stored(nvars_model_file))

    ! now read in the names of the variables
    do i = 1, nvars_model_file
       read (99, '(a256)') header_line
       ipos = index(header_line, '#') + 1
       varnames_stored(i) = trim(adjustl(header_line(ipos:)))
    enddo
    
    ! allocate storage for the model data
    allocate (base_state(npts_model, nvars_model))
    allocate (base_r(npts_model))

    do i = 1, npts_model
       read(99,*) base_r(i), (vars_stored(j), j = 1, nvars_model_file)

       base_state(i,:) = ZERO

       do j = 1, nvars_model_file

          found = .false.
       
          if (trim(varnames_stored(j)) == "density") then
             base_state(i,idens_model) = vars_stored(j)
             found = .true.

          else if (trim(varnames_stored(j)) == "temperature") then
             base_state(i,itemp_model) = vars_stored(j)
             found = .true.

          else if (trim(varnames_stored(j)) == "pressure") then
             base_state(i,ipres_model) = vars_stored(j)
             found = .true.

          else
             do n = 1, nspec
                if (trim(varnames_stored(j)) == spec_names(n)) then
                   base_state(i,ispec_model-1+n) = vars_stored(j)
                   found = .true.
                   exit
                endif
             enddo
          endif

          if (.NOT. found) then
             print *, 'ERROR: variable not found: ', varnames_stored(j)
          endif
          
       enddo

    end do

    close(99)

    call helmeos_init

    dr_in = base_r(2) - base_r(1)
    dr = dr_in * dble(npts_model)  / dble(n_base)
    if ( parallel_IOProcessor() ) then
      print *,'DR_IN  ',dr_in
      print *,'NEW DR ',dr
    end if

    if (spherical .eq. 1) then
      if (dble(npts_model)*dr_in .lt. HALF*prob_hi(1)*dsqrt(3.0_dp_t)) then
         print *,'OOPS - NEED LARGER RADIAL ARRAY OR SMALLER DOMAIN '
         print *,'-- base state goes to radius ',dble(npts_model)*dr_in
         print *,'-- problem domain goes from ',prob_lo(1),' to ',prob_hi(1)
         print *,'-- which means corners are at a radius of ', &
                  HALF*prob_hi(1)*dsqrt(3.0_dp_t)
         stop
      end if
    end if

    j_cutoff = n_base-1
    do j = 0,n_base-1

       ! compute the coordinate height at this level
       ! NOTE: we are assuming that the basestate is in the y-direction
       ! and that ymin = 0.0
       r = prob_lo(2) + (dble(j) + HALF)*dr

       d_ambient = interpolate(r, npts_model, base_r, base_state(:,idens_model))
       t_ambient = interpolate(r, npts_model, base_r, base_state(:,itemp_model))
       p_ambient = interpolate(r, npts_model, base_r, base_state(:,ipres_model))

       do n = 1, nspec
          xn_ambient(n) = interpolate(r, npts_model, base_r, base_state(:,ispec_model-1+n))
       enddo

       ! use the EOS to make the state consistent
       temp_row(1) = t_ambient
       den_row(1)  = d_ambient
       p_row(1)    = p_ambient

       ! (rho,T) --> p,h
       input_flag = 1

       call eos(input_flag, den_row, temp_row, &
                npts, nspec, &
                xn_ambient, aion, zion, &
                p_row, h_row, e_row, &
                cv_row, cp_row, xne_row, eta_row, pele_row, &
                dpdt_row, dpdr_row, dedt_row, dedr_row, &
                dpdX_row, dhdX_row, &
                gam1_row, cs_row, s_row, &
                dsdt_row, dsdr_row, &
                do_diag)
       
       s0(j, rho_comp ) = d_ambient
       s0(j,rhoh_comp ) = d_ambient * h_row(1)
       s0(j,spec_comp:spec_comp+nspec-1) = d_ambient * xn_ambient(1:nspec)
       p0(j)    = p_row(1)
       temp0(j) = t_ambient

       ! keep track of the height where we drop below the cutoff density
       if (s0(j,rho_comp) .lt. cutoff_density .and. j_cutoff .eq. n_base-1) j_cutoff = j

    end do
 

    ! KEEP RHO0 CONSTANT ABOVE J_CUTOFF
    if ( parallel_IOProcessor() ) &
      print *,'JCUT_comp ',j_cutoff

    do j = j_cutoff,n_base-1
       s0(j, rho_comp ) = s0(j_cutoff, rho_comp )
       s0(j,rhoh_comp ) = s0(j_cutoff,rhoh_comp )
       s0(j,spec_comp:spec_comp+nspec-1) = s0(j_cutoff,spec_comp:spec_comp+nspec-1)
       p0(j)            = p0(j_cutoff)
       temp0(j)         = temp0(j_cutoff)
    end do

      ! RECALCULATE T, RHO_H

    do j = 0,n_base-1

       den_row(1)  = s0(j,rho_comp)
       xn_ambient(1:nspec) = s0(j,spec_comp:spec_comp-1+nspec)/den_row(1)
       temp_row(1) = temp0(j)
       p_row(1)    =   p0(j)

!       (rho, p) --> T, h
       input_flag = 4

       call eos(input_flag, den_row, temp_row, &
                npts, nspec, &
                xn_ambient, aion, zion, &
                p_row, h_row, e_row, &
                cv_row, cp_row, xne_row, eta_row, pele_row, &
                dpdt_row, dpdr_row, dedt_row, dedr_row, &
                dpdX_row, dhdX_row, &
                gam1_row, cs_row, s_row, &
                dsdt_row, dsdr_row, &
                do_diag)
       
       temp0(j) = temp_row(1)
       gam1(j)  = gam1_row(1)
       s0(j,rhoh_comp) = h_row(1) * s0(j,rho_comp)

       s0(j,trac_comp:) = ZERO
       
    end do

  end subroutine init_base_state

  function interpolate(r, npts, model_r, model_var)

    ! given the array of model coordinates (model_r), and variable (model_var),
    ! find the value of model_var at point r using linear interpolation.
    ! Eventually, we can do something fancier here.

    real(kind=dp_t) :: interpolate
    real(kind=dp_t), intent(in) :: r
    integer :: npts
    real(kind=dp_t), dimension(npts) :: model_r, model_var

    real(kind=dp_t) :: val, slope

    integer :: i, id

    ! find the location in the coordinate array where we want to interpolate
    do i = 1, npts
       if (model_r(i) >= r) exit
    enddo

    id = i

    if (id == 1) then
       slope = (model_var(id+1) - model_var(id))/(model_r(id+1) - model_r(id))
    else if (id == npts) then
       slope = (model_var(id) - model_var(id-1))/(model_r(id) - model_r(id-1))
    else
       slope = (model_var(id+1) - model_var(id-1))/(model_r(id+1) - model_r(id-1))
    endif

    interpolate = slope*(r - model_r(id)) + model_var(id)

    return

  end function interpolate


  subroutine write_base_state(sd_name,s0,temp0,p0,div_coeff)

    character(len=9), intent(in) :: sd_name
    real(kind=dp_t) , intent(in) :: s0(:,:),p0(:),temp0(:),div_coeff(:)

    integer :: i, n

    if (parallel_IOProcessor()) &
       write(6,*) 'Writing base state to ',sd_name

    open(unit=99,file=sd_name,form = "formatted", access = "sequential",action="write")
    do i = 1, size(s0,dim=1)
       write(99,1000)  s0(i,rho_comp), p0(i), s0(i,rhoh_comp), &
                       (s0(i,n), n=spec_comp,spec_comp+nspec-1), temp0(i), div_coeff(i)
    end do
    close(99)

1000 format(16(e18.12,1x))

  end subroutine write_base_state


  subroutine read_base_state(sd_name,s0,temp0,p0,div_coeff)
    
    character(len=9), intent(in   ) :: sd_name
    real(kind=dp_t) , intent(inout) :: s0(:,:),p0(:),temp0(:),div_coeff(:)

    integer :: i, n

    print *,'Reading base state from ',sd_name

    open(unit=99,file=sd_name)
    do i = 1, size(s0,dim=1)
       read(99,*)  s0(i,rho_comp), p0(i), s0(i,rhoh_comp), &
                   (s0(i,n), n=spec_comp,spec_comp+nspec-1), temp0(i), div_coeff(i)
    end do
    close(99)

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
