module init_module

  use bl_types
  use bl_constants_module
  use bc_module
  use define_bc_module
  use multifab_module
  use fill_3d_module
  use eos_module
  use variables
  use network
  use geometry
  use ml_layout_module
  use ml_restriction_module
  use multifab_fill_ghost_module

  implicit none
  public :: initscalardata, initveldata

contains

  subroutine initscalardata(nlevs,s,s0,p0,dx,bc,mla)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: s0(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(s(1)%dim),hi(s(1)%dim),ng,dm
    integer :: i,n
    
    ng = s(1)%ng
    dm = s(1)%dim

    do n=1,nlevs

       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          sop => dataptr(s(n), i)
          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))
          select case (dm)
          case (2)
             call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx(n,:), s0(n,:,:), p0(n,:))
          case (3)

          end select
       end do

    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(s(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(s(nlevs),rho_comp,dm+rho_comp,nscal,bc(nlevs))

    else

       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(s(n-1),s(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(s(n),s(n-1),ng,mla%mba%rr(n-1,:), &
                                         bc(n-1),bc(n),1,dm+rho_comp,nscal)

       enddo

    end if

  end subroutine initscalardata

  subroutine initscalardata_2d(s,lo,hi,ng,dx,s0,p0)

    use probin_module, only: prob_lo_x, prob_lo_y, perturb_model

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0(0:)

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
          y = prob_lo_y + (dble(j)+HALF) * dx(2)
       
          do i = lo(1), hi(1)
             x = prob_lo_x + (dble(i)+HALF) * dx(1)
          
             call perturb_2d(x, y, s0(j,temp_comp), p0(j), s0(j,:), &
                             dens_pert, rhoh_pert, rhoX_pert, trac_pert)

             s(i,j,rho_comp) = dens_pert
             s(i,j,rhoh_comp) = rhoh_pert
             s(i,j,spec_comp:spec_comp+nspec-1) = rhoX_pert(1:)
             s(i,j,trac_comp:trac_comp+ntrac-1) = trac_pert(:)
          enddo
       enddo
    endif
    
  end subroutine initscalardata_2d

  subroutine initveldata(nlevs,u,s0,p0,dx,bc,mla)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: u(:)
    real(kind=dp_t), intent(in   ) :: s0(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: uop(:,:,:,:)
    integer :: lo(u(1)%dim),hi(u(1)%dim),ng,dm
    integer :: i,n
    
    ng = u(1)%ng
    dm = u(1)%dim

    do n=1,nlevs

       do i = 1, u(n)%nboxes
          if ( multifab_remote(u(n), i) ) cycle
          uop => dataptr(u(n), i)
          lo =  lwb(get_box(u(n), i))
          hi =  upb(get_box(u(n), i))
          select case (dm)
          case (2)
             call initveldata_2d(uop(:,:,1,:), lo, hi, ng, dx(n,:), s0(n,:,:), p0(n,:))
          case (3)
          end select
       end do
       
    end do

    if (nlevs .eq. 1) then

       ! fill ghost cells for two adjacent grids at the same level
       ! this includes periodic domain boundary ghost cells
       call multifab_fill_boundary(u(nlevs))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(u(nlevs),1,1,dm,bc(nlevs))
    else
    
       ! the loop over nlevs must count backwards to make sure the finer grids are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering it
          call ml_cc_restriction(u(n-1),u(n),mla%mba%rr(n-1,:))

          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(u(n),u(n-1),ng,mla%mba%rr(n-1,:), &
                                         bc(n-1),bc(n),1,1,dm)
       enddo
       
    end if

  end subroutine initveldata

  subroutine initveldata_2d(u,lo,hi,ng,dx,s0,p0)

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0(0:)

    ! initial the velocity
    u(:,:,1) = ZERO
    u(:,:,2) = ZERO

  end subroutine initveldata_2d


  subroutine perturb_2d(x, y, t0, p0, s0, dens_pert, rhoh_pert, rhoX_pert, trac_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y
    real(kind=dp_t), intent(in ) :: t0, p0, s0(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp
    real(kind=dp_t) :: x1, y1, r1
    integer :: i, j

    
!   x1 = 1.0d8
    x1 = 0.72d8
    y1 = 8.5d7

    ! Tanh bubbles
    r1 = sqrt( (x-x1)**2 +(y-y1)**2 ) / 2.5e6
    
    ! One bubble for convergence test
    temp = t0 * (ONE + TWO * ( &
         .300_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r1))) ) )
          
    ! Use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    temp_eos(1) = temp
    p_eos(1) = p0
    den_eos(1) = s0(rho_comp)
    xn_eos(1,:) = s0(spec_comp:)/s0(rho_comp)

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

    dens_pert = den_eos(1)
    rhoh_pert = den_eos(1)*h_eos(1)
    rhoX_pert(:) = dens_pert*xn_eos(1,:)
    
    if ( (r1 .lt. 2.0) ) then
      trac_pert(:) = ONE
    else
      trac_pert(:) = ZERO
    end if

  end subroutine perturb_2d

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

    integer :: iratio, offset

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

    print *, npts_model, '    points found in the initial model file'

    ! now read in the number of variables
    read (99, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) nvars_model_file

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

    print *,'nBASE ',n_base
    print *,'nptsmodel ',npts_model

    dr_in = base_r(2) - base_r(1)

    dr = (prob_hi(2)-prob_lo(2)) / dble(n_base)

    print *,'DR OF MODEL.HSE   ',dr_in
    print *,'DR OF CALCULATION ',dr

    iratio = dr / dr_in
    offset = prob_lo(2) / dr_in + 1
    print *,'IRATIO ',iratio
    print *,'OFFSET ',offset

    print *,'RADIUS AT OFFSET ',base_r(offset), prob_lo(2)

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

       if (iratio .eq. 1) then
         d_ambient = base_state(j+offset,idens_model)
         t_ambient = base_state(j+offset,itemp_model)
         p_ambient = base_state(j+offset,ipres_model)
       else if (iratio .eq. 2) then
         d_ambient = HALF * (base_state(2*j+offset,idens_model) + base_state(2*j+1+offset,idens_model))
         t_ambient = HALF * (base_state(2*j+offset,itemp_model) + base_state(2*j+1+offset,itemp_model))
         p_ambient = HALF * (base_state(2*j+offset,ipres_model) + base_state(2*j+1+offset,ipres_model))
       else if (iratio .eq. 4) then
         d_ambient = FOURTH * (base_state(4*j  +offset,idens_model) + base_state(4*j+1+offset,idens_model) &
                              +base_state(4*j+2+offset,idens_model) + base_state(4*j+3+offset,idens_model))
         t_ambient = FOURTH * (base_state(4*j  +offset,itemp_model) + base_state(4*j+1+offset,itemp_model) &
                              +base_state(4*j+2+offset,itemp_model) + base_state(4*j+3+offset,itemp_model))
         p_ambient = FOURTH * (base_state(4*j  +offset,ipres_model) + base_state(4*j+1+offset,ipres_model) &
                              +base_state(4*j+2+offset,ipres_model) + base_state(4*j+3+offset,ipres_model))
!        print *,'R IS ',j,r
!        print *,'AVERAGING RS: ',base_r(4*j+offset)
!        print *,'AVERAGING RS: ',base_r(4*j+1+offset)
!        print *,'AVERAGING RS: ',base_r(4*j+2+offset)
!        print *,'AVERAGING RS: ',base_r(4*j+3+offset)
       else 
         d_ambient = interpolate(r, npts_model, base_r, base_state(:,idens_model))
         t_ambient = interpolate(r, npts_model, base_r, base_state(:,itemp_model))
         p_ambient = interpolate(r, npts_model, base_r, base_state(:,ipres_model))
       end if

       do n = 1, nspec
          xn_ambient(n) = interpolate(r, npts_model, base_r, base_state(:,ispec_model-1+n))
       enddo

       ! use the EOS to make the state consistent
       temp_eos(1) = t_ambient
       den_eos(1)  = d_ambient
       p_eos(1)    = p_ambient

       call eos(eos_input_rt, den_eos, temp_eos, &
                npts, nspec, &
                xn_ambient, &
                p_eos, h_eos, e_eos, &
                cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                dpdX_eos, dhdX_eos, &
                gam1_eos, cs_eos, s_eos, &
                dsdt_eos, dsdr_eos, &
                do_diag)
       
       s0(j, rho_comp ) = d_ambient
       s0(j,rhoh_comp ) = d_ambient * h_eos(1)
       s0(j,spec_comp:spec_comp+nspec-1) = d_ambient * xn_ambient(1:nspec)
       p0(j)    = p_eos(1)
       temp0(j) = t_ambient

       ! keep track of the height where we drop below the cutoff density
       if (s0(j,rho_comp) .lt. cutoff_density .and. j_cutoff .eq. n_base-1) j_cutoff = j

    end do
 

    ! KEEP RHO0 CONSTANT ABOVE J_CUTOFF
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

       den_eos(1)  = s0(j,rho_comp)
       xn_ambient(1:nspec) = s0(j,spec_comp:spec_comp-1+nspec)/den_eos(1)
       temp_eos(1) = temp0(j)
       p_eos(1)    =   p0(j)

       call eos(eos_input_rp, den_eos, temp_eos, &
                npts, nspec, &
                xn_ambient, &
                p_eos, h_eos, e_eos, &
                cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                dpdX_eos, dhdX_eos, &
                gam1_eos, cs_eos, s_eos, &
                dsdt_eos, dsdr_eos, &
                do_diag)
       
       temp0(j) = temp_eos(1)
       gam1(j)  = gam1_eos(1)
       s0(j,rhoh_comp) = h_eos(1) * s0(j,rho_comp)

       s0(j,trac_comp:) = ZERO
       
    end do

  end subroutine init_base_state

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
