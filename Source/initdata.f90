module init_module

  use bl_types
  use bl_constants_module
  use bc_module
  use inlet_bc
  use setbc_module
  use define_bc_module
  use multifab_module
  use make_div_coeff_module
  use eos_module
  use variables
  use network

  implicit none

contains

  subroutine initdata (u,s,rho0,rhoh0,p0,temp0,dx,prob_lo,prob_hi,bc,nscal)

    type(multifab) , intent(inout) :: u,s
    real(kind=dp_t), intent(in   ) ::  rho0(:)
    real(kind=dp_t), intent(in   ) :: rhoh0(:)
    real(kind=dp_t), intent(in   ) ::    p0(:)
    real(kind=dp_t), intent(in   ) :: temp0(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: prob_lo(:)
    real(kind=dp_t), intent(in   ) :: prob_hi(:)
    type(bc_level) , intent(in   ) :: bc
    integer        , intent(in   ) :: nscal

    real(kind=dp_t), pointer:: uop(:,:,:,:), sop(:,:,:,:)
    integer :: lo(u%dim),hi(u%dim),ng,dm
    integer :: i,n
    
    ng = u%ng
    dm = u%dim

    do i = 1, u%nboxes
       if ( multifab_remote(u, i) ) cycle
       uop => dataptr(u, i)
       sop => dataptr(s, i)
       lo =  lwb(get_box(u, i))
       hi =  upb(get_box(u, i))

       select case (dm)
       case (2)
          call initdata_2d(uop(:,:,1,:), sop(:,:,1,:), lo, hi, ng, dx, &
                           prob_lo, prob_hi, rho0, rhoh0, p0, temp0)
          do n = 1,dm
             call setbc_2d(uop(:,:,1,n), lo, ng, bc%adv_bc_level_array(i,:,:,   n),dx,   n)
          end do
          do n = 1,nscal
             call setbc_2d(sop(:,:,1,n), lo, ng, bc%adv_bc_level_array(i,:,:,dm+n),dx,dm+n)
          end do

       case (3)
          call initdata_3d(uop(:,:,:,:), sop(:,:,:,:), lo, hi, ng, dx, &
                           prob_lo, prob_hi, rho0, rhoh0, p0, temp0)
          do n = 1, dm
             call setbc_3d(uop(:,:,:,n), lo, ng, bc%adv_bc_level_array(i,:,:,   n),dx,   n)
          end do
          do n = 1, nscal
             call setbc_3d(sop(:,:,:,n), lo, ng, bc%adv_bc_level_array(i,:,:,dm+n),dx,dm+n)
          end do
       end select
    end do

    call multifab_fill_boundary(u)
    call multifab_fill_boundary(s)

  end subroutine initdata

  subroutine initdata_2d (u,s,lo,hi,ng,dx,prob_lo,prob_hi,rho0,rhoh0,p0,temp0)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real (kind = dp_t), intent(in ) :: prob_lo(:)
    real (kind = dp_t), intent(in ) :: prob_hi(:)
    real(kind=dp_t), intent(in   ) ::  rho0(lo(2):)
    real(kind=dp_t), intent(in   ) :: rhoh0(lo(2):)
    real(kind=dp_t), intent(in   ) ::    p0(lo(2):)
    real(kind=dp_t), intent(in   ) :: temp0(lo(2):)

    !     Local variables
    integer :: i, j, n
    real(kind=dp_t) :: x,y,r,r0,r1,r2,temp
    real(kind=dp_t) :: dens_pert, rhoh_pert, xn_pert(nspec)

    logical, parameter :: perturbModel = .true.

    ! initial the domain with the base state
    u = ZERO
    s = ZERO

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          s(i,j,1) = rho0(j)
          s(i,j,2) = rhoh0(j)
       enddo
    enddo

    
    ! add an optional perturbation
    do j = lo(2), hi(2)
       y = prob_lo(2) + (dble(j)+HALF) * dx(2)
       
       do i = lo(1), hi(1)
          x = prob_lo(1) + (dble(i)+HALF) * dx(1)
          
          if (perturbModel) then
             call perturb_2d(x, y, temp0(j), p0(j), rho0(j), dens_pert, rhoh_pert, xn_pert)
             s(i,j,1) = dens_pert
             s(i,j,2) = rhoh_pert
             s(i,j,spec_comp:spec_comp+nspec-1) = dens_pert*xn_pert(:)
          endif

       enddo
    enddo

    
    if (size(s,dim=3).gt.2) then
       do n = 3, size(s,dim=3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                s(i,j,n) = ONE
             end do
          end do
       end do
    end if
    
  end subroutine initdata_2d

  subroutine initdata_3d (u,s,lo,hi,ng,dx,prob_lo,prob_hi,rho0,rhoh0,p0,temp0)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(out) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    real (kind = dp_t), intent(in ) :: prob_lo(:)
    real (kind = dp_t), intent(in ) :: prob_hi(:)
    real(kind=dp_t), intent(in   ) ::  rho0(lo(3):)
    real(kind=dp_t), intent(in   ) :: rhoh0(lo(3):)
    real(kind=dp_t), intent(in   ) ::    p0(lo(3):)
    real(kind=dp_t), intent(in   ) :: temp0(lo(3):)

    !     Local variables
    integer :: i, j, k, n
    real(kind=dp_t) :: x,y,z,r,r0,r1,r2,temp
    real(kind=dp_t) :: dens_pert, rhoh_pert, xn_pert(nspec)
    logical, parameter :: perturbModel = .true.

    ! initial the domain with the base state
    u = ZERO
    s = ZERO

    do k = lo(3), hi(3)
    do j = lo(2), hi(2)
    do i = lo(1), hi(1)
       s(i,j,k,1) = rho0(k)
       s(i,j,k,2) = rhoh0(k)
    enddo
    enddo
    enddo

    
    ! add an optional perturbation
    do k = lo(3), hi(3)
       z = prob_lo(3) + (dble(k)+HALF) * dx(3)
       
       do j = lo(2), hi(2)
        y = prob_lo(2) + (dble(j)+HALF) * dx(2)
        do i = lo(1), hi(1)
          x = prob_lo(1) + (dble(i)+HALF) * dx(1)
          
          if (perturbModel) then
             call perturb_3d(x, y, z, temp0(k), p0(k), rho0(k), dens_pert, rhoh_pert, xn_pert)
             s(i,j,k,1) = dens_pert
             s(i,j,k,2) = rhoh_pert
             s(i,j,k,spec_comp:spec_comp+nspec-1) = dens_pert*xn_pert(:)
          endif
        enddo
       enddo
    enddo
    
    if (size(s,dim=4).gt.2) then
       do n = 4, size(s,dim=4)
          do k = lo(3), hi(3)
          do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             s(i,j,k,n) = ONE
          end do
          end do
          end do
       end do
    end if
    
  end subroutine initdata_3d

  subroutine perturb_2d(x, y, t0, p0, rho0, dens_pert, rhoh_pert, xn_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y
    real(kind=dp_t), intent(in ) :: t0, p0, rho0
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert
    real(kind=dp_t), dimension(nspec), intent(out) :: xn_pert

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
    
    temp = t0 * (ONE + TWO * ( &
         .0625_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r0))) + &
         .1875_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r1))) + &
         .1250_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r2)))  ) )
          
    ! Use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    temp_row(1) = temp
    p_row(1) = p0
    den_row(1) = rho0
    input_flag = 3      ! (t, p) -> (rho, h)

    call eos(input_flag, den_row, temp_row, npts, nspec, &
         xmass, aion, zion, &
         p_row, h_row, e_row, &
         cv_row, cp_row, xne_row, eta_row, &
         pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
         s_row,do_diag)

    dens_pert = den_row(1)
    rhoh_pert = den_row(1)*h_row(1)

  end subroutine perturb_2d

  subroutine perturb_3d(x, y, z, t0, p0, rho0, dens_pert, rhoh_pert, xn_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y, z
    real(kind=dp_t), intent(in ) :: t0, p0, rho0
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert
    real(kind=dp_t), dimension(nspec), intent(out) :: xn_pert

    real(kind=dp_t) :: temp
    real(kind=dp_t) :: x0, y0, z0, x1, y1, z1, x2, y2, z2
    integer :: i, j, k
    real(kind=dp_t) :: r0, r1, r2

    y0 = 5.0d7
    z0 = 6.5d7
    
    y1 = 1.2d8
    z1 = 8.5d7
    
    y2 = 2.0d8
    z2 = 7.5d7

    ! Tanh bubbles
    r0 = sqrt( (y-y0)**2 +(z-z0)**2 ) / 2.5e6
    r1 = sqrt( (y-y1)**2 +(z-z1)**2 ) / 2.5e6
    r2 = sqrt( (y-y2)**2 +(z-z2)**2 ) / 2.5e6
    
    temp = t0 * (ONE + TWO * ( &
         .0625_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r0))) + &
         .1875_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r1))) + &
         .1250_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0-r2)))  ) )
          
    ! Use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    temp_row(1) = temp
    p_row(1) = p0
    den_row(1) = rho0
    input_flag = 3      ! (t, p) -> (rho, h)

    call eos(input_flag, den_row, temp_row, npts, nspec, &
         xmass, aion, zion, &
         p_row, h_row, e_row, &
         cv_row, cp_row, xne_row, eta_row, &
         pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
         s_row,do_diag)

    dens_pert = den_row(1)
    rhoh_pert = den_row(1)*h_row(1)

  end subroutine perturb_3d

  subroutine init_base_state (div_coef_type,div_coeff,div_coeff_half,gam1, &
                              rho0,temp0,rhoh0,p0,dx,prob_lo,prob_hi,grav,anelastic_cutoff)

    integer        , intent(in   ) :: div_coef_type
    real(kind=dp_t), intent(inout) :: div_coeff(0:)
    real(kind=dp_t), intent(inout) :: div_coeff_half(0:)
    real(kind=dp_t), intent(inout) ::  gam1(0:)
    real(kind=dp_t), intent(inout) ::  rho0(0:)
    real(kind=dp_t), intent(inout) :: temp0(0:)
    real(kind=dp_t), intent(inout) :: rhoh0(0:)
    real(kind=dp_t), intent(inout) ::    p0(0:)
    real(kind=dp_t), intent(in   ) :: prob_lo(:)
    real(kind=dp_t), intent(in   ) :: prob_hi(:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), intent(in   ) :: grav
    real(kind=dp_t), intent(in   ) :: anelastic_cutoff
    
    integer :: i,j,nx,j_cutoff

    real(kind=dp_t) :: r
    real(kind=dp_t) :: d_ambient,t_ambient,p_ambient
    real(kind=dp_t) :: integral, temp_term_lo, temp_term_hi
    real(kind=dp_t) :: temp_min,p0_lo,p0_hi

    integer, parameter :: nvars_model = 6
    integer :: npts_model
    real(kind=dp_t), allocatable :: base_state(:,:), base_r(:)
    integer :: ipos
    character (len=256) :: header_line

    logical :: do_diag

    real(kind=dp_t), parameter :: cutoff_density = 2.5d6

    nx = size(div_coeff,dim=1)

    do_diag = .false.

    ! open the model file and read in the header
    ! the first line has the number of points in the model
    open(99,file="model.hse")
    read (99, '(a256)') header_line
    ipos = index(header_line, '=') + 1
    read (header_line(ipos:),*) npts_model

    ! the base state from the model input file is contained in two arrays.
    ! base_r(:) holds the coordinate information and base_state(:,var) holds
    ! variable var as a function of height.
    !     base_state(i,1) =  rho(i)
    !     base_state(i,2) = temp(i)
    !     base_state(i,3) = pres(i)
    !     base_state(i,4) =  c12(i)
    !     base_state(i,5) =  o16(i)

    allocate (base_state(npts_model, nvars_model))
    allocate (base_r(npts_model))

    print *, '<<< npts_model = ', npts_model, ' >>>'

    do i = 1, npts_model
       read(99,*) base_r(i),base_state(i,1),base_state(i,2),base_state(i,3), &
            base_state(i,4),base_state(i,5)
    end do
    close(99)

    call helmeos_init


    print *,'DIV_COEF_TYPE ',div_coef_type
    j_cutoff = nx-1
    do j = 0,nx-1

       ! compute the coordinate height at this level
       ! NOTE: we are assuming that the basestate is in the y-direction
       ! and that ymin = 0.0
       r = prob_lo(2) + (dble(j) + HALF)*dx(2)

       d_ambient = interpolate(r, npts_model, base_r, base_state(:,1))
       t_ambient = interpolate(r, npts_model, base_r, base_state(:,2))
       p_ambient = interpolate(r, npts_model, base_r, base_state(:,3))


       ! use the EOS to make the state consistent
       temp_row(1) = t_ambient
       den_row(1)  = d_ambient
       p_row(1)    = p_ambient

       ! (rho,T) --> p,h
       input_flag = 1

       call eos(input_flag, den_row, temp_row, npts, nspec, &
            xmass, aion, zion, &
            p_row, h_row, e_row, &
            cv_row, cp_row, xne_row, eta_row, &
            pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
            s_row,do_diag)

       
       rho0(j)  = d_ambient
       p0(j)    = p_row(1)
       temp0(j) = t_ambient
       rhoh0(j) = d_ambient * h_row(1)

       ! keep track of the height where we drop below the cutoff density
       if (rho0(j) .lt. cutoff_density .and. j_cutoff .eq. nx-1) j_cutoff = j

    end do
 

    ! KEEP RHO0 CONSTANT ABOVE J_CUTOFF
    print *,'JCUT_comp ',j_cutoff

    do j = j_cutoff,nx-1
       rho0(j)  =  rho0(j_cutoff)
       p0(j)    =    p0(j_cutoff)
       temp0(j) = temp0(j_cutoff)
       rhoh0(j) = rhoh0(j_cutoff)
    end do

    ! RECALCULATE P0 ABOVE J_CUTOFF
!     p0_lo = HALF * (p0(j_cutoff) + p0(j_cutoff+1))
!     do j = j_cutoff+1,nx-1
!       p0_hi = p0_lo - abs(grav) * rho0(j) * dx(2)
!       p0(j) = HALF * (p0_lo + p0_hi)
!       p0_lo = p0_hi
!     end do

      ! RECALCULATE T, RHO_H

    do j = 0,nx-1

       den_row(1)  = rho0(j)
       temp_row(1) = temp0(j)
       p_row(1)    =   p0(j)

!       (rho, p) --> T, h
       input_flag = 4

       call eos(input_flag, den_row, temp_row, npts, nspec, &
            xmass, aion, zion, &
            p_row, h_row, e_row, &
            cv_row, cp_row, xne_row, eta_row, &
            pele_row, dpdt_row, dpdr_row, dedt_row, dedr_row, gam1_row, cs_row, &
            s_row,do_diag)
       
       temp0(j) = temp_row(1)
       gam1(j)  = gam1_row(1)
       rhoh0(j) = h_row(1) * rho0(j)

    end do


    ! compute the coefficient inside the divergence constraint
    !     div_coef_type = 1   INCOMPRESSIBLE
    !     div_coef_type = 2     ANELASTIC 
    !     div_coef_type = 3        P-I
    
    if (div_coef_type .eq. 2) then
       div_coeff = rho0
       div_coeff_half(   0) = rho0( 0)
       div_coeff_half(nx+1) = rho0(nx)
       do j = 1,nx
          div_coeff_half(j) = HALF * (rho0(j) + rho0(j-1))
       end do
    else
       call make_div_coeff(div_coeff,div_coeff_half,rho0,p0,gam1,grav,dx(2),anelastic_cutoff)
    end if

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


  subroutine write_base_state(sd_name,rho0,rhoh0,temp0,p0,div_coeff)

    character(len=9), intent(in) :: sd_name
    real(kind=dp_t) , intent(in) :: rho0(:),p0(:),rhoh0(:),temp0(:),div_coeff(:)

    integer :: i

    print *,'Writing base state to ',sd_name

    open(unit=99,file=sd_name,form = "formatted", access = "sequential",action="write")
    do i = 1, size(rho0,dim=1)
       write(99,1000)  rho0(i), p0(i), rhoh0(i), temp0(i), div_coeff(i)
    end do
    close(99)

1000 format(e18.12,1x,e18.12,1x,e18.12,1x,e18.12,1x,e18.12)

  end subroutine write_base_state


  subroutine read_base_state(sd_name,rho0,rhoh0,temp0,p0,div_coeff)
    
    character(len=9), intent(in   ) :: sd_name
    real(kind=dp_t) , intent(inout) :: rho0(:),p0(:),rhoh0(:),temp0(:),div_coeff(:)

    integer :: i

    print *,'Reading base state from ',sd_name

    open(unit=99,file=sd_name)
    do i = 1, size(rho0,dim=1)
       read(99,*)  rho0(i), p0(i), rhoh0(i), temp0(i), div_coeff(i)
    end do
    close(99)

  end subroutine read_base_state


  subroutine impose_pressure_bcs(p,mla,mult)

    type(multifab ), intent(inout) :: p(:)
    type(ml_layout), intent(in   ) :: mla
    real(kind=dp_t), intent(in   ) :: mult
    
    type(box)           :: bx,pd
    integer             :: i,n,nlevs
     
    nlevs = size(p,dim=1)

    do n = 1,nlevs
       pd = layout_get_pd(mla%la(n))
       do i = 1, p(n)%nboxes; if ( remote(p(n),i) ) cycle
          bx = get_ibox(p(n),i)
          if (bx%lo(2) == pd%lo(2)) then
             bx%hi(2) = bx%lo(2)
             call setval(p(n),mult,bx)
          end if
       end do
    end do
    
  end subroutine impose_pressure_bcs

end module init_module
