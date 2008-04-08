module init_module

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_physbc_module
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

  private
  public :: initscalardata, initveldata, scalar_diags

contains

  subroutine initscalardata(nlevs,s,s_background,p0_background,dx,bc,mla)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: s_background(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_background(:,0:)
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
          if ( multifab_remote(s(n),i) ) cycle
          sop => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx(n,:), s_background(n,:,:), p0_background(n,:))
          case (3)
             call initscalardata_3d(n,sop(:,:,:,:), lo, hi, ng, dx(n,:), s_background(n,:,:), p0_background(n,:))
          end select
       end do
    enddo

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

  subroutine initscalardata_2d(s,lo,hi,ng,dx,s_background,p0_background)

    use probin_module, only: prob_lo_x, prob_lo_y, perturb_model

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s_background(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_background(0:)

    ! Local variables
    integer         :: i,j
    real(kind=dp_t) :: x,y
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)

    ! initial the domain with the base state
    s = ZERO

    ! initialize the scalars
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          s(i,j,rho_comp)  = s_background(j,rho_comp)
          s(i,j,rhoh_comp) = s_background(j,rhoh_comp)
          s(i,j,temp_comp) = s_background(j,temp_comp)
          s(i,j,spec_comp:spec_comp+nspec-1) = &
               s_background(j,spec_comp:spec_comp+nspec-1)
          s(i,j,trac_comp:trac_comp+ntrac-1) = &
               s_background(j,trac_comp:trac_comp+ntrac-1)
       enddo
    enddo
    
    ! add an optional perturbation
    if (perturb_model) then
       do j = lo(2), hi(2)
          y = prob_lo_y + (dble(j)+HALF) * dx(2)
          
          do i = lo(1), hi(1)
             x = prob_lo_x + (dble(i)+HALF) * dx(1)
          
             call perturb_2d(x, y, p0_background(j), s_background(j,:), &
                             dens_pert, rhoh_pert, rhoX_pert, temp_pert, trac_pert)

             s(i,j,rho_comp) = dens_pert
             s(i,j,rhoh_comp) = rhoh_pert
             s(i,j,temp_comp) = temp_pert
             s(i,j,spec_comp:spec_comp+nspec-1) = rhoX_pert(1:)
             s(i,j,trac_comp:trac_comp+ntrac-1) = trac_pert(:)
          enddo
       enddo
    endif
    
  end subroutine initscalardata_2d

  subroutine initscalardata_3d(n,s,lo,hi,ng,dx,s_background,p0_background)

    use probin_module, only: prob_lo_x, prob_lo_y, prob_lo_z, perturb_model
    
    integer           , intent(in   ) :: n,lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s_background(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_background(0:)

    !     Local variables
    integer         :: i,j,k,comp
    real(kind=dp_t) :: x,y,z
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)

    ! initial the domain with the base state
    s = ZERO
  
    if (spherical .eq. 1) then

       ! initialize the scalars
       call fill_3d_data(n,s(:,:,:,rho_comp), s_background(:,rho_comp), lo,hi,dx,ng)
       call fill_3d_data(n,s(:,:,:,rhoh_comp),s_background(:,rhoh_comp),lo,hi,dx,ng)
       call fill_3d_data(n,s(:,:,:,temp_comp),s_background(:,temp_comp),lo,hi,dx,ng)

       do comp = spec_comp, spec_comp+nspec-1
          call fill_3d_data(n,s(:,:,:,comp),s_background(:,comp),lo,hi,dx,ng)
       end do

       do comp = trac_comp, trac_comp+ntrac-1
          call fill_3d_data(n,s(:,:,:,comp),s_background(:,comp),lo,hi,dx,ng)
       end do

    else 

       ! initialize the scalars
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                s(i,j,k,rho_comp)  = s_background(k,rho_comp)
                s(i,j,k,rhoh_comp) = s_background(k,rhoh_comp)
                s(i,j,k,temp_comp) = s_background(k,temp_comp)
                s(i,j,k,spec_comp:spec_comp+nspec-1) = &
                     s_background(k,spec_comp:spec_comp+nspec-1)
                s(i,j,k,trac_comp:trac_comp+ntrac-1) = &
                     s_background(k,trac_comp:trac_comp+ntrac-1)
             enddo
          enddo
       enddo
       
       if (perturb_model) then

          ! add an optional perturbation
          do k = lo(3), hi(3)
             z = prob_lo_z + (dble(k)+HALF) * dx(3)
             
             do j = lo(2), hi(2)
                y = prob_lo_y + (dble(j)+HALF) * dx(2)
                
                do i = lo(1), hi(1)
                   x = prob_lo_x + (dble(i)+HALF) * dx(1)
                   
                   call perturb_3d(x, y, z, p0_background(k), s_background(k,:), &
                                   dens_pert, rhoh_pert, rhoX_pert, temp_pert, trac_pert)

                   s(i,j,k,rho_comp) = dens_pert
                   s(i,j,k,rhoh_comp) = rhoh_pert
                   s(i,j,k,temp_comp) = temp_pert
                   s(i,j,k,spec_comp:spec_comp+nspec-1) = rhoX_pert(:)
                   s(i,j,k,trac_comp:trac_comp+ntrac-1) = trac_pert(:)
                enddo
             enddo
          enddo
       endif

    end if
    
  end subroutine initscalardata_3d

  subroutine initveldata(nlevs,u,s_background,p0_background,dx,bc,mla)

    integer        , intent(in   ) :: nlevs
    type(multifab) , intent(inout) :: u(:)
    real(kind=dp_t), intent(in   ) :: s_background(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_background(:,0:)
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
          if ( multifab_remote(u(n),i) ) cycle
          uop => dataptr(u(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call initveldata_2d(uop(:,:,1,:), lo, hi, ng, dx(n,:), &
                                 s_background(n,:,:), p0_background(n,:))
          case (3) 
             call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx(n,:), &
                                 s_background(n,:,:), p0_background(n,:))
          end select
       end do
    enddo

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

  subroutine initveldata_2d(u,lo,hi,ng,dx,s_background,p0_background)

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s_background(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_background(0:)

    ! Local variables

    ! initial the velocity
    u = ZERO

  end subroutine initveldata_2d

  subroutine initveldata_3d(u,lo,hi,ng,dx,s_background,p0_background)

    integer           , intent(in   ) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s_background(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_background(0:)

    ! Local variables

    ! initial the velocity
    u = ZERO
    
  end subroutine initveldata_3d


  subroutine perturb_2d(x, y, p0_background, s_background, dens_pert, rhoh_pert, rhoX_pert, temp_pert, trac_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y
    real(kind=dp_t), intent(in ) :: p0_background, s_background(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp,t0
    real(kind=dp_t) :: x0, y0, x1, y1, x2, y2
    real(kind=dp_t) :: r0, r1, r2

    t0 = s_background(temp_comp)

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
          
    ! Use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    temp_eos(1) = temp
    p_eos(1) = p0_background
    den_eos(1) = s_background(rho_comp)
    xn_eos(1,:) = s_background(spec_comp:spec_comp+nspec-1)/s_background(rho_comp)

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
    rhoX_pert = dens_pert*xn_eos(1,:)

    temp_pert = temp
    
!   if ( (r0 .lt. 2.0) .or. (r1 .lt. 2.0) .or. (r2 .lt. 2.0) ) then
!     trac_pert = ONE
!   else
      trac_pert = ZERO
!   end if

  end subroutine perturb_2d

  subroutine perturb_3d(x, y, z, p0_background, s_background, dens_pert, rhoh_pert, &
                        rhoX_pert, temp_pert, trac_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y, z
    real(kind=dp_t), intent(in ) :: p0_background, s_background(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp, t0
    real(kind=dp_t) :: x0, y0, z0, x1, y1, z1, x2, y2, z2
    real(kind=dp_t) :: r0, r1, r2

    t0 = s_background(temp_comp)

    x0 = 5.0d7
    y0 = 5.0d7
    z0 = 6.5d7
    
    x1 = 1.2d8
    y1 = 1.2d8
    z1 = 8.5d7
    
    x2 = 2.0d8
    y2 = 2.0d8
    z2 = 7.5d7

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
    temp_eos(1) = temp
    p_eos(1) = p0_background
    den_eos(1) = s_background(rho_comp)
    xn_eos(1,:) = s_background(spec_comp:spec_comp+nspec-1)/s_background(rho_comp)

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
    rhoX_pert = dens_pert*xn_eos(1,:)

    temp_pert = temp
    
!   if (r1 .lt. 2.0) then
!     trac_pert = ONE
!   else
      trac_pert = ZERO
!   end if

  end subroutine perturb_3d

  subroutine scalar_diags (istep,s,s_background,p0_background,dx)

    integer        , intent(in   ) :: istep
    type(multifab) , intent(inout) :: s
    real(kind=dp_t), intent(in)    :: s_background(:,:)
    real(kind=dp_t), intent(in)    :: p0_background(:)
    real(kind=dp_t), intent(in)    :: dx(:)

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng,dm
    integer :: i
    
    ng = s%ng
    dm = s%dim

    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       sop => dataptr(s, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))

       select case (dm)
       case (2)
          call scalar_diags_2d(istep, sop(:,:,1,:), lo, hi, ng, dx, s_background, p0_background)
       case (3)
!         call scalar_diags_3d(istep, sop(:,:,:,:), lo, hi, ng, dx, s_background)
       end select
    end do

  end subroutine scalar_diags

  subroutine scalar_diags_2d (istep, s,lo,hi,ng,dx,s_background,p0_background)

    use probin_module, only: grav_const

    integer, intent(in) :: istep, lo(:), hi(:), ng
    real (kind = dp_t), intent(in) ::  s(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in) :: dx(:)
    real(kind=dp_t)   , intent(in) :: s_background(0:,:)
    real(kind=dp_t)   , intent(in) :: p0_background(0:)

    ! Local variables
    integer :: i, j
    real(kind=dp_t) :: fac, mass, mass0
    real(kind=dp_t), allocatable ::  rhoavg(:)
    real(kind=dp_t), allocatable :: rhopert(:)
    real(kind=dp_t), allocatable ::    pavg(:)
    character(len=11) :: file_name
    character(len=10) :: file_name2
    character(len= 8) :: file_name3

    allocate(rhopert(lo(2):hi(2)))
    allocate( rhoavg(lo(2):hi(2)))
    allocate(   pavg(lo(2):hi(2)))

    write(unit=file_name ,fmt='("rhopert",i4.4)') istep
    write(unit=file_name2,fmt='("rhoavg",i4.4)') istep
    write(unit=file_name3,fmt='("pavg",i4.4)') istep
    open(90,file=file_name)
    open(91,file=file_name2)
    open(92,file=file_name3)

    fac = ONE / dble(hi(1)-lo(1)+1)
    mass  = ZERO
    mass0 = ZERO
    do j = lo(2), hi(2)
      rhoavg(j) = ZERO
      rhopert(j) = ZERO
      do i = lo(1), hi(1)
         rhopert(j) = rhopert(j) + (s(i,j,rho_comp) - s_background(j,rho_comp))
         rhoavg(j) = rhoavg(j) +  s(i,j,rho_comp)
      enddo
      rhoavg(j)  = rhoavg(j) * fac
      rhopert(j)  = rhopert(j) * fac
      write(90,*) (dble(j)+HALF)*dx(2),rhopert(j)
      write(91,*) (dble(j)+HALF)*dx(2),rhoavg(j)
      mass  = mass  + rhoavg(j)
      mass0 = mass0 + s_background(j,rho_comp)
    enddo

!   print *,'TOTAL MASS ',istep, mass, mass0

    pavg(hi(2)) = p0_background(hi(2))
    do j = hi(2)-1,lo(2),-1
      pavg(j) = pavg(j+1) + 0.5d0 * (rhoavg(j+1)+rhoavg(j))*abs(grav_const)*dx(2)
    enddo
    do j = lo(2),hi(2)
      write(92,*) (dble(j)+HALF)*dx(2),p0_background(j),pavg(j)
    enddo

    deallocate(rhoavg,rhopert,pavg)

  end subroutine scalar_diags_2d


end module init_module
