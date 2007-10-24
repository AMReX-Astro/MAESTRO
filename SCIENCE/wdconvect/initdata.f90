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

  implicit none

contains

  subroutine initscalardata (s,s0,p0,dx,perturb_model, &
                             prob_lo,prob_hi,bc)

    type(multifab) , intent(inout) :: s
    real(kind=dp_t), intent(in   ) ::    s0(0:,:)
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
       select case (dm)
       case (2)
          do n = 1,nscal
             call setbc_2d(sop(:,:,1,n), lo, ng, &
                           bc%adv_bc_level_array(i,:,:,dm+n),dx,dm+n)
          end do
       case (3)
          do n = 1, nscal
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
    real(kind=dp_t), intent(in   ) ::    s0(0:,:)
    real(kind=dp_t), intent(in   ) ::    p0(0:)

    !     Local variables
    integer :: i, j, n
    real(kind=dp_t) :: x,y,r,r0,r1,r2,temp
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)

    ! initial the domain with the base state
    s = ZERO

    ! initialize the scalars
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          s(i,j,rho_comp)  = s0(j,rho_comp)
          s(i,j,rhoh_comp) = s0(j,rhoh_comp)
          s(i,j,temp_comp) = s0(j,temp_comp)

          s(i,j,spec_comp:spec_comp+nspec-1) = &
               s0(j,spec_comp:spec_comp+nspec-1)
       enddo
    enddo
    
    ! add an optional perturbation
    if (perturb_model) then
       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j)+HALF) * dx(2)
       
          do i = lo(1), hi(1)
             x = prob_lo(1) + (dble(i)+HALF) * dx(1)
          
             call perturb_2d(x, y, p0(j), s0(j,:), &
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

  subroutine initscalardata_3d (s,lo,hi,ng,dx, perturb_model, &
                                prob_lo,prob_hi,s0,p0)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in ) :: dx(:)
    logical,            intent(in ) :: perturb_model
    real (kind = dp_t), intent(in ) :: prob_lo(:)
    real (kind = dp_t), intent(in ) :: prob_hi(:)
    real(kind=dp_t), intent(in   ) ::    s0(0:,:)
    real(kind=dp_t), intent(in   ) ::    p0(0:)

    !     Local variables
    integer :: i, j, k, n
    real(kind=dp_t) :: x,y,z,r,r0,r1,r2,temp
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)

    ! initial the domain with the base state
    s = ZERO
  
    if (spherical .eq. 1) then

       ! initialize the scalars
       call fill_3d_data (s(:,:,:,rho_comp), s0(:,rho_comp), lo,hi,dx,ng)
       call fill_3d_data (s(:,:,:,rhoh_comp),s0(:,rhoh_comp),lo,hi,dx,ng)
       call fill_3d_data (s(:,:,:,temp_comp),s0(:,temp_comp),lo,hi,dx,ng)

       do n = spec_comp, spec_comp+nspec-1
          call fill_3d_data (s(:,:,:,n),s0(:,n),lo,hi,dx,ng)
       end do

    else 

       ! initialize the scalars
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                s(i,j,k,rho_comp)  = s0(k,rho_comp)
                s(i,j,k,rhoh_comp) = s0(k,rhoh_comp)
                s(i,j,k,temp_comp) = s0(k,temp_comp)

                s(i,j,k,spec_comp:spec_comp+nspec-1) = &
                     s0(k,spec_comp:spec_comp+nspec-1)
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
                   
                   call perturb_3d(x, y, z, p0(k), s0(k,:), &
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

    !     Local variables
    integer :: i, j, n
    real(kind=dp_t) :: x,y,r,r0,r1,r2,temp
    real(kind=dp_t) :: dens_pert, rhoh_pert, rhoX_pert(nspec), trac_pert(ntrac)

    ! initial the velocity
    u = ZERO

  end subroutine initveldata_2d

  ! the velocity is initialized to zero plus a perturbation which is a
  ! summation of 27 fourier modes with random amplitudes and phase
  ! shifts over a square of length "velpert_scale".  The parameter
  ! "velpert_radius" is the cutoff radius, such that the actual
  ! perturbation that gets added to the initial velocity decays to
  ! zero quickly if r > "velpert_radius", via the tanh function.
  ! The steepness of the cutoff is controlled by "velpert_steep".  The
  ! relative amplitude of the modes is controlled by
  ! "velpert_amplitude".
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

    ! Local variables
    integer :: i, j, k
    integer :: iloc, jloc, kloc

    ! L2 norm of k
    real(kind=dp_t) :: normk(3,3,3)

    ! random numbers between -1 and 1
    real(kind=dp_t) :: alpha(3,3,3), beta(3,3,3), gamma(3,3,3)

    ! random numbers between 0 and 2*pi
    real(kind=dp_t) :: phix(3,3,3), phiy(3,3,3), phiz(3,3,3)

    ! cos and sin of (2*pi*kx/L + phix), etc
    real(kind=dp_t) :: cx(3,3,3), cy(3,3,3), cz(3,3,3)
    real(kind=dp_t) :: sx(3,3,3), sy(3,3,3), sz(3,3,3)

    ! location of center of star
    real(kind=dp_t) :: xc(3)

    ! radius, or distance, to center of star
    real(kind=dp_t) :: rloc

    ! the point we're at
    real(kind=dp_t) :: xloc(3)

    ! perturbational velocity to add
    real(kind=dp_t) :: upert(3)

    ! random number
    real(kind=dp_t) :: rand

    ! initialize the velocity to zero everywhere
    u = ZERO

    ! load in random numbers alpha, beta, gamma, phix, phiy, and phiz
    do i=1,3
       do j=1,3
          do k=1,3
             call random_number(rand)
             rand = 2.0d0*rand - 1.0d0
             alpha(i,j,k) = rand
             call random_number(rand)
             rand = 2.0d0*rand - 1.0d0
             beta(i,j,k) = rand
             call random_number(rand)
             rand = 2.0d0*rand - 1.0d0
             gamma(i,j,k) = rand
             call random_number(rand)
             rand = 2.0d0*M_PI*rand
             phix(i,j,k) = rand
             call random_number(rand)
             rand = 2.0d0*M_PI*rand
             phiy(i,j,k) = rand
             call random_number(rand)
             rand = 2.0d0*M_PI*rand
             phiz(i,j,k) = rand
          enddo
       enddo
    enddo

    ! compute the norm of k
    do i=1,3
       do j=1,3
          do k=1,3
             normk(i,j,k) = sqrt(dble(i)**2+dble(j)**2+dble(k)**2)
          enddo
       enddo
    enddo

    ! define where center of star is
    ! this currently assumes the star is at the center of the domain
    do i=1,3
       xc(i) = 0.5d0*(prob_lo(i)+prob_hi(i))
    enddo

    ! now do the big loop over all points in the domain
    do iloc = lo(1),hi(1)
       do jloc = lo(2),hi(2)
          do kloc = lo(3),hi(3)

             ! set perturbational velocity to zero
             upert = ZERO

             ! compute where we physically are
             xloc(1) = prob_lo(1) + (dble(iloc)+0.5d0)*dx(1)
             xloc(2) = prob_lo(2) + (dble(jloc)+0.5d0)*dx(2)
             xloc(3) = prob_lo(3) + (dble(kloc)+0.5d0)*dx(3)

             ! compute distance to the center of the star
             rloc = ZERO
             do i=1,3
                rloc = rloc + (xloc(i) - xc(i))**2
             enddo
             rloc = sqrt(rloc)

             ! loop over the 27 combinations of fourier components
             do i=1,3
                do j=1,3
                   do k=1,3
                      ! compute cosines and sines
                      cx(i,j,k) = cos(2.0d0*M_PI*dble(i)*xloc(1)/velpert_scale + phix(i,j,k))
                      cy(i,j,k) = cos(2.0d0*M_PI*dble(j)*xloc(2)/velpert_scale + phiy(i,j,k))
                      cz(i,j,k) = cos(2.0d0*M_PI*dble(k)*xloc(3)/velpert_scale + phiz(i,j,k))
                      sx(i,j,k) = sin(2.0d0*M_PI*dble(i)*xloc(1)/velpert_scale + phix(i,j,k))
                      sy(i,j,k) = sin(2.0d0*M_PI*dble(j)*xloc(2)/velpert_scale + phiy(i,j,k))
                      sz(i,j,k) = sin(2.0d0*M_PI*dble(k)*xloc(3)/velpert_scale + phiz(i,j,k))
                   enddo
                enddo
             enddo

             ! loop over the 27 combinations of fourier components
             do i=1,3
                do j=1,3
                   do k=1,3
                      ! compute contribution from perturbation velocity from each mode
                      upert(1) = upert(1) + &
                           (-gamma(i,j,k)*dble(j)*cx(i,j,k)*cz(i,j,k)*sy(i,j,k) &
                            +beta(i,j,k)*dble(k)*cx(i,j,k)*cy(i,j,k)*sz(i,j,k)) &
                            / normk(i,j,k)

                      upert(2) = upert(2) + &
                           (gamma(i,j,k)*dble(i)*cy(i,j,k)*cz(i,j,k)*sx(i,j,k) &
                            -alpha(i,j,k)*dble(k)*cx(i,j,k)*cy(i,j,k)*sz(i,j,k)) &
                            / normk(i,j,k)

                      upert(3) = upert(3) + &
                           (-beta(i,j,k)*dble(i)*cy(i,j,k)*cz(i,j,k)*sx(i,j,k) &
                            +alpha(i,j,k)*dble(j)*cx(i,j,k)*cz(i,j,k)*sy(i,j,k)) &
                            / normk(i,j,k)
                   enddo
                enddo
             enddo

             ! apply the cutoff function to the perturbational velocity
             do i=1,3
                upert(i) = velpert_amplitude*upert(i)*(0.5d0+0.5d0*tanh((velpert_radius-rloc)/velpert_steep))
             enddo

             ! add perturbational velocity to background velocity
             do i=1,3
                u(iloc,jloc,kloc,i) = u(iloc,jloc,kloc,i) + upert(i)
             enddo

          enddo
       enddo
    enddo
      
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
    xn_zone(:) = s0(spec_comp:spec_comp+nspec-1)/s0(rho_comp)

    ! (t, p) -> (rho, h)

    call eos(eos_input_tp, den_row, temp_row, &
             npts, nspec, &
             xn_zone, &
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

    real(kind=dp_t) :: temp, t0
    real(kind=dp_t) :: x0, y0, z0, x1, y1, z1, x2, y2, z2
    integer :: i, j, k
    real(kind=dp_t) :: r0, r1, r2

    t0 = s0(temp_comp)

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
    temp_row(1) = temp
    p_row(1) = p0
    den_row(1) = s0(rho_comp)
    xn_zone(:) = s0(spec_comp:spec_comp+nspec-1)/s0(rho_comp)

    ! (t, p) -> (rho, h)

    call eos(eos_input_tp, den_row, temp_row, &
             npts, nspec, &
             xn_zone, &
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
