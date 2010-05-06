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
  public :: initscalardata, initscalardata_on_level, initveldata

contains

  subroutine initscalardata(s,s0_init,p0_init,dx,bc,mla)

    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(dm),hi(dm),ng
    integer :: i,n
    
    ng = s(1)%ng

    do n=1,nlevs
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n),i) ) cycle
          sop => dataptr(s(n),i)
          lo =  lwb(get_box(s(n),i))
          hi =  upb(get_box(s(n),i))
          select case (dm)
          case (2)
             call initscalardata_2d(sop(:,:,1,:), lo, hi, ng, dx(n,:), s0_init(n,:,:), &
                                    p0_init(n,:))
          case (3)
             if (spherical .eq. 1) then
                call initscalardata_3d_sphr(sop(:,:,:,:), lo, hi, ng, dx(n,:), &
                                            s0_init(1,:,:), p0_init(1,:))
             else
                call initscalardata_3d(sop(:,:,:,:), lo, hi, ng, dx(n,:), s0_init(n,:,:), &
                                       p0_init(n,:))
             end if
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
                                         bc(n-1),bc(n),rho_comp,dm+rho_comp,nscal, &
                                         fill_crse_input=.false.)

       enddo

    end if

  end subroutine initscalardata

  subroutine initscalardata_on_level(n,s,s0_init,p0_init,dx,bc)

    integer        , intent(in   ) :: n
    type(multifab) , intent(inout) :: s
    real(kind=dp_t), intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_level) , intent(in   ) :: bc

    ! local
    integer                  :: ng,i
    integer                  :: lo(dm),hi(dm)
    real(kind=dp_t), pointer :: sop(:,:,:,:)

    ng = s%ng

    do i = 1, s%nboxes
       if ( multifab_remote(s,i) ) cycle
       sop => dataptr(s,i)
       lo =  lwb(get_box(s,i))
       hi =  upb(get_box(s,i))
       select case (dm)
       case (2)
          call initscalardata_2d(sop(:,:,1,:),lo,hi,ng,dx,s0_init,p0_init)
       case (3)
          if (spherical .eq. 1) then
                call initscalardata_3d_sphr(sop(:,:,:,:),lo,hi,ng,dx,s0_init,p0_init)
             else
                call initscalardata_3d(sop(:,:,:,:),lo,hi,ng,dx,s0_init,p0_init)
             end if
       end select
    end do

    call multifab_fill_boundary(s)

    call multifab_physbc(s,rho_comp,dm+rho_comp,nscal,bc)

  end subroutine initscalardata_on_level

  subroutine initscalardata_2d(s,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: prob_lo, perturb_model

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

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
          s(i,j,rho_comp)  = s0_init(j,rho_comp)
          s(i,j,rhoh_comp) = s0_init(j,rhoh_comp)
          s(i,j,temp_comp) = s0_init(j,temp_comp)
          s(i,j,spec_comp:spec_comp+nspec-1) = &
               s0_init(j,spec_comp:spec_comp+nspec-1)
          s(i,j,trac_comp:trac_comp+ntrac-1) = &
               s0_init(j,trac_comp:trac_comp+ntrac-1)
       enddo
    enddo
    
    ! add an optional perturbation
    if (perturb_model) then
       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j)+HALF) * dx(2)
          
          do i = lo(1), hi(1)
             x = prob_lo(1) + (dble(i)+HALF) * dx(1)
          
             call perturb_2d(x, y, p0_init(j), s0_init(j,:), &
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

  subroutine initscalardata_3d(s,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: prob_lo, perturb_model
    
    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    !     Local variables
    integer         :: i,j,k
    real(kind=dp_t) :: x,y,z
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)

    ! initial the domain with the base state
    s = ZERO
    
    ! initialize the scalars
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             s(i,j,k,rho_comp)  = s0_init(k,rho_comp)
             s(i,j,k,rhoh_comp) = s0_init(k,rhoh_comp)
             s(i,j,k,temp_comp) = s0_init(k,temp_comp)
             s(i,j,k,spec_comp:spec_comp+nspec-1) = &
                  s0_init(k,spec_comp:spec_comp+nspec-1)
             s(i,j,k,trac_comp:trac_comp+ntrac-1) = &
                  s0_init(k,trac_comp:trac_comp+ntrac-1)
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

                call perturb_3d(x, y, z, p0_init(k), s0_init(k,:), &
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

  end subroutine initscalardata_3d

  subroutine initscalardata_3d_sphr(s,lo,hi,ng,dx,s0_init,p0_init)

    use probin_module, only: prob_lo, perturb_model
    
    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    !     Local variables
    integer         :: i,j,k,comp
    real(kind=dp_t) :: x,y,z
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)
    real(kind=dp_t), allocatable :: p0_cart(:,:,:,:)

    ! initial the domain with the base state
    s = ZERO
    
    ! if we are spherical, we want to make sure that p0 is good, since that is
    ! what is needed for HSE.  Therefore, we will put p0 onto a cart array and
    ! then initialize h from rho, X, and p0.
    allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    ! initialize the scalars
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,s0_init(:,rho_comp), &
                                      s(:,:,:,rho_comp:),lo,hi,dx,ng)

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,s0_init(:,temp_comp), &
                                      s(:,:,:,temp_comp:),lo,hi,dx,ng)

    ! initialize p0_cart
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,p0_init(:), &
                                      p0_cart(:,:,:,1:),lo,hi,dx,0)

    ! initialize species
    do comp = spec_comp, spec_comp+nspec-1
       call put_1d_array_on_cart_3d_sphr(.false.,.false.,s0_init(:,comp), &
                                         s(:,:,:,comp:),lo,hi,dx,ng)
    end do

    ! initialize tracers
    do comp = trac_comp, trac_comp+ntrac-1
       call put_1d_array_on_cart_3d_sphr(.false.,.false.,s0_init(:,comp), &
                                         s(:,:,:,comp:),lo,hi,dx,ng)
    end do

    ! initialize (rho h) using the EOS
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             temp_eos(1) = s(i,j,k,temp_comp)
             p_eos(1) = p0_cart(i,j,k,1)
             den_eos(1) = s(i,j,k,rho_comp)
             xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

             call eos(eos_input_rp, den_eos, temp_eos, &
                      npts, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      .false.)

             s(i,j,k,rhoh_comp) = den_eos(1)*h_eos(1)
             s(i,j,k,temp_comp) = temp_eos(1)

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

                call perturb_3d_sphr(x, y, z, p0_cart(i,j,k,1), s(i,j,k,:), &
                                     dens_pert, rhoh_pert, rhoX_pert, temp_pert, trac_pert)

                s(i,j,k,rho_comp) = dens_pert
                s(i,j,k,rhoh_comp) = rhoh_pert
                s(i,j,k,temp_comp) = temp_pert
                s(i,j,k,spec_comp:spec_comp+nspec-1) = rhoX_pert(:)
                s(i,j,k,trac_comp:trac_comp+ntrac-1) = trac_pert(:)
             enddo
          enddo
       enddo

    end if

    deallocate(p0_cart)

  end subroutine initscalardata_3d_sphr

  subroutine initveldata(u,s0_init,p0_init,dx,bc,mla)

    use geometry, only: nlevs

    type(multifab) , intent(inout) :: u(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: uop(:,:,:,:)
    integer :: lo(dm),hi(dm),ng
    integer :: i,n
    
    ng = u(1)%ng

    do n=1,nlevs
       do i = 1, u(n)%nboxes
          if ( multifab_remote(u(n),i) ) cycle
          uop => dataptr(u(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call initveldata_2d(uop(:,:,1,:), lo, hi, ng, dx(n,:), &
                                 s0_init(n,:,:), p0_init(n,:))
          case (3) 
             if (spherical .eq. 1) then
                call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx(n,:), &
                                    s0_init(1,:,:), p0_init(1,:))
             else
                call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx(n,:), &
                                    s0_init(n,:,:), p0_init(n,:))
             end if
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
                                         bc(n-1),bc(n),1,1,dm,fill_crse_input=.false.)
       enddo
       
    end if

  end subroutine initveldata

  subroutine initveldata_2d(u,lo,hi,ng,dx,s0_init,p0_init)

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    ! Local variables

    ! initial the velocity
    u = ZERO

  end subroutine initveldata_2d

  subroutine initveldata_3d(u,lo,hi,ng,dx,s0_init,p0_init)

    integer           , intent(in   ) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    ! Local variables

    ! initial the velocity
    u = ZERO
    
  end subroutine initveldata_3d

  subroutine perturb_2d(x, y, p0_init, s0_init, dens_pert, rhoh_pert, rhoX_pert, &
                        temp_pert, trac_pert)

    use probin_module, only: pert_temp_factor, pert_rad_factor

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp,t0
    real(kind=dp_t) :: x1, y1, r1, x2, y2, r2, x3, y3, r3

    t0 = s0_init(temp_comp)

    x1 = 5.0d7
    y1 = 6.5d7
    r1 = pert_rad_factor*sqrt( (x-x1)**2 +(y-y1)**2 ) / 2.5d6
    
    x2 = 1.2d8
    y2 = 8.5d7
    r2 = pert_rad_factor*sqrt( (x-x2)**2 +(y-y2)**2 ) / 2.5d6
    
    x3 = 2.0d8
    y3 = 7.5d7
    r3 = pert_rad_factor*sqrt( (x-x3)**2 +(y-y3)**2 ) / 2.5d6

    temp = t0 * (1.d0 + pert_temp_factor* &
                 (0.150d0 * (1.d0 + tanh(2.d0-r1)) + &
                  0.300d0 * (1.d0 + tanh(2.d0-r2)) + &
                  0.225d0 * (1.d0 + tanh(2.d0-r3))))
          
    ! Use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    temp_eos(1) = temp
    p_eos(1) = p0_init
    den_eos(1) = s0_init(rho_comp)
    xn_eos(1,:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_tp, den_eos, temp_eos, &
             npts, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             .false.)

    dens_pert = den_eos(1)
    rhoh_pert = den_eos(1)*h_eos(1)
    rhoX_pert = dens_pert*xn_eos(1,:)

    temp_pert = temp

    trac_pert = ZERO

  end subroutine perturb_2d

  subroutine perturb_3d(x, y, z, p0_init, s0_init, dens_pert, rhoh_pert, &
                        rhoX_pert, temp_pert, trac_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles

    real(kind=dp_t), intent(in ) :: x, y, z
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp, t0
    real(kind=dp_t) :: x0, y0, z0, r0

    t0 = s0_init(temp_comp)

    x0 = 3.6d7
    y0 = 3.6d7
    z0 = 8.5d7

    r0 = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 ) / 2.5d6

    temp = t0 * (ONE + TWO * (0.15d0 * (1.d0 + tanh((2.d0-r0)))))

    ! Use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    temp_eos(1) = temp
    p_eos(1) = p0_init
    den_eos(1) = s0_init(rho_comp)
    xn_eos(1,:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_tp, den_eos, temp_eos, &
             npts, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             .false.)

    dens_pert = den_eos(1)
    rhoh_pert = den_eos(1)*h_eos(1)
    rhoX_pert = dens_pert*xn_eos(1,:)

    temp_pert = temp

    trac_pert = ZERO

  end subroutine perturb_3d

  subroutine perturb_3d_sphr(x, y, z, p0_init, s0_init, dens_pert, rhoh_pert, &
                             rhoX_pert, temp_pert, trac_pert)

    ! apply an optional perturbation to the initial temperature field
    ! to see some bubbles
    real(kind=dp_t), intent(in ) :: x, y, z
    real(kind=dp_t), intent(in ) :: p0_init, s0_init(:)
    real(kind=dp_t), intent(out) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t), intent(out) :: rhoX_pert(:)
    real(kind=dp_t), intent(out) :: trac_pert(:)

    real(kind=dp_t) :: temp, t0
    real(kind=dp_t) :: x0, y0, z0, r0

    t0 = s0_init(temp_comp)

    ! center of star is at 2.5d8
    x0 = 2.5d8
    y0 = 2.5d8
    z0 = 3.0d8

    ! Tanh bubbles
    r0 = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 ) / 2.5d6
    
    ! This case works
    temp = t0 * (ONE + TWO*(.150_dp_t * 0.5_dp_t * (1.0_dp_t + tanh((2.0_dp_t-r0)))))

    ! Use the EOS to make this temperature perturbation occur at constant 
    ! pressure
    temp_eos(1) = temp
    p_eos(1) = p0_init
    den_eos(1) = s0_init(rho_comp)
    xn_eos(1,:) = s0_init(spec_comp:spec_comp+nspec-1)/s0_init(rho_comp)

    call eos(eos_input_tp, den_eos, temp_eos, &
             npts, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             .false.)

    dens_pert = den_eos(1)
    rhoh_pert = den_eos(1)*h_eos(1)
    rhoX_pert = dens_pert*xn_eos(1,:)

    temp_pert = temp

    trac_pert = ZERO

  end subroutine perturb_3d_sphr

end module init_module
