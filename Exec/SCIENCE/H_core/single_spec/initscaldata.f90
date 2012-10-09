module init_scalar_module

  use multifab_module
  use ml_layout_module
  use bl_constants_module
  use multifab_physbc_module
  use ml_restriction_module
  use multifab_fill_ghost_module
  use eos_module
  use variables
  use network
  use fill_3d_module, only: put_1d_array_on_cart_3d_sphr

  implicit none

  private
  public :: initscalardata, initscalardata_on_level

contains

  subroutine initscalardata(s,s0_init,p0_init,dx,bc,mla)

    use geometry, only: spherical

    type(multifab) , intent(inout) :: s(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: sop(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim),ng
    integer :: i,n,dm,nlevs
    
    dm = mla%dim
    nlevs = mla%nlevel

    ng = nghost(s(1))

    do n=1,nlevs
       do i = 1, nboxes(s(n))
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

    use geometry, only: spherical

    integer        , intent(in   ) :: n
    type(multifab) , intent(inout) :: s
    real(kind=dp_t), intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_level) , intent(in   ) :: bc

    ! local
    integer                  :: ng,i,dm
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    real(kind=dp_t), pointer :: sop(:,:,:,:)

    dm = get_dim(s)

    ng = nghost(s)

    do i = 1, nboxes(s)
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
    use init_perturb_module

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
    use init_perturb_module
    
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
    use init_perturb_module
    use geometry, only: center
    use probin_module, only: prob_lo, perturb_model, velpert_amplitude, &
         velpert_radius, velpert_steep, velpert_scale
    use mt19937_module


    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(inout) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    !     Local variables
    integer         :: i,j,k,comp
    real(kind=dp_t) :: x,y,z, r
    real(kind=dp_t) :: dens_pert, rhoh_pert, temp_pert
    real(kind=dp_t) :: rhoX_pert(nspec), trac_pert(ntrac)
    real(kind=dp_t), allocatable :: p0_cart(:,:,:,:)

    ! random numbers between -1 and 1
    real(kind=dp_t) :: alpha(3,3,3), beta(3,3,3), gamma(3,3,3)

    ! random numbers between 0 and 2*pi
    real(kind=dp_t) :: phix(3,3,3), phiy(3,3,3), phiz(3,3,3)

    ! L2 norm of k
    real(kind=dp_t) :: normk(3,3,3)

    ! random number
    real(kind=dp_t) :: rand
    
    integer :: ii, jj, kk

    ! cos and sin of (2*pi*kx/L + phix), etc
    real(kind=dp_t) :: cx(3,3,3), cy(3,3,3), cz(3,3,3)
    real(kind=dp_t) :: sx(3,3,3), sy(3,3,3), sz(3,3,3)

    real(kind=dp_t) :: theta,phi


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
       z = prob_lo(3) + (dble(k)+HALF) * dx(3)

       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j)+HALF) * dx(2)
          
          do i = lo(1), hi(1)
             x = prob_lo(1) + (dble(i)+HALF) * dx(1)

             temp_eos(1) = s(i,j,k,temp_comp)
             p_eos(1) = p0_cart(i,j,k,1)
             den_eos(1) = s(i,j,k,rho_comp)
             xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

             r = sqrt( (x -center(1))**2 + (y -center(2))**2 + &
                  (z -center(3))**2)

             call eos(r, eos_input_rp, den_eos, temp_eos, &
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
       if (.false.) then

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
          
       else
          ! random temperature fluctuations          

          ! load in random numbers alpha, beta, gamma, phix, phiy, and phiz
          call init_genrand(20908)
          do i=1,3
             do j=1,3
                do k=1,3
                   rand = genrand_real1()
                   rand = 2.0d0*rand - 1.0d0
                   alpha(i,j,k) = rand
                   rand = genrand_real1()
                   rand = 2.0d0*rand - 1.0d0
                   beta(i,j,k) = rand
                   rand = genrand_real1()
                   rand = 2.0d0*rand - 1.0d0
                   gamma(i,j,k) = rand
                   rand = genrand_real1()
                   rand = 2.0d0*M_PI*rand
                   phix(i,j,k) = rand
                   rand = genrand_real1()
                   rand = 2.0d0*M_PI*rand
                   phiy(i,j,k) = rand
                   rand = genrand_real1()
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
          
          do k = lo(3), hi(3)
             z = prob_lo(3) + (dble(k)+HALF) * dx(3)
             
             do j = lo(2), hi(2)
                y = prob_lo(2) + (dble(j)+HALF) * dx(2)
                
                do i = lo(1), hi(1)
                   x = prob_lo(1) + (dble(i)+HALF) * dx(1)

                   temp_pert = ZERO
                   
                   ! loop over the 27 combinations of fourier components
                   do ii=1,3
                      do jj=1,3
                         do kk=1,3
                            ! compute cosiines and siines
                            cx(ii,jj,kk) = cos(2.0d0*M_PI*dble(ii)*x/velpert_scale + phix(ii,jj,kk))
                            cy(ii,jj,kk) = cos(2.0d0*M_PI*dble(jj)*y/velpert_scale + phiy(ii,jj,kk))
                            cz(ii,jj,kk) = cos(2.0d0*M_PI*dble(kk)*z/velpert_scale + phiz(ii,jj,kk))
                            sx(ii,jj,kk) = sin(2.0d0*M_PI*dble(ii)*x/velpert_scale + phix(ii,jj,kk))
                            sy(ii,jj,kk) = sin(2.0d0*M_PI*dble(jj)*y/velpert_scale + phiy(ii,jj,kk))
                            sz(ii,jj,kk) = sin(2.0d0*M_PI*dble(kk)*z/velpert_scale + phiz(ii,jj,kk))
                         enddo
                      enddo
                   enddo

                   ! loop over the 27 combiinatiions of fouriier components
                   do ii=1,3
                      do jj=1,3
                         do kk=1,3
                            ! compute contriibutiion from perturbatiion velociity from each mode
                            temp_pert = temp_pert + &
                                 (alpha(ii,jj,kk)*dble(ii)*cy(ii,jj,kk)*cz(ii,jj,kk)*sx(ii,jj,kk) &
                                 -gamma(ii,jj,kk)*dble(jj)*cx(ii,jj,kk)*cz(ii,jj,kk)*sy(ii,jj,kk) &
                                 +beta(ii,jj,kk)*dble(kk)*cx(ii,jj,kk)*cy(ii,jj,kk)*sz(ii,jj,kk) ) &
                                 / normk(ii,jj,kk)
                         enddo
                      enddo
                   enddo
                   
                   r = sqrt( (x-center(1))**2 + (y-center(2))**2 + (z-center(3))**2 ) 

                   ! apply the cutoff function to the perturbational velocity
                   ! with 2D hack y is like radius
                   temp_pert = velpert_amplitude * temp_pert &
                        *(0.5d0+0.5d0*tanh((velpert_radius - r)/velpert_steep))
                   
                   ! add perturbational velocity to background velocity
                   temp_pert = temp_pert + s(i,j,k,temp_comp)

                   ! use the EOS to make this temperature perturbation occur at
                   ! constant  pressure -- this adjusts the base state density (rho0)
                   ! In calling function will make sure HSE and EOS are satisfied by: 
                   !   1. set rho0 to the the average of the 3D state s(rho_comp)
                   !   2. adjust p0 so that HSE is satisfied (dp0/dr = g0*rho0)
                   !   3. use 3D rho and p to recompute 3D T, h via the EOS
                   !   4. set T0 to the average of the 3D T, h0 to average of 3D h
                   temp_eos(1) = temp_pert
                   p_eos(1) = p0_cart(i,j,k,1)
                   den_eos(1) = s(i,j,k,rho_comp)
                   xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1) / &
                        s(i,j,k,rho_comp)

                   call eos(r, eos_input_tp, den_eos, temp_eos, &
                        npts, &
                        xn_eos, &
                        p_eos, h_eos, e_eos, &
                        cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                        dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                        dpdX_eos, dhdX_eos, &
                        gam1_eos, cs_eos, s_eos, &
                        dsdt_eos, dsdr_eos, &
                        .false.)
                   
                   s(i,j,k,rho_comp) = den_eos(1)
                   s(i,j,k,rhoh_comp) = den_eos(1)*h_eos(1)
                   s(i,j,k,temp_comp) = temp_pert
                   s(i,j,k,spec_comp:spec_comp+nspec-1) = den_eos(1)*xn_eos(1,:)
                   s(i,j,k,trac_comp:trac_comp+ntrac-1) = ZERO

                enddo
             enddo
          enddo

       endif
    end if

    deallocate(p0_cart)

  end subroutine initscalardata_3d_sphr

end module init_scalar_module
