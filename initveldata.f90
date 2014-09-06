module init_vel_module

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_physbc_module
  use define_bc_module
  use multifab_module
  use eos_module
  use variables
  use network
  use geometry, only: nr, spherical
  use ml_layout_module
  use ml_cc_restriction_module
  use multifab_fill_ghost_module

  implicit none

  real(dp_t), save :: pert_height

  private
  public :: initveldata

contains

  subroutine initveldata(u,s0_init,p0_init,dx,bc,mla)
    
    use probin_module, only: prob_lo, prob_hi, num_vortices
    use mt19937_module


    type(multifab) , intent(inout) :: u(:)
    real(kind=dp_t), intent(in   ) :: s0_init(:,0:,:)
    real(kind=dp_t), intent(in   ) :: p0_init(:,0:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: bc(:)
    type(ml_layout), intent(inout) :: mla

    real(kind=dp_t), pointer:: uop(:,:,:,:)
    integer :: lo(mla%dim),hi(mla%dim),ng,dm,nlevs
    integer :: i,j,k,n

    real(kind=dp_t) :: xloc_vortices(num_vortices)
    real(kind=dp_t) :: offset

    ! random numbers between -1 and 1                                                     
    real(kind=dp_t) :: alpha(3,3,3), beta(3,3,3), gamma(3,3,3)

    ! random numbers between 0 and 2*pi                                                   
    real(kind=dp_t) :: phix(3,3,3), phiy(3,3,3), phiz(3,3,3)

    ! L2 norm of k                                                                        
    real(kind=dp_t) :: normk(3,3,3)

    ! random number                                                                       
    real(kind=dp_t) :: rand


    dm = mla%dim
    nlevs = mla%nlevel

    if (dm == 2) then
       ! for now, this is calculated even if we don't use the velocity field in
       ! the initveldata_2d routine below
       offset = (prob_hi(1) - prob_lo(1)) / (num_vortices + 1)

       do i = 1, num_vortices
          xloc_vortices(i) = dble(i) * offset
       enddo

    else
       ! for 3-d we'll do some random velocity field stuff
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
    end if


    ng = nghost(u(1))

    do n=1,nlevs
       do i = 1, nfabs(u(n))
          uop => dataptr(u(n),i)
          lo =  lwb(get_box(u(n),i))
          hi =  upb(get_box(u(n),i))
          select case (dm)
          case (2)
             call initveldata_2d(uop(:,:,1,:), lo, hi, ng, dx(n,:), &
                                 s0_init(n,:,:), p0_init(n,:), xloc_vortices)
          ! 3d doesn't currently have vortice information coded !
          case (3) 
             call initveldata_3d(uop(:,:,:,:), lo, hi, ng, dx(n,:), &
                                 s0_init(n,:,:), p0_init(n,:), &
                                 alpha, beta, gamma, phix, phiy, phiz, normk)
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
    
       ! the loop over nlevs must count backwards to make sure the finer grids
       ! are done first
       do n=nlevs,2,-1

          ! set level n-1 data to be the average of the level n data covering 
          ! it
          call ml_cc_restriction(u(n-1),u(n),mla%mba%rr(n-1,:))
          
          ! fill level n ghost cells using interpolation from level n-1 data
          ! note that multifab_fill_boundary and multifab_physbc are called for
          ! both levels n-1 and n
          call multifab_fill_ghost_cells(u(n),u(n-1),ng,mla%mba%rr(n-1,:), &
                                         bc(n-1),bc(n),1,1,dm)

       enddo

    end if

  end subroutine initveldata


  subroutine initveldata_2d(u,lo,hi,ng,dx,s0_init,p0_init,xloc_vortices)

    use probin_module, only: prob_lo, apply_vel_field, velpert_scale, &
                             velpert_amplitude, velpert_height_loc

    integer           , intent(in   ) :: lo(:),hi(:),ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)
    real(kind=dp_t)   , intent(in   ) :: xloc_vortices(:)

    ! local variables
    real(kind=dp_t) :: xloc(2), upert(2)
    integer :: i, j, vortex
    real(kind=dp_t) :: xdist, ydist, r


    u = ZERO

    if (apply_vel_field) then

       do j = lo(2), hi(2)

          xloc(2) = prob_lo(2) + (dble(j)+HALF)*dx(2)

          ydist = xloc(2) - velpert_height_loc
          
          do i = lo(1), hi(1)

             upert = ZERO

             xloc(1) = prob_lo(1) + (dble(i)+HALF)*dx(1)

             ! loop over each vortex
             do vortex = 1, size(xloc_vortices, dim=1)

                xdist = xloc(1) - xloc_vortices(vortex)

                r = xdist**2 + ydist**2
                r = sqrt(r)

                ! e.g. Calder et al. ApJSS 143, 201-229 (2002)
                ! we set things up so that every other vortex has the same
                ! orientation
                upert(1) = upert(1) - ydist * &
                     velpert_amplitude * exp( -r**2/(TWO*velpert_scale)) &
                     * (-ONE)**vortex

                upert(2) = upert(2) + xdist * &
                     velpert_amplitude * exp(-r**2/(TWO*velpert_scale)) &
                     * (-ONE)**vortex
             enddo

             u(i,j,:) = u(i,j,:) + upert(:)

          enddo

       enddo
       
                
    endif

  end subroutine initveldata_2d


  subroutine initveldata_3d(u,lo,hi,ng,dx,s0_init,p0_init, &
                            alpha, beta, gamma, phix, phiy, phiz, normk)

    use probin_module, only: apply_vel_field, prob_lo, velpert_height_loc, &
                             velpert_steep, velpert_scale, velpert_amplitude

    integer           , intent(in   ) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(in   ) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(in   ) :: p0_init(0:)

    ! random numbers between -1 and 1                                                     
    real(kind=dp_t), intent(in) :: alpha(3,3,3), beta(3,3,3), gamma(3,3,3)

    ! random numbers between 0 and 2*pi                                                   
    real(kind=dp_t), intent(in) :: phix(3,3,3), phiy(3,3,3), phiz(3,3,3)

    ! L2 norm of k                                                                        
    real(kind=dp_t), intent(in) :: normk(3,3,3)

    ! Local variables                                                                     
    integer :: i, j, k
    integer :: iloc, jloc, kloc

    ! cos and sin of (2*pi*kx/L + phix), etc                                              
    real(kind=dp_t) :: cx(3,3,3), cy(3,3,3), cz(3,3,3)
    real(kind=dp_t) :: sx(3,3,3), sy(3,3,3), sz(3,3,3)

    ! the point we're at                                                                  
    real(kind=dp_t) :: xloc(3)

    ! perturbational velocity to add                                                      
    real(kind=dp_t) :: upert(3)


    ! initial the velocity
    u = ZERO

    ! this is the same velocity field that was used in the first sub-Chandra
    ! paper (Zingale et al. 2013)
    if (apply_vel_field) then

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
                           ( -beta(i,j,k)*dble(i)*cy(i,j,k)*cz(i,j,k)*sx(i,j,k) &
                            +alpha(i,j,k)*dble(j)*cx(i,j,k)*cz(i,j,k)*sy(i,j,k)) &
                            / normk(i,j,k)
                   enddo
                enddo
             enddo

             ! apply the cutoff function to the perturbational velocity                             
             do i=1,3
                upert(i) = velpert_amplitude*upert(i) * &
                     HALF*(ONE + tanh((velpert_height_loc + velpert_scale - velpert_steep - xloc(3))/velpert_steep))* &
                     HALF*(ONE + tanh((xloc(3) - velpert_height_loc - velpert_steep)/velpert_steep))
             enddo

             ! add perturbational velocity to background velocity                                   
             do i=1,3
                u(iloc,jloc,kloc,i) = u(iloc,jloc,kloc,i) + upert(i)
             enddo

          enddo
       enddo
    enddo

       

    endif

    
  end subroutine initveldata_3d

end module init_vel_module
