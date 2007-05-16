module advect_base_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use mkflux_module
  use eos_module
  use variables
  use geometry
  use make_grav_module
  use make_div_coeff_module

  implicit none

contains

   subroutine advect_base(vel,Sbar_in,p0_old,p0_new, &
                          s0_old,s0_new,temp0, &
                          gam1,div_coeff, &
                          dz,dt,anelastic_cutoff)

      real(kind=dp_t), intent(in   ) :: vel(:)
      real(kind=dp_t), intent(in   ) :: Sbar_in(:)
      real(kind=dp_t), intent(in   ) :: p0_old(:), s0_old(:,:)
      real(kind=dp_t), intent(  out) :: p0_new(:), s0_new(:,:)
      real(kind=dp_t), intent(inout) :: temp0(:),gam1(:)
      real(kind=dp_t), intent(in   ) :: div_coeff(:)
      real(kind=dp_t), intent(in   ) :: dz,dt,anelastic_cutoff

      integer :: i

      if (spherical .eq. 0) then

        call advect_base_state_planar(vel,p0_old,p0_new,s0_old,s0_new,temp0, &
                                      gam1,dz,dt,anelastic_cutoff)

      else

        call advect_base_state_spherical(vel,Sbar_in,p0_old,p0_new,s0_old,s0_new,temp0, &
                                         gam1,div_coeff,&
                                         dt,anelastic_cutoff)
      end if


   end subroutine advect_base

   subroutine advect_base_state_planar (vel,p0_old,p0_new,s0_old,s0_new,temp0, &
                                        gam1,dz,dt,anelastic_cutoff)

      real(kind=dp_t), intent(in   ) :: vel(:)
      real(kind=dp_t), intent(in   ) :: p0_old(:), s0_old(:,:)
      real(kind=dp_t), intent(  out) :: p0_new(:), s0_new(:,:)
      real(kind=dp_t), intent(inout) :: temp0(:)
      real(kind=dp_t), intent(inout) :: gam1(:)
      real(kind=dp_t), intent(in   ) :: dz,dt,anelastic_cutoff

!     Local variables
      integer :: i, j, n, nz
      real(kind=dp_t) :: coeff

      real (kind = dp_t), allocatable :: H(:,:)
      real (kind = dp_t), allocatable :: force(:)
      real (kind = dp_t), allocatable :: edge(:)

      nz = size(p0_new,dim=1)

      allocate(    force(nz))
      allocate(     edge(nz+1))

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE P0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      force = ZERO
      call mkflux_1d(p0_old,edge,vel,force,1,dz,dt)
      do j = 1,nz
        p0_new(j) = p0_old(j) - &
             dt / dz * HALF * (vel(j) + vel(j+1)) * (edge(j+1) - edge(j))
      end do

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE RHOX0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do n = spec_comp,spec_comp+nspec-1
         do j = 1,nz
            force(j) = -s0_old(j,n) * (vel(j+1) - vel(j)) / dz
         end do

         call mkflux_1d(s0_old(:,n),edge,vel,force,1,dz,dt)

         do j = 1,nz
            s0_new(j,n) = s0_old(j,n) - &
                 dt / dz * (edge(j+1) * vel(j+1) - edge(j) * vel(j))
         end do

      enddo

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE RHO0 FROM RHOX0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j = 1,nz
        s0_new(j,rho_comp) =  zero
        do n = spec_comp,spec_comp+nspec-1
          s0_new(j,rho_comp) =  s0_new(j,rho_comp) + s0_new(j,n)
        end do
      end do


!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE RHOH0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j = 1,nz
         force(j) = -s0_old(j,rhoh_comp) * (vel(j+1) - vel(j)) / dz
      end do

      call mkflux_1d(s0_old(:,rhoh_comp),edge,vel,force,1,dz,dt)

      do j = 1,nz
         s0_new(j,rhoh_comp) = s0_old(j,rhoh_comp) - &
              dt / dz * (edge(j+1) * vel(j+1) - edge(j) * vel(j))
      end do


!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     MAKE TEMP0 AND GAM1 FROM P0 AND RHO0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j = 1,nz

         den_row(1)  = s0_new(j,rho_comp)
         temp_row(1) = temp0(j)
         p_row(1)    = p0_new(j)
         xn_zone(1:) = s0_new(j,spec_comp:)/s0_new(j,rho_comp)

         ! (rho,P) --> T, h
         input_flag = 4

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

         temp0(j) = temp_row(1)
         gam1(j) = gam1_row(1)


      end do

      deallocate(force,edge)

   end subroutine advect_base_state_planar

   subroutine advect_base_state_spherical (vel,Sbar_in,p0_old,p0_new,s0_old,s0_new,temp0, &
                                           gam1,div_coeff_old,& 
                                           dt,anelastic_cutoff)

      real(kind=dp_t), intent(in   ) :: vel(:),Sbar_in(:)
      real(kind=dp_t), intent(in   ) :: p0_old(:), s0_old(:,:)
      real(kind=dp_t), intent(  out) :: p0_new(:), s0_new(:,:)
      real(kind=dp_t), intent(inout) :: temp0(:), gam1(:)
      real(kind=dp_t), intent(in   ) :: div_coeff_old(:)
      real(kind=dp_t), intent(in   ) :: dt,anelastic_cutoff

!     Local variables
      integer :: i, j, k, n, nz
      real(kind=dp_t) :: dtdr,divbetaw,betahalf,factor,integral
      real(kind=dp_t) :: div_w0

      real (kind = dp_t), allocatable :: m(:)
      real (kind = dp_t), allocatable :: force(:)
      real (kind = dp_t), allocatable :: eta(:)
      real (kind = dp_t), allocatable :: edge(:)
      real (kind = dp_t), allocatable :: div_coeff_new(:)
      real (kind = dp_t), allocatable :: beta(:),beta_new(:),beta_nh(:)
      real (kind = dp_t), allocatable :: gam1_old(:)
      real (kind = dp_t), allocatable :: grav_cell(:)

      dtdr = dt / dr

      nz = size(p0_new,dim=1)

      ! Cell-centered
      allocate(force(nz))
      allocate(eta(nz))
      allocate(m(nz))
      allocate(gam1_old(nz))
      allocate(grav_cell(nz))
      allocate(div_coeff_new(nz))

      ! Edge-centered
      allocate(edge(nz+1))
      allocate(beta(nz+1),beta_new(nz+1),beta_nh(nz+1))

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE RHOX0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do n = spec_comp,spec_comp+nspec-1

         ! compute the force -- include the geometric source term that
         ! results from expanding out the spherical divergence
         do j = 1,nz
            force(j) = -s0_old(j,n) * (vel(j+1) - vel(j)) / dr - &
                       2.0_dp_t*s0_old(j,n)*HALF*(vel(j) + vel(j+1))/z(j)
         end do

         call mkflux_1d(s0_old(:,n),edge,vel,force,1,dr,dt)

         do j = 1,nz
            s0_new(j,n) = s0_old(j,n) - &
                 dtdr / z(j)**2 * ( zl(j+1)**2 * edge(j+1) * vel(j+1) &
                                   -zl(j  )**2 * edge(j  ) * vel(j  ))
         end do

      enddo

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE RHO0 FROM RHOX0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j = 1,nz
        s0_new(j,rho_comp) =  zero
        do n = spec_comp,spec_comp+nspec-1
          s0_new(j,rho_comp) =  s0_new(j,rho_comp) + s0_new(j,n)
        end do
      end do

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE P0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Put beta_old on edges
      call put_1d_beta_on_edges(div_coeff_old,beta)
 
      ! Update p0 -- predictor
      do j = 1,nz
         divbetaw = one / (z(j)**2) * (zl(j+1)**2 * beta(j+1) * vel(j+1) - &
                                       zl(j  )**2 * beta(j  ) * vel(j  ) ) / dr
         betahalf = half * (beta(j+1) + beta(j))
         factor = half * dt * gam1(j) * (Sbar_in(j) - divbetaw / betahalf)
         p0_new(j) = p0_old(j) * (one + factor ) / (one - factor)
 
      end do
 
      do j = 1,nz
         ! (rho, p) --> T,h, etc
         input_flag = 4
         den_row(1)  = s0_new(j,rho_comp)
           p_row(1) =  p0_new(j)
         gam1_old(j) = gam1(j)
 
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
         gam1(j) = gam1_row(1)
      end do
 
      call make_grav_cell(grav_cell,s0_new(:,rho_comp))
 
      ! Define beta^n+1 at cell edges using the new gravity above
      call make_div_coeff(div_coeff_new,s0_new(:,rho_comp),p0_new,gam1,grav_cell,anelastic_cutoff)
      call put_1d_beta_on_edges(div_coeff_new,beta_new)

      ! time-centered beta
      beta_nh = HALF*(beta + beta_new)
 
      ! Update p0 -- corrector
      do j = 1,nz
         divbetaw = one / (z(j)**2) * (zl(j+1)**2 * beta_nh(j+1) * vel(j+1) - &
                                       zl(j  )**2 * beta_nh(j  ) * vel(j  ) ) / dr
         betahalf = half * (beta_nh(j+1) + beta_nh(j))
         factor = half * dt * (Sbar_in(j) - divbetaw / betahalf)
         p0_new(j) = p0_old(j) * (one + factor * gam1_old(j)) / (one - factor * gam1(j))
 
      end do

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE RHOH0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j = 1,nz

         div_w0 = (vel(j+1) - vel(j)) / dr 

         force(j) = -s0_old(j,rhoh_comp) * div_w0 - &
              2.0_dp_t*s0_old(j,rhoh_comp)*HALF*(vel(j) + vel(j+1))/z(j)

         eta(j) = gam1_old(j) * p0_old(j) * (Sbar_in(j) - div_w0)

      end do

      call mkflux_1d(s0_old(:,rhoh_comp),edge,vel,force,1,dr,dt)

      do j = 1,nz

         s0_new(j,rhoh_comp) = s0_old(j,rhoh_comp) - &
              dtdr / z(j)**2 * ( zl(j+1)**2 * edge(j+1) * vel(j+1) &
                                -zl(j  )**2 * edge(j  ) * vel(j  ))

         s0_new(j,rhoh_comp) = s0_new(j,rhoh_comp) + dt * eta(j)

      end do


!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     MAKE TEMP0 AND GAM1 FROM P0 AND RHO0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j = 1,nz

         den_row(1)  = s0_new(j,rho_comp)
         temp_row(1) = temp0(j)
         p_row(1)    = p0_new(j)
         xn_zone(1:) = s0_new(j,spec_comp:)/s0_new(j,rho_comp)

         ! (rho,P) --> T, h
         input_flag = 4

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

         temp0(j) = temp_row(1)
         gam1(j) = gam1_row(1)

      end do

      deallocate(force,eta,edge,beta,beta_new,beta_nh,div_coeff_new,gam1_old,grav_cell)

   end subroutine advect_base_state_spherical

end module advect_base_module
