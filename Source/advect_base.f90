module advect_base_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use mkflux_module
  use eos_module
  use variables
  use geometry
  use make_grav_module
  use cell_to_edge_module
  use make_div_coeff_module

  implicit none

contains

   subroutine advect_base(vel,Sbar_in,p0_old,p0_new, &
                          s0_old,s0_new,&
                          gam1,div_coeff, &
                          dz,dt,anelastic_cutoff)

      real(kind=dp_t), intent(in   ) :: vel(0:)
      real(kind=dp_t), intent(in   ) :: Sbar_in(0:)
      real(kind=dp_t), intent(in   ) :: p0_old(0:), s0_old(0:,:)
      real(kind=dp_t), intent(  out) :: p0_new(0:), s0_new(0:,:)
      real(kind=dp_t), intent(inout) :: gam1(0:)
      real(kind=dp_t), intent(in   ) :: div_coeff(0:)
      real(kind=dp_t), intent(in   ) :: dz,dt,anelastic_cutoff

      if (spherical .eq. 0) then

        call advect_base_state_planar(vel,p0_old,p0_new,s0_old,s0_new,&
                                      gam1,dz,dt)

      else

        call advect_base_state_spherical(vel,Sbar_in,p0_old,p0_new,s0_old,s0_new,&
                                         gam1,div_coeff,&
                                         dt,anelastic_cutoff)
      end if


   end subroutine advect_base

   subroutine advect_base_state_planar (vel,p0_old,p0_new,s0_old,s0_new,&
                                        gam1,dz,dt)

      real(kind=dp_t), intent(in   ) :: vel(0:)
      real(kind=dp_t), intent(in   ) :: p0_old(0:), s0_old(0:,:)
      real(kind=dp_t), intent(  out) :: p0_new(0:), s0_new(0:,:)
      real(kind=dp_t), intent(inout) :: gam1(0:)
      real(kind=dp_t), intent(in   ) :: dz,dt

!     Local variables
      integer :: j, n, nz

      real (kind = dp_t), allocatable :: force(:)
      real (kind = dp_t), allocatable :: edge(:)

      ! nz is the size of a cell-centered quantity
      nz = size(p0_new,dim=1)

      allocate(    force(0:nz-1))
      allocate(     edge(0:nz))

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE P0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      force = ZERO
      call mkflux_1d(p0_old,edge,vel,force,1,dz,dt)
      do j = 0,nz-1
        p0_new(j) = p0_old(j) - &
             dt / dz * HALF * (vel(j) + vel(j+1)) * (edge(j+1) - edge(j))
      end do

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE RHOX0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do n = spec_comp,spec_comp+nspec-1
         do j = 0,nz-1
            force(j) = -s0_old(j,n) * (vel(j+1) - vel(j)) / dz
         end do

         call mkflux_1d(s0_old(:,n),edge,vel,force,1,dz,dt)

         do j = 0,nz-1
            s0_new(j,n) = s0_old(j,n) - &
                 dt / dz * (edge(j+1) * vel(j+1) - edge(j) * vel(j))
         end do

      enddo

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE RHO0 FROM RHOX0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j = 0,nz-1
        s0_new(j,rho_comp) =  s0_old(j,rho_comp)
        do n = spec_comp,spec_comp+nspec-1
          s0_new(j,rho_comp) =  s0_new(j,rho_comp) + (s0_new(j,n)-s0_old(j,n))
        end do
      end do


!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE RHOH0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j = 0,nz-1
         force(j) = -s0_old(j,rhoh_comp) * (vel(j+1) - vel(j)) / dz
      end do

      call mkflux_1d(s0_old(:,rhoh_comp),edge,vel,force,1,dz,dt)

      do j = 0,nz-1
         s0_new(j,rhoh_comp) = s0_old(j,rhoh_comp) - &
              dt / dz * (edge(j+1) * vel(j+1) - edge(j) * vel(j))
      end do


!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     MAKE TEMP0 AND GAM1 FROM P0 AND RHO0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j = 0,nz-1

         den_row(1)  = s0_new(j,rho_comp)
         temp_row(1) = s0_old(j,temp_comp)
         p_row(1)    = p0_new(j)
         xn_row(1,:) = s0_new(j,spec_comp:spec_comp+nspec-1)/s0_new(j,rho_comp)

         ! (rho,P) --> T, h
         call eos(eos_input_rp, den_row, temp_row, &
                  npts, nspec, &
                  xn_row, &
                  p_row, h_row, e_row, &
                  cv_row, cp_row, xne_row, eta_row, pele_row, &
                  dpdt_row, dpdr_row, dedt_row, dedr_row, &
                  dpdX_row, dhdX_row, &
                  gam1_row, cs_row, s_row, &
                  dsdt_row, dsdr_row, &
                  do_diag)

         s0_new(j,temp_comp) = temp_row(1)
         gam1(j) = gam1_row(1)

      end do

      deallocate(force,edge)

   end subroutine advect_base_state_planar

   subroutine advect_base_state_spherical (vel,Sbar_in,p0_old,p0_new,s0_old,s0_new,&
                                           gam1,div_coeff_old,& 
                                           dt,anelastic_cutoff)

      real(kind=dp_t), intent(in   ) :: vel(0:),Sbar_in(0:)
      real(kind=dp_t), intent(in   ) :: p0_old(0:), s0_old(0:,:)
      real(kind=dp_t), intent(  out) :: p0_new(0:), s0_new(0:,:)
      real(kind=dp_t), intent(inout) :: gam1(0:)
      real(kind=dp_t), intent(in   ) :: div_coeff_old(0:)
      real(kind=dp_t), intent(in   ) :: dt,anelastic_cutoff

!     Local variables
      integer :: j, n, nz
      real(kind=dp_t) :: dtdr,divbetaw,betahalf,factor
      real(kind=dp_t) :: div_w0

      real (kind = dp_t), allocatable :: force(:)
      real (kind = dp_t), allocatable :: eta(:)
      real (kind = dp_t), allocatable :: edge(:)
      real (kind = dp_t), allocatable :: div_coeff_new(:)
      real (kind = dp_t), allocatable :: beta(:),beta_new(:),beta_nh(:)
      real (kind = dp_t), allocatable :: gam1_old(:)
      real (kind = dp_t), allocatable :: grav_cell(:)

      dtdr = dt / dr

      ! nz is the size of a cell-centered quantity
      nz = size(p0_new,dim=1)

      ! Cell-centered
      allocate(force(0:nz-1))
      allocate(eta(0:nz-1))
      allocate(gam1_old(0:nz-1))
      allocate(grav_cell(0:nz-1))
      allocate(div_coeff_new(0:nz-1))

      ! Edge-centered
      allocate(edge(0:nz))
      allocate(beta(0:nz),beta_new(0:nz),beta_nh(0:nz))

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE RHOX0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do n = spec_comp,spec_comp+nspec-1

         ! compute the force -- include the geometric source term that
         ! results from expanding out the spherical divergence
         do j = 0,nz-1
            force(j) = -s0_old(j,n) * (vel(j+1) - vel(j)) / dr - &
                       2.0_dp_t*s0_old(j,n)*HALF*(vel(j) + vel(j+1))/z(j)
         end do

         call mkflux_1d(s0_old(:,n),edge,vel,force,1,dr,dt)

         do j = 0,nz-1
            s0_new(j,n) = s0_old(j,n) - &
                 dtdr / z(j)**2 * ( zl(j+1)**2 * edge(j+1) * vel(j+1) &
                                   -zl(j  )**2 * edge(j  ) * vel(j  ))
         end do

      enddo

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE RHO0 FROM RHOX0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j = 0,nz-1
        s0_new(j,rho_comp) =  s0_old(j,rho_comp)
        do n = spec_comp,spec_comp+nspec-1
          s0_new(j,rho_comp) =  s0_new(j,rho_comp) + (s0_new(j,n)-s0_old(j,n))
        end do
      end do

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE P0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! Put beta_old on edges
      call cell_to_edge(div_coeff_old,beta)
 
      ! Update p0 -- predictor
      do j = 0,nz-1
         divbetaw = one / (z(j)**2) * (zl(j+1)**2 * beta(j+1) * vel(j+1) - &
                                       zl(j  )**2 * beta(j  ) * vel(j  ) ) / dr
         betahalf = div_coeff_old(j)
         factor = half * dt * gam1(j) * (Sbar_in(j) - divbetaw / betahalf)
         p0_new(j) = p0_old(j) * (one + factor ) / (one - factor)
 
      end do
 
      do j = 0,nz-1
         ! (rho, p) --> T,h, etc

         den_row(1)  = s0_new(j,rho_comp)
        temp_row(1)  = s0_old(j,temp_comp) 
           p_row(1)  = p0_new(j)
         xn_row(1,:) = s0_new(j,spec_comp:spec_comp+nspec-1)/s0_new(j,rho_comp)

         gam1_old(j) = gam1(j)
 
         call eos(eos_input_rp, den_row, temp_row, & 
                  npts, nspec, & 
                  xn_row, & 
                  p_row, h_row, e_row, & 
                  cv_row, cp_row, xne_row, eta_row, pele_row, &
                  dpdt_row, dpdr_row, dedt_row, dedr_row, &
                  dpdX_row, dhdX_row, &
                  gam1_row, cs_row, s_row, & 
                  dsdt_row, dsdr_row, &
                  do_diag) 

        gam1(j) = gam1_row(1)
        s0_new(j,temp_comp) = temp_row(1)
      end do
 
      call make_grav_cell(grav_cell,s0_new(:,rho_comp))
 
      ! Define beta^n+1 at cell edges using the new gravity above
      call make_div_coeff(div_coeff_new,s0_new(:,rho_comp),p0_new,gam1,grav_cell,anelastic_cutoff)
      call cell_to_edge(div_coeff_new,beta_new)

      ! time-centered beta
      beta_nh = HALF*(beta + beta_new)
 
      ! Update p0 -- corrector
      do j = 0,nz-1
         divbetaw = one / (z(j)**2) * (zl(j+1)**2 * beta_nh(j+1) * vel(j+1) - &
                                       zl(j  )**2 * beta_nh(j  ) * vel(j  ) ) / dr
         betahalf = HALF*(div_coeff_old(j) + div_coeff_new(j))
         factor = half * dt * (Sbar_in(j) - divbetaw / betahalf)
         p0_new(j) = p0_old(j) * (one + factor * gam1_old(j)) / (one - factor * gam1(j))
 
      end do

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE RHOH0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j = 0,nz-1

         div_w0 = (vel(j+1) - vel(j)) / dr 

         force(j) = -s0_old(j,rhoh_comp) * div_w0 - &
              2.0_dp_t*s0_old(j,rhoh_comp)*HALF*(vel(j) + vel(j+1))/z(j)

         ! add eta at time-level n to the force for the prediction
         eta(j) = gam1_old(j) * p0_old(j) * (Sbar_in(j) - div_w0)
         force(j) = force(j) + eta(j)

         ! construct a new, time-centered eta for the final update
         eta(j) = HALF*(gam1(j)*p0_new(j) + gam1_old(j)*p0_old(j))* &
              (Sbar_in(j) - div_w0)
      end do

      call mkflux_1d(s0_old(:,rhoh_comp),edge,vel,force,1,dr,dt)

      do j = 0,nz-1

         s0_new(j,rhoh_comp) = s0_old(j,rhoh_comp) - &
              dtdr / z(j)**2 * ( zl(j+1)**2 * edge(j+1) * vel(j+1) &
                                -zl(j  )**2 * edge(j  ) * vel(j  ))

         s0_new(j,rhoh_comp) = s0_new(j,rhoh_comp) + dt * eta(j)

      end do


!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     MAKE TEMP0 AND GAM1 FROM P0 AND RHO0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j = 0,nz-1

         den_row(1)  = s0_new(j,rho_comp)
         temp_row(1) = s0_new(j,temp_comp)
         p_row(1)    = p0_new(j)
         xn_row(1,:) = s0_new(j,spec_comp:spec_comp+nspec-1)/s0_new(j,rho_comp)

         ! (rho,P) --> T, h
         call eos(eos_input_rp, den_row, temp_row, &
                  npts, nspec, &
                  xn_row, &
                  p_row, h_row, e_row, &
                  cv_row, cp_row, xne_row, eta_row, pele_row, &
                  dpdt_row, dpdr_row, dedt_row, dedr_row, &
                  dpdX_row, dhdX_row, &
                  gam1_row, cs_row, s_row, &
                  dsdt_row, dsdr_row, &
                  do_diag)

         s0_new(j,temp_comp) = temp_row(1)
         gam1(j) = gam1_row(1)

      end do

      deallocate(force,eta,edge,beta,beta_new,beta_nh,div_coeff_new,gam1_old,grav_cell)

   end subroutine advect_base_state_spherical

end module advect_base_module
