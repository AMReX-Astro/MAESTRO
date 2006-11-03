module advect_base_module

  use bl_types
  use bl_constants_module
  use bc_module
  use multifab_module
  use heating_module
  use mkflux_module
  use make_div_coeff_module
  use variables
  use eos_module
  use make_grav_module

  implicit none

contains

   subroutine advect_base(vel,Sbar_in,p0_old,p0_new, &
                          s0_old,s0_new,temp0, &
                          gam1,div_coeff_edge, &
                          dr,dt,anelastic_cutoff,spherical)

      real(kind=dp_t), intent(in   ) :: vel(:)
      real(kind=dp_t), intent(in   ) :: Sbar_in(:)
      real(kind=dp_t), intent(in   ) :: p0_old(:), s0_old(:,:)
      real(kind=dp_t), intent(  out) :: p0_new(:), s0_new(:,:)
      real(kind=dp_t), intent(inout) :: temp0(:),gam1(:),div_coeff_edge(:)
      real(kind=dp_t), intent(in   ) :: dr,dt,anelastic_cutoff
      integer        , intent(in   ) :: spherical

      integer :: i

      print *, '<<< advect base >>>'

      if (spherical == 0) then

        call advect_base_state_planar(vel,p0_old,p0_new,s0_old,s0_new,temp0, &
                                      gam1,dr,dt,anelastic_cutoff)

      else

        call advect_base_state_spherical(vel,Sbar_in,p0_old,p0_new,s0_old,s0_new,temp0, &
                                         gam1,div_coeff_edge,&
                                         dr,dt,anelastic_cutoff)
      end if


   end subroutine advect_base

   subroutine advect_base_state_planar (vel,p0_old,p0_new,s0_old,s0_new,temp0, &
                                        gam1,dr,dt,anelastic_cutoff)

      real(kind=dp_t), intent(in   ) :: vel(:)
      real(kind=dp_t), intent(in   ) :: p0_old(:), s0_old(:,:)
      real(kind=dp_t), intent(  out) :: p0_new(:), s0_new(:,:)
      real(kind=dp_t), intent(inout) :: temp0(:)
      real(kind=dp_t), intent(inout) :: gam1(:)
      real(kind=dp_t), intent(in   ) :: dr,dt,anelastic_cutoff

!     Local variables
      integer :: i, j, n, nz
      real(kind=dp_t) :: coeff

      real (kind = dp_t), allocatable :: H(:,:)
      real (kind = dp_t), allocatable :: force(:)
      real (kind = dp_t), allocatable :: edge(:)
      real (kind = dp_t), allocatable :: grav_edge(:)

      nz = size(p0_new,dim=1)

      allocate(    force(nz))
      allocate(     edge(nz+1))
      allocate(grav_edge(nz+1))

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE P0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      force = ZERO
      call mkflux_1d(p0_old,edge,vel,force,1,dr,dt)
      do j = 1,nz
        p0_new(j) = p0_old(j) - dt / dr * HALF * (vel(j) + vel(j+1)) * (edge(j+1) - edge(j))
      end do

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE RHOX0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do n = spec_comp,spec_comp+nspec-1
         do j = 1,nz
            force(j) = s0_old(j,n) * (vel(j+1) - vel(j)) / dr
         end do
         call mkflux_1d(s0_old(:,n),edge,vel,force,1,dr,dt)
         do j = 1,nz
            s0_new(j,n) = s0_old(j,n) - dt / dr * (edge(j+1) * vel(j+1) - edge(j) * vel(j))
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
!     MAKE TEMP0, RHOH0 AND GAM1 FROM P0 AND RHO0
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
                  do_diag)

         temp0(j) = temp_row(1)
         gam1(j) = gam1_row(1)

         s0_new(j,rhoh_comp) = s0_new(j,rho_comp) * h_row(1)

      end do

      deallocate(force,edge,grav_edge)

   end subroutine advect_base_state_planar

   subroutine advect_base_state_spherical (vel,Sbar_in,p0_old,p0_new,s0_old,s0_new,temp0, &
                                           gam1,div_coeff_edge,& 
                                           dr,dt,anelastic_cutoff)

      real(kind=dp_t), intent(in   ) :: vel(:),Sbar_in(:)
      real(kind=dp_t), intent(in   ) :: p0_old(:), s0_old(:,:)
      real(kind=dp_t), intent(  out) :: p0_new(:), s0_new(:,:)
      real(kind=dp_t), intent(inout) :: temp0(:)
      real(kind=dp_t), intent(inout) :: gam1(:)
      real(kind=dp_t), intent(inout) :: div_coeff_edge(:)
      real(kind=dp_t), intent(in   ) :: dr,dt,anelastic_cutoff

!     Local variables
      integer :: i, j, k, n, nz
      real(kind=dp_t) :: dtdr,divbetaw,betahalf,factor,integral

      real (kind = dp_t), allocatable :: z(:),zl(:)
      real (kind = dp_t), allocatable :: m(:)
      real (kind = dp_t), allocatable :: force(:)
      real (kind = dp_t), allocatable :: edge(:)
      real (kind = dp_t), allocatable :: beta(:),beta_new(:),beta_nh(:)
      real (kind = dp_t), allocatable :: gam1_old(:)
      real (kind = dp_t), allocatable :: grav_cell(:)
      real (kind = dp_t), allocatable :: grav_edge(:)

      dtdr = dt / dr

      nz = size(p0_new,dim=1)

      ! Cell-centered
      allocate(force(nz),z(nz))
      allocate(m(nz))
      allocate(gam1_old(nz))
      allocate(grav_cell(nz))

      ! Edge-centered
      allocate(edge(nz+1),zl(nz+1))
      allocate(beta(nz+1),beta_nh(nz+1))
      allocate(grav_edge(nz+1))

      ! z(j)  is the location of the cell center of cell (j)
      ! zl(j) is the location of the lower edge of cell (j)
      z(1) = 0.5_dp_t*dr
      do j = 2,nz
         z(j) = z(j-1) + dr
      enddo
      zl(1) = zero
      do j = 2,nz+1
         zl(j) = zl(j-1) + dr
      enddo

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE RHOX0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do n = spec_comp,spec_comp+nspec-1
         do j = 1,nz
            force(j) = s0_old(j,n) * (vel(j+1) - vel(j)) / dr
         end do
         call mkflux_1d(s0_old(:,n),edge,vel,force,1,dr,dt)
         do j = 1,nz
            s0_new(j,n) = s0_old(j,n) - dtdr / zl(j)**2 * ( z(j+1)**2 * edge(j+1) * vel(j+1) &
                                                           -z(  j)**2 * edge(j  ) * vel(j  ))
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

     do j = 1,nz+1
       beta(j) = div_coeff_edge(j)
     end do

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
                 do_diag) 
        gam1(j) = gam1_row(1)
     end do

     call make_grav_edge(grav_edge,s0_new(:,rho_comp),dr,1)
     call make_grav_cell(grav_cell,s0_new(:,rho_comp),dr,1)

     ! Define beta^n+1 at cell edges using the new gravity above
     beta_new(1) = 1.5d0 * s0_new(1,rho_comp) - 0.5d0 * s0_new(2,rho_comp)
     do j = 2,nz+1
        integral  = s0_new(j-1,rho_comp) * grav_cell(j) * dr / (gam1(j-1) * p0_new(j-1))
        beta_new(j) = beta_new(j-1) * exp(-integral)
     end do 

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
!     MAKE TEMP0, RHOH0 AND GAM1 FROM P0 AND RHO0
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
                  do_diag)

         temp0(j) = temp_row(1)
         gam1(j) = gam1_row(1)

         s0_new(j,rhoh_comp) = s0_new(j,rho_comp) * h_row(1)

      end do

      deallocate(z,zl,force,edge,beta,beta_nh,gam1_old,grav_edge,grav_cell)

   end subroutine advect_base_state_spherical

end module advect_base_module
