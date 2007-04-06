module correct_base_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use mkflux_module
  use eos_module
  use variables
  use geometry
  use make_grav_module

  implicit none

contains

   subroutine correct_base(rhopertbar,p0,s0,temp0,gam1, &
                          grav_cell,div_coeff_edge, &
                          dz,dt,anelastic_cutoff)

      real(kind=dp_t), intent(in   ) :: rhopertbar(:)
      real(kind=dp_t), intent(inout) :: p0(:), s0(:,:)
      real(kind=dp_t), intent(inout) :: temp0(:),gam1(:),div_coeff_edge(:)
      real(kind=dp_t), intent(in   ) :: grav_cell(:)
      real(kind=dp_t), intent(in   ) :: dz,dt,anelastic_cutoff

      integer :: i

      if (spherical .eq. 0) then

        call correct_base_state_planar(rhopertbar,p0,s0,temp0,gam1,grav_cell,dz, &
                                       dt,anelastic_cutoff)

      else

        print *,'NO ROUTINE YET FOR SPHERICAL CORRECT_BASE'
        stop
!       call correct_base_state_spherical(rhopertbar,p0,s0,temp0,gam1,div_coeff_edge,&
!                                        dt,anelastic_cutoff)
      end if


   end subroutine correct_base

   subroutine correct_base_state_planar (rhopertbar,p0,s0,temp0,gam1,grav_cell,dz,dt,anelastic_cutoff)

      real(kind=dp_t), intent(in   ) :: rhopertbar(0:)
      real(kind=dp_t), intent(inout) :: p0(0:), s0(0:,:)
      real(kind=dp_t), intent(inout) :: temp0(0:)
      real(kind=dp_t), intent(inout) :: gam1(0:)
      real(kind=dp_t), intent(in   ) :: grav_cell(0:)
      real(kind=dp_t), intent(in   ) :: dz,dt,anelastic_cutoff


!     Local variables
      integer :: i, j, n, nz
      real(kind=dp_t) :: xfrac(nspec), enthalpy
      real(kind=dp_t), allocatable :: deltap0(:)

      nz = size(p0,dim=1)
      allocate(deltap0(0:nz-1))

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     CORRECT RHOX0 AND RHOH0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j = 0,nz-1
        ! Compute X_i and h
        do n = spec_comp,spec_comp+nspec-1
          xfrac(n-spec_comp+1) = s0(j,n) / s0(j,rho_comp)
        end do
        enthalpy = s0(j,rhoh_comp) / s0(j,rho_comp)

        ! Update rho
        s0(j,rho_comp) = s0(j,rho_comp) + rhopertbar(j)

        ! Compute new (rho X)_i and (rho h)
        do n = spec_comp,spec_comp+nspec-1
          s0(j,n) = s0(j,rho_comp) * xfrac(n-spec_comp+1)
        end do
        s0(j,rhoh_comp) = enthalpy * s0(j,rho_comp)
      end do

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     UPDATE P0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j = 0,nz-1
        deltap0(j) = -rhopertbar(j) * grav_cell(j) * dz
      end do

      do j = 1,nz-1
        p0(j) = p0(j) + HALF * (deltap0(j-1) + deltap0(j))
      end do

!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     MAKE TEMP0 AND GAM1 FROM P0 AND RHO0
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j = 0,nz-1

         den_row(1)  = s0(j,rho_comp)
         temp_row(1) = temp0(j)
         p_row(1)    = p0(j)
         xn_zone(1:) = s0(j,spec_comp:)/s0(j,rho_comp)

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

      deallocate(deltap0)

   end subroutine correct_base_state_planar

end module correct_base_module
