module react_base_module

  ! 
  ! This is React Base from paper III.
  ! 

  use bl_types
  use bl_constants_module
  use multifab_module
  use variables
  use eos_module
  use network

  implicit none

contains

   subroutine react_base(p0_in,s0_in,temp0_in,rho_omegadotbar,rho_Hextbar,dr,dt_in,p0_out,s0_out,gam1_out)

      real(kind=dp_t), intent(in   ) :: p0_in( :), s0_in( :,:), temp0_in(:)
      real(kind=dp_t), intent(in   ) :: rho_omegadotbar(:,:)
      real(kind=dp_t), intent(in   ) :: rho_Hextbar(:)
      real(kind=dp_t), intent(in   ) :: dr, dt_in
      real(kind=dp_t), intent(  out) :: p0_out(:), s0_out(:,:)
      real(kind=dp_t), intent(inout) :: gam1_out(:)

      integer :: j,n,nz

      nz    = size(rho_omegadotbar,dim=1)

      do j = 1,nz

         ! (rho X)_out = (rho X)_in + dt_in * (rho omegadotbar)_in
         do n = spec_comp,spec_comp+nspec-1
           s0_out(j,n) = s0_in(j,n) + dt_in * rho_omegadotbar(j,n-spec_comp+1) 
         end do

         ! p_out = p_in
         p0_out(j) = p0_in(j)

         ! rho_out = rho_in
         s0_out(j,rho_comp) = s0_in(j,rho_comp)

         den_row(1)  = s0_in(j,rho_comp)
         temp_row(1) = temp0_in(j)
         p_row(1)    = p0_in(j)

         do n = spec_comp,spec_comp+nspec-1
           xn_zone(n-spec_comp+1) = s0_out(j,n)/s0_out(j,rho_comp)
         end do

         ! (rho,P,X) --> T, h
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

         s0_out(j,rhoh_comp) = s0_in(j,rhoh_comp)
         do n = spec_comp,spec_comp+nspec-1
           s0_out(j,rhoh_comp) = s0_out(j,rhoh_comp) &
             - dt_in * rho_omegadotbar(j,n-spec_comp+1) * ebin(n-spec_comp+1)
         end do
         s0_out(j,rhoh_comp) = s0_out(j,rhoh_comp) + dt_in * rho_Hextbar(j)

         gam1_out(j) = gam1_row(1)

      end do

   end subroutine react_base

end module react_base_module
