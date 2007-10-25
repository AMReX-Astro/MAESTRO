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
  use probin_module

  implicit none

contains

   subroutine react_base(p0_in,s0_in,rho_omegadotbar,rho_Hextbar,dt_in,p0_out,s0_out,gam1_out)

      real(kind=dp_t), intent(in   ) :: p0_in(0:), s0_in(0:,:)
      real(kind=dp_t), intent(in   ) :: rho_omegadotbar(0:,:)
      real(kind=dp_t), intent(in   ) :: rho_Hextbar(0:)
      real(kind=dp_t), intent(in   ) :: dt_in
      real(kind=dp_t), intent(  out) :: p0_out(0:), s0_out(0:,:)
      real(kind=dp_t), intent(inout) :: gam1_out(0:)

      integer :: j,n,nz

      nz    = size(rho_omegadotbar,dim=1)

      do j = 0,nz-1

         ! (rho X)_out = (rho X)_in + dt_in * (rho omegadotbar)_in
         do n = spec_comp,spec_comp+nspec-1
           s0_out(j,n) = s0_in(j,n) + dt_in * rho_omegadotbar(j,n-spec_comp+1) 
         end do

         ! p_out = p_in
         p0_out(j) = p0_in(j)

         ! rho_out = rho_in
         s0_out(j,rho_comp) = s0_in(j,rho_comp)

         den_row(1)  = s0_in(j,rho_comp)
         temp_row(1) = s0_in(j,temp_comp)
         p_row(1)    = p0_in(j)

         do n = spec_comp,spec_comp+nspec-1
           xn_row(1,n-spec_comp+1) = s0_out(j,n)/s0_out(j,rho_comp)
         end do

         s0_out(j,rhoh_comp) = s0_in(j,rhoh_comp)
         if(.not. use_big_h) then
            do n = spec_comp,spec_comp+nspec-1
               s0_out(j,rhoh_comp) = s0_out(j,rhoh_comp) &
                    -dt_in*rho_omegadotbar(j,n-spec_comp+1)*ebin(n-spec_comp+1)
            end do
         endif
         s0_out(j,rhoh_comp) = s0_out(j,rhoh_comp) + dt_in * rho_Hextbar(j)

         ! Only do this to evaluate a new Gamma.
         ! (rho,P,X) --> T, h
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

!        We shouldn't update temp here since we don't update it in react-state.
         s0_out(j,temp_comp) = s0_in(j,temp_comp)

         gam1_out(j) = gam1_row(1)

      end do

   end subroutine react_base

end module react_base_module
