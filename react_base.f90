module react_base_module
  ! 
  ! This is React Base from paper III.
  ! 
  use bl_types

  implicit none
  
  private

  public :: react_base
  
contains
  
  subroutine react_base(nlevs,p0_in,s0_in,rho_omegadotbar,rho_Hextbar,dt_in,p0_out, &
                        s0_out,gamma1bar_out)

    use geometry, only: nr
    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp
    use eos_module
    use bl_prof_module
    use probin_module, ONLY: use_big_h, enthalpy_pred_type
     
    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(in   ) :: p0_in(:,0:), s0_in(:,0:,:)
    real(kind=dp_t), intent(in   ) :: rho_omegadotbar(:,0:,:)
    real(kind=dp_t), intent(in   ) :: rho_Hextbar(:,0:)
    real(kind=dp_t), intent(in   ) :: dt_in
    real(kind=dp_t), intent(  out) :: p0_out(:,0:), s0_out(:,0:,:)
    real(kind=dp_t), intent(inout) :: gamma1bar_out(:,0:)
    
    integer :: n,r,comp

    type(bl_prof_timer), save :: bpt

    call build(bpt, "react_base")
    
    do n=1,nlevs

       do r = 0,nr(n)-1
          
          ! (rho X)_out = (rho X)_in + dt_in * (rho omegadotbar)_in
          do comp = spec_comp,spec_comp+nspec-1
             s0_out(n,r,comp) = s0_in(n,r,comp) &
                  + dt_in * rho_omegadotbar(n,r,comp-spec_comp+1) 
          end do
          
          ! p_out = p_in
          p0_out(n,r) = p0_in(n,r)
          
          ! rho_out = rho_in
          s0_out(n,r,rho_comp) = s0_in(n,r,rho_comp)
          
          den_eos(1)  = s0_in(n,r,rho_comp)
          temp_eos(1) = s0_in(n,r,temp_comp)
          p_eos(1)    = p0_in(n,r)
          
          do comp = spec_comp,spec_comp+nspec-1
             xn_eos(1,comp-spec_comp+1) = s0_out(n,r,comp)/s0_out(n,r,rho_comp)
          end do
          
          s0_out(n,r,rhoh_comp) = s0_in(n,r,rhoh_comp)
          if(.not. use_big_h) then
             do comp = spec_comp,spec_comp+nspec-1
                s0_out(n,r,rhoh_comp) = s0_out(n,r,rhoh_comp) &
                     -dt_in*rho_omegadotbar(n,r,comp-spec_comp+1)*ebin(comp-spec_comp+1)
             end do
          endif
          s0_out(n,r,rhoh_comp) = s0_out(n,r,rhoh_comp) + dt_in * rho_Hextbar(n,r)
          
          ! We shouldn't update temp here since we don't update it in react-state.
          s0_out(n,r,temp_comp) = s0_in(n,r,temp_comp)
          
          if (enthalpy_pred_type .eq. 1) then
             ! Only do this to evaluate a new Gamma.
             ! (rho,P,X) --> T, h
             call eos(eos_input_rp, den_eos, temp_eos, &
                      npts, nspec, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      do_diag)

             gamma1bar_out(n,r) = gam1_eos(1)

          end if
          
       end do
       
    enddo ! end loop over levels

    call destroy(bpt)

  end subroutine react_base

end module react_base_module
