module react_base_module
  ! 
  ! This is React Base from paper III.
  ! 
  use bl_types

  implicit none
  
  private

  public :: react_base
  
contains
  
  subroutine react_base(nlevs,s0_in,rho_omegadotbar,rho_Hextbar,dt_in, &
                        s0_out,gamma1bar_out)

    use geometry, only: nr
    use variables, only: rho_comp, spec_comp, rhoh_comp
    use eos_module
    use bl_prof_module
     
    integer        , intent(in   ) :: nlevs
    real(kind=dp_t), intent(in   ) :: s0_in(:,0:,:)
    real(kind=dp_t), intent(in   ) :: rho_omegadotbar(:,0:,:)
    real(kind=dp_t), intent(in   ) :: rho_Hextbar(:,0:)
    real(kind=dp_t), intent(in   ) :: dt_in
    real(kind=dp_t), intent(  out) :: s0_out(:,0:,:)
    real(kind=dp_t), intent(inout) :: gamma1bar_out(:,0:)
    
    integer :: n,r,comp

    type(bl_prof_timer), save :: bpt

    call build(bpt, "react_base")
    
    do n=1,nlevs

       do r = 0,nr(n)-1
          
          ! rho_out = rho_in
          s0_out(n,r,rho_comp) = s0_in(n,r,rho_comp)
          
          ! update enthalpy
          s0_out(n,r,rhoh_comp) = s0_in(n,r,rhoh_comp)
          do comp = spec_comp,spec_comp+nspec-1
             s0_out(n,r,rhoh_comp) = s0_out(n,r,rhoh_comp) &
                  -dt_in*rho_omegadotbar(n,r,comp-spec_comp+1)*ebin(comp-spec_comp+1)
          end do
          s0_out(n,r,rhoh_comp) = s0_out(n,r,rhoh_comp) + dt_in * rho_Hextbar(n,r)

       end do
       
    enddo ! end loop over levels

    call destroy(bpt)

  end subroutine react_base

end module react_base_module
