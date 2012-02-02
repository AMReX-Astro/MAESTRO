module base_state_module

  use bl_types
  use network, only: nspec
  
  implicit none

  private

  public :: init_base_state

contains

  subroutine init_base_state(n,model_file,s0_init,p0_init,dx,&
       diffusion_coefficient)

    use bl_prof_module
    use parallel
    use bl_error_module
    use bl_constants_module
    use eos_module
    use network, only: network_species_index
    use probin_module, only: ambient_h, ambient_dens, &
                             ambient_he4, ambient_c12, ambient_fe56, &
                             thermal_conductivity
    use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp,  &
                         trac_comp, ntrac
    use geometry, only: nr
    use inlet_bc_module, only: set_inlet_bcs

    
    integer           , intent(in   ) :: n
    character(len=256), intent(in   ) :: model_file
    real(kind=dp_t)   , intent(inout) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(inout) :: p0_init(0:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)
    real(kind=dp_t)   , intent(  out) :: diffusion_coefficient

    ! local
    integer :: ihe4, ic12, ife56
    integer :: r

    type(bl_prof_timer), save :: bpt
    

    call build(bpt, "init_base_state")

    ihe4  = network_species_index("He4")
    ic12  = network_species_index("C12")
    ife56 = network_species_index("Fe56")

    if (ihe4 < 0 .or. ic12 < 0 .or. ife56 < 0) &
       call bl_error("Invalid species in init_base_state.")

    h_eos = ambient_h
    den_eos  = ambient_dens
    
    xn_eos(ihe4)  = ambient_he4
    xn_eos(ic12)  = ambient_c12
    xn_eos(ife56) = ambient_fe56

    call eos(eos_input_rh, den_eos, temp_eos, &
         xn_eos, &
         p_eos, h_eos, e_eos, &
         cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
         dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
         dpdX_eos, dhdX_eos, &
         gam1_eos, cs_eos, s_eos, &
         dsdt_eos, dsdr_eos, &
         .false.)       

    diffusion_coefficient = thermal_conductivity / (cp_eos * ambient_dens)

    do r = 0, nr(n) - 1

       s0_init(r,rho_comp ) = den_eos
       s0_init(r,rhoh_comp) = den_eos * h_eos
       s0_init(r,temp_comp) = temp_eos

       s0_init(r,spec_comp:spec_comp+nspec-1) = den_eos * xn_eos(1:nspec)

       p0_init(r) = p_eos

       if (ntrac .gt. 0) then
          s0_init(r,trac_comp:trac_comp+ntrac-1) = ZERO
       end if

    enddo

    ! initialize any inlet BC parameters
    call set_inlet_bcs()

    call destroy(bpt)

  end subroutine init_base_state

end module base_state_module
