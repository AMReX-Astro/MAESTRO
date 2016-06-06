module base_state_module

  use bl_types
  use network, only: nspec, network_init
  
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
    use eos_module, only: eos_input_rh, eos
    use eos_type_module
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
    type (eos_t) :: eos_state

    type(bl_prof_timer), save :: bpt
    
    call network_init()

    call build(bpt, "init_base_state")

    ihe4  = network_species_index("helium-4")
    ic12  = network_species_index("carbon-12")
    ife56 = network_species_index("iron-56")

    if (ihe4 < 0 .or. ic12 < 0 .or. ife56 < 0) then
       print *, ihe4, ic12, ife56
       call bl_error("Invalid species in init_base_state.")
    endif

    eos_state%h         = ambient_h
    eos_state%rho       = ambient_dens

    eos_state%xn(:) = ZERO
    eos_state%xn(ihe4)  = ambient_he4
    eos_state%xn(ic12)  = ambient_c12
    eos_state%xn(ife56) = ambient_fe56

    call eos(eos_input_rh, eos_state)

    diffusion_coefficient = thermal_conductivity / (eos_state%cp * ambient_dens)

    do r = 0, nr(n) - 1

       s0_init(r,rho_comp ) = eos_state%rho
       s0_init(r,rhoh_comp) = eos_state%rho * eos_state%h
       s0_init(r,temp_comp) = eos_state%T

       s0_init(r,spec_comp:spec_comp+nspec-1) = eos_state%rho * eos_state%xn(1:nspec)

       p0_init(r) = eos_state%p

       if (ntrac .gt. 0) then
          s0_init(r,trac_comp:trac_comp+ntrac-1) = ZERO
       end if

    enddo

    ! initialize any inlet BC parameters
    call set_inlet_bcs()

    call destroy(bpt)

  end subroutine init_base_state

end module base_state_module
