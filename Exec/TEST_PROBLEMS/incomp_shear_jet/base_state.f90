!Incompressible shear jet base state
!Commenting is heavy to facilitate understanding for a MAESTRO new-comer (me).

!MAESTRO models fluids with a 1D base state that is always in hydrostatic
!equilibrium (HSE) overlaid by a state that captures deviations from HSE.

!Two geometries are currently supported:
!   1) Plane-parallel: In this case the base state grid is parallel to ("locked
!   to") the vertical dimension of the non-HSE state.
!   2) Spherical: In this case the base state runs radially outwards from the
!   center while the non-HSE state is represented with a Cartesian grid.  In
!   this case, the base-state grid is not locked to the non-HSE grid.

module base_state_module
  !TODO: Delete any unused modules/variables
  use bl_types              !BoxLib types
  use network, only: nspec  !Number of species
  
  implicit none

  private

  public :: init_base_state

contains

  subroutine init_base_state(n,model_file,s0_init,p0_init,dx)
    !--- Variable declaration / initialization ---
    !-- Modules --
    use bl_prof_module          !BoxLib profiling module. Contains utilities for
                                !creating timers to facilitate optimization,
                                !i.e. profiling. Alternatively, one can use
                                !GNU's profiler 'gprof.' 
    use parallel                !BoxLib parallel module.
    use bl_error_module         
    use bl_constants_module
    use eos_module, only: eos, eos_input_rp    !Equation of state module
    use eos_type_module
    use network,       only: spec_names, network_species_index
    !Probin module gives global access to all problem parameters.
    use probin_module, only: base_cutoff_density, anelastic_cutoff, &   
                             buoyancy_cutoff_factor, &
                             prob_lo, prob_hi, &
                             small_temp, small_dens, &
                             rho_base, p_base !, do_SNe
    !Variables is a module containing indices used for accessing variables in 
    !storage arrays (e.g. the 's' multifab).
        !rho_comp  - density component of fluid state
        !rhoh_comp - density*enthalpy
        !temp_comp - temperature
        !spec_comp - species mass fractions (stored as X*rho)
        !trac_comp - tracers
        !ntrac     - number of tracers
    use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp, trac_comp, ntrac
    use geometry,  only: dr, spherical, nr       !Module w/ info about chosen
                                                 !geometry
                                                    !nr - number of base state
                                                    !cells
                                                    !dr - base state grid size
    !-- Args --
    integer           , intent(in   ) :: n              !Number of base state cells
    character(len=256), intent(in   ) :: model_file     !Not used for this problem
    real(kind=dp_t)   , intent(inout) :: s0_init(0:,:)  !Initial base state fluid state
    real(kind=dp_t)   , intent(inout) :: p0_init(0:)    !Initial base state pressure
    real(kind=dp_t)   , intent(in   ) :: dx(:)          !Grid size array

    !-- Local variables --
    integer         :: i,j,r,comp
    integer         :: ia, ib
    real(kind=dp_t) :: rloc,rshr1,rshr2
    real(kind=dp_t) :: d_ambient,t_ambient,p_ambient,xn_ambient(nspec)
    real(kind=dp_t) :: t_guess
    real(kind=dp_t) :: xn_jet(nspec), xn_still(nspec)
    real(kind=dp_t) :: min_dens, max_dens, min_temp, max_temp
    real(kind=dp_t) :: dpdr, rhog
    real(kind=dp_t) :: max_hse_error
    type(bl_prof_timer), save      :: bpt
    real(kind=dp_t),     parameter :: SMALL   = 1.d-12
    real(kind=dp_t),     parameter :: TINY    = 1.0d-10
    character(len=*),    parameter :: FMT_SEP = "(78('-'))"
    character(len=*),    parameter :: FMT_MSG = "(a60,g18.10)"
    type (eos_t) :: eos_state


    !--- Execution ---
    call build(bpt, "init_base_state")  !Build the timer

    !-- Error checking --
    !If the spherical geometry is chosen by the
    !user we raise an error since isj is not compatible
    if (spherical .eq. 1) then
       call bl_error("ERROR: Incompressible shear jet is not valid for spherical geometry")
    endif

    !Here's an example of using the parallel module.  The module assigns one
    !processor for I/O and parallel_IOProcessor() returns true if the currently
    !executing processor is the I/O processor.
    if ( parallel_IOProcessor() .and. n == 1) then
       ! Output block for cutoff densities
       write (*,FMT_SEP)
       write (*,*)   'cutoff densities:'
       write (*,FMT_MSG) '    low density cutoff (for mapping the model) =      ', &
            base_cutoff_density
       write (*,FMT_MSG) '    buoyancy cutoff density                           '
       write (*,FMT_MSG) '        (for zeroing rho - rho_0, centrifugal term) = ', &
            buoyancy_cutoff_factor*base_cutoff_density
       write (*,FMT_MSG) '    anelastic cutoff =                                ', &
            anelastic_cutoff
       write (*,FMT_MSG) ' '
    end if

    !Check min density
    min_dens = rho_base
    if (min_dens < small_dens) then
       if ( parallel_IOProcessor() .and. n == 1) then
          print *, ' '
          print *, 'WARNING: minimum model density is lower than the EOS cutoff'
          print *, '         density, small_dens'
       endif
    endif

    if ( parallel_IOProcessor() .and. n == 1) then
       ! Close the cutoff density output block
       write (*,FMT_SEP)
       write (*,*)   ' '
    end if

    !-- Initialize s0_init and p0_init --
    !Location of the two thin shear layers
    rshr1 = 0.25d0*(prob_lo(size(dx)) + prob_hi(size(dx)))
    rshr2 = 0.75d0*(prob_lo(size(dx)) + prob_hi(size(dx)))

    ! A and B act as tags.  
    ! A tags the part of the fluid initially in the "jet zone," i.e. the part of
    ! the fluid with positive initial velocity.
    ! B tags the part of the fluid with negative initial velocity
    ia = network_species_index("A")
    ib = network_species_index("B")

    ! xn is an array where the ith element is the mass fraction (percent) of the
    ! ith species.  Initially, the middle half of the domain is composed entirely of
    ! "A" particles and negligible "B" particles, while the reverse is true of
    ! the outer regions of the domain.

    ! As we evolve in time these arrays track how the two fluids mix.
    xn_jet(:)      = SMALL
    xn_jet(ia)     = ONE - (nspec-1)*SMALL

    xn_still(:)  = SMALL
    xn_still(ib) = ONE - (nspec-1)*SMALL

    ! set a guess for the temperature of the EOS calls
    t_guess = 1.e-8

    !Iterate through each base state cell and initialize
    !   -The components of the fluid state 's':
    !       density, enthalpy, species mass fractions, and temperature
    !   -The pressure (note, pressure is NOT a component of the 's' multifab)
    do r=0,nr(n)-1
        ! Current height above the bottom of the domain
        rloc = (dble(r) + HALF)*dr(n)

        !Init density, pressure, and temp
        d_ambient = rho_base
        p_ambient = p_base
        t_ambient = t_guess
        
        ! Depending on location, initialize the mass fraction
        if (rloc > rshr1 .and. rloc < rshr2) then
           ! Middle -- jet
           xn_ambient(:) = xn_jet(:)
        else
           ! Outer regions -- still fluid
           xn_ambient(:) = xn_still(:)
        endif

        ! use the EOS to make the state consistent
        ! We set density and pressure, and from this the EoS yields many
        ! thermodynamic quantities (temperature and enthalpy being the two we
        ! care about in this problem).
        eos_state%T     = t_ambient
        eos_state%rho   = d_ambient
        eos_state%p     = p_ambient
        eos_state%xn(:) = xn_ambient(:)

        ! (rho,p) --> T, h
        call eos(eos_input_rp, eos_state, .false.)

        !Now that we've calculated all of the ambient values and churned them
        !through the EoS we can finally initialize the fluid state.
        s0_init(r, rho_comp)                    = d_ambient
        s0_init(r, rhoh_comp)                   = d_ambient * eos_state%h
        s0_init(r, spec_comp:spec_comp+nspec-1) = d_ambient * xn_ambient(1:nspec)
        s0_init(r, temp_comp)                   = eos_state%T
        p0_init(r)                              = p_base

        !We don't use tracers in this problem, so ntrac is always 0.
        if (ntrac .gt. 0) then
           s0_init(r,trac_comp:trac_comp+ntrac-1) = ZERO
        end if
    end do

    !-- Post State Initialization --
    !Check that the temperature is consistent with the EoS 
    min_temp = minval(s0_init(:,temp_comp))
    if (min_temp < small_temp) then
       if ( parallel_IOProcessor() .and. n == 1) then
          print *, ' '
          print *, 'WARNING: minimum model temperature is lower than the EOS cutoff'
          print *, '         temperature, small_temp'
       endif
    endif

    call destroy(bpt)
  end subroutine init_base_state
end module base_state_module
