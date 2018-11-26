module base_state_module

  use bl_types
  use network, only: nspec
  
  implicit none

  private

  public :: init_base_state

contains

  subroutine init_base_state(n,model_file,s0_init,p0_init,dx)

    use bl_prof_module
    use fundamental_constants_module
    use parallel
    use bl_error_module
    use bl_constants_module
    use eos_module, only: eos_input_rp, eos
    use extern_probin_module, only: eos_gamma
    use eos_type_module
    use network, only: spec_names, network_species_index
    use probin_module, only: base_cutoff_density, anelastic_cutoff, &
                             buoyancy_cutoff_factor, &
                             prob_lo, prob_hi, &
                             small_temp, small_dens, grav_const, &
                             T_0,rho_0

    use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp, trac_comp, ntrac
    use geometry, only: dr, spherical, polar, nr
    use inlet_bc_module, only: set_inlet_bcs
    
    integer           , intent(in   ) :: n
    character(len=256), intent(in   ) :: model_file
    real(kind=dp_t)   , intent(inout) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(inout) :: p0_init(0:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    ! local
    integer         :: i,j,r,comp
    real(kind=dp_t) :: rloc,rmid
    real(kind=dp_t) :: d_ambient,t_ambient,p_ambient,xn_ambient(nspec)

    real(kind=dp_t) :: min_dens, max_dens, min_temp, max_temp,t_guess
    real(kind=dp_t) :: hp, N0, gas_const

    type(bl_prof_timer), save :: bpt
    
    real(kind=dp_t), parameter :: TINY = 1.0d-10

    real(kind=dp_t) :: dpdr, rhog
    real(kind=dp_t) :: max_hse_error

    real(kind=dp_t), parameter :: SMALL = 1.d-12

    type (eos_t) :: eos_state

    call build(bpt, "init_base_state")

887 format(78('-'))
888 format(a60,g18.10)
889 format(a60)

    if (spherical .eq. 1 .or. polar .eq. 1)  then
       call bl_error("ERROR: rt base_state is not valid for spherical or polar")
    endif


    if ( parallel_IOProcessor() .and. n == 1) then
       ! output block for cutoff density information
       write (*,887)
       write (*,*)   'cutoff densities:'
       write (*,888) '    low density cutoff (for mapping the model) =      ', &
            base_cutoff_density
       write (*,888) '    buoyancy cutoff density                           '
       write (*,888) '        (for zeroing rho - rho_0, centrifugal term) = ', &
            buoyancy_cutoff_factor*base_cutoff_density
       write (*,888) '    anelastic cutoff =                                ', &
            anelastic_cutoff
       write (*,888) ' '
    end if

    min_dens = rho_0

    if (anelastic_cutoff > min_dens .or. base_cutoff_density > min_dens) then
       call bl_error("ERROR: for the RT problem, the anelastic and base cutoff densities > min(rho)")
    endif

    if (min_dens < small_dens) then
       if ( parallel_IOProcessor() .and. n == 1) then
          print *, ' '
          print *, 'WARNING: minimum model density is lower than the EOS cutoff'
          print *, '         density, small_dens'
       endif
    endif

    if ( parallel_IOProcessor() .and. n == 1) then
       ! close the cutoff density output block
       write (*,887)
       write (*,*)   ' '
    end if
    
    !compute the gas constant R
    gas_const = k_B * n_A
    ! the pressure scale height is, assuming a molecular weight of 1
    hp = -gas_const*T_0/(grav_const)
    ! N_0 is then
    N0 = sqrt(-grav_const/hp * (eos_gamma - ONE)/eos_gamma)
    
    ! set a guess for the temperature for the EOS calls
    t_guess = 1.e-8

    ! fill the base state arrays
    do r=0,nr(n)-1

       ! height above the bottom of the domain
       rloc = prob_lo(size(dx)) + (dble(r) + HALF)*dr(n)
       d_ambient = rho_0 * exp(-rloc/hp)
       p_ambient = rho_0*gas_const*T_0 * exp(-rloc/hp)
       xn_ambient(:) = ONE
       t_ambient = t_guess
       
       ! use the EOS to make the state consistent
       eos_state%T     = t_ambient
       eos_state%rho   = d_ambient
       eos_state%p     = p_ambient
       eos_state%xn(:) = xn_ambient(:)

       ! (rho,p) --> T, h
       call eos(eos_input_rp, eos_state)

       s0_init(r, rho_comp) = d_ambient
       s0_init(r,rhoh_comp) = d_ambient * eos_state%h
       s0_init(r,spec_comp:spec_comp+nspec-1) = d_ambient * xn_ambient(1:nspec)
       p0_init(r) = eos_state%p
       s0_init(r,temp_comp) = eos_state%T

       if (ntrac .gt. 0) then
          s0_init(r,trac_comp:trac_comp+ntrac-1) = ZERO
       end if

    end do

    min_temp = minval(s0_init(:,temp_comp))

    if (min_temp < small_temp) then
       if ( parallel_IOProcessor() .and. n == 1) then
          print *, ' '
          print *, 'WARNING: minimum model temperature is lower than the EOS cutoff'
          print *, '         temperature, small_temp'
       endif
    endif



    max_hse_error = -1.d30

    do r=1,nr(n)-1

       rloc = prob_lo(size(dx)) + (dble(r) + HALF)*dr(n)

       dpdr = (p0_init(r) - p0_init(r-1))/dr(n)
       rhog = HALF*(s0_init(r,rho_comp) + s0_init(r-1,rho_comp))*grav_const

       max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(dpdr))

    enddo

    if ( parallel_IOProcessor() ) then
       write (*,*) ' '
       write (*,*) 'Maximum HSE Error = ', max_hse_error
       write (*,*) '   (after putting initial model into base state arrays, and'
       write (*,*) '    for density < base_cutoff_density)'
       write (*,887)
       write (*,*) ' '
    endif




    ! initialize any inlet BC parameters
    call set_inlet_bcs()


    call destroy(bpt)

  end subroutine init_base_state

end module base_state_module
