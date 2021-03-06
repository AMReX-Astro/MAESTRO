! init_base_state is used to initialize the base state arrays from the
! model file.  The actual reading of the model file is handled by the
! model_parser_module in Util/
!
! Note: The initial base state quantities returned from this routine
! are only a temporary base state.  These quantities are mapped onto
! the full 2- or 3-d state in initscaldata.f90 and a new base state is
! created after initialization by averaging the density and calling
! enforce_HSE in initialize.f90.

module base_state_module

  use bl_types, only: dp_t
  use network, only: nspec
  
  implicit none

  real(dp_t), save :: base_cutoff_density_loc
  real(dp_t), save :: rho_above_cutoff, rhoh_above_cutoff
  real(dp_t), save :: spec_above_cutoff(nspec), p_above_cutoff
  real(dp_t), save :: temp_above_cutoff
  real(dp_t), save :: trac_above_cutoff(10) ! if ntrac=0 we'd have a problem if we used ntrac

  private

  public :: init_base_state

contains

  subroutine init_base_state(n,model_file,s0_init,p0_init,dx)

    use bl_prof_module, only: bl_prof_timer, build, destroy
    use parallel, only: parallel_IOProcessor
    use bl_constants_module, only: ZERO, HALF, ONE, FOUR3RD, M_PI
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use network, only: spec_names
    use probin_module, only: base_cutoff_density, prob_lo, &
                             grav_const, planar_invsq_mass, &
                             do_planar_invsq_grav, do_2d_planar_octant, &
                             print_init_hse_diag
    use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp, trac_comp, ntrac
    use geometry, only: dr, spherical, nr
    use inlet_bc_module, only: set_inlet_bcs
    use fundamental_constants_module, only: Gconst
    use model_parser_module, only: read_model_file, interpolate, npts_model, &
                                   model_r, &
                                   idens_model, itemp_model, ipres_model, &
                                   ispec_model
    use simple_log_module

    integer           , intent(in   ) :: n
    character(len=256), intent(in   ) :: model_file
    real(kind=dp_t)   , intent(inout) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(inout) :: p0_init(0:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    ! local
    integer         :: r,comp
    real(kind=dp_t) :: rloc,dr_in,rmax,starting_rad,mod_dr
    real(kind=dp_t) :: d_ambient,t_ambient,p_ambient,xn_ambient(nspec)
    real(kind=dp_t) :: sumX

    type(bl_prof_timer), save :: bpt
    
    real(kind=dp_t), parameter :: TINY = 1.0d-10

    real(kind=dp_t) :: mencl, g, r_l, r_r, dpdr, rhog
    real(kind=dp_t) :: max_hse_error

    logical, save :: firstCall = .true.

    type (eos_t) :: eos_state

    call build(bpt, "init_base_state")

887 format(78('-'))
888 format(a60,g18.10)

    if (firstCall) then

       base_cutoff_density_loc = 1.d99

       ! only need to read in the initial model once -- model_parser_module
       ! stores the model data
       call read_model_file(model_file)

       firstCall = .false.

    endif


    dr_in = (model_r(npts_model) - model_r(1)) / dble(npts_model-1)
    rmax = model_r(npts_model)

    if ( parallel_IOProcessor() ) then
       write (*,887)
       if (spherical .ne. 1) then
          call log('model file mapping, level: ', n)
       else
          call log('model file mapping (spherical base state)')
       endif

       call log('dr of MAESTRO base state =                            ', dr(n))
       call log('dr of input file data =                               ', dr_in)
       call log(' ')
       call log('maximum radius (cell-centered) of input model =       ', rmax)
       
       if (dr(n) .lt. dr_in) then
          mod_dr = mod(dr_in,dr(n))
       else
          mod_dr = mod(dr(n),dr_in)
       endif

       if (mod_dr .gt. TINY) then
          call log(' ')
          call log("WARNING: resolution of base state array is not an integer")
          call log("         multiple of the initial model's resolution.     ")
          call log("         make sure this is a desired property as this    ")
          call log("         could lead to aliasing when performing the      ")
          call log("         interpolation.                                  ")
          call log(" ")
          call log("modulus = ", mod_dr)
       endif
    end if

    if (spherical .eq. 0) then
       starting_rad = prob_lo(size(dx))
    else
       starting_rad = ZERO
    endif
    
    do r=0,nr(n)-1

       rloc = starting_rad + (dble(r) + HALF)*dr(n)

       ! here we account for r > rmax of the model.hse array, assuming
       ! that the state stays constant beyond rmax
       rloc = min(rloc, rmax)

       ! also, if we've falled below the cutoff density, just keep the
       ! model constant
       if (rloc .gt. base_cutoff_density_loc) then

          s0_init(r,rho_comp) = rho_above_cutoff
          s0_init(r,rhoh_comp) = rhoh_above_cutoff
          s0_init(r,spec_comp:spec_comp+nspec-1) = spec_above_cutoff(1:nspec)
          p0_init(r) = p_above_cutoff
          s0_init(r,temp_comp) = temp_above_cutoff
          if(ntrac .gt. 0) then
             s0_init(r,trac_comp:trac_comp+ntrac-1) = trac_above_cutoff(1:ntrac)
          end if

       else

          d_ambient = interpolate(rloc, idens_model)
          t_ambient = interpolate(rloc, itemp_model)
          p_ambient = interpolate(rloc, ipres_model)

          sumX = ZERO
          do comp = 1, nspec
             xn_ambient(comp) = max(ZERO, min(ONE, &
                  interpolate(rloc, ispec_model-1+comp)))
             sumX = sumX + xn_ambient(comp)
          enddo
          xn_ambient = xn_ambient/sumX

          ! use the EOS to make the state consistent
          eos_state%T     = t_ambient
          eos_state%rho   = d_ambient
          eos_state%p     = p_ambient
          eos_state%xn(:) = xn_ambient(:)

          ! (rho,T) --> p,h
          call eos(eos_input_rt, eos_state)

          s0_init(r, rho_comp ) = d_ambient
          s0_init(r,rhoh_comp ) = d_ambient * eos_state%h
          s0_init(r,spec_comp:spec_comp+nspec-1) = d_ambient * xn_ambient(1:nspec)
          p0_init(r) = eos_state%p ! p_ambient !

          s0_init(r,temp_comp) = t_ambient
          
          if (ntrac .gt. 0) then
             s0_init(r,trac_comp:trac_comp+ntrac-1) = ZERO
          end if

          ! keep track of the height where we drop below the cutoff density
          if (s0_init(r,rho_comp) .le. base_cutoff_density .and. &
               base_cutoff_density_loc .eq. 1.d99) then

             if ( parallel_IOProcessor() ) then
                write (*,*) ' '
                write (*,*) 'setting r_cutoff to ', r
             end if

             base_cutoff_density_loc = rloc

             rho_above_cutoff = s0_init(r,rho_comp)
             rhoh_above_cutoff = s0_init(r,rhoh_comp)
             spec_above_cutoff(1:nspec) = s0_init(r,spec_comp:spec_comp+nspec-1)
             p_above_cutoff = p0_init(r)
             temp_above_cutoff = s0_init(r,temp_comp)
             if (ntrac .gt. 0) then
                trac_above_cutoff(1:ntrac) = s0_init(r,trac_comp:trac_comp+ntrac-1)
             end if

          end if

       end if

    end do

    ! check whether we are in HSE

    mencl = zero
    
    if (spherical .eq. 1 .OR. do_2d_planar_octant .eq. 1) then
       mencl = four3rd*m_pi*dr(n)**3*s0_init(0,rho_comp)
    endif

    max_hse_error = -1.d30

    do r=1,nr(n)-1
       
       rloc = starting_rad + (dble(r) + HALF)*dr(n)
       rloc = min(rloc, rmax)

       if (rloc .lt. base_cutoff_density_loc) then

          r_r = starting_rad + dble(r+1)*dr(n)
          r_l = starting_rad + dble(r)*dr(n)

          if (spherical .eq. 1 .OR. do_2d_planar_octant .eq. 1) then
             g = -Gconst*mencl/r_l**2
             mencl = mencl &
                  + four3rd*m_pi*dr(n)*(r_l**2+r_l*r_r+r_r**2)*s0_init(r,rho_comp)
          else
             if (.not. do_planar_invsq_grav) then
                g = grav_const
             else
                g = -Gconst*planar_invsq_mass / r_l**2
             endif
          endif

          dpdr = (p0_init(r) - p0_init(r-1))/dr(n)
          rhog = HALF*(s0_init(r,rho_comp) + s0_init(r-1,rho_comp))*g

          if (print_init_hse_diag) then
             if ( parallel_IOProcessor() ) then
                print *, 'r, dpdr, rhog, err: ', rloc, dpdr, rhog, &
                     abs(dpdr - rhog)/abs(rhog)
             endif
          endif
          
          max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(rhog))

       end if

    enddo

    if ( parallel_IOProcessor() ) then
       call log(' ')
       call log('Maximum HSE Error = ', max_hse_error)
       call log('   (after putting initial model into base state arrays, and')
       call log('    for density < base_cutoff_density)')
       write (*,887)
       call log(' ')
    endif


    ! initialize any inlet BC parameters
    call set_inlet_bcs()

    call destroy(bpt)

  end subroutine init_base_state

end module base_state_module
