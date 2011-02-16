module base_state_module

  use bl_types
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

    use bl_prof_module
    use parallel
    use bl_error_module
    use bl_constants_module
    use eos_module
    use network, only: spec_names
    use probin_module, only: base_cutoff_density, anelastic_cutoff, prob_lo, prob_hi, &
                             buoyancy_cutoff_factor, sponge_start_factor, &
                             sponge_center_density, &
                             small_temp, small_dens, grav_const, planar_invsq_mass, &
                             do_planar_invsq_grav
    use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp, trac_comp, ntrac
    use geometry, only: dr, spherical, nr
    use inlet_bc_module, only: set_inlet_bcs
    use fundamental_constants_module, only: Gconst
    use model_parser_module
    
    integer           , intent(in   ) :: n
    character(len=256), intent(in   ) :: model_file
    real(kind=dp_t)   , intent(inout) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(inout) :: p0_init(0:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    ! local
    logical, save   :: firstCall = .true.
    integer         :: r,comp
    real(kind=dp_t) :: rloc,dr_in,rmax,starting_rad,mod_dr
    real(kind=dp_t) :: d_ambient,t_ambient,p_ambient,xn_ambient(nspec)
    real(kind=dp_t) :: sumX

    real(kind=dp_t) :: min_dens, max_dens, min_temp, max_temp, eps

    type(bl_prof_timer), save :: bpt
    
    real(kind=dp_t), parameter :: TINY = 1.0d-10

    real(kind=dp_t) :: mencl, g, r_l, r_r, dpdr, rhog
    real(kind=dp_t) :: max_hse_error, max_hse_err_r

    call build(bpt, "init_base_state")

887 format(78('-'))
888 format(a60,g18.10)
889 format(a60)


    if (firstCall) then

       ! only need to read in the initial model once -- model_parser_module
       ! stores the model data
       call read_model_file(model_file)

       firstCall = .false.
    endif

    eps = 1.d-8

    max_dens = maxval(model_state(:,idens_model))
    min_dens = minval(model_state(:,idens_model))

    max_temp = maxval(model_state(:,itemp_model))
    min_temp = minval(model_state(:,itemp_model))

    if ( parallel_IOProcessor() .and. n == 1) then
       write (*,889) ' '
       write (*,888) '    minimum density of model =                        ', &
            min_dens
       write (*,888) '    maximum density of model =                        ', &
            max_dens
       write (*,*)   ' '
       write (*,888) '    minimum temperature of model =                    ', &
            min_temp
       write (*,888) '    maximum temperature of model =                    ', &
            max_temp
       write (*,887)
       write (*,889) ' '
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
       write (*,888) '    sponge start density =                            ', &
            sponge_start_factor*sponge_center_density
       write (*,888) '    sponge center density =                           ', &
            sponge_center_density
       write (*,888) ' '
    end if

    if (min_dens < base_cutoff_density .OR. min_dens < anelastic_cutoff) then
       if ( parallel_IOProcessor() .and. n == 1) then
          print *, ' '
          print *, 'WARNING: minimum model density is lower than one of the cutoff densities'
          print *, '         make sure that the cutoff densities are lower than any density'
          print *, '         of dynamical interest'
       end if
    endif

    if (min_dens + eps > base_cutoff_density .or. min_dens + eps > anelastic_cutoff) then
       if ( parallel_IOProcessor() .and. n == 1) then
          print *, ' '
          print *, 'WARNING: minimum model density is larger than, or very close to '
          print *,'          the cutoff density.' 
       end if
    end if

    if (min_temp < small_temp) then
       if ( parallel_IOProcessor() .and. n == 1) then
          print *, ' '
          print *, 'WARNING: minimum model temperature is lower than the EOS cutoff'
          print *, '         temperature, small_temp'
       endif
    endif

    if (min_dens < small_dens) then
       if ( parallel_IOProcessor() .and. n == 1) then
          print *, ' '
          print *, 'WARNING: minimum model density is lower than the EOS cutoff'
          print *, '         density, small_dens'
       endif
    endif


    if (anelastic_cutoff < base_cutoff_density) then
       print *, 'ERROR: anelastic cutoff should be at a higher density than the base state'
       print *, '       cutoff density.'
       call bl_error("anelastic cutoff < base_cutoff_density")
    endif

    if ( parallel_IOProcessor() .and. n == 1) then
       ! close the cutoff density output block
       write (*,887)
       write (*,*)   ' '
    end if

    dr_in = (model_r(npts_model) - model_r(1)) / dble(npts_model-1)
    rmax = model_r(npts_model)

    if ( parallel_IOProcessor() ) then
       write (*,887)
       if (spherical .ne. 1) then
          write (*,*)   'model file mapping, level:', n
       else
          write (*,*)   'model file mapping (spherical base state)'
       endif

       write (*,888) 'dr of MAESTRO base state =                            ', &
            dr(n)
       write (*,888) 'dr of input file data =                               ', &
            dr_in
       write (*,*) ' '
       write (*,888) 'maximum radius (cell-centered) of input model =       ', &
            rmax
       
       if (dr(n) .lt. dr_in) then
          mod_dr = mod(dr_in,dr(n))
       else
          mod_dr = mod(dr(n),dr_in)
       endif

       if (mod_dr .gt. TINY) then
          print *, ''
          print *, "WARNING: resolution of base state array is not an integer"
          print *, "         multiple of the initial model's resolution.     "
          print *, "         make sure this is a desired property as this    "
          print *, "         could lead to aliasing when performing the      "
          print *, "         interpolation.                                  "
          print *, " "
          print *, "modulus = ", mod_dr
       endif
    end if

    if (spherical .eq. 0) then
       starting_rad = prob_lo(size(dx))
    else
       starting_rad = ZERO
    endif

    if (n .eq. 1) then
       base_cutoff_density_loc = 1.d99
    end if
    
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

          d_ambient = interpolate(rloc, npts_model, model_r, model_state(:,idens_model))
          t_ambient = interpolate(rloc, npts_model, model_r, model_state(:,itemp_model))
          p_ambient = interpolate(rloc, npts_model, model_r, model_state(:,ipres_model))

          sumX = ZERO
          do comp = 1, nspec
             xn_ambient(comp) = max(ZERO,min(ONE, &
                  interpolate(rloc, npts_model, model_r, model_state(:,ispec_model-1+comp))))
             sumX = sumX + xn_ambient(comp)
          enddo
          xn_ambient = xn_ambient/sumX

          ! use the EOS to make the state consistent
          temp_eos(1) = t_ambient
          den_eos(1)  = d_ambient
          p_eos(1)    = p_ambient
          xn_eos(1,:) = xn_ambient(:)

          write(10,*) rloc, d_ambient, t_ambient

          ! (rho,T) --> p,h
          call eos(eos_input_rt, den_eos, temp_eos, &
                   npts, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, &
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   .false.)

          s0_init(r, rho_comp ) = d_ambient
          s0_init(r,rhoh_comp ) = d_ambient * h_eos(1)
          s0_init(r,spec_comp:spec_comp+nspec-1) = d_ambient * xn_ambient(1:nspec)
          p0_init(r) = p_eos(1)

          s0_init(r,temp_comp) = t_ambient

          write(11,*) rloc, s0_init(r,rhoh_comp)

          if (ntrac .gt. 0) then
             s0_init(r,trac_comp:trac_comp+ntrac-1) = ZERO
          end if

          ! keep track of the height where we drop below the cutoff density
          if (s0_init(r,rho_comp) .le. base_cutoff_density .and. &
               base_cutoff_density_loc .eq. 1.d99 .and. n .eq. 1) then

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

! FIXME!!!!!
       open(unit=11,file='hse_error.dat',form = "formatted", access = "sequential",action="write")
!!!!!!!!!!!!

    mencl = zero

! HACK   
!    if (spherical .eq. 1) then
       mencl = four3rd*m_pi*dr(n)**3*s0_init(0,rho_comp)
!    endif

    max_hse_error = -1.d30

    do r=1,nr(n)-1
       
       rloc = starting_rad + (dble(r) + HALF)*dr(n)
       rloc = min(rloc, rmax)

       if (rloc .lt. base_cutoff_density_loc) then

          r_r = starting_rad + dble(r+1)*dr(n)
          r_l = starting_rad + dble(r)*dr(n)

          if (spherical .eq. 1) then
             g = -Gconst*mencl/r_l**2
             mencl = mencl &
                  + four3rd*m_pi*dr(n)*(r_l**2+r_l*r_r+r_r**2)*s0_init(r,rho_comp)
          else
! FIXME! HACK to map octant into a planar geometry 
             g = -Gconst*mencl/r_l**2
             mencl = mencl &
                  + four3rd*m_pi*dr(n)*(r_l**2+r_l*r_r+r_r**2)*s0_init(r,rho_comp)
          endif

          dpdr = (p0_init(r) - p0_init(r-1))/dr(n)
          rhog = HALF*(s0_init(r,rho_comp) + s0_init(r-1,rho_comp))*g

! FIXME !!!!!!!!!!!!
          write(11,1000)r,(r+HALF)*dr(n), abs(dpdr - rhog)/abs(rhog), &
                        dpdr, rhog 
1000 format(i8,32(e30.20,1x))
!!!!!!!!!!!!!!!!!

!          max_hse_error = max(max_hse_error, abs(dpdr - rhog)/abs(rhog))
          if ( abs(dpdr - rhog)/abs(rhog) .gt. max_hse_error ) then !&
!               .and. rloc .lt. 8.15d10) then
             max_hse_error = abs(dpdr - rhog)/abs(rhog)
             max_hse_err_r = rloc
          end if

       end if

    enddo

! FIXME
    close(11)
!!!!!

    if ( parallel_IOProcessor() ) then
       write (*,*) ' '
       write (*,*) 'Maximum HSE Error = ', max_hse_error
       write (*,*) '        at radius = ', max_hse_err_r
       write (*,*) '   (after putting initial model into base state arrays, and'
       write (*,*) '    for density < base_cutoff_density)'
       write (*,887)
       write (*,*) ' '
    endif


    ! initialize any inlet BC parameters
    call set_inlet_bcs()

    call destroy(bpt)

  end subroutine init_base_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function interpolate(r, npts, model_r, model_var)

    use bl_constants_module

    ! given the array of model coordinates (model_r), and variable (model_var),
    ! find the value of model_var at point r using linear interpolation.
    ! Eventually, we can do something fancier here.

    real(kind=dp_t) :: interpolate
    real(kind=dp_t), intent(in) :: r
    integer :: npts
    real(kind=dp_t), dimension(npts) :: model_r, model_var

    real(kind=dp_t) :: slope
    real(kind=dp_t) :: minvar, maxvar

    integer :: i, id

    ! find the location in the coordinate array where we want to interpolate
    do i = 1, npts
       if (model_r(i) >= r) exit
    enddo
    if(i .gt. 1 .and. i .lt. npts+1) then
       if(abs(r-model_r(i-1)) .lt. abs(r-model_r(i))) then
          i = i-1
       end if
    end if
    if (i == npts+1) then
       i = npts
    end if

    id = i

    if (id == 1) then

       slope = (model_var(id+1) - model_var(id))/(model_r(id+1) - model_r(id))
       interpolate = slope*(r - model_r(id)) + model_var(id)

       ! safety check to make sure interpolate lies within the bounding points
       minvar = min(model_var(id+1),model_var(id))
       maxvar = max(model_var(id+1),model_var(id))
       interpolate = max(interpolate,minvar)
       interpolate = min(interpolate,maxvar)

    else if (id == npts) then

       slope = (model_var(id) - model_var(id-1))/(model_r(id) - model_r(id-1))
       interpolate = slope*(r - model_r(id)) + model_var(id)

       ! safety check to make sure interpolate lies within the bounding points
       minvar = min(model_var(id),model_var(id-1))
       maxvar = max(model_var(id),model_var(id-1))
       interpolate = max(interpolate,minvar)
       interpolate = min(interpolate,maxvar)

    else

       if (r .ge. model_r(id)) then

          ! we should not wind up in here

          slope = (model_var(id+1) - model_var(id))/(model_r(id+1) - model_r(id))
          interpolate = slope*(r - model_r(id)) + model_var(id)
          
          ! safety check to make sure interpolate lies within the bounding points
          minvar = min(model_var(id+1),model_var(id))
          maxvar = max(model_var(id+1),model_var(id))
          interpolate = max(interpolate,minvar)
          interpolate = min(interpolate,maxvar)
          
       else

          slope = (model_var(id) - model_var(id-1))/(model_r(id) - model_r(id-1))
          interpolate = slope*(r - model_r(id)) + model_var(id)
          
          ! safety check to make sure interpolate lies within the bounding points
          minvar = min(model_var(id),model_var(id-1))
          maxvar = max(model_var(id),model_var(id-1))
          interpolate = max(interpolate,minvar)
          interpolate = min(interpolate,maxvar)
          
       end if

    end if

    return

  end function interpolate

end module base_state_module
