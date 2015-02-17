! check_cutoff_values performs some simple checks at the start of a
! run to ensure that the various cutoff parameters
! (base_cutoff_density, anelastic_cutoff, buoyancy_cutoff_factor, ...)
! are reasonable.

module check_cutoff_module

  use bl_types, only: dp_t
  implicit none

contains

  subroutine check_cutoff_values()

    use parallel, only: parallel_IOProcessor
    use sponge_module, only: sponge_start_density
    use probin_module, only: buoyancy_cutoff_factor, base_cutoff_density, &
                             anelastic_cutoff, small_dens, small_temp, do_sponge
    use model_parser_module, only: model_initialized, model_state, &
                                   idens_model, itemp_model
    use bl_error_module, only: bl_error
    use simple_log_module

    real (kind=dp_t) :: model_min_dens, model_max_dens
    real (kind=dp_t) :: model_min_temp, model_max_temp

    real(kind=dp_t), parameter :: eps = 1.d-8


887 format(78('-'))
888 format(a60,g18.10)
889 format(a60)


    ! compute the extrema of the initial model -- if we read one in
    if (model_initialized) then
       model_max_dens = maxval(model_state(:,idens_model))
       model_min_dens = minval(model_state(:,idens_model))
       
       model_max_temp = maxval(model_state(:,itemp_model))
       model_min_temp = minval(model_state(:,itemp_model))

       if ( parallel_IOProcessor() ) then
          write (*,889) ' '
          write (*,887)
          write (*,*)   'initial model extrema: '
          write (*,888) '    minimum density of model =                        ', &
               model_min_dens
          write (*,888) '    maximum density of model =                        ', &
               model_max_dens
          write (*,*)   ' '
          write (*,888) '    minimum temperature of model =                    ', &
               model_min_temp
          write (*,888) '    maximum temperature of model =                    ', &
               model_max_temp
          write (*,887)
          write (*,889) ' '
       endif
    endif


    ! report on the various cutoffs
    if ( parallel_IOProcessor() ) then
          ! output block for cutoff density information
       write (*,887)
       call log('physical cutoff densities:')
       call log('    low density cutoff (for mapping the model) =      ', base_cutoff_density)
       call log('    buoyancy cutoff density                           ')
       call log('        (for zeroing rho - rho_0, centrifugal term) = ', &
            buoyancy_cutoff_factor*base_cutoff_density)
       call log('    anelastic cutoff =                                ', anelastic_cutoff)
       if (do_sponge) then
          call log('    sponge start density =                            ', sponge_start_density)
       endif
       call log(' ')
       call log('thermodynamics cutoffs:')
       call log('    EOS temperature floor =                           ', small_temp)
       call log('    EOS density floor =                               ', small_dens)
    end if

    
    if (model_initialized) then

       ! some cutoff sanity checks
       if ( model_min_dens < base_cutoff_density .OR. &
            model_min_dens < anelastic_cutoff) then
          if ( parallel_IOProcessor() ) then
             call log(' ')
             call log('WARNING: minimum model density is lower than one of the cutoff densities')
             call log('         make sure that the cutoff densities are lower than any density')
             call log('         of dynamical interest')
          end if
       endif
       
       if ( model_min_dens + eps > base_cutoff_density .or. &
            model_min_dens + eps > anelastic_cutoff) then
          if ( parallel_IOProcessor() ) then
             call log(' ')
             call log('WARNING: minimum model density is larger than, or very close to ')
             call log('          the cutoff density.')
          end if
       end if

       if (model_min_temp < small_temp) then
          if ( parallel_IOProcessor() ) then
             call log(' ')
             call log('WARNING: minimum model temperature is lower than the EOS cutoff')
             call log('         temperature, small_temp')
          endif
       endif

       if (model_min_dens < small_dens) then
          if ( parallel_IOProcessor() ) then
             call log(' ')
             call log('WARNING: minimum model density is lower than the EOS cutoff')
             call log('         density, small_dens')
          endif
       endif

    endif

       
    if (anelastic_cutoff < base_cutoff_density) then
       print *, 'ERROR: anelastic cutoff should be at a higher density than the base state'
       print *, '       cutoff density.'
       call bl_error("anelastic cutoff < base_cutoff_density")
    endif


    if (do_sponge) then
       if (buoyancy_cutoff_factor*base_cutoff_density > sponge_start_density) then
          if ( parallel_IOProcessor() ) then
             call log(' ')
             call log('WARNING: buoyancy cutoff occurs at densities greater than those sponged')
             call log('         The buoyancy_cutoff_factor should be lowered to ensure the')
             call log('         buoyancy cutoff occurs within the sponged region')
          endif
       endif
    endif
    

    if ( parallel_IOProcessor() ) then
       ! close the cutoff density output block
       write (*,887)
       call log(' ')
    end if


  end subroutine check_cutoff_values

end module check_cutoff_module
