module check_cutoff_module

  implicit none

contains

  subroutine check_cutoff_values()

    use parallel
    use sponge_module, only: sponge_start_density
    use probin_module, only: buoyancy_cutoff_factor, base_cutoff_density, anelastic_cutoff, &
         small_dens, small_temp
    use model_parser_module
    use bl_error_module

    real (kind=dp_t) :: min_dens, max_dens
    real (kind=dp_t) :: min_temp, max_temp

    real(kind=dp_t), parameter :: eps = 1.d-8


887 format(78('-'))
888 format(a60,g18.10)
889 format(a60)


    ! compute the extrema of the initial model
    max_dens = maxval(model_state(:,idens_model))
    min_dens = minval(model_state(:,idens_model))
    
    max_temp = maxval(model_state(:,itemp_model))
    min_temp = minval(model_state(:,itemp_model))

    if ( parallel_IOProcessor() ) then
       write (*,889) ' '
       write (*,887)
       write (*,*)   'initial model extrema: '
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


    ! report on the various cutoffs
    if ( parallel_IOProcessor() ) then
          ! output block for cutoff density information
       write (*,887)
       write (*,*)   'physical cutoff densities:'
       write (*,888) '    low density cutoff (for mapping the model) =      ', &
            base_cutoff_density
       write (*,888) '    buoyancy cutoff density                           '
       write (*,888) '        (for zeroing rho - rho_0, centrifugal term) = ', &
            buoyancy_cutoff_factor*base_cutoff_density
       write (*,888) '    anelastic cutoff =                                ', &
            anelastic_cutoff
       write (*,888) '    sponge start density =                            ', &
            sponge_start_density
       write (*,888) ' '

       write (*,*)   'thermodynamics cutoffs:'
       write (*,888) '    EOS temperature floor =                           ', &
            small_temp
       write (*,888) '    EOS density floor =                               ', &
            small_dens
    end if


    ! some cutoff sanity checks
    if (min_dens < base_cutoff_density .OR. min_dens < anelastic_cutoff) then
       if ( parallel_IOProcessor() ) then
          print *, ' '
          print *, 'WARNING: minimum model density is lower than one of the cutoff densities'
          print *, '         make sure that the cutoff densities are lower than any density'
          print *, '         of dynamical interest'
       end if
    endif
       
    if (min_dens + eps > base_cutoff_density .or. min_dens + eps > anelastic_cutoff) then
       if ( parallel_IOProcessor() ) then
          print *, ' '
          print *, 'WARNING: minimum model density is larger than, or very close to '
          print *,'          the cutoff density.' 
       end if
    end if

    if (buoyancy_cutoff_factor*base_cutoff_density > sponge_start_density) then
       if ( parallel_IOProcessor() ) then
          print *, ' '
          print *, 'WARNING: buoyancy cutoff occurs at densities greater than those sponged'
          print *, '         The buoyancy_cutoff_factor should be lowered to ensure the'
          print *, '         buoyancy cutoff occurs within the sponged region'
       endif
    endif

    if (min_temp < small_temp) then
       if ( parallel_IOProcessor() ) then
          print *, ' '
          print *, 'WARNING: minimum model temperature is lower than the EOS cutoff'
          print *, '         temperature, small_temp'
       endif
    endif

    if (min_dens < small_dens) then
       if ( parallel_IOProcessor() ) then
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

    if ( parallel_IOProcessor() ) then
       ! close the cutoff density output block
       write (*,887)
       write (*,*)   ' '
    end if

  end subroutine check_cutoff_values

end module check_cutoff_module
