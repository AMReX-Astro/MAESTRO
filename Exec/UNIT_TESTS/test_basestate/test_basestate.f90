module test_basestate_module

  use bl_types
  use bl_error_module
  use bl_constants_module
  implicit none

  public get_heating

contains

  subroutine get_heating(Hbar,s0,time,dt)
  
    use geometry, ONLY : nr, spherical, r_cc_loc
    use probin_module, ONLY : heating_time, heating_rad, heating_peak, &
         heating_sigma, prob_type, &
         cooling_rad, cooling_peak, cooling_sigma
    use variables
    use network

    real(dp_t), intent(inout) :: Hbar(0:)
    real(dp_t), intent(in   ) :: s0(0:,:)
    real(dp_t), intent(in   ) :: time,dt

    real(dp_t) :: fac
    integer :: r

    integer         :: h1_comp
    integer         :: he4_comp               
    integer         :: c12_comp
    integer         :: n14_comp
    integer         :: o16_comp
    real(kind=dp_t) :: rho, T_6_third, X_CNO, X_1, g14
    real(kind=dp_t) :: tmp1, tmp2, tmp3


    Hbar(:) = 0.d0
    
    if (prob_type .eq. 1) then

       if (time .le. heating_time) then

          if ( (time+dt) .gt. heating_time ) then
             fac = (heating_time - time) / dt
          else
             fac = 1.d0
          end if

          do r = 0, nr(1)-1
             if (spherical .eq. 0) then
                ! plane-parallel -- do the heating term in paper II (section 4)
                Hbar(r) = fac * heating_peak * &
                     exp(-((r_cc_loc(1,r) - heating_rad)**2)/ heating_sigma)
             else
                ! spherical -- lower amplitude heating term
                Hbar(r) = fac * heating_peak * &
                     exp(-((r_cc_loc(1,r) - heating_rad)**2)/ heating_sigma)
             endif
          enddo
       end if

    elseif (prob_type .eq. 2) then

       ! analytic heating modeling CNO cycle

       h1_comp = spec_comp - 1 + network_species_index("hydrogen-1")
       c12_comp = spec_comp - 1 + network_species_index("carbon-12")
       n14_comp = spec_comp - 1 + network_species_index("nitrogen-14")
       o16_comp = spec_comp - 1 + network_species_index("oxygen-16")

       do r = 0, nr(1)-1
          rho = s0(r,rho_comp)
          T_6_third = (s0(r,temp_comp) / 1.0d6) ** THIRD
          tmp1 = s0(r,c12_comp)
          tmp2 = s0(r,n14_comp)
          tmp3 = s0(r,o16_comp)
          X_CNO = (tmp1 + tmp2 + tmp3) / rho
          X_1 = s0(r,h1_comp) / rho
          tmp1 =   2.7d-3 * T_6_third
          tmp2 = -7.78d-3 * T_6_third**2
          tmp3 = -1.49d-4 * T_6_third**3
          g14 = 1.0_dp_t + tmp1 + tmp2 + tmp3
          tmp1 = 8.67d27 * g14 * X_CNO * X_1 * rho / T_6_third**2
          tmp2 = dexp(-1.5228d2 / T_6_third)
          Hbar(r) = tmp1 * tmp2
       enddo

    elseif (prob_type .eq. 3) then

       if (.not. (network_name == "triple_alpha" .or. &
                  network_name == "triple_alpha_plus_cago" .or. &
                  network_name == "general-triple_alpha_plus_o.net")) then
          call bl_error("ERROR: invalid network for prob_type == 3")
       endif

       he4_comp = spec_comp - 1 + network_species_index("helium-4")


       ! off-center heating for sub_chandra
       if (time .le. heating_time) then

          if ( (time+dt) .gt. heating_time ) then
             fac = (heating_time - time) / dt
          else
             fac = 1.d0
          end if

          do r = 0, nr(1)-1
             if (spherical .eq. 0) then
                call bl_error("ERROR: heating not supported")
             else
                ! spherical -- lower amplitude heating term
                Hbar(r) = fac * heating_peak * &
                     exp(-((r_cc_loc(1,r) - heating_rad)**2)/ heating_sigma**2)

                ! only heat if there is He-4
                Hbar(r) = (s0(r,he4_comp)/s0(r,rho_comp)) * Hbar(r)
             endif
          enddo
       end if

    elseif (prob_type .eq. 4) then
       ! Apply both heating and cooling for an Urca process

       if (time .le. heating_time) then

          if ( (time+dt) .gt. heating_time ) then
             fac = (heating_time - time) / dt
          else
             fac = 1.d0
          end if

          do r = 0, nr(1)-1
             if (spherical .eq. 0) then
                ! plane-parallel -- do the heating term in paper II (section 4)
                ! plus a similar cooling term for Urca
                Hbar(r) = fac * (&
                     heating_peak * exp(-((r_cc_loc(1,r) - heating_rad)**2)/ heating_sigma) + &
                     cooling_peak * exp(-((r_cc_loc(1,r) - cooling_rad)**2)/ cooling_sigma))
                
             else
                ! spherical -- lower amplitude heating/cooling term
                Hbar(r) = fac * (&
                     heating_peak * exp(-((r_cc_loc(1,r) - heating_rad)**2)/ heating_sigma) + &
                     cooling_peak * exp(-((r_cc_loc(1,r) - cooling_rad)**2)/ cooling_sigma))
             endif
          enddo
       end if
       
    else

       write(*,*) prob_type
       call bl_error("prob_type not yet supported.")       

    endif

    return
  end subroutine get_heating

  subroutine make_Sbar(Sbar, s0, Hbar)

    use variables
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    use network, only: nspec
    use geometry, only : nr

    real(dp_t), intent(inout) :: Sbar(0:)
    real(dp_t), intent(in   ) :: s0(0:,:)
    real(dp_t), intent(in   ) :: Hbar(0:)

    integer :: r
    type (eos_t) :: eos_state

     do r=0,nr(1)-1

        ! (rho, T) --> p,h, etc
        eos_state%rho   = s0(r,rho_comp)
        eos_state%T     = s0(r,temp_comp)
        eos_state%xn(:) = s0(r,spec_comp:spec_comp-1+nspec)/s0(r,rho_comp)

        call eos(eos_input_rt, eos_state)

        Sbar(r) = Hbar(r) * eos_state%dpdt / (eos_state%rho * eos_state%cp * eos_state%dpdr)

     enddo

  end subroutine make_Sbar


end module test_basestate_module

