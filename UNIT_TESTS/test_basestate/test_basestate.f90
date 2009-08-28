module test_basestate_module

  use bl_types
  implicit none

  public get_heating

contains

  subroutine get_heating(Hbar,time,dt)
  
    use geometry, ONLY : nr, spherical, r_cc_loc

    real(dp_t), intent(inout) :: Hbar(0:)
    real(dp_t), intent(in   ) :: time,dt

    real(dp_t) :: y_0, fac, t_stop
    integer :: r

!   For heating away from the center
    y_0 = 5.d7

!   For heating at the center
!   y_0 = 0.d0

    Hbar(:) = 0.d0

    t_stop = 0.5d0
    
    if (time .le. t_stop) then

       if ( (time+dt) .gt. t_stop ) then
           fac = (t_stop - time) / dt
       else
           fac = 1.d0
       end if
 
       do r = 0, nr(1)-1
          if (spherical .eq. 0) then
             ! plane-parallel -- do the heating term in paper II (section 4)
             Hbar(r) = 1.d17 * exp(-((r_cc_loc(1,r) - y_0)**2)/ 1.d14)
          else
             ! spherical -- lower amplitude heating term
             Hbar(r) = fac * 1.0d16 * exp(-((r_cc_loc(1,r) - y_0)**2)/ 1.d14)
          endif
       enddo
    end if

    return
  end subroutine get_heating

  subroutine make_Sbar(Sbar, s0, Hbar)

    use variables
    use eos_module
    use network, only: nspec
    use geometry, only : nr

    real(dp_t), intent(inout) :: Sbar(0:)
    real(dp_t), intent(in   ) :: s0(0:,:)
    real(dp_t), intent(in   ) :: Hbar(0:)

    integer :: r

     do r=0,nr(1)-1

        ! (rho, T) --> p,h, etc
        den_eos(1)  = s0(r,rho_comp)
        temp_eos(1) = s0(r,temp_comp)
        xn_eos(1,:) = s0(r,spec_comp:spec_comp-1+nspec)/s0(r,rho_comp)

        call eos(eos_input_rt, den_eos, temp_eos, NP, nspec, &
                 xn_eos, &
                 p_eos, h_eos, e_eos, &
                 cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                 dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                 dpdX_eos, dhdX_eos, &
                 gam1_eos, cs_eos, s_eos, &
                 dsdt_eos, dsdr_eos, &
                 do_diag)

        Sbar(r) = Hbar(r) * dpdt_eos(1) / (den_eos(1) * cp_eos(1) * dpdr_eos(1))

     enddo

  end subroutine make_Sbar


end module test_basestate_module

