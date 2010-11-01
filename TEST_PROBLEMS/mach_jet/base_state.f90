module base_state_module

  use bl_types
  use network, only: nspec
  
  implicit none

  private

  public :: init_base_state

contains

  subroutine init_base_state(n,model_file,s0_init,p0_init,dx)

    use bl_prof_module
    use parallel
    use bl_error_module
    use bl_constants_module
    use eos_module
    use network, only: spec_names, network_species_index
    use probin_module, only: grav_const, do_stratified, do_isentropic
    use variables, only: rho_comp, rhoh_comp, temp_comp, spec_comp, trac_comp, ntrac
    use geometry, only: dr, spherical, nr
    use inlet_bc_module
    
    integer           , intent(in   ) :: n
    character(len=256), intent(in   ) :: model_file
    real(kind=dp_t)   , intent(inout) :: s0_init(0:,:)
    real(kind=dp_t)   , intent(inout) :: p0_init(0:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    ! local
    integer         :: r,j
    real(kind=dp_t) :: z,H

    type(bl_prof_timer), save :: bpt

    call build(bpt, "init_base_state")

    if (spherical .eq. 1) then
       call bl_error("ERROR: base_state.f90 is not valid for spherical")
    endif

    if (do_stratified) then

       ! use the EOS to make the state consistent
       den_eos(1) = 1.d-3
       p_eos(1) = 1.d6

       ! only initialize the first species
       xn_eos(1,:) = ZERO
       xn_eos(1,1) = 1.d0

       p0_init(0) = p_eos(1)

       ! H = pres_base / dens_base / abs(grav_const)
       H = 1.d6 / 1.d-3 / abs(grav_const)

       ! set an initial guess for the temperature -- this will be reset
       ! by the EOS
       temp_eos(1) = 10.d0

       do j=0,nr(n)-1

          z = (dble(j)+HALF) * dr(1)

          if (do_isentropic) then
             den_eos(1) = 1.d-3*(grav_const*1.d-3*(gamma_const - 1.0)*z/ &
               (gamma_const*1.d6) + 1.d0)**(1.d0/(gamma_const - 1.d0))
          else
             den_eos(1) = 1.d-3*exp(-z/H)
          end if

          s0_init(j, rho_comp) = den_eos(1)

          if (j.eq.0) then
             p0_init(j) = p0_init(j) - &
                  dr(1) * HALF * s0_init(j,rho_comp) * &
                  abs(grav_const)
          else if (j.gt.0) then
             p0_init(j) = p0_init(j-1) - &
                  dr(1) * HALF * (s0_init(j,rho_comp) + s0_init(j-1,rho_comp)) * &
                  abs(grav_const)
          end if
             
          p_eos(1) = p0_init(j)

          ! use the EOS to make the state consistent

          ! (rho,p) --> T, h
          call eos(eos_input_rp, den_eos, temp_eos, &
                   npts, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, &
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   .false.)

          s0_init(j, rho_comp) = den_eos(1)
          s0_init(j,rhoh_comp) = den_eos(1)*h_eos(1)
          
          s0_init(j,spec_comp:spec_comp-1+nspec) = ZERO
          s0_init(j,spec_comp) = den_eos(1)
          
          s0_init(j,temp_comp) = temp_eos(1)
          s0_init(j,trac_comp) = ZERO
          
       end do

    else

       ! use the EOS to make the state consistent
       temp_eos(1) = 10.d0
       den_eos(1)  = 1.d-3
       p_eos(1)    = 1.d6
       xn_eos(1,:) = 1.d0

       ! (rho,p) --> T, h
       call eos(eos_input_rp, den_eos, temp_eos, &
                npts, &
                xn_eos, &
                p_eos, h_eos, e_eos, &
                cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                dpdX_eos, dhdX_eos, &
                gam1_eos, cs_eos, s_eos, &
                dsdt_eos, dsdr_eos, &
                .false.)

       s0_init(0:nr(n)-1, rho_comp) = den_eos(1)
       s0_init(0:nr(n)-1,rhoh_comp) = den_eos(1)*h_eos(1)
       s0_init(0:nr(n)-1,spec_comp) = den_eos(1)
       s0_init(0:nr(n)-1,temp_comp) = temp_eos(1)
       s0_init(0:nr(n)-1,trac_comp) = ZERO
    
       p0_init(0:nr(n)-1) = p_eos(1)

    end if

    call set_inlet_bcs()

    call destroy(bpt)

  end subroutine init_base_state

end module base_state_module
