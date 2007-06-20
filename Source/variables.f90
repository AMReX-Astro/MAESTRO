! a module to provide integer indices into the various storage arrays
! for accessing the different variables by name

module variables

  implicit none

  integer, save :: rho_comp, rhoh_comp, spec_comp, trac_comp, press_comp, derive_comp, derive_spec_comp
  integer, save :: icomp_vel, icomp_rho, icomp_rhoh, icomp_spec, icomp_trac
  integer, save :: icomp_magvel, icomp_vort,icomp_enthalpy,icomp_tfromrho,icomp_tpert,icomp_rhopert
  integer, save :: icomp_machno,icomp_dg,icomp_gp
  integer, save :: icomp_tfromH,icomp_dp,icomp_dT
  integer, save :: icomp_omegadot,icomp_enuc,icomp_sponge
  integer, save :: n_plot_comps

contains

  subroutine init_variables(dm, nscal, nspec)

    integer, intent(in) :: dm, nscal, nspec

    rho_comp    = 1
    rhoh_comp   = 2
    spec_comp   = rhoh_comp + 1
    trac_comp = spec_comp + nspec
    press_comp  = dm + nscal + 1

    icomp_vel      = 1
    icomp_rho      = dm+1
    icomp_rhoh     = icomp_rho +1
    icomp_spec     = icomp_rhoh+1
    icomp_trac     = icomp_spec+nspec

    derive_comp    = icomp_vel + dm + nscal

    icomp_magvel   = derive_comp
    icomp_vort     = derive_comp+1
    icomp_enthalpy = derive_comp+2
    icomp_rhopert  = derive_comp+3
    icomp_tfromrho = derive_comp+4
    icomp_tfromH   = derive_comp+5
    icomp_tpert    = derive_comp+6
    icomp_machno   = derive_comp+7
    icomp_dp       = derive_comp+8
    icomp_dg       = derive_comp+9

    icomp_gp       = derive_comp+10

    derive_spec_comp = icomp_gp + dm

    icomp_omegadot = derive_spec_comp
    icomp_enuc     = derive_spec_comp+nspec
    icomp_dT       = derive_spec_comp+nspec+1
    icomp_sponge   = icomp_dT+1

    n_plot_comps = icomp_sponge - icomp_vel + 1

  end subroutine init_variables

end module variables
