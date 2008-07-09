!
! A module to provide integer indices into the various storage arrays
! for accessing the different variables by name.
!
module variables

  implicit none

  integer, save :: rho_comp, rhoh_comp, spec_comp, temp_comp, trac_comp, press_comp
  integer, save :: foextrap_comp, hoextrap_comp
  integer, save :: icomp_vel, icomp_rho, icomp_rhoh, icomp_spec, icomp_trac
  integer, save :: icomp_w0, icomp_rho0, icomp_rhoh0, icomp_p0, icomp_velr
  integer, save :: icomp_magvel, icomp_velplusw0, icomp_mom, icomp_vort, icomp_divu
  integer, save :: icomp_enthalpy,icomp_tfromp,icomp_tpert,icomp_rhopert
  integer, save :: icomp_machno,icomp_dg,icomp_gp,icomp_entropy
  integer, save :: icomp_tfromH,icomp_dp,icomp_dT
  integer, save :: icomp_omegadot,icomp_enuc,icomp_sponge
  integer, save :: n_plot_comps
  integer, save :: ntrac,nscal

contains

  subroutine init_variables(dm, nspec)

    integer, intent(in) :: dm, nspec

    rho_comp    = 1
    rhoh_comp   = 2
    spec_comp   = rhoh_comp + 1
    temp_comp   = spec_comp + nspec
    trac_comp   = temp_comp + 1

    ntrac = 1

    ! The "3" here refers to rho, rhoh, and temp
    nscal = nspec + ntrac + 3

    press_comp  = dm + nscal + 1

    foextrap_comp = press_comp + 1
    hoextrap_comp = foextrap_comp + 1

  end subroutine init_variables

  subroutine init_plot_variables(dm)

    use network, only: nspec
    use probin_module, only: plot_spec, plot_trac, plot_base
    use geometry, only: spherical

    integer, intent(in) :: dm

    integer :: first_derive_comp

    icomp_vel      = 1
    icomp_rho      = dm+1
    icomp_rhoh     = icomp_rho +1
    first_derive_comp = icomp_rhoh + 1

    if (plot_spec) then
      icomp_spec     = icomp_rhoh+1
      first_derive_comp = first_derive_comp + nspec
    end if

    if (plot_trac) then
      icomp_trac  = first_derive_comp
      first_derive_comp = first_derive_comp + ntrac
    end if

    if (plot_base) then
       icomp_w0 = first_derive_comp
       first_derive_comp = first_derive_comp + dm

       icomp_rho0 = first_derive_comp
       first_derive_comp = first_derive_comp + 1

       icomp_rhoh0 = first_derive_comp
       first_derive_comp = first_derive_comp + 1

       icomp_p0 = first_derive_comp
       first_derive_comp = first_derive_comp + 1
    end if

    if (spherical .eq. 1) then
       icomp_velr = first_derive_comp
       first_derive_comp = first_derive_comp + 1
    endif

    icomp_velplusw0 = first_derive_comp
    icomp_magvel    = first_derive_comp+1
    icomp_mom       = first_derive_comp+2
    icomp_vort      = first_derive_comp+3
    icomp_divu      = first_derive_comp+4
    icomp_enthalpy  = first_derive_comp+5
    icomp_rhopert   = first_derive_comp+6
    icomp_tfromp    = first_derive_comp+7
    icomp_tfromH    = first_derive_comp+8
    icomp_tpert     = first_derive_comp+9
    icomp_machno    = first_derive_comp+10
    icomp_dp        = first_derive_comp+11
    icomp_dg        = first_derive_comp+12
    icomp_entropy   = first_derive_comp+13
    icomp_dT        = first_derive_comp+14
    icomp_sponge    = first_derive_comp+15
    icomp_gp        = first_derive_comp+16

    if (plot_spec) then
      icomp_omegadot = icomp_gp + dm
      icomp_enuc     = icomp_omegadot + nspec
      n_plot_comps = icomp_enuc
    else
      n_plot_comps = icomp_gp + dm - 1
    end if

  end subroutine init_plot_variables

end module variables
