!
! A module to provide integer indices into the various storage arrays
! for accessing the different variables by name.
!
module variables

  use bl_types

  implicit none

  integer, save :: rho_comp, rhoh_comp, spec_comp, temp_comp, trac_comp, press_comp
  integer, save :: foextrap_comp, hoextrap_comp
  integer, save :: icomp_vel, icomp_rho, icomp_rhoh, icomp_h, icomp_spec, icomp_trac
  integer, save :: icomp_w0, icomp_divw0, icomp_rho0, icomp_rhoh0, icomp_h0
  integer, save :: icomp_p0, icomp_velr, icomp_velc
  integer, save :: icomp_magvel, icomp_mom, icomp_vort, icomp_src
  integer, save :: icomp_tfromp,icomp_tpert,icomp_rhopert,icomp_rhohpert
  integer, save :: icomp_machno,icomp_cs
  integer, save :: icomp_dg,icomp_pi,icomp_gpi,icomp_pioverp0,icomp_p0pluspi
  integer, save :: icomp_entropy,icomp_entropypert
  integer, save :: icomp_tfromH,icomp_dp,icomp_dT
  integer, save :: icomp_omegadot,icomp_enuc,icomp_Hext, icomp_eta, icomp_sponge
  integer, save :: icomp_thermal, icomp_conductivity
  integer, save :: icomp_ad_excess
  integer, save :: icomp_part
  integer, save :: icomp_proc
  integer, save :: icomp_pidivu

  ! the total number of plot components
  integer, save :: n_plot_comps = 0

  integer, save :: ntrac,nscal
  real(kind=dp_t), save :: rel_eps

contains

  function get_next_plot_index(num) result (next)

    ! return the next starting index for a plotfile quantity,
    ! and increment the counter of plotfile quantities by num
    integer :: num, next

    next = n_plot_comps + 1
    n_plot_comps = n_plot_comps + num

    return
  end function get_next_plot_index

  subroutine init_variables()

    use probin_module, only: dm_in
    use network, only: nspec

    rho_comp    = 1
    rhoh_comp   = 2
    spec_comp   = rhoh_comp + 1
    temp_comp   = spec_comp + nspec
    trac_comp   = temp_comp + 1

    ntrac = 1

    ! The "3" here refers to rho, rhoh, and temp
    nscal = nspec + ntrac + 3

    press_comp  = dm_in + nscal + 1

    foextrap_comp = press_comp + 1
    hoextrap_comp = foextrap_comp + 1

  end subroutine init_variables

  subroutine init_plot_variables()

    use network, only: nspec
    use probin_module, only: plot_spec, plot_trac, plot_base, use_thermal_diffusion, &
         plot_omegadot, plot_Hnuc, plot_Hext, plot_eta, plot_ad_excess, &
         use_tfromp, plot_h_with_use_tfromp, plot_gpi, plot_cs, dm_in, use_particles, &
         plot_processors, plot_pidivu
    use geometry, only: spherical

    icomp_vel      = get_next_plot_index(dm_in)
    icomp_rho      = get_next_plot_index(1)
    if (.not. use_tfromp .or. (use_tfromp .and. plot_h_with_use_tfromp)) then
       icomp_rhoh     = get_next_plot_index(1)
       icomp_h        = get_next_plot_index(1)
    end if

    if (plot_spec) icomp_spec = get_next_plot_index(nspec)
    if (plot_trac) icomp_trac = get_next_plot_index(ntrac)

    if (plot_base) then
       icomp_w0    = get_next_plot_index(dm_in)
       icomp_divw0 = get_next_plot_index(1)
       icomp_rho0  = get_next_plot_index(1)
       icomp_rhoh0 = get_next_plot_index(1)
       icomp_h0    = get_next_plot_index(1)
       icomp_p0    = get_next_plot_index(1)
    end if

    if (spherical .eq. 1) then
       icomp_velr = get_next_plot_index(1)
       icomp_velc = get_next_plot_index(1)
    end if

    icomp_magvel      = get_next_plot_index(1)
    icomp_mom         = get_next_plot_index(1)
    icomp_vort        = get_next_plot_index(1)
    icomp_src         = get_next_plot_index(1)
    icomp_rhopert     = get_next_plot_index(1)

    if (.not. use_tfromp .or. (use_tfromp .and. plot_h_with_use_tfromp)) then
       icomp_rhohpert    = get_next_plot_index(1)
    endif

    icomp_tfromp      = get_next_plot_index(1)
    if (.not. use_tfromp .or. (use_tfromp .and. plot_h_with_use_tfromp)) then
       icomp_tfromH      = get_next_plot_index(1)
       icomp_dT          = get_next_plot_index(1)
       icomp_dp          = get_next_plot_index(1)
    endif

    icomp_tpert       = get_next_plot_index(1)
    icomp_machno      = get_next_plot_index(1)
    if (plot_cs) then
       icomp_cs          = get_next_plot_index(1)
    end if
    icomp_dg          = get_next_plot_index(1)
    icomp_entropy     = get_next_plot_index(1)
    icomp_entropypert = get_next_plot_index(1)
    icomp_sponge      = get_next_plot_index(1)

    icomp_pi          = get_next_plot_index(1)
    if (plot_gpi) then
       icomp_gpi         = get_next_plot_index(dm_in)
    endif

    if (plot_base) then
       icomp_pioverp0    = get_next_plot_index(1)
       icomp_p0pluspi    = get_next_plot_index(1)
    end if

    if (plot_omegadot) then
       icomp_omegadot = get_next_plot_index(nspec)
    end if

    if (plot_Hnuc) then
       icomp_enuc     = get_next_plot_index(1)
    end if

    if (plot_Hext) then
      icomp_Hext     = get_next_plot_index(1)
    end if

    if (plot_eta) then
      icomp_eta     = get_next_plot_index(1)
    end if

    if (use_thermal_diffusion) then
       icomp_thermal = get_next_plot_index(1)
       icomp_conductivity = get_next_plot_index(1)
    endif

    if (plot_ad_excess) then
       icomp_ad_excess = get_next_plot_index(1)
    endif

    if (use_particles) then
       icomp_part = get_next_plot_index(1)
    endif

    if (plot_processors) then
       icomp_proc = get_next_plot_index(1)
    endif

    if (plot_pidivu) then
       icomp_pidivu = get_next_plot_index(1)
    endif

  end subroutine init_plot_variables

end module variables
