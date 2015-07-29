!
! A module to provide integer indices into the various storage arrays
! for accessing the different variables by name.
!
module variables

  use bl_types

  implicit none

  integer, save :: rho_comp, rhoh_comp, spec_comp, temp_comp, pi_comp
  integer, save :: trac_comp, press_comp
  integer, save :: foextrap_comp, hoextrap_comp

  type plot_t
     integer :: icomp_vel = -1
     integer :: icomp_rho = -1
     integer :: icomp_rhoh = -1
     integer :: icomp_h = -1
     integer :: icomp_spec = -1
     integer :: icomp_trac = -1
     integer :: icomp_w0 = -1
     integer :: icomp_divw0 = -1
     integer :: icomp_rho0 = -1
     integer :: icomp_rhoh0 = -1
     integer :: icomp_h0 = -1
     integer :: icomp_p0 = -1
     integer :: icomp_velr = -1
     integer :: icomp_velc = -1
     integer :: icomp_magvel = -1
     integer :: icomp_mom = -1
     integer :: icomp_vort = -1
     integer :: icomp_src = -1
     integer :: icomp_tfromp = -1
     integer :: icomp_tpert = -1
     integer :: icomp_rhopert = -1
     integer :: icomp_rhohpert = -1
     integer :: icomp_machno = -1
     integer :: icomp_cs = -1
     integer :: icomp_dg = -1
     integer :: icomp_pi = -1
     integer :: icomp_gpi = -1
     integer :: icomp_pioverp0 = -1
     integer :: icomp_p0pluspi = -1
     integer :: icomp_entropy = -1
     integer :: icomp_entropypert = -1
     integer :: icomp_tfromH = -1
     integer :: icomp_dp = -1
     integer :: icomp_dT = -1
     integer :: icomp_omegadot = -1
     integer :: icomp_enuc = -1
     integer :: icomp_Hext = -1
     integer :: icomp_eta = -1
     integer :: icomp_sponge = -1
     integer :: icomp_thermal = -1
     integer :: icomp_conductivity = -1
     integer :: icomp_ad_excess = -1
     integer :: icomp_part = -1
     integer :: icomp_proc = -1
     integer :: icomp_pidivu = -1
     integer :: n_plot_comps = 0
  end type plot_t

  integer, save :: ntrac,nscal
  real(kind=dp_t), save :: rel_eps

contains

  function get_next_plot_index(num, n_plot_comps) result (next)

    ! return the next starting index for a plotfile quantity, and
    ! increment the counter of plotfile quantities, n_plot_comps, by
    ! num
    integer :: num, next
    integer, intent(inout) :: n_plot_comps

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
    pi_comp     = temp_comp + 1
    trac_comp   = pi_comp + 1

    ntrac = 1

    ! The "4" here refers to rho, rhoh, temp, and pi (the perturbation pressure)
    nscal = nspec + ntrac + 4

    ! press_comp here is used in the elliptic solves.  This slot is here
    ! for the bc tower
    press_comp  = dm_in + nscal + 1

    foextrap_comp = press_comp + 1
    hoextrap_comp = foextrap_comp + 1

  end subroutine init_variables

  subroutine init_plot_variables(plotidx)

    use network, only: nspec
    use probin_module, only: plot_spec, plot_trac, plot_base, use_thermal_diffusion, &
         plot_omegadot, plot_Hnuc, plot_Hext, plot_eta, plot_ad_excess, &
         use_tfromp, plot_h_with_use_tfromp, plot_gpi, plot_cs, dm_in, use_particles, &
         plot_processors, plot_pidivu
    use geometry, only: spherical

    type(plot_t), intent(inout) :: plotidx

    np = 0

    plotidx%icomp_vel      = get_next_plot_index(dm_in, np)
    plotidx%icomp_rho      = get_next_plot_index(1, np)
    if (.not. use_tfromp .or. (use_tfromp .and. plot_h_with_use_tfromp)) then
       plotidx%icomp_rhoh     = get_next_plot_index(1, np)
       plotidx%icomp_h        = get_next_plot_index(1, np)
    end if

    if (plot_spec) plotidx%icomp_spec = get_next_plot_index(nspec, np)
    if (plot_trac) plotidx%icomp_trac = get_next_plot_index(ntrac, np)

    if (plot_base) then
       plotidx%icomp_w0    = get_next_plot_index(dm_in, np)
       plotidx%icomp_divw0 = get_next_plot_index(1, np)
       plotidx%icomp_rho0  = get_next_plot_index(1, np)
       plotidx%icomp_rhoh0 = get_next_plot_index(1, np)
       plotidx%icomp_h0    = get_next_plot_index(1, np)
       plotidx%icomp_p0    = get_next_plot_index(1, np)
    end if

    if (spherical .eq. 1) then
       plotidx%icomp_velr = get_next_plot_index(1, np)
       plotidx%icomp_velc = get_next_plot_index(1, np)
    end if

    plotidx%icomp_magvel      = get_next_plot_index(1, np)
    plotidx%icomp_mom         = get_next_plot_index(1, np)
    plotidx%icomp_vort        = get_next_plot_index(1, np)
    plotidx%icomp_src         = get_next_plot_index(1, np)
    plotidx%icomp_rhopert     = get_next_plot_index(1, np)

    if (.not. use_tfromp .or. (use_tfromp .and. plot_h_with_use_tfromp)) then
       plotidx%icomp_rhohpert    = get_next_plot_index(1, np)
    endif

    icomp_tfromp      = get_next_plot_index(1, np)
    if (.not. use_tfromp .or. (use_tfromp .and. plot_h_with_use_tfromp)) then
       plotidx%icomp_tfromH      = get_next_plot_index(1, np)
       plotidx%icomp_dT          = get_next_plot_index(1, np)
       plotidx%icomp_dp          = get_next_plot_index(1, np)
    endif

    plotidx%icomp_tpert       = get_next_plot_index(1, np)
    plotidx%icomp_machno      = get_next_plot_index(1, np)
    if (plot_cs) then
       plotidx%icomp_cs          = get_next_plot_index(1, np)
    end if
    plotidx%icomp_dg          = get_next_plot_index(1, np)
    plotidx%icomp_entropy     = get_next_plot_index(1, np)
    plotidx%icomp_entropypert = get_next_plot_index(1, np)
    plotidx%icomp_sponge      = get_next_plot_index(1, np)

    plotidx%icomp_pi          = get_next_plot_index(1, np)
    if (plot_gpi) then
       plotidx%icomp_gpi         = get_next_plot_index(dm_in, np)
    endif

    if (plot_base) then
       plotidx%icomp_pioverp0    = get_next_plot_index(1, np)
       plotidx%icomp_p0pluspi    = get_next_plot_index(1, np)
    end if

    if (plot_omegadot) then
       plotidx%icomp_omegadot = get_next_plot_index(nspec, np)
    end if

    if (plot_Hnuc) then
       plotidx%icomp_enuc     = get_next_plot_index(1, np)
    end if

    if (plot_Hext) then
      plotidx%icomp_Hext     = get_next_plot_index(1, np)
    end if

    if (plot_eta) then
      plotidx%icomp_eta     = get_next_plot_index(1, np)
    end if

    if (use_thermal_diffusion) then
       plotidx%icomp_thermal = get_next_plot_index(1, np)
       plotidx%icomp_conductivity = get_next_plot_index(1, np)
    endif

    if (plot_ad_excess) then
       plotidx%icomp_ad_excess = get_next_plot_index(1, np)
    endif

    if (use_particles) then
       plotidx%icomp_part = get_next_plot_index(1, np)
    endif

    if (plot_processors) then
       plotidx%icomp_proc = get_next_plot_index(1, np)
    endif

    if (plot_pidivu) then
       plotidx%icomp_pidivu = get_next_plot_index(1, np)
    endif

    plotidx%n_plot_comps = np

  end subroutine init_plot_variables

end module variables
