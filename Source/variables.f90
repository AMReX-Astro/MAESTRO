!
! A module to provide integer indices into the various storage arrays
! for accessing the different variables by name.
!
module variables

  use bl_types
  use network, only: nspec

  implicit none

  integer, save :: rho_comp, rhoh_comp, spec_comp, temp_comp, pi_comp
  integer, save :: trac_comp, press_comp
  integer, save :: foextrap_comp, hoextrap_comp

  integer, save :: ntrac, nscal
  real(kind=dp_t), save :: rel_eps

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
     integer :: icomp_s_cc = -1
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
     integer :: icomp_brunt = -1
     integer :: icomp_hp = -1
     integer :: icomp_grav = -1
     
     integer :: n_plot_comps = 0

     integer :: n_species = 0
     integer :: spec(nspec)      ! index in the list of species in the network
     integer :: ipf_spec(nspec)  ! index in the plotfile for a given species

     character(len=20), allocatable :: names(:)

     integer :: plot_int
     real(dp_t) :: plot_dt
     character(len=256) :: base_name

   contains
     procedure :: next_index => get_next_plot_index

  end type plot_t

contains

  function get_next_plot_index(this, num) result (next)

    ! return the next starting index for a plotfile quantity, and
    ! increment the counter of plotfile quantities, n_plot_comps, by
    ! num

    class(plot_t), intent(inout) :: this
    integer, intent(in) :: num
    integer :: next

    next = this%n_plot_comps + 1
    this%n_plot_comps = this%n_plot_comps + num

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


  subroutine init_plot_variables(p, plot_int, plot_dt, base_name)

    use network, only: nspec
    use probin_module, only: plot_spec, plot_trac, plot_base, use_thermal_diffusion, &
         plot_omegadot, plot_Hnuc, plot_Hext, plot_eta, plot_ad_excess, plot_brunt_freq, &
         use_tfromp, plot_h_with_use_tfromp, plot_gpi, plot_cs, dm_in, use_particles, &
         plot_processors, plot_pidivu, plot_hp, plot_grav
    use geometry, only: spherical, polar

    type(plot_t), intent(inout) :: p
    integer, intent(in) :: plot_int
    real(dp_t), intent(in) :: plot_dt
    character (len=*), intent(in) :: base_name

    ! general plotfile stuff
    p%plot_int = plot_int
    p%plot_dt = plot_dt
    p%base_name = base_name

    ! variable information
    p%icomp_vel      = p%next_index(dm_in)
    p%icomp_rho      = p%next_index(1)
    if (.not. use_tfromp .or. (use_tfromp .and. plot_h_with_use_tfromp)) then
       p%icomp_rhoh     = p%next_index(1)
       p%icomp_h        = p%next_index(1)
    end if

    if (plot_spec) p%icomp_spec = p%next_index(nspec)
    if (plot_trac) p%icomp_trac = p%next_index(ntrac)

    if (plot_base) then
       p%icomp_w0    = p%next_index(dm_in)
       p%icomp_divw0 = p%next_index(1)
       p%icomp_rho0  = p%next_index(1)
       p%icomp_rhoh0 = p%next_index(1)
       p%icomp_h0    = p%next_index(1)
       p%icomp_p0    = p%next_index(1)
    end if

    if (spherical .eq. 1 .or. polar .eq. 1) then
       p%icomp_velr = p%next_index(1)
       p%icomp_velc = p%next_index(1)
    end if

    p%icomp_magvel      = p%next_index(1)
    p%icomp_mom         = p%next_index(1)
    p%icomp_vort        = p%next_index(1)
    p%icomp_s_cc        = p%next_index(1)
    p%icomp_rhopert     = p%next_index(1)

    if (.not. use_tfromp .or. (use_tfromp .and. plot_h_with_use_tfromp)) then
       p%icomp_rhohpert    = p%next_index(1)
    endif

    p%icomp_tfromp      = p%next_index(1)
    if (.not. use_tfromp .or. (use_tfromp .and. plot_h_with_use_tfromp)) then
       p%icomp_tfromH      = p%next_index(1)
       p%icomp_dT          = p%next_index(1)
       p%icomp_dp          = p%next_index(1)
    endif

    p%icomp_tpert       = p%next_index(1)
    p%icomp_machno      = p%next_index(1)
    if (plot_cs) then
       p%icomp_cs          = p%next_index(1)
    end if
    p%icomp_dg          = p%next_index(1)
    p%icomp_entropy     = p%next_index(1)
    p%icomp_entropypert = p%next_index(1)
    p%icomp_sponge      = p%next_index(1)

    p%icomp_pi          = p%next_index(1)
    if (plot_gpi) then
       p%icomp_gpi         = p%next_index(dm_in)
    endif

    if (plot_base) then
       p%icomp_pioverp0    = p%next_index(1)
       p%icomp_p0pluspi    = p%next_index(1)
    end if

    if (plot_omegadot) then
       p%icomp_omegadot = p%next_index(nspec)
    end if

    if (plot_Hnuc) then
       p%icomp_enuc     = p%next_index(1)
    end if

    if (plot_Hext) then
      p%icomp_Hext     = p%next_index(1)
    end if

    if (plot_eta) then
      p%icomp_eta     = p%next_index(1)
    end if

    if (use_thermal_diffusion) then
       p%icomp_thermal = p%next_index(1)
       p%icomp_conductivity = p%next_index(1)
    endif

    if (plot_ad_excess) then
       p%icomp_ad_excess = p%next_index(1)
    endif

    if (plot_brunt_freq) then
       p%icomp_brunt = p%next_index(1)
    endif

    if (plot_grav) then
       p%icomp_grav = p%next_index(1)
    endif
    
    if (plot_hp) then
       p%icomp_hp = p%next_index(1)
    endif
    
    if (use_particles) then
       p%icomp_part = p%next_index(1)
    endif

    if (plot_processors) then
       p%icomp_proc = p%next_index(1)
    endif

    if (plot_pidivu) then
       p%icomp_pidivu = p%next_index(1)
    endif

  end subroutine init_plot_variables


  subroutine init_miniplot_variables(p, plot_int, plot_dt, base_name)

    use bl_error_module
    use network, only: nspec, network_species_index
    use probin_module, only: mini_plot_var1, mini_plot_var2, mini_plot_var3, &
                             mini_plot_var4, mini_plot_var5, mini_plot_var6, &
                             mini_plot_var7, mini_plot_var8, mini_plot_var9, &
                             dm_in, use_tfromp

    type(plot_t), intent(inout) :: p
    integer, intent(in) :: plot_int
    real(dp_t), intent(in) :: plot_dt
    character (len=*), intent(in) :: base_name

    integer :: nvar, n, idx

    character (len=256) :: vars(9)

    ! general plotfile stuff
    p%plot_int = plot_int
    p%plot_dt = plot_dt
    p%base_name = base_name

    ! variable information
    
    nvar = 0

    ! hacky, but we don't want to deal with character array namelist
    ! variables
    if (.not. mini_plot_var1 == "") then
       nvar = nvar + 1
       vars(nvar) = trim(mini_plot_var1)
    endif

    if (.not. mini_plot_var2 == "") then
       nvar = nvar + 1
       vars(nvar) = trim(mini_plot_var2)
    endif

    if (.not. mini_plot_var3 == "") then
       nvar = nvar + 1
       vars(nvar) = trim(mini_plot_var3)
    endif

    if (.not. mini_plot_var4 == "") then
       nvar = nvar + 1
       vars(nvar) = trim(mini_plot_var4)
    endif

    if (.not. mini_plot_var5 == "") then
       nvar = nvar + 1
       vars(nvar) = trim(mini_plot_var5)
    endif

    if (.not. mini_plot_var6 == "") then
       nvar = nvar + 1
       vars(nvar) = trim(mini_plot_var6)
    endif

    if (.not. mini_plot_var7 == "") then
       nvar = nvar + 1
       vars(nvar) = trim(mini_plot_var7)
    endif

    if (.not. mini_plot_var8 == "") then
       nvar = nvar + 1
       vars(nvar) = trim(mini_plot_var8)
    endif

    if (.not. mini_plot_var9 == "") then
       nvar = nvar + 1
       vars(nvar) = trim(mini_plot_var9)
    endif

    do n = 1, nvar
       select case (vars(n))

       case ("density")
          p%icomp_rho = p%next_index(1)

       case ("species")
          p%icomp_spec = p%next_index(nspec)

       case ("radvel")
          ! note: we assume in make_plot_variables.f90 that if
          ! velr > 0 then velc > 0 too
          p%icomp_velr = p%next_index(nspec)
          p%icomp_velc = p%next_index(nspec)

       case ("velocity")
          p%icomp_vel = p%next_index(dm_in)

       case ("vorticity")
          p%icomp_vort = p%next_index(1)

       case ("temperature")
          if (use_tfromp) then
             p%icomp_tfromp = p%next_index(1)
          else
             p%icomp_tfromh = p%next_index(1)
          endif

       case ("enuc")
          p%icomp_enuc = p%next_index(1)

       case ("mach")
          p%icomp_machno = p%next_index(1)

       case default

          ! did we specify a particular species?
          idx = network_species_index(vars(n))
          if (idx > 0) then
             p%n_species = p%n_species + 1
             p%spec(p%n_species) = idx
             p%ipf_spec(p%n_species) = p%next_index(1)
          else
             print *, n, vars(n)
             call bl_error("ERROR, plot variable undefined", vars(n))
          endif
       end select

    enddo

  end subroutine init_miniplot_variables

end module variables
