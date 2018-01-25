module urca_composition_module

  use bl_types
  use network, only: nspec

  implicit none

  private xn_in, xn_out
  public set_urca_composition, c12_in, c12_out, o16_in, o16_out, &
         ne23_in, ne23_out, na23_in, na23_out, urca_23_dens, &
         urca_shell_type, shell_atan_kappa

  real (kind=dp_t) :: c12_in  = 0.0d0, c12_out  = 0.0d0, &
                      o16_in  = 0.0d0, o16_out  = 0.0d0, &
                      ne23_in = 0.0d0, ne23_out = 0.0d0, &
                      na23_in = 0.0d0, na23_out = 0.0d0, &
                      urca_23_dens = 0.0d0, shell_atan_kappa = 0.0d0, &
                      na_ne_23 = 0.0d0

  real (kind=dp_t) :: xn_in(nspec), xn_out(nspec)

  character (len=128) :: urca_shell_type = ""

contains

  subroutine init_urca_composition()
    use bl_error_module
    use network, only: network_species_index

    implicit none

    integer :: ihe4, ic12, io16, ine20, ine23, ina23, img23

    ! get the species indices
    ihe4  = network_species_index("helium-4")
    ic12  = network_species_index("carbon-12")
    io16  = network_species_index("oxygen-16")
    ine20 = network_species_index("neon-20")
    ine23 = network_species_index("neon-23")
    ina23 = network_species_index("sodium-23")
    img23 = network_species_index("magnesium-23")

    ! check the required species exist
    if (ihe4 < 0 .or. &
        ic12 < 0 .or. io16 < 0 .or. ine20 < 0 .or. &
        ine23 < 0 .or. ina23 < 0 .or. img23 < 0) then
       call bl_error("ERROR: species not defined")
    endif

    ! check the species mass fractions
    if (c12_in < 0.0_dp_t .or. c12_in > 1.0_dp_t) then
       call bl_error("ERROR: c12_in must be between 0 and 1")
    endif
    if (c12_out < 0.0_dp_t .or. c12_out > 1.0_dp_t) then
       call bl_error("ERROR: c12_out must be between 0 and 1")
    endif
    if (o16_in < 0.0_dp_t .or. o16_in > 1.0_dp_t) then
       call bl_error("ERROR: o16_in must be between 0 and 1")
    endif
    if (o16_out < 0.0_dp_t .or. o16_out > 1.0_dp_t) then
       call bl_error("ERROR: o16_out must be between 0 and 1")
    endif
    if (ne23_in < 0.0_dp_t .or. ne23_in > 1.0_dp_t) then
       call bl_error("ERROR: ne23_in must be between 0 and 1")
    endif
    if (ne23_out < 0.0_dp_t .or. ne23_out > 1.0_dp_t) then
       call bl_error("ERROR: ne23_out must be between 0 and 1")
    endif
    if (na23_in < 0.0_dp_t .or. na23_in > 1.0_dp_t) then
       call bl_error("ERROR: na23_in must be between 0 and 1")
    endif
    if (na23_out < 0.0_dp_t .or. na23_out > 1.0_dp_t) then
       call bl_error("ERROR: na23_out must be between 0 and 1")
    endif

    ! Set the composition within and outside the urca shell
    ! If a non-jump profile is used at the shell, these are asymptotic
    xn_in(:)    = 0.0d0
    xn_in(ic12) = c12_in
    xn_in(io16) = o16_in
    xn_in(ine23) = ne23_in
    xn_in(ina23) = na23_in

    xn_out(:)    = 0.0d0
    xn_out(ic12) = c12_out
    xn_out(io16) = o16_out
    xn_out(ine23) = ne23_out
    xn_out(ina23) = na23_out
  end subroutine init_urca_composition


  subroutine set_urca_composition(eos_state, xn)

    ! Construct composition profiles given a choice of
    ! the type of profile denoted by urca_shell_type.
    ! "jump": impose sharp species discontinuity
    ! "atan": use arctan to smooth composition profile
    use bl_error_module
    use bl_types, only: dp_t
    use network
    use eos_type_module, only: eos_t

    type (eos_t), intent(in) :: eos_state
    real (kind=dp_t), intent(out), DIMENSION(nspec) :: xn

    if (urca_shell_type .eq. "jump") then
       call composition_jump(eos_state, urca_23_dens, xn, xn_in, xn_out)
    else if (urca_shell_type .eq. "atan") then
       call composition_atan(eos_state, urca_23_dens, xn, xn_in, xn_out, shell_atan_kappa)
    else if (urca_shell_type .eq. "equilibrium") then
       call composition_equilibrium(eos_state, xn)
    else
       call bl_error("ERROR: invalid urca_shell_type")
    end if

    ! Renormalize species
    call renormalize_species(xn)

  end subroutine set_urca_composition


  subroutine renormalize_species(xn)

    ! Renormalize the mass fractions so they sum to 1
    use bl_types, only: dp_t

    real (kind=dp_t) :: xn(nspec)
    real (kind=dp_t) :: sumx

    sumx = sum(xn)
    xn(:) = xn(:)/sumx

  end subroutine renormalize_species


  subroutine composition_jump(eos_state, shell_rho, xn, xn_left, xn_right)
    use bl_types
    use bl_constants_module, only: HALF
    use network
    use eos_type_module, only: eos_t

    type (eos_t), intent(in) :: eos_state
    real (kind=dp_t), intent( in) :: shell_rho
    real (kind=dp_t), intent(out), dimension(nspec) :: xn
    real (kind=dp_t), intent( in), dimension(nspec) :: xn_left, xn_right

    if (eos_state % rho < shell_rho) then
       xn = xn_right
    else if (eos_state % rho > shell_rho) then
       xn = xn_left
    else
       xn = HALF*(xn_right + xn_left)
    end if
  end subroutine composition_jump


  subroutine composition_atan(eos_state, shell_rho, xn, xn_left, xn_right, &
                              shell_atan_kappa)
    use bl_types
    use bl_constants_module, only: TWO, HALF, M_PI
    use network
    use eos_type_module, only: eos_t

    type (eos_t), intent(in) :: eos_state
    real (kind=dp_t), intent(in) :: shell_rho, shell_atan_kappa
    real (kind=dp_t), intent(out), dimension(nspec) :: xn
    real (kind=dp_t), intent(in),  dimension(nspec) :: xn_left, xn_right
    real (kind=dp_t), dimension(nspec) :: A, B

    B = HALF * (xn_left + xn_right)
    A = xn_right - B

    xn = (TWO/M_PI) * A * atan(shell_atan_kappa * (shell_rho - eos_state % rho)) + B
  end subroutine composition_atan


  subroutine composition_equilibrium(eos_state, xn)
    use bl_types
    use network
    use eos_type_module, only: eos_t
    use burn_type_module, only: burn_t, eos_to_burn

    type (eos_t), intent(in) :: eos_state
    real (kind=dp_t), intent(out), dimension(nspec) :: xn
    type (burn_t) :: burn_state

    call eos_to_burn(eos_state, burn_state)

    ! Keep the mass fraction sum X(ne23) + X(na23) = na_ne_23
    ! Find the A=23 mass fractions such that the A=23 Urca rate equilibrium is maintained
    ! STUB

  end subroutine composition_equilibrium
    

end module urca_composition_module
