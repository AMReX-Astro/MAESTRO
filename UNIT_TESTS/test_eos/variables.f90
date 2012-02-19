!
! A module to provide integer indices into the various storage arrays
! for accessing the different variables by name.
!
module variables

  use bl_types

  implicit none

  integer, save :: rho_comp, temp_comp, spec_comp, &
                   h_comp, p_comp, e_comp, s_comp

  integer, save :: press_comp, foextrap_comp, hoextrap_comp

  ! the total number of plot components
  integer, save :: n_plot_comps = 0

  character (len=20), save, allocatable :: varnames(:)

  integer, save :: nscal

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
    use network, only: nspec, short_spec_names

    integer :: n

    rho_comp    = 1
    temp_comp   = 2
    spec_comp   = 3
    h_comp      = spec_comp + nspec
    p_comp      = h_comp + 1
    e_comp      = p_comp + 1
    s_comp      = e_comp + 1
    
    nscal = 6 + nspec

    allocate (varnames(nscal))

    varnames(rho_comp)  = "density"
    varnames(temp_comp) = "temperature"
    do n = 1, nspec
       varnames(spec_comp-1+n) = "X(" // trim(short_spec_names(n)) // ")"
    enddo
    varnames(h_comp)    = "enthalpy"
    varnames(p_comp)    = "pressure"
    varnames(e_comp)    = "internal energy"
    varnames(s_comp)    = "entropy"

    press_comp  = dm_in + nscal + 1

    foextrap_comp = press_comp + 1
    hoextrap_comp = foextrap_comp + 1

  end subroutine init_variables

  subroutine init_plot_variables()

    ! not used in this problem

  end subroutine init_plot_variables

end module variables
