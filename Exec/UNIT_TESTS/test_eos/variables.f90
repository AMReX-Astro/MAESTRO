!
! A module to provide integer indices into the various storage arrays
! for accessing the different variables by name.
!
module variables

  use bl_types

  implicit none

  integer, save :: rho_comp, rhoh_comp, temp_comp, spec_comp, &
                   h_comp, p_comp, e_comp, s_comp
  integer, save :: tfromrh_comp, rfromtp_comp, tfromrp_comp, &
                   tfromre_comp, tfromps_comp
  integer, save :: tfromrh_err_comp, rfromtp_err_comp, tfromrp_err_comp, &
                   tfromre_err_comp, tfromps_err_comp

  ! this is only needed to compile -- we don't have any ghost cells, so the
  ! value here doesn't matter
  integer, save :: foextrap_comp = -1  

  ! the total number of state quantities we will deal with
  integer, save :: nscal = 0

  character (len=20), save, allocatable :: varnames(:)

contains

  function get_next_scal_index(num) result (next)

    ! return the next starting index for a state quantity,
    ! and increment the counter of state quantities by num
    integer :: num, next

    next = nscal + 1
    nscal = nscal + num

    return
  end function get_next_scal_index

  subroutine init_variables()

    use probin_module, only: dm_in
    use network, only: nspec, short_spec_names

    integer :: n

    rho_comp    = get_next_scal_index(1)
    rhoh_comp   = get_next_scal_index(1)
    temp_comp   = get_next_scal_index(1)
    spec_comp   = get_next_scal_index(nspec)
    h_comp      = get_next_scal_index(1)
    p_comp      = get_next_scal_index(1)
    e_comp      = get_next_scal_index(1)
    s_comp      = get_next_scal_index(1)

    tfromrh_comp = get_next_scal_index(1)
    rfromtp_comp = get_next_scal_index(1)
    tfromrp_comp = get_next_scal_index(1)
    tfromre_comp = get_next_scal_index(1)
    tfromps_comp = get_next_scal_index(1)

    tfromrh_err_comp = get_next_scal_index(1)
    rfromtp_err_comp = get_next_scal_index(1)
    tfromrp_err_comp = get_next_scal_index(1)
    tfromre_err_comp = get_next_scal_index(1)
    tfromps_err_comp = get_next_scal_index(1)

    allocate (varnames(nscal))

    varnames(rho_comp)  = "density"
    varnames(rhoh_comp)  = "rho * h"    
    varnames(temp_comp) = "temperature"
    do n = 1, nspec
       varnames(spec_comp-1+n) = "X(" // trim(short_spec_names(n)) // ")"
    enddo
    varnames(h_comp)    = "enthalpy"
    varnames(p_comp)    = "pressure"
    varnames(e_comp)    = "internal energy"
    varnames(s_comp)    = "entropy"

    varnames(tfromrh_comp) = "T(rho,h)"
    varnames(rfromtp_comp) = "rho(T,p)"
    varnames(tfromrp_comp) = "T(rho,p)"
    varnames(tfromre_comp) = "T(rho,e)"
    varnames(tfromps_comp) = "T(p,s)"

    varnames(tfromrh_err_comp) = "[T(rho,h) - T]/T"
    varnames(rfromtp_err_comp) = "[rho(T,p) - rho]/rho"
    varnames(tfromrp_err_comp) = "[T(rho,p) - T]/T"
    varnames(tfromre_err_comp) = "[T(rho,e) - T]/T"
    varnames(tfromps_err_comp) = "[TT(p,s) - T]/T"

  end subroutine init_variables

  subroutine init_plot_variables()

    ! not used in this problem

  end subroutine init_plot_variables

end module variables
