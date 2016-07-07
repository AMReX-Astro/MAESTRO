module rpar_indices

  use network, only: nspec, nrat
  implicit none

  !The number of equations the burner wants to evolve.  For rprox, we
  !evolve the species and temperature.
  integer, parameter :: neqs  = nspec + 1

  !The maximum order the burner will allow for variable-order integrators
  integer, parameter :: burn_max_order = 6

  !If the burner can be operate on a vector, 
  !this is the number of points to use in the vector.
  integer, parameter :: burn_npts = 1

  !This tells the burner if we want to reset data before burning.  In this case,
  !we always want to.  In some use cases, you may want to call the burner
  !several times, using data from previous burns.
  logical, parameter :: reset = .true.

  !This tells the burner if we want it to reuse the Jacobian.
  logical, parameter :: reuse = .false.

  !The total number of components in the user's real parameter array.  This is
  !for storing any real data the network's RHS or Jacobian wants.
  integer, parameter :: n_rpar_comps = 11 + 3*nspec + 2*nrat
 
  integer, parameter :: irp_dens          = 1
  integer, parameter :: irp_cp            = irp_dens   + 1
  integer, parameter :: irp_dhdX          = irp_cp     + 1
  integer, parameter :: irp_T9_eos        = irp_dhdX   + nspec
  integer, parameter :: irp_dTcrit        = irp_T9_eos + nspec

  ! rate-related 
  integer, parameter :: irp_rates         = irp_dTcrit + nspec
  integer, parameter :: irp_drtdt         = irp_rates  + nrat
  integer, parameter :: irp_dlambCNOdh1   = irp_drtdt  + nrat
  integer, parameter :: irp_drs1dhe4      = irp_dlambCNOdh1  + 1
  integer, parameter :: irp_drr1dh1       = irp_drs1dhe4     + 1
  integer, parameter :: irp_dlambda1dhe4  = irp_drr1dh1      + 1
  integer, parameter :: irp_dlambda2dhe4  = irp_dlambda1dhe4 + 1
  integer, parameter :: irp_delta1        = irp_dlambda2dhe4 + 1
  integer, parameter :: irp_delta2        = irp_delta1 + 1
  integer, parameter :: irp_r56eff        = irp_delta2 + 1
  integer, parameter :: irp_dr56effdt     = irp_r56eff + 1

contains

  !function get_next_rpar_index(num) result (next)

  !  ! return the next starting index for a plotfile quantity,
  !  ! and increment the counter of plotfile quantities by num
  !  integer :: num, next

  !  next = n_rpar_comps + 1
  !  n_rpar_comps = n_rpar_comps + num

  !  return
  !end function get_next_rpar_index


  subroutine init_rpar_indices(nrat, nspec)

    integer, intent(in) :: nrat, nspec

    !irp_dens  = get_next_rpar_index(1)
    !irp_cp    = get_next_rpar_index(1)
    !irp_dhdX  = get_next_rpar_index(nspec)
    !irp_T9_eos = get_next_rpar_index(nspec)
    !irp_dTcrit = get_next_rpar_index(nspec)

    !!============================================================
    !! only rate-related things below this point
    !irp_rates = get_next_rpar_index(nrat)
    !irp_drtdt = get_next_rpar_index(nrat)

    !irp_dlambCNOdh1 = get_next_rpar_index(1)
    !irp_drs1dhe4 = get_next_rpar_index(1)
    !irp_drr1dh1 = get_next_rpar_index(1)
    !irp_dlambda1dhe4 = get_next_rpar_index(1)
    !irp_dlambda2dhe4 = get_next_rpar_index(1)

    !irp_delta1 = get_next_rpar_index(1)
    !irp_delta2 = get_next_rpar_index(1)

    !irp_r56eff = get_next_rpar_index(1)
    !irp_dr56effdt = get_next_rpar_index(1)

  end subroutine init_rpar_indices

end module rpar_indices
