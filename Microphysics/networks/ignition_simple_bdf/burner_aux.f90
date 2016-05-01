!module network_indices
!
!  ! this module is for use only within this network -- these
!  ! quantities should not be accessed in general MAESTRO routines.
!  ! Instead the species indices should be queried via
!  ! network_species_index()
!
!  implicit none
!
!  integer, parameter :: ic12_ = 1
!  integer, parameter :: io16_ = 2
!  integer, parameter :: img24_ = 3
!
!end module network_indices

module burner_data
  use network, only: nspec, nspec_advance
  implicit none

  !The number of equations the burner wants to evolve.  For ignition_simple, we
  !evolve the advancing species as well as temperature.
  integer, parameter :: neqs  = nspec_advance + 1

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
  integer, parameter :: n_rpar_comps  = 8 + nspec
  !!$acc declare copyin(n_rpar_comps)

  !Indices for the data in rpar.
  integer, parameter :: irp_dens      = 1
  integer, parameter :: irp_cp        = irp_dens + 1
  integer, parameter :: irp_dhdx      = irp_cp + 1
  integer, parameter :: irp_o16       = irp_dhdx + nspec
  integer, parameter :: irp_rate      = irp_o16 + 1
  integer, parameter :: irp_dratedt   = irp_rate + 1
  integer, parameter :: irp_sc1212    = irp_dratedt + 1
  integer, parameter :: irp_dsc1212dt = irp_sc1212 + 1
  integer, parameter :: irp_xc12tmp   = irp_dsc1212dt + 1
  !!$acc declare copyin(irp_dens, irp_cp, irp_dhdx, irp_o16, irp_rate, &
  !!$acc    irp_dratedt, irp_sc1212, irp_dsc1212dt, irp_xc12tmp)

contains

  !function get_next_rpar_index(num) result (next)

  !  ! return the next starting index for a plotfile quantity,
  !  ! and increment the counter of plotfile quantities by num
  !  integer :: num, next

  !  next = n_rpar_comps + 1
  !  n_rpar_comps = n_rpar_comps + num

  !  return
  !end function get_next_rpar_index


  subroutine init_rpar_indices(nspec)

    integer, intent(in) :: nspec

    !irp_dens  = get_next_rpar_index(1)
    !!$acc update device(irp_dens)
    !irp_cp    = get_next_rpar_index(1)
    !!$acc update device(irp_cp)
    !irp_dhdX  = get_next_rpar_index(nspec)
    !!$acc update device(irp_dhdX)
    !irp_o16   = get_next_rpar_index(1)
    !!$acc update device(irp_o16)

    !irp_rate      = get_next_rpar_index(1)
    !!$acc update device(irp_rate)
    !irp_dratedt   = get_next_rpar_index(1)
    !!$acc update device(irp_dratedt)
    !irp_sc1212    = get_next_rpar_index(1)
    !!$acc update device(irp_sc1212)
    !irp_dsc1212dt = get_next_rpar_index(1)
    !!$acc update device(irp_dsc1212dt)
    !irp_xc12tmp   = get_next_rpar_index(1)
    !!$acc update device(irp_xc12tmp)

  end subroutine init_rpar_indices

end module burner_data
!end module rpar_indices
