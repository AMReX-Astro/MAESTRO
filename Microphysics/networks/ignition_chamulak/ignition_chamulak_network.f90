! the network module provides the information about the species we are
! advecting:
!
! nspec      -- the number of species
!
! aion       -- atomic number
! zion       -- proton number
!
! spec_names -- the name of the isotope
! short_spec_names -- the abbreviated name of the isotope
!
! This module contains two routines:
!
!  network_init()        -- initialize the isotope properties
!
!  network_species_index -- return the index of the species given its name
!

module actual_network

  use bl_types

  implicit none

  character (len=*), parameter :: network_name = "ignition_chamulak"

  ! nspec = number of species this network carries
  ! nspec_advance = the number of species that are explicitly integrated
  !                 in the ODE solve (the others are solved for
  !                 algebraically).
  integer, parameter :: nspec = 3
  integer, parameter :: nspec_advance = 1
  integer, parameter :: naux  = 0

  character (len=16), save :: spec_names(nspec)
  character (len= 5), save :: short_spec_names(nspec)
  character (len= 5), save :: short_aux_names(naux)

  real(kind=dp_t), save :: aion(nspec), zion(nspec)

  ! M12_chamulak is the effective number of C12 nuclei destroyed per
  ! reaction
  real(kind=dp_t), parameter :: M12_chamulak = 2.93_dp_t

!  logical, save :: network_initialized = .false.

contains

  subroutine actual_network_init()

    use network_indices
    use rpar_indices

    spec_names(ic12_)  = "carbon-12"
    spec_names(io16_)  = "oxygen-16"
    spec_names(iash_)  = "ash"

    short_spec_names(ic12_)  = "C12"
    short_spec_names(io16_)  = "O16"
    short_spec_names(iash_) = "ash"

    ! the ash from C12 burning according to Chamulak et al. is a mixture
    ! of C13, O16, Ne20, and Na23.   Ne20 + alpha results 60% of the time,
    ! while Na23 + p is the other 40%.   Fusing 6 C12 will result in
    ! 1.2 Na23, 1.2 O16 (from the alpha), 1.8 Ne20, and 1.8 C13.
    ! The ash state will have an A and Z corresponding to this mixture.
    aion(ic12_)  = 12.0_dp_t
    aion(io16_)  = 16.0_dp_t
    aion(iash_)  = 18.0_dp_t

    zion(ic12_)  = 6.0_dp_t
    zion(io16_)  = 8.0_dp_t
    zion(iash_)  = 8.8_dp_t

    ! rpar is VODE's way of passing information into the RHS and
    ! jacobian routines.  Here we initialize some indices to make
    ! sense of what is stored in the rpar() array.
    call init_rpar_indices(nspec)

!    network_initialized = .true.

  end subroutine actual_network_init

  function get_ebin_value(density) result (ebin)

    use fundamental_constants_module

    real (kind=dp_t) :: density, rho9, q_eff, ebin

    ! Chamulak et al. provide the q-value resulting from C12 burning,
    ! given as 3 different values (corresponding to 3 different densities).
    ! Here we do a simple quadratic fit to the 3 values provided (see
    ! Chamulak et al., p. 164, column 2).

    ! our convention is that the binding energies are negative.  We convert
    ! from the MeV values that are traditionally written in astrophysics
    ! papers by multiplying by 1.e6 eV/MeV * 1.60217646e-12 erg/eV.  The
    ! MeV values are per nucleus, so we divide by aion to make it per
    ! nucleon and we multiple by Avogardo's # (6.0221415e23) to get the
    ! value in erg/g
    rho9 = density/1.0e9_dp_t

    ! q_eff is effective heat evolved per reaction (given in MeV)
    q_eff = 0.06_dp_t*rho9**2 + 0.02_dp_t*rho9 + 8.83_dp_t

    ! convert from MeV to ergs / gram.  Here M12_chamulak is the
    ! number of C12 nuclei destroyed in a single reaction and 12.0 is
    ! the mass of a single C12 nuclei.  Also note that our convention
    ! is that binding energies are negative.
    ebin = -q_eff*MeV2eV*eV2erg*n_A/(M12_chamulak*12.0_dp_t)

    return

  end function get_ebin_value

  ! function network_species_index(name) result(r)

  !   character(len=*) :: name
  !   integer :: r, n

  !   r = -1

  !   do n = 1, nspec
  !      if (name == spec_names(n) .or. name == short_spec_names(n)) then
  !         r = n
  !         exit
  !      endif
  !   enddo
  !   return
  ! end function network_species_index

  subroutine network_finalize()

  end subroutine network_finalize

end module actual_network
