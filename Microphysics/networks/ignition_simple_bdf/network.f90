! the network module provides the information about the species we are
! advecting: 
!
! nspec      -- the number of species
!
! aion       -- atomic number
! zion       -- proton number
! eion       -- nuclear binding energy (in erg/g)
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

module network

  use bl_types

  implicit none

  character (len=*), parameter :: network_name = "ignition_simple"

  ! nspec = number of species this network carries
  ! nspec_advance = the number of species that are explicitly integrated
  !                 in the ODE solve (the others are solved for 
  !                 algebraically).
  integer, parameter :: nspec = 3
  integer, parameter :: nspec_advance = 1
  integer, parameter :: naux  = 0

  ! These indices are intended only for internal network usage.  
  ! External program units should use network_species_index()
  integer, parameter :: ic12_  = 1
  integer, parameter :: io16_  = 2
  integer, parameter :: img24_ = 3

  character (len=16) :: spec_names(nspec) 
  character (len= 5) :: short_spec_names(nspec)
  character (len= 5) :: short_aux_names(naux)

  real(kind=dp_t) :: aion(nspec), zion(nspec), ebin(nspec)
  !$acc declare create(aion(:), zion(:), ebin(:))

  logical :: network_initialized = .false.

contains
  
  subroutine network_init()

    spec_names(ic12_)  = "carbon-12"
    spec_names(io16_)  = "oxygen-16"
    spec_names(img24_) = "magnesium-24"

    short_spec_names(ic12_)  = "C12"
    short_spec_names(io16_)  = "O16"
    short_spec_names(img24_) = "Mg24"

    aion(ic12_)  = 12.0_dp_t
    aion(io16_)  = 16.0_dp_t
    aion(img24_) = 24.0_dp_t
    
    zion(ic12_)  = 6.0_dp_t
    zion(io16_)  = 8.0_dp_t
    zion(img24_) = 12.0_dp_t

    ! our convention is that the binding energies are negative.  We convert
    ! from the MeV values that are traditionally written in astrophysics 
    ! papers by multiplying by 1.e6 eV/MeV * 1.60217646e-12 erg/eV.  The
    ! MeV values are per nucleus, so we divide by aion to make it per
    ! nucleon and we multiple by Avogardo's # (6.0221415e23) to get the 
    ! value in erg/g
    ebin(ic12_)  = -7.4103097e18_dp_t     !  92.16294 MeV
    ebin(io16_)  = -7.6959672e18_dp_t     ! 127.62093 MeV
    ebin(img24_) = -7.9704080e18_dp_t     ! 198.2579  MeV

    !$acc update device(aion(:), zion(:), ebin(:))

    ! rpar is VODE's way of passing information into the RHS and
    ! jacobian routines.  Here we initialize some indices to make
    ! sense of what is stored in the rpar() array.
    !call init_rpar_indices(nspec)

    network_initialized = .true.

  end subroutine network_init

  function network_species_index(name) result(r)

    character(len=*) :: name
    integer :: r, n

    r = -1

    do n = 1, nspec
       if (name == spec_names(n) .or. name == short_spec_names(n)) then
          r = n
          exit
       endif
    enddo
    return
  end function network_species_index


  subroutine network_finalize()

  end subroutine network_finalize

end module network
