! the network module provides the information about the species we are 
! advecting:
!
! nspec      -- the number of species in the network
! nrat       -- the number of reactions in this network
!
! aion       -- atomic number
! zion       -- proton number
! ebin       -- nuclear binding energy
!
! spec_names -- the name of the isotope
! short_spec_names -- abbreviated names
!
! reac_names -- the name of the reaction
!
!
! This module contains three routines:
!
!  network_init()         -- initialize the isotope properties
!
!  network_species_index  -- return the index of the species given its name
!
!  network_reaction_index -- return the index of the reaction given its name
!

module network

  use bl_types

  implicit none

  integer, parameter :: nspec = 3, nrat = 2

  character (len=16), save :: spec_names(nspec)
  character (len=5),  save :: short_spec_names(nspec)

  character (len=10), save :: reac_names(nrat)

  real(kind=dp_t), save :: aion(nspec), zion(nspec), ebin(nspec)

  logical, save :: network_initialized = .false.

contains

  subroutine network_init()

    ! set the names
    spec_names(1) = "helium-4"
    spec_names(2) = "carbon-12"
    spec_names(3) = "nickel-56"

    short_spec_names(1)  = "He4"
    short_spec_names(2)  = "C12"
    short_spec_names(3) = "Ni56"

    reac_names(1) = "forward"   ! 3 He4 --> C12
    reac_names(2) = "backward"  !   C12 --> 3 He4 

    ! set the species properties
    aion(1) =  4.0_dp_t
    aion(2) = 12.0_dp_t
    aion(3) = 56.0_dp_t

    zion(1) =  2.0_dp_t
    zion(2) =  6.0_dp_t
    zion(3) = 28.0_dp_t

    ! our convention is that binding energy is negative.  The following are 
    ! the binding energies per unit mass (erg / g) obtained by converting
    ! the energies in MeV to erg then multiplying by (N_A / aion) where 
    ! N_A = 6.0221415e23 is Avogadro's number
    ebin(1) = -6.8253797e18_dp_t    !  28.39603 MeV / nucleon
    ebin(2) = -7.4103097e18_dp_t    !  92.16294 MeV / nucleon
    ebin(3) = -8.3391412e18_dp_t    ! 484.00300 MeV / nucleon

    ! done initializing
    network_initialized = .true.

  end subroutine network_init

  
  function network_species_index(name)

    character (len=*) :: name
    integer :: network_species_index, n

    network_species_index = -1

    do n = 1, nspec
       if (trim(name) == trim(spec_names(n))) then
          network_species_index = n
          exit
       endif
    enddo

    return
  end function network_species_index


  function network_reaction_index(name)
    
    character(len=*) :: name
    integer :: network_reaction_index, n

    network_reaction_index = -1

    do n = 1, nrat
       if (trim(name) == trim(reac_names(n))) then
          network_reaction_index = n
          exit
       endif
    enddo

    return
  end function network_reaction_index


end module network
