! This module provides information about the species in the network:
!
! nspec       -- the number of species
!
! spec_names  -- the names of the species
! rates       -- rate jacobian 
!
!
! This module contains 2 routines:
!
!  network_init()         -- initialize the isotope properties
!
!  network_species_index  -- return the index of the species given its name
!
!

module network

  use bl_types

  implicit none

  integer, parameter :: nspec = 2
  integer, parameter :: naux  = 0

  character, save :: spec_names(nspec)

  real(kind=dp_t), save :: rates(nspec,nspec)

  logical, save :: network_initialized = .false.

contains

  subroutine network_init()

    use bl_constants_module
    
    ! set the names of the species
    spec_names(1) = "A"
    spec_names(2) = "B"

    rates(:,:) = ZERO

    network_initialized = .true.

    return

  end subroutine network_init



  function network_species_index(name)

    character (len=*) :: name
    
    integer :: network_species_index, n

    do n = 1, nspec

       if (name == spec_names(n)) then
          network_species_index = n
          exit
       endif

    enddo

    return

  end function network_species_index


end module network
