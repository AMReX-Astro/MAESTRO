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
!
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

  integer, parameter :: nspec = 3

  character (len=16), save :: spec_names(nspec) 
  real(kind=dp_t), save :: aion(nspec), zion(nspec), ebin(nspec)

  logical, save :: network_initialized = .false.

contains
  
  subroutine network_init()

    spec_names(1) = "carbon-12"
    spec_names(2) = "oxygen-16"
    spec_names(3) = "magnesium-24"
    
    aion(1) = 12.0_dp_t
    aion(2) = 16.0_dp_t
    aion(3) = 24.0_dp_t
    
    zion(1) = 6.0_dp_t
    zion(2) = 8.0_dp_t
    zion(3) = 12.0_dp_t

    ! our convention is that the binding energies are negative.  We convert
    ! from the MeV values that are traditionally written in astrophysics 
    ! papers by multiplying by 1.e6 eV/MeV * 1.60217646e-12 erg/eV.  The
    ! MeV values are per nucleus, so we divide by aion to make it per
    ! nucleon and we multiple by Avogardo's # (6.0221415e23) to get the 
    ! value in erg/g
    ebin(1) = -7.4103097e18_dp_t     !  92.16294 MeV
    ebin(2) = -7.6959672e18_dp_t     ! 127.62093 MeV
    ebin(3) = -7.9704080e18_dp_t     ! 198.2579  MeV

    network_initialized = .true.

  end subroutine network_init

  
  function network_species_index(name)

    character(len=*) :: name
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

end module network
