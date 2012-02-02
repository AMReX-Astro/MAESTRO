module conductivity_module

  use bl_types
  implicit none

  real (kind=dp_t), save, private :: conductivity_constant

contains

  subroutine conductivity_init(cond_const)

    real (kind=dp_t), optional :: cond_const

    if (present(cond_const)) then
       conductivity_constant = cond_const
    else
       conductivity_constant = 1.0_dp_t
    endif

  end subroutine conductivity_init


  subroutine conducteos(input, dens, temp, &
                        nspecies, & 
                        xmass, &
                        pres, enthalpy, eint, &
                        c_v, c_p, ne, eta, pele, &
                        dPdT, dPdR, dEdT, dEdR, &
                        dPdX, dhdX, &
                        gam1, cs, entropy, &
                        dsdT, dsdR, &
                        do_eos_diag, &
                        conductivity)

    use eos_module

    implicit none

    ! arguments
    integer input,nspecies
    logical do_eos_diag
    double precision dens, temp
    double precision xmass(nspecies)
    double precision pres, enthalpy, eint
    double precision c_v, c_p
    double precision ne, eta, pele
    double precision dPdT, dPdR
    double precision dEdT, dEdR
    double precision gam1, entropy, cs
    double precision dPdX(nspecies), dhdX(nspecies)
    double precision dsdT, dsdR
    double precision conductivity

    ! first things first, call the eos
    call eos(input, dens, temp, &
             xmass, &
             pres, enthalpy, eint, &
             c_v, c_p, ne, eta, pele, &
             dPdT, dPdR, dEdT, dEdR, &
             dPdX, dhdX, &
             gam1, cs, entropy, &
             dsdT, dsdR, &
             do_eos_diag)

    ! fill the conductivity
    conductivity = conductivity_constant

  end subroutine conducteos

end module conductivity_module

