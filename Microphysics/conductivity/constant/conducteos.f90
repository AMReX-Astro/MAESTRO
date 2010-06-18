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
                        npoints, nspecies, & 
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
    integer input,npoints,nspecies
    logical do_eos_diag
    double precision dens(npoints), temp(npoints)
    double precision xmass(npoints,nspecies)
    double precision pres(npoints), enthalpy(npoints), eint(npoints)
    double precision c_v(npoints), c_p(npoints)
    double precision ne(npoints), eta(npoints), pele(npoints)
    double precision dPdT(npoints), dPdR(npoints)
    double precision dEdT(npoints), dEdR(npoints)
    double precision gam1(npoints), entropy(npoints), cs(npoints)
    double precision dPdX(npoints,nspecies), dhdX(npoints,nspecies)
    double precision dsdT(npoints), dsdR(npoints)
    double precision conductivity(npoints)
    double precision xmass_temp(nspecies)

    ! first things first, call the eos
    call eos(input, dens, temp, &
             npoints, &
             xmass, &
             pres, enthalpy, eint, &
             c_v, c_p, ne, eta, pele, &
             dPdT, dPdR, dEdT, dEdR, &
             dPdX, dhdX, &
             gam1, cs, entropy, &
             dsdT, dsdR, &
             do_eos_diag)

    ! fill the conductivity
    conductivity(1) = conductivity_constant

  end subroutine conducteos

end module conductivity_module

