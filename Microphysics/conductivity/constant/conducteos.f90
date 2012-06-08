module conductivity_module

  use bl_types
  implicit none

  real (kind=dp_t), save, private :: conductivity_constant

  interface conducteos
     module procedure conducteos_old
     module procedure conducteos_new
  end interface conducteos

contains

  subroutine conductivity_init(cond_const)

    real (kind=dp_t), optional :: cond_const

    if (present(cond_const)) then
       conductivity_constant = cond_const
    else
       conductivity_constant = 1.0_dp_t
    endif

  end subroutine conductivity_init


  subroutine conducteos_new(input, eos_state, do_diag, conductivity)

    use eos_type_module
    use network, only: nspec

    integer         , intent(in   ) :: input
    type (eos_t)    , intent(inout) :: eos_state
    logical         , intent(in   ) :: do_diag
    real (kind=dp_t), intent(inout) :: conductivity

    call conducteos_old(input, eos_state%rho, eos_state%T, &
                        nspec, &
                        eos_state%xn, &
                        eos_state%p, eos_state%h, eos_state%e, &
                        eos_state%cv, eos_state%cp, eos_state%xne, &
                        eos_state%eta, eos_state%pele, &
                        eos_state%dpdT, eos_state%dpdr, &
                        eos_state%dedT, eos_state%dedr, &
                        eos_state%dpdX, eos_state%dhdX, &
                        eos_state%gam1, eos_state%cs, eos_state%s, &
                        eos_state%dsdT, eos_state%dsdr, &
                        do_diag, &
                        conductivity)
    
  end subroutine conducteos_new

  subroutine conducteos_old(input, dens, temp, &
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

  end subroutine conducteos_old

end module conductivity_module

