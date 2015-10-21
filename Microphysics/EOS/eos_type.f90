module eos_type_module
  
  use bl_types
  use network
  use eos_data_module
  use mempool_module
  use bl_constants_module
  
  implicit none

  ! A generic structure holding thermodynamic quantities and their derivatives,
  ! plus some other quantities of interest.

  ! rho      -- mass density (g/cm**3)
  ! T        -- temperature (K)
  ! xn       -- the mass fractions of the individual isotopes
  ! p        -- the pressure (dyn/cm**2)
  ! h        -- the enthalpy (erg/g)
  ! e        -- the internal energy (erg/g)
  ! s        -- the entropy (erg/g/K)
  ! c_v      -- specific heat at constant volume
  ! c_p      -- specific heat at constant pressure
  ! ne       -- number density of electrons + positrons
  ! np       -- number density of positrons only
  ! eta      -- degeneracy parameter
  ! pele     -- electron pressure + positron pressure
  ! ppos     -- position pressure only
  ! mu       -- mean molecular weight
  ! mu_e     -- mean number of nucleons per electron
  ! y_e      -- electron fraction == 1 / mu_e
  ! dPdT     -- d pressure/ d temperature
  ! dPdr     -- d pressure/ d density
  ! dedT     -- d energy/ d temperature
  ! dedr     -- d energy/ d density
  ! dsdT     -- d entropy/ d temperature
  ! dsdr     -- d entropy/ d density
  ! dhdT     -- d enthalpy/ d temperature
  ! dhdr     -- d enthalpy/ d density
  ! dPdX     -- d pressure / d xmass
  ! dhdX     -- d enthalpy / d xmass at constant pressure
  ! gam1     -- first adiabatic index (d log P/ d log rho) |_s
  ! cs       -- sound speed
  ! abar     -- average atomic number ( sum_k {X_k} ) / ( sum_k {X_k/A_k} )
  ! zbar     -- average proton number ( sum_k {Z_k X_k/ A_k} ) / ( sum_k {X_k/A_k} )
  ! dpdA     -- d pressure/ d abar
  ! dpdZ     -- d pressure/ d zbar
  ! dedA     -- d energy/ d abar
  ! dedZ     -- d energy/ d zbar

  ! Initialize the main quantities to an unphysical number
  ! so that we know if the user forgot to initialize them
  ! when calling the EOS in a particular mode.

  double precision, parameter :: init_num  = -1.0d200
  double precision, parameter :: init_test = -1.0d199

  type :: eos_t

    double precision :: rho         = init_num
    double precision :: T           = init_num
    double precision :: p           = init_num
    double precision :: e           = init_num
    double precision :: h           = init_num
    double precision :: s           = init_num
    double precision :: dpdT
    double precision :: dpdr
    double precision :: dedT
    double precision :: dedr
    double precision :: dhdT
    double precision :: dhdr
    double precision :: dsdT
    double precision :: dsdr
    double precision :: dpde
    double precision :: dpdr_e

    double precision :: xn(nspec)   = init_num
    double precision :: aux(naux)   = init_num
    double precision :: cv
    double precision :: cp
    double precision :: xne
    double precision :: xnp
    double precision :: eta
    double precision :: pele
    double precision :: ppos
    double precision :: mu
    double precision :: mu_e
    double precision :: y_e
    double precision :: dedX(nspec)
    double precision :: dpdX(nspec)
    double precision :: dhdX(nspec)
    double precision :: gam1
    double precision :: cs

    double precision :: abar
    double precision :: zbar
    double precision :: dpdA

    double precision :: dpdZ
    double precision :: dedA
    double precision :: dedZ

    logical :: reset
    
  end type eos_t


contains

  
  ! Given a set of mass fractions, calculate quantities that depend
  ! on the composition like abar and zbar.

  subroutine composition(state)

    use bl_constants_module
    use network

    implicit none

    type (eos_t), intent(inout) :: state

    ! Calculate abar, the mean nucleon number,
    ! zbar, the mean proton number,
    ! mu, the mean molecular weight,
    ! mu_e, the mean number of nucleons per electron, and
    ! y_e, the electron fraction.

    state % mu_e = ONE / (sum(state % xn(:) * zion(:) / aion(:)))       
    state % y_e = ONE / state % mu_e

    state % abar = ONE / (sum(state % xn(:) / aion(:)))
    state % zbar = state % abar / state % mu_e

  end subroutine composition

  
  
  ! Compute thermodynamic derivatives with respect to xn(:)

  subroutine composition_derivatives(state)

    use bl_constants_module
    use network

    implicit none

    type (eos_t), intent(inout) :: state

    
    state % dpdX(:) = state % dpdA * (state % abar/aion(:)) &
         * (aion(:) - state % abar)             &
         + state % dpdZ * (state % abar/aion(:)) &
         * (zion(:) - state % zbar)

    state % dEdX(:) = state % dedA * (state % abar/aion(:)) &
         * (aion(:) - state % abar)             &
         + state % dedZ * (state % abar/aion(:)) &
         * (zion(:) - state % zbar)

    state % dhdX(:) = state % dedX(:) &
         + (state % p / state % rho**2 - state % dedr) &
         *  state % dPdX(:) / state % dPdr

  end subroutine composition_derivatives



  ! Normalize the mass fractions: they must be individually positive and less than one,
  ! and they must all sum to unity.

  subroutine normalize_abundances(state)

    use bl_constants_module
    use network

    implicit none

    type (eos_t), intent(inout) :: state

    integer :: i    

    state % xn(:) = max(smallx, min(ONE, state % xn(:)))
    
    state % xn(:) = state % xn(:) / sum(state % xn(:))

  end subroutine normalize_abundances



  subroutine eos_copy(state_in, state_out)

    type(eos_t) :: state_in, state_out

    state_out % rho = state_in % rho
    state_out % T   = state_in % T
    state_out % p   = state_in % p
    state_out % e   = state_in % e
    state_out % h   = state_in % h
    state_out % s   = state_in % s
    state_out % dpdT = state_in % dpdT
    state_out % dpdr = state_in % dpdr
    state_out % dedT = state_in % dedT
    state_out % dedr = state_in % dedr
    state_out % dhdT = state_in % dhdT
    state_out % dhdr = state_in % dhdr
    state_out % dsdT = state_in % dsdT
    state_out % dsdr = state_in % dsdr
    state_out % dpde = state_in % dpde
    state_out % dpdr_e = state_in % dpdr_e
    state_out % xn = state_in % xn
    state_out % aux = state_in % aux
    state_out % cv = state_in % cv
    state_out % cp = state_in % cp
    state_out % xne = state_in % xne
    state_out % xnp = state_in % xnp
    state_out % eta = state_in % eta
    state_out % pele = state_in % pele
    state_out % ppos = state_in % ppos
    state_out % mu = state_in % mu
    state_out % mu_e = state_in % mu
    state_out % y_e = state_in % y_e
    state_out % dedX = state_in % dedX
    state_out % dpdX = state_in % dpdX
    state_out % dhdX = state_in % dhdX
    state_out % gam1 = state_in % gam1
    state_out % cs = state_in % cs
    state_out % abar = state_in % abar
    state_out % zbar = state_in % zbar
    state_out % dpdA = state_in % dpdA
    state_out % dpdZ = state_in % dpdZ
    state_out % dedA = state_in % dedA
    state_out % dedZ = state_in % dedZ

    state_out % reset = state_in % reset


  end subroutine eos_copy



  subroutine eos_type_error(err, input, pt_index)

    use bl_error_module
    
    implicit none

    integer,           intent(in) :: err
    integer,           intent(in) :: input
    integer, optional, intent(in) :: pt_index(3)

    integer :: dim_ptindex

    character (len=64) :: err_string, zone_string, eos_input_str

    write(eos_input_str, '(A13, I1)') ' EOS input = ', input

    if (err .eq. ierr_input) then

      err_string = 'EOS: invalid input.'

    elseif (err .eq. ierr_init) then

      err_string = 'EOS: input variables were not initialized.'

    elseif (err .eq. ierr_init_xn) then

      err_string = 'EOS: species abundances were not initialized.'

    else

      err_string = 'EOS: invalid input to error handler.'

    endif

    err_string = err_string // eos_input_str

    ! this format statement is for writing into zone_string -- make sure that
    ! the len of z_err can accomodate this format specifier
1001 format(1x,"zone index info: i = ", i5)
1002 format(1x,"zone index info: i = ", i5, '  j = ', i5)
1003 format(1x,"zone index info: i = ", i5, '  j = ', i5, '  k = ', i5)

    if (present(pt_index)) then

       dim_ptindex = 3

       if (pt_index(3) .eq. -99) then
          dim_ptindex = 2
          if (pt_index(2) .eq. -99) then
             dim_ptindex = 1
             if (pt_index(1) .eq. -99) then
                dim_ptindex = 0
             endif
          endif
       endif

       if (dim_ptindex .eq. 1) then 
          write (zone_string,1001) pt_index(1)
       else if (dim_ptindex .eq. 2) then 
          write (zone_string,1002) pt_index(1), pt_index(2)
       else if (dim_ptindex .eq. 3) then 
          write (zone_string,1003) pt_index(1), pt_index(2), pt_index(3)
       end if

    else

      zone_string = ''

    endif

    err_string = err_string // zone_string
    
    call bl_error(err_string)

  end subroutine eos_type_error

  
end module eos_type_module
