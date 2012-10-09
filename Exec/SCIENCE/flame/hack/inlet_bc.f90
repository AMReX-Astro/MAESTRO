! inlet_bc_module serves as a container to hold the inflow boundary 
! condition information.
!
! These quantities are initialized through a call to set_inlet_bcs(),
! which should be done on initialization and restart.

module inlet_bc_module

  use bl_types
  use bl_constants_module
  use bl_space
  use network, only: nspec
  use multifab_module
  use ml_layout_module

  implicit none

  real(dp_t), save    :: INLET_VEL         ! normal velocity through boundary
  real(dp_t), save    :: INLET_RHO
  real(dp_t), save    :: INLET_RHOH
  real(dp_t), save    :: INLET_TEMP
  real(dp_t), save    :: INLET_RHOX(nspec)
  real(dp_t), save    :: INLET_TRA

  logical, save :: inlet_bc_initialized = .false.

  real(dp_t), save, private :: INLET_VEL_OLD = ZERO

contains

  subroutine set_initial_inlet_bcs()

    ! initialize the inflow boundary condition variables
    ! this is called at the top of each timestep, to allow for
    ! time-varying inlet BCs based on the current fluid state.

    use eos_module
    use probin_module, ONLY: dens_fuel, temp_fuel, xc12_fuel, vel_fuel
    use network, only: network_species_index

    integer, save :: ic12, io16
    logical, save :: firstCall_params = .true.

    if (firstCall_params) then

       ic12 = network_species_index("carbon-12")
       io16 = network_species_index("oxygen-16")

       firstCall_params = .false.
    endif


    ! determine the enthalpy that is consistent with our
    ! inlet rho, T, and X
    den_eos(1)  = dens_fuel
    temp_eos(1) = temp_fuel

    xn_eos(1,:) = ZERO
    xn_eos(1,ic12) = xc12_fuel
    xn_eos(1,io16) = 1.d0 - xc12_fuel

    call eos(eos_input_rt, den_eos, temp_eos, &
             npts, &
             xn_eos, &
             p_eos, h_eos, e_eos, &
             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
             dpdX_eos, dhdX_eos, &
             gam1_eos, cs_eos, s_eos, &
             dsdt_eos, dsdr_eos, &
             .false.)

    INLET_RHO     = dens_fuel
    INLET_RHOH    = dens_fuel*h_eos(1)
    INLET_TEMP    = temp_fuel
    INLET_RHOX(:) = dens_fuel*xn_eos(1,:)
    INLET_TRA     = ZERO

    INLET_VEL     = ZERO

    inlet_bc_initialized = .true.

  end subroutine set_initial_inlet_bcs

  subroutine update_inlet_bcs(time,dx,mla,s,rho_Hnuc,rho_omegadot,halftime)

    ! initialize the inflow boundary condition variables
    ! this is called at the top of each timestep, to allow for
    ! time-varying inlet BCs based on the current fluid state.

    use probin_module, ONLY: dens_fuel, temp_fuel, xc12_fuel, vel_fuel, &
         prob_hi_x, prob_lo_x
    use geometry, only: dm, nlevs
    use network, only: network_species_index

    real(kind=dp_t) :: time, dx(:,:)
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: rho_Hnuc(:)
    type(multifab) , intent(in   ) :: rho_omegadot(:)
    logical        , intent(in   ) :: halftime

    ! local variables
    real(kind=dp_t), pointer::  sp(:,:,:,:)
    real(kind=dp_t), pointer::  rhnp(:,:,:,:)
    real(kind=dp_t), pointer::  rwp(:,:,:,:)
    logical        , pointer::  mp(:,:,:,:)

    real(kind=dp_t) :: C12_flame_speed, C12_flame_speed_level, C12_flame_speed_local
    integer :: lo(dm),hi(dm)
    integer :: ng_s,ng_rhn,ng_rw
    integer :: n, i

    integer, save :: ic12, io16
    logical, save :: firstCall_params = .true.

    if (firstCall_params) then

       ic12 = network_species_index("carbon-12")
       io16 = network_species_index("oxygen-16")

       firstCall_params = .false.
    endif


    ! now determine the inlet velocity
    ! here, we compute an estimate of the flame speed
    ng_s = s(1)%ng
    ng_rhn = rho_Hnuc(1)%ng
    ng_rw = rho_omegadot(1)%ng


    ! initialize the flame speed
    C12_flame_speed = ZERO

    ! loop over levels and compute the flame speed
    do n = 1, nlevs

       ! initialize the local (processor's version) and level quantities to 0   
       C12_flame_speed_level = ZERO
       C12_flame_speed_local = ZERO

       ! loop over boxes in a given level                                       
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle
          sp => dataptr(s(n) , i)
          rhnp => dataptr(rho_Hnuc(n), i)
          rwp => dataptr(rho_omegadot(n), i)
          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))
          select case (dm)
          case (2)
             if (n .eq. nlevs) then
                call compute_inlet_vel_2d(n,time,dx(n,:), &
                                          sp(:,:,1,:),ng_s, &
                                          rhnp(:,:,1,1),ng_rhn, &
                                          rwp(:,:,1,:),ng_rw, &
                                          lo,hi, &
                                          C12_flame_speed_local)
             else
                mp => dataptr(mla%mask(n), i)
                call compute_inlet_vel_2d(n,time,dx(n,:), &
                                          sp(:,:,1,:),ng_s, &
                                          rhnp(:,:,1,1),ng_rhn, &
                                          rwp(:,:,1,:),ng_rw, &
                                          lo,hi, &
                                          C12_flame_speed_local, &
                                          mp(:,:,1,1))
             endif
          case (3)
             call bl_error("ERROR: 3D compute_inlet_vel not implemented")
          end select
       end do

       ! do a parallel reduce for the current level
       call parallel_reduce(C12_flame_speed_level, C12_flame_speed_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())


       ! reduce the current level's data with the global data                                
       if (parallel_IOProcessor()) then
          C12_flame_speed = C12_flame_speed + C12_flame_speed_level
       endif

    enddo

    ! normalize                                      

    ! our flame speed estimate is based on the carbon destruction rate  
    ! V_eff = - int { rho omegadot dx } / W (rho X)^in                      
    C12_flame_speed = -C12_flame_speed*dx(1,1)*dx(1,2)/ &
         ((prob_hi_x - prob_lo_x)*inlet_rhox(ic12))


    ! broadcast the flame_speed to all processors as the new
    ! inlet velocity

    if (halftime) then
       INLET_VEL = HALF*(C12_flame_speed + INLET_VEL_OLD)
    else
       INLET_VEL_OLD = INLET_VEL
       INLET_VEL = C12_flame_speed
    endif


!    INLET_VEL = vel_fuel

    if (parallel_IOProcessor()) then
       print *, 'inlet velocity = ', INLET_VEL
    endif

  end subroutine update_inlet_bcs

  subroutine compute_inlet_vel_2d(n,time,dx, &
                                  s,ng_s, &
                                  rho_Hnuc,ng_rhn, &
                                  rho_omegadot,ng_rw, &
                                  lo,hi, &
                                  C12_flame_speed, &
                                  mask)

    use bl_constants_module
    use network, only: nspec, network_species_index
    use probin_module, only: prob_lo

    integer, intent(in) :: n, lo(:), hi(:), ng_s, ng_rhn, ng_rw
    real (kind=dp_t), intent(in   ) :: time, dx(:)
    real (kind=dp_t), intent(in   ) ::        s(lo(1)-ng_s:,  lo(2)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho_Hnuc(lo(1)-ng_rhn:,lo(2)-ng_rhn:)
    real (kind=dp_t), intent(in   ) :: rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,:)
    real (kind=dp_t), intent(inout) :: C12_flame_speed
    logical,          intent(in   ), optional :: mask(lo(1):,lo(2):)

    ! local variables
    integer            :: i, j
    real (kind=dp_t)   :: weight
    logical            :: cell_valid
    real (kind=dp_t)   :: x, y

    integer, save :: ic12, io16

    logical, save :: firstCall = .true.

    if (firstCall) then

       ic12 = network_species_index("carbon-12")
       io16 = network_species_index("oxygen-16")

       firstCall = .false.
    endif

    ! weight is the factor by which the volume of a cell at the current level                
    ! relates to the volume of a cell at the coarsest level of refinement.                   
    weight = 1.d0 / 4.d0**(n-1)

    do j = lo(2), hi(2)
       y = prob_lo(2) + (dble(j) + HALF) * dx(2)

       do i = lo(1), hi(1)
          x = prob_lo(1) + (dble(i) + HALF) * dx(1)

          cell_valid = .true.
          if (present(mask)) then
             if ( (.not. mask(i,j)) ) cell_valid = .false.
          endif

          if (cell_valid) then

             ! compute the flame speed by integrating the carbon                             
             ! consumption rate.                                                             
             C12_flame_speed = C12_flame_speed + weight*rho_omegadot(i,j,ic12)
          endif

       enddo
    enddo

  end subroutine compute_inlet_vel_2d

end module inlet_bc_module
