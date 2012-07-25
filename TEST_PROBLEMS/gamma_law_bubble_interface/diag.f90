! a generic diag module for MAESTRO.  This simply computes the 
! maximum Mach number on the domain and outputs it each timestep
! to maestro_diag.out

module diag_module

  use parallel, only: parallel_IOProcessor, parallel_IOProcessorNode, &
                      parallel_reduce, MPI_MAX, MPI_SUM
  use bl_types, only: dp_t
  use bl_IO_module, only: unit_new
  use bl_error_module, only: bl_error
  use box_module, only: lwb, upb
  use multifab_module, only: multifab, multifab_build, multifab_build_edge, &
                             dataptr, get_box, multifab_remote, &
                             nghost, setval, destroy
  use ml_layout_module, only: ml_layout
  use define_bc_module, only: bc_tower

  implicit none

  private

  public :: diag, flush_diag

contains

  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  subroutine diag(newtime,dt,dx,s,rho_Hnuc,rho_Hext,thermal,rho_omegadot, &
                  rho0,rhoh0,p0,tempbar, &
                  gamma1bar,div_coeff, &
                  u,w0,normal, &
                  mla,the_bc_tower)

    use bl_prof_module, only: bl_prof_timer, build, destroy
    use geometry, only: spherical
    use bl_constants_module, only: ZERO
    use probin_module, only: prob_lo_x, prob_lo_y, prob_lo_z, &
                             prob_hi_x, prob_hi_y, prob_hi_z, &
                             job_name
    use network, only: network_species_index
    use fill_3d_module, only: put_1d_array_on_cart, make_w0mac
    use variables, only: foextrap_comp

    real(kind=dp_t), intent(in   ) :: dt,dx(:,:),newtime
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: rho_Hnuc(:)
    type(multifab) , intent(in   ) :: rho_Hext(:)
    type(multifab) , intent(in   ) :: thermal(:)
    type(multifab) , intent(in   ) :: rho_omegadot(:)
    type(multifab) , intent(in   ) :: u(:)
    type(multifab) , intent(in   ) :: normal(:)
    real(kind=dp_t), intent(in   ) ::      rho0(:,0:)
    real(kind=dp_t), intent(in   ) ::     rhoh0(:,0:)
    real(kind=dp_t), intent(in   ) ::        p0(:,0:)
    real(kind=dp_t), intent(in   ) ::   tempbar(:,0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(:,0:)
    real(kind=dp_t), intent(in   ) :: div_coeff(:,0:)
    real(kind=dp_t), intent(in   ) ::        w0(:,0:)
    type(ml_layout), intent(in   ) :: mla
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! Local
    real(kind=dp_t), pointer::  sp(:,:,:,:)
    real(kind=dp_t), pointer::  rhnp(:,:,:,:)
    real(kind=dp_t), pointer::  rhep(:,:,:,:)
    real(kind=dp_t), pointer::  rwp(:,:,:,:)
    real(kind=dp_t), pointer::  up(:,:,:,:)
    real(kind=dp_t), pointer :: nop(:,:,:,:)
    real(kind=dp_t), pointer :: w0rp(:,:,:,:)
    real(kind=dp_t), pointer :: w0xp(:,:,:,:)
    real(kind=dp_t), pointer :: w0yp(:,:,:,:)
    real(kind=dp_t), pointer :: w0zp(:,:,:,:)
    logical        , pointer::  mp(:,:,:,:)

    type(multifab) :: w0r_cart(mla%nlevel)
    type(multifab) ::    w0mac(mla%nlevel,mla%dim)

    real(kind=dp_t) :: Mach_max, Mach_max_level, Mach_max_local
    real(kind=dp_t) :: temp_max, temp_max_level, temp_max_local
    real(kind=dp_t) :: enuc_max, enuc_max_level, enuc_max_local
    real(kind=dp_t) :: Hext_max, Hext_max_level, Hext_max_local

    real(kind=dp_t) :: kin_ener, kin_ener_level, kin_ener_local
    real(kind=dp_t) :: int_ener, int_ener_level, int_ener_local
    real(kind=dp_t) :: pot_ener, pot_ener_level, pot_ener_local

    real(kind=dp_t) :: sum_data_level(3), sum_data_local(3)


    integer :: lo(mla%dim),hi(mla%dim),dm,nlevs
    integer :: ng_s,ng_u,ng_n,ng_rhn,ng_rhe,ng_rw,ng_w,ng_wm
    integer :: i,n,comp
    integer :: un
    logical :: lexist

    logical, save :: firstCall_io = .true.

    type(bl_prof_timer), save :: bpt

    call build(bpt, "diagnostics")

    dm = mla%dim
    nlevs = mla%nlevel

    if (spherical .eq. 1) then

       do n=1,nlevs
          
          do comp=1,dm
             ! w0mac will contain an edge-centered w0 on a Cartesian
             ! grid, for use in computing divergences.
             call multifab_build_edge(w0mac(n,comp), mla%la(n),1,1,comp)
             call setval(w0mac(n,comp), ZERO, all=.true.)
          enddo

          ! w0r_cart is w0 but onto a Cartesian grid in cell-centered
          ! as a scalar.  Since w0 is the radial expansion velocity,
          ! w0r_cart is the radial w0 in a zone
          call multifab_build(w0r_cart(n), mla%la(n),1,0)
          call setval(w0r_cart(n), ZERO, all=.true.)
       end do

       ! put w0 on Cartesian edges as a vector
       call make_w0mac(mla,w0,w0mac,dx,the_bc_tower%bc_tower_array)


       ! put w0 in Cartesian cell-centers as a scalar (the radial
       ! expansion velocity)
       call put_1d_array_on_cart(w0,w0r_cart,1,.true.,.false.,dx, &
                                 the_bc_tower%bc_tower_array,mla)
    endif

    ng_s   = nghost(s(1))
    ng_u   = nghost(u(1))
    ng_n   = nghost(normal(1))
    ng_rhn = nghost(rho_Hnuc(1))
    ng_rhe = nghost(rho_Hext(1))
    ng_rw  = nghost(rho_omegadot(1))
    if (spherical == 1) then
       ng_n   = nghost(normal(1))
       ng_w   = nghost(w0r_cart(1))
       ng_wm  = nghost(w0mac(1,1))
    endif


    !=========================================================================
    ! initialize
    !=========================================================================
    Mach_max = ZERO
    temp_max = ZERO
    enuc_max = ZERO
    Hext_max = ZERO

    kin_ener = ZERO
    int_ener = ZERO
    pot_ener = ZERO


    !=========================================================================
    ! loop over the levels and compute the global quantities
    !=========================================================================
    do n = 1, nlevs

       ! initialize the local (processor's version) and level quantities to 0
       Mach_max_level = ZERO
       Mach_max_local = ZERO

       temp_max_level = ZERO
       temp_max_local = ZERO

       enuc_max_level = ZERO
       enuc_max_local = ZERO

       Hext_max_level = ZERO
       Hext_max_local = ZERO

       kin_ener_level = ZERO
       kin_ener_local = ZERO

       int_ener_level = ZERO
       int_ener_local = ZERO

       pot_ener_level = ZERO
       pot_ener_local = ZERO


       !----------------------------------------------------------------------
       ! loop over boxes in a given level
       !----------------------------------------------------------------------
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle

          sp => dataptr(s(n) , i)
          rhnp => dataptr(rho_Hnuc(n), i)
          rhep => dataptr(rho_Hext(n), i)
          rwp => dataptr(rho_omegadot(n), i)
          up => dataptr(u(n) , i)

          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))

          select case (dm)

          case (1)
             call bl_error("diag not implemented in 1-d")

          case (2)
             if (n .eq. nlevs) then
                call diag_2d(n,newtime,dt,dx(n,:), &
                             sp(:,:,1,:),ng_s, &
                             rhnp(:,:,1,1),ng_rhn, &
                             rhep(:,:,1,1),ng_rhe, &
                             rwp(:,:,1,:),ng_rw, &
                             rho0(n,:),rhoh0(n,:), &
                             p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                             up(:,:,1,:),ng_u, &
                             w0(n,:), &
                             lo,hi, &
                             Mach_max_local,temp_max_local, &
                             enuc_max_local,Hext_max_local, &
                             kin_ener_local,int_ener_local,pot_ener_local)
             else
                mp => dataptr(mla%mask(n), i)
                call diag_2d(n,newtime,dt,dx(n,:), &
                             sp(:,:,1,:),ng_s, &
                             rhnp(:,:,1,1),ng_rhn, &
                             rhep(:,:,1,1),ng_rhe, &
                             rwp(:,:,1,:),ng_rw, &
                             rho0(n,:),rhoh0(n,:), &
                             p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                             up(:,:,1,:),ng_u, &
                             w0(n,:), &
                             lo,hi, &
                             Mach_max_local,temp_max_local, &
                             enuc_max_local,Hext_max_local, &
                             kin_ener_local,int_ener_local,pot_ener_local, &
                             mp(:,:,1,1))
             endif

          case (3)
             call bl_error("diag not implemented in 3-d")

          end select
       end do

       !----------------------------------------------------------------------
       ! do the appropriate parallel reduction for the current level
       !----------------------------------------------------------------------

       ! NOTE: only the I/O Processor will have the correct reduced value
       call parallel_reduce(Mach_max_level, Mach_max_local, MPI_MAX, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(temp_max_level, temp_max_local, MPI_MAX, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(enuc_max_level, enuc_max_local, MPI_MAX, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(Hext_max_level, Hext_max_local, MPI_MAX, &
                            proc = parallel_IOProcessorNode())


       sum_data_local(1) = kin_ener_local
       sum_data_local(2) = int_ener_local
       sum_data_local(3) = pot_ener_local

       call parallel_reduce(sum_data_level, sum_data_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       kin_ener_level = sum_data_level(1)
       int_ener_level = sum_data_level(2)
       pot_ener_level = sum_data_level(3)



       !----------------------------------------------------------------------
       ! reduce the current level's data with the global data
       !----------------------------------------------------------------------
       if (parallel_IOProcessor()) then
          Mach_max = max(Mach_max, Mach_max_level)
          temp_max = max(temp_max, temp_max_level)
          enuc_max = max(enuc_max, enuc_max_level)
          Hext_max = max(Hext_max, Hext_max_level)

          kin_ener = kin_ener + kin_ener_level
          int_ener = int_ener + int_ener_level
          pot_ener = pot_ener + pot_ener_level

       endif

    end do


    !=========================================================================
    ! normalize 
    !=========================================================================

    ! normalize any integral quantities here

    ! the volume we normalize with is that of a single coarse-level          
    ! zone.  This is because the weight used in the loop over cells          
    ! was with reference to the coarse level                                 

    if (dm == 2) then
       kin_ener = kin_ener*dx(1,1)*dx(1,2)
       int_ener = int_ener*dx(1,1)*dx(1,2)
       pot_ener = pot_ener*dx(1,1)*dx(1,2)
    else
       call bl_error("only 2-d supported in diag.f90")
    endif


    !=========================================================================
    ! output
    !=========================================================================
 999 format("# job name: ",a)
1000 format(1x,10(g24.10,1x))
1001 format("#",10(a24,1x))

    if (parallel_IOProcessor()) then

       ! open the diagnostic files for output, taking care not to overwrite
       ! an existing file
       un = unit_new()
       inquire(file="maestro_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un, file="maestro_diag.out", &
               status="old", position="append")
       else
          open(unit=un, file="maestro_diag.out", status="new")
       endif


       ! write out the headers
       if (firstCall_io) then

          ! radvel
          write (un, *) " "
          write (un, 999) trim(job_name)
          write (un, 1001) "time", "max Mach #", "max T (K)", "max enuc (erg/g/s)", "max Hext (erg/g/s)", &
               "kin energy (erg)", "int energy (erg)", "pot energy (erg)"

          firstCall_io = .false.
       endif

       ! write out the data
       write (un,1000) newtime, Mach_max, temp_max, enuc_max, Hext_max, &
            kin_ener, int_ener, pot_ener

       close(un)

    endif
    
    !=========================================================================
    ! clean-up
    !=========================================================================
    call destroy(bpt)

    if (spherical.eq.1) then
       do n = 1, nlevs
          call destroy(w0r_cart(n))
          do comp=1,dm
             call destroy(w0mac(n,comp))
          enddo
       end do
    endif

  end subroutine diag


  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  subroutine flush_diag()
    ! flush_diag is called immediately before checkpointing.  If an
    ! implementation of these diagnostic routines wants to buffer the
    ! data a write out a lot of timestep's worth of information all 
    ! at once, flush_diag() is the routine that should do the writing.

  end subroutine flush_diag



  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  subroutine diag_2d(n,newtime,dt,dx, &
                     s,ng_s, &
                     rho_Hnuc,ng_rhn, &
                     rho_Hext,ng_rhe, &
                     rho_omegadot,ng_rw, &
                     rho0,rhoh0,p0,tempbar,gamma1bar, &
                     u,ng_u, &
                     w0, &
                     lo,hi, &
                     Mach_max,temp_max, &
                     enuc_max,Hext_max, &
                     kin_ener,int_ener,pot_ener, &
                     mask)

    use variables, only: rho_comp, spec_comp, temp_comp
    use bl_constants_module, only: HALF
    use network, only: nspec
    use probin_module, only: prob_lo, grav_const
    use eos_module, only: eos_input_rt, eos
    use eos_type_module

    integer, intent(in) :: n, lo(:), hi(:), ng_s, ng_u, ng_rhn, ng_rhe, ng_rw
    real (kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho_Hnuc(lo(1)-ng_rhn:,lo(2)-ng_rhn:)
    real (kind=dp_t), intent(in   ) :: rho_Hext(lo(1)-ng_rhe:,lo(2)-ng_rhe:)
    real (kind=dp_t), intent(in   ) :: rho_omegadot(lo(1)-ng_rw:,lo(2)-ng_rw:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:), rhoh0(0:), &
                                         p0(0:),tempbar(0:),gamma1bar(0:)
    real (kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: w0(0:)
    real (kind=dp_t), intent(in   ) :: newtime, dt, dx(:)
    real (kind=dp_t), intent(inout) :: Mach_max, temp_max, enuc_max, Hext_max
    real (kind=dp_t), intent(inout) :: kin_ener, int_ener, pot_ener
    logical,          intent(in   ), optional :: mask(lo(1):,lo(2):)

    !     Local variables
    integer            :: i, j
    real (kind=dp_t)   :: weight
    logical            :: cell_valid
    real (kind=dp_t)   :: x, y
    real (kind=dp_t)   :: vel

    type (eos_t) :: eos_state


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

             ! vel is the magnitude of the velocity, including w0
             vel = sqrt(  u(i,j,1)**2 + &
                        ( u(i,j,2) + HALF*(w0(j) + w0(j+1)) )**2 )

             
             ! call the EOS to get the sound speed and internal energy       
             eos_state%T     = s(i,j,temp_comp)
             eos_state%rho   = s(i,j,rho_comp)
             eos_state%xn(:) = s(i,j,spec_comp:spec_comp+nspec-1)/eos_state%rho

             call eos(eos_input_rt, eos_state, .false.)


             ! max Mach number                                       
             Mach_max = max(Mach_max,vel/eos_state%cs)

             ! max temp and enuc
             temp_max = max(temp_max,s(i,j,temp_comp))
             enuc_max = max(enuc_max,rho_Hnuc(i,j)/s(i,j,rho_comp))
             Hext_max = max(Hext_max,rho_Hext(i,j)/s(i,j,rho_comp))

             ! energies
             kin_ener = kin_ener + weight*s(i,j,rho_comp)*vel**2
             int_ener = int_ener + weight*s(i,j,rho_comp)*eos_state%e
             pot_ener = pot_ener + weight*s(i,j,rho_comp)*grav_const*y  ! zero-point is arbitrary


          endif  ! cell valid

       enddo
    enddo

  end subroutine diag_2d

end module diag_module
