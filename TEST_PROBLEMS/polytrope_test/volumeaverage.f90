! computes volume averaged quantities, such as energies

module volave_module

  use parallel, only: parallel_IOProcessor, parallel_IOProcessorNode, &
                      parallel_reduce, MPI_SUM
  use bl_types, only: dp_t
  use bl_IO_module, only: unit_new
  use box_module, only: lwb, upb
  use multifab_module, only: multifab, multifab_build, multifab_build_edge, &
                             dataptr, get_box, multifab_remote, &
                             nghost, setval, destroy
  use ml_layout_module, only: ml_layout
  use define_bc_module, only: bc_tower

  implicit none

  private

  public :: volave

contains

  subroutine volave(newtime,dx,s, &
                  u,w0,rho0,div_coeff, &
                  mla)

    use bl_prof_module, only: bl_prof_timer, build, destroy
    use geometry, only: spherical
    use bl_constants_module, only: ZERO
    use probin_module, only: prob_lo_x, prob_lo_y, prob_lo_z, &
                             prob_hi_x, prob_hi_y, prob_hi_z, &
                             job_name
    use network, only: network_species_index
    use fill_3d_module, only: put_1d_array_on_cart, make_w0mac
    use variables, only: foextrap_comp

    !**** inputs
    real(kind=dp_t), intent(in   ) :: dx(:,:),newtime
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: u(:)
    real(kind=dp_t), intent(in   ) :: w0(:,0:)
    real(kind=dp_t), intent(in   ) :: rho0(:,0:)
    real(kind=dp_t), intent(in   ) :: div_coeff(:,0:)

    !**** local variables
    real(kind=dp_t) :: KE_tot, KE_level, KE_local
    real(kind=dp_t) :: PE_tot, PE_level, PE_local
    real(kind=dp_t) :: PKE_tot, PKE_level, PKE_local
    real(kind=dp_t) :: PPE_tot, PPE_level, PPE_local
    real(kind=dp_t) :: KEnl_tot, KEnl_level, KEnl_local
    real(kind=dp_t) :: PEnl_tot, PEnl_level, PEnl_local
    real(kind=dp_t) :: IEnl_tot, IEnl_level, IEnl_local
    real(kind=dp_t), pointer::  sp(:,:,:,:)
    real(kind=dp_t), pointer::  up(:,:,:,:)
    logical        , pointer::  mp(:,:,:,:)
    real(kind=dp_t) :: volume
    integer :: dm, nlevs, n, i, lo(mla%dim), hi(mla%dim)
    integer :: ng_s, ng_u
    integer :: un
    logical :: lexist

    logical, save :: firstCall_io = .true.

    type(bl_prof_timer), save :: bpt

    call build(bpt, "volume_average")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_s   = nghost(s(1))
    ng_u   = nghost(u(1))

    !===========================================================================
    ! INITIALIZE
    !===========================================================================

    KE_tot = ZERO
    PE_tot = ZERO
    PKE_tot = ZERO
    PPE_tot = ZERO
    KEnl_tot = ZERO
    PEnl_tot = ZERO
    IEnl_tot = ZERO

    !===========================================================================
    ! Loop over levels
    !===========================================================================

    do n = 1, nlevs

       ! Set level and local totals to zero
       KE_level = ZERO
       KE_local = ZERO
       PE_level = ZERO
       PE_local = ZERO
       PKE_level = ZERO
       PKE_local = ZERO
       PPE_level = ZERO
       PPE_local = ZERO
       KEnl_level = ZERO
       KEnl_local = ZERO
       PEnl_level = ZERO
       PEnl_local = ZERO
       IEnl_level = ZERO
       IEnl_local = ZERO

       !========================================================================
       ! Loop over boxes in a level
       !========================================================================
       do i = 1, s(n)%nboxes

          if ( multifab_remote(s(n), i) ) cycle

          sp => dataptr(s(n), i)
          up => dataptr(u(n), i)

          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))

          select case (dm)

          case (2)
             if (n .eq. nlevs) then
                call volave_2d(n,dx(n,:), &
                             sp(:,:,1,:),ng_s, &
                             up(:,:,1,:),ng_u, &
                             w0(n,:),rho0(n,:),div_coeff(n,:), &
                             lo,hi, &
                             KE_local,PE_local,PKE_local,PPE_local, &
                             KEnl_local,PEnl_local,IEnl_local)
             else
                mp => dataptr(mla%mask(n), i)
                call volave_2d(n,dx(n,:), &
                             sp(:,:,1,:),ng_s, &
                             up(:,:,1,:),ng_u, &
                             w0(n,:),rho0(n,:),div_coeff(n,:), &
                             lo,hi, &
                             KE_local,PE_local,PKE_local,PPE_local, &
                             KEnl_local,PEnl_local,IEnl_local, &
                             mp(:,:,1,1))
             endif

          end select ! number of dimensions

       end do ! loop over boxes

       !----------------------------------------------------------------------
       ! sum up the quantities for the current level
       !----------------------------------------------------------------------

       ! NOTE: only the I/O Processor will have the correct reduced value
       call parallel_reduce(KE_level, KE_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(PE_level, PE_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(PKE_level, PKE_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(PPE_level, PPE_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(KEnl_level, KEnl_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(PEnl_level, PEnl_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       call parallel_reduce(IEnl_level, IEnl_local, MPI_SUM, &
                            proc = parallel_IOProcessorNode())

       !----------------------------------------------------------------------
       ! sum the current level quantity and the total quantity
       !----------------------------------------------------------------------
       if (parallel_IOProcessor()) then
          KE_tot = KE_tot + KE_level
          PE_tot = PE_tot + PE_level
          PKE_tot = PKE_tot + PKE_level
          PPE_tot = PPE_tot + PPE_level
          KEnl_tot = KEnl_tot + KEnl_level
          PEnl_tot = PEnl_tot + PEnl_level
          IEnl_tot = IEnl_tot + IEnl_level
       endif

    end do ! loop over levels


    !=========================================================================
    ! normalize
    !=========================================================================

    ! find volume

    select case (dm)

    case (2)

       volume = (prob_hi_x - prob_lo_x)*(prob_hi_y - prob_lo_y)

    end select ! number of dimensions

    KE_tot = KE_tot/volume
    PE_tot = PE_tot/volume
    PKE_tot = PKE_tot/volume
    PPE_tot = PPE_tot/volume
    KEnl_tot = KEnl_tot/volume
    PEnl_tot = PEnl_tot/volume
    IEnl_tot = IEnl_tot/volume

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
       inquire(file="volumeaverages.dat", exist=lexist)
       if (lexist) then
          open(unit=un, file="volumeaverages.dat", &
               status="old", position="append")
       else
          open(unit=un, file="volumeaverages.dat", status="new")
       endif


       ! write out the headers
       if (firstCall_io) then

          ! radvel
          write (un, *) " "
          write (un, 999) trim(job_name)
          write (un, 1001) "time", "Average pseudo-KE", "Average pseudo-PE", "Average KE", "Average PE", "Average Total KE", &
                            "Average Total PE", "Average Total IE"

          firstCall_io = .false.
       endif

       ! write out the data
       write (un,1000) newtime, PKE_tot, PPE_tot, KE_tot, PE_tot, KEnl_tot, PEnl_tot, IEnl_tot

       close(un)

    endif

    !=========================================================================
    ! clean-up
    !=========================================================================
    call destroy(bpt)

  end subroutine volave

  subroutine volave_2d(n,dx, &
                     s,ng_s, &
                     u,ng_u, &
                     w0,rho0,div_coeff, &
                     lo,hi, &
                     KE_local,PE_local,PKE_local,PPE_local, &
                     KEnl_local,PEnl_local,IEnl_local, &
                     mask)

    use variables, only: rho_comp, temp_comp
    use bl_constants_module, only: HALF, ZERO
    use probin_module, only: prob_lo, grav_const

    !**** inputs

    integer, intent(in) :: n, lo(:), hi(:), ng_s, ng_u
    real (kind=dp_t), intent(in   ) :: s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: u(lo(1)-ng_u:,lo(2)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: w0(0:)
    real (kind=dp_t), intent(in   ) :: rho0(0:)
    real (kind=dp_t), intent(in   ) :: div_coeff(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)
    logical,          intent(in   ), optional :: mask(lo(1):,lo(2):)

    !**** outputs

    real (kind=dp_t), intent(  out) :: KE_local
    real (kind=dp_t), intent(  out) :: PE_local
    real (kind=dp_t), intent(  out) :: PKE_local
    real (kind=dp_t), intent(  out) :: PPE_local
    real (kind=dp_t), intent(  out) :: KEnl_local
    real (kind=dp_t), intent(  out) :: PEnl_local
    real (kind=dp_t), intent(  out) :: IEnl_local

    !**** local variables

    real (kind=dp_t)   :: velsqr, velsqrlin
    integer            :: i, j
    real (kind=dp_t)   :: x, y
    logical            :: cell_valid
    real (kind=dp_t)   :: rho1
    real (kind=dp_t)   :: N2(lo(2):hi(2))
    real (kind=dp_t)   :: rho_over_beta(lo(2):hi(2))
    real (kind=dp_t)   :: entropy(lo(2):hi(2))

    ! initializing KE and PE

    KE_local = ZERO
    PE_local = ZERO
    PKE_local = ZERO
    PPE_local = ZERO
    KEnl_local = ZERO
    PEnl_local = ZERO
    IEnl_local = ZERO

    ! calculate N2

    do j = lo(2), hi(2)
       rho_over_beta(j) = rho0(j) / div_coeff(j)
       entropy(j) = log(rho_over_beta(j))
    enddo

    do j = lo(2), hi(2)-1
       N2(j) = grav_const*(entropy(j+1)-entropy(j))/dx(2)
    enddo
    N2(hi(2)) = N2(hi(2)-1)

    do j = lo(2), hi(2)
       y = prob_lo(2) + (dble(j) + HALF) * dx(2)

       do i = lo(1), hi(1)
          x = prob_lo(1) + (dble(i) + HALF) * dx(1)

          cell_valid = .true.
          if (present(mask)) then
             if ( (.not. mask(i,j)) ) cell_valid = .false.
          endif

          if (cell_valid) then

             velsqr = u(i,j,1)**2 + ( u(i,j,2) + HALF*(w0(j) + w0(j+1)) )**2

             rho1 = s(i,j,rho_comp) - rho0(j)

!**** nonlinear energies
             KEnl_local = KEnl_local + HALF*velsqr*s(i,j,rho_comp)*dx(1)*dx(2)
!             IEnl_local = IEnl_local + s(i,j,rho_comp)*s(i,j,temp_comp)*2.5d0*dx(1)*dx(2)
!             KEnl_local = KEnl_local - grav_const*y*rho0(j)*dx(1)*dx(2)
             IEnl_local = IEnl_local - grav_const*y*rho1*dx(1)*dx(2)
             PEnl_local = PEnl_local - grav_const*y*s(i,j,rho_comp)*dx(1)*dx(2)

!**** linear energies
             KE_local = KE_local + HALF*velsqr*rho0(j)*dx(1)*dx(2)
             PE_local = PE_local + HALF*grav_const**2*rho0(j)/(N2(j)) &
                        * (rho1/rho0(j))**2*dx(1)*dx(2)

!**** linear pseudo-energies
             PKE_local = PKE_local + HALF*velsqr*rho0(j)*div_coeff(j)*dx(1)*dx(2)
             PPE_local = PPE_local + HALF*rho0(j)*div_coeff(j)/(N2(j))*(rho1/rho0(j))**2*dx(1)*dx(2)

          endif ! cell valid

       enddo
    enddo

  end subroutine volave_2d

end module volave_module
