! test_he_burn diagnostic routine
!
! output files:
!
!    test_he_burn_enuc_diag.out:
!        peak nuclear energy generation rate (erg / g / s)
!        x/y/z location of peak 
!        max ratio of cooling / nuc energy generation
!        x/y/z location of max
!
!    test_he_burn_temp_diag.out:
!        peak temperature (K)
!        x/y/z location of peak
!



module diag_module

  use bl_types
  use bl_IO_module
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: diag, flush_diag

contains

  subroutine diag(time,dt,dx,s,rho_Hnuc,rho_Hext,thermal,rho_omegadot, &
                  rho0,rhoh0,p0,tempbar, &
                  gamma1bar,div_coeff, &
                  u,w0,normal, &
                  mla,the_bc_tower)

    use bl_prof_module
    use bl_constants_module

    real(kind=dp_t), intent(in   ) :: dt,dx(:,:),time
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(in   ) :: rho_Hnuc(:)
    type(multifab) , intent(in   ) :: rho_Hext(:)
    type(multifab) , intent(in   ) :: rho_omegadot(:)    
    type(multifab) , intent(in   ) :: thermal(:)
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
    real(kind=dp_t), pointer::  thp(:,:,:,:)
    real(kind=dp_t), pointer::  up(:,:,:,:)
    logical,         pointer :: mp(:,:,:,:)

    real(kind=dp_t) :: T_max, T_max_local, T_max_level
    real(kind=dp_t) :: x_Tmax_local, y_Tmax_local, z_Tmax_local
    real(kind=dp_t) :: x_Tmax_level, y_Tmax_level, z_Tmax_level
    real(kind=dp_t) :: x_Tmax, y_Tmax, z_Tmax
    ! buffers
    real(kind=dp_t) :: T_max_data_local(1), T_max_coords_local(mla%dim)
    real(kind=dp_t), allocatable :: T_max_data(:), T_max_coords(:)

    real(kind=dp_t) :: enuc_max, enuc_max_local, enuc_max_level
    real(kind=dp_t) :: x_enucmax_local, y_enucmax_local, z_enucmax_local
    real(kind=dp_t) :: x_enucmax_level, y_enucmax_level, z_enucmax_level
    real(kind=dp_t) :: x_enucmax, y_enucmax, z_enucmax
    ! buffers
    real(kind=dp_t) :: enuc_max_data_local(1), enuc_max_coords_local(mla%dim)
    real(kind=dp_t), allocatable :: enuc_max_data(:), enuc_max_coords(:)

    real(kind=dp_t) :: ratio_max, ratio_max_local, ratio_max_level
    real(kind=dp_t) :: x_ratiomax_local, y_ratiomax_local, z_ratiomax_local
    real(kind=dp_t) :: x_ratiomax_level, y_ratiomax_level, z_ratiomax_level
    real(kind=dp_t) :: x_ratiomax, y_ratiomax, z_ratiomax
    ! buffers
    real(kind=dp_t) :: ratio_max_data_local(1), ratio_max_coords_local(mla%dim)
    real(kind=dp_t), allocatable :: ratio_max_data(:), ratio_max_coords(:)

    integer :: lo(mla%dim),hi(mla%dim),ng_s,ng_u,ng_rhn,ng_rhe, ng_th, dm,nlevs
    integer :: i,n, index_max
    integer :: un, un2
    logical :: lexist

    logical, save :: firstCall = .true.

    type(bl_prof_timer), save :: bpt

    call build(bpt, "diagnostics")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_s = s(1)%ng
    ng_u = u(1)%ng
    ng_rhn = rho_Hnuc(1)%ng
    ng_rhe = rho_Hext(1)%ng
    ng_th  = thermal(1)%ng

    ! initialize
    T_max       = ZERO

    x_Tmax      = ZERO
    y_Tmax      = ZERO
    z_Tmax      = ZERO

    enuc_max    = ZERO

    x_enucmax   = ZERO
    y_enucmax   = ZERO
    z_enucmax   = ZERO

    ratio_max   = ZERO
    
    x_ratiomax  = ZERO
    y_ratiomax  = ZERO
    z_ratiomax  = ZERO
    

    ! loop over the levels and calculate global quantities
    do n = 1, nlevs

       ! initialize local and level quantities
       T_max_local = ZERO
       T_max_level = ZERO

       x_Tmax_local = ZERO
       y_Tmax_local = ZERO
       z_Tmax_local = ZERO
       x_Tmax_level = ZERO
       y_Tmax_level = ZERO
       z_Tmax_level = ZERO       

       enuc_max_local = ZERO
       enuc_max_level = ZERO

       x_enucmax_local = ZERO
       y_enucmax_local = ZERO
       z_enucmax_local = ZERO
       x_enucmax_level = ZERO
       y_enucmax_level = ZERO
       z_enucmax_level = ZERO

       ratio_max_local = ZERO
       ratio_max_level = ZERO

       x_ratiomax_local = ZERO
       y_ratiomax_local = ZERO
       z_ratiomax_local = ZERO
       x_ratiomax_level = ZERO
       y_ratiomax_level = ZERO
       z_ratiomax_level = ZERO


       ! loop over the boxes at the current level
       do i = 1, s(n)%nboxes

          if ( multifab_remote(s(n), i) ) cycle

          sp => dataptr(s(n) , i)
          rhnp => dataptr(rho_Hnuc(n), i)
          rhep => dataptr(rho_Hext(n), i)
          thp => dataptr(thermal(n), i)
          up => dataptr(u(n) , i)
          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))

          ! build the local data
          select case (dm)
          case (2)
             ! only do those boxes that aren't masked
             if (n .eq. nlevs) then
                
                call diag_2d(n,time,dt,dx(n,:), &
                             sp(:,:,1,:),ng_s, &
                             rhnp(:,:,1,1),ng_rhn, &
                             rhep(:,:,1,1),ng_rhe, &
                             thp(:,:,1,1), ng_th, &
                             rho0(n,:),rhoh0(n,:), &
                             p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                             up(:,:,1,:),ng_u, &
                             w0(n,:), &
                             lo,hi, &
                             T_max_local, &
                             x_Tmax_local, y_Tmax_local, &
                             enuc_max_local, &
                             x_enucmax_local, y_enucmax_local, &
                             ratio_max_local, &
                             x_ratiomax_local, y_ratiomax_local)
             else
                mp => dataptr(mla%mask(n),i)
                call diag_2d(n,time,dt,dx(n,:), &
                             sp(:,:,1,:),ng_s, &
                             rhnp(:,:,1,1),ng_rhn, &
                             rhep(:,:,1,1),ng_rhe, &
                             thp(:,:,1,1), ng_th, &
                             rho0(n,:),rhoh0(n,:), &
                             p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                             up(:,:,1,:),ng_u, &
                             w0(n,:), &
                             lo,hi, &
                             T_max_local, &
                             x_Tmax_local, y_Tmax_local, &
                             enuc_max_local, &
                             x_enucmax_local, y_enucmax_local, &
                             ratio_max_local, &
                             x_ratiomax_local, y_ratiomax_local, &
                             mp(:,:,1,1))
             endif
          case (3)
             ! only do those boxes that aren't masked
             if (n .eq. nlevs) then

                call diag_3d(n,time,dt,dx(n,:), &
                             sp(:,:,:,:),ng_s, &
                             rhnp(:,:,:,1),ng_rhn, &
                             rhep(:,:,:,1),ng_rhe, &
                             thp(:,:,:,1), ng_th, &
                             rho0(n,:),rhoh0(n,:), &
                             p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                             up(:,:,:,:),ng_u, &
                             w0(n,:), &
                             lo,hi, &
                             T_max_local, &
                             x_Tmax_local, y_Tmax_local, z_Tmax_local, &
                             enuc_max_local, &
                             x_enucmax_local, y_enucmax_local, z_enucmax_local, &
                             ratio_max_local, &
                             x_ratiomax_local, y_ratiomax_local, z_ratiomax_local)

             else
                mp => dataptr(mla%mask(n),i)
                call diag_3d(n,time,dt,dx(n,:), &
                             sp(:,:,:,:),ng_s, &
                             rhnp(:,:,:,1),ng_rhn, &
                             rhep(:,:,:,1),ng_rhe, &
                             thp(:,:,:,1), ng_th, &
                             rho0(n,:),rhoh0(n,:), &
                             p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                             up(:,:,:,:),ng_u, &
                             w0(n,:), &
                             lo,hi, &
                             T_max_local, &
                             x_Tmax_local, y_Tmax_local, z_Tmax_local, &
                             enuc_max_local, &
                             x_enucmax_local, y_enucmax_local, z_enucmax_local,&
                             ratio_max_local, &
                             x_ratiomax_local, y_ratiomax_local, z_ratiomax_local, &
                             mp(:,:,:,1))
             endif
          end select
       end do

       ! build the level data
       
       ! get the maximum temperature and its location
       ! gather all the T_max data into an array; find the index of the maximum
       allocate(T_max_data(parallel_nprocs()))
       T_max_data_local(1) = T_max_local

       call parallel_gather(T_max_data_local, T_max_data, 1, &
                            root = parallel_IOProcessorNode())
       
       index_max = maxloc(T_max_data, dim=1)

       ! gather all the T_max_coords into an array and use index_max to 
       ! get the correct location
       allocate(T_max_coords(dm*parallel_nprocs()))
       T_max_coords_local(1) = x_Tmax_local
       T_max_coords_local(2) = y_Tmax_local
       if (dm > 2) T_max_coords_local(3) = z_Tmax_local

       call parallel_gather(T_max_coords_local , T_max_coords, dm, &
                            root = parallel_IOProcessorNode())

       T_max_level = T_max_data(index_max)

       x_Tmax_level = T_max_coords(dm*(index_max-1) + 1)
       y_Tmax_level = T_max_coords(dm*(index_max-1) + 2)
       if (dm > 2) z_Tmax_level = T_max_coords(dm*(index_max-1) + 3)

       deallocate(T_max_data, T_max_coords)

       ! get the maximum enuc and its location
       allocate(enuc_max_data(parallel_nprocs()))
       enuc_max_data_local(1) = enuc_max_local

       call parallel_gather(enuc_max_data_local, enuc_max_data, 1, &
                            root = parallel_IOProcessorNode())

       index_max = maxloc(enuc_max_data, dim=1)

       allocate(enuc_max_coords(dm*parallel_nprocs()))
       enuc_max_coords_local(1) = x_enucmax_local
       enuc_max_coords_local(2) = y_enucmax_local
       if (dm > 2) enuc_max_coords_local(3) = z_enucmax_local

       call parallel_gather(enuc_max_coords_local, enuc_max_coords, dm, &
                            root = parallel_IOProcessorNode())

       enuc_max_level = enuc_max_data(index_max)

       x_enucmax_level = enuc_max_coords(dm*(index_max - 1) + 1)
       y_enucmax_level = enuc_max_coords(dm*(index_max - 1) + 2)
       if (dm > 2) z_enucmax_level = enuc_max_coords(dm*(index_max - 1) + 3)

       deallocate(enuc_max_data, enuc_max_coords)

       ! get the maximum cooling to energy generation ratio and its location
       allocate(ratio_max_data(parallel_nprocs()))
       ratio_max_data_local(1) = ratio_max_local

       call parallel_gather(ratio_max_data_local, ratio_max_data, 1, &
                            root = parallel_IOProcessorNode())

       index_max = maxloc(ratio_max_data, dim=1)

       allocate(ratio_max_coords(dm*parallel_nprocs()))
       ratio_max_coords_local(1) = x_ratiomax_local
       ratio_max_coords_local(2) = y_ratiomax_local
       if (dm > 2) ratio_max_coords_local(3) = z_ratiomax_local
       
       call parallel_gather(ratio_max_coords_local,ratio_max_coords, dm, &
                            root = parallel_IOProcessorNode())

       ratio_max_level = ratio_max_data(index_max)

       x_ratiomax_level = ratio_max_coords(dm*(index_max - 1) + 1)
       y_ratiomax_level = ratio_max_coords(dm*(index_max - 1) + 2)
       if (dm > 2) z_ratiomax_level = ratio_max_coords(dm*(index_max - 1) + 3)

       deallocate(ratio_max_data, ratio_max_coords)


       ! reduce the current level's data with the global data
       if (parallel_IOProcessor()) then

          if (T_max_level > T_max) then
             T_max = T_max_level

             x_Tmax = x_Tmax_level
             y_Tmax = y_Tmax_level
             if (dm > 2) z_Tmax = z_Tmax_level

          endif

          if (enuc_max_level > enuc_max) then
             enuc_max = enuc_max_level

             x_enucmax = x_enucmax_level
             y_enucmax = y_enucmax_level
             if (dm > 2) z_enucmax = z_enucmax_level

          endif

          if (ratio_max_level > ratio_max) then
             ratio_max = ratio_max_level

             x_ratiomax = x_ratiomax_level
             y_ratiomax = y_ratiomax_level
             if (dm > 2) z_ratiomax = z_ratiomax_level
          endif

       endif

    end do


1000 format(1x,10(g20.10,1x))
1001 format("#",10(a18,1x))

    if (parallel_IOProcessor()) then

       ! open the diagnostic files for output, taking care not to overwrite
       ! an existing file
       un = unit_new()
       inquire(file="test_he_burn_temp_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un, file="test_he_burn_temp_diag.out", &
               status="old", position="append")
       else
          open(unit=un, file="test_he_burn_temp_diag.out", status="new")
       endif

       un2 = unit_new()
       inquire(file="test_he_burn_enuc_diag.out", exist=lexist)
       if (lexist) then
          open(unit=un2, file="test_he_burn_enuc_diag.out", &
               status="old", position="append")
       else
          open(unit=un2, file="test_he_burn_enuc_diag.out", status="new")
       endif

       ! print the headers
       if(firstCall) then
          if (dm > 2) then
             write(un , 1001) "time", "max{T}", "x_loc", "y_loc", "z_loc"
             write(un2, 1001) "time", "max{enuc}", "x_loc", "y_loc", "z_loc", &
                              "max{cool/nuc}", "x_loc", "y_loc", "z_loc"
          else
             write(un , 1001) "time", "max{T}", "x_loc", "y_loc"
             write(un2, 1001) "time", "max{enuc}", "x_loc", "y_loc", &
                              "max{cool/nuc}", "x_loc", "y_loc"
          endif

          firstCall = .false.
       endif

       ! print the data
       if (dm > 2) then
          write(un , 1000) time, T_max, x_Tmax, y_Tmax, z_Tmax
          write(un2, 1000) time, enuc_max, x_enucmax, y_enucmax, z_enucmax, &
                           ratio_max, x_ratiomax, y_ratiomax, z_ratiomax
       else
          write(un , 1000) time, T_max, x_Tmax, y_Tmax
          write(un2, 1000) time, enuc_max, x_enucmax, y_enucmax, &
                           ratio_max, x_ratiomax, y_ratiomax
       endif

       close(un )
       close(un2)
    endif

    call destroy(bpt)

  end subroutine diag

  subroutine diag_2d(n,time,dt,dx, &
                     s,ng_s, &
                     rho_Hnuc,ng_rhn, &
                     rho_Hext,ng_rhe, &
                     thermal, ng_th, &
                     rho0,rhoh0,p0,tempbar,gamma1bar, &
                     u,ng_u, &
                     w0, &
                     lo,hi, &
                     T_max, &
                     x_Tmax, y_Tmax, &
                     enuc_max, &
                     x_enucmax, y_enucmax, &
                     ratio_max, &
                     x_ratiomax, y_ratiomax, &
                     mask)

    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp
    use network, only: nspec
    use probin_module, only: prob_lo, base_cutoff_density
    use bl_constants_module, only: HALF, ZERO

    integer, intent(in) :: n, lo(:), hi(:), ng_s, ng_u, ng_rhn, ng_rhe, ng_th
    real (kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho_Hnuc(lo(1)-ng_rhn:,lo(2)-ng_rhn:)
    real (kind=dp_t), intent(in   ) :: rho_Hext(lo(1)-ng_rhe:,lo(2)-ng_rhe:)
    real (kind=dp_t), intent(in   ) :: thermal(lo(1)-ng_th:,lo(2)-ng_th:)
    real (kind=dp_t), intent(in   ) :: rho0(0:), rhoh0(0:), &
                                         p0(0:),tempbar(0:),gamma1bar(0:)
    real (kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: w0(0:)
    real (kind=dp_t), intent(in   ) :: time, dt, dx(:)
    real (kind=dp_t), intent(inout) :: T_max, enuc_max
    real (kind=dp_t), intent(inout) :: x_Tmax, y_Tmax
    real (kind=dp_t), intent(inout) :: x_enucmax, y_enucmax
    real (kind=dp_t), intent(inout) :: ratio_max
    real (kind=dp_t), intent(inout) :: x_ratiomax, y_ratiomax
    logical,          intent(in   ), optional :: mask(lo(1):,lo(2):)

    !     Local variables
    integer                     :: i, j
    real (kind=dp_t) :: x, y 
    real (kind=dp_t) :: enuc_local, ratio_local
    logical :: cell_valid

    do j = lo(2), hi(2)
       y = prob_lo(2) + (dble(j) + HALF) * dx(2)

       do i = lo(1), hi(1)
          x = prob_lo(1) + (dble(i) + HALF) * dx(1)

          cell_valid = .true.
          if (present(mask)) then
             if ( (.not. mask(i,j)) ) cell_valid = .false.
          endif

          if (cell_valid .and. s(i,j,rho_comp) > base_cutoff_density) then

             ! temperature diagnostic
             if (s(i,j,temp_comp) > T_max) then
                
                T_max = s(i,j,temp_comp)
                
                x_Tmax = x
                y_Tmax = y

             endif

             ! enuc diagnostic
             enuc_local = rho_Hnuc(i,j)/s(i,j,rho_comp)
             if (enuc_local > enuc_max) then

                enuc_max = enuc_local
             
                x_enucmax = x
                y_enucmax = y

             endif

             ! ratio cooling / nuc. energy generation
             ! cooling is where thermal < ZERO
             ratio_local = ZERO

             if (thermal(i,j) < ZERO) then
                ratio_local = abs(thermal(i,j)) / rho_Hnuc(i,j)
             endif
             
             if (ratio_local > ratio_max) then

                ratio_max = ratio_local

                x_ratiomax = x
                y_ratiomax = y

             endif

          endif

       enddo
    enddo

  end subroutine diag_2d

  subroutine diag_3d(n,time,dt,dx, &
                     s,ng_s, &
                     rho_Hnuc,ng_rhn, &
                     rho_Hext,ng_rhe, &
                     thermal, ng_th, &
                     rho0,rhoh0,p0,tempbar,gamma1bar, &
                     u,ng_u,w0, &
                     lo,hi, &
                     T_max, &
                     x_Tmax, y_Tmax, z_Tmax, &
                     enuc_max, &
                     x_enucmax, y_enucmax, z_enucmax, &
                     ratio_max, &
                     x_ratiomax, y_ratiomax, z_ratiomax, &
                     mask)

    use variables, only: rho_comp, spec_comp, temp_comp, rhoh_comp
    use network, only: nspec
    use bl_constants_module, only: HALF, ZERO
    use probin_module, only: prob_lo, base_cutoff_density

    integer, intent(in) :: n,lo(:), hi(:), ng_s, ng_u, ng_rhn, ng_rhe, ng_th
    real (kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho_Hnuc(lo(1)-ng_rhn:,lo(2)-ng_rhn:,lo(3)-ng_rhn:)
    real (kind=dp_t), intent(in   ) :: rho_Hext(lo(1)-ng_rhe:,lo(2)-ng_rhe:,lo(3)-ng_rhe:)
    real (kind=dp_t), intent(in   ) :: thermal(lo(1)-ng_th:,lo(2)-ng_th:,lo(3)-ng_th:)
    real (kind=dp_t), intent(in   ) :: rho0(0:), rhoh0(0:), &
                                         p0(0:),tempbar(0:),gamma1bar(0:)
    real (kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: w0(0:)
    real (kind=dp_t), intent(in   ) :: time, dt, dx(:)
    real (kind=dp_t), intent(inout) :: T_max, enuc_max
    real (kind=dp_t), intent(inout) :: x_Tmax, y_Tmax, z_Tmax
    real (kind=dp_t), intent(inout) :: x_enucmax, y_enucmax, z_enucmax
    real (kind=dp_t), intent(inout) :: ratio_max
    real (kind=dp_t), intent(inout) :: x_ratiomax, y_ratiomax, z_ratiomax
    logical,          intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    !     Local variables
    integer                     :: i, j, k
    real (kind=dp_t) :: x, y , z
    real (kind=dp_t) :: enuc_local, ratio_local
    logical :: cell_valid

    do k = lo(3), hi(3)
       z = prob_lo(3) + (dble(k) + HALF) * dx(3)       
       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j) + HALF) * dx(2)
          do i = lo(1), hi(1)
             x = prob_lo(1) + (dble(i) + HALF) * dx(1)

             cell_valid = .true.
             if (present(mask)) then
                if ( (.not. mask(i,j,k)) ) cell_valid = .false.
             endif

             if (cell_valid .and. s(i,j,k,rho_comp) > base_cutoff_density) then

                ! temperature diagnostic
                if (s(i,j,k,temp_comp) > T_max) then
                   
                   T_max = s(i,j,k,temp_comp)
                   
                   x_Tmax = x
                   y_Tmax = y
                   z_Tmax = z
                   
                endif

                ! enuc diagnostic
                enuc_local = rho_Hnuc(i,j,k)/s(i,j,k,rho_comp)
                if (enuc_local > enuc_max) then
                
                   enuc_max = enuc_local

                   x_enucmax = x
                   y_enucmax = y
                   z_enucmax = z

                endif

                ! ratio cooling / nuc. energy generation
                ! cooling is where thermal < ZERO
                ratio_local = ZERO
                if (thermal(i,j,k) < ZERO) then
                   ratio_local = abs(thermal(i,j,k)) / rho_Hnuc(i,j,k)
                endif
             
                if (ratio_local > ratio_max) then

                   ratio_max = ratio_local

                   x_ratiomax = x
                   y_ratiomax = y
                   z_ratiomax = z

                endif

             endif

          enddo
       enddo
    enddo

  end subroutine diag_3d

  subroutine flush_diag()
    ! flush_diag is called before checkpointing.  If you want to buffer the 
    ! results of the diag routines and write out a lot of timesteps at once,
    ! flush_diag() is the routine that should do the writing

  end subroutine flush_diag

end module diag_module
