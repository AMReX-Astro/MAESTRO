! a simply sanity checker to ensure that MAESTRO is still operating in
! a valid regime.  This is based, initially, off of diag.f90

module sanity_module

  use bl_types
  use bl_IO_module
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none

  private

  public :: sanity_check

contains

  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  subroutine sanity_check(newtime,dx,s,rho0,rhoh0,p0,tempbar, &
                          gamma1bar,div_coeff, &
                          u,w0,normal, &
                          mla,the_bc_tower)

    use bl_prof_module
    use geometry, only: spherical
    use bl_constants_module
    use network, only: network_species_index
    use inlet_bc_module
    use fill_3d_module, only: put_1d_array_on_cart, make_w0mac
    use variables, only: foextrap_comp
    use probin_module, only: mach_max_abort

    real(kind=dp_t), intent(in   ) :: dx(:,:),newtime
    type(multifab) , intent(in   ) :: s(:)
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

    integer :: lo(mla%dim),hi(mla%dim),dm,nlevs
    integer :: ng_s,ng_u,ng_n,ng_rhn,ng_rhe,ng_rw,ng_w,ng_wm
    integer :: i,n,comp
    logical :: lexist

    logical, save :: firstCall_io = .true.

    type(bl_prof_timer), save :: bpt

    call build(bpt, "sanity check")

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
    if (spherical == 1) then
       ng_n   = nghost(normal(1))
       ng_w   = nghost(w0r_cart(1))
       ng_wm  = nghost(w0mac(1,1))
    endif


    !=========================================================================
    ! initialize
    !=========================================================================
    Mach_max = ZERO


    !=========================================================================
    ! loop over the levels and compute the global quantities
    !=========================================================================
    do n = 1, nlevs

       ! initialize the local (processor's version) and level quantities to 0
       Mach_max_level = ZERO
       Mach_max_local = ZERO


       !----------------------------------------------------------------------
       ! loop over boxes in a given level
       !----------------------------------------------------------------------
       do i = 1, s(n)%nboxes
          if ( multifab_remote(s(n), i) ) cycle

          sp => dataptr(s(n) , i)
          up => dataptr(u(n) , i)

          lo =  lwb(get_box(s(n), i))
          hi =  upb(get_box(s(n), i))

          select case (dm)

          case (1)
             if (n .eq. nlevs) then
                call sanity_1d(n,newtime,dx(n,:), &
                               sp(:,1,1,:),ng_s, &
                               rho0(n,:),rhoh0(n,:), &
                               p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                               up(:,1,1,:),ng_u, &
                               w0(n,:), &
                               lo,hi, &
                               Mach_max_local)
             else
                mp => dataptr(mla%mask(n), i)
                call sanity_1d(n,newtime,dx(n,:), &
                               sp(:,1,1,:),ng_s, &
                               rho0(n,:),rhoh0(n,:), &
                               p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                               up(:,1,1,:),ng_u, &
                               w0(n,:), &
                               lo,hi, &
                               Mach_max_local, &
                               mp(:,1,1,1))
             endif

          case (2)
             if (n .eq. nlevs) then
                call sanity_2d(n,newtime,dx(n,:), &
                               sp(:,:,1,:),ng_s, &
                               rho0(n,:),rhoh0(n,:), &
                               p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                               up(:,:,1,:),ng_u, &
                               w0(n,:), &
                               lo,hi, &
                               Mach_max_local)
             else
                mp => dataptr(mla%mask(n), i)
                call sanity_2d(n,newtime,dx(n,:), &
                               sp(:,:,1,:),ng_s, &
                               rho0(n,:),rhoh0(n,:), &
                               p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                               up(:,:,1,:),ng_u, &
                               w0(n,:), &
                               lo,hi, &
                               Mach_max_local, &
                               mp(:,:,1,1))
             endif

          case (3)
             if (spherical == 1) then

                nop => dataptr(normal(n) , i)
                w0rp => dataptr(w0r_cart(n), i)
                w0xp => dataptr(w0mac(n,1), i)
                w0yp => dataptr(w0mac(n,2), i)
                w0zp => dataptr(w0mac(n,3), i)

                if (n .eq. nlevs) then
                   call sanity_3d_sph(n,newtime,dx(n,:), &
                                      sp(:,:,:,:),ng_s, &
                                      rho0(n,:),rhoh0(n,:), &
                                      p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                                      up(:,:,:,:),ng_u, &
                                      w0rp(:,:,:,1), ng_w, &
                                      w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1),ng_wm, &
                                      nop(:,:,:,:),ng_n, &
                                      lo,hi, &
                                      Mach_max_local)
                else
                   mp => dataptr(mla%mask(n), i)
                   call sanity_3d_sph(n,newtime,dx(n,:), &
                                      sp(:,:,:,:),ng_s, &
                                      rho0(n,:),rhoh0(n,:), &
                                      p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                                      up(:,:,:,:),ng_u, &
                                      w0rp(:,:,:,1), ng_w, &
                                      w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1),ng_wm, &
                                      nop(:,:,:,:),ng_n, &
                                      lo,hi, &
                                      Mach_max_local, &
                                      mp(:,:,:,1))
                endif

             else
                if (n .eq. nlevs) then
                   call sanity_3d(n,newtime,dx(n,:), &
                                  sp(:,:,:,:),ng_s, &
                                  rho0(n,:),rhoh0(n,:), &
                                  p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                                  up(:,:,:,:),ng_u, &
                                  w0(n,:), &
                                  lo,hi, &
                                  Mach_max_local)
                else
                   mp => dataptr(mla%mask(n), i)
                   call sanity_3d(n,newtime,dx(n,:), &
                                  sp(:,:,:,:),ng_s, &
                                  rho0(n,:),rhoh0(n,:), &
                                  p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                                  up(:,:,:,:),ng_u, &
                                  w0(n,:), &
                                  lo,hi, &
                                  Mach_max_local, &
                                  mp(:,:,:,1))
                endif

             endif
          end select
       end do

       !----------------------------------------------------------------------
       ! do the appropriate parallel reduction for the current level
       !----------------------------------------------------------------------

       ! NOTE: only the I/O Processor will have the correct reduced value
       call parallel_reduce(Mach_max_level, Mach_max_local, MPI_MAX, &
                            proc = parallel_IOProcessorNode())


       !----------------------------------------------------------------------
       ! reduce the current level's data with the global data
       !----------------------------------------------------------------------
       if (parallel_IOProcessor()) then
          Mach_max = max(Mach_max, Mach_max_level)
       endif

    end do


    !=========================================================================
    ! normalize 
    !=========================================================================

    ! normalize any integral quantities here


    !=========================================================================
    ! perform sanity tests
    !=========================================================================
    if (Mach_max >= mach_max_abort) then
       call bl_error("ERROR: Mach number too high")
    end if


    !=========================================================================
    ! clean-up
    !=========================================================================
    call destroy(bpt)

    if (spherical == 1) then
       do n = 1, nlevs
          call destroy(w0r_cart(n))
          do comp=1, dm
             call destroy(w0mac(n,comp))
          enddo
       enddo
    endif

  end subroutine sanity_check



  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  subroutine sanity_1d(n,newtime,dx, &
                       s,ng_s, &
                       rho0,rhoh0,p0,tempbar,gamma1bar, &
                       u,ng_u, &
                       w0, &
                       lo,hi, &
                       Mach_max, &
                       mask)

    use variables, only: rho_comp, spec_comp, temp_comp
    use bl_constants_module
    use network, only: nspec
    use probin_module, only: prob_lo
    use eos_module

    integer, intent(in) :: n, lo(:), hi(:), ng_s, ng_u
    real (kind=dp_t), intent(in   ) ::         s(lo(1)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:), rhoh0(0:), &
                                         p0(0:),tempbar(0:),gamma1bar(0:)
    real (kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: w0(0:)
    real (kind=dp_t), intent(in   ) :: newtime, dx(:)
    real (kind=dp_t), intent(inout) :: Mach_max
    logical,          intent(in   ), optional :: mask(lo(1):)

    !     Local variables
    integer            :: i
    logical            :: cell_valid
    real (kind=dp_t)   :: vel, x


    do i = lo(1), hi(1)
       x = prob_lo(1) + (dble(i)+HALF) * dx(1)             

       cell_valid = .true.
       if (present(mask)) then
          if ( (.not. mask(i)) ) cell_valid = .false.
       endif
       
       if (cell_valid) then

          ! vel is the magnitude of the velocity, including w0
          vel = sqrt( (u(i,1) + HALF*(w0(i) + w0(i+1)) )**2 )

             
          ! call the EOS to get the sound speed and internal energy       
          temp_eos(1) = s(i,temp_comp)
          den_eos(1)  = s(i,rho_comp)
          xn_eos(1,:) = s(i,spec_comp:spec_comp+nspec-1)/den_eos(1)

          call eos(x, eos_input_rt, den_eos, temp_eos, &
                   npts, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, &
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   .false.)


          ! max Mach number                                       
          Mach_max = max(Mach_max,vel/cs_eos(1))

       endif  ! cell valid

    enddo

  end subroutine sanity_1d
  

  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  subroutine sanity_2d(n,newtime,dx, &
                       s,ng_s, &
                       rho0,rhoh0,p0,tempbar,gamma1bar, &
                       u,ng_u, &
                       w0, &
                       lo,hi, &
                       Mach_max, &
                       mask)

    use variables, only: rho_comp, spec_comp, temp_comp
    use bl_constants_module
    use network, only: nspec
    use probin_module, only: prob_lo
    use eos_module

    integer, intent(in) :: n, lo(:), hi(:), ng_s, ng_u
    real (kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:), rhoh0(0:), &
                                         p0(0:),tempbar(0:),gamma1bar(0:)
    real (kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: w0(0:)
    real (kind=dp_t), intent(in   ) :: newtime, dx(:)
    real (kind=dp_t), intent(inout) :: Mach_max
    logical,          intent(in   ), optional :: mask(lo(1):,lo(2):)

    !     Local variables
    integer            :: i, j
    logical            :: cell_valid
    real (kind=dp_t)   :: vel, y


    do j = lo(2), hi(2)
       y = prob_lo(2) + (dble(j)+HALF) * dx(2)

       do i = lo(1), hi(1)

          cell_valid = .true.
          if (present(mask)) then
             if ( (.not. mask(i,j)) ) cell_valid = .false.
          endif

          if (cell_valid) then

             ! vel is the magnitude of the velocity, including w0
             vel = sqrt(  u(i,j,1)**2 + &
                        ( u(i,j,2) + HALF*(w0(j) + w0(j+1)) )**2 )

             
             ! call the EOS to get the sound speed and internal energy       
             temp_eos(1) = s(i,j,temp_comp)
             den_eos(1)  = s(i,j,rho_comp)
             xn_eos(1,:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)

             pt_index_eos(:) = (/i, j, -1/)       

             call eos(y, eos_input_rt, den_eos, temp_eos, &
                      npts, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      .false., pt_index_eos)


             ! max Mach number                                       
             Mach_max = max(Mach_max,vel/cs_eos(1))

          endif  ! cell valid

       enddo
    enddo

  end subroutine sanity_2d


  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  subroutine sanity_3d(n,newtime,dx, &
                       s,ng_s, &
                       rho0,rhoh0,p0,tempbar,gamma1bar, &
                       u,ng_u, &
                       w0, &
                       lo,hi, &
                       Mach_max,&
                       mask)

    use variables, only: rho_comp, spec_comp, temp_comp
    use bl_constants_module
    use network, only: nspec
    use probin_module, only: prob_lo
    use eos_module

    integer, intent(in) :: n, lo(:), hi(:), ng_s, ng_u
    real (kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:), rhoh0(0:), &
                                         p0(0:),tempbar(0:),gamma1bar(0:)
    real (kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: w0(0:)
    real (kind=dp_t), intent(in   ) :: newtime, dx(:)
    real (kind=dp_t), intent(inout) :: Mach_max
    logical,          intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    !     Local variables
    integer            :: i, j, k
    logical            :: cell_valid
    real (kind=dp_t)   :: vel, z


    do k = lo(3), hi(3)
       z = prob_lo(3) + (dble(k)+HALF) * dx(3)             

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             cell_valid = .true.
             if (present(mask)) then
                if ( (.not. mask(i,j,k)) ) cell_valid = .false.
             endif

             if (cell_valid) then

                ! vel is the magnitude of the velocity, including w0
                vel = sqrt(  u(i,j,k,1)**2 + &
                             u(i,j,k,2)**2 + &
                           ( u(i,j,k,3) + HALF*(w0(k) + w0(k+1)) )**2 )

             
                ! call the EOS to get the sound speed and internal energy       
                temp_eos(1) = s(i,j,k,temp_comp)
                den_eos(1)  = s(i,j,k,rho_comp)
                xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

                call eos(z, eos_input_rt, den_eos, temp_eos, &
                         npts, &
                         xn_eos, &
                         p_eos, h_eos, e_eos, &
                         cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                         dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                         dpdX_eos, dhdX_eos, &
                         gam1_eos, cs_eos, s_eos, &
                         dsdt_eos, dsdr_eos, &
                         .false.)


                ! max Mach number                                       
                Mach_max = max(Mach_max,vel/cs_eos(1))

             endif  ! cell valid

          enddo
       enddo
    enddo

  end subroutine sanity_3d


  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  subroutine sanity_3d_sph(n,newtime,dx, &
                           s,ng_s, &
                           rho0,rhoh0,p0,tempbar,gamma1bar, &
                           u,ng_u, &
                           w0r,ng_w, &
                           w0macx,w0macy,w0macz,ng_wm, &
                           normal,ng_n, &
                           lo,hi, &                         
                           Mach_max, &
                           mask)

    use variables, only: rho_comp, spec_comp, temp_comp
    use bl_constants_module
    use network, only: nspec
    use probin_module, only: prob_lo
    use eos_module
    use geometry, only: center

    integer, intent(in) :: n, lo(:), hi(:), ng_s, ng_u, ng_w, ng_wm, ng_n
    real (kind=dp_t), intent(in   ) ::      s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:), rhoh0(0:), &
                                         p0(0:),tempbar(0:),gamma1bar(0:)
    real (kind=dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
    real (kind=dp_t), intent(in   ) ::      w0r(lo(1)-ng_w:  ,lo(2)-ng_w:  ,lo(3)-ng_w:)
    real (kind=dp_t), intent(in   ) ::   w0macx(lo(1)-ng_wm: ,lo(2)-ng_wm: ,lo(3)-ng_wm:)
    real (kind=dp_t), intent(in   ) ::   w0macy(lo(1)-ng_wm: ,lo(2)-ng_wm: ,lo(3)-ng_wm:)
    real (kind=dp_t), intent(in   ) ::   w0macz(lo(1)-ng_wm: ,lo(2)-ng_wm: ,lo(3)-ng_wm:)
    real (kind=dp_t), intent(in   ) ::   normal(lo(1)-ng_n:  ,lo(2)-ng_n:  ,lo(3)-ng_n:,:)
    real (kind=dp_t), intent(in   ) :: newtime, dx(:)
    real (kind=dp_t), intent(inout) :: Mach_max
    logical,          intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    !     Local variables
    integer            :: i, j, k
    logical            :: cell_valid
    real (kind=dp_t)   :: vel,x,y,z,r


    do k = lo(3), hi(3)
       z = prob_lo(3) + (dble(k)+HALF) * dx(3)             
       
       do j = lo(2), hi(2)
          y = prob_lo(2) + (dble(j)+HALF) * dx(2)             
       
          do i = lo(1), hi(1)
             x = prob_lo(1) + (dble(i)+HALF) * dx(1)             

             cell_valid = .true.
             if (present(mask)) then
                if ( (.not. mask(i,j,k)) ) cell_valid = .false.
             endif

             if (cell_valid) then

                ! vel is the magnitude of the velocity, including w0
                vel = sqrt( (u(i,j,k,1)+HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)))**2 + &
                            (u(i,j,k,2)+HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)))**2 + &
                            (u(i,j,k,3)+HALF*(w0macz(i,j,k)+w0macz(i,j,k+1)))**2)

             
                ! call the EOS to get the sound speed and internal energy       
                temp_eos(1) = s(i,j,k,temp_comp)
                den_eos(1)  = s(i,j,k,rho_comp)
                xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

                r = sqrt( (x -center(1))**2 + (y -center(2))**2 + &
                     (z -center(3))**2)

                call eos(r, eos_input_rt, den_eos, temp_eos, &
                         npts, &
                         xn_eos, &
                         p_eos, h_eos, e_eos, &
                         cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                         dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                         dpdX_eos, dhdX_eos, &
                         gam1_eos, cs_eos, s_eos, &
                         dsdt_eos, dsdr_eos, &
                         .false.)


                ! max Mach number                                       
                Mach_max = max(Mach_max,vel/cs_eos(1))

             endif  ! cell valid

          enddo
       enddo
    enddo

  end subroutine sanity_3d_sph

end module sanity_module
