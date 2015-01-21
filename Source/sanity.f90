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

    integer :: lo(3),hi(3),dm,nlevs
    integer :: ng_n,ng_w,ng_wm
    integer :: i,n,comp

    type(bl_prof_timer), save :: bpt

    call build(bpt, "sanity check")

    dm = mla%dim
    nlevs = mla%nlevel

    print *, "*** in sanity ***"
    
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
       lo(:) = 1; hi(:) = 1

       do i = 1, nfabs(s(n))

          sp => dataptr(s(n) , i)
          up => dataptr(u(n) , i)
          
          lo(1:dm) =  lwb(get_box(s(n), i))
          hi(1:dm) =  upb(get_box(s(n), i))

          if (spherical == 1) then
             
             nop => dataptr(normal(n) , i)
             w0rp => dataptr(w0r_cart(n), i)
             w0xp => dataptr(w0mac(n,1), i)
             w0yp => dataptr(w0mac(n,2), i)
             w0zp => dataptr(w0mac(n,3), i)

             if (n .eq. nlevs) then
                call sanity_3d_sph(n,newtime,dx(n,:), &
                                   sp, lbound(sp),ubound(sp), &
                                   rho0(n,:),rhoh0(n,:), &
                                   p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                                   up, lbound(up),ubound(up), &
                                   w0rp(:,:,:,1), ng_w, &
                                   w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1),ng_wm, &
                                   nop(:,:,:,:),ng_n, &
                                   lo,hi, &
                                   Mach_max_local)
             else
                mp => dataptr(mla%mask(n), i)
                call sanity_3d_sph(n,newtime,dx(n,:), &
                                   sp, lbound(sp),ubound(sp), &
                                   rho0(n,:),rhoh0(n,:), &
                                   p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                                   up, lbound(up),ubound(up), &
                                   w0rp(:,:,:,1), ng_w, &
                                   w0xp(:,:,:,1),w0yp(:,:,:,1),w0zp(:,:,:,1),ng_wm, &
                                   nop(:,:,:,:),ng_n, &
                                   lo,hi, &
                                   Mach_max_local, &
                                   mp(:,:,:,1))
             endif

          else
             if (n .eq. nlevs) then
                call sanity_cart(dm, n,newtime,dx(n,:), &
                                 sp, lbound(sp), ubound(sp), &
                                 rho0(n,:),rhoh0(n,:), &
                                 p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                                 up, lbound(up), ubound(up), &
                                 w0(n,:), &
                                 lo,hi, &
                                 Mach_max_local)
             else
                mp => dataptr(mla%mask(n), i)
                call sanity_cart(dm, n,newtime,dx(n,:), &
                                 sp, lbound(sp), ubound(sp), &
                                 rho0(n,:),rhoh0(n,:), &
                                 p0(n,:),tempbar(n,:),gamma1bar(n,:), &
                                 up, lbound(up), ubound(up), &
                                 w0(n,:), &
                                 lo,hi, &
                                 Mach_max_local, &
                                 mp(:,:,:,1))
             endif
                
          endif
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
  subroutine sanity_cart(dm, n,newtime,dx, &
                         s, slo, shi, &
                         rho0,rhoh0,p0,tempbar,gamma1bar, &
                         u, ulo, uhi, &
                         w0, &
                         lo,hi, &
                         Mach_max,&
                         mask)

    use variables, only: rho_comp, spec_comp, temp_comp
    use bl_constants_module
    use network, only: nspec
    use eos_module, only: eos_input_rt, eos
    use eos_type_module
    
    integer, intent(in) :: dm, n, lo(:), hi(:)
    integer, intent(in) :: slo(4), shi(4), ulo(4), uhi(4)
    real (kind=dp_t), intent(in   ) :: s(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),slo(4):shi(4))
    real (kind=dp_t), intent(in   ) :: rho0(0:), rhoh0(0:), &
                                         p0(0:),tempbar(0:),gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),ulo(4):uhi(4))
    real (kind=dp_t), intent(in   ) :: w0(0:)
    real (kind=dp_t), intent(in   ) :: newtime, dx(:)
    real (kind=dp_t), intent(inout) :: Mach_max
    logical,          intent(in   ), optional :: mask(lo(1):,lo(2):,lo(3):)

    !     Local variables
    integer            :: i, j, k
    logical            :: cell_valid
    real (kind=dp_t)   :: vel

    type (eos_t) :: eos_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             cell_valid = .true.
             if (present(mask)) then
                if ( (.not. mask(i,j,k)) ) cell_valid = .false.
             endif

             if (cell_valid) then

                ! vel is the magnitude of the velocity, including w0
                select case (dm)
                case (1)
                   vel =  u(i,j,k,1) + HALF*(w0(i) + w0(i+1))

                case (2)
                   vel = sqrt(  u(i,j,k,1)**2 + &
                              ( u(i,j,k,2) + HALF*(w0(j) + w0(j+1)) )**2 )

                case (3)
                   vel = sqrt(  u(i,j,k,1)**2 + &
                                u(i,j,k,2)**2 + &
                              ( u(i,j,k,3) + HALF*(w0(k) + w0(k+1)) )**2 )                   
                end select
             
                ! call the EOS to get the sound speed and internal energy       
                eos_state%T     = s(i,j,k,temp_comp)
                eos_state%rho   = s(i,j,k,rho_comp)
                eos_state%xn(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                call eos(eos_input_rt, eos_state, .false.)


                ! max Mach number                                       
                Mach_max = max(Mach_max,vel/eos_state%cs)

             endif  ! cell valid

          enddo
       enddo
    enddo

  end subroutine sanity_cart


  !xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  subroutine sanity_3d_sph(n,newtime,dx, &
                           s,slo,shi, &
                           rho0,rhoh0,p0,tempbar,gamma1bar, &
                           u,ulo,uhi, &
                           w0r,ng_w, &
                           w0macx,w0macy,w0macz,ng_wm, &
                           normal,ng_n, &
                           lo,hi, &                         
                           Mach_max, &
                           mask)

    use variables, only: rho_comp, spec_comp, temp_comp
    use bl_constants_module
    use network, only: nspec
    use eos_module, only: eos_input_rt, eos
    use eos_type_module

    integer, intent(in) :: n, lo(:), hi(:), ng_w, ng_wm, ng_n
    integer, intent(in) :: slo(4), shi(4), ulo(4), uhi(4)
    real (kind=dp_t), intent(in   ) :: s(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),slo(4):shi(4))    
    real (kind=dp_t), intent(in   ) :: rho0(0:), rhoh0(0:), &
                                         p0(0:),tempbar(0:),gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: u(ulo(1):uhi(1),ulo(2):uhi(2),ulo(3):uhi(3),ulo(4):uhi(4))    
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
    real (kind=dp_t)   :: vel
 
    type (eos_t) :: eos_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

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
                eos_state%T     = s(i,j,k,temp_comp)
                eos_state%rho   = s(i,j,k,rho_comp)
                eos_state%xn(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/eos_state%rho

                call eos(eos_input_rt, eos_state, .false.)


                ! max Mach number                                       
                Mach_max = max(Mach_max,vel/eos_state%cs)

             endif  ! cell valid

          enddo
       enddo
    enddo

  end subroutine sanity_3d_sph

end module sanity_module
