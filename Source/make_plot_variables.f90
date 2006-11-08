
module plot_variables_module

  use bl_types
  use multifab_module
  use eos_module
  use fill_3d_module
  use network
  use variables
  use geometry

  implicit none

contains

  subroutine make_enthalpy (enthalpy,comp,s)

    integer        , intent(in   ) :: comp
    type(multifab) , intent(inout) :: enthalpy
    type(multifab) , intent(in   ) :: s

    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: rp(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng,dm
    integer :: i

    ng = s%ng
    dm = s%dim

    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       sp => dataptr(s, i)
       rp => dataptr(enthalpy, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (2)
          call make_enthalpy_2d(rp(:,:,1,comp),sp(:,:,1,:), lo, hi, ng)
       case (3)
          call make_enthalpy_3d(rp(:,:,:,comp),sp(:,:,:,:), lo, hi, ng)
       end select
    end do

  end subroutine make_enthalpy

  subroutine make_enthalpy_2d (enthalpy,s,lo,hi,ng)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: enthalpy(lo(1):,lo(2):)  
    real (kind = dp_t), intent(in   ) ::    s(lo(1)-ng:,lo(2)-ng:,:)

!     Local variables
    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          enthalpy(i,j) = s(i,j,rhoh_comp) / s(i,j,rho_comp)
       enddo
    enddo

  end subroutine make_enthalpy_2d

  subroutine make_enthalpy_3d (enthalpy,s,lo,hi,ng)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: enthalpy(lo(1):,lo(2):, lo(3):)  
    real (kind = dp_t), intent(in   ) ::    s(lo(1)-ng:,lo(2)-ng:, lo(3)-ng:,:)
    
!     Local variables
    integer :: i, j, k
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             enthalpy(i,j,k) = s(i,j,k,rhoh_comp) / s(i,j,k,rho_comp)
          enddo
       enddo
    end do
    
  end subroutine make_enthalpy_3d

  subroutine make_tfromH (plotdata,comp_t,comp_dp,state,p0,temp0, dx)

    integer        , intent(in   ) :: comp_t,comp_dp
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: state
    real(kind=dp_t), intent(in   ) ::    p0(:)
    real(kind=dp_t), intent(in   ) :: temp0(:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: tp(:,:,:,:)
    integer :: lo(state%dim),hi(state%dim),ng,dm
    integer :: i

    ng = state%ng
    dm = state%dim

    do i = 1, state%nboxes
       if ( multifab_remote(state, i) ) cycle
       sp => dataptr(state, i)
       tp => dataptr(plotdata, i)
       lo =  lwb(get_box(state, i))
       hi =  upb(get_box(state, i))
       select case (dm)
       case (2)
          call maketfromH_2d(tp(:,:,1,comp_t),tp(:,:,1,comp_dp),sp(:,:,1,:), lo, hi, ng, p0, temp0)
       case (3)
          if (spherical .eq. 1) then
            call maketfromH_3d_sphr(tp(:,:,:,comp_t),tp(:,:,:,comp_dp),sp(:,:,:,:), &
                                    lo, hi, ng, p0, temp0, dx)
          else
            call maketfromH_3d_cart(tp(:,:,:,comp_t),tp(:,:,:,comp_dp),sp(:,:,:,:), &
                                    lo, hi, ng, p0, temp0)
          end if
       end select
    end do

  end subroutine make_tfromH

  subroutine maketfromH_2d (T,deltaP,state,lo,hi,ng,p0,temp0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) ::      T(lo(1)   :,lo(2):)  
    real (kind = dp_t), intent(  out) :: deltaP(lo(1)   :,lo(2):)  
    real (kind = dp_t), intent(in   ) ::  state(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::     p0(lo(2):)
    real (kind = dp_t), intent(in   ) ::  temp0(lo(2):)

    !     Local variables
    integer :: i, j

    integer :: nx

    do_diag = .false.

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          ! (rho, H) --> T, p
            
          den_row(1)  = state(i,j,rho_comp)
          h_row(1)    = state(i,j,rhoh_comp) / state(i,j,rho_comp)
          p_row(1)    = p0(j)
          temp_row(1) = temp0(j)
          xn_zone(:) = state(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

          input_flag = 2
  
          call eos(input_flag, den_row, temp_row, &
                   npts, nspec, &
                   xn_zone, aion, zion, &
                   p_row, h_row, e_row, &
                   cv_row, cp_row, xne_row, eta_row, pele_row, &
                   dpdt_row, dpdr_row, dedt_row, dedr_row, &
                   dpdX_row, dhdX_row, &
                   gam1_row, cs_row, s_row, &
                   do_diag)
          
          T(i,j) = log(temp_row(1))/log(10.)

          deltaP(i,j) = (p_row(1)-p0(j))/ p0(j)
          
       enddo
    enddo

  end subroutine maketfromH_2d

  subroutine maketfromH_3d_cart (T,deltaP,state,lo,hi,ng,p0,temp0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) ::      T(lo(1)   :,lo(2):   ,lo(3):     )  
    real (kind = dp_t), intent(  out) :: deltaP(lo(1)   :,lo(2):   ,lo(3):     )  
    real (kind = dp_t), intent(in   ) :: state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) ::    p0(lo(3):)
    real (kind = dp_t), intent(in   ) :: temp0(lo(3):)

    !     Local variables
    integer :: i, j, k

    integer :: nx

    do_diag = .false.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! (rho, H) --> T, p
            
             den_row(1)  = state(i,j,k,rho_comp)
             h_row(1)    = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)
             p_row(1)    = p0(k)
             temp_row(1) = temp0(k)
             xn_zone(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

             input_flag = 2
  
             call eos(input_flag, den_row, temp_row, &
                      npts, nspec, &
                      xn_zone, aion, zion, &
                      p_row, h_row, e_row, &
                      cv_row, cp_row, xne_row, eta_row, pele_row, &
                      dpdt_row, dpdr_row, dedt_row, dedr_row, &
                      dpdX_row, dhdX_row, &
                      gam1_row, cs_row, s_row, &
                      do_diag)

             T(i,j,k) = log(temp_row(1))/log(10.)

             deltaP(i,j,k) = (p_row(1)-p0(k))/ p0(k)

          enddo
       enddo
    enddo

  end subroutine maketfromH_3d_cart

  subroutine maketfromH_3d_sphr (T,deltaP,state,lo,hi,ng,p0,t0,dx)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) ::      T(lo(1)   :,lo(2):   ,lo(3):     )  
    real (kind = dp_t), intent(  out) :: deltaP(lo(1)   :,lo(2):   ,lo(3):     )  
    real (kind = dp_t), intent(in   ) :: state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) ::    p0(lo(3):)
    real (kind = dp_t), intent(in   ) ::    t0(lo(3):)
    real (kind = dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer :: i, j, k, nx
    real (kind=dp_t), allocatable :: t0_cart(:,:,:)
    real (kind=dp_t), allocatable :: p0_cart(:,:,:)

    allocate(t0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(t0_cart,t0,dx,0)

    allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(p0_cart,p0,dx,0)

    do_diag = .false.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            
             den_row(1)  = state(i,j,k,rho_comp)
             h_row(1)    = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)
             p_row(1)    = p0_cart(i,j,k)
             temp_row(1) = t0_cart(i,j,k)
             xn_zone(:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

             ! (rho, H) --> T, p
             input_flag = 2
  
             call eos(input_flag, den_row, temp_row, &
                      npts, nspec, &
                      xn_zone, aion, zion, &
                      p_row, h_row, e_row, &
                      cv_row, cp_row, xne_row, eta_row, pele_row, &
                      dpdt_row, dpdr_row, dedt_row, dedr_row, &
                      dpdX_row, dhdX_row, &
                      gam1_row, cs_row, s_row, &
                      do_diag)

             T(i,j,k) = log(temp_row(1))/log(10.)

             deltaP(i,j,k) = (p_row(1)-p0_cart(i,j,k))/ p0_cart(i,j,k)

          enddo
       enddo
    enddo
   
    deallocate(t0_cart,p0_cart)

  end subroutine maketfromH_3d_sphr

  subroutine make_tfromrho (plotdata,comp_tfromrho,comp_tpert,comp_rhopert, &
                            comp_machno,comp_deltag, &
                            s,u,s0,t0,p0,time,dx)

    integer        , intent(in   ) :: comp_tfromrho,comp_tpert
    integer        , intent(in   ) :: comp_rhopert, comp_machno, comp_deltag
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: s
    type(multifab) , intent(in   ) :: u
    real(kind=dp_t), intent(in   ) :: s0(:,:)
    real(kind=dp_t), intent(inout) :: t0(:)
    real(kind=dp_t), intent(in   ) :: p0(:)
    real(kind=dp_t), intent(in   ) :: time,dx(:)

    real(kind=dp_t), pointer:: sp(:,:,:,:),tp(:,:,:,:),up(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng,dm
    integer :: i

    ng = s%ng
    dm = s%dim

    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       tp => dataptr(plotdata, i)
       sp => dataptr(s, i)
       up => dataptr(u, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (2)
          call maketfromrho_2d(tp(:,:,1,comp_tfromrho),tp(:,:,1,comp_tpert), &
                               tp(:,:,1,comp_rhopert ), &
                               tp(:,:,1,comp_machno  ),tp(:,:,1,comp_deltag), &
                               sp(:,:,1,:), up(:,:,1,:), &
                               lo, hi, ng, s0, t0, p0, time, dx)
       case (3)
          if (spherical .eq. 1) then
            call maketfromrho_3d_sphr(tp(:,:,:,comp_tfromrho),tp(:,:,:,comp_tpert), &
                                      tp(:,:,:,comp_rhopert ), &
                                      tp(:,:,:,comp_machno  ),tp(:,:,:,comp_deltag), &
                                      sp(:,:,:,:), up(:,:,:,:), &
                                      lo, hi, ng, s0, t0, p0, time, dx)
          else
            call maketfromrho_3d_cart(tp(:,:,:,comp_tfromrho),tp(:,:,:,comp_tpert), &
                                      tp(:,:,:,comp_rhopert ), &
                                      tp(:,:,:,comp_machno  ),tp(:,:,:,comp_deltag), &
                                      sp(:,:,:,:), up(:,:,:,:), &
                                      lo, hi, ng, s0, t0, p0, time, dx)
          endif
       end select
    end do

  end subroutine make_tfromrho

  subroutine maketfromrho_2d (t,tpert,rhopert,machno,deltagamma,s,u,lo,hi,ng,s0,t0,p0,time,dx)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(  out)     ::      t(lo(1):,lo(2):)  
    real (kind=dp_t), intent(  out) ::      tpert(lo(1):,lo(2):)  
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1):,lo(2):)  
    real (kind=dp_t), intent(  out) ::     machno(lo(1):,lo(2):)  
    real (kind=dp_t), intent(  out) :: deltagamma(lo(1):,lo(2):)  
    real (kind=dp_t), intent(in   ) ::  s(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(in   ) ::  u(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(in   ) :: s0(lo(2):,:)
    real (kind=dp_t), intent(inout) :: t0(lo(2):)
    real (kind=dp_t), intent(in   ) :: p0(lo(2):)
    real (kind=dp_t), intent(in   ) :: time,dx(:)

    !     Local variables
    integer          :: i, j
    real (kind=dp_t) :: vel
    real (kind=dp_t), allocatable :: gam10(:)

    allocate(gam10(lo(2):hi(2)))

    do_diag = .false.

    ! First make sure we have t0 right.
    do j = lo(2), hi(2)
        den_row(1) = s0(j,rho_comp)
       temp_row(1) = t0(j)
          p_row(1) = p0(j)
       xn_zone(:) = s0(j,spec_comp:spec_comp+nspec-1)/den_row(1)

       ! (rho,P) --> T,h
       input_flag = 4

       call eos(input_flag, den_row, temp_row, &
                npts, nspec, &
                xn_zone, aion, zion, &
                p_row, h_row, e_row, & 
                cv_row, cp_row, xne_row, eta_row, pele_row, &
                dpdt_row, dpdr_row, dedt_row, dedr_row, &
                dpdX_row, dhdX_row, &
                gam1_row, cs_row, s_row, &
                do_diag)
        t0(j) = temp_row(1)
     gam10(j) = gam1_row(1)
    end do

    ! Then compute the perturbation
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          den_row(1) = s(i,j,rho_comp)
          temp_row(1) = t0(j)
          p_row(1) = p0(j)
          xn_zone(:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_row(1)

          ! (rho,P) --> T,h
          input_flag = 4

          call eos(input_flag, den_row, temp_row, &
                   npts, nspec, &
                   xn_zone, aion, zion, &
                   p_row, h_row, e_row, & 
                   cv_row, cp_row, xne_row, eta_row, pele_row, &
                   dpdt_row, dpdr_row, dedt_row, dedr_row, &
                   dpdX_row, dhdX_row, &
                   gam1_row, cs_row, s_row, &
                   do_diag)
          
          t(i,j) = log(temp_row(1))/log(10.)
          tpert(i,j) = t(i,j) - temp_row(1)

          rhopert(i,j) = s(i,j,rho_comp) - s0(j,rho_comp)

          vel = sqrt(u(i,j,1)*u(i,j,1) + u(i,j,2)*u(i,j,2))
          machno(i,j) = vel / cs_row(1)

          deltagamma(i,j) = gam1_row(1) - gam10(j)
       enddo
    enddo

    deallocate(gam10)

  end subroutine maketfromrho_2d

  subroutine maketfromrho_3d_cart (t,tpert,rhopert,machno,deltagamma, &
                                   s,u,lo,hi,ng,s0,t0,p0,time,dx)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(  out) ::          t(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::      tpert(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::     machno(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) :: deltagamma(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(in   ) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(in   ) ::  u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(in   ) :: s0(lo(3):,:)
    real (kind=dp_t), intent(inout) :: t0(lo(3):)
    real (kind=dp_t), intent(in   ) :: p0(lo(3):)
    real (kind=dp_t), intent(in   ) :: time,dx(:)

    !     Local variables
    integer          :: i, j, k
    real (kind=dp_t) :: vel
    real (kind=dp_t), allocatable :: gam10(:)

    allocate(gam10(lo(3):hi(3)))

    do_diag = .false.

    ! First make sure we have t0 right.
    do k = lo(3), hi(3)
        den_row(1) = s0(k,rho_comp)
       temp_row(1) = t0(k)
          p_row(1) = p0(k)
       xn_zone(:) = s0(k,spec_comp:spec_comp+nspec-1)/den_row(1)

       ! (rho,P) --> T,h
       input_flag = 4

       call eos(input_flag, den_row, temp_row, &
                npts, nspec, &
                xn_zone, aion, zion, &
                p_row, h_row, e_row, & 
                cv_row, cp_row, xne_row, eta_row, pele_row, &
                dpdt_row, dpdr_row, dedt_row, dedr_row, &
                dpdX_row, dhdX_row, &
                gam1_row, cs_row, s_row, &
                do_diag)
        t0(k) = temp_row(1)
     gam10(k) = gam1_row(1)
    end do

    ! Then compute the perturbation and Mach number
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_row(1) = s(i,j,k,rho_comp)
             temp_row(1) = t0(k)
             p_row(1) = p0(k)
             xn_zone(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

             ! (rho,P) --> T,h
             input_flag = 4

             call eos(input_flag, den_row, temp_row, &
                      npts, nspec, &
                      xn_zone, aion, zion, &
                      p_row, h_row, e_row, & 
                      cv_row, cp_row, xne_row, eta_row, pele_row, &
                      dpdt_row, dpdr_row, dedt_row, dedr_row, &
                      dpdX_row, dhdX_row, &
                      gam1_row, cs_row, s_row, &
                      do_diag)

             t(i,j,k) = log(temp_row(1))/log(10.)
             tpert(i,j,k) = temp_row(1) - t0(k)

             rhopert(i,j,k) = s(i,j,k,rho_comp) - s0(k,rho_comp)

             vel = sqrt(u(i,j,k,1)*u(i,j,k,1) + u(i,j,k,2)*u(i,j,k,2) + u(i,j,k,3)*u(i,j,k,3))
             machno(i,j,k) = vel / cs_row(1)

             deltagamma(i,j,k) = gam1_row(1) - gam10(k)
          enddo
       enddo
    enddo

    deallocate(gam10)

   end subroutine maketfromrho_3d_cart

  subroutine maketfromrho_3d_sphr (t,tpert,rhopert,machno,deltagamma, &
                                   s,u,lo,hi,ng,s0,t0,p0,time,dx)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(  out) ::          t(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::      tpert(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::     machno(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) :: deltagamma(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(in   ) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(in   ) ::  u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(in   ) :: s0(lo(3):,:)
    real (kind=dp_t), intent(inout) :: t0(lo(3):)
    real (kind=dp_t), intent(in   ) :: p0(lo(3):)
    real (kind=dp_t), intent(in   ) :: time,dx(:)

    !     Local variables
    integer          :: i, j, k
    real (kind=dp_t) :: vel
    real (kind=dp_t), allocatable :: gam10(:)
    real (kind=dp_t), allocatable :: rho0_cart(:,:,:)
    real (kind=dp_t), allocatable ::   t0_cart(:,:,:)
    real (kind=dp_t), allocatable ::   p0_cart(:,:,:)
    real (kind=dp_t), allocatable :: gam0_cart(:,:,:)

    allocate(gam10(lo(3):hi(3)))

    do_diag = .false.

    ! First make sure we have t0 right.
    do k = lo(3), hi(3)
        den_row(1) = s0(k,rho_comp)
       temp_row(1) = t0(k)
          p_row(1) = p0(k)
       xn_zone(:) = s0(k,spec_comp:spec_comp+nspec-1)/den_row(1)

       ! (rho,P) --> T,h
       input_flag = 4

       call eos(input_flag, den_row, temp_row, &
                npts, nspec, &
                xn_zone, aion, zion, &
                p_row, h_row, e_row, & 
                cv_row, cp_row, xne_row, eta_row, pele_row, &
                dpdt_row, dpdr_row, dedt_row, dedr_row, &
                dpdX_row, dhdX_row, &
                gam1_row, cs_row, s_row, &
                do_diag)
        t0(k) = temp_row(1)
     gam10(k) = gam1_row(1)
    end do

    allocate(rho0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(rho0_cart,s0(:,rho_comp),dx,0)

    allocate(t0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(t0_cart,t0,dx,0)

    allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(p0_cart,p0,dx,0)

    allocate(gam0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(gam0_cart,gam10,dx,0)

    ! Then compute the perturbation and Mach number
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_row(1) = s(i,j,k,rho_comp)
             temp_row(1) = t0_cart(i,j,k)
             p_row(1) = p0_cart(i,j,k)
             xn_zone(:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_row(1)

             ! (rho,P) --> T,h
             input_flag = 4

             call eos(input_flag, den_row, temp_row, &
                      npts, nspec, &
                      xn_zone, aion, zion, &
                      p_row, h_row, e_row, & 
                      cv_row, cp_row, xne_row, eta_row, pele_row, &
                      dpdt_row, dpdr_row, dedt_row, dedr_row, &
                      dpdX_row, dhdX_row, &
                      gam1_row, cs_row, s_row, &
                      do_diag)

             t(i,j,k) = log(temp_row(1))/log(10.)
             tpert(i,j,k) = temp_row(1) - t0_cart(i,j,k)

             rhopert(i,j,k) = s(i,j,k,rho_comp) - rho0_cart(i,j,k)

             vel = sqrt(u(i,j,k,1)*u(i,j,k,1) + u(i,j,k,2)*u(i,j,k,2) + u(i,j,k,3)*u(i,j,k,3))
             machno(i,j,k) = vel / cs_row(1)

             deltagamma(i,j,k) = gam1_row(1) - gam0_cart(i,j,k)
          enddo
       enddo
    enddo

    deallocate(gam10)
    deallocate(rho0_cart,t0_cart,p0_cart,gam0_cart)

   end subroutine maketfromrho_3d_sphr

   subroutine make_XfromrhoX (plotdata,comp,s)

      integer        , intent(in   ) :: comp
      type(multifab) , intent(in   ) :: s
      type(multifab) , intent(inout) :: plotdata

      real(kind=dp_t), pointer:: sp(:,:,:,:)
      real(kind=dp_t), pointer:: pp(:,:,:,:)
      integer :: lo(s%dim),hi(s%dim),ng,dm
      integer :: i

      ng = s%ng
      dm = s%dim

      do i = 1, s%nboxes
         if ( multifab_remote(s, i) ) cycle
         sp => dataptr(s, i)
         pp => dataptr(plotdata, i)
         lo =  lwb(get_box(s, i))
         hi =  upb(get_box(s, i))
         select case (dm)
            case (2)
              call makeXfromrhoX_2d(pp(:,:,1,comp:),sp(:,:,1,1),sp(:,:,1,spec_comp:), lo, hi, ng)
            case (3)
              call makeXfromrhoX_3d(pp(:,:,:,comp:),sp(:,:,:,1),sp(:,:,:,spec_comp:), lo, hi, ng)
         end select
      end do

   end subroutine make_XfromrhoX

   subroutine makeXfromrhoX_2d (X,rho,rhoX,lo,hi,ng)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) ::    X(lo(1)   :,lo(2)   :,:)  
      real (kind = dp_t), intent(in   ) :: rhoX(lo(1)-ng:,lo(2)-ng:,:)
      real (kind = dp_t), intent(in   ) ::  rho(lo(1)-ng:,lo(2)-ng:  )

!     Local variables
      integer :: i, j, n

      do n = 1, nspec
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
         X(i,j,n) = rhoX(i,j,n) / rho(i,j)
      enddo
      enddo
      enddo

   end subroutine makeXfromrhoX_2d

   subroutine makeXfromrhoX_3d (X,rho,rhoX,lo,hi,ng)

      implicit none

      integer, intent(in) :: lo(:), hi(:), ng
      real (kind = dp_t), intent(  out) ::    X(lo(1)   :,lo(2)   :,lo(3)   :,:)  
      real (kind = dp_t), intent(in   ) :: rhoX(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
      real (kind = dp_t), intent(in   ) ::  rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:  )

!     Local variables
      integer :: i, j, k, n

      do n = 1, nspec
      do k = lo(3), hi(3)
      do j = lo(2), hi(2)
      do i = lo(1), hi(1)
         X(i,j,k,n) = rhoX(i,j,k,n) / rho(i,j,k)
      enddo
      enddo
      enddo
      enddo

   end subroutine makeXfromrhoX_3d

end module plot_variables_module
