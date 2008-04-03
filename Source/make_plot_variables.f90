module plot_variables_module

  use bl_types
  use multifab_module

  implicit none

  private

  public :: make_enthalpy, make_tfromH, make_tfromrho, make_XfromrhoX
  public :: make_omegadot, make_deltaT

contains

  subroutine make_enthalpy(enthalpy,comp,s)

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

  subroutine make_enthalpy_2d(enthalpy,s,lo,hi,ng)

    use variables, only: rho_comp, rhoh_comp

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: enthalpy(lo(1):,lo(2):)  
    real (kind = dp_t), intent(in   ) ::    s(lo(1)-ng:,lo(2)-ng:,:)

    ! Local variables
    integer :: i, j, comp

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          enthalpy(i,j) = s(i,j,rhoh_comp)/s(i,j,rho_comp)
       enddo
    enddo

  end subroutine make_enthalpy_2d

  subroutine make_enthalpy_3d(enthalpy,s,lo,hi,ng)

    use variables, only: rho_comp, rhoh_comp

    integer, intent(in)               :: lo(:),hi(:),ng
    real (kind = dp_t), intent(  out) :: enthalpy(lo(1):,lo(2):,lo(3):)  
    real (kind = dp_t), intent(in   ) :: s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)

    ! Local variables
    integer :: i, j, k, comp

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             enthalpy(i,j,k) = s(i,j,k,rhoh_comp)/s(i,j,k,rho_comp)
          enddo
       enddo
    end do

  end subroutine make_enthalpy_3d

  subroutine make_tfromH(n,plotdata,comp_t,comp_dp,state,p0,tempbar,dx)

    use geometry, only: spherical

    integer        , intent(in   ) :: n,comp_t,comp_dp
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: state
    real(kind=dp_t), intent(in   ) :: p0(0:)
    real(kind=dp_t), intent(in   ) :: tempbar(0:)
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
          call make_tfromH_2d(tp(:,:,1,comp_t),tp(:,:,1,comp_dp),sp(:,:,1,:), lo, hi, &
                              ng, p0, tempbar)
       case (3)
          if (spherical .eq. 1) then
             call make_tfromH_3d_sphr(n,tp(:,:,:,comp_t),tp(:,:,:,comp_dp),sp(:,:,:,:), &
                                      lo, hi, ng, p0, tempbar, dx)
          else
             call make_tfromH_3d_cart(tp(:,:,:,comp_t),tp(:,:,:,comp_dp),sp(:,:,:,:), &
                                      lo, hi, ng, p0, tempbar)
          end if
       end select
    end do

  end subroutine make_tfromH

  subroutine make_tfromH_2d(T,deltaP,state,lo,hi,ng,p0,tempbar)

    use variables, only: rho_comp, spec_comp, rhoh_comp
    use eos_module
    use bl_constants_module

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) ::      T(lo(1)   :,lo(2):)  
    real (kind = dp_t), intent(  out) :: deltaP(lo(1)   :,lo(2):)  
    real (kind = dp_t), intent(in   ) ::  state(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::  p0(0:)
    real (kind = dp_t), intent(in   ) ::  tempbar(0:)

    ! Local variables
    integer :: i, j, comp

    do_diag = .false.

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          ! (rho, H) --> T, p

          den_eos(1)  = state(i,j,rho_comp)
          p_eos(1)    = p0(j)
          temp_eos(1) = tempbar(j)
          xn_eos(1,:) = state(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)
          h_eos(1) = state(i,j,rhoh_comp) / state(i,j,rho_comp)

          call eos(eos_input_rh, den_eos, temp_eos, &
                   npts, nspec, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, &
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   do_diag)

          T(i,j) = temp_eos(1)

          deltaP(i,j) = abs(p_eos(1)-p0(j))/ p0(j)

       enddo
    enddo

  end subroutine make_tfromH_2d

  subroutine make_tfromH_3d_cart(T,deltaP,state,lo,hi,ng,p0,tempbar)

    use variables, only: rho_comp, spec_comp, rhoh_comp
    use eos_module

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) ::      T(lo(1)   :,lo(2):   ,lo(3):     )  
    real (kind = dp_t), intent(  out) :: deltaP(lo(1)   :,lo(2):   ,lo(3):     )  
    real (kind = dp_t), intent(in   ) :: state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) :: p0(0:)
    real (kind = dp_t), intent(in   ) :: tempbar(0:)

    ! Local variables
    integer :: i, j, k, comp

    do_diag = .false.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! (rho, H) --> T, p
             den_eos(1)  = state(i,j,k,rho_comp)
             p_eos(1)    = p0(k)
             temp_eos(1) = tempbar(k)
             xn_eos(1,:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)
             h_eos(1) = state(i,j,k,rhoh_comp)/state(i,j,k,rho_comp)

             call eos(eos_input_rh, den_eos, temp_eos, &
                      npts, nspec, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      do_diag)

             ! T(i,j,k) = log(temp_eos(1))/log(10.)
             T(i,j,k) = temp_eos(1)

             deltaP(i,j,k) = (p_eos(1)-p0(k))/ p0(k)

          enddo
       enddo
    enddo

  end subroutine make_tfromH_3d_cart

  subroutine make_tfromH_3d_sphr(n,T,deltaP,state,lo,hi,ng,p0,tempbar,dx)

    use variables, only: rho_comp, rhoh_comp, spec_comp
    use eos_module
    use fill_3d_module

    integer, intent(in) :: n, lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) ::      T(lo(1)   :,lo(2):   ,lo(3):     )  
    real (kind = dp_t), intent(  out) :: deltaP(lo(1)   :,lo(2):   ,lo(3):     )  
    real (kind = dp_t), intent(in   ) :: state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) ::    p0(0:)
    real (kind = dp_t), intent(in   ) ::    tempbar(0:)
    real (kind = dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i, j, k
    real (kind=dp_t), allocatable :: tempbar_cart(:,:,:)
    real (kind=dp_t), allocatable :: p0_cart(:,:,:)

    allocate(tempbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(n,tempbar_cart,tempbar,lo,hi,dx,0)

    allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(n,p0_cart,p0,lo,hi,dx,0)

    do_diag = .false.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_eos(1)  = state(i,j,k,rho_comp)
             h_eos(1)    = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)
             p_eos(1)    = p0_cart(i,j,k)
             temp_eos(1) = tempbar_cart(i,j,k)
             xn_eos(1,:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

             ! (rho, H) --> T, p
             call eos(eos_input_rh, den_eos, temp_eos, &
                      npts, nspec, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      do_diag)

             T(i,j,k) = temp_eos(1)

             deltaP(i,j,k) = (p_eos(1)-p0_cart(i,j,k))/ p0_cart(i,j,k)

          enddo
       enddo
    enddo

    deallocate(tempbar_cart,p0_cart)

  end subroutine make_tfromH_3d_sphr

  subroutine make_tfromrho(n,plotdata,comp_tfromrho,comp_tpert,comp_rhopert, &
                           comp_machno,comp_deltag,s,u,rho0,tempbar,gamma1bar,p0,dx)

    use geometry, only: spherical

    integer        , intent(in   ) :: n,comp_tfromrho,comp_tpert
    integer        , intent(in   ) :: comp_rhopert, comp_machno
    integer        , intent(in   ) :: comp_deltag
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: s
    type(multifab) , intent(in   ) :: u
    real(kind=dp_t), intent(in   ) :: rho0(0:)
    real(kind=dp_t), intent(in   ) :: tempbar(0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real(kind=dp_t), intent(in   ) :: p0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
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
          call make_tfromrho_2d(tp(:,:,1,comp_tfromrho),tp(:,:,1,comp_tpert), &
                                tp(:,:,1,comp_rhopert ), &
                                tp(:,:,1,comp_machno  ),tp(:,:,1,comp_deltag), &
                                sp(:,:,1,:), up(:,:,1,:), &
                                lo, hi, ng, rho0, tempbar, gamma1bar, p0)
       case (3)
          if (spherical .eq. 1) then
             call make_tfromrho_3d_sphr(n,tp(:,:,:,comp_tfromrho),tp(:,:,:,comp_tpert), &
                                        tp(:,:,:,comp_rhopert ), &
                                        tp(:,:,:,comp_machno  ),tp(:,:,:,comp_deltag), &
                                        sp(:,:,:,:), up(:,:,:,:), &
                                        lo, hi, ng, rho0, tempbar, gamma1bar, p0, dx)
          else
             call make_tfromrho_3d_cart(tp(:,:,:,comp_tfromrho),tp(:,:,:,comp_tpert), &
                                        tp(:,:,:,comp_rhopert ), &
                                        tp(:,:,:,comp_machno  ),tp(:,:,:,comp_deltag), &
                                        sp(:,:,:,:), up(:,:,:,:), &
                                        lo, hi, ng, rho0, tempbar, gamma1bar, p0)
          endif
       end select
    end do

  end subroutine make_tfromrho

  subroutine make_tfromrho_2d(t,tpert,rhopert,machno,deltagamma,s,u,lo,hi,ng,rho0,tempbar, &
                              gamma1bar,p0)

    use eos_module
    use variables, only: rho_comp, spec_comp

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(  out)     ::      t(lo(1):,lo(2):)  
    real (kind=dp_t), intent(  out) ::      tpert(lo(1):,lo(2):)  
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1):,lo(2):)  
    real (kind=dp_t), intent(  out) ::     machno(lo(1):,lo(2):)  
    real (kind=dp_t), intent(  out) :: deltagamma(lo(1):,lo(2):)  
    real (kind=dp_t), intent(in   ) ::  s(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(in   ) ::  u(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:)
    real (kind=dp_t), intent(in   ) :: tempbar(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)

    !     Local variables
    integer          :: i, j
    real (kind=dp_t) :: vel

    do_diag = .false.

    ! Then compute the perturbation
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          den_eos(1) = s(i,j,rho_comp)
          temp_eos(1) = tempbar(j)
          p_eos(1) = p0(j)
          xn_eos(1,:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)

          ! (rho,P) --> T,h
          call eos(eos_input_rp, den_eos, temp_eos, &
                   npts, nspec, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, & 
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   do_diag)

          t(i,j) = temp_eos(1)
          tpert(i,j) = temp_eos(1) - tempbar(j)

          rhopert(i,j) = s(i,j,rho_comp) - rho0(j)

          vel = sqrt(u(i,j,1)*u(i,j,1) + u(i,j,2)*u(i,j,2))
          machno(i,j) = vel / cs_eos(1)

          deltagamma(i,j) = gam1_eos(1) - gamma1bar(j)
       enddo
    enddo

  end subroutine make_tfromrho_2d

  subroutine make_tfromrho_3d_cart(t,tpert,rhopert,machno,deltagamma,s,u,lo,hi, &
                                   ng,rho0,tempbar,gamma1bar,p0)

    use variables, only: rho_comp, spec_comp
    use eos_module

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(  out) ::          t(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::      tpert(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::     machno(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) :: deltagamma(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(in   ) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(in   ) ::  u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:)
    real (kind=dp_t), intent(in   ) :: tempbar(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)

    ! Local variables
    integer          :: i, j, k
    real (kind=dp_t) :: vel

    do_diag = .false.

    ! Then compute the perturbation and Mach number
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_eos(1) = s(i,j,k,rho_comp)
             temp_eos(1) = tempbar(k)
             p_eos(1) = p0(k)
             xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

             ! (rho,P) --> T,h
             call eos(eos_input_rp, den_eos, temp_eos, &
                      npts, nspec, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, & 
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      do_diag)

             t(i,j,k) = temp_eos(1)
             tpert(i,j,k) = temp_eos(1) - tempbar(k)

             rhopert(i,j,k) = s(i,j,k,rho_comp) - rho0(k)

             vel = sqrt(u(i,j,k,1)*u(i,j,k,1) + u(i,j,k,2)*u(i,j,k,2) + u(i,j,k,3)*u(i,j,k,3))
             machno(i,j,k) = vel / cs_eos(1)

             deltagamma(i,j,k) = gam1_eos(1) - gamma1bar(k)
          enddo
       enddo
    enddo

  end subroutine make_tfromrho_3d_cart

  subroutine make_tfromrho_3d_sphr(n,t,tpert,rhopert,machno,deltagamma, &
                                   s,u,lo,hi,ng,rho0,tempbar,gamma1bar,p0,dx)

    use geometry, only: nr
    use variables, only: rho_comp, spec_comp
    use eos_module
    use fill_3d_module

    integer, intent(in)             :: n,lo(:),hi(:),ng
    real (kind=dp_t), intent(  out) ::          t(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::      tpert(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::     machno(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) :: deltagamma(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(in   ) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(in   ) ::  u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:)
    real (kind=dp_t), intent(in   ) :: tempbar(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer          :: i, j, k, r
    real (kind=dp_t) :: vel
    real (kind=dp_t), allocatable ::  rho0_cart(:,:,:)
    real (kind=dp_t), allocatable ::    tempbar_cart(:,:,:)
    real (kind=dp_t), allocatable ::    p0_cart(:,:,:)
    real (kind=dp_t), allocatable ::  gam0_cart(:,:,:)

    do_diag = .false.

    allocate(rho0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(n,rho0_cart,rho0,lo,hi,dx,0)

    allocate(tempbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(n,tempbar_cart,tempbar,lo,hi,dx,0)

    allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(n,p0_cart,p0,lo,hi,dx,0)

    allocate(gam0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(n,gam0_cart,gamma1bar,lo,hi,dx,0)

    ! Then compute the perturbation and Mach number
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_eos(1) = s(i,j,k,rho_comp)
             temp_eos(1) = tempbar_cart(i,j,k)
             p_eos(1) = p0_cart(i,j,k)
             xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

             ! (rho,P) --> T,h
             call eos(eos_input_rp, den_eos, temp_eos, &
                      npts, nspec, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, & 
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      do_diag)

             t(i,j,k) = temp_eos(1)
             tpert(i,j,k) = temp_eos(1) - tempbar_cart(i,j,k)

             rhopert(i,j,k) = s(i,j,k,rho_comp) - rho0_cart(i,j,k)

             vel = sqrt(u(i,j,k,1)*u(i,j,k,1)+u(i,j,k,2)*u(i,j,k,2)+u(i,j,k,3)*u(i,j,k,3))
             machno(i,j,k) = vel / cs_eos(1)

             deltagamma(i,j,k) = gam1_eos(1) - gam0_cart(i,j,k)
          enddo
       enddo
    enddo

    deallocate(rho0_cart,tempbar_cart,p0_cart,gam0_cart)

  end subroutine make_tfromrho_3d_sphr

  subroutine make_XfromrhoX(plotdata,comp,s)

    use network, only: nspec
    use variables, only: spec_comp

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
          call make_XfromrhoX_2d(pp(:,:,1,comp:),sp(:,:,1,1), &
                                 sp(:,:,1,spec_comp:spec_comp+nspec-1), lo, hi, ng)
       case (3)
          call make_XfromrhoX_3d(pp(:,:,:,comp:),sp(:,:,:,1), &
                                 sp(:,:,:,spec_comp:spec_comp+nspec-1), lo, hi, ng)
       end select
    end do

  end subroutine make_XfromrhoX

  subroutine make_XfromrhoX_2d(X,rho,rhoX,lo,hi,ng)

    use network, only: nspec

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) ::    X(lo(1)   :,lo(2)   :,:)  
    real (kind = dp_t), intent(in   ) :: rhoX(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::  rho(lo(1)-ng:,lo(2)-ng:  )

    !     Local variables
    integer :: i, j, comp

    do comp = 1, nspec
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             X(i,j,comp) = rhoX(i,j,comp) / rho(i,j)
          enddo
       enddo
    enddo

  end subroutine make_XfromrhoX_2d

  subroutine make_XfromrhoX_3d(X,rho,rhoX,lo,hi,ng)

    use network, only: nspec

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) ::    X(lo(1)   :,lo(2)   :,lo(3)   :,:)  
    real (kind = dp_t), intent(in   ) :: rhoX(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) ::  rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:  )

    !     Local variables
    integer :: i, j, k, comp

    do comp = 1, nspec
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                X(i,j,k,comp) = rhoX(i,j,k,comp) / rho(i,j,k)
             enddo
          enddo
       enddo
    enddo

  end subroutine make_XfromrhoX_3d


  subroutine make_omegadot(plotdata,comp,comp_enuc,s,rho_omegadot)

    use network, only: nspec
    use variables, only: rho_comp

    integer        , intent(in   ) :: comp, comp_enuc
    type(multifab) , intent(in   ) :: s, rho_omegadot
    type(multifab) , intent(inout) :: plotdata

    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: rop(:,:,:,:)
    real(kind=dp_t), pointer:: pp(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng_s,ng_o,dm
    integer :: i

    ng_s = s%ng
    ng_o = rho_omegadot%ng
    dm = s%dim

    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       sp => dataptr(s, i)
       rop => dataptr(rho_omegadot, i)
       pp => dataptr(plotdata, i)

       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))

       select case (dm)
       case (2)
          call make_omega_2d(pp(:,:,1,comp:comp-1+nspec),pp(:,:,1,comp_enuc), &
                             rop(:,:,1,:),sp(:,:,1,rho_comp), lo, hi, ng_s, ng_o)
       case (3)
          call make_omega_3d(pp(:,:,:,comp:comp-1+nspec),pp(:,:,:,comp_enuc), &
                             rop(:,:,:,:),sp(:,:,:,rho_comp), lo, hi, ng_s, ng_o)
       end select
    end do

  end subroutine make_omegadot


  subroutine make_omega_2d(omegadot,enuc,rho_omegadot,rho,lo,hi,ng_s,ng_o)

    use network, only: nspec, ebin
    use bl_constants_module

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_o
    real (kind = dp_t), intent(  out) :: omegadot(lo(1):,lo(2):,:)
    real (kind = dp_t), intent(  out) :: enuc(lo(1):,lo(2):)
    real (kind = dp_t), intent(in   ) :: rho_omegadot(lo(1)-ng_o:,lo(2)-ng_o:,:)
    real (kind = dp_t), intent(in   ) :: rho(lo(1)-ng_s:,lo(2)-ng_s:)

    ! Local variables
    integer :: i, j, comp

    enuc(:,:) = ZERO

    do comp = 1, nspec
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             omegadot(i,j,comp) = rho_omegadot(i,j,comp) / rho(i,j)
             enuc(i,j) = enuc(i,j) - omegadot(i,j,comp)*ebin(comp)
          enddo
       enddo
    enddo

  end subroutine make_omega_2d

  subroutine make_omega_3d(omegadot,enuc,rho_omegadot,rho,lo,hi,ng_s,ng_o)

    use network, only: nspec, ebin
    use bl_constants_module

    integer, intent(in) :: lo(:), hi(:), ng_s, ng_o
    real (kind = dp_t), intent(  out) :: omegadot(lo(1):,lo(2):,lo(3):,:)  
    real (kind = dp_t), intent(  out) :: enuc(lo(1):,lo(2):,lo(3):)
    real (kind = dp_t), intent(in   ) :: rho_omegadot(lo(1)-ng_o:,lo(2)-ng_o:,lo(3)-ng_o:,:)
    real (kind = dp_t), intent(in   ) :: rho(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:  )

    ! Local variables
    integer :: i, j, k, comp

    enuc(:,:,:) = ZERO

    do comp = 1, nspec
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                omegadot(i,j,k,comp) = rho_omegadot(i,j,k,comp) / rho(i,j,k)
                enuc(i,j,k) = enuc(i,j,k) - omegadot(i,j,k,comp)*ebin(comp)
             enddo
          enddo
       enddo
    enddo

  end subroutine make_omega_3d

  subroutine make_deltaT(plotdata,comp_dT,comp_tfromH,comp_tfromrho)

    integer        , intent(in   ) :: comp_dT, comp_tfromH, comp_tfromrho
    type(multifab) , intent(inout) :: plotdata

    real(kind=dp_t), pointer:: tp(:,:,:,:)
    integer :: lo(plotdata%dim),hi(plotdata%dim)
    integer :: i,dm

    dm = plotdata%dim

    do i = 1, plotdata%nboxes
       if ( multifab_remote(plotdata, i) ) cycle
       tp => dataptr(plotdata, i)
       lo =  lwb(get_box(plotdata, i))
       hi =  upb(get_box(plotdata, i))
       select case (dm)
       case (2)
          call make_deltaT_2d(tp(:,:,1,comp_dT),tp(:,:,1,comp_tfromH), &
                              tp(:,:,1,comp_tfromrho), lo, hi)
       case (3)
          call make_deltaT_3d(tp(:,:,:,comp_dT),tp(:,:,:,comp_tfromH), &
                              tp(:,:,:,comp_tfromrho), lo, hi)
       end select
    end do

  end subroutine make_deltaT

  subroutine make_deltaT_2d(dT,tfromH,tfromrho,lo,hi)

    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(  out) ::       dT(lo(1):,lo(2):)
    real (kind = dp_t), intent(in   ) ::   tfromH(lo(1):,lo(2):)
    real (kind = dp_t), intent(in   ) :: tfromrho(lo(1):,lo(2):)

    !     Local variables
    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          dT(i,j) = tfromH(i,j) - tfromrho(i,j)
       end do
    end do

  end subroutine make_deltaT_2d

  subroutine make_deltaT_3d(dT,tfromH,tfromrho,lo,hi)

    integer, intent(in) :: lo(:), hi(:)
    real (kind = dp_t), intent(  out) ::       dT(lo(1):,lo(2):,lo(3):)
    real (kind = dp_t), intent(in   ) ::   tfromH(lo(1):,lo(2):,lo(3):)
    real (kind = dp_t), intent(in   ) :: tfromrho(lo(1):,lo(2):,lo(3):)

    !     Local variables
    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dT(i,j,k) = tfromH(i,j,k) - tfromrho(i,j,k)
          end do
       end do
    end do

  end subroutine make_deltaT_3d

end module plot_variables_module

