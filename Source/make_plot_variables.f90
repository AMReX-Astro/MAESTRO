module plot_variables_module

  use bl_types
  use multifab_module

  implicit none

  private

  public :: make_enthalpy, make_tfromH, make_tfromp, make_XfromrhoX, make_entropypert
  public :: make_omegadot, make_deltaT, make_divw0

contains

  subroutine make_enthalpy(enthalpy,comp,s)

    integer        , intent(in   ) :: comp
    type(multifab) , intent(inout) :: enthalpy
    type(multifab) , intent(in   ) :: s

    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: rp(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng_s,ng_e,dm
    integer :: i

    ng_s = s%ng
    ng_e = enthalpy%ng
    dm = s%dim

    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       sp => dataptr(s, i)
       rp => dataptr(enthalpy, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (2)
          call make_enthalpy_2d(rp(:,:,1,comp),ng_e,sp(:,:,1,:),ng_s,lo,hi)
       case (3)
          call make_enthalpy_3d(rp(:,:,:,comp),ng_e,sp(:,:,:,:),ng_s,lo,hi)
       end select
    end do

  end subroutine make_enthalpy

  subroutine make_enthalpy_2d(enthalpy,ng_e,s,ng_s,lo,hi)

    use variables, only: rho_comp, rhoh_comp

    integer, intent(in) :: lo(:), hi(:), ng_e, ng_s
    real (kind = dp_t), intent(  out) :: enthalpy(lo(1)-ng_e:,lo(2)-ng_e:)  
    real (kind = dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,:)

    ! Local variables
    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          enthalpy(i,j) = s(i,j,rhoh_comp)/s(i,j,rho_comp)
       enddo
    enddo

  end subroutine make_enthalpy_2d

  subroutine make_enthalpy_3d(enthalpy,ng_e,s,ng_s,lo,hi)

    use variables, only: rho_comp, rhoh_comp

    integer, intent(in)               :: lo(:),hi(:),ng_e,ng_s
    real (kind = dp_t), intent(  out) :: enthalpy(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)  
    real (kind = dp_t), intent(in   ) ::        s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)

    ! Local variables
    integer :: i, j, k

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
    integer :: lo(state%dim),hi(state%dim),ng_s,ng_p,dm
    integer :: i

    ng_s = state%ng
    ng_p = plotdata%ng
    dm = state%dim

    do i = 1, state%nboxes
       if ( multifab_remote(state, i) ) cycle
       sp => dataptr(state, i)
       tp => dataptr(plotdata, i)
       lo =  lwb(get_box(state, i))
       hi =  upb(get_box(state, i))
       select case (dm)
       case (2)
          call make_tfromH_2d(tp(:,:,1,comp_t),tp(:,:,1,comp_dp),ng_p,sp(:,:,1,:),ng_s, &
                              lo,hi,p0,tempbar)
       case (3)
          if (spherical .eq. 1) then
             call make_tfromH_3d_sphr(n,tp(:,:,:,comp_t),tp(:,:,:,comp_dp),ng_p, &
                                      sp(:,:,:,:),ng_s,lo,hi,p0,tempbar,dx)
          else
             call make_tfromH_3d_cart(tp(:,:,:,comp_t),tp(:,:,:,comp_dp),ng_p, &
                                      sp(:,:,:,:),ng_s,lo,hi,p0,tempbar)
          end if
       end select
    end do

  end subroutine make_tfromH

  subroutine make_tfromH_2d(T,deltaP,ng_p,state,ng_s,lo,hi,p0,tempbar)

    use variables, only: rho_comp, spec_comp, rhoh_comp
    use eos_module
    use bl_constants_module

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind = dp_t), intent(  out) ::      T(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind = dp_t), intent(  out) :: deltaP(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind = dp_t), intent(in   ) ::  state(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind = dp_t), intent(in   ) ::  p0(0:)
    real (kind = dp_t), intent(in   ) ::  tempbar(0:)

    ! Local variables
    integer :: i, j

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

  subroutine make_tfromH_3d_cart(T,deltaP,ng_p,state,ng_s,lo,hi,p0,tempbar)

    use variables, only: rho_comp, spec_comp, rhoh_comp
    use eos_module

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind = dp_t), intent(  out) ::      T(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind = dp_t), intent(  out) :: deltaP(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind = dp_t), intent(in   ) ::  state(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind = dp_t), intent(in   ) :: p0(0:)
    real (kind = dp_t), intent(in   ) :: tempbar(0:)

    ! Local variables
    integer :: i, j, k

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

  subroutine make_tfromH_3d_sphr(n,T,deltaP,ng_p,state,ng_s,lo,hi,p0,tempbar,dx)

    use variables, only: rho_comp, rhoh_comp, spec_comp
    use eos_module
    use fill_3d_module

    integer, intent(in) :: n, lo(:), hi(:), ng_p, ng_s
    real (kind = dp_t), intent(  out) ::      T(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind = dp_t), intent(  out) :: deltaP(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind = dp_t), intent(in   ) ::  state(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind = dp_t), intent(in   ) ::    p0(0:)
    real (kind = dp_t), intent(in   ) ::    tempbar(0:)
    real (kind = dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i, j, k
    real (kind=dp_t), allocatable :: tempbar_cart(:,:,:,:)
    real (kind=dp_t), allocatable :: p0_cart(:,:,:,:)

    allocate(tempbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,tempbar,tempbar_cart,lo,hi,dx,0,0)

    allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,p0,p0_cart,lo,hi,dx,0,0)
    do_diag = .false.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_eos(1)  = state(i,j,k,rho_comp)
             h_eos(1)    = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)
             p_eos(1)    = p0_cart(i,j,k,1)
             temp_eos(1) = tempbar_cart(i,j,k,1)
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

             deltaP(i,j,k) = (p_eos(1)-p0_cart(i,j,k,1))/ p0_cart(i,j,k,1)

          enddo
       enddo
    enddo

    deallocate(tempbar_cart,p0_cart)

  end subroutine make_tfromH_3d_sphr


  subroutine make_tfromp(n,plotdata, &
                         comp_tfromp,comp_tpert,comp_rhopert,comp_rhohpert, &
                         comp_machno,comp_deltag,comp_entropy, &
                         s,u,rho0,rhoh0,tempbar,gamma1bar,p0,dx)

    use geometry, only: spherical

    integer        , intent(in   ) :: n,comp_tfromp,comp_tpert
    integer        , intent(in   ) :: comp_rhopert, comp_rhohpert, comp_machno
    integer        , intent(in   ) :: comp_deltag, comp_entropy
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: s
    type(multifab) , intent(in   ) :: u
    real(kind=dp_t), intent(in   ) :: rho0(0:)
    real(kind=dp_t), intent(in   ) :: rhoh0(0:)
    real(kind=dp_t), intent(in   ) :: tempbar(0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real(kind=dp_t), intent(in   ) :: p0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), pointer:: sp(:,:,:,:),tp(:,:,:,:),up(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),i,dm
    integer :: ng_p,ng_s,ng_u

    ng_p = plotdata%ng
    ng_s = s%ng
    ng_u = u%ng
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
          call make_tfromp_2d(tp(:,:,1,comp_tfromp),tp(:,:,1,comp_tpert), &
                              tp(:,:,1,comp_rhopert ),tp(:,:,1,comp_rhohpert), &
                              tp(:,:,1,comp_machno  ),tp(:,:,1,comp_deltag), &
                              tp(:,:,1,comp_entropy ), ng_p, &
                              sp(:,:,1,:), ng_s, up(:,:,1,:), ng_u, &
                              lo, hi, rho0, rhoh0, tempbar, gamma1bar, p0)
       case (3)
          if (spherical .eq. 1) then
             call make_tfromp_3d_sphr(n,tp(:,:,:,comp_tfromp),tp(:,:,:,comp_tpert), &
                                      tp(:,:,:,comp_rhopert ),tp(:,:,:,comp_rhohpert), &
                                      tp(:,:,:,comp_machno  ),tp(:,:,:,comp_deltag), &
                                      tp(:,:,:,comp_entropy ), ng_p, &
                                      sp(:,:,:,:), ng_s, up(:,:,:,:), ng_u, &
                                      lo, hi, rho0, rhoh0, tempbar, gamma1bar, p0, dx)
          else
             call make_tfromp_3d_cart(tp(:,:,:,comp_tfromp),tp(:,:,:,comp_tpert), &
                                      tp(:,:,:,comp_rhopert ),tp(:,:,:,comp_rhohpert), &
                                      tp(:,:,:,comp_machno  ),tp(:,:,:,comp_deltag), &
                                      tp(:,:,:,comp_entropy ), ng_p, &
                                      sp(:,:,:,:), ng_s, up(:,:,:,:), ng_u, &
                                      lo, hi, rho0, rhoh0, tempbar, gamma1bar, p0)
          endif
       end select
    end do

  end subroutine make_tfromp

  subroutine make_tfromp_2d(t,tpert,rhopert,rhohpert,machno,deltagamma,entropy, &
                            ng_p,s,ng_s,u,ng_u, &
                            lo,hi,rho0,rhoh0,tempbar,gamma1bar,p0)

    use eos_module
    use variables, only: rho_comp, rhoh_comp, spec_comp

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s, ng_u
    real (kind=dp_t), intent(  out) ::          t(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(  out) ::      tpert(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(  out) ::   rhohpert(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(  out) ::     machno(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(  out) :: deltagamma(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(  out) ::    entropy(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind=dp_t), intent(in   ) ::          u(lo(1)-ng_u:,lo(2)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:)
    real (kind=dp_t), intent(in   ) :: rhoh0(0:)
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

          rhopert(i,j)  = s(i,j,rho_comp)  - rho0(j)
          rhohpert(i,j) = s(i,j,rhoh_comp) - rhoh0(j)

          vel = sqrt(u(i,j,1)*u(i,j,1) + u(i,j,2)*u(i,j,2))
          machno(i,j) = vel / cs_eos(1)

          deltagamma(i,j) = gam1_eos(1) - gamma1bar(j)

          entropy(i,j) = s_eos(1)
       enddo
    enddo

  end subroutine make_tfromp_2d

  subroutine make_tfromp_3d_cart(t,tpert,rhopert,rhohpert,machno,deltagamma,entropy, &
                                 ng_p,s,ng_s,u,ng_u, &
                                 lo,hi,rho0,rhoh0,tempbar,gamma1bar,p0)

    use variables, only: rho_comp, rhoh_comp, spec_comp
    use eos_module

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s, ng_u
    real (kind=dp_t), intent(  out) ::          t(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::      tpert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::   rhohpert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::     machno(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) :: deltagamma(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::    entropy(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) ::          u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:)
    real (kind=dp_t), intent(in   ) :: rhoh0(0:)
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

             rhopert(i,j,k)  = s(i,j,k,rho_comp)  - rho0(k)
             rhohpert(i,j,k) = s(i,j,k,rhoh_comp) - rhoh0(k)

             vel = sqrt(u(i,j,k,1)*u(i,j,k,1)+u(i,j,k,2)*u(i,j,k,2)+u(i,j,k,3)*u(i,j,k,3))
             machno(i,j,k) = vel / cs_eos(1)

             deltagamma(i,j,k) = gam1_eos(1) - gamma1bar(k)

             entropy(i,j,k) = s_eos(1)
          enddo
       enddo
    enddo

  end subroutine make_tfromp_3d_cart

  subroutine make_tfromp_3d_sphr(n,t,tpert,rhopert,rhohpert,machno,deltagamma,entropy, &
                                 ng_p,s,ng_s,u,ng_u, &
                                 lo,hi,rho0,rhoh0,tempbar,gamma1bar,p0,dx)

    use variables, only: rho_comp, rhoh_comp, spec_comp
    use eos_module
    use fill_3d_module

    integer, intent(in)             :: n,lo(:),hi(:),ng_p,ng_s,ng_u
    real (kind=dp_t), intent(  out) ::          t(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::      tpert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::   rhohpert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::     machno(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) :: deltagamma(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::    entropy(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) ::          u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:)
    real (kind=dp_t), intent(in   ) :: rhoh0(0:)
    real (kind=dp_t), intent(in   ) :: tempbar(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer          :: i, j, k
    real (kind=dp_t) :: vel
    real (kind=dp_t), allocatable ::  rho0_cart(:,:,:,:)
    real (kind=dp_t), allocatable ::  rhoh0_cart(:,:,:,:)
    real (kind=dp_t), allocatable ::  tempbar_cart(:,:,:,:)
    real (kind=dp_t), allocatable ::  p0_cart(:,:,:,:)
    real (kind=dp_t), allocatable ::  gamma1bar_cart(:,:,:,:)

    do_diag = .false.

    allocate(rho0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,rho0,rho0_cart,lo,hi,dx,0,0)

    allocate(rhoh0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,rhoh0,rhoh0_cart,lo,hi,dx,0,0)

    allocate(tempbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,tempbar,tempbar_cart,lo,hi,dx,0,0)

    allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)) 
    call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,p0,p0_cart,lo,hi,dx,0,0)

    allocate(gamma1bar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)) 
    call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,gamma1bar,gamma1bar_cart,lo,hi, &
                                      dx,0,0)

    ! Then compute the perturbation and Mach number
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_eos(1) = s(i,j,k,rho_comp)
             temp_eos(1) = tempbar_cart(i,j,k,1)
             p_eos(1) = p0_cart(i,j,k,1)
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
             tpert(i,j,k) = temp_eos(1) - tempbar_cart(i,j,k,1)

             rhopert(i,j,k)  = s(i,j,k,rho_comp)  -  rho0_cart(i,j,k,1)
             rhohpert(i,j,k) = s(i,j,k,rhoh_comp) - rhoh0_cart(i,j,k,1)

             vel = sqrt(u(i,j,k,1)*u(i,j,k,1)+u(i,j,k,2)*u(i,j,k,2)+u(i,j,k,3)*u(i,j,k,3))
             machno(i,j,k) = vel / cs_eos(1)

             deltagamma(i,j,k) = gam1_eos(1) - gamma1bar_cart(i,j,k,1)

             entropy(i,j,k) = s_eos(1)
          enddo
       enddo
    enddo

    deallocate(rho0_cart,rhoh0_cart,tempbar_cart,p0_cart,gamma1bar_cart)

  end subroutine make_tfromp_3d_sphr


  subroutine make_entropypert(n,plotdata,comp_entropy,comp_entropypert,entropybar,dx)

    use geometry, only: spherical

    integer        , intent(in   ) :: n,comp_entropy,comp_entropypert
    type(multifab) , intent(inout) :: plotdata
    real(kind=dp_t), intent(in   ) :: entropybar(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    real(kind=dp_t), pointer:: tp(:,:,:,:)
    integer :: lo(plotdata%dim),hi(plotdata%dim),dm,ng_p
    integer :: i

    dm = plotdata%dim
    ng_p = plotdata%ng

    do i = 1, plotdata%nboxes
       if ( multifab_remote(plotdata, i) ) cycle
       tp => dataptr(plotdata, i)
       lo =  lwb(get_box(plotdata, i))
       hi =  upb(get_box(plotdata, i))
       select case (dm)
       case (2)
          call make_entropypert_2d(tp(:,:,1,comp_entropy), &
                                   tp(:,:,1,comp_entropypert),ng_p, &
                                   lo, hi, entropybar)
       case (3)
          if (spherical .eq. 1) then
             call make_entropypert_3d_sphr(n,tp(:,:,:,comp_entropy), &
                                             tp(:,:,:,comp_entropypert),ng_p, &
                                           lo, hi, entropybar, dx)
          else
             call make_entropypert_3d_cart(tp(:,:,:,comp_entropy), &
                                           tp(:,:,:,comp_entropypert),ng_p, &
                                           lo, hi, entropybar)
          endif
       end select
    end do

  end subroutine make_entropypert

  subroutine make_entropypert_2d(entropy,entropypert,ng_p,lo,hi,entropybar)

    integer, intent(in) :: lo(:), hi(:), ng_p
    real (kind=dp_t), intent(inout) ::     entropy(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(  out) :: entropypert(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(in   ) :: entropybar(0:)

    !     Local variables
    integer          :: i, j

    ! Compute entropy - entropybar
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          entropypert(i,j) = (entropy(i,j) - entropybar(j))/entropybar(j)
       enddo
    enddo

  end subroutine make_entropypert_2d

  subroutine make_entropypert_3d_cart(entropy,entropypert,ng_p,lo,hi,entropybar)

    integer, intent(in) :: lo(:), hi(:), ng_p
    real (kind=dp_t), intent(inout) ::     entropy(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) :: entropypert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(in   ) :: entropybar(0:)

    ! Local variables
    integer          :: i, j, k


    ! Compute entropy - entropybar
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             entropypert(i,j,k) = (entropy(i,j,k) - entropybar(k))/entropybar(k)
          enddo
       enddo
    enddo

  end subroutine make_entropypert_3d_cart

  subroutine make_entropypert_3d_sphr(n,entropy,entropypert,ng_p,lo,hi,entropybar,dx)

    use fill_3d_module

    integer, intent(in)             :: n,lo(:),hi(:),ng_p
    real (kind=dp_t), intent(inout) ::     entropy(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) :: entropypert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(in   ) :: entropybar(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer          :: i, j, k
    real (kind=dp_t), allocatable ::  entropybar_cart(:,:,:,:)

    allocate(entropybar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1)) 
    call put_1d_array_on_cart_3d_sphr(n,.false.,.false.,entropybar,entropybar_cart,lo,hi, &
                                      dx,0,0)

    ! Compute entropy-entropybar
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             entropypert(i,j,k) = (entropy(i,j,k) - entropybar_cart(i,j,k,1))/entropybar_cart(i,j,k,1)
          enddo
       enddo
    enddo

    deallocate(entropybar_cart)

  end subroutine make_entropypert_3d_sphr

  subroutine make_XfromrhoX(plotdata,comp,s)

    use network, only: nspec
    use variables, only: spec_comp

    integer        , intent(in   ) :: comp
    type(multifab) , intent(in   ) :: s
    type(multifab) , intent(inout) :: plotdata

    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: pp(:,:,:,:)
    integer :: lo(s%dim),hi(s%dim),ng_s,ng_p,dm
    integer :: i

    ng_s = s%ng
    ng_p = plotdata%ng
    dm = s%dim

    do i = 1, s%nboxes
       if ( multifab_remote(s, i) ) cycle
       sp => dataptr(s, i)
       pp => dataptr(plotdata, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (2)
          call make_XfromrhoX_2d(pp(:,:,1,comp:),ng_p,sp(:,:,1,1), &
                                 sp(:,:,1,spec_comp:spec_comp+nspec-1),ng_s,lo,hi)
       case (3)
          call make_XfromrhoX_3d(pp(:,:,:,comp:),ng_p,sp(:,:,:,1), &
                                 sp(:,:,:,spec_comp:spec_comp+nspec-1),ng_s,lo,hi)
       end select
    end do

  end subroutine make_XfromrhoX

  subroutine make_XfromrhoX_2d(X,ng_p,rho,rhoX,ng_s,lo,hi)

    use network, only: nspec

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind = dp_t), intent(  out) ::    X(lo(1)-ng_p:,lo(2)-ng_p:,:)  
    real (kind = dp_t), intent(in   ) :: rhoX(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind = dp_t), intent(in   ) ::  rho(lo(1)-ng_s:,lo(2)-ng_s:)

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

  subroutine make_XfromrhoX_3d(X,ng_p,rho,rhoX,ng_s,lo,hi)

    use network, only: nspec

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind = dp_t), intent(  out) ::    X(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,:)  
    real (kind = dp_t), intent(in   ) :: rhoX(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind = dp_t), intent(in   ) ::  rho(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)

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
    integer :: lo(s%dim),hi(s%dim),ng_s,ng_o,ng_p,dm
    integer :: i

    ng_s = s%ng
    ng_o = rho_omegadot%ng
    ng_p = plotdata%ng
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
          call make_omega_2d(pp(:,:,1,comp:comp-1+nspec),pp(:,:,1,comp_enuc),ng_p, &
                             rop(:,:,1,:),ng_o,sp(:,:,1,rho_comp),ng_s,lo,hi)
       case (3)
          call make_omega_3d(pp(:,:,:,comp:comp-1+nspec),pp(:,:,:,comp_enuc),ng_p, &
                             rop(:,:,:,:),ng_o,sp(:,:,:,rho_comp),ng_s,lo,hi)
       end select
    end do

  end subroutine make_omegadot


  subroutine make_omega_2d(omegadot,enuc,ng_p,rho_omegadot,ng_o,rho,ng_s,lo,hi)

    use network, only: nspec, ebin
    use bl_constants_module

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_o, ng_s
    real (kind = dp_t), intent(  out) ::     omegadot(lo(1)-ng_p:,lo(2)-ng_p:,:)
    real (kind = dp_t), intent(  out) ::         enuc(lo(1)-ng_p:,lo(2)-ng_p:)
    real (kind = dp_t), intent(in   ) :: rho_omegadot(lo(1)-ng_o:,lo(2)-ng_o:,:)
    real (kind = dp_t), intent(in   ) ::          rho(lo(1)-ng_s:,lo(2)-ng_s:)

    ! Local variables
    integer :: i, j, comp

    enuc = ZERO

    do comp = 1, nspec
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             omegadot(i,j,comp) = rho_omegadot(i,j,comp) / rho(i,j)
             enuc(i,j) = enuc(i,j) - omegadot(i,j,comp)*ebin(comp)
          enddo
       enddo
    enddo

  end subroutine make_omega_2d

  subroutine make_omega_3d(omegadot,enuc,ng_p,rho_omegadot,ng_o,rho,ng_s,lo,hi)

    use network, only: nspec, ebin
    use bl_constants_module

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_o, ng_s
    real (kind = dp_t), intent(  out) ::     omegadot(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,:)
    real (kind = dp_t), intent(  out) ::         enuc(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real (kind = dp_t), intent(in   ) :: rho_omegadot(lo(1)-ng_o:,lo(2)-ng_o:,lo(3)-ng_o:,:)
    real (kind = dp_t), intent(in   ) ::          rho(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)

    ! Local variables
    integer :: i, j, k, comp

    enuc = ZERO

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

  subroutine make_deltaT(plotdata,comp_dT,comp_tfromH,comp_tfromp)

    integer        , intent(in   ) :: comp_dT, comp_tfromH, comp_tfromp
    type(multifab) , intent(inout) :: plotdata

    real(kind=dp_t), pointer:: tp(:,:,:,:)
    integer :: lo(plotdata%dim),hi(plotdata%dim)
    integer :: i,dm,ng_p

    ng_p = plotdata%ng
    dm = plotdata%dim

    do i = 1, plotdata%nboxes
       if ( multifab_remote(plotdata, i) ) cycle
       tp => dataptr(plotdata, i)
       lo =  lwb(get_box(plotdata, i))
       hi =  upb(get_box(plotdata, i))
       select case (dm)
       case (2)
          call make_deltaT_2d(tp(:,:,1,comp_dT),tp(:,:,1,comp_tfromH), &
                              tp(:,:,1,comp_tfromp), ng_p, lo, hi)
       case (3)
          call make_deltaT_3d(tp(:,:,:,comp_dT),tp(:,:,:,comp_tfromH), &
                              tp(:,:,:,comp_tfromp), ng_p, lo, hi)
       end select
    end do

  end subroutine make_deltaT

  subroutine make_deltaT_2d(dT,tfromH,tfromp,ng_p,lo,hi)

    integer, intent(in) :: lo(:), hi(:), ng_p
    real (kind = dp_t), intent(  out) ::     dT(lo(1)-ng_p:,lo(2)-ng_p:)
    real (kind = dp_t), intent(in   ) :: tfromH(lo(1)-ng_p:,lo(2)-ng_p:)
    real (kind = dp_t), intent(in   ) :: tfromp(lo(1)-ng_p:,lo(2)-ng_p:)

    !     Local variables
    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          dT(i,j) = tfromH(i,j) - tfromp(i,j)
       end do
    end do

  end subroutine make_deltaT_2d

  subroutine make_deltaT_3d(dT,tfromH,tfromp,ng_p,lo,hi)

    integer, intent(in) :: lo(:), hi(:), ng_p
    real (kind = dp_t), intent(  out) ::     dT(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real (kind = dp_t), intent(in   ) :: tfromH(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real (kind = dp_t), intent(in   ) :: tfromp(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)

    !     Local variables
    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dT(i,j,k) = tfromH(i,j,k) - tfromp(i,j,k)
          end do
       end do
    end do

  end subroutine make_deltaT_3d

  subroutine make_divw0(w0_cart,divw0,dx)

    type(multifab) , intent(in   ) :: w0_cart
    type(multifab) , intent(inout) :: divw0
    real(kind=dp_t), intent(in   ) :: dx(:)

    real(kind=dp_t), pointer :: w0p(:,:,:,:)
    real(kind=dp_t), pointer :: dwp(:,:,:,:)
    integer                  :: lo(w0_cart%dim),hi(w0_cart%dim)
    integer                  :: i,dm,ng_w0,ng_dw

    dm = w0_cart%dim
    ng_w0 = w0_cart%ng
    ng_dw = divw0%ng

    do i=1,w0_cart%nboxes
       if ( multifab_remote(w0_cart, i) ) cycle
       w0p => dataptr(w0_cart, i)
       dwp => dataptr(divw0, i)
       lo = lwb(get_box(w0_cart, i))
       hi = upb(get_box(w0_cart, i))
       select case (dm)
       case (2)
          call make_divw0_2d(w0p(:,:,1,:), ng_w0, dwp(:,:,1,1), ng_dw, lo, hi, dx)
       case (3)
          call make_divw0_3d(w0p(:,:,:,:), ng_w0, dwp(:,:,:,1), ng_dw, lo, hi, dx)
       end select
    end do

  end subroutine make_divw0

  subroutine make_divw0_2d(w0_cart,ng_w0,divw0,ng_dw,lo,hi,dx)

    integer, intent(in)               :: lo(:), hi(:), ng_w0, ng_dw
    real (kind = dp_t), intent(in   ) :: w0_cart(lo(1)-ng_w0:,lo(2)-ng_w0:,:)
    real (kind = dp_t), intent(inout) ::   divw0(lo(1)-ng_dw:,lo(2)-ng_dw:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    !     Local variables
    integer :: i,j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          divw0(i,j) = (w0_cart(i+1,j,1)-w0_cart(i-1,j,1)) / (2.d0 * dx(1)) &
                      +(w0_cart(i,j+1,2)-w0_cart(i,j-1,2)) / (2.d0 * dx(2))

       end do
    end do

  end subroutine make_divw0_2d

  subroutine make_divw0_3d(w0_cart,ng_w0,divw0,ng_dw,lo,hi,dx)

    integer, intent(in)               :: lo(:), hi(:), ng_w0, ng_dw
    real (kind = dp_t), intent(in   ) :: w0_cart(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:,:)
    real (kind = dp_t), intent(inout) ::   divw0(lo(1)-ng_dw:,lo(2)-ng_dw:,lo(3)-ng_dw:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    !     Local variables
    integer :: i,j,k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             divw0(i,j,k) = (w0_cart(i+1,j,k,1)-w0_cart(i-1,j,k,1)) / (2.d0 * dx(1)) &
                           +(w0_cart(i,j+1,k,2)-w0_cart(i,j-1,k,2)) / (2.d0 * dx(2)) &
                           +(w0_cart(i,j,k+1,3)-w0_cart(i,j,k-1,3)) / (2.d0 * dx(3))
          end do
       end do
    end do

  end subroutine make_divw0_3d

end module plot_variables_module

