
module plot_variables_module

  use bl_types
  use multifab_module
  use eos_module
  use fill_3d_module
  use network
  use variables
  use geometry
  use probin_module

  implicit none

  private
  public :: make_enthalpy, make_tfromH, make_tfromrho, make_XfromrhoX, &
            make_omegadot, make_deltaT

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
    integer :: i, j, n
    real(kind=dp_t) :: qreact

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          qreact = 0.0d0
          if(use_big_h) then
             do n=1,nspec
                qreact = qreact + ebin(n)*s(i,j,spec_comp+n-1)/s(i,j,rho_comp)
             enddo
             enthalpy(i,j) = s(i,j,rhoh_comp)/s(i,j,rho_comp) - qreact
          else
             enthalpy(i,j) = s(i,j,rhoh_comp)/s(i,j,rho_comp)
          endif

       enddo
    enddo

  end subroutine make_enthalpy_2d

  subroutine make_enthalpy_3d (enthalpy,s,lo,hi,ng)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) :: enthalpy(lo(1):,lo(2):, lo(3):)  
    real (kind = dp_t), intent(in   ) ::    s(lo(1)-ng:,lo(2)-ng:, lo(3)-ng:,:)
    
!     Local variables
    integer :: i, j, k, n
    real(kind=dp_t) :: qreact
    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             qreact = 0.0d0
             if(use_big_h) then
                do n=1,nspec
                   qreact = qreact + ebin(n)*s(i,j,k,spec_comp+n-1)/s(i,j,k,rho_comp)
                enddo
                enthalpy(i,j,k) = s(i,j,k,rhoh_comp)/s(i,j,k,rho_comp) - qreact
             else
                enthalpy(i,j,k) = s(i,j,k,rhoh_comp)/s(i,j,k,rho_comp)
             endif

          enddo
       enddo
    end do
    
  end subroutine make_enthalpy_3d

  subroutine make_tfromH (plotdata,comp_t,comp_dp,state,p0,t0, dx)

    integer        , intent(in   ) :: comp_t,comp_dp
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: state
    real(kind=dp_t), intent(in   ) :: p0(0:)
    real(kind=dp_t), intent(in   ) :: t0(0:)
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
          call maketfromH_2d(tp(:,:,1,comp_t),tp(:,:,1,comp_dp),sp(:,:,1,:), lo, hi, ng, p0, t0)
       case (3)
          if (spherical .eq. 1) then
            call maketfromH_3d_sphr(tp(:,:,:,comp_t),tp(:,:,:,comp_dp),sp(:,:,:,:), &
                                    lo, hi, ng, p0, t0, dx)
          else
            call maketfromH_3d_cart(tp(:,:,:,comp_t),tp(:,:,:,comp_dp),sp(:,:,:,:), &
                                    lo, hi, ng, p0, t0)
          end if
       end select
    end do

  end subroutine make_tfromH

  subroutine maketfromH_2d (T,deltaP,state,lo,hi,ng,p0,t0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) ::      T(lo(1)   :,lo(2):)  
    real (kind = dp_t), intent(  out) :: deltaP(lo(1)   :,lo(2):)  
    real (kind = dp_t), intent(in   ) ::  state(lo(1)-ng:,lo(2)-ng:,:)
    real (kind = dp_t), intent(in   ) ::  p0(0:)
    real (kind = dp_t), intent(in   ) ::  t0(0:)

    !     Local variables
    integer :: i, j, n
    real(kind=dp_t) qreact

    do_diag = .false.

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          ! (rho, H) --> T, p
            
          den_eos(1)  = state(i,j,rho_comp)
          p_eos(1)    = p0(j)
          temp_eos(1) = t0(j)
          xn_eos(1,:) = state(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)

          qreact = 0.0d0
          if(use_big_h) then
             do n=1,nspec
                qreact = qreact + ebin(n)*xn_eos(1,n)
             enddo
             h_eos(1) = state(i,j,rhoh_comp) / state(i,j,rho_comp) - qreact
          else
             h_eos(1) = state(i,j,rhoh_comp) / state(i,j,rho_comp)
          endif

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
          
!         T(i,j) = log(temp_eos(1))/log(10.)
          T(i,j) = temp_eos(1)

          deltaP(i,j) = (p_eos(1)-p0(j))/ p0(j)
          
       enddo
    enddo

  end subroutine maketfromH_2d

  subroutine maketfromH_3d_cart (T,deltaP,state,lo,hi,ng,p0,t0)

    implicit none
    integer, intent(in) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(  out) ::      T(lo(1)   :,lo(2):   ,lo(3):     )  
    real (kind = dp_t), intent(  out) :: deltaP(lo(1)   :,lo(2):   ,lo(3):     )  
    real (kind = dp_t), intent(in   ) :: state(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind = dp_t), intent(in   ) :: p0(0:)
    real (kind = dp_t), intent(in   ) :: t0(0:)

    !     Local variables
    integer :: i, j, k, n
    real(kind=dp_t) qreact

    do_diag = .false.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! (rho, H) --> T, p
              den_eos(1)  = state(i,j,k,rho_comp)
             p_eos(1)    = p0(k)
             temp_eos(1) = t0(k)
             xn_eos(1,:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

             qreact = 0.0d0
             if(use_big_h) then
                do n=1,nspec
                   qreact = qreact + ebin(n)*xn_eos(1,n)
                enddo
                h_eos(1) = state(i,j,k,rhoh_comp)/state(i,j,k,rho_comp) - qreact
             else
                h_eos(1) = state(i,j,k,rhoh_comp)/state(i,j,k,rho_comp)
             endif

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

!            T(i,j,k) = log(temp_eos(1))/log(10.)
             T(i,j,k) = temp_eos(1)

             deltaP(i,j,k) = (p_eos(1)-p0(k))/ p0(k)

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
    real (kind = dp_t), intent(in   ) ::    p0(0:)
    real (kind = dp_t), intent(in   ) ::    t0(0:)
    real (kind = dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer :: i, j, k
    real (kind=dp_t), allocatable :: t0_cart(:,:,:)
    real (kind=dp_t), allocatable :: p0_cart(:,:,:)

    allocate(t0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(t0_cart,t0,lo,hi,dx,0)

    allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(p0_cart,p0,lo,hi,dx,0)

    do_diag = .false.

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
            
             den_eos(1)  = state(i,j,k,rho_comp)
             h_eos(1)    = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)
             if(use_big_h) then
               print*,"WARNING: H conversion not defined in maketfromH_3d_sphr"
             endif
             p_eos(1)    = p0_cart(i,j,k)
             temp_eos(1) = t0_cart(i,j,k)
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

!            T(i,j,k) = log(temp_eos(1))/log(10.)
             T(i,j,k) = temp_eos(1)

             deltaP(i,j,k) = (p_eos(1)-p0_cart(i,j,k))/ p0_cart(i,j,k)

          enddo
       enddo
    enddo
   
    deallocate(t0_cart,p0_cart)

  end subroutine maketfromH_3d_sphr

  subroutine make_tfromrho (plotdata,comp_tfromrho,comp_tpert,comp_rhopert, &
                            comp_machno,comp_deltag,comp_spert, &
                            s,u,s0,p0,dx)

    integer        , intent(in   ) :: comp_tfromrho,comp_tpert
    integer        , intent(in   ) :: comp_rhopert, comp_machno
    integer        , intent(in   ) :: comp_deltag, comp_spert
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: s
    type(multifab) , intent(in   ) :: u
    real(kind=dp_t), intent(in   ) :: s0(0:,:)
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
          call maketfromrho_2d(tp(:,:,1,comp_tfromrho),tp(:,:,1,comp_tpert), &
                               tp(:,:,1,comp_rhopert ), &
                               tp(:,:,1,comp_machno  ),tp(:,:,1,comp_deltag), &
                               tp(:,:,1,comp_spert   ), &
                               sp(:,:,1,:), up(:,:,1,:), &
                               lo, hi, ng, s0, p0)
       case (3)
          if (spherical .eq. 1) then
            call maketfromrho_3d_sphr(tp(:,:,:,comp_tfromrho),tp(:,:,:,comp_tpert), &
                                      tp(:,:,:,comp_rhopert ), &
                                      tp(:,:,:,comp_machno  ),tp(:,:,:,comp_deltag), &
                                      tp(:,:,:,comp_spert   ), &
                                      sp(:,:,:,:), up(:,:,:,:), &
                                      lo, hi, ng, s0, p0, dx)
          else
            call maketfromrho_3d_cart(tp(:,:,:,comp_tfromrho),tp(:,:,:,comp_tpert), &
                                      tp(:,:,:,comp_rhopert ), &
                                      tp(:,:,:,comp_machno  ),tp(:,:,:,comp_deltag), &
                                      tp(:,:,:,comp_spert   ), &
                                      sp(:,:,:,:), up(:,:,:,:), &
                                      lo, hi, ng, s0, p0)
          endif
       end select
    end do

  end subroutine make_tfromrho

  subroutine maketfromrho_2d (t,tpert,rhopert,machno,deltagamma,spert, &
                              s,u,lo,hi,ng,s0,p0)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(  out)     ::      t(lo(1):,lo(2):)  
    real (kind=dp_t), intent(  out) ::      tpert(lo(1):,lo(2):)  
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1):,lo(2):)  
    real (kind=dp_t), intent(  out) ::     machno(lo(1):,lo(2):)  
    real (kind=dp_t), intent(  out) :: deltagamma(lo(1):,lo(2):)  
    real (kind=dp_t), intent(  out) ::      spert(lo(1):,lo(2):)  
    real (kind=dp_t), intent(in   ) ::  s(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(in   ) ::  u(lo(1)-ng:,lo(2)-ng:,:)
    real (kind=dp_t), intent(in   ) :: s0(0:,:)
    real (kind=dp_t), intent(in   ) :: p0(0:)

    !     Local variables
    integer          :: i, j
    real (kind=dp_t) :: vel
!   real (kind=dp_t), allocatable :: gam10(:), entr0(:)

!   allocate(gam10(lo(2):hi(2)))
!   allocate(entr0(lo(2):hi(2)))

    do_diag = .false.

    ! We now assume that the temperature coming in in the base state is correct, but
    !   we do this eos call to get gam10 and entr0.
!   do j = lo(2), hi(2)
!       den_eos(1) = s0(j,rho_comp)
!      temp_eos(1) = s0(j,temp_comp)
!         p_eos(1) = p0(j)
!      xn_eos(1,:) = s0(j,spec_comp:spec_comp+nspec-1)/den_eos(1)

!      ! (rho,P) --> T,h

!      call eos(eos_input_rp, den_eos, temp_eos, &
!               npts, nspec, &
!               xn_eos, &
!               p_eos, h_eos, e_eos, & 
!               cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
!               dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
!               dpdX_eos, dhdX_eos, &
!               gam1_eos, cs_eos, s_eos, &
!               dsdt_eos, dsdr_eos, &
!               do_diag)
!    gam10(j) = gam1_eos(1)
!    entr0(j) = s_eos(1)

!   end do

    ! Then compute the perturbation
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          den_eos(1)  = s(i,j,rho_comp)
          temp_eos(1) = s0(j,temp_comp)
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
          tpert(i,j) = temp_eos(1) - s0(j,temp_comp)

          rhopert(i,j) = s(i,j,rho_comp) - s0(j,rho_comp)

          vel = sqrt(u(i,j,1)*u(i,j,1) + u(i,j,2)*u(i,j,2))
          machno(i,j) = vel / cs_eos(1)

!         deltagamma(i,j) = gam1_eos(1) - gam10(j)
!         spert(i,j) = s_eos(1) - entr0(j)

       enddo
    enddo

    deltagamma = 0.d0
         spert = 0.d0

!   deallocate(gam10)
!   deallocate(entr0)

  end subroutine maketfromrho_2d

  subroutine maketfromrho_3d_cart (t,tpert,rhopert,machno,deltagamma,spert, &
                                   s,u,lo,hi,ng,s0,p0)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(  out) ::          t(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::      tpert(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::     machno(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) :: deltagamma(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::      spert(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(in   ) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(in   ) ::  u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(in   ) :: s0(0:,:)
    real (kind=dp_t), intent(in   ) :: p0(0:)

    !     Local variables
    integer          :: i, j, k
    real (kind=dp_t) :: vel
!    real (kind=dp_t), allocatable :: gam10(:), entr0(:)

!    allocate(gam10(lo(3):hi(3)))
!    allocate(entr0(lo(3):hi(3)))

    do_diag = .false.

    ! We now assume that the temperature coming in in the base state is correct, but
    !   we do this eos call to get gam10 and entr0.
!    do k = lo(3), hi(3)
!        den_eos(1) = s0(k,rho_comp)
!       temp_eos(1) = s0(k,temp_comp)
!          p_eos(1) = p0(k)
!       xn_eos(1,:) = s0(k,spec_comp:spec_comp+nspec-1)/den_eos(1)

       ! (rho,P) --> T,h
!       call eos(eos_input_rp, den_eos, temp_eos, &
!                npts, nspec, &
!                xn_eos, &
!                p_eos, h_eos, e_eos, & 
!                cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
!                dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
!                dpdX_eos, dhdX_eos, &
!                gam1_eos, cs_eos, s_eos, &
!                dsdt_eos, dsdr_eos, &
!                do_diag)
!     gam10(k) = gam1_eos(1)
!     entr0(k) = s_eos(1)
!    end do

    ! Then compute the perturbation and Mach number
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_eos(1) = s(i,j,k,rho_comp)
             temp_eos(1) = s0(k,temp_comp)
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

!            t(i,j,k) = log(temp_eos(1))/log(10.)
             t(i,j,k) = temp_eos(1)
             tpert(i,j,k) = temp_eos(1) - s0(k,temp_comp)

             rhopert(i,j,k) = s(i,j,k,rho_comp) - s0(k,rho_comp)

             vel = sqrt(u(i,j,k,1)*u(i,j,k,1) + u(i,j,k,2)*u(i,j,k,2) + u(i,j,k,3)*u(i,j,k,3))
             machno(i,j,k) = vel / cs_eos(1)

!             deltagamma(i,j,k) = gam1_eos(1) - gam10(k)
!             spert(i,j,k) = s_eos(1) - entr0(k)
          enddo
       enddo
    enddo

    deltagamma = 0.d0
    spert = 0.d0

!    deallocate(gam10)

   end subroutine maketfromrho_3d_cart

  subroutine maketfromrho_3d_sphr (t,tpert,rhopert,machno,deltagamma,spert, &
                                   s,u,lo,hi,ng,s0,p0,dx)

    implicit none

    integer, intent(in) :: lo(:), hi(:), ng
    real (kind=dp_t), intent(  out) ::          t(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::      tpert(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::     machno(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) :: deltagamma(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(  out) ::      spert(lo(1):,lo(2):,lo(3):)  
    real (kind=dp_t), intent(in   ) ::  s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(in   ) ::  u(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real (kind=dp_t), intent(in   ) :: s0(0:,:)
    real (kind=dp_t), intent(in   ) :: p0(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer          :: i, j, k
    real (kind=dp_t) :: vel
    real (kind=dp_t), allocatable :: gam10(:), entr0(:)
    real (kind=dp_t), allocatable ::  rho0_cart(:,:,:)
    real (kind=dp_t), allocatable ::    t0_cart(:,:,:)
    real (kind=dp_t), allocatable ::    p0_cart(:,:,:)
    real (kind=dp_t), allocatable ::  gam0_cart(:,:,:)
    real (kind=dp_t), allocatable :: entr0_cart(:,:,:)

    integer :: nr 

    nr = size(s0,dim=1)

    allocate(gam10(0:nr-1))
    allocate(entr0(0:nr-1))

    do_diag = .false.

    ! We now assume that the temperature coming in in the base state is correct, but
    !   we do this eos call to get gam10 and entr0.
    do k = 0, nr-1
        den_eos(1) = s0(k,rho_comp)
       temp_eos(1) = s0(k,temp_comp)
          p_eos(1) = p0(k)
        xn_eos(1,:) = s0(k,spec_comp:spec_comp+nspec-1)/den_eos(1)

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
     gam10(k) = gam1_eos(1)
     entr0(k) = s_eos(1)
    end do

    allocate(rho0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(rho0_cart,s0(:,rho_comp),lo,hi,dx,0)

    allocate(t0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(t0_cart,s0(:,temp_comp),lo,hi,dx,0)

    allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(p0_cart,p0,lo,hi,dx,0)

    allocate(gam0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(gam0_cart,gam10,lo,hi,dx,0)

    allocate(entr0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    call fill_3d_data(entr0_cart,entr0,lo,hi,dx,0)

    ! Then compute the perturbation and Mach number
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

              den_eos(1) = s(i,j,k,rho_comp)
             temp_eos(1) = t0_cart(i,j,k)
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

!            t(i,j,k) = log(temp_eos(1))/log(10.)
             t(i,j,k) = temp_eos(1)
             tpert(i,j,k) = temp_eos(1) - t0_cart(i,j,k)

             rhopert(i,j,k) = s(i,j,k,rho_comp) - rho0_cart(i,j,k)

             vel = sqrt(u(i,j,k,1)*u(i,j,k,1) + u(i,j,k,2)*u(i,j,k,2) + u(i,j,k,3)*u(i,j,k,3))
             machno(i,j,k) = vel / cs_eos(1)

             deltagamma(i,j,k) = gam1_eos(1) - gam0_cart(i,j,k)
             
             spert(i,j,k) = s_eos(1) - entr0_cart(i,j,k)
          enddo
       enddo
    enddo

    deallocate(gam10,entr0)
    deallocate(rho0_cart,t0_cart,p0_cart,gam0_cart,entr0_cart)

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


   subroutine make_omegadot(plotdata,comp,comp_enuc,s,rho_omegadot)

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
           call makeomega_2d(pp(:,:,1,comp:comp-1+nspec),pp(:,:,1,comp_enuc), &
                             rop(:,:,1,:),sp(:,:,1,rho_comp), lo, hi, ng_s, ng_o)
        case (3)
           call makeomega_3d(pp(:,:,:,comp:comp-1+nspec),pp(:,:,:,comp_enuc), &
                             rop(:,:,:,:),sp(:,:,:,rho_comp), lo, hi, ng_s, ng_o)
        end select
     end do
     
   end subroutine make_omegadot


   subroutine makeomega_2d (omegadot,enuc,rho_omegadot,rho,lo,hi,ng_s,ng_o)

     implicit none
     
     integer, intent(in) :: lo(:), hi(:), ng_s, ng_o
     real (kind = dp_t), intent(  out) :: omegadot(lo(1):,lo(2):,:)
     real (kind = dp_t), intent(  out) :: enuc(lo(1):,lo(2):)
     real (kind = dp_t), intent(in   ) :: rho_omegadot(lo(1)-ng_o:,lo(2)-ng_o:,:)
     real (kind = dp_t), intent(in   ) :: rho(lo(1)-ng_s:,lo(2)-ng_s:)
     
     ! Local variables
     integer :: i, j, n

     enuc(:,:) = ZERO

     do n = 1, nspec
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              omegadot(i,j,n) = rho_omegadot(i,j,n) / rho(i,j)
              enuc(i,j) = enuc(i,j) - omegadot(i,j,n)*ebin(n)
           enddo
        enddo
     enddo
     
   end subroutine makeomega_2d
   
   subroutine makeomega_3d (omegadot,enuc,rho_omegadot,rho,lo,hi,ng_s,ng_o)

     implicit none
     
     integer, intent(in) :: lo(:), hi(:), ng_s, ng_o
     real (kind = dp_t), intent(  out) :: omegadot(lo(1):,lo(2):,lo(3):,:)  
     real (kind = dp_t), intent(  out) :: enuc(lo(1):,lo(2):,lo(3):)
     real (kind = dp_t), intent(in   ) :: rho_omegadot(lo(1)-ng_o:,lo(2)-ng_o:,lo(3)-ng_o:,:)
     real (kind = dp_t), intent(in   ) :: rho(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:  )

     ! Local variables
     integer :: i, j, k, n

     enuc(:,:,:) = ZERO

     do n = 1, nspec
        do k = lo(3), hi(3)
           do j = lo(2), hi(2)
              do i = lo(1), hi(1)
                 omegadot(i,j,k,n) = rho_omegadot(i,j,k,n) / rho(i,j,k)
                 enuc(i,j,k) = enuc(i,j,k) - omegadot(i,j,k,n)*ebin(n)
              enddo
           enddo
        enddo
     enddo
     
   end subroutine makeomega_3d

  subroutine make_deltaT (plotdata,comp_dT,comp_tfromH,comp_tfromrho)

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
          call makedeltaT_2d(tp(:,:,1,comp_dT),tp(:,:,1,comp_tfromH), tp(:,:,1,comp_tfromrho), lo, hi)
       case (3)
          call makedeltaT_3d(tp(:,:,:,comp_dT),tp(:,:,:,comp_tfromH), tp(:,:,:,comp_tfromrho), lo, hi)
       end select
    end do

  end subroutine make_deltaT

  subroutine makedeltaT_2d (dT,tfromH,tfromrho,lo,hi)

    implicit none
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

  end subroutine makedeltaT_2d

  subroutine makedeltaT_3d (dT,tfromH,tfromrho,lo,hi)

    implicit none
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

  end subroutine makedeltaT_3d
   
 end module plot_variables_module
