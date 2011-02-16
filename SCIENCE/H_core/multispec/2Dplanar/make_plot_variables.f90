module plot_variables_module

  use bl_types
  use multifab_module
  use define_bc_module
  use probin_module, only: use_tfromp

  implicit none

  private

  public :: make_abar_change, make_stability, make_ad_excess
  public :: make_conductivity, make_pi_cc
  public :: make_tfromH, make_tfromp, make_entropypert
  public :: make_deltaT, make_divw0, make_vorticity, make_magvel, make_velrc
  public :: make_rhopert, make_rhohpert

contains


  !---------------------------------------------------------------------------
  ! check parcel stability
  !---------------------------------------------------------------------------
  subroutine make_abar_change(plotdata,comp_abar_change,state,normal)

    use geometry, only: spherical

    type(multifab), intent(inout) :: plotdata
    integer,        intent(in   ) :: comp_abar_change
    type(multifab), intent(in   ) :: state
    type(multifab), intent(in   ) :: normal
    
    real(kind=dp_t), pointer :: sp(:,:,:,:), cp(:,:,:,:), np(:,:,:,:)
    real(kind=dp_t), pointer :: ep(:,:,:,:)
    integer :: lo(get_dim(plotdata)), hi(get_dim(plotdata)), ng_s, ng_c, ng_n
    integer :: i, dm

    dm = get_dim(plotdata)

    ng_c = nghost(plotdata)
    ng_s = nghost(state)
    ng_n = nghost(normal)

    do i = 1, nboxes(state)
       if (multifab_remote(state, i)) cycle
       sp => dataptr(state, i)
       cp => dataptr(plotdata, i)
       lo = lwb(get_box(state, i))
       hi = upb(get_box(state, i))
       select case (dm)
       case (1)
          call bl_error('abar/abar_initial not written for 3d')
       case (2)
          call make_abar_change_2d(cp(:,:,1,comp_abar_change), ng_c, &
                                 sp(:,:,1,:), ng_s, &
                                 lo, hi)
       case (3)
          call bl_error('abar/abar_initial not written for 3d')
      end select

    enddo
  end subroutine make_abar_change
  
  subroutine make_abar_change_2d(abar_change, ng_ac, state, ng_s, lo, hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use model_parser_module
    use network, only: nspec, aion
    use probin_module, only: base_cutoff_density
    use bl_constants_module

    integer,         intent(in   ) :: lo(:), hi(:), ng_ac, ng_s
    real(kind=dp_t), intent(  out) :: abar_change(lo(1)-ng_ac:,lo(2)-ng_ac:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,:)

    real(kind=dp_t) :: abar_local

    integer :: i, j, n
    

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          abar_local = 0.d0

          do n=1,nspec
             abar_local = abar_local + state(i,j,spec_comp+n-1) / &
                  (state(i,j,rho_comp) * aion(n))
          enddo
       
          abar_local = 1.0d0/abar_local
    
          ! meanA indexing starts at 1 
          abar_change(i,j) = (meanA(j) - abar_local)/meanA(j) 

       enddo

    enddo

  end subroutine make_abar_change_2d

  !---------------------------------------------------------------------------
  ! check parcel stability
  !---------------------------------------------------------------------------
  subroutine make_stability(plotdata,comp_stability,state,normal)

    use geometry, only: spherical

    type(multifab), intent(inout) :: plotdata
    integer,        intent(in   ) :: comp_stability
    type(multifab), intent(in   ) :: state
    type(multifab), intent(in   ) :: normal
    
    real(kind=dp_t), pointer :: sp(:,:,:,:), cp(:,:,:,:), np(:,:,:,:)
    real(kind=dp_t), pointer :: ep(:,:,:,:)
    integer :: lo(get_dim(plotdata)), hi(get_dim(plotdata)), ng_s, ng_c, ng_n
    integer :: i, dm

    dm = get_dim(plotdata)

    ng_c = nghost(plotdata)
    ng_s = nghost(state)
    ng_n = nghost(normal)

    do i = 1, nboxes(state)
       if (multifab_remote(state, i)) cycle
       sp => dataptr(state, i)
       cp => dataptr(plotdata, i)
       lo = lwb(get_box(state, i))
       hi = upb(get_box(state, i))
       select case (dm)
       case (1)
          call make_stability_1d(cp(:,1,1,comp_stability), ng_c, &
                                 sp(:,1,1,:), ng_s, &
                                 lo, hi)
       case (2)
          call make_stability_2d(cp(:,:,1,comp_stability), ng_c, &
                                 sp(:,:,1,:), ng_s, &
                                 lo, hi)
       case (3)
          if (spherical .eq. 1) then
             np => dataptr(normal, i)       
             call make_stability_3d_sphr(cp(:,:,:,comp_stability), ng_c, &
                                         sp(:,:,:,:), ng_s, &
                                         np(:,:,:,:), ng_n, lo, hi)
          else
             call make_stability_3d(cp(:,:,:,comp_stability), ng_c, &
                                    sp(:,:,:,:), ng_s, &
                                    lo, hi)
          endif
       end select

    enddo
  end subroutine make_stability
  
  subroutine make_stability_1d(stability, ng_ad, state, ng_s, lo, hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use network, only: nspec
    use probin_module, only: base_cutoff_density
    use bl_constants_module

    integer,         intent(in   ) :: lo(:), hi(:), ng_ad, ng_s
    real(kind=dp_t), intent(  out) :: stability(lo(1)-ng_ad:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,:)

    real(kind=dp_t) :: pres(lo(1):hi(1)), nabla_ad(lo(1):hi(1))
    real(kind=dp_t) :: chi_rho, chi_t, dt, dp, nabla

    integer :: i
    
    call bl_error("make_plot_variables.f90:: parcel stability not written for 1D")

! use this a base to write the function

    do i = lo(1), hi(1)

       den_eos(1) = state(i,rho_comp)
       temp_eos(1) = state(i,temp_comp)
       xn_eos(1,:) = state(i,spec_comp:spec_comp+nspec-1)/den_eos(1)

       pt_index_eos(:) = (/i, -1, -1/)       

       call eos(eos_input_rt, den_eos, temp_eos, &
                npts, &
                xn_eos, &
                p_eos, h_eos, e_eos, &
                cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                dpdX_eos, dhdX_eos, &
                gam1_eos, cs_eos, s_eos, &
                dsdt_eos, dsdr_eos, &
                .false., &
                pt_index_eos)       
       
       pres(i) = p_eos(1)

       chi_rho = den_eos(1) * dpdr_eos(1) / p_eos(1)
       chi_t = temp_eos(1) * dpdt_eos(1) / p_eos(1)
       nabla_ad(i) = (gam1_eos(1) - chi_rho) / (chi_t * gam1_eos(1))

    enddo

    do i = lo(1), hi(1)
       if (state(i,rho_comp) <= base_cutoff_density) then
          nabla = ZERO
       else
          ! forward difference
          if (i == lo(1)) then
             dt = state(i+1,temp_comp) - state(i,temp_comp)
             dp = pres(i+1) - pres(i)
             ! backward difference
          else if (i == hi(1)) then
             dt = state(i,temp_comp) - state(i-1,temp_comp)
             dp = pres(i) - pres(i-1)
             ! centered difference
          else
             dt = state(i+1,temp_comp) - state(i-1,temp_comp)
             dp = pres(i+1) - pres(i-1)
          endif

          nabla = pres(i) * dt / (dp * state(i,temp_comp))
       endif

       stability(i) = nabla - nabla_ad(i)
    enddo

  end subroutine make_stability_1d

  subroutine make_stability_2d(stability, ng_ad, state, ng_s, lo, hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use network, only: nspec
    use probin_module, only: base_cutoff_density
    use bl_constants_module

    integer,         intent(in   ) :: lo(:), hi(:), ng_ad, ng_s
    real(kind=dp_t), intent(  out) :: stability(lo(1)-ng_ad:,lo(2)-ng_ad:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,:)

    real(kind=dp_t) :: pres(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)
    real(kind=dp_t) :: entr(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)

    integer :: i, j
    

    do j = lo(2), hi(2)+1
       do i = lo(1), hi(1)

          !get entropy and pressure given rho and T
          den_eos(1) = state(i,j,rho_comp)
          temp_eos(1) = state(i,j,temp_comp)
          xn_eos(1,:) = state(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)

          pt_index_eos(:) = (/i, j, -1/)       

          call eos(eos_input_rt, den_eos, temp_eos, &
                   npts, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, &
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   .false., &
                   pt_index_eos)       


          pres(i,j) = p_eos(1)
          entr(i,j) = s_eos(1)

       enddo
    enddo

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          ! don't know if really need to fill these, but put a reasonable 
          !   guess in here
          den_eos(1) = state(i,j,rho_comp)
          temp_eos(1) = state(i,j,temp_comp)

          ! retain composition and entropy at i
          xn_eos(1,:) = state(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)
          s_eos(1) = entr(i,j)
          ! enforce pressure equilibration at j+1
          p_eos(1) = pres(i,j+1)

          pt_index_eos(:) = (/i, j, -1/)       
         
          call eos(eos_input_ps, den_eos, temp_eos, &
                   npts, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, &
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   .false., pt_index_eos)

          ! compare density of risen parcel to local density
          stability(i,j) = (den_eos(1) - state(i,j+1,rho_comp))

       enddo
    enddo

  end subroutine make_stability_2d

  subroutine make_stability_3d(stability, ng_ad, state, ng_s, lo, hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use network, only: nspec
    use probin_module, only: base_cutoff_density
    use bl_constants_module

    integer,         intent(in   ) :: lo(:), hi(:), ng_ad, ng_s
    real(kind=dp_t), intent(  out) :: stability(lo(1)-ng_ad:,lo(2)-ng_ad:,lo(3)-ng_ad:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)

    real(kind=dp_t) :: pres(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: nabla_ad(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: chi_rho, chi_t, dt, dp, nabla

    integer :: i, j, k

    call bl_error("make_plot_variables.f90:: parcel stability not written for 3D")


! use this a base to write the function

    !$OMP PARALLEL DO PRIVATE(i,j,k,chi_rho,chi_t)    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_eos(1) = state(i,j,k,rho_comp)
             temp_eos(1) = state(i,j,k,temp_comp)
             xn_eos(1,:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

             pt_index_eos(:) = (/i, j, k/)       

             call eos(eos_input_rt, den_eos, temp_eos, &
                      npts, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      .false., &
                      pt_index_eos)       
       
             pres(i,j,k) = p_eos(1)

             chi_rho = den_eos(1) * dpdr_eos(1) / p_eos(1)
             chi_t = temp_eos(1) * dpdt_eos(1) / p_eos(1)
             nabla_ad(i,j,k) = (gam1_eos(1) - chi_rho) / (chi_t * gam1_eos(1))

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,dt,dp,nabla)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,rho_comp) <= base_cutoff_density) then
                nabla = ZERO
             else
                ! forward difference
                if (k == lo(3)) then
                   dt = state(i,j,k+1,temp_comp) - state(i,j,k,temp_comp)
                   dp = pres(i,j,k+1) - pres(i,j,k)
                   ! backward difference
                else if (k == hi(3)) then
                   dt = state(i,j,k,temp_comp) - state(i,j,k-1,temp_comp)
                   dp = pres(i,j,k) - pres(i,j,k-1)
                   ! centered difference
                else
                   dt = state(i,j,k+1,temp_comp) - state(i,j,k-1,temp_comp)
                   dp = pres(i,j,k+1) - pres(i,j,k-1)
                endif

                nabla = pres(i,j,k) * dt / (dp * state(i,j,k,temp_comp))
             endif

             stability(i,j,k) = nabla - nabla_ad(i,j,k)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine make_stability_3d

  subroutine make_stability_3d_sphr(stability, ng_ad, state, ng_s, &
                                    normal, ng_n, lo, hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use network, only: nspec
    use probin_module, only: base_cutoff_density
    use bl_constants_module

    integer,         intent(in   ) :: lo(:), hi(:), ng_ad, ng_s, ng_n
    real(kind=dp_t), intent(  out) :: stability(lo(1)-ng_ad:,lo(2)-ng_ad:,lo(3)-ng_ad:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind = dp_t), intent(in   ) :: normal(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)  

    real(kind=dp_t) :: pres(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: nabla_ad(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: chi_rho, chi_t, nabla
    real(kind=dp_t) :: dp(4), dt(4)

    integer :: i, j, k, c

    call bl_error("make_plot_variables.f90:: parcel stability not written for spherical")

! use this a base to write the function

    !$OMP PARALLEL DO PRIVATE(i,j,k,chi_rho,chi_t)    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_eos(1) = state(i,j,k,rho_comp)
             temp_eos(1) = state(i,j,k,temp_comp)
             xn_eos(1,:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

             pt_index_eos(:) = (/i, j, k/)       

             call eos(eos_input_rt, den_eos, temp_eos, &
                      npts, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      .false., &
                      pt_index_eos)       
       
             pres(i,j,k) = p_eos(1)

             chi_rho = den_eos(1) * dpdr_eos(1) / p_eos(1)
             chi_t = temp_eos(1) * dpdt_eos(1) / p_eos(1)
             nabla_ad(i,j,k) = (gam1_eos(1) - chi_rho) / (chi_t * gam1_eos(1))

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,dt,dp,nabla)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,rho_comp) <= base_cutoff_density) then
                nabla = ZERO
             else
                ! compute gradient

                ! forward difference
                if (k == lo(3)) then
                   dt(3) = state(i,j,k+1,temp_comp) - state(i,j,k,temp_comp)
                   dp(3) = pres(i,j,k+1) - pres(i,j,k)
                   ! backward difference
                else if (k == hi(3)) then
                   dt(3) = state(i,j,k,temp_comp) - state(i,j,k-1,temp_comp)
                   dp(3) = pres(i,j,k) - pres(i,j,k-1)
                   ! centered difference
                else
                   dt(3) = state(i,j,k+1,temp_comp) - state(i,j,k-1,temp_comp)
                   dp(3) = pres(i,j,k+1) - pres(i,j,k-1)
                endif

                if (j == lo(2)) then
                   dt(2) = state(i,j+1,k,temp_comp) - state(i,j,k,temp_comp)
                   dp(2) = pres(i,j+1,k) - pres(i,j,k)
                   ! backward difference
                else if (j == hi(2)) then
                   dt(2) = state(i,j,k,temp_comp) - state(i,j-1,k,temp_comp)
                   dp(2) = pres(i,j,k) - pres(i,j-1,k)
                   ! centered difference
                else
                   dt(2) = state(i,j+1,k,temp_comp) - state(i,j-1,k,temp_comp)
                   dp(2) = pres(i,j+1,k) - pres(i,j-1,k)
                endif

                if (i == lo(1)) then
                   dt(1) = state(i+1,j,k,temp_comp) - state(i,j,k,temp_comp)
                   dp(1) = pres(i+1,j,k) - pres(i,j,k)
                   ! backward difference
                else if (i == hi(1)) then
                   dt(1) = state(i,j,k,temp_comp) - state(i-1,j,k,temp_comp)
                   dp(1) = pres(i,j,k) - pres(i-1,j,k)
                   ! centered difference
                else
                   dt(1) = state(i+1,j,k,temp_comp) - state(i-1,j,k,temp_comp)
                   dp(1) = pres(i+1,j,k) - pres(i-1,j,k)
                endif

                ! dot into normal to get d/dr
                dp(4) = 0.d0
                dt(4) = 0.d0
                do c = 1,3
                   dp(4) = dp(4) + dp(c)*normal(i,j,k,c) 
                   dt(4) = dt(4) + dt(c)*normal(i,j,k,c) 
                enddo

                nabla = pres(i,j,k)*dt(4) / (dp(4)*state(i,j,k,temp_comp))
             endif

             stability(i,j,k) = (nabla - nabla_ad(i,j,k)) / nabla
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine make_stability_3d_sphr

  !---------------------------------------------------------------------------
  ! make_ad_excess
  !---------------------------------------------------------------------------
  subroutine make_ad_excess(plotdata,comp_ad_excess,state,normal)

    use geometry, only: spherical

    type(multifab), intent(inout) :: plotdata
    integer,        intent(in   ) :: comp_ad_excess
    type(multifab), intent(in   ) :: state
    type(multifab), intent(in   ) :: normal
    
    real(kind=dp_t), pointer :: sp(:,:,:,:), cp(:,:,:,:), np(:,:,:,:)
    integer :: lo(get_dim(plotdata)), hi(get_dim(plotdata)), ng_s, ng_c, ng_n
    integer :: i, dm

    dm = get_dim(plotdata)

    ng_c = nghost(plotdata)
    ng_s = nghost(state)
    ng_n = nghost(normal)

    do i = 1, nboxes(state)
       if (multifab_remote(state, i)) cycle
       sp => dataptr(state, i)
       cp => dataptr(plotdata, i)
       lo = lwb(get_box(state, i))
       hi = upb(get_box(state, i))
       select case (dm)
       case (1)
          call make_ad_excess_1d(cp(:,1,1,comp_ad_excess), ng_c, &
                                 sp(:,1,1,:), ng_s, &
                                 lo, hi)
       case (2)
          call make_ad_excess_2d(cp(:,:,1,comp_ad_excess), ng_c, &
                                 sp(:,:,1,:), ng_s, &
                                 lo, hi)
       case (3)
          if (spherical .eq. 1) then
             np => dataptr(normal, i)       
             call make_ad_excess_3d_sphr(cp(:,:,:,comp_ad_excess), ng_c, &
                                         sp(:,:,:,:), ng_s, &
                                         np(:,:,:,:), ng_n, lo, hi)
          else
             call make_ad_excess_3d(cp(:,:,:,comp_ad_excess), ng_c, &
                                    sp(:,:,:,:), ng_s, &
                                    lo, hi)
          endif
       end select

    enddo
  end subroutine make_ad_excess
  
  subroutine make_ad_excess_1d(ad_excess, ng_ad, state, ng_s, lo, hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use network, only: nspec
    use probin_module, only: base_cutoff_density
    use bl_constants_module

    integer,         intent(in   ) :: lo(:), hi(:), ng_ad, ng_s
    real(kind=dp_t), intent(  out) :: ad_excess(lo(1)-ng_ad:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,:)

    real(kind=dp_t) :: pres(lo(1):hi(1)), nabla_ad(lo(1):hi(1))
    real(kind=dp_t) :: chi_rho, chi_t, dt, dp, nabla

    integer :: i
    
    do i = lo(1), hi(1)

       den_eos(1) = state(i,rho_comp)
       temp_eos(1) = state(i,temp_comp)
       xn_eos(1,:) = state(i,spec_comp:spec_comp+nspec-1)/den_eos(1)

       pt_index_eos(:) = (/i, -1, -1/)       

       call eos(eos_input_rt, den_eos, temp_eos, &
                npts, &
                xn_eos, &
                p_eos, h_eos, e_eos, &
                cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                dpdX_eos, dhdX_eos, &
                gam1_eos, cs_eos, s_eos, &
                dsdt_eos, dsdr_eos, &
                .false., &
                pt_index_eos)       
       
       pres(i) = p_eos(1)

       chi_rho = den_eos(1) * dpdr_eos(1) / p_eos(1)
       chi_t = temp_eos(1) * dpdt_eos(1) / p_eos(1)
       nabla_ad(i) = (gam1_eos(1) - chi_rho) / (chi_t * gam1_eos(1))

    enddo

    do i = lo(1), hi(1)
       if (state(i,rho_comp) <= base_cutoff_density) then
          nabla = ZERO
       else
          ! forward difference
          if (i == lo(1)) then
             dt = state(i+1,temp_comp) - state(i,temp_comp)
             dp = pres(i+1) - pres(i)
             ! backward difference
          else if (i == hi(1)) then
             dt = state(i,temp_comp) - state(i-1,temp_comp)
             dp = pres(i) - pres(i-1)
             ! centered difference
          else
             dt = state(i+1,temp_comp) - state(i-1,temp_comp)
             dp = pres(i+1) - pres(i-1)
          endif

          nabla = pres(i) * dt / (dp * state(i,temp_comp))
       endif

       ad_excess(i) = nabla - nabla_ad(i)
    enddo

  end subroutine make_ad_excess_1d

  subroutine make_ad_excess_2d(ad_excess, ng_ad, state, ng_s, lo, hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use network, only: nspec
    use probin_module, only: base_cutoff_density
    use bl_constants_module

    integer,         intent(in   ) :: lo(:), hi(:), ng_ad, ng_s
    real(kind=dp_t), intent(  out) :: ad_excess(lo(1)-ng_ad:,lo(2)-ng_ad:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,:)

    real(kind=dp_t) :: pres(lo(1):hi(1),lo(2):hi(2))
    real(kind=dp_t) :: nabla_ad(lo(1):hi(1),lo(2):hi(2))
    real(kind=dp_t) :: chi_rho, chi_t, dt, dp, nabla

    integer :: i, j
    
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          den_eos(1) = state(i,j,rho_comp)
          temp_eos(1) = state(i,j,temp_comp)
          xn_eos(1,:) = state(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)

          pt_index_eos(:) = (/i, j, -1/)       

          call eos(eos_input_rt, den_eos, temp_eos, &
                   npts, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, &
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   .false., &
                   pt_index_eos)       
       
          pres(i,j) = p_eos(1)

          chi_rho = den_eos(1) * dpdr_eos(1) / p_eos(1)
          chi_t = temp_eos(1) * dpdt_eos(1) / p_eos(1)
          nabla_ad(i,j) = (gam1_eos(1) - chi_rho) / (chi_t * gam1_eos(1))

       enddo
    enddo

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          if (state(i,j,rho_comp) <= base_cutoff_density) then
             nabla = ZERO
          else
             ! forward difference
             if (j == lo(2)) then
                dt = state(i,j+1,temp_comp) - state(i,j,temp_comp)
                dp = pres(i,j+1) - pres(i,j)
                ! backward difference
             else if (j == hi(2)) then
                dt = state(i,j,temp_comp) - state(i,j-1,temp_comp)
                dp = pres(i,j) - pres(i,j-1)
                ! centered difference
             else
                dt = state(i,j+1,temp_comp) - state(i,j-1,temp_comp)
                dp = pres(i,j+1) - pres(i,j-1)
             endif

             nabla = pres(i,j) * dt / (dp * state(i,j,temp_comp))
          endif

          ad_excess(i,j) = (nabla - nabla_ad(i,j))/abs(nabla_ad(i,j))
       enddo
    enddo

  end subroutine make_ad_excess_2d

  subroutine make_ad_excess_3d(ad_excess, ng_ad, state, ng_s, lo, hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use network, only: nspec
    use probin_module, only: base_cutoff_density
    use bl_constants_module

    integer,         intent(in   ) :: lo(:), hi(:), ng_ad, ng_s
    real(kind=dp_t), intent(  out) :: ad_excess(lo(1)-ng_ad:,lo(2)-ng_ad:,lo(3)-ng_ad:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)

    real(kind=dp_t) :: pres(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: nabla_ad(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: chi_rho, chi_t, dt, dp, nabla

    integer :: i, j, k

    !$OMP PARALLEL DO PRIVATE(i,j,k,chi_rho,chi_t)    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_eos(1) = state(i,j,k,rho_comp)
             temp_eos(1) = state(i,j,k,temp_comp)
             xn_eos(1,:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

             pt_index_eos(:) = (/i, j, k/)       

             call eos(eos_input_rt, den_eos, temp_eos, &
                      npts, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      .false., &
                      pt_index_eos)       
       
             pres(i,j,k) = p_eos(1)

             chi_rho = den_eos(1) * dpdr_eos(1) / p_eos(1)
             chi_t = temp_eos(1) * dpdt_eos(1) / p_eos(1)
             nabla_ad(i,j,k) = (gam1_eos(1) - chi_rho) / (chi_t * gam1_eos(1))

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,dt,dp,nabla)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,rho_comp) <= base_cutoff_density) then
                nabla = ZERO
             else
                ! forward difference
                if (k == lo(3)) then
                   dt = state(i,j,k+1,temp_comp) - state(i,j,k,temp_comp)
                   dp = pres(i,j,k+1) - pres(i,j,k)
                   ! backward difference
                else if (k == hi(3)) then
                   dt = state(i,j,k,temp_comp) - state(i,j,k-1,temp_comp)
                   dp = pres(i,j,k) - pres(i,j,k-1)
                   ! centered difference
                else
                   dt = state(i,j,k+1,temp_comp) - state(i,j,k-1,temp_comp)
                   dp = pres(i,j,k+1) - pres(i,j,k-1)
                endif

                nabla = pres(i,j,k) * dt / (dp * state(i,j,k,temp_comp))
             endif

             ad_excess(i,j,k) = nabla - nabla_ad(i,j,k)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine make_ad_excess_3d

  subroutine make_ad_excess_3d_sphr(ad_excess, ng_ad, state, ng_s, &
                                    normal, ng_n, lo, hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use network, only: nspec
    use probin_module, only: base_cutoff_density
    use bl_constants_module

    integer,         intent(in   ) :: lo(:), hi(:), ng_ad, ng_s, ng_n
    real(kind=dp_t), intent(  out) :: ad_excess(lo(1)-ng_ad:,lo(2)-ng_ad:,lo(3)-ng_ad:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind = dp_t), intent(in   ) :: normal(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)  

    real(kind=dp_t) :: pres(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: nabla_ad(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    real(kind=dp_t) :: chi_rho, chi_t, nabla
    real(kind=dp_t) :: dp(4), dt(4)

    integer :: i, j, k, c

    !$OMP PARALLEL DO PRIVATE(i,j,k,chi_rho,chi_t)    
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_eos(1) = state(i,j,k,rho_comp)
             temp_eos(1) = state(i,j,k,temp_comp)
             xn_eos(1,:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

             pt_index_eos(:) = (/i, j, k/)       

             call eos(eos_input_rt, den_eos, temp_eos, &
                      npts, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      .false., &
                      pt_index_eos)       
       
             pres(i,j,k) = p_eos(1)

             chi_rho = den_eos(1) * dpdr_eos(1) / p_eos(1)
             chi_t = temp_eos(1) * dpdt_eos(1) / p_eos(1)
             nabla_ad(i,j,k) = (gam1_eos(1) - chi_rho) / (chi_t * gam1_eos(1))

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j,k,dt,dp,nabla)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (state(i,j,k,rho_comp) <= base_cutoff_density) then
                nabla = ZERO
             else
                ! compute gradient

                ! forward difference
                if (k == lo(3)) then
                   dt(3) = state(i,j,k+1,temp_comp) - state(i,j,k,temp_comp)
                   dp(3) = pres(i,j,k+1) - pres(i,j,k)
                   ! backward difference
                else if (k == hi(3)) then
                   dt(3) = state(i,j,k,temp_comp) - state(i,j,k-1,temp_comp)
                   dp(3) = pres(i,j,k) - pres(i,j,k-1)
                   ! centered difference
                else
                   dt(3) = state(i,j,k+1,temp_comp) - state(i,j,k-1,temp_comp)
                   dp(3) = pres(i,j,k+1) - pres(i,j,k-1)
                endif

                if (j == lo(2)) then
                   dt(2) = state(i,j+1,k,temp_comp) - state(i,j,k,temp_comp)
                   dp(2) = pres(i,j+1,k) - pres(i,j,k)
                   ! backward difference
                else if (j == hi(2)) then
                   dt(2) = state(i,j,k,temp_comp) - state(i,j-1,k,temp_comp)
                   dp(2) = pres(i,j,k) - pres(i,j-1,k)
                   ! centered difference
                else
                   dt(2) = state(i,j+1,k,temp_comp) - state(i,j-1,k,temp_comp)
                   dp(2) = pres(i,j+1,k) - pres(i,j-1,k)
                endif

                if (i == lo(1)) then
                   dt(1) = state(i+1,j,k,temp_comp) - state(i,j,k,temp_comp)
                   dp(1) = pres(i+1,j,k) - pres(i,j,k)
                   ! backward difference
                else if (i == hi(1)) then
                   dt(1) = state(i,j,k,temp_comp) - state(i-1,j,k,temp_comp)
                   dp(1) = pres(i,j,k) - pres(i-1,j,k)
                   ! centered difference
                else
                   dt(1) = state(i+1,j,k,temp_comp) - state(i-1,j,k,temp_comp)
                   dp(1) = pres(i+1,j,k) - pres(i-1,j,k)
                endif

                ! dot into normal to get d/dr
                dp(4) = 0.d0
                dt(4) = 0.d0
                do c = 1,3
                   dp(4) = dp(4) + dp(c)*normal(i,j,k,c) 
                   dt(4) = dt(4) + dt(c)*normal(i,j,k,c) 
                enddo

                nabla = pres(i,j,k)*dt(4) / (dp(4)*state(i,j,k,temp_comp))
             endif

             ad_excess(i,j,k) = (nabla - nabla_ad(i,j,k)) / nabla
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine make_ad_excess_3d_sphr


  !---------------------------------------------------------------------------
  ! make_conductivity
  !---------------------------------------------------------------------------
  subroutine make_conductivity(plotdata,comp_cond,state)

    use geometry, only: spherical

    type(multifab),  intent(inout) :: plotdata
    integer,         intent(in   ) :: comp_cond
    type(multifab),  intent(in   ) :: state

    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: cp(:,:,:,:)
    integer :: lo(get_dim(plotdata)), hi(get_dim(plotdata)), ng_s, ng_c
    integer :: i,dm

    dm = get_dim(plotdata)

    ng_s = nghost(state)
    ng_c = nghost(plotdata)

    do i = 1, nboxes(state)
       if (multifab_remote(state, i)) cycle
       sp => dataptr(state,i)
       cp => dataptr(plotdata,i)
       lo = lwb(get_box(state,i))
       hi = upb(get_box(state,i))
       select case (dm)
       case (1)
          call make_conductivity_1d(cp(:,1,1,comp_cond), ng_c, &
                                    sp(:,1,1,:), ng_s, lo, hi)
       case (2)
          call make_conductivity_2d(cp(:,:,1,comp_cond), ng_c, &
                                    sp(:,:,1,:), ng_s, lo, hi)
       case (3)
          call make_conductivity_3d(cp(:,:,:,comp_cond), ng_c, &
                                    sp(:,:,:,:), ng_s, lo, hi)
       end select

    enddo

  end subroutine make_conductivity

  subroutine make_conductivity_1d(cond,ng_c,state,ng_s,lo,hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use conductivity_module
    use network, only: nspec

    integer,         intent(in   ) :: lo(:), hi(:), ng_c, ng_s
    real(kind=dp_t), intent(  out) :: cond(lo(1)-ng_c:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,:)

    ! local
    integer :: i

    do i = lo(1), hi(1)

       den_eos(1)  = state(i,rho_comp)
       temp_eos(1) = state(i,temp_comp)
       xn_eos(1,:) = state(i,spec_comp:spec_comp+nspec-1) / den_eos(1)

       call conducteos(eos_input_rt, den_eos, temp_eos, &
                       npts, nspec, &
                       xn_eos, &
                       p_eos, h_eos, e_eos, &
                       cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                       dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                       dpdX_eos, dhdX_eos, &
                       gam1_eos, cs_eos, s_eos, &
                       dsdt_eos, dsdr_eos, &
                       .false., conduct_eos)

       cond(i) = conduct_eos(1)

    enddo

  end subroutine make_conductivity_1d

  subroutine make_conductivity_2d(cond,ng_c,state,ng_s,lo,hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use conductivity_module
    use network, only: nspec

    integer,         intent(in   ) :: lo(:), hi(:), ng_c, ng_s
    real(kind=dp_t), intent(  out) :: cond(lo(1)-ng_c:,lo(2)-ng_c:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,:)

    ! local
    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          den_eos(1) = state(i,j,rho_comp)
          temp_eos(1) = state(i,j,temp_comp)
          xn_eos(1,:) = state(i,j,spec_comp:spec_comp+nspec-1) / den_eos(1)

          call conducteos(eos_input_rt, den_eos, temp_eos, &
                          npts, nspec, &
                          xn_eos, &
                          p_eos, h_eos, e_eos, &
                          cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                          dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                          dpdX_eos, dhdX_eos, &
                          gam1_eos, cs_eos, s_eos, &
                          dsdt_eos, dsdr_eos, &
                          .false., conduct_eos)

          cond(i,j) = conduct_eos(1)

       enddo
    enddo

  end subroutine make_conductivity_2d

  subroutine make_conductivity_3d(cond,ng_c,state,ng_s,lo,hi)

    use variables, only: rho_comp, temp_comp, spec_comp
    use eos_module
    use conductivity_module
    use network, only: nspec

    integer,         intent(in   ) :: lo(:), hi(:), ng_c, ng_s
    real(kind=dp_t), intent(  out) :: cond(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
    real(kind=dp_t), intent(in   ) :: state(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)

    ! local
    integer :: i, j, k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_eos(1)  = state(i,j,k,rho_comp)
             temp_eos(1) = state(i,j,k,temp_comp)
             xn_eos(1,:) = state(i,j,k,spec_comp:spec_comp+nspec-1) / &
                           den_eos(1)

             call conducteos(eos_input_rt, den_eos, temp_eos, &
                             npts, nspec, &
                             xn_eos, &
                             p_eos, h_eos, e_eos, &
                             cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                             dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                             dpdX_eos, dhdX_eos, &
                             gam1_eos, cs_eos, s_eos, &
                             dsdt_eos, dsdr_eos, &
                             .false., conduct_eos)

             cond(i,j,k) = conduct_eos(1)

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine make_conductivity_3d


  !---------------------------------------------------------------------------
  ! make_pi_cc
  !---------------------------------------------------------------------------
  subroutine make_pi_cc(mla,pi,pi_cc,the_bc_level)

    use ml_layout_module
    use bc_module

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: pi(:)
    type(multifab) , intent(inout) :: pi_cc(:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    real(kind=dp_t), pointer :: ppn(:,:,:,:)
    real(kind=dp_t), pointer :: ppc(:,:,:,:)
    logical,         pointer ::  mp(:,:,:,:)

    integer :: i,n,ng_pn,ng_pc,dm,nlevs
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t) :: ncell_proc(mla%nlevel), ncell(mla%nlevel)
    real(kind=dp_t) :: pisum_proc(mla%nlevel), pisum(mla%nlevel)

    real(kind=dp_t) :: weight,avg

    dm = mla%dim
    nlevs = mla%nlevel

    ncell      = 0.d0
    pisum      = 0.d0
    ncell_proc = 0.d0
    pisum_proc = 0.d0
    ng_pn      = nghost(pi(1))
    ng_pc      = nghost(pi_cc(1))

    do n=1,nlevs
       weight = 2.d0**(dm*(n-1))
       do i=1,nboxes(pi_cc(n))
          if ( multifab_remote(pi_cc(n), i) ) cycle
          ppn => dataptr(pi(n), i)
          ppc => dataptr(pi_cc(n), i)
          lo  =  lwb(get_box(pi_cc(n), i))
          hi  =  upb(get_box(pi_cc(n), i))
          select case (dm)
          case (1)
             if (n .eq. nlevs) then
                call make_pi_cc_1d(weight,ppn(:,1,1,1),ng_pn,ppc(:,1,1,1),ng_pc, &
                                   lo,hi,ncell_proc(n),pisum_proc(n))
             else
                mp => dataptr(mla%mask(n), i)
                call make_pi_cc_1d(weight,ppn(:,1,1,1),ng_pn,ppc(:,1,1,1),ng_pc, &
                                   lo,hi,ncell_proc(n),pisum_proc(n),mp(:,1,1,1))
             end if
          case (2)
             if (n .eq. nlevs) then
                call make_pi_cc_2d(weight,ppn(:,:,1,1),ng_pn,ppc(:,:,1,1),ng_pc, &
                                   lo,hi,ncell_proc(n),pisum_proc(n))
             else
                mp => dataptr(mla%mask(n), i)
                call make_pi_cc_2d(weight,ppn(:,:,1,1),ng_pn,ppc(:,:,1,1),ng_pc, &
                                   lo,hi,ncell_proc(n),pisum_proc(n),mp(:,:,1,1))
             end if
          case (3)
             if (n .eq. nlevs) then
                call make_pi_cc_3d(weight,ppn(:,:,:,1),ng_pn,ppc(:,:,:,1),ng_pc, &
                                   lo,hi,ncell_proc(n),pisum_proc(n))
             else
                mp => dataptr(mla%mask(n), i)
                call make_pi_cc_3d(weight,ppn(:,:,:,1),ng_pn,ppc(:,:,:,1),ng_pc, &
                                   lo,hi,ncell_proc(n),pisum_proc(n),mp(:,:,:,1))
             end if
          end select
       end do
    end do

    call parallel_reduce(ncell, ncell_proc, MPI_SUM)
    call parallel_reduce(pisum, pisum_proc, MPI_SUM)

    ! now ncell will contain the total number of cells over all levels
    ! now picum will contain the sum of (volume weighted) pi over all levels
    do n=2,nlevs
       ncell(1) = ncell(1) + ncell(n)
       pisum(1) = pisum(1) + pisum(n)
    end do

    ! divide the sum by the number of cells
    avg = pisum(1)/ncell(1)

    ! if there are no outlet boundary conditions, normalize pi_cc so the
    ! sum over the domain is zero
    if (.not.(any(the_bc_level(1)%phys_bc_level_array(:,:,:) .eq. OUTLET))) then
       do n=1,nlevs
          call multifab_sub_sub_s(pi_cc(n),avg,ng_pc)
       end do
    end if

  end subroutine make_pi_cc

  subroutine make_pi_cc_1d(weight,pi,ng_pn,pi_cc,ng_pc,lo,hi,ncell,pisum,mask)

    real (kind=dp_t), intent(in   )           :: weight
    integer         , intent(in   )           :: lo(:), hi(:), ng_pn, ng_pc
    real (kind=dp_t), intent(in   )           ::    pi(lo(1)-ng_pn:)
    real (kind=dp_t), intent(inout)           :: pi_cc(lo(1)-ng_pc:)
    real (kind=dp_t), intent(inout)           :: ncell,pisum
    logical         , intent(in   ), optional ::  mask(lo(1):      )

    ! local
    integer :: i

    logical :: cell_valid

    do i=lo(1),hi(1)

       pi_cc(i) = (pi(i) + pi(i+1)) / 2.d0

       ! make sure the cell isn't covered by finer cells
       cell_valid = .true.
       if ( present(mask) ) then
          cell_valid = mask(i)
       end if
       
       if (cell_valid) then
          pisum = pisum + weight*pi_cc(i)
          ncell = ncell + weight
       end if
       
    end do

  end subroutine make_pi_cc_1d

  subroutine make_pi_cc_2d(weight,pi,ng_pn,pi_cc,ng_pc,lo,hi,ncell,pisum,mask)

    real (kind=dp_t), intent(in   )           :: weight
    integer         , intent(in   )           :: lo(:), hi(:), ng_pn, ng_pc
    real (kind=dp_t), intent(in   )           ::    pi(lo(1)-ng_pn:,lo(2)-ng_pn:)
    real (kind=dp_t), intent(inout)           :: pi_cc(lo(1)-ng_pc:,lo(2)-ng_pc:)
    real (kind=dp_t), intent(inout)           :: ncell,pisum
    logical         , intent(in   ), optional ::  mask(lo(1):      ,lo(2):      )

    ! local
    integer :: i,j

    logical :: cell_valid

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          pi_cc(i,j) = (pi(i,j) + pi(i+1,j) + pi(i,j+1) + pi(i+1,j+1)) / 4.d0

          ! make sure the cell isn't covered by finer cells
          cell_valid = .true.
          if ( present(mask) ) then
             cell_valid = mask(i,j)
          end if

          if (cell_valid) then
             pisum = pisum + weight*pi_cc(i,j)
             ncell = ncell + weight
          end if
             
       end do
    end do

  end subroutine make_pi_cc_2d

  subroutine make_pi_cc_3d(weight,pi,ng_pn,pi_cc,ng_pc,lo,hi,ncell,pisum,mask)

    real(kind=dp_t), intent(in   )           :: weight
    integer        , intent(in   )           :: lo(:), hi(:), ng_pn, ng_pc
    real(kind=dp_t), intent(in   )           ::    pi(lo(1)-ng_pn:,lo(2)-ng_pn:,lo(3)-ng_pn:)
    real(kind=dp_t), intent(inout)           :: pi_cc(lo(1)-ng_pc:,lo(2)-ng_pc:,lo(3)-ng_pc:)
    real(kind=dp_t), intent(inout)           :: ncell,pisum
    logical        , intent(in   ), optional ::  mask(lo(1):      ,lo(2):      ,lo(3):      )

    ! local
    integer :: i,j,k

    logical :: cell_valid

    !$OMP PARALLEL DO PRIVATE(i,j,k,cell_valid) reduction(+:pisum,ncell)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             
             pi_cc(i,j,k) = (pi(i,j,k) + pi(i+1,j,k) + pi(i,j+1,k) + pi(i,j,k+1) &
                  + pi(i+1,j+1,k) + pi(i+1,j,k+1) + pi(i,j+1,k+1) + pi(i+1,j+1,k+1)) / 8.d0
             
             ! make sure the cell isn't covered by finer cells
             cell_valid = .true.
             if ( present(mask) ) then
                cell_valid = mask(i,j,k)
             end if
             
             if (cell_valid) then
                pisum = pisum + weight*pi_cc(i,j,k)
                ncell = ncell + weight
             end if
             
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine make_pi_cc_3d


  !---------------------------------------------------------------------------
  ! make_tfromH
  !---------------------------------------------------------------------------
  subroutine make_tfromH(plotdata,comp_t,comp_tpert,comp_dp,state,p0,tempbar,dx)

    use geometry, only: spherical

    integer        , intent(in   ) :: comp_t,comp_tpert,comp_dp
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: state
    real(kind=dp_t), intent(in   ) :: p0(0:),tempbar(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: tp(:,:,:,:)
    integer :: lo(get_dim(plotdata)),hi(get_dim(plotdata)),ng_s,ng_p
    integer :: i,dm

    dm = get_dim(plotdata)

    ng_s = nghost(state)
    ng_p = nghost(plotdata)

    do i = 1, nboxes(state)
       if ( multifab_remote(state, i) ) cycle
       sp => dataptr(state, i)
       tp => dataptr(plotdata, i)
       lo =  lwb(get_box(state, i))
       hi =  upb(get_box(state, i))
       select case (dm)
       case (1)
          call make_tfromH_1d(tp(:,1,1,comp_t),tp(:,1,1,comp_tpert), &
                              tp(:,1,1,comp_dp),ng_p,sp(:,1,1,:),ng_s, &
                              lo,hi,p0,tempbar)
       case (2)
          call make_tfromH_2d(tp(:,:,1,comp_t),tp(:,:,1,comp_tpert), &
                              tp(:,:,1,comp_dp),ng_p,sp(:,:,1,:),ng_s, &
                              lo,hi,p0,tempbar)
       case (3)
          if (spherical .eq. 1) then
             call make_tfromH_3d_sphr(tp(:,:,:,comp_t),tp(:,:,:,comp_tpert), &
                                      tp(:,:,:,comp_dp),ng_p, &
                                      sp(:,:,:,:),ng_s,lo,hi,p0,tempbar,dx)
          else
             call make_tfromH_3d_cart(tp(:,:,:,comp_t),tp(:,:,:,comp_tpert), &
                                      tp(:,:,:,comp_dp),ng_p, &
                                      sp(:,:,:,:),ng_s,lo,hi,p0,tempbar)
          end if
       end select
    end do

  end subroutine make_tfromH

  subroutine make_tfromH_1d(T,tpert,deltaP,ng_p,state,ng_s,lo,hi,p0,tempbar)

    use variables, only: rho_comp, spec_comp, rhoh_comp, temp_comp
    use eos_module
    use network, only: nspec
    use bl_constants_module

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind = dp_t), intent(  out) ::      T(lo(1)-ng_p:)
    real (kind = dp_t), intent(  out) ::  tpert(lo(1)-ng_p:)
    real (kind = dp_t), intent(  out) :: deltaP(lo(1)-ng_p:)
    real (kind = dp_t), intent(in   ) ::  state(lo(1)-ng_s:,:)
    real (kind = dp_t), intent(in   ) ::  p0(0:),tempbar(0:)

    ! Local variables
    integer :: i

    do i = lo(1), hi(1)

       ! (rho, H) --> T, p
       den_eos(1)  = state(i,rho_comp)
       p_eos(1)    = p0(i)
       temp_eos(1) = state(i,temp_comp)
       xn_eos(1,:) = state(i,spec_comp:spec_comp+nspec-1)/den_eos(1)
       h_eos(1) = state(i,rhoh_comp) / state(i,rho_comp)

       pt_index_eos(:) = (/i, -1, -1/)       

       call eos(eos_input_rh, den_eos, temp_eos, &
                npts, &
                xn_eos, &
                p_eos, h_eos, e_eos, &
                cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                dpdX_eos, dhdX_eos, &
                gam1_eos, cs_eos, s_eos, &
                dsdt_eos, dsdr_eos, &
                .false., &
                pt_index_eos)

       T(i) = temp_eos(1)
       if (.not. use_tfromp) tpert(i) = temp_eos(1) - tempbar(i)

       deltaP(i) = abs(p_eos(1)-p0(i))/ p0(i)

    enddo

  end subroutine make_tfromH_1d

  subroutine make_tfromH_2d(T,tpert,deltaP,ng_p,state,ng_s,lo,hi,p0,tempbar)

    use variables, only: rho_comp, spec_comp, rhoh_comp, temp_comp
    use eos_module
    use network, only: nspec
    use bl_constants_module

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind = dp_t), intent(  out) ::      T(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind = dp_t), intent(  out) ::  tpert(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind = dp_t), intent(  out) :: deltaP(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind = dp_t), intent(in   ) ::  state(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind = dp_t), intent(in   ) ::  p0(0:),tempbar(0:)

    ! Local variables
    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          ! (rho, H) --> T, p
          den_eos(1)  = state(i,j,rho_comp)
          p_eos(1)    = p0(j)
          temp_eos(1) = state(i,j,temp_comp)
          xn_eos(1,:) = state(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)
          h_eos(1) = state(i,j,rhoh_comp) / state(i,j,rho_comp)

          pt_index_eos(:) = (/i, j, -1/)

          call eos(eos_input_rh, den_eos, temp_eos, &
                   npts, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, &
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   .false., &
                   pt_index_eos)

          T(i,j) = temp_eos(1)
          if (.not. use_tfromp) tpert(i,j) = temp_eos(1) - tempbar(j)

          deltaP(i,j) = abs(p_eos(1)-p0(j))/ p0(j)

       enddo
    enddo

  end subroutine make_tfromH_2d

  subroutine make_tfromH_3d_cart(T,tpert,deltaP,ng_p,state,ng_s,lo,hi,p0,tempbar)

    use variables, only: rho_comp, spec_comp, rhoh_comp, temp_comp
    use eos_module
    use network, only: nspec

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind = dp_t), intent(  out) ::      T(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind = dp_t), intent(  out) ::  tpert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind = dp_t), intent(  out) :: deltaP(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind = dp_t), intent(in   ) ::  state(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind = dp_t), intent(in   ) :: p0(0:),tempbar(0:)

    ! Local variables
    integer :: i, j, k

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! (rho, H) --> T, p
             den_eos(1)  = state(i,j,k,rho_comp)
             p_eos(1)    = p0(k)
             temp_eos(1) = state(i,j,k,temp_comp)
             xn_eos(1,:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)
             h_eos(1) = state(i,j,k,rhoh_comp)/state(i,j,k,rho_comp)

             pt_index_eos(:) = (/i, j, k/)

             call eos(eos_input_rh, den_eos, temp_eos, &
                      npts, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      .false., &
                      pt_index_eos)

             T(i,j,k) = temp_eos(1)
             if (.not. use_tfromp) tpert(i,j,k) = temp_eos(1) - tempbar(k)

             deltaP(i,j,k) = (p_eos(1)-p0(k))/ p0(k)

          enddo
       enddo
    enddo

  end subroutine make_tfromH_3d_cart

  subroutine make_tfromH_3d_sphr(T,tpert,deltaP,ng_p,state,ng_s,lo,hi,p0,tempbar,dx)

    use variables, only: rho_comp, rhoh_comp, spec_comp, temp_comp
    use eos_module
    use network, only: nspec
    use fill_3d_module

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind = dp_t), intent(  out) ::      T(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind = dp_t), intent(  out) ::  tpert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind = dp_t), intent(  out) :: deltaP(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind = dp_t), intent(in   ) ::  state(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind = dp_t), intent(in   ) :: p0(0:),tempbar(0:)
    real (kind = dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer          :: i, j, k
    real (kind=dp_t), allocatable :: p0_cart(:,:,:,:)
    real (kind=dp_t), allocatable :: tempbar_cart(:,:,:,:)

    allocate(p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    allocate(  tempbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,p0,p0_cart,lo,hi,dx,0)
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,tempbar,tempbar_cart,lo,hi,dx,0)

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_eos(1)  = state(i,j,k,rho_comp)
             h_eos(1)    = state(i,j,k,rhoh_comp) / state(i,j,k,rho_comp)
             temp_eos(1) = state(i,j,k,temp_comp)
             xn_eos(1,:) = state(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

             pt_index_eos(:) = (/i, j, k/)

             ! (rho, H) --> T, p
             call eos(eos_input_rh, den_eos, temp_eos, &
                      npts, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, &
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      .false., &
                      pt_index_eos)

             T(i,j,k) = temp_eos(1)
             if (.not. use_tfromp) tpert(i,j,k) = temp_eos(1) - tempbar_cart(i,j,k,1)
             
             deltaP(i,j,k) = (p_eos(1)-p0_cart(i,j,k,1))/ p0_cart(i,j,k,1)

          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(p0_cart,tempbar_cart)

  end subroutine make_tfromH_3d_sphr


  !---------------------------------------------------------------------------
  ! make_tfromp
  !---------------------------------------------------------------------------
  subroutine make_tfromp(plotdata,comp_tfromp,comp_tpert, &
                         comp_machno,comp_cs,comp_deltag,comp_entropy, comp_magvel, &
                         s,tempbar,gamma1bar,p0,dx)

    use geometry, only: spherical

    integer        , intent(in   ) :: comp_tfromp,comp_tpert
    integer        , intent(in   ) :: comp_machno,comp_cs
    integer        , intent(in   ) :: comp_deltag, comp_entropy, comp_magvel
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: tempbar(0:)
    real(kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real(kind=dp_t), intent(in   ) :: p0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), pointer:: sp(:,:,:,:),tp(:,:,:,:)
    integer :: lo(get_dim(plotdata)),hi(get_dim(plotdata)),i
    integer :: ng_p,ng_s,dm

    dm = get_dim(plotdata)

    ng_p = nghost(plotdata)
    ng_s = nghost(s)

    do i = 1, nboxes(s)
       if ( multifab_remote(s, i) ) cycle
       tp => dataptr(plotdata, i)
       sp => dataptr(s, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (1)
          call make_tfromp_1d(tp(:,1,1,comp_tfromp  ),tp(:,1,1,comp_tpert), &
                              tp(:,1,1,comp_machno  ), tp(:,1,1,comp_cs  ), &
                              tp(:,1,1,comp_deltag), &
                              tp(:,1,1,comp_entropy ),tp(:,1,1,comp_magvel), &
                              ng_p, &
                              sp(:,1,1,:), ng_s, &
                              lo, hi, tempbar, gamma1bar, p0)
       case (2)
          call make_tfromp_2d(tp(:,:,1,comp_tfromp),tp(:,:,1,comp_tpert), &
                              tp(:,:,1,comp_machno  ), tp(:,:,1,comp_cs), &
                              tp(:,:,1,comp_deltag), &
                              tp(:,:,1,comp_entropy ),tp(:,:,1,comp_magvel), &
                              ng_p, &
                              sp(:,:,1,:), ng_s, &
                              lo, hi, tempbar, gamma1bar, p0)
       case (3)
          if (spherical .eq. 1) then
             call make_tfromp_3d_sphr(tp(:,:,:,comp_tfromp),tp(:,:,:,comp_tpert), &
                                      tp(:,:,:,comp_machno  ), tp(:,:,:,comp_cs), &
                                      tp(:,:,:,comp_deltag), &
                                      tp(:,:,:,comp_entropy ),tp(:,:,:,comp_magvel), &
                                      ng_p, &
                                      sp(:,:,:,:), ng_s, &
                                      lo, hi, tempbar, gamma1bar, p0, dx)
          else
             call make_tfromp_3d_cart(tp(:,:,:,comp_tfromp),tp(:,:,:,comp_tpert), &
                                      tp(:,:,:,comp_machno  ), tp(:,:,:,comp_cs), &
                                      tp(:,:,:,comp_deltag), &
                                      tp(:,:,:,comp_entropy ),tp(:,:,:,comp_magvel), &
                                      ng_p, &
                                      sp(:,:,:,:), ng_s, &
                                      lo, hi, tempbar, gamma1bar, p0)
          endif
       end select
    end do

  end subroutine make_tfromp

  subroutine make_tfromp_1d(t,tpert,machno,cs,deltagamma,entropy,magvel, &
                            ng_p,s,ng_s,lo,hi,tempbar,gamma1bar,p0)

    use eos_module
    use network, only: nspec
    use variables, only: rho_comp, spec_comp, temp_comp
    use probin_module, only: plot_cs

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind=dp_t), intent(  out) ::          t(lo(1)-ng_p:)  
    real (kind=dp_t), intent(  out) ::      tpert(lo(1)-ng_p:)  
    real (kind=dp_t), intent(  out) ::     machno(lo(1)-ng_p:)  
    real (kind=dp_t), intent(  out) ::         cs(lo(1)-ng_p:)  
    real (kind=dp_t), intent(  out) :: deltagamma(lo(1)-ng_p:)  
    real (kind=dp_t), intent(  out) ::    entropy(lo(1)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::     magvel(lo(1)-ng_p:)
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: tempbar(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)

    !     Local variables
    integer          :: i

    ! Then compute the perturbation
    do i = lo(1), hi(1)

       den_eos(1) = s(i,rho_comp)
       temp_eos(1) = s(i,temp_comp)
       p_eos(1) = p0(i)
       xn_eos(1,:) = s(i,spec_comp:spec_comp+nspec-1)/den_eos(1)

       pt_index_eos(:) = (/i, -1, -1/)

       ! (rho,P) --> T,h
       call eos(eos_input_rp, den_eos, temp_eos, &
                npts, &
                xn_eos, &
                p_eos, h_eos, e_eos, & 
                cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                dpdX_eos, dhdX_eos, &
                gam1_eos, cs_eos, s_eos, &
                dsdt_eos, dsdr_eos, &
                .false., &
                pt_index_eos)

       t(i) = temp_eos(1)
       if (use_tfromp) tpert(i) = temp_eos(1) - tempbar(i)

       if (plot_cs) cs(i) = cs_eos(1)

       machno(i) = magvel(i) / cs_eos(1)
       deltagamma(i) = gam1_eos(1) - gamma1bar(i)

       entropy(i) = s_eos(1)
    enddo

  end subroutine make_tfromp_1d

  subroutine make_tfromp_2d(t,tpert,machno,cs,deltagamma,entropy,magvel, &
                            ng_p,s,ng_s,lo,hi,tempbar,gamma1bar,p0)

    use eos_module
    use network, only: nspec
    use variables, only: rho_comp, spec_comp, temp_comp
    use probin_module, only: plot_cs

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind=dp_t), intent(  out) ::          t(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(  out) ::      tpert(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(  out) ::     machno(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(  out) ::         cs(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(  out) :: deltagamma(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(  out) ::    entropy(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::     magvel(lo(1)-ng_p:,lo(2)-ng_p:)
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: tempbar(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)

    !     Local variables
    integer          :: i, j

    ! Then compute the perturbation
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)

          den_eos(1) = s(i,j,rho_comp)
          temp_eos(1) = s(i,j,temp_comp)
          p_eos(1) = p0(j)
          xn_eos(1,:) = s(i,j,spec_comp:spec_comp+nspec-1)/den_eos(1)

          pt_index_eos(:) = (/i, j, -1/)

          ! (rho,P) --> T,h
          call eos(eos_input_rp, den_eos, temp_eos, &
                   npts, &
                   xn_eos, &
                   p_eos, h_eos, e_eos, & 
                   cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                   dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                   dpdX_eos, dhdX_eos, &
                   gam1_eos, cs_eos, s_eos, &
                   dsdt_eos, dsdr_eos, &
                   .false., &
                   pt_index_eos)

          t(i,j) = temp_eos(1)
          if (use_tfromp) tpert(i,j) = temp_eos(1) - tempbar(j)

          if (plot_cs) cs(i,j) = cs_eos(1)

          machno(i,j) = magvel(i,j) / cs_eos(1)

          deltagamma(i,j) = gam1_eos(1) - gamma1bar(j)

          entropy(i,j) = s_eos(1)
       enddo
    enddo

  end subroutine make_tfromp_2d

  subroutine make_tfromp_3d_cart(t,tpert,machno,cs,deltagamma,entropy,magvel, &
                                 ng_p,s,ng_s,lo,hi,tempbar,gamma1bar,p0)

    use variables, only: rho_comp, spec_comp, temp_comp
    use eos_module
    use network, only: nspec
    use probin_module, only: plot_cs

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind=dp_t), intent(  out) ::          t(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::      tpert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::     machno(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::         cs(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) :: deltagamma(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::    entropy(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::     magvel(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: tempbar(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)

    ! Local variables
    integer          :: i, j, k

    ! Then compute the perturbation and Mach number
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             den_eos(1) = s(i,j,k,rho_comp)
             temp_eos(1) = s(i,j,k,temp_comp)
             p_eos(1) = p0(k)
             xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

             pt_index_eos(:) = (/i, j, k/)

             ! (rho,P) --> T,h
             call eos(eos_input_rp, den_eos, temp_eos, &
                      npts, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, & 
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      .false., &
                      pt_index_eos)

             t(i,j,k) = temp_eos(1)
             if (use_tfromp) tpert(i,j,k) = temp_eos(1) - tempbar(k)

             if (plot_cs) cs(i,j,k) = cs_eos(1)

             machno(i,j,k) = magvel(i,j,k) / cs_eos(1)

             deltagamma(i,j,k) = gam1_eos(1) - gamma1bar(k)

             entropy(i,j,k) = s_eos(1)
          enddo
       enddo
    enddo

  end subroutine make_tfromp_3d_cart

  subroutine make_tfromp_3d_sphr(t,tpert,machno,cs,deltagamma,entropy,magvel, &
                                 ng_p,s,ng_s,lo,hi,tempbar,gamma1bar,p0,dx)

    use variables, only: rho_comp, spec_comp, temp_comp
    use eos_module
    use network, only: nspec
    use fill_3d_module
    use probin_module, only: plot_cs

    integer         , intent(in   ) :: lo(:),hi(:),ng_p,ng_s
    real (kind=dp_t), intent(  out) ::          t(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::      tpert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::     machno(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::         cs(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) :: deltagamma(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) ::    entropy(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::     magvel(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: tempbar(0:)
    real (kind=dp_t), intent(in   ) :: gamma1bar(0:)
    real (kind=dp_t), intent(in   ) :: p0(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer          :: i, j, k


    real (kind=dp_t), allocatable ::   tempbar_cart(:,:,:,:)
    real (kind=dp_t), allocatable ::        p0_cart(:,:,:,:)
    real (kind=dp_t), allocatable :: gamma1bar_cart(:,:,:,:)

    allocate(  tempbar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    allocate(       p0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))
    allocate(gamma1bar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,tempbar,tempbar_cart,lo,hi,dx,0)
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,p0,p0_cart,lo,hi,dx,0)
    call put_1d_array_on_cart_3d_sphr(.false.,.false.,gamma1bar,gamma1bar_cart,lo,hi,dx,0)

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! Then compute the perturbation and Mach number
             den_eos(1) = s(i,j,k,rho_comp)
             temp_eos(1) = s(i,j,k,temp_comp)
             p_eos(1) = p0_cart(i,j,k,1)
             xn_eos(1,:) = s(i,j,k,spec_comp:spec_comp+nspec-1)/den_eos(1)

             pt_index_eos(:) = (/i, j, k/)

             ! (rho,P) --> T,h
             call eos(eos_input_rp, den_eos, temp_eos, &
                      npts, &
                      xn_eos, &
                      p_eos, h_eos, e_eos, & 
                      cv_eos, cp_eos, xne_eos, eta_eos, pele_eos, &
                      dpdt_eos, dpdr_eos, dedt_eos, dedr_eos, &
                      dpdX_eos, dhdX_eos, &
                      gam1_eos, cs_eos, s_eos, &
                      dsdt_eos, dsdr_eos, &
                      .false., &
                      pt_index_eos)

             t(i,j,k) = temp_eos(1)
             if (use_tfromp) tpert(i,j,k) = temp_eos(1) - tempbar_cart(i,j,k,1)

             if (plot_cs) cs(i,j,k) = cs_eos(1)

             machno(i,j,k) = magvel(i,j,k) / cs_eos(1)

             deltagamma(i,j,k) = gam1_eos(1) - gamma1bar_cart(i,j,k,1)

             entropy(i,j,k) = s_eos(1)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(tempbar_cart,p0_cart,gamma1bar_cart)

  end subroutine make_tfromp_3d_sphr


  !---------------------------------------------------------------------------
  ! make_rhopert
  !---------------------------------------------------------------------------
  subroutine make_rhopert(plotdata,comp_rhopert,s,rho0,dx)

    use geometry, only: spherical

    integer        , intent(in   ) :: comp_rhopert
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: rho0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), pointer:: sp(:,:,:,:),tp(:,:,:,:)
    integer :: lo(get_dim(plotdata)),hi(get_dim(plotdata)),i
    integer :: ng_p,ng_s,dm

    dm = get_dim(plotdata)

    ng_p = nghost(plotdata)
    ng_s = nghost(s)

    do i = 1, nboxes(s)
       if ( multifab_remote(s, i) ) cycle
       tp => dataptr(plotdata, i)
       sp => dataptr(s, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (1)
          call make_rhopert_1d(tp(:,1,1,comp_rhopert ),ng_p, &
                               sp(:,1,1,:), ng_s, &
                               lo, hi, rho0)
       case (2)
          call make_rhopert_2d(tp(:,:,1,comp_rhopert ),ng_p, &
                               sp(:,:,1,:), ng_s, &
                               lo, hi, rho0)
       case (3)
          if (spherical .eq. 1) then
             call make_rhopert_3d_sphr(tp(:,:,:,comp_rhopert ),ng_p, &
                                       sp(:,:,:,:), ng_s, &
                                       lo, hi, rho0, dx)
          else
             call make_rhopert_3d_cart(tp(:,:,:,comp_rhopert ),ng_p, &
                                       sp(:,:,:,:), ng_s, &
                                       lo, hi, rho0)
          endif
       end select
    end do

  end subroutine make_rhopert

  subroutine make_rhopert_1d(rhopert,ng_p,s,ng_s,lo,hi,rho0)

    use variables, only: rho_comp

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:)

    !     Local variables
    integer          :: i

    ! Then compute the perturbation
    do i = lo(1), hi(1)
       rhopert(i)  = s(i,rho_comp)  - rho0(i)
    enddo

  end subroutine make_rhopert_1d

  subroutine make_rhopert_2d(rhopert,ng_p,s,ng_s,lo,hi,rho0)

    use variables, only: rho_comp

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:)

    !     Local variables
    integer          :: i, j

    ! Then compute the perturbation
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          rhopert(i,j)  = s(i,j,rho_comp)  - rho0(j)
       enddo
    enddo

  end subroutine make_rhopert_2d

  subroutine make_rhopert_3d_cart(rhopert,ng_p,s,ng_s,lo,hi,rho0)

    use variables, only: rho_comp

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:)

    ! Local variables
    integer          :: i, j, k

    ! Then compute the perturbation and Mach number
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhopert(i,j,k)  = s(i,j,k,rho_comp)  - rho0(k)
          enddo
       enddo
    enddo

  end subroutine make_rhopert_3d_cart

  subroutine make_rhopert_3d_sphr(rhopert,ng_p,s,ng_s,lo,hi,rho0,dx)

    use variables, only: rho_comp
    use fill_3d_module

    integer         , intent(in   ) :: lo(:),hi(:),ng_p,ng_s
    real (kind=dp_t), intent(  out) ::    rhopert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rho0(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer          :: i, j, k

    real (kind=dp_t), allocatable ::      rho0_cart(:,:,:,:)

    allocate(     rho0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,rho0,rho0_cart,lo,hi,dx,0)


    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhopert(i,j,k)  = s(i,j,k,rho_comp)  -  rho0_cart(i,j,k,1)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(rho0_cart)

  end subroutine make_rhopert_3d_sphr


  !---------------------------------------------------------------------------
  ! make_rhohpert
  !---------------------------------------------------------------------------
  subroutine make_rhohpert(plotdata,comp_rhohpert,s,rhoh0,dx)

    use geometry, only: spherical

    integer        , intent(in   ) :: comp_rhohpert
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: rhoh0(0:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    real(kind=dp_t), pointer:: sp(:,:,:,:),tp(:,:,:,:)
    integer :: lo(get_dim(plotdata)),hi(get_dim(plotdata)),i
    integer :: ng_p,ng_s,dm

    dm = get_dim(plotdata)

    ng_p = nghost(plotdata)
    ng_s = nghost(s)

    do i = 1, nboxes(s)
       if ( multifab_remote(s, i) ) cycle
       tp => dataptr(plotdata, i)
       sp => dataptr(s, i)
       lo =  lwb(get_box(s, i))
       hi =  upb(get_box(s, i))
       select case (dm)
       case (1)
          call make_rhohpert_1d(tp(:,1,1,comp_rhohpert ),ng_p, &
                                sp(:,1,1,:), ng_s, &
                                lo, hi, rhoh0)
       case (2)
          call make_rhohpert_2d(tp(:,:,1,comp_rhohpert ),ng_p, &
                                sp(:,:,1,:), ng_s, &
                                lo, hi, rhoh0)
       case (3)
          if (spherical .eq. 1) then
             call make_rhohpert_3d_sphr(tp(:,:,:,comp_rhohpert ),ng_p, &
                                        sp(:,:,:,:), ng_s, &
                                        lo, hi, rhoh0, dx)
          else
             call make_rhohpert_3d_cart(tp(:,:,:,comp_rhohpert ),ng_p, &
                                        sp(:,:,:,:), ng_s, &
                                        lo, hi, rhoh0)
          endif
       end select
    end do

  end subroutine make_rhohpert

  subroutine make_rhohpert_1d(rhohpert,ng_p,s,ng_s,lo,hi,rhoh0)

    use variables, only: rhoh_comp

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind=dp_t), intent(  out) ::   rhohpert(lo(1)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rhoh0(0:)

    !     Local variables
    integer          :: i

    ! Then compute the perturbation
    do i = lo(1), hi(1)
       rhohpert(i)  = s(i,rhoh_comp)  - rhoh0(i)
    enddo

  end subroutine make_rhohpert_1d

  subroutine make_rhohpert_2d(rhohpert,ng_p,s,ng_s,lo,hi,rhoh0)

    use variables, only: rhoh_comp

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind=dp_t), intent(  out) ::   rhohpert(lo(1)-ng_p:,lo(2)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rhoh0(0:)

    !     Local variables
    integer          :: i, j

    ! Then compute the perturbation
    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          rhohpert(i,j)  = s(i,j,rhoh_comp)  - rhoh0(j)
       enddo
    enddo

  end subroutine make_rhohpert_2d

  subroutine make_rhohpert_3d_cart(rhohpert,ng_p,s,ng_s,lo,hi,rhoh0)

    use variables, only: rhoh_comp

    integer, intent(in) :: lo(:), hi(:), ng_p, ng_s
    real (kind=dp_t), intent(  out) ::   rhohpert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rhoh0(0:)

    ! Local variables
    integer          :: i, j, k

    ! Then compute the perturbation and Mach number
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhohpert(i,j,k)  = s(i,j,k,rhoh_comp)  - rhoh0(k)
          enddo
       enddo
    enddo

  end subroutine make_rhohpert_3d_cart

  subroutine make_rhohpert_3d_sphr(rhohpert,ng_p,s,ng_s,lo,hi,rhoh0,dx)

    use variables, only: rhoh_comp
    use fill_3d_module

    integer         , intent(in   ) :: lo(:),hi(:),ng_p,ng_s
    real (kind=dp_t), intent(  out) ::    rhohpert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(in   ) ::          s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real (kind=dp_t), intent(in   ) :: rhoh0(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer          :: i, j, k

    real (kind=dp_t), allocatable ::      rhoh0_cart(:,:,:,:)

    allocate(     rhoh0_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,rhoh0,rhoh0_cart,lo,hi,dx,0)


    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rhohpert(i,j,k)  = s(i,j,k,rhoh_comp)  -  rhoh0_cart(i,j,k,1)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(rhoh0_cart)

  end subroutine make_rhohpert_3d_sphr


  !---------------------------------------------------------------------------
  ! make_entropypert
  !---------------------------------------------------------------------------
  subroutine make_entropypert(plotdata,comp_entropy,comp_entropypert,entropybar,dx)

    use geometry, only: spherical

    integer        , intent(in   ) :: comp_entropy,comp_entropypert
    type(multifab) , intent(inout) :: plotdata(:)
    real(kind=dp_t), intent(in   ) :: entropybar(:,0:)
    real(kind=dp_t), intent(in   ) ::         dx(:,:)

    ! local
    real(kind=dp_t), pointer :: tp(:,:,:,:)
    integer                  :: lo(get_dim(plotdata(1))),hi(get_dim(plotdata(1)))
    integer                  :: ng_p,n,i,dm,nlevs

    dm = get_dim(plotdata(1))
    nlevs = size(plotdata)

    ng_p = nghost(plotdata(1))

    do n=1,nlevs
       do i = 1, nboxes(plotdata(n))
          if ( multifab_remote(plotdata(n), i) ) cycle
          tp => dataptr(plotdata(n), i)
          lo =  lwb(get_box(plotdata(n), i))
          hi =  upb(get_box(plotdata(n), i))
          select case (dm)
          case (1)
             call make_entropypert_1d(tp(:,1,1,comp_entropy), &
                                      tp(:,1,1,comp_entropypert),ng_p, &
                                      lo, hi, entropybar(n,:))
          case (2)
             call make_entropypert_2d(tp(:,:,1,comp_entropy), &
                                      tp(:,:,1,comp_entropypert),ng_p, &
                                      lo, hi, entropybar(n,:))
          case (3)
             if (spherical .eq. 1) then
                call make_entropypert_3d_sphr(tp(:,:,:,comp_entropy), &
                                              tp(:,:,:,comp_entropypert),ng_p, &
                                              lo, hi, entropybar(1,:), dx(n,:))
             else
                call make_entropypert_3d_cart(tp(:,:,:,comp_entropy), &
                                              tp(:,:,:,comp_entropypert),ng_p, &
                                              lo, hi, entropybar(n,:))
             endif
          end select
       end do
    end do

  end subroutine make_entropypert

  subroutine make_entropypert_1d(entropy,entropypert,ng_p,lo,hi,entropybar)

    integer, intent(in) :: lo(:), hi(:), ng_p
    real (kind=dp_t), intent(inout) ::     entropy(lo(1)-ng_p:)
    real (kind=dp_t), intent(  out) :: entropypert(lo(1)-ng_p:)
    real (kind=dp_t), intent(in   ) :: entropybar(0:)

    !     Local variables
    integer          :: i

    ! Compute entropy - entropybar
    do i = lo(1), hi(1)
       entropypert(i) = (entropy(i) - entropybar(i))/entropybar(i)
    enddo

  end subroutine make_entropypert_1d

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

  subroutine make_entropypert_3d_sphr(entropy,entropypert,ng_p,lo,hi,entropybar,dx)

    use fill_3d_module

    integer         , intent(in)    :: lo(:),hi(:),ng_p
    real (kind=dp_t), intent(inout) ::     entropy(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(  out) :: entropypert(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)  
    real (kind=dp_t), intent(in   ) :: entropybar(0:)
    real (kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer          :: i, j, k
    real (kind=dp_t), allocatable :: entropybar_cart(:,:,:,:)

    allocate(entropybar_cart(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),1))

    call put_1d_array_on_cart_3d_sphr(.false.,.false.,entropybar,entropybar_cart, &
                                      lo,hi,dx,0)

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             ! Compute entropy-entropybar
             entropypert(i,j,k) = &
                  (entropy(i,j,k) - entropybar_cart(i,j,k,1))/entropybar_cart(i,j,k,1)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    deallocate(entropybar_cart)

  end subroutine make_entropypert_3d_sphr


  !---------------------------------------------------------------------------
  ! make_deltaT
  !---------------------------------------------------------------------------
  subroutine make_deltaT(plotdata,comp_dT,comp_tfromH,comp_tfromp)

    integer        , intent(in   ) :: comp_dT, comp_tfromH, comp_tfromp
    type(multifab) , intent(inout) :: plotdata

    real(kind=dp_t), pointer:: tp(:,:,:,:)
    integer :: lo(get_dim(plotdata)),hi(get_dim(plotdata))
    integer :: i,ng_p,dm

    dm = get_dim(plotdata)

    ng_p = nghost(plotdata)

    do i = 1, nboxes(plotdata)
       if ( multifab_remote(plotdata, i) ) cycle
       tp => dataptr(plotdata, i)
       lo =  lwb(get_box(plotdata, i))
       hi =  upb(get_box(plotdata, i))
       select case (dm)
       case (1)
          call make_deltaT_1d(tp(:,1,1,comp_dT),tp(:,1,1,comp_tfromH), &
                              tp(:,1,1,comp_tfromp), ng_p, lo, hi)
       case (2)
          call make_deltaT_2d(tp(:,:,1,comp_dT),tp(:,:,1,comp_tfromH), &
                              tp(:,:,1,comp_tfromp), ng_p, lo, hi)
       case (3)
          call make_deltaT_3d(tp(:,:,:,comp_dT),tp(:,:,:,comp_tfromH), &
                              tp(:,:,:,comp_tfromp), ng_p, lo, hi)
       end select
    end do

  end subroutine make_deltaT

  subroutine make_deltaT_1d(dT,tfromH,tfromp,ng_p,lo,hi)

    integer, intent(in) :: lo(:), hi(:), ng_p
    real (kind = dp_t), intent(  out) ::     dT(lo(1)-ng_p:)
    real (kind = dp_t), intent(in   ) :: tfromH(lo(1)-ng_p:)
    real (kind = dp_t), intent(in   ) :: tfromp(lo(1)-ng_p:)

    !     Local variables
    integer :: i

    do i = lo(1), hi(1)
       dT(i) = (tfromH(i) - tfromp(i)) / tfromH(i)
    end do

  end subroutine make_deltaT_1d

  subroutine make_deltaT_2d(dT,tfromH,tfromp,ng_p,lo,hi)

    integer, intent(in) :: lo(:), hi(:), ng_p
    real (kind = dp_t), intent(  out) ::     dT(lo(1)-ng_p:,lo(2)-ng_p:)
    real (kind = dp_t), intent(in   ) :: tfromH(lo(1)-ng_p:,lo(2)-ng_p:)
    real (kind = dp_t), intent(in   ) :: tfromp(lo(1)-ng_p:,lo(2)-ng_p:)

    !     Local variables
    integer :: i, j

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          dT(i,j) = (tfromH(i,j) - tfromp(i,j)) / tfromH(i,j)
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

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dT(i,j,k) = (tfromH(i,j,k) - tfromp(i,j,k)) / tfromH(i,j,k)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine make_deltaT_3d


  !---------------------------------------------------------------------------
  ! make_divw0
  !---------------------------------------------------------------------------
  subroutine make_divw0(divw0,comp_divw0,w0,w0mac,dx)

    use geometry, only: spherical

    type(multifab) , intent(inout) :: divw0
    integer        , intent(in   ) :: comp_divw0
    real(kind=dp_t), intent(in   ) :: w0(0:)
    type(multifab) , intent(in   ) :: w0mac(:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    real(kind=dp_t), pointer :: w0xp(:,:,:,:)
    real(kind=dp_t), pointer :: w0yp(:,:,:,:)
    real(kind=dp_t), pointer :: w0zp(:,:,:,:)
    real(kind=dp_t), pointer :: dwp(:,:,:,:)
    integer                  :: lo(get_dim(divw0)),hi(get_dim(divw0))
    integer                  :: i,ng_w0,ng_dw,dm

    dm = get_dim(divw0)

    ng_dw = nghost(divw0)

    do i=1,nboxes(divw0)
       if ( multifab_remote(divw0, i) ) cycle
       dwp => dataptr(divw0, i)
       lo = lwb(get_box(divw0, i))
       hi = upb(get_box(divw0, i))
       select case (dm)
       case (1)
          call make_divw0_1d(w0, dwp(:,1,1,comp_divw0), ng_dw, lo, hi, dx)
       case (2)
          call make_divw0_2d(w0, dwp(:,:,1,comp_divw0), ng_dw, lo, hi, dx)
       case (3)
          if(spherical .eq. 1) then
             ng_w0 = nghost(w0mac(1))
             w0xp => dataptr(w0mac(1), i)
             w0yp => dataptr(w0mac(2), i)
             w0zp => dataptr(w0mac(3), i)
             call make_divw0_3d_sphr(w0xp(:,:,:,1), w0yp(:,:,:,1), w0zp(:,:,:,1), ng_w0, &
                                     dwp(:,:,:,comp_divw0), ng_dw, lo, hi, dx)
          else
             call make_divw0_3d(w0, dwp(:,:,:,comp_divw0), ng_dw, lo, hi, dx)
          end if
       end select
    end do

  end subroutine make_divw0

  subroutine make_divw0_1d(w0,divw0,ng_dw,lo,hi,dx)

    integer, intent(in)            :: lo(:), hi(:), ng_dw
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(inout) :: divw0(lo(1)-ng_dw:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer :: i

    do i = lo(1), hi(1)
       divw0(i) = ( w0(i+1) - w0(i) )/dx(1)
    end do

  end subroutine make_divw0_1d

  subroutine make_divw0_2d(w0,divw0,ng_dw,lo,hi,dx)

    integer, intent(in)            :: lo(:), hi(:), ng_dw
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(inout) :: divw0(lo(1)-ng_dw:,lo(2)-ng_dw:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer :: j

    do j = lo(2), hi(2)
       divw0(:,j) = ( w0(j+1) - w0(j) )/dx(2)
    end do

  end subroutine make_divw0_2d

  subroutine make_divw0_3d(w0,divw0,ng_dw,lo,hi,dx)

    integer, intent(in)            :: lo(:), hi(:), ng_dw
    real(kind=dp_t), intent(in   ) :: w0(0:)
    real(kind=dp_t), intent(inout) :: divw0(lo(1)-ng_dw:,lo(2)-ng_dw:,lo(3)-ng_dw:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    !     Local variables
    integer :: k

    do k = lo(3), hi(3)
       divw0(:,:,k) = ( w0(k+1) - w0(k) )/dx(3)
    end do

  end subroutine make_divw0_3d

  subroutine make_divw0_3d_sphr(w0macx,w0macy,w0macz,ng_w0,divw0,ng_dw,lo,hi,dx)

    integer, intent(in)               :: lo(:), hi(:), ng_w0, ng_dw
    real (kind = dp_t), intent(inout) :: w0macx(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real (kind = dp_t), intent(inout) :: w0macy(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real (kind = dp_t), intent(inout) :: w0macz(lo(1)-ng_w0:,lo(2)-ng_w0:,lo(3)-ng_w0:)
    real (kind = dp_t), intent(inout) ::  divw0(lo(1)-ng_dw:,lo(2)-ng_dw:,lo(3)-ng_dw:)
    real(kind=dp_t)   , intent(in   ) :: dx(:)

    ! Local variables
    integer :: i,j,k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             divw0(i,j,k) = (w0macx(i+1,j,k)-w0macx(i,j,k)) / dx(1) &
                           +(w0macy(i,j+1,k)-w0macy(i,j,k)) / dx(2) &
                           +(w0macz(i,j,k+1)-w0macz(i,j,k)) / dx(3)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

  end subroutine make_divw0_3d_sphr


  !---------------------------------------------------------------------------
  ! make_vorticity
  !---------------------------------------------------------------------------
  subroutine make_vorticity(vort,comp,u,dx,bc)

    use bl_prof_module

    integer        , intent(in   ) :: comp
    type(multifab) , intent(inout) :: vort
    type(multifab) , intent(in   ) :: u
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(bc_level) , intent(in   ) :: bc

    real(kind=dp_t), pointer:: up(:,:,:,:)
    real(kind=dp_t), pointer:: vp(:,:,:,:)
    integer :: lo(get_dim(vort)),hi(get_dim(vort))
    integer :: i,ng_u,ng_v,dm

    type(bl_prof_timer), save :: bpt

    call build(bpt, "make_vort")

    dm = get_dim(vort)

    ng_u = nghost(u)
    ng_v = nghost(vort)

    do i = 1, nboxes(u)
       if ( multifab_remote(u, i) ) cycle
       up => dataptr(u, i)
       vp => dataptr(vort, i)
       lo =  lwb(get_box(u, i))
       hi =  upb(get_box(u, i))
       select case (dm)
       case (1)
          call setval(vort,0.d0,comp,1,all=.true.)
       case (2)
          call makevort_2d(vp(:,:,1,comp),ng_v,up(:,:,1,:),ng_u,lo,hi,dx, &
                           bc%phys_bc_level_array(i,:,:))
       case (3)
          call makevort_3d(vp(:,:,:,comp),ng_v,up(:,:,:,:),ng_u,lo,hi,dx, &
                           bc%phys_bc_level_array(i,:,:))
       end select
    end do

    call destroy(bpt)

  end subroutine make_vorticity

  subroutine makevort_2d(vort,ng_v,u,ng_u,lo,hi,dx,bc)

    use bc_module
    use bl_constants_module

    integer           , intent(in   ) :: lo(:), hi(:), ng_v, ng_u
    real (kind = dp_t), intent(  out) :: vort(lo(1)-ng_v:,lo(2)-ng_v:)  
    real (kind = dp_t), intent(in   ) ::    u(lo(1)-ng_u:,lo(2)-ng_u:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    integer           , intent(in   ) :: bc(:,:)

    !     Local variables
    integer :: i, j
    real (kind = dp_t) :: vx,uy

    do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          vx = (u(i+1,j,2) - u(i-1,j,2)) / (2.d0*dx(1)) 
          uy = (u(i,j+1,1) - u(i,j-1,1)) / (2.d0*dx(2))
          vort(i,j) = vx - uy
       enddo
    enddo

    if (bc(1,1) .eq. INLET .or. bc(1,1) .eq. SLIP_WALL .or. bc(1,1) .eq. NO_SLIP_WALL) then
       i = lo(1)
       do j = lo(2), hi(2)
          vx = (u(i+1,j,2) + 3.d0*u(i,j,2) - 4.d0*u(i-1,j,2)) / dx(1)
          uy = (u(i,j+1,1) - u(i,j-1,1)) / (2.d0*dx(2))
          vort(i,j) = vx - uy
       end do
    end if

    if (bc(1,2) .eq. INLET .or. bc(1,2) .eq. SLIP_WALL .or. bc(1,2) .eq. NO_SLIP_WALL) then
       i = hi(1)
       do j = lo(2), hi(2)
          vx = -(u(i-1,j,2) + 3.d0*u(i,j,2) - 4.d0*u(i+1,j,2)) / dx(1)
          uy = (u(i,j+1,1) - u(i,j-1,1)) / (2.d0*dx(2))
          vort(i,j) = vx - uy
       end do
    end if

    if (bc(2,1) .eq. INLET .or. bc(2,1) .eq. SLIP_WALL .or. bc(2,1) .eq. NO_SLIP_WALL) then
       j = lo(2)
       do i = lo(1), hi(1)
          vx = (u(i+1,j,2) - u(i-1,j,2)) / (2.d0*dx(1)) 
          uy = (u(i,j+1,1) + 3.d0*u(i,j,1) - 4.d0*u(i,j-1,1)) / dx(2)
          vort(i,j) = vx - uy
       end do
    end if

    if (bc(2,2) .eq. INLET .or. bc(2,2) .eq. SLIP_WALL .or. bc(2,2) .eq. NO_SLIP_WALL) then
       j = hi(2)
       do i = lo(1), hi(1)
          vx =  (u(i+1,j,2) - u(i-1,j,2)) / (2.d0*dx(1)) 
          uy = -(u(i,j-1,1) + 3.d0*u(i,j,1) - 4.d0*u(i,j+1,1)) / dx(2)
          vort(i,j) = vx - uy
       end do
    end if

  end subroutine makevort_2d

  subroutine makevort_3d(vort,ng_v,u,ng_u,lo,hi,dx,bc)

    use bc_module
    use bl_constants_module

    integer           , intent(in   ) :: lo(:), hi(:), ng_v, ng_u
    real (kind = dp_t), intent(  out) :: vort(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:)
    real (kind = dp_t), intent(in   ) ::    u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)  
    real (kind = dp_t), intent(in   ) :: dx(:)
    integer           , intent(in   ) :: bc(:,:)

    !     Local variables
    integer :: i, j, k
    logical :: fix_lo_x,fix_hi_x,fix_lo_y,fix_hi_y,fix_lo_z,fix_hi_z
    real (kind = dp_t) :: wy,vz,uz,wx,vx,uy

    !$OMP PARALLEL DO PRIVATE(i,j,k,uy,uz,vx,vz,wx,wy)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             uy = uycen(i,j,k)
             uz = uzcen(i,j,k)
             vx = vxcen(i,j,k)
             vz = vzcen(i,j,k)
             wx = wxcen(i,j,k)
             wy = wycen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

    fix_lo_x = ( bc(1,1) .eq. INLET .or. bc(1,1) .eq. NO_SLIP_WALL )
    fix_hi_x = ( bc(1,2) .eq. INLET .or. bc(1,2) .eq. NO_SLIP_WALL )

    fix_lo_y = ( bc(2,1) .eq. INLET .or. bc(2,1) .eq. NO_SLIP_WALL )
    fix_hi_y = ( bc(2,2) .eq. INLET .or. bc(2,2) .eq. NO_SLIP_WALL )

    fix_lo_z = ( bc(3,1) .eq. INLET .or. bc(3,1) .eq. NO_SLIP_WALL )
    fix_hi_z = ( bc(3,2) .eq. INLET .or. bc(3,2) .eq. NO_SLIP_WALL )

    !
    !     First do all the faces
    !
    if (fix_lo_x) then
       i = lo(1)
       !$OMP PARALLEL DO PRIVATE(j,k,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(i)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             vx = vxlo(i,j,k)
             wx = wxlo(i,j,k)
             uy = uycen(i,j,k)
             wy = wycen(i,j,k)
             uz = uzcen(i,j,k)
             vz = vzcen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_hi_x) then
       i = hi(1)
       !$OMP PARALLEL DO PRIVATE(j,k,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(i)
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             vx = vxhi(i,j,k)
             wx = wxhi(i,j,k)
             uy = uycen(i,j,k)
             wy = wycen(i,j,k)
             uz = uzcen(i,j,k)
             vz = vzcen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_lo_y) then
       j = lo(2)
       !$OMP PARALLEL DO PRIVATE(i,k,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(j)
       do k = lo(3),hi(3)
          do i = lo(1),hi(1)
             vx = vxcen(i,j,k)
             wx = wxcen(i,j,k)
             uy = uylo(i,j,k)
             wy = wylo(i,j,k)
             uz = uzcen(i,j,k)
             vz = vzcen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_hi_y) then
       j = hi(2)
       !$OMP PARALLEL DO PRIVATE(i,k,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(j)
       do k = lo(3),hi(3)
          do i = lo(1),hi(1)
             vx = vxcen(i,j,k)
             wx = wxcen(i,j,k)
             uy = uyhi(i,j,k)
             wy = wyhi(i,j,k)
             uz = uzcen(i,j,k)
             vz = vzcen(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_lo_z) then
       k = lo(3)
       !$OMP PARALLEL DO PRIVATE(i,j,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(k)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             vx = vxcen(i,j,k)
             wx = wxcen(i,j,k)
             uy = uycen(i,j,k)
             wy = wycen(i,j,k)
             uz = uzlo(i,j,k)
             vz = vzlo(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_hi_z) then
       k = hi(3)
       !$OMP PARALLEL DO PRIVATE(i,j,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(k)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             vx = vxcen(i,j,k)
             wx = wxcen(i,j,k)
             uy = uycen(i,j,k)
             wy = wycen(i,j,k)
             uz = uzhi(i,j,k)
             vz = vzhi(i,j,k)
             vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
          end do
       end do
       !$OMP END PARALLEL DO
    end if
    !
    !     Next do all the edges
    !
    if (fix_lo_x .and. fix_lo_y) then
       i = lo(1)
       j = lo(2)
       !$OMP PARALLEL DO PRIVATE(k,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(i,j)
       do k = lo(3),hi(3)
          vx = vxlo(i,j,k)
          wx = wxlo(i,j,k)
          uy = uylo(i,j,k)
          wy = wylo(i,j,k)
          uz = uzcen(i,j,k)
          vz = vzcen(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_hi_x .and. fix_lo_y) then
       i = hi(1)
       j = lo(2)
       !$OMP PARALLEL DO PRIVATE(k,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(i,j)
       do k = lo(3),hi(3)
          vx = vxhi(i,j,k)
          wx = wxhi(i,j,k)
          uy = uylo(i,j,k)
          wy = wylo(i,j,k)
          uz = uzcen(i,j,k)
          vz = vzcen(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_lo_x .and. fix_hi_y) then
       i = lo(1)
       j = hi(2)
       !$OMP PARALLEL DO PRIVATE(k,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(i,j)
       do k = lo(3),hi(3)
          vx = vxlo(i,j,k)
          wx = wxlo(i,j,k)
          uy = uyhi(i,j,k)
          wy = wyhi(i,j,k)
          uz = uzcen(i,j,k)
          vz = vzcen(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_hi_x .and. fix_hi_y) then
       i = hi(1)
       j = hi(2)
       !$OMP PARALLEL DO PRIVATE(k,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(i,j)
       do k = lo(3),hi(3)
          vx = vxhi(i,j,k)
          wx = wxhi(i,j,k)
          uy = uyhi(i,j,k)
          wy = wyhi(i,j,k)
          uz = uzcen(i,j,k)
          vz = vzcen(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_lo_x .and. fix_lo_z) then
       i = lo(1)
       k = lo(3)
       !$OMP PARALLEL DO PRIVATE(j,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(i,k)
       do j = lo(2),hi(2)
          vx = vxlo(i,j,k)
          wx = wxlo(i,j,k)
          uy = uycen(i,j,k)
          wy = wycen(i,j,k)
          uz = uzlo(i,j,k)
          vz = vzlo(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_hi_x .and. fix_lo_z) then
       i = hi(1)
       k = lo(3)
       !$OMP PARALLEL DO PRIVATE(j,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(i,k)
       do j = lo(2),hi(2)
          vx = vxhi(i,j,k)
          wx = wxhi(i,j,k)
          uy = uycen(i,j,k)
          wy = wycen(i,j,k)
          uz = uzlo(i,j,k)
          vz = vzlo(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_lo_x .and. fix_hi_z) then
       i = lo(1)
       k = hi(3)
       !$OMP PARALLEL DO PRIVATE(j,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(i,k)
       do j = lo(2),hi(2)
          vx = vxlo(i,j,k)
          wx = wxlo(i,j,k)
          uy = uycen(i,j,k)
          wy = wycen(i,j,k)
          uz = uzhi(i,j,k)
          vz = vzhi(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_hi_x .and. fix_hi_z) then
       i = hi(1)
       k = hi(3)
       !$OMP PARALLEL DO PRIVATE(j,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(i,k)
       do j = lo(2),hi(2)
          vx = vxhi(i,j,k)
          wx = wxhi(i,j,k)
          uy = uycen(i,j,k)
          wy = wycen(i,j,k)
          uz = uzhi(i,j,k)
          vz = vzhi(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_lo_y .and. fix_lo_z) then
       j = lo(2)
       k = lo(3)
       !$OMP PARALLEL DO PRIVATE(i,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(j,k)
       do i = lo(1),hi(1)
          vx = vxcen(i,j,k)
          wx = wxcen(i,j,k)
          uy = uylo(i,j,k)
          wy = wylo(i,j,k)
          uz = uzlo(i,j,k)
          vz = vzlo(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_hi_y .and. fix_lo_z) then
       j = hi(2)
       k = lo(3)
       !$OMP PARALLEL DO PRIVATE(i,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(j,k)
       do i = lo(1),hi(1)
          vx = vxcen(i,j,k)
          wx = wxcen(i,j,k)
          uy = uyhi(i,j,k)
          wy = wyhi(i,j,k)
          uz = uzlo(i,j,k)
          vz = vzlo(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_lo_y .and. fix_hi_z) then
       j = lo(2)
       k = hi(3)
       !$OMP PARALLEL DO PRIVATE(i,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(j,k)
       do i = lo(1),hi(1)
          vx = vxcen(i,j,k)
          wx = wxcen(i,j,k)
          uy = uylo(i,j,k)
          wy = wylo(i,j,k)
          uz = uzhi(i,j,k)
          vz = vzhi(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
       !$OMP END PARALLEL DO
    end if

    if (fix_hi_y .and. fix_hi_z) then
       j = hi(2)
       k = hi(3)
       !$OMP PARALLEL DO PRIVATE(i,uy,uz,vx,vz,wx,wy) FIRSTPRIVATE(j,k)
       do i = lo(1),hi(1)
          vx = vxcen(i,j,k)
          wx = wxcen(i,j,k)
          uy = uyhi(i,j,k)
          wy = wyhi(i,j,k)
          uz = uzhi(i,j,k)
          vz = vzhi(i,j,k)
          vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
       end do
       !$OMP END PARALLEL DO
    end if
    !
    !     Finally do all the corners
    !
    if (fix_lo_x .and. fix_lo_y .and. fix_lo_z) then
       i = lo(1)
       j = lo(2)
       k = lo(3)
       vx = vxlo(i,j,k)
       wx = wxlo(i,j,k)
       uy = uylo(i,j,k)
       wy = wylo(i,j,k)
       uz = uzlo(i,j,k)
       vz = vzlo(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_hi_x .and. fix_lo_y .and. fix_lo_z) then
       i = hi(1)
       j = lo(2)
       k = lo(3)
       vx = vxhi(i,j,k)
       wx = wxhi(i,j,k)
       uy = uylo(i,j,k)
       wy = wylo(i,j,k)
       uz = uzlo(i,j,k)
       vz = vzlo(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_lo_x .and. fix_hi_y .and. fix_lo_z) then
       i = lo(1)
       j = hi(2)
       k = lo(3)
       vx = vxlo(i,j,k)
       wx = wxlo(i,j,k)
       uy = uyhi(i,j,k)
       wy = wyhi(i,j,k)
       uz = uzlo(i,j,k)
       vz = vzlo(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_hi_x .and. fix_hi_y .and. fix_lo_z) then
       i = hi(1)
       j = hi(2)
       k = lo(3)
       vx = vxhi(i,j,k)
       wx = wxhi(i,j,k)
       uy = uyhi(i,j,k)
       wy = wyhi(i,j,k)
       uz = uzlo(i,j,k)
       vz = vzlo(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_lo_x .and. fix_lo_y .and. fix_hi_z) then
       i = lo(1)
       j = lo(2)
       k = hi(3)
       vx = vxlo(i,j,k)
       wx = wxlo(i,j,k)
       uy = uylo(i,j,k)
       wy = wylo(i,j,k)
       uz = uzhi(i,j,k)
       vz = vzhi(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_hi_x .and. fix_lo_y .and. fix_hi_z) then
       i = hi(1)
       j = lo(2)
       k = hi(3)
       vx = vxhi(i,j,k)
       wx = wxhi(i,j,k)
       uy = uylo(i,j,k)
       wy = wylo(i,j,k)
       uz = uzhi(i,j,k)
       vz = vzhi(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_lo_x .and. fix_hi_y .and. fix_hi_z) then
       i = lo(1)
       j = hi(2)
       k = hi(3)
       vx = vxlo(i,j,k)
       wx = wxlo(i,j,k)
       uy = uyhi(i,j,k)
       wy = wyhi(i,j,k)
       uz = uzhi(i,j,k)
       vz = vzhi(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

    if (fix_hi_x .and. fix_hi_y .and. fix_hi_z) then
       i = hi(1)
       j = hi(2)
       k = hi(3)
       vx = vxhi(i,j,k)
       wx = wxhi(i,j,k)
       uy = uyhi(i,j,k)
       wy = wyhi(i,j,k)
       uz = uzhi(i,j,k)
       vz = vzhi(i,j,k)
       vort(i,j,k) = vorfun(uy,uz,vx,vz,wx,wy)
    end if

  contains

    function uycen(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i,j+1,k,1)-u(i,j-1,k,1))/dx(2)
    end function uycen

    function uylo(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i,j+1,k,1)+THREE*u(i,j,k,1)-FOUR*u(i,j-1,k,1))/(THREE*dx(2))
    end function uylo

    function uyhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = -(u(i,j-1,k,1)+THREE*u(i,j,k,1)-FOUR*u(i,j+1,k,1))/(THREE*dx(2))
    end function uyhi

    function uzcen(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i,j,k+1,1)-u(i,j,k-1,1))/dx(3)
    end function uzcen

    function uzlo(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i,j,k+1,1)+THREE*u(i,j,k,1)-FOUR*u(i,j,k-1,1))/(THREE*dx(3))
    end function uzlo

    function uzhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r =-(u(i,j,k-1,1)+THREE*u(i,j,k,1)-FOUR*u(i,j,k+1,1))/(THREE*dx(3))
    end function uzhi

    function vxcen(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i+1,j,k,2)-u(i-1,j,k,2))/dx(1)
    end function vxcen

    function vxlo(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i+1,j,k,2)+THREE*u(i,j,k,2)-FOUR*u(i-1,j,k,2))/(THREE*dx(1))
    end function vxlo

    function vxhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r =-(u(i-1,j,k,2)+THREE*u(i,j,k,2)-FOUR*u(i+1,j,k,2))/(THREE*dx(1))
    end function vxhi

    function vzcen(i,j,k) result(r) 
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i,j,k+1,2)-u(i,j,k-1,2))/dx(3)
    end function vzcen

    function vzlo(i,j,k) result(r) 
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i,j,k+1,2)+THREE*u(i,j,k,2)-FOUR*u(i,j,k-1,2))/(THREE*dx(3))
    end function vzlo

    function vzhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r =-(u(i,j,k-1,2)+THREE*u(i,j,k,2)-FOUR*u(i,j,k+1,2))/(THREE*dx(3))
    end function vzhi

    function wxcen(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i+1,j,k,3)-u(i-1,j,k,3))/dx(1)
    end function wxcen

    function wxlo(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i+1,j,k,3)+THREE*u(i,j,k,3)-FOUR*u(i-1,j,k,3))/(THREE*dx(1))
    end function wxlo

    function wxhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r =-(u(i-1,j,k,3)+THREE*u(i,j,k,3)-FOUR*u(i+1,j,k,3))/(THREE*dx(1))
    end function wxhi

    function wycen(i,j,k) result(r) 
      integer :: i,j,k
      real(dp_t) :: r
      r = HALF*(u(i,j+1,k,3)-u(i,j-1,k,3))/dx(2)
    end function wycen

    function wylo(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r = (u(i,j+1,k,3)+THREE*u(i,j,k,3)-FOUR*u(i,j-1,k,3))/(THREE*dx(2))
    end function wylo

    function wyhi(i,j,k) result(r)
      integer :: i,j,k
      real(dp_t) :: r
      r =-(u(i,j-1,k,3)+THREE*u(i,j,k,3)-FOUR*u(i,j+1,k,3))/(THREE*dx(2))
    end function wyhi

    function vorfun(uy,uz,vx,vz,wx,wy) result(r)
      real(dp_t) :: uy,uz,vx,vz,wx,wy
      real(dp_t) :: r
      r = sqrt((wy-vz)**2+(uz-wx)**2+(vx-uy)**2)
    end function vorfun

  end subroutine makevort_3d


  !---------------------------------------------------------------------------
  ! make_magvel
  !---------------------------------------------------------------------------
  subroutine make_magvel(plotdata,comp_magvel,comp_mom,s,u,w0,w0mac)

    use bc_module
    use bl_constants_module
    use geometry, only : spherical
    use variables, only : rho_comp

    integer        , intent(in   ) :: comp_magvel, comp_mom
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: s, u
    real(kind=dp_t), intent(in   ) :: w0(0:)
    type(multifab) , intent(in   ) :: w0mac(:)

    real(kind=dp_t), pointer:: pp(:,:,:,:)
    real(kind=dp_t), pointer:: sp(:,:,:,:)
    real(kind=dp_t), pointer:: up(:,:,:,:)
    real(kind=dp_t), pointer:: wxp(:,:,:,:)
    real(kind=dp_t), pointer:: wyp(:,:,:,:)
    real(kind=dp_t), pointer:: wzp(:,:,:,:)

    integer :: lo(get_dim(plotdata)),hi(get_dim(plotdata)),ng_p,ng_s,ng_u,ng_w
    integer :: i,dm

    dm = get_dim(plotdata)

    ng_s = nghost(s)
    ng_u = nghost(u)
    ng_p = nghost(plotdata)

    do i = 1, nboxes(u)
       if ( multifab_remote(u, i) ) cycle
       pp => dataptr(plotdata, i)
       sp => dataptr(s, i)
       up => dataptr(u, i)
       lo =  lwb(get_box(u, i))
       hi =  upb(get_box(u, i))
       select case (dm)
       case (1)
          call makemagvel_1d(pp(:,1,1,comp_magvel),pp(:,1,1,comp_mom),ng_p, &
                             sp(:,1,1,rho_comp),ng_s,up(:,1,1,1),ng_u,w0,lo,hi)
       case (2)
          call makemagvel_2d(pp(:,:,1,comp_magvel),pp(:,:,1,comp_mom),ng_p, &
                             sp(:,:,1,rho_comp),ng_s,up(:,:,1,:),ng_u,w0,lo,hi)
       case (3)
          if (spherical .eq. 1) then
             wxp => dataptr(w0mac(1), i)
             wyp => dataptr(w0mac(2), i)
             wzp => dataptr(w0mac(3), i)
             ng_w = nghost(w0mac(1))
             call makemagvel_3d_sphr(pp(:,:,:,comp_magvel),pp(:,:,:,comp_mom),ng_p, &
                                     sp(:,:,:,rho_comp),ng_s,up(:,:,:,:),ng_u, &
                                     wxp(:,:,:,1),wyp(:,:,:,1),wzp(:,:,:,1),ng_w,lo,hi)
          else
             call makemagvel_3d_cart(pp(:,:,:,comp_magvel),pp(:,:,:,comp_mom),ng_p, &
                                     sp(:,:,:,rho_comp),ng_s,up(:,:,:,:),ng_u,w0,lo,hi)
          end if
       end select
    end do

  end subroutine make_magvel

  subroutine makemagvel_1d(magvel,mom,ng_p,rho,ng_s,u,ng_u,w0,lo,hi)

    integer           , intent(in   ) :: lo(:), hi(:), ng_p, ng_u, ng_s
    real (kind = dp_t), intent(  out) :: magvel(lo(1)-ng_p:)
    real (kind = dp_t), intent(  out) ::    mom(lo(1)-ng_p:)
    real (kind = dp_t), intent(in   ) ::    rho(lo(1)-ng_s:)
    real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u:)
    real (kind = dp_t), intent(in   ) :: w0(0:)

    !     Local variables
    integer :: i
    real (kind = dp_t) :: w0_cent

    ! Recall w0 is edge-centered
    do i = lo(1), hi(1)
       w0_cent = 0.5d0 * (w0(i) + w0(i+1))
       magvel(i) = abs(u(i)+w0_cent)
       mom(i) = rho(i)*magvel(i)
    enddo

  end subroutine makemagvel_1d

  subroutine makemagvel_2d(magvel,mom,ng_p,rho,ng_s,u,ng_u,w0,lo,hi)

    integer           , intent(in   ) :: lo(:), hi(:), ng_p, ng_u, ng_s
    real (kind = dp_t), intent(  out) :: magvel(lo(1)-ng_p:,lo(2)-ng_p:)
    real (kind = dp_t), intent(  out) ::    mom(lo(1)-ng_p:,lo(2)-ng_p:)
    real (kind = dp_t), intent(in   ) ::    rho(lo(1)-ng_s:,lo(2)-ng_s:)  
    real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,:)  
    real (kind = dp_t), intent(in   ) :: w0(0:)

    !     Local variables
    integer :: i, j
    real (kind = dp_t) :: w0_cent

    ! Recall w0 is edge-centered
    do j = lo(2), hi(2)
       w0_cent = 0.5d0 * (w0(j) + w0(j+1))
       do i = lo(1), hi(1)
          magvel(i,j) = sqrt( u(i,j,1)**2 + (u(i,j,2)+w0_cent)**2 )
          mom(i,j) = rho(i,j)*magvel(i,j)
       enddo
    enddo

  end subroutine makemagvel_2d

  subroutine makemagvel_3d_cart(magvel,mom,ng_p,rho,ng_s,u,ng_u,w0,lo,hi)

    use geometry, only : spherical

    integer           , intent(in   ) :: lo(:), hi(:), ng_p, ng_u, ng_s
    real (kind = dp_t), intent(  out) :: magvel(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real (kind = dp_t), intent(  out) ::    mom(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real (kind = dp_t), intent(in   ) ::    rho(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:) 
    real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:) 
    real (kind = dp_t), intent(in   ) :: w0(0:)

    !     Local variables
    integer :: i, j, k
    real (kind = dp_t) :: w0_cent

    ! Recall w0 is edge-centered
    do k = lo(3), hi(3)
       w0_cent = 0.5d0 * (w0(k) + w0(k+1))
       do j = lo(2), hi(2)
       do i = lo(1), hi(1)
          magvel(i,j,k) = sqrt(u(i,j,k,1)**2 + u(i,j,k,2)**2 + (u(i,j,k,3)+w0_cent)**2)
          mom(i,j,k) = rho(i,j,k)*magvel(i,j,k)
       enddo
       enddo
    enddo

  end subroutine makemagvel_3d_cart

  subroutine makemagvel_3d_sphr(magvel,mom,ng_p,rho,ng_s,u,ng_u, &
                                w0macx,w0macy,w0macz,ng_w,lo,hi)


    use bl_constants_module

    integer           , intent(in   ) :: lo(:), hi(:), ng_p, ng_u, ng_w, ng_s
    real (kind = dp_t), intent(  out) :: magvel(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real (kind = dp_t), intent(  out) ::    mom(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real (kind = dp_t), intent(in   ) ::    rho(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:) 
    real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:) 
    real (kind = dp_t), intent(in   ) :: w0macx(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:)
    real (kind = dp_t), intent(in   ) :: w0macy(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:)
    real (kind = dp_t), intent(in   ) :: w0macz(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:)

    !     Local variables
    integer :: i, j, k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             magvel(i,j,k) = sqrt( (u(i,j,k,1)+HALF*(w0macx(i,j,k)+w0macx(i+1,j,k)))**2 + &
                                   (u(i,j,k,2)+HALF*(w0macy(i,j,k)+w0macy(i,j+1,k)))**2 + &
                                   (u(i,j,k,3)+HALF*(w0macz(i,j,k)+w0macz(i,j,k+1)))**2)
             mom(i,j,k) = rho(i,j,k)*magvel(i,j,k)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine makemagvel_3d_sphr


  !---------------------------------------------------------------------------
  ! make_velrc
  !---------------------------------------------------------------------------
  subroutine make_velrc(plotdata,comp_velr,comp_velc,u,w0r_cart,normal)

    use bc_module
    use bl_constants_module
    use geometry, only: spherical

    integer        , intent(in   ) :: comp_velr, comp_velc
    type(multifab) , intent(inout) :: plotdata
    type(multifab) , intent(in   ) :: u
    type(multifab) , intent(in   ) :: w0r_cart
    type(multifab) , intent(in   ) :: normal

    ! local
    real(kind=dp_t), pointer:: pp(:,:,:,:)
    real(kind=dp_t), pointer:: up(:,:,:,:)
    real(kind=dp_t), pointer:: nop(:,:,:,:)
    real(kind=dp_t), pointer:: w0rp(:,:,:,:)
    integer :: lo(get_dim(plotdata)),hi(get_dim(plotdata)),ng_p,ng_u,ng_n,ng_w
    integer :: i,dm

    dm = get_dim(plotdata)

    ng_u = nghost(u)
    ng_p = nghost(plotdata)
    ng_n = nghost(normal)
    ng_w = nghost(w0r_cart)

    if (spherical .ne. 1) then
       call bl_error("unable to create radial and circumferential velocity -- not spherical geometry")
    endif

    do i = 1, nboxes(u)

       if ( multifab_remote(u, i) ) cycle

       pp => dataptr(plotdata, i)
       up => dataptr(u, i)
       nop => dataptr(normal, i)
       w0rp => dataptr(w0r_cart, i)
       lo =  lwb(get_box(u, i))
       hi =  upb(get_box(u, i))

       call makevelrc_3d_sphr(pp(:,:,:,comp_velr),pp(:,:,:,comp_velc),&
                              ng_p,up(:,:,:,:),ng_u, &
                              w0rp(:,:,:,1),ng_w,nop(:,:,:,:),ng_n,lo,hi)
    end do

  end subroutine make_velrc

  subroutine makevelrc_3d_sphr(velr,velc,ng_p,u,ng_u,w0r,ng_w,normal,ng_n,lo,hi)

    integer           , intent(in   ) :: lo(:), hi(:), ng_p, ng_u, ng_n, ng_w
    real (kind = dp_t), intent(  out) ::   velr(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real (kind = dp_t), intent(  out) ::   velc(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
    real (kind = dp_t), intent(in   ) ::      u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
    real (kind = dp_t), intent(in   ) ::    w0r(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:)
    real (kind = dp_t), intent(in   ) :: normal(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:,:)  

    !     Local variables
    integer :: i, j, k

    !$OMP PARALLEL DO PRIVATE(i,j,k)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             velr(i,j,k) = u(i,j,k,1)*normal(i,j,k,1) + &
                           u(i,j,k,2)*normal(i,j,k,2) + &
                           u(i,j,k,3)*normal(i,j,k,3) 

             velc(i,j,k) = (u(i,j,k,1)-velr(i,j,k)*normal(i,j,k,1)) * &
                           (u(i,j,k,1)-velr(i,j,k)*normal(i,j,k,1))
             velc(i,j,k) = velc(i,j,k) + &
                           (u(i,j,k,2)-velr(i,j,k)*normal(i,j,k,2)) * &
                           (u(i,j,k,2)-velr(i,j,k)*normal(i,j,k,2))
             velc(i,j,k) = velc(i,j,k) + &
                           (u(i,j,k,3)-velr(i,j,k)*normal(i,j,k,3)) * &
                           (u(i,j,k,3)-velr(i,j,k)*normal(i,j,k,3))
             velc(i,j,k) = sqrt(velc(i,j,k))

             velr(i,j,k) = velr(i,j,k) + w0r(i,j,k)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

  end subroutine makevelrc_3d_sphr

end module plot_variables_module

